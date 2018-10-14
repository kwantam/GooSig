#!/usr/bin/python
#
# (C) 2018 Riad S. Wahby <rsw@cs.stanford.edu>

import itertools
from itertools import islice
import sys

from libGooPy.defs import Defs
import libGooPy.util as lutil

# python 2/3 hack
if sys.version_info[0] == 2:
    zip = itertools.izip    # pylint: disable=redefined-builtin,no-member
    range = xrange          # pylint: disable=redefined-builtin,undefined-variable

class RSAGroupOps(object):
    # NOTE you should use an RSA modulus whose factorization is unknown!
    #      In other words, *don't* just generate a modulus! defs.py provides
    #      a few candidates for you to try.
    def __init__(self, Gdesc, modbits=None, prng=None):
        self.n = Gdesc.modulus
        self.nOver2 = self.n // 2
        self.g = Gdesc.g
        self.h = Gdesc.h

        # random scalars can cut off a bit without any danger; this makes rand_scalar faster
        self.nbits_rand = lutil.clog2(self.n) - 1

        # combs for g and h
        # NOTE: you really want to store all combs on disk!
        #       I'd recommend having a telescope of combs supporting up to (say)
        #       8192-bit RSA keys (i.e., a ~(2 * 8192 + chalbits + 1) sized comb)
        #       sized for the group of unknown order (see big_prod_size below)
        #
        # Only P needs "big" and "small" combs; V just needs "tiny"
        #
        combsize = Defs.combsize
        self.combsize = combsize
        # figure out comb sizes
        if modbits is not None:
            # largest exponent P handles is the greater of
            #       chalbits + 2 * log2(P's RSA modulus) + 1
            # and
            #       chalbits + log2(P's RSA modulus) + log2(n) + 1
            big_prod_size = max(2 * modbits, modbits + lutil.clog2(self.n))
            self.big_comb = combsize * ((big_prod_size + Defs.chalbits + combsize) // combsize)
            # P has to do a bunch of arithmetic with exponents of size log2(|G|) for G this group
            self.small_comb = combsize * ((lutil.clog2(self.n) + combsize - 1) // combsize)
            self.tiny_comb = -1
        else:
            self.big_comb = self.small_comb = -1
            # V's work is over tiny moduli, just chalbits long
            self.tiny_comb = combsize * ((Defs.chalbits + 1 + combsize) // combsize)

        # compute the combs
        if modbits is not None:
            self.gcomb_big = self._precomp_comb(self.g, self.big_comb)
            self.hcomb_big = self._precomp_comb(self.h, self.big_comb)
            self.gcomb_small = self._precomp_comb(self.g, self.small_comb)
            self.hcomb_small = self._precomp_comb(self.h, self.small_comb)
            self.gcomb_tiny = self.hcomb_tiny = None
        else:
            self.gcomb_big = self.gcomb_small = self.hcomb_big = self.hcomb_small = None
            self.gcomb_tiny = self._precomp_comb(self.g, self.tiny_comb)
            self.hcomb_tiny = self._precomp_comb(self.h, self.tiny_comb)

        if prng is None:
            # /dev/urandom
            self.prng = lutil.rand
        else:
            # provided RNG
            self.prng = prng

    def quot(self, b):
        # compute the representative of (Z/n)/{1,-1}, i.e., min(|b|, n-|b|)
        if b > self.nOver2:
            return self.n - b
        return b

    def is_quot(self, b):
        return b <= self.nOver2

    def pow(self, b, e):
        # native pow is pretty fast already
        return self.quot(pow(b, e, self.n))

    def _precomp_wind2(self, b1, b2, winsize):
        ret = [1] * (2 ** (2 * winsize))
        winlen = 2 ** winsize
        for i in range(0, winlen):
            if i != 0:
                ret[i] = (ret[i-1] * b1) % self.n
            for j in range(1, winlen):
                ret[i + winlen * j] = (ret[i + winlen * (j - 1)] * b2) % self.n
        return ret

    def pow2(self, b1, e1, b2, e2):
        pctable = self._precomp_wind2(b1, b2, Defs.winsize)

        # exponents as bits, padded to multiple of winsize
        nwins = (max(e1.bit_length(), e2.bit_length()) + (Defs.winsize - 1)) // Defs.winsize
        e1bits = lutil.num_to_bits(e1, nwins * Defs.winsize)
        e2bits = lutil.num_to_bits(e2, nwins * Defs.winsize)

        winlen = 2 ** Defs.winsize
        ret = 1
        for win in range(0, nwins):
            if win > 0:
                ret = pow(ret, winlen, self.n)

            e1val = self._from_win(islice(e1bits, Defs.winsize))
            e2val = self._from_win(islice(e2bits, Defs.winsize))

            ret *= pctable[e1val + winlen * e2val]
            ret %= self.n

        return self.quot(ret)

    @staticmethod
    def _from_win(vs):
        ret = 0
        for v in vs:
            ret <<= 1
            ret += 1 if v else 0
        return ret

    def _precomp_comb(self, b, max_comb):
        # first, make a "basic" comb
        pow_bpc = 2 ** (max_comb // self.combsize)
        comb_basic = [None] * self.combsize
        comb_basic[0] = b
        for i in range(1, self.combsize):
            comb_basic[i] = pow(comb_basic[i-1], pow_bpc, self.n)

        # next, make a "windowed" comb
        comb_windowed = [None] * (2 ** self.combsize)
        comb_windowed[0] = 1
        for i in range(0, self.combsize):
            p_off = 1 << i
            comb_windowed[p_off] = comb_basic[i]
            for j in range(p_off + 1, 2 * p_off):
                comb_windowed[j] = (comb_windowed[j - p_off] * comb_windowed[p_off]) % self.n

        return comb_windowed

    def _to_comb_bits(self, e, max_comb):
        ebits = list(lutil.num_to_bits(e, max_comb))
        bpc = max_comb // self.combsize
        comb_words_iter = zip(*( ebits[j*bpc:(j+1)*bpc] for j in range(0, self.combsize) ))
        return ( self._from_win(x) for x in comb_words_iter )

    def powgh(self, e1, e2):
        loge = max(int(e1).bit_length(), int(e2).bit_length())
        if loge <= self.tiny_comb:
            max_comb = self.tiny_comb
            gcomb = self.gcomb_tiny
            hcomb = self.hcomb_tiny
        elif loge <= self.small_comb:
            max_comb = self.small_comb
            gcomb = self.gcomb_small
            hcomb = self.hcomb_small
        elif loge <= self.big_comb:
            max_comb = self.big_comb
            gcomb = self.gcomb_big
            hcomb = self.hcomb_big
        else:
            raise ValueError("got unexpectedly large exponent in powgh (%d > %d)" % (loge, self.big_comb))

        e1bits = self._to_comb_bits(e1, max_comb)
        e2bits = self._to_comb_bits(e2, max_comb)

        ret = 1
        for (e1v, e2v) in zip(e1bits, e2bits):
            ret = pow(ret, 2, self.n)
            ret = (ret * gcomb[e1v]) % self.n
            ret = (ret * hcomb[e2v]) % self.n

        return self.quot(ret)

    def mul(self, m1, m2):
        return self.quot((m1 * m2) % self.n)

    def inv2(self, b1, b2):
        b12Inv = lutil.invert_modp(b1 * b2, self.n)
        return (self.quot((b2 * b12Inv) % self.n), self.quot((b1 * b12Inv) % self.n))

    def rand_scalar(self):
        return self.prng.getrandbits(self.nbits_rand)

def main(nreps):
    import libGooPy.test_util as tu

    # test on random RSA modulus
    (p, q) = lutil.rand.sample(Defs.primes_1024, 2)
    n = p * q
    Grandom = Defs.gen_group_obj(n, 5, 7)

    t1 = RSAGroupOps(Defs.Grsa2048, 2048)
    t2 = RSAGroupOps(Grandom, 2048)

    def test_pow2():
        "pow2,RSA_chal,RSA_rand"

        (b1, b2, e1, e2) = ( lutil.rand.getrandbits(2048) for _ in range(0, 4) )
        out1 = t1.quot((pow(b1, e1, Defs.Grsa2048.modulus) * pow(b2, e2, Defs.Grsa2048.modulus)) % Defs.Grsa2048.modulus)
        t1o = t1.pow2(b1, e1, b2, e2)

        out2 = t2.quot((pow(b1, e1, n) * pow(b2, e2, n)) % n)
        t2o = t2.pow2(b1, e1, b2, e2)

        return (out1 == t1o, out2 == t2o)

    def test_powgh():
        "powgh,RSA_chal,RSA_rand"

        (e1, e2) = ( lutil.rand.getrandbits(2 * 2048 + Defs.chalbits + 2) for _ in range(0, 2) )

        out1 = t1.quot((pow(2, e1, Defs.Grsa2048.modulus) * pow(3, e2, Defs.Grsa2048.modulus)) % Defs.Grsa2048.modulus)
        t1o = t1.powgh(e1, e2)

        (e1_s, e2_s) = ( x >> (2048 + Defs.chalbits) for x in (e1, e2) )
        out2 = t2.quot((pow(5, e1_s, n) * pow(7, e2_s, n)) % n)
        t2o = t2.powgh(e1_s, e2_s)

        return (out1 == t1o, out2 == t2o)

    def test_inv2():
        "inv2,RSA_chal,RSA_rand"

        (e1, e2) = ( lutil.rand.getrandbits(2048) for _ in range(0, 2) )
        (e1Inv, e2Inv) = t1.inv2(e1, e2)
        t1pass = t1.quot((e1 * e1Inv) % Defs.Grsa2048.modulus) == 1 and t1.quot((e2 * e2Inv) % Defs.Grsa2048.modulus) == 1

        (e1_s, e2_s) = ( x >> 1536 for x in (e1, e2) )
        (e1_sInv, e2_sInv) = t2.inv2(e1_s, e2_s)
        t2pass = t2.quot((e1_s * e1_sInv) % n) == 1 and t2.quot((e2_s * e2_sInv) % n) == 1

        return (t1pass, t2pass)

    tu.run_all_tests(nreps, "group_ops", test_pow2, test_powgh, test_inv2)

if __name__ == "__main__":
    try:
        nr = int(sys.argv[1])
    except:
        nr = 32
    main(nr)
