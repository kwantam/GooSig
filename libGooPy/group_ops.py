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

def _from_win(vs):
    ret = 0
    for v in vs:
        ret <<= 1
        ret += 1 if v else 0
    return ret

def _wnaf(r, w, bitlen=None):
    if bitlen is None:
        bitlen = r.bit_length() + 1
    out = [None] * bitlen
    for i in reversed(range(0, bitlen)):
        val = 0
        if r % 2:
            val = r & ((1<<w)-1)
            if val & (1<<(w-1)):
                val -= (1<<w)
            r -= val
        out[i] = val
        r = r >> 1
    assert r == 0
    return out

class _CombPrecomp(list):
    def __init__(self, g, combspec, gops):
        (ppa, aps, nshifts, bpw, _, size) = combspec
        self.gops = gops
        self.points_per_add = ppa
        self.adds_per_shift = aps
        self.nshifts = nshifts
        self.bits_per_window = bpw
        self.nbits = bpw * ppa
        nskip = 2 ** ppa - 1
        self.points_per_subcomb = nskip

        # allocate space
        list.__init__(self, [None] * size)

        # compute bottom comb
        self[0] = g
        for i in range(1, ppa):
            oval = 1 << i
            ival = oval // 2
            self[oval - 1] = gops.pow(self[ival - 1], 1 << bpw)
            for j in range(oval + 1, 2 * oval):
                self[j - 1] = gops.mul(self[j - oval - 1], self[oval - 1])

        # compute upper combs from bottom
        powval = 1 << nshifts
        for i in range(1, aps):
            for j in range(0, nskip):
                self[i * nskip + j] = gops.pow(self[i * nskip + j - nskip], powval)

    def to_comb_exp(self, e):
        ebits = list(lutil.num_to_bits(e, self.nbits))
        nsh = self.nshifts
        nsh_tot = self.adds_per_shift * self.points_per_add
        ebits_split = [ ebits[j*nsh:(j+1)*nsh] for j in range(0, nsh_tot) ]
        iters = []
        for i in reversed(range(0, self.adds_per_shift)):
            iters.append( _from_win(x) for x in zip(*( ebits_split[i + j*self.adds_per_shift] for j in range(0, self.points_per_add) )) )
        return zip(*iters)

    @staticmethod
    def gen_opt_combs(nbits, maxsize=None):
        # an "optimal" comb for a given #ops is the one that uses the least storage
        #                   for a given storage size is the one that uses the least ops
        opt_combs = {}

        def _gen_comb_result(nshifts, adds_per_shift, points_per_add, bits_per_window):
            nops = nshifts * (adds_per_shift + 1) - 1
            size = (2 ** points_per_add - 1) * adds_per_shift
            result = (points_per_add, adds_per_shift, nshifts, bits_per_window, nops, size)

            best_so_far = opt_combs.get(nops, None)
            if best_so_far is None or best_so_far[5] > size:
                opt_combs[nops] = result

        for points_per_add in range(2, 18):
            bits_per_window = (nbits + points_per_add - 1) // points_per_add
            for adds_per_shift in range(1, lutil.isqrt(bits_per_window) + 2):
                if bits_per_window % adds_per_shift != 0:
                    # only factorizations of bits_per_window are useful
                    continue
                nshifts = bits_per_window // adds_per_shift

                _gen_comb_result(nshifts, adds_per_shift, points_per_add, bits_per_window)
                _gen_comb_result(adds_per_shift, nshifts, points_per_add, bits_per_window)

        ret_all = []
        ret = None
        sm = None
        for nops in sorted(opt_combs.keys()):
            opt_comb_val = opt_combs[nops]

            if sm is not None and sm <= opt_comb_val[5]:
                continue

            sm = opt_comb_val[5]
            ret_all.append(opt_comb_val)

            if ret is None and maxsize is not None and opt_comb_val[5] <= maxsize:
                ret = opt_comb_val
                break

        if maxsize is None:
            return ret_all
        return ret

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
        # P needs two comb sizes, V needs one
        #
        # figure out comb sizes
        if modbits is not None:
            # largest exponent P handles is the greater of
            #       chalbits + 2 * log2(P's RSA modulus) + 1
            #       chalbits + log2(P's RSA modulus) + log2(n) + 1
            big_nbits = max(2 * modbits, modbits + lutil.clog2(self.n)) + Defs.chalbits + 1
            big_combspec = _CombPrecomp.gen_opt_combs(big_nbits, Defs.max_comb_size)
            small_nbits = lutil.clog2(self.n)
            small_combspec = _CombPrecomp.gen_opt_combs(small_nbits, Defs.max_comb_size)
            self.combs = [(_CombPrecomp(self.g, small_combspec, self), _CombPrecomp(self.h, small_combspec, self))
                         ,(_CombPrecomp(self.g, big_combspec, self), _CombPrecomp(self.h, big_combspec, self))
                         ]
        else:
            tiny_nbits = Defs.chalbits
            tiny_combspec = _CombPrecomp.gen_opt_combs(tiny_nbits, Defs.max_comb_size)
            self.combs = [(_CombPrecomp(self.g, tiny_combspec, self), _CombPrecomp(self.h, tiny_combspec, self))]

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

    def _precomp_wnaf(self, b, bInv, winsize):
        tablen = 1 << (winsize - 2)
        pctabP = [None] * tablen
        pctabN = [None] * tablen
        bSq = pow(b, 2, self.n)
        bInvSq = pow(bInv, 2, self.n)
        pctabP[0] = b
        pctabN[0] = bInv
        for i in range(1, tablen):
            pctabP[i] = (pctabP[i-1] * bSq) % self.n
            pctabN[i] = (pctabN[i-1] * bInvSq) % self.n
        return (pctabP, pctabN)

    def pow2_wnaf(self, b1, b1Inv, e1, b2, b2Inv, e2):
        (pctabP1, pctabN1) = self._precomp_wnaf(b1, b1Inv, Defs.winsize)
        (pctabP2, pctabN2) = self._precomp_wnaf(b2, b2Inv, Defs.winsize)

        totlen = max(e1.bit_length(), e2.bit_length()) + 1
        e1bits = _wnaf(e1, Defs.winsize, totlen)
        e2bits = _wnaf(e2, Defs.winsize, totlen)

        ret = 1
        for (w1, w2) in zip(e1bits, e2bits):
            if ret != 1:
                ret = pow(ret, 2, self.n)

            if w1 > 0:
                ret = (ret * pctabP1[(w1-1)//2]) % self.n
            elif w1 < 0:
                ret = (ret * pctabN1[(-1-w1)//2]) % self.n

            if w2 > 0:
                ret = (ret * pctabP2[(w2-1)//2]) % self.n
            elif w2 < 0:
                ret = (ret * pctabN2[(-1-w2)//2]) % self.n

        return self.quot(ret)

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

            e1val = _from_win(islice(e1bits, Defs.winsize))
            e2val = _from_win(islice(e2bits, Defs.winsize))

            ret *= pctable[e1val + winlen * e2val]
            ret %= self.n

        return self.quot(ret)

    def powgh(self, e1, e2):
        loge = max(int(e1).bit_length(), int(e2).bit_length())
        gcomb = hcomb = None
        for c in self.combs:
            if loge <= c[0].nbits:
                (gcomb, hcomb) = c
                break
        if gcomb is None or hcomb is None:
            raise ValueError("got unexpectedly large exponent in powgh")

        ret = 1
        for (e1vs, e2vs) in zip(gcomb.to_comb_exp(e1), hcomb.to_comb_exp(e2)):
            if ret != 1:
                ret = pow(ret, 2, self.n)

            for (idx, (e1v, e2v)) in enumerate(zip(e1vs, e2vs)):
                if e1v != 0:
                    ret = (ret * gcomb[idx * gcomb.points_per_subcomb + e1v - 1]) % self.n
                if e2v != 0:
                    ret = (ret * hcomb[idx * hcomb.points_per_subcomb + e2v - 1]) % self.n

        return self.quot(ret)

    def mul(self, m1, m2):
        return self.quot((m1 * m2) % self.n)

    def inv2(self, b1, b2):
        b12Inv = lutil.invert_modp(b1 * b2, self.n)
        return (self.quot((b2 * b12Inv) % self.n), self.quot((b1 * b12Inv) % self.n))

    def inv5(self, b1, b2, b3, b4, b5):
        b12 = (b1 * b2) % self.n
        b34 = (b3 * b4) % self.n
        b1234 = (b12 * b34) % self.n
        b12345 = (b1234 * b5) % self.n

        b12345Inv = lutil.invert_modp(b12345, self.n)
        b1234Inv = (b12345Inv * b5) % self.n
        b34Inv = (b1234Inv * b12) % self.n
        b12Inv = (b1234Inv * b34) % self.n

        return ((b12Inv * b2) % self.n, (b12Inv * b1) % self.n, (b34Inv * b4) % self.n, (b34Inv * b3) % self.n, (b12345Inv * b1234) % self.n)

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

    def test_pow2_wnaf():
        "pow2_wnaf,RSA_chal,RSA_rand"

        (b1, b2, e1, e2) = ( lutil.rand.getrandbits(2048) for _ in range(0, 4) )
        (b1Inv, b2Inv)= t1.inv2(b1, b2)
        out1 = t1.quot((pow(b1, e1, t1.n) * pow(b2, e2, t1.n)) % t1.n)
        t1o = t1.pow2_wnaf(b1, b1Inv, e1, b2, b2Inv, e2)

        (b1Inv, b2Inv)= t2.inv2(b1, b2)
        out2 = t2.quot((pow(b1, e1, t2.n) * pow(b2, e2, t2.n)) % t2.n)
        t2o = t2.pow2_wnaf(b1, b1Inv, e1, b2, b2Inv, e2)

        return (out1 == t1o, out2 == t2o)

    def test_pow2():
        "pow2,RSA_chal,RSA_rand"

        (b1, b2, e1, e2) = ( lutil.rand.getrandbits(2048) for _ in range(0, 4) )
        out1 = t1.quot((pow(b1, e1, t1.n) * pow(b2, e2, t1.n)) % t1.n)
        t1o = t1.pow2(b1, e1, b2, e2)

        out2 = t2.quot((pow(b1, e1, t2.n) * pow(b2, e2, t2.n)) % t2.n)
        t2o = t2.pow2(b1, e1, b2, e2)

        return (out1 == t1o, out2 == t2o)

    def test_powgh():
        "powgh,RSA_chal,RSA_rand"

        (e1, e2) = ( lutil.rand.getrandbits(2 * 2048 + Defs.chalbits + 2) for _ in range(0, 2) )

        out1 = t1.quot((pow(2, e1, t1.n) * pow(3, e2, t1.n)) % t1.n)
        t1o = t1.powgh(e1, e2)

        (e1_s, e2_s) = ( x >> (2048 + Defs.chalbits) for x in (e1, e2) )
        out2 = t2.quot((pow(5, e1_s, t2.n) * pow(7, e2_s, t2.n)) % t2.n)
        t2o = t2.powgh(e1_s, e2_s)

        return (out1 == t1o, out2 == t2o)

    def test_inv2():
        "inv2,RSA_chal,RSA_rand"

        (e1, e2) = ( lutil.rand.getrandbits(2048) for _ in range(0, 2) )
        (e1Inv, e2Inv) = t1.inv2(e1, e2)
        t1pass = t1.quot((e1 * e1Inv) % t1.n) == 1 and t1.quot((e2 * e2Inv) % t1.n) == 1

        (e1_s, e2_s) = ( x >> 1536 for x in (e1, e2) )
        (e1_sInv, e2_sInv) = t2.inv2(e1_s, e2_s)
        t2pass = t2.quot((e1_s * e1_sInv) % t2.n) == 1 and t2.quot((e2_s * e2_sInv) % t2.n) == 1

        return (t1pass, t2pass)

    def test_inv5():
        "inv5,RSA_chal,RSA_rand"

        eVals = tuple( lutil.rand.getrandbits(2048) for _ in range(0, 5) )
        eInvs = t1.inv5(*eVals)
        t1pass = all( t1.quot((e * eInv) % t1.n) == 1 for (e, eInv) in zip(eVals, eInvs) )

        eInvs = t2.inv5(*eVals)
        t2pass = all( t2.quot((e * eInv) % t2.n) == 1 for (e, eInv) in zip(eVals, eInvs) )

        return (t1pass, t2pass)

    tu.run_all_tests(nreps, "group_ops", test_pow2_wnaf, test_pow2, test_powgh, test_inv2, test_inv5)

if __name__ == "__main__":
    try:
        nr = int(sys.argv[1])
    except:
        nr = 32
    main(nr)
