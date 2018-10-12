#!/usr/bin/python
#
# (C) 2018 Riad S. Wahby <rsw@cs.stanford.edu>

import itertools
import sys

from libsig.defs import Defs
import libsig.util as lutil

# python 2/3 hack
if sys.version_info[0] == 2:
    zip = itertools.izip    # pylint: disable=redefined-builtin,no-member
    range = xrange          # pylint: disable=redefined-builtin,undefined-variable

class RSAGroupOps(object):
    def __init__(self, n, g, h, prng=None):
        self.n = n
        self.g = g
        self.h = h

        # random scalars can cut off a bit without any danger; this makes rand_scalar faster
        self.nbits_rand = lutil.clog2(n) - 1

        # largest exponent we will to handle is something like 128 + 2 * log2(n) + 1
        combsize = Defs.combsize
        self.big_comb = combsize * ((2 * lutil.clog2(self.n) + Defs.chalbits + 1 + combsize) // combsize)
        self.small_comb = combsize * ((lutil.clog2(self.n) + combsize - 1) // combsize)
        self.combsize = combsize

        # precompute combs for g, h
        # NOTE: you really want to store these on disk!
        self.gcomb_big = self._precomp_comb(self.g, self.big_comb)
        self.gcomb_small = self._precomp_comb(self.g, self.small_comb)
        self.hcomb_big = self._precomp_comb(self.h, self.big_comb)
        self.hcomb_small = self._precomp_comb(self.h, self.small_comb)

        if prng is None:
            # /dev/urandom
            self.prng = lutil.rand
        else:
            # provided RNG
            self.prng = prng

    def pow(self, b, e):
        # native pow is pretty fast already
        return pow(b, e, self.n)

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
        e1bits = lutil.num_to_bits(e1)
        e2bits = lutil.num_to_bits(e2)
        nwins = (max(len(e1bits), len(e2bits)) + (Defs.winsize - 1)) // Defs.winsize
        e1bits = ([False] * (Defs.winsize * nwins - len(e1bits))) + e1bits
        e2bits = ([False] * (Defs.winsize * nwins - len(e2bits))) + e2bits

        winlen = 2 ** Defs.winsize
        ret = 1
        for win in range(0, nwins):
            if win > 0:
                ret = pow(ret, winlen, self.n)

            e1val = self._from_win(e1bits[Defs.winsize * win : Defs.winsize * (win + 1)])
            e2val = self._from_win(e2bits[Defs.winsize * win : Defs.winsize * (win + 1)])

            ret *= pctable[e1val + winlen * e2val]
            ret %= self.n

        return ret

    @staticmethod
    def _from_win(ls):
        ret = 0
        for v in ls[:-1]:
            ret += 1 if v else 0
            ret <<= 1
        ret += 1 if ls[-1] else 0
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
        ebits = lutil.num_to_bits(e)
        ebits = ([False] * (max_comb - int(e).bit_length())) + lutil.num_to_bits(e)
        assert len(ebits) == max_comb
        bpc = max_comb // self.combsize
        comb_words_iter = zip(*( ebits[j*bpc:(j+1)*bpc] for j in range(0, self.combsize) ))
        return ( self._from_win(x) for x in comb_words_iter )

    def powgh(self, e1, e2):
        loge = max(int(e1).bit_length(), int(e2).bit_length())
        if loge <= self.small_comb:
            max_comb = self.small_comb
            gcomb = self.gcomb_small
            hcomb = self.hcomb_small
        elif loge <= self.big_comb:
            max_comb = self.big_comb
            gcomb = self.gcomb_big
            hcomb = self.hcomb_big
        else:
            raise ValueError("got unexpectedly large exponent in powgh")

        e1bits = self._to_comb_bits(e1, max_comb)
        e2bits = self._to_comb_bits(e2, max_comb)

        ret = 1
        for (e1v, e2v) in zip(e1bits, e2bits):
            ret = pow(ret, 2, self.n)
            ret = (ret * gcomb[e1v]) % self.n
            ret = (ret * hcomb[e2v]) % self.n

        return ret

    def mul(self, m1, m2):
        return (m1 * m2) % self.n

    def inv2(self, b1, b2):
        b12Inv = lutil.invert_modp(b1 * b2, self.n)
        return ((b2 * b12Inv) % self.n, (b1 * b12Inv) % self.n)

    def rand_scalar(self):
        return self.prng.getrandbits(self.nbits_rand)
