#!/usr/bin/python
#
# (C) 2018 Riad S. Wahby <rsw@cs.stanford.edu>

import sys

from libsig.defs import Defs
import libsig.util as lutil

# python 2/3 hack
if sys.version_info[0] == 2:
    range = xrange      # pylint: disable=redefined-builtin,undefined-variable

class _RSAGroupOps(object):
    pctable = None

    def __init__(self, n, g, h, prng=None):
        self.n = n
        self.g = g
        self.h = h

        # random scalars can cut off a bit without any danger; this makes rand_scalar faster
        self.nbits_rand = lutil.clog2(n) - 1

        # precompute windows of g, h
        self._precompute_windows()

        if prng is None:
            self.prng = lutil.rand
        else:
            self.prng = prng

    def _precompute_windows(self):
        self.pctable = [None] * (2 ** (2 * Defs.winsize))
        winlen = 2 ** Defs.winsize
        for i in range(0, winlen):
            if i == 0:
                self.pctable[0] = 1
            else:
                self.pctable[i] = (self.pctable[i-1] * self.g) % self.n
            for j in range(1, winlen):
                self.pctable[i + winlen * j] = (self.pctable[i + winlen * (j-1)] * self.h) % self.n

    def rand_scalar(self):
        return self.prng.getrandbits(self.nbits_rand)

    def pow(self, b, e):
        return pow(b, e, self.n)

    @staticmethod
    def _from_win(ls):
        ret = 0
        for v in ls[:-1]:
            ret += 1 if v else 0
            ret <<= 1
        ret += 1 if ls[-1] else 0
        return ret

    def powgh(self, e1, e2):
        e1bits = lutil.num_to_bits(e1)
        e2bits = lutil.num_to_bits(e2)
        nwins = (max(len(e1bits), len(e2bits)) + (Defs.winsize - 1)) // Defs.winsize

        # pad out exponents to multiple of winsize (for ease of windowing)
        e1bits = ([False] * (Defs.winsize * nwins - len(e1bits))) + e1bits
        e2bits = ([False] * (Defs.winsize * nwins - len(e2bits))) + e2bits

        winlen = 2 ** Defs.winsize
        ret = 1
        for win in range(0, nwins):
            ret = pow(ret, winlen, self.n)

            e1val = self._from_win(e1bits[Defs.winsize*win:Defs.winsize*(win+1)])
            e2val = self._from_win(e2bits[Defs.winsize*win:Defs.winsize*(win+1)])

            ret *= self.pctable[e1val + winlen * e2val]
            ret %= self.n

        return ret

    def div(self, n, d):
        dInv = lutil.invert_modp(d, self.n)
        return (n * dInv) % self.n

class Commit(object):
    def __init__(self, gops):
        self.gops = gops

    def commit(self, v):
        s = self.gops.rand_scalar()
        com = self.gops.powgh(v, s)
        return (com, s)
