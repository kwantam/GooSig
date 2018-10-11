#!/usr/bin/python
#
# (C) 2018 Riad S. Wahby <rsw@cs.stanford.edu>

import hashlib

import libsig.util as lutil

# NOTE AES-CTR is almost certainly much faster, but this sticks to the Python stdlib
class SHA512PRNG(object):
    def __init__(self, prng_key):
        self.prng = hashlib.sha512()
        self.rnum = 0
        self.prng.update(b"libsig_prng,%s" % str(prng_key).encode("ascii"))

    def _next_rand(self):
        self.prng.update(b",%d" % self.rnum)
        self.rnum += 1
        return self.prng.hexdigest()

    def getrandbits(self, nbits):
        r = 0
        b = 0
        while b < nbits:
            r <<= 512
            r += int(self._next_rand(), 16)
            b += 512
        if b > nbits:
            r >>= (b - nbits)
        return r

    def _randrange(self, maxval):
        nbits = lutil.clog2(maxval)
        ret = maxval
        while ret >= maxval:
            ret = self.getrandbits(nbits)
        return ret

    def randrange(self, start, stop=None):
        if stop is None:
            return self._randrange(start)

        if stop <= start:
            raise ValueError("require stop > start in randrange(start, stop)")

        return start + self._randrange(stop - start)
