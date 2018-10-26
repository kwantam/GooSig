#!/usr/bin/python
#
# (C) 2018 Dan Boneh, Riad S. Wahby <rsw@cs.stanford.edu>

from libGooPy.defs import Defs
import libGooPy.primes as lprimes
import libGooPy.util as lutil

# NOTE AES-CTR is almost certainly much faster, but this sticks to the Python stdlib
class HashPRNG(object):
    def __init__(self, prng_key):
        self.prng = Defs.hashfn(b"libGooPy_prng,%s" % str(prng_key).encode("utf-8"))
        self.rnum = 0
        self.r = 0

    def _next_rand(self):
        self.prng.update(b"%016x" % self.rnum)
        self.rnum += 1
        return int(self.prng.hexdigest(), 16)

    def getrandbits(self, nbits):
        r = self.r
        b = r.bit_length()
        hashbits = self.prng.digest_size * 8
        while b < nbits:
            r <<= hashbits
            r += self._next_rand()
            b += hashbits
        self.r = r & ((1 << (b - nbits)) - 1)
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

def fs_chal(ver_only, *args):
    fs_hash = Defs.hashfn(b"libGooPy:")

    # hash the inputs
    for arg in args:
        fs_hash.update(str(arg).encode("utf-8"))

    # build PRNG instance and compute chal and ell
    prng_key = fs_hash.hexdigest()
    fs_hash_prng = HashPRNG(prng_key)
    chal = fs_hash_prng.getrandbits(Defs.chalbits)
    if ver_only:
        ell = fs_hash_prng.getrandbits(Defs.chalbits)
    else:
        ell = lprimes.next_prime(fs_hash_prng.getrandbits(Defs.chalbits))

    return (chal, ell)
