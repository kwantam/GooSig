#!/usr/bin/python
#
# (C) 2018 Dan Boneh, Riad S. Wahby <rsw@cs.stanford.edu>

from libGooPy.defs import Defs
import libGooPy.primes as lprimes
import libGooPy.util as lutil

# NOTE AES-CTR is almost certainly faster for most machines, but this approach sticks to the Python stdlib
class _HashPRNG(object):
    def __init__(self, prng):
        self.prng = prng
        self.rnum = 0
        self.r_save = 0

    def _next_rand(self):
        self.prng.update(b"%016x" % self.rnum)
        self.rnum += 1
        return int(self.prng.hexdigest(), 16)

    def getrandbits(self, nbits):
        r = self.r_save
        b = r.bit_length()
        hashbits = self.prng.digest_size * 8
        while b < nbits:
            r <<= hashbits
            r += self._next_rand()
            b += hashbits
        self.r_save = r & ((1 << (b - nbits)) - 1)
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

    @classmethod
    def new(cls, *seed):
        new_prng = Defs.hashfn(b"libGooPy:")
        for s in seed:
            new_prng.update(str(s).encode("utf-8"))
        return cls(new_prng)

def fs_chal(ver_only, *args):
    fs_hash_prng = _HashPRNG.new(*args)

    chal = fs_hash_prng.getrandbits(Defs.chalbits)

    if ver_only:
        ell = fs_hash_prng.getrandbits(Defs.chalbits)
    if not ver_only:
        # for prover, call next_prime on ell_r to get ell
        ell = lprimes.next_prime(fs_hash_prng.getrandbits(Defs.chalbits), Defs.elldiff_max)

    return (chal, ell)

def expand_sprime(s_prime):
    return _HashPRNG.new(s_prime).getrandbits(Defs.rand_exponent_size)
