#!/usr/bin/python
#
# (C) 2018 Dan Boneh, Riad S. Wahby <rsw@cs.stanford.edu>

import libGooPy.primes as lprimes
import libGooPy.util as lutil

class RSAKey(object):
    def __init__(self, p, q):
        self.p = p
        self.q = q
        assert lprimes.is_prime(p)
        assert lprimes.is_prime(q)
        self.n = self.p * self.q

        # find a decryption exponent
        self.lam = (p - 1) * (q - 1) // lutil.gcd(p - 1, q - 1)
        for e in lprimes.primes_skip(1):
            if e > 1000:
                raise RuntimeError("could find a suitable exponent!")
            d = lutil.invert_modp(e, self.lam)
            if d is not None:
                self.e = e
                self.d = d
                break

        assert (self.d * self.e) % self.lam == 1

    def encrypt(self, m):
        # NOTE this is not real RSA encryption! You should use RSA-OAEP or the like.
        return pow(m, self.e, self.n)

    def decrypt(self, c):
        # NOTE this is not real RSA decryption! You should use RSA-OAEP or the like.
        return pow(c, self.d, self.n)
