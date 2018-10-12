#!/usr/bin/python
#
# (C) 2018 Riad S. Wahby <rsw@cs.stanford.edu>

import sys

import libsig.group_ops as lgops
from libsig.defs import Defs
import libsig.prng as lprng
import libsig.util as lutil

# python 2/3 hack
if sys.version_info[0] == 2:
    range = xrange      # pylint: disable=redefined-builtin,undefined-variable

class HSSigProver(object):
    def __init__(self, p, q, gops=None):
        self.p = p
        self.q = q
        if gops is None:
            gops = lgops.RSAGroupOps(Defs.Grsa.modulus, Defs.Grsa.g, Defs.Grsa.h)
        self.gops = gops

        assert lutil.is_prime(p)
        assert lutil.is_prime(q)

    def sign(self, C1, s, msg):
        # NOTE one assumes that s will have been encrypted to our public key.
        #      This function expects that s has already been decrypted.
        n = self.p * self.q
        assert C1 == self.gops.powgh(n, s), "C1 does not appear to commit to our RSA modulus with opening s"

        ###
        ### Preliminaries: compute values P needs to run the ZKPOK
        ###
        # find t
        t = None
        for t in Defs.primes:
            w = lutil.sqrt_modn(t, self.p, self.q)
            if w is not None:
                break
        if w is None or t is None:
            RuntimeError("did not find a prime quadratic residue less than 1000 mod N!")

        a = (w**2 - t) // n
        assert a * n == w**2 - t, "w^2 - t was not divisible by N!"

        # commitment to w
        s1 = self.gops.rand_scalar()
        (C2, s1) = self.gops.powgh(w, s1)

        # inverses of C1 and C2
        (C1Inv, C2Inv) = self.gops.inv2(C1, C2)

        ###
        ### P's first message: commit to randomness
        ###
        # P's randomness
        (r_w, r_w2, r_s1, r_a, r_an, r_s1w, r_sa) = ( self.gops.rand_scalar() for _ in range(0, 7) )

        # P's first message
        A = self.gops.powgh(r_w, r_s1)
        B = self.gops.mul(self.gops.pow(C2Inv, r_w), self.gops.powgh(r_w2, r_s1w))
        C = self.gops.mul(self.gops.pow(C1Inv, r_a), self.gops.powgh(r_an, r_sa))
        D = r_w2 - r_an

        ###
        ### V's message: random challenge and random prime
        ###
        (chal, ell) = lprng.fs_chal(self.gops.g, self.gops.h, self.gops.n, C1, C2, t, A, B, C, D, msg)

        ###
        ### P's second message: compute quotient message
        ###
        # compute z' = c*(w, w2, s1, a, an, s1w, sa) + (r_w, r_w2, r_s1, r_a, r_an, r_s1w, r_sa)
        z_w = chal * w + r_w
        z_w2 = chal * w * w + r_w2
        z_s1 = chal * s1 + r_s1
        z_a = chal * a + r_a
        z_an = chal * a * n + r_an
        z_s1w = chal * s1 * w + r_s1w
        z_sa = chal * s * a + r_sa

        # compute quotient commitments
        Aq = self.gops.powgh(z_w // ell, z_s1 // ell)
        Bq = self.gops.mul(self.gops.pow(C2Inv, z_w // ell), self.gops.powgh(z_w2 // ell, z_s1w // ell))
        Cq = self.gops.mul(self.gops.pow(C1Inv, z_a // ell), self.gops.powgh(z_an // ell, z_sa // ell))
        Dq = (z_w2 - z_an) // ell

        # compute z'
        z_prime = tuple( z_foo % ell for z_foo in (z_w, z_w2, z_s1, z_a, z_an, z_s1w, z_sa) )

        ###
        ### signature: (chal, ell, Aq, Bq, Cq, Dq, z_prime)
        ###
        sigma = (chal, ell, Aq, Bq, Cq, Dq, z_prime)
        return (C2, t, sigma)
