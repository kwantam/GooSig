#!/usr/bin/python
#
# (C) 2018 Dan Boneh, Riad S. Wahby <rsw@cs.stanford.edu>

import sys

import libGooPy.group_ops as lgops
from libGooPy.consts import Grsa2048
from libGooPy.defs import Defs
import libGooPy.prng as lprng
import libGooPy.util as lutil

# python 2/3 hack
if sys.version_info[0] == 2:
    range = xrange      # pylint: disable=redefined-builtin,undefined-variable

class GooSigSigner(object):
    def __init__(self, rsakey, gops=None):
        self.rsakey = rsakey
        if gops is None:
            modbits = lutil.clog2(self.rsakey.n)
            gops = lgops.RSAGroupOps(Grsa2048, modbits)
        self.gops = gops

    def sign(self, C0, C1, msg):
        s_prime = self.rsakey.decrypt(C0)
        s = lprng.expand_sprime(s_prime)
        assert C1 == self.gops.reduce(self.gops.powgh(self.rsakey.n, s)), "C1 does not appear to commit to our RSA modulus with opening s"

        ###
        ### Preliminaries: compute values P needs to run the ZKPOK
        ###
        # find t
        t = None
        for t in Defs.primes:
            w = lutil.sqrt_modn(t, self.rsakey.p, self.rsakey.q)
            if w is not None:
                break
        if w is None or t is None:
            RuntimeError("did not find a prime quadratic residue less than 1000 mod N!")

        a = (w**2 - t) // self.rsakey.n
        assert a * self.rsakey.n == w**2 - t, "w^2 - t was not divisible by N!"

        # commitment to w
        s1 = self.gops.rand_scalar()
        C2 = self.gops.reduce(self.gops.powgh(w, s1))

        # inverses of C1 and C2
        (C1Inv, C2Inv) = self.gops.inv2(C1, C2)

        ###
        ### P's first message: commit to randomness
        ###
        # P's randomness
        (r_w, r_w2, r_s1, r_a, r_an, r_s1w, r_sa) = ( self.gops.rand_scalar() for _ in range(0, 7) )

        # P's first message
        A = self.gops.reduce(self.gops.powgh(r_w, r_s1))
        B = self.gops.reduce(self.gops.mul(self.gops.pow(C2Inv, C2, r_w), self.gops.powgh(r_w2, r_s1w)))
        C = self.gops.reduce(self.gops.mul(self.gops.pow(C1Inv, C1, r_a), self.gops.powgh(r_an, r_sa)))
        D = r_w2 - r_an

        ###
        ### V's message: random challenge and random prime
        ###
        (chal, ell) = lprng.fs_chal(False, self.gops.desc, C1, C2, t, A, B, C, D, msg)

        ###
        ### P's second message: compute quotient message
        ###
        # compute z' = c*(w, w2, s1, a, an, s1w, sa) + (r_w, r_w2, r_s1, r_a, r_an, r_s1w, r_sa)
        z_w = chal * w + r_w
        z_w2 = chal * w * w + r_w2
        z_s1 = chal * s1 + r_s1
        z_a = chal * a + r_a
        z_an = chal * a * self.rsakey.n + r_an
        z_s1w = chal * s1 * w + r_s1w
        z_sa = chal * s * a + r_sa

        # compute quotient commitments
        Aq = self.gops.reduce(self.gops.powgh(z_w // ell, z_s1 // ell))
        Bq = self.gops.reduce(self.gops.mul(self.gops.pow(C2Inv, C2, z_w // ell), self.gops.powgh(z_w2 // ell, z_s1w // ell)))
        Cq = self.gops.reduce(self.gops.mul(self.gops.pow(C1Inv, C2, z_a // ell), self.gops.powgh(z_an // ell, z_sa // ell)))
        Dq = (z_w2 - z_an) // ell

        # compute z'
        z_prime = tuple( z_foo % ell for z_foo in (z_w, z_w2, z_s1, z_a, z_an, z_s1w, z_sa) )

        ###
        ### signature: (chal, ell, Aq, Bq, Cq, Dq, z_prime)
        ###
        sigma = (chal, ell, Aq, Bq, Cq, Dq, z_prime)
        return (C2, t, sigma)
