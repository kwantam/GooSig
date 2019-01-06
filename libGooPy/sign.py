#!/usr/bin/python
#
# (C) 2018 Dan Boneh, Riad S. Wahby <rsw@cs.stanford.edu>

import hashlib
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
        C0dec = self.rsakey.decrypt(C0)
        s_prime = C0dec & ((1 << 256) - 1)
        hC1 = C0dec >> 256
        s = lprng.expand_sprime(s_prime)
        assert C1 == self.gops.reduce(self.gops.powgh(self.rsakey.n, s)), "C1 does not appear to commit to our RSA modulus with opening s"
        assert hC1 == int(hashlib.sha256(str(C1).encode("utf-8")).hexdigest(), 16)

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

        ## UPDATE 2019 Jan 06: commit to value of a with randomness s2
        # commitment to a
        s2 = self.gops.rand_scalar()
        C3 = self.gops.reduce(self.gops.powgh(a, s2))

        # inverses of C1 and C2
        (C1Inv, C2Inv) = self.gops.inv2(C1, C2)

        ###
        ### P's first message: commit to randomness
        ###
        ## UPDATE 2019 Jan 06: added r_s2 to P's randomness
        # P's randomness (except for r_s1; see "V's message", below)
        (r_w, r_w2, r_a, r_an, r_s1w, r_sa, r_s2) = ( self.gops.rand_scalar() for _ in range(0, 7) )

        ## UPDATE 2019 Jan 06: (B, C, D) renamed to (C, D, E); new group element B added
        # P's first message (except for A; see "V's message", below)
        B = self.gops.reduce(self.gops.powgh(r_a, r_s2))
        C = self.gops.reduce(self.gops.mul(self.gops.pow(C2Inv, C2, r_w), self.gops.powgh(r_w2, r_s1w)))
        D = self.gops.reduce(self.gops.mul(self.gops.pow(C1Inv, C1, r_a), self.gops.powgh(r_an, r_sa)))
        E = r_w2 - r_an

        ###
        ### V's message: random challenge and random prime
        ###
        ell = None
        while ell is None or ell.bit_length() != 128:
            # randomize the signature until Fiat-Shamir returns an admissable ell
            # NOTE it's not necessary to re-start the whole signature!
            #      Just pick a new r_s1, which only requires re-computing A.
            r_s1 = self.gops.rand_scalar()
            A = self.gops.reduce(self.gops.powgh(r_w, r_s1))
            ## UPDATE 2019 Jan 06: C3, E added as inputs to the hash function
            (chal, ell) = lprng.fs_chal(False, self.gops.desc, C1, C2, C3, t, A, B, C, D, E, msg)

        ###
        ### P's second message: compute quotient message
        ###
        ## UPDATE 2019 Jan 06: z_s2 added to z
        # compute z = c*(w, w2, s1, a, an, s1w, sa, s2) + (r_w, r_w2, r_s1, r_a, r_an, r_s1w, r_sa, r_s2)
        z_w = chal * w + r_w
        z_w2 = chal * w * w + r_w2
        z_s1 = chal * s1 + r_s1
        z_a = chal * a + r_a
        z_an = chal * a * self.rsakey.n + r_an
        z_s1w = chal * s1 * w + r_s1w
        z_sa = chal * s * a + r_sa
        z_s2 = chal * s2 + r_s2

        ## UPDATE 2019 Jan 06: (Bq, Cq, Dq) renamed to (Cq, Dq, Eq); new group element Bq added
        # compute quotient commitments
        Aq = self.gops.reduce(self.gops.powgh(z_w // ell, z_s1 // ell))
        Bq = self.gops.reduce(self.gops.powgh(z_a // ell, z_s2 // ell))
        Cq = self.gops.reduce(self.gops.mul(self.gops.pow(C2Inv, C2, z_w // ell), self.gops.powgh(z_w2 // ell, z_s1w // ell)))
        Dq = self.gops.reduce(self.gops.mul(self.gops.pow(C1Inv, C2, z_a // ell), self.gops.powgh(z_an // ell, z_sa // ell)))
        Eq = (z_w2 - z_an) // ell

        ## UPDATE 2019 Jan 06: z'_s2 added to z'
        # compute z'
        z_prime = tuple( z_foo % ell for z_foo in (z_w, z_w2, z_s1, z_a, z_an, z_s1w, z_sa, z_s2) )

        ###
        ### signature: (chal, ell, Aq, Bq, Cq, Dq, z_prime)
        ###
        ## UPDATE 2019 Jan 06: C3, Eq added to the outputs from the sign function
        sigma = (chal, ell, Aq, Bq, Cq, Dq, Eq, z_prime)
        return (C2, C3, t, sigma)
