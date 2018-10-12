#!/usr/bin/python
#
# (C) 2018 Riad S. Wahby <rsw@cs.stanford.edu>

import sys

import libsig.group_ops as lgops
from libsig.defs import Defs
import libsig.prng as lprng

# python 2/3 hack
if sys.version_info[0] == 2:
    range = xrange      # pylint: disable=redefined-builtin,undefined-variable

class GooSigVerifier(object):
    def __init__(self, gops=None):
        if gops is None:
            gops = lgops.RSAGroupOps(Defs.Grsa, modbits=None)
        self.gops = gops

    def verify(self, pubkey, msg, sigma):
        (C1, C2, t) = pubkey
        (chal, ell, Aq, Bq, Cq, Dq, z_prime) = sigma
        (zp_w, zp_w2, zp_s1, zp_a, zp_an, zp_s1w, zp_sa) = z_prime

        # make sure that the public key is valid!
        if t not in Defs.primes:
            raise ValueError("public key element t is not a prime less than 1000")

        # compute inverses of C1 and C2
        (C1Inv, C2Inv) = self.gops.inv2(C1, C2)

        ###
        ### Step 1: reconstruct A, B, C, and D from signature
        ###
        A = self.gops.mul(self.gops.pow2(Aq, ell, C2Inv, chal), self.gops.powgh(zp_w, zp_s1))
        B = self.gops.mul(self.gops.pow2(Bq, ell, C2Inv, zp_w), self.gops.powgh(zp_w2, zp_s1w))
        C = self.gops.mul(self.gops.pow2(Cq, ell, C1Inv, zp_a), self.gops.powgh(zp_an, zp_sa))

        # make sure sign of (zp_w2 - zp_an) is positive
        zp_w2_m_an = zp_w2 - zp_an
        if zp_w2_m_an < 0:
            zp_w2_m_an += ell
        D = Dq * ell + zp_w2_m_an - t * chal

        ###
        ### Step 2: recompute implicitly claimed V message, viz., chal and ell
        ###
        (chal_out, ell_out) = lprng.fs_chal(self.gops.g, self.gops.h, self.gops.n, C1, C2, t, A, B, C, D, msg)

        # final check
        if chal != chal_out or ell != ell_out:
            return False
        return True
