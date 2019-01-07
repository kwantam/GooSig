#!/usr/bin/python
#
# (C) 2018 Dan Boneh, Riad S. Wahby <rsw@cs.stanford.edu>

import sys

import libGooPy.group_ops as lgops
from libGooPy.consts import Grsa2048
from libGooPy.defs import Defs
import libGooPy.primes as lprimes
import libGooPy.prng as lprng

# python 2/3 hack
if sys.version_info[0] == 2:
    range = xrange      # pylint: disable=redefined-builtin,undefined-variable

class GooSigVerifier(object):
    def __init__(self, gops=None):
        if gops is None:
            gops = lgops.RSAGroupOps(Grsa2048, modbits=None)
        self.gops = gops

    def verify(self, pubkey, msg, sigma):
        ## UPDATE 2019 Jan 06: pubkey adds element C3; sigma adds element Eq; z' adds element z'_s2
        (C1, C2, C3, t) = pubkey
        (chal, ell, Aq, Bq, Cq, Dq, Eq, z_prime) = sigma
        (zp_w, zp_w2, zp_s1, zp_a, zp_an, zp_s1w, zp_sa, zp_s2) = z_prime

        # make sure that the public key is valid
        if t not in Defs.primes:
            # t must be one of the small primes in our list
            return False
        ## UPDATE 2019 Jan 06: check that C3 and Dq are reduced, too
        if not all( self.gops.is_reduced(b) for b in (C1, C2, C3, Aq, Bq, Cq, Dq) ):
            # all group elements must be the "canonical" element of the quotient group (Z/n)/{1,-1}
            return False

        ## UPDATE 2019 Jan 06: invert C3 and Dq, too
        # compute inverses of C1, C2, C3, Aq, Bq, Cq, Dq
        # NOTE: Since we're inverting C1 and C2, we can get inverses of Aq, Bq, Cq for ~free.
        #       This lets us use signed-digit exponentiation below, which is much faster.
        (C1Inv, C2Inv, C3Inv, AqInv, BqInv, CqInv, DqInv) = self.gops.inv7(C1, C2, C3, Aq, Bq, Cq, Dq)

        ###
        ### Step 1: reconstruct A, B, C, and D from signature
        ###
        ## UPDATE 2019 Jan 06: (B, C) renamed to (C, D); compute new group element B
        A = self.gops.reduce(self.gops.mul(self.gops.pow2(Aq, AqInv, ell, C2Inv, C2, chal), self.gops.powgh(zp_w, zp_s1)))
        B = self.gops.reduce(self.gops.mul(self.gops.pow2(Bq, BqInv, ell, C3Inv, C3, chal), self.gops.powgh(zp_a, zp_s2)))
        C = self.gops.reduce(self.gops.mul(self.gops.pow2(Cq, CqInv, ell, C2Inv, C2, zp_w), self.gops.powgh(zp_w2, zp_s1w)))
        D = self.gops.reduce(self.gops.mul(self.gops.pow2(Dq, DqInv, ell, C1Inv, C1, zp_a), self.gops.powgh(zp_an, zp_sa)))

        ## UPDATE 2019 Jan 06: D renamed to E
        # make sure sign of (zp_w2 - zp_an) is positive
        zp_w2_m_an = zp_w2 - zp_an
        E = Eq * ell + zp_w2_m_an - t * chal
        if zp_w2_m_an < 0:
            E += ell

        ###
        ### Step 2: recompute implicitly claimed V message, viz., chal and ell
        ###
        ## UPDATE 2019 Jan 06: C3, E added as inputs to the hash function
        (chal_out, ell_r_out) = lprng.fs_chal(True, self.gops.desc, C1, C2, C3, t, A, B, C, D, E, msg)

        # final checks
        # chal has to match AND 0 <= (ell_r_out - ell) <= elldiff_max AND ell is prime
        elldiff = ell - ell_r_out
        if chal != chal_out or elldiff < 0 or elldiff > Defs.elldiff_max or not lprimes.is_prime(ell):
            return False

        return True
