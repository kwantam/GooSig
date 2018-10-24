#!/usr/bin/python
#
# (C) 2018 Riad S. Wahby <rsw@cs.stanford.edu>

import sys

import libGooPy.group_ops as lgops
from libGooPy.consts import Grsa2048
from libGooPy.defs import Defs
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
        (C1, C2, t) = pubkey
        (chal, ell, Aq, Bq, Cq, Dq, z_prime) = sigma
        (zp_w, zp_w2, zp_s1, zp_a, zp_an, zp_s1w, zp_sa) = z_prime

        # make sure that the public key is valid
        if t not in Defs.primes:
            # t must be one of the small primes in our list
            return False
        if not all( self.gops.is_reduced(b) for b in (C1, C2, Aq, Bq, Cq) ):
            # all group elements must be the "canonical" element of the quotient group (Z/n)/{1,-1}
            return False

        # compute inverses of C1, C2, Aq, Bq, Cq
        # NOTE: Since we're inverting C1 and C2, we can get inverses of Aq, Bq, Cq for ~free.
        #       This lets us use signed-digit exponentiation below, which is much faster.
        (C1Inv, C2Inv, AqInv, BqInv, CqInv) = self.gops.inv5(C1, C2, Aq, Bq, Cq)

        ###
        ### Step 1: reconstruct A, B, C, and D from signature
        ###
        A = self.gops.reduce(self.gops.mul(self.gops.pow2(Aq, AqInv, ell, C2Inv, C2, chal), self.gops.powgh(zp_w, zp_s1)))
        B = self.gops.reduce(self.gops.mul(self.gops.pow2(Bq, BqInv, ell, C2Inv, C2, zp_w), self.gops.powgh(zp_w2, zp_s1w)))
        C = self.gops.reduce(self.gops.mul(self.gops.pow2(Cq, CqInv, ell, C1Inv, C1, zp_a), self.gops.powgh(zp_an, zp_sa)))

        # make sure sign of (zp_w2 - zp_an) is positive
        zp_w2_m_an = zp_w2 - zp_an
        D = Dq * ell + zp_w2_m_an - t * chal
        if zp_w2_m_an < 0:
            D += ell

        ###
        ### Step 2: recompute implicitly claimed V message, viz., chal and ell
        ###
        (chal_out, ell_out) = lprng.fs_chal(self.gops.desc, C1, C2, t, A, B, C, D, msg)

        # final check
        if chal != chal_out or ell != ell_out:
            return False

        return True

    def verify_simple(self, pubkey, msg, sigma):
        (C1, C2) = pubkey
        (chal, ell, Aq, z_prime) = sigma
        (zp_n, zp_s) = z_prime

        # make sure that the public key and signature include valid group elements
        if not all( self.gops.is_reduced(b) for b in (C1, Aq) ):
            # all group elements must be the "canonical" element of the quotient group (Z/n)/{1,-1}
            return False

        # compute inverses of C1 and Aq
        # NOTE: As above, can get inverse of Aq for free from inverse of C1 and then
        #       use signed-digit exponentiation.
        (C1Inv, AqInv) = self.gops.inv2(C1, Aq)

        ###
        ### Step 1: reconstruct A from signature
        ###
        A = self.gops.reduce(self.gops.mul(self.gops.pow2(Aq, AqInv, ell, C1Inv, C1, chal), self.gops.powgh(zp_n, zp_s)))

        ###
        ### Step 2: recompute implicitly claimed V message, viz., chal and ell
        ###
        (chal_out, ell_out) = lprng.fs_chal(self.gops.desc, C1, C2, A, msg)

        # final check
        if chal != chal_out or ell != ell_out:
            return False

        return True
