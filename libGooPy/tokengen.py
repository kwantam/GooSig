#!/usr/bin/python
#
# (C) 2018 Dan Boneh, Riad S. Wahby <rsw@cs.stanford.edu>

from libGooPy.consts import Grsa2048
from libGooPy.defs import Defs
import libGooPy.group_ops as lgops
import libGooPy.util as lu

class GooSigTokenGen(object):
    def __init__(self, gops):
        self.gops = gops if gops is not None else lgops.RSAGroupOps(Grsa2048, Defs.max_rsa_keysize)

    def send_tokens(self, rsakey):
        # NOTE: in the real protocol, select a 256-bit s' and expand to 2048-bit s,
        #       e.g., with AESEnc(s', 0), ..., AESEnc(s', 3)
        s = lu.rand.randrange(rsakey.n)

        # the challenge: a commitment to the RSA modulus
        C1 = self.gops.reduce(self.gops.powgh(rsakey.n, s))

        # s, encrypted to the pubkey
        C0_pre = rsakey.encrypt(s)

        # make a ciphertext C0 indistinguishable from random (max_rsa_keysize + 8)-bit integer
        ct_lim = 1 << (Defs.max_rsa_keysize + 8)
        r_lim = 1 + ct_lim // rsakey.n
        C0 = ct_lim
        while C0 >= ct_lim:
            c0_r = lu.rand.randrange(r_lim)
            C0 = C0_pre + c0_r * rsakey.n

        return (C0, C1)
