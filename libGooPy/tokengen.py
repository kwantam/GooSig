#!/usr/bin/python
#
# (C) 2018 Dan Boneh, Riad S. Wahby <rsw@cs.stanford.edu>

import hashlib

from libGooPy.consts import Grsa2048
from libGooPy.defs import Defs
import libGooPy.group_ops as lgops
import libGooPy.prng as lprng
import libGooPy.util as lu

class GooSigTokGen(object):
    def __init__(self, gops):
        self.gops = gops if gops is not None else lgops.RSAGroupOps(Grsa2048, Defs.max_rsa_keysize)

    def send_tokens(self, rsapubkey):
        # NOTE: in the real protocol, select a 256-bit s' and expand to 2048-bit s,
        #       e.g., with AESEnc(s', 0), ..., AESEnc(s', 3)
        s_prime = lu.rand.getrandbits(256)
        s = lprng.expand_sprime(s_prime)

        # the challenge: a commitment to the RSA modulus
        C1 = self.gops.reduce(self.gops.powgh(rsapubkey.n, s))
        hC1 = int(hashlib.sha256(str(C1).encode("utf-8")).hexdigest(), 16)

        # (Hash(C1) || s_prime), encrypted to the pubkey
        C0_pre = rsapubkey.encrypt((hC1 << 256) | s_prime)

        # make a ciphertext C0 indistinguishable from a random (max_rsa_keysize + 8)-bit integer
        ct_lim = 1 << (Defs.max_rsa_keysize + 8)
        # ceiling of (ct_lim - C0_pre) / rsapubkey.n, ensuring C0_pre + rsapubkey.n * r_lim >= ct_lim
        r_lim = (ct_lim - C0_pre + rsapubkey.n - 1) // rsapubkey.n
        C0 = ct_lim
        while C0 >= ct_lim:
            c0_rand = lu.rand.randint(0, r_lim)     # NOTE randint returns in [0, r_lim]
            C0 = C0_pre + c0_rand * rsapubkey.n

        assert C0 % rsapubkey.n == C0_pre
        return (C0, C1)
