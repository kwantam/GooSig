#!/usr/bin/python
#
# (C) 2018 Dan Boneh, Riad S. Wahby <rsw@cs.stanford.edu>

import hashlib
import sys

import libGooPy.primes as lprimes
import libGooPy.prng as lprng
import libGooPy.util as lutil

# RSA-OAEP enc/dec using SHA-256
# NOTE this is a non-standard implementation that you should probably not use except for benchmarking
class RSAPubKey(object):
    HASHLEN = 32    # SHA256 digest size

    def __init__(self, n, e):
        self.n = n
        self.e = e

        n_octets = (self.n.bit_length() + 7) // 8
        if n_octets < 128:
            raise ValueError("RSAKey does not support <1024-bit moduli")
        self.max_mlen = n_octets - 2 * self.HASHLEN - 2
        self.dblen = n_octets - 1 - self.HASHLEN
        self.dbmask = (1 << (8 * self.dblen)) - 1
        self.lhash = int(hashlib.sha256(b"libGooPy_RSA_OAEP_LABEL").hexdigest(), 16) << (8 * (self.dblen - self.HASHLEN))

    @staticmethod
    def mask_gen(seed, length):
        return lprng.HashPRNG.new(seed).getrandbits(8 * length)

    def encrypt(self, m):
        mlen = (m.bit_length() + 7) // 8    # round up to some number of bytes
        if mlen > self.max_mlen:
            raise ValueError("message is too long")

        data = self.lhash | (1 << (8 * mlen)) | m
        seed = lutil.rand.getrandbits(8 * self.HASHLEN)

        dbMask = self.mask_gen(seed, self.dblen)
        maskedDb = dbMask ^ data

        sMask = self.mask_gen(maskedDb, self.HASHLEN)
        maskedSeed = sMask ^ seed

        enc_msg = (maskedSeed << (8 * self.dblen)) | maskedDb

        return pow(enc_msg, self.e, self.n)

class RSAKey(RSAPubKey):
    def __init__(self, p, q):
        self.p = p
        self.q = q
        assert p != q
        assert lprimes.is_prime(p)
        assert lprimes.is_prime(q)

        # find a decryption exponent (must be coprime to carmichael(p * q) == lcm(p-1, q-1)
        lam = (p - 1) * (q - 1) // lutil.gcd(p - 1, q - 1)
        for e in lprimes.primes_skip(1):
            if e > 1000:
                raise RuntimeError("could find a suitable exponent!")
            d = lutil.invert_modp(e, lam)
            if d is not None:
                break
        self.d = d

        # now that we have n and e, initialize the pubkey
        RSAPubKey.__init__(self, p * q, e)  # pylint: disable=undefined-loop-variable

        assert (self.d * self.e) % lam == 1

    def get_public_key(self):
        return RSAPubKey(self.n, self.e)

    def decrypt(self, c):
        enc_msg = pow(c, self.d, self.n)

        maskedDb = enc_msg & self.dbmask
        maskedSeed = enc_msg >> (8 * self.dblen)
        if maskedSeed.bit_length() > 8 * self.HASHLEN:
            raise ValueError("invalid ciphertext")

        sMask = self.mask_gen(maskedDb, self.HASHLEN)
        seed = maskedSeed ^ sMask

        dbMask = self.mask_gen(seed, self.dblen)
        data = dbMask ^ maskedDb

        data ^= self.lhash
        if data >> (8 * (self.dblen - self.HASHLEN)) != 0:
            raise ValueError("invalid ciphertext")

        dlen = (7 + data.bit_length()) // 8
        data ^= 1 << (8 * (dlen - 1))
        if data >> (8 * (dlen - 1)) != 0:
            raise ValueError("invalid padding")

        if data.bit_length() > 8 * self.max_mlen:
            raise ValueError("invalid message")

        return data

def main(nreps):
    import libGooPy.test_util as tu     # pylint: disable=bad-option-value,import-outside-toplevel

    def test_endec():
        "RSA endec,RSA2048,RSA4096"

        (p1, q1) = lutil.rand.sample(tu.primes_1024, 2)
        (p2, q2) = lutil.rand.sample(tu.primes_2048, 2)

        r1 = RSAKey(p1, q1)
        rp1 = r1.get_public_key()
        m1 = lutil.rand.getrandbits(512)
        c1 = rp1.encrypt(m1)
        d1 = r1.decrypt(c1)

        r2 = RSAKey(p2, q2)
        rp2 = r2.get_public_key()
        m2 = lutil.rand.getrandbits(512)
        c2 = rp2.encrypt(m2)
        d2 = r2.decrypt(c2)

        return (m1 == d1, m2 == d2)

    tu.run_all_tests(nreps, "RSA", test_endec)

if __name__ == "__main__":
    try:
        nr = int(sys.argv[1])
    except:
        nr = 32
    main(nr)
