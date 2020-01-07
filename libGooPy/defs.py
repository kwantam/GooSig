#!/usr/bin/python
#
# (C) 2018 Dan Boneh, Riad S. Wahby <rsw@cs.stanford.edu>

from itertools import takewhile
import hashlib

import libGooPy.primes as lprimes

_lambda = 128       # security parameter; used to compute chalbits and ellbits below.
class Defs(object):
    max_rsa_keysize = 4096
    rand_exponent_size = 2048

    winsize = 6
    max_rsa_comb_size = 512
    max_bqf_comb_size = 64

    hashfn = hashlib.sha256
    ## UPDATE 2020 Jan 07: need ell to be in the first 2^(2 * Defs.chal_bits) primes
    ## because of a sqrt attack on Fiat-Shamir-ized Wesolowski proofs. See S3.3 of
    ##     Boneh, Buenz, Fisch. "A Survey of Two Verifiable Delay Functions."
    ##     ePrint # 2018/712, https://eprint.iacr.org/2018/712
    chalbits = _lambda
    ellbits = 1 + 2 * _lambda + int(_lambda - 1).bit_length()
    elldiff_max = 8 * _lambda

    # this is the list of primes from which P can choose to prove ability to take sqrt mod her RSA modulus
    # each one is a QR mod N with probability 1/2, so this list suffices except with vanishing probability
    primes = list( takewhile(lambda x: x < 1000, lprimes.primes()) )
