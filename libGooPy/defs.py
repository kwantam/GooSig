#!/usr/bin/python
#
# (C) 2018 Dan Boneh, Riad S. Wahby <rsw@cs.stanford.edu>

from itertools import takewhile
import hashlib

import libGooPy.primes as lprimes

class Defs(object):
    max_rsa_keysize = 4096
    rand_exponent_size = 2048

    winsize = 6
    max_rsa_comb_size = 512
    max_bqf_comb_size = 64

    hashfn = hashlib.sha256
    chalbits = 128
    pdelta = 18
    ft_prime_opts = None    # gets set below

    # this is the list of primes from which P can choose to prove ability to take sqrt mod her RSA modulus
    # each one is a QR mod N with probability 1/2, so this list suffices except with vanishing probability
    primes = list( takewhile(lambda x: x < 1000, lprimes.primes()) )

Defs.ft_prime_opts = lprimes.gen_ft_prime_opts(Defs.chalbits, Defs.pdelta)
