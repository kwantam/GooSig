#!/usr/bin/python
#
# (C) 2018 Riad S. Wahby <rsw@cs.stanford.edu>


#from libsig.defs import Defs
import libsig.util as lutil

class RSAGroupOps(object):
    def __init__(self, n, g, h, prng_key=None):
        self.n = n
        self.g = g
        self.h = h
        if prng_key is None:
            prng_key = lutil.rand.getrandbits(1024)
        self.prng = lutil.SHA512PRNG(prng_key)

    def rand_scalar(self):
        pass

    def pow(self, b, e):
        pass

    def pow2(self, b1, e1, b2, e2):
        pass

    def div(self, n, d):
        pass
