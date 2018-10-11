#!/usr/bin/python
#
# (C) 2017-18 Riad S. Wahby <rsw@cs.stanford.edu>

import bz2
from collections import deque
import hashlib
import pickle
import sys

import libsig.util as lutil

# python 2/3 hack
if sys.version_info[0] == 2:
    range = xrange      # pylint: disable=redefined-builtin,undefined-variable
else:
    long = int          # pylint: disable=redefined-builtin

class FiatShamir(object):
    VERSION = 500
    ndb = None
    rvstart = None
    rvend = None
    wDiv = None
    prng = None

    def __init__(self, q):
        self.q = q

        nbits = lutil.clog2(q)
        self.mask = (1 << nbits) - 1

        self.rnum = 0
        self.io = deque()
        self.trans = deque()
        self.can_put = True

        self.words_per_val = (nbits + 511) // 512
        self.hash = hashlib.sha512()
        self.hash.update(b"libsig,")

    def hash_update(self, vals):
        for elm in lutil.flatten_iter(vals):
            assert isinstance(elm, (int, long, type(None)))
            self.hash_append(elm)

    def hash_append(self, elm):
        self.hash.update(b"%s," % str(elm).encode("ascii"))
        self.prng = None

    def rand_scalar(self):
        ret = self.q
        while ret >= self.q:
            ret = self._get_rand()
        return ret

    # NOTE you could use AES-CTR instead, but this sticks to Python's stdlib
    def _get_rand(self):
        if self.prng is None:
            self.prng = hashlib.sha512()
            self.prng.update(b"libsig_prng,%s" % self.hash.hexdigest().encode("ascii"))

        val = 0
        for _ in range(0, self.words_per_val):
            val <<= 512
            self.prng.update(b",%d" % self.rnum)
            val += int(self.prng.hexdigest(), 16)
            self.rnum += 1

        val &= self.mask
        return val

    def put(self, vals, is_io=False):
        if not self.can_put:
            assert False, "Cannot put values into an input FS transcript"

        call_items = False
        if isinstance(vals, dict):
            call_items = True
        elif not isinstance(vals, list):
            vals = [vals]

        if is_io:
            self.io.append(vals)
        else:
            self.trans.append(vals)

        if call_items:
            self.hash_update(vals.items())
        else:
            self.hash_update(vals)

    def take(self, take_io=False):
        if self.can_put:
            assert False, "Cannot take values from an output FS transcript"

        if take_io:
            vals = self.io.popleft()
        else:
            vals = self.trans.popleft()

        if isinstance(vals, dict):
            self.hash_update(vals.items())
        else:
            self.hash_update(vals)

        return vals

    @classmethod
    def from_string(cls, string):
        (q, iovals, tvals, ndb, rvstart, rvend, wDiv) = cls.unpack_proof(string)

        ret = cls(q)
        ret.can_put = False
        ret.io.extend(iovals)
        ret.trans.extend(tvals)
        ret.ndb = ndb
        ret.rvstart = rvstart
        ret.rvend = rvend
        ret.wDiv = wDiv

        return ret

    @classmethod
    def proof_size(cls, string):
        (_, _, tvals, _, _, _, _) = cls.unpack_proof(string)

        nelems = len(list(lutil.flatten_iter(tvals)))
        size = len(bz2.compress(pickle.dumps(tvals, -1)))

        return (nelems, size)

    @classmethod
    def unpack_proof(cls, string):
        (q, iovals, tvals, ndb, rvstart, rvend, wDiv, version) = pickle.loads(string)
        assert version == cls.VERSION, "Version: expected %d, got %d" % (cls.VERSION, version)
        return (q, iovals, tvals, ndb, rvstart, rvend, wDiv)

    def to_string(self, full=True):
        if full:
            return pickle.dumps((self.q, self.io, self.trans, self.ndb, self.rvstart, self.rvend, self.wDiv, self.VERSION), -1)
        return pickle.dumps(self.trans, -1)
