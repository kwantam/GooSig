#!/usr/bin/python
#
# (C) 2018 Riad S. Wahby <rsw@cs.stanford.edu>

import itertools
import sys

from libGooPy.defs import Defs
from libGooPy.group_mixins import _CombMixin, _WnafMixin, _RandMixin
import libGooPy.util as lutil

# python 2/3 hack
if sys.version_info[0] == 2:
    zip = itertools.izip    # pylint: disable=redefined-builtin,no-member
    range = xrange          # pylint: disable=redefined-builtin,undefined-variable

class RSAGroupOps(_RandMixin, _WnafMixin, _CombMixin):
    # NOTE you should use an RSA modulus whose factorization is unknown!
    #      In other words, *don't* just generate a modulus! defs.py provides
    #      a few candidates for you to try.
    def __init__(self, Gdesc, modbits=None, prng=None):
        self.n = Gdesc.modulus
        self.nOver2 = self.n // 2
        self.g = Gdesc.g
        self.h = Gdesc.h
        self.id = 1
        self.desc = (self.n, self.g, self.h)

        _RandMixin.__init__(self, lutil.clog2(self.n) - 1, prng)
        _WnafMixin.__init__(self, False)
        _CombMixin.__init__(self, modbits)

    def reduce(self, b):
        # compute the representative of (Z/n)/{1,-1}, i.e., min(|b|, n-|b|)
        if b > self.nOver2:
            return self.n - b
        return b

    def is_reduced(self, b):
        return b <= self.nOver2

    def sqr(self, b):
        return pow(b, 2, self.n)

    def mul(self, m1, m2):
        return (m1 * m2) % self.n

    @lutil.overrides(_WnafMixin)
    def pow(self, b, bInv, e):
        return pow(b, e, self.n)

    def inv(self, b):
        return lutil.invert_modp(b, self.n)

    def inv2(self, b1, b2):
        b12Inv = self.inv(b1 * b2)
        return ((b2 * b12Inv) % self.n, (b1 * b12Inv) % self.n)

    def inv5(self, b1, b2, b3, b4, b5):
        b12 = (b1 * b2) % self.n
        b34 = (b3 * b4) % self.n
        b1234 = (b12 * b34) % self.n
        b12345 = (b1234 * b5) % self.n

        b12345Inv = self.inv(b12345)
        b1234Inv = (b12345Inv * b5) % self.n
        b34Inv = (b1234Inv * b12) % self.n
        b12Inv = (b1234Inv * b34) % self.n

        return ((b12Inv * b2) % self.n, (b12Inv * b1) % self.n, (b34Inv * b4) % self.n, (b34Inv * b3) % self.n, (b12345Inv * b1234) % self.n)

class ClassGroupOps(_RandMixin, _WnafMixin, _CombMixin):
    def __init__(self, Gdesc, modbits=None, prng=None):
        self.D = Gdesc.disc
        assert self.D < 0 and self.D % 4 == 1 and lutil.is_prime(-self.D)
        self.g = Gdesc.g
        self.h = Gdesc.h
        self.id = Gdesc.id
        self.L = Gdesc.L
        assert self.L ** 4 <= -self.D and (self.L + 1) ** 4 > -self.D
        self.desc = (self.D, self.g, self.h)

        # number of random bits for exponents; NOTE shouldn't we halve this?
        _RandMixin.__init__(self, lutil.clog2(self.D) - 1, prng)
        _WnafMixin.__init__(self, True)
        # NOTE hack
        (old_comb_size, Defs.max_comb_size) = (Defs.max_comb_size, 32)
        _CombMixin.__init__(self, modbits)
        Defs.max_comb_size = old_comb_size

    # Algorithm 5.4.2 of Cohen's "A Course in Computational Algebraic Number Theory"
    @staticmethod
    def reduce(f):
        (a, b, c) = f

        while True:
            if -a < b <= a:
                if a > c:
                    b = -b
                    (a, c) = (c, a)
                else:
                    if a == c and b < 0:
                        b = -b
                    return (a, b, c)

            (q, r) = divmod(b, 2*a)
            if r > a:
                r -= 2 * a
                q += 1

            c = c - ((b + r) * q) // 2
            b = r

    @staticmethod
    def is_reduced(f):
        (a, b, c) = f
        return (-a < b <= a < c) or (0 <= b <= a == c)

    # NUCOMP of Daniel Shanks
    # Adapted from
    #   Jacobson, M.J. and van der Poorten, A.J., "Computational Aspects of NUCOMP." Proc. ANTS 2002.
    def mul(self, m1, m2):
        if m1[0] == 1:
            return m2
        if m2[0] == 1:
            return m1

        # unpack, swapping m1 and m2 if w1 < w2
        ((u1, v1, w1), (u2, v2, w2)) = (m2, m1) if m1[2] < m2[2] else (m1, m2)

        # Step 1
        s = (v1 + v2) // 2
        m = v2 - s

        # Step 2
        (c, b, F) = lutil.ext_euclid(u1, u2)
        assert u1 * c + u2 * b == F

        # Steps 2--4
        if s % F == 0:
            G = F
            Bx = m * b
            By = u1 // G
        else:
            (y, G) = lutil.ext_euclid(s, F, left_only=True)
            assert (G - y *s) % F == 0
            H = F // G
            By = u1 // G
            # Step 4
            l1 = (b * (w1 % H)) % H
            l2 = (c * (w2 % H)) % H
            l = (y * (l1 + l2)) % H
            Bx = b * m // H + l * By // H
        Cy = u2 // G
        Dy = s // G

        # Step 5 (truncated Euclidean division)
        bx = Bx % By
        by = By
        x = 1
        y = z = 0
        while abs(by) > self.L and bx != 0:
            ((q, bx), by) = (divmod(by, bx), bx)
            (y, x) = (x, y - q * x)
            z += 1
        (by, y) = (-by, -y) if z % 2 == 1 else (by, y)
        (ax, ay) = (G * x, G * y)

        # Steps 6--7
        if z == 0:
            # Step 6
            Q1 = Cy * bx
            (cx, dx) = ((Q1 - m) // By, (bx * Dy - w2) // By)
            ret = (by * Cy, v2 - 2 * Q1, bx * cx - G * dx)
        else:
            # Step 7
            (cx, dx) = ((Cy * bx - m * x) // By, (Dy * bx - w2 * x) // By)
            (Q1, Q3) = (by * cx, y * dx)
            (Q2, Q4) = (Q1 + m, Q3 + Dy)
            dy = Q4 // x
            cy = Q2 // bx if bx != 0 else (cx * dy - w1) // dx
            ret = (by * cy - ay * dy, G * (Q3 + Q4) - Q1 - Q2, bx * cx - ax * dx)

        assert self.discrim(ret) == self.D
        return self.reduce(ret)

    # NUCOMP of Daniel Shanks
    # Adapted from
    #   Jacobson, M.J. and van der Poorten, A.J., "Computational Aspects of NUCOMP." Proc. ANTS 2002.
    def sqr(self, m):
        if m[0] == 1:
            return m
        (u, v, w) = m

        # Step 1
        (y, G) = lutil.ext_euclid(v, u, left_only=True)
        (By, Dy) = (u // G, v // G)

        # Step 2
        Bx = (y * w) % By

        # Step 3
        bx = Bx
        by = By
        x = 1
        y = z = 0
        while abs(by) > self.L and bx != 0:
            ((q, bx), by) = (divmod(by, bx), bx)
            (y, x) = (x, y - q * x)
            z += 1
        (by, y) = (-by, -y) if z % 2 == 1 else (by, y)
        (ax, ay) = (G * x, G * y)

        # Steps 4--5
        if z == 0:
            # Step 4
            dx = (bx * Dy - w) // By
            (u3, w3) = (by ** 2, bx ** 2)
            ret = (u3, v - (bx + by) ** 2 + u3 + w3, w3 - G * dx)
        else:
            # Step 5
            dx = (bx * Dy - w * x) // By
            Q1 = dx * y
            dy = Q1 + Dy
            v3 = G * (dy + Q1)
            dy = dy // x
            (u3, w3) = (by ** 2, bx ** 2)
            ret = (u3 - ay * dy, v3 - (bx + by) ** 2 + u3 + w3, w3 - ax * dx)

        assert self.discrim(ret) == self.D
        return self.reduce(ret)

    @staticmethod
    def discrim(m):
        (a, b, c) = m
        return b * b - 4 * a * c

    @staticmethod
    def inv(m):
        (a, b, c) = m
        return (a, -b, c)

    @classmethod
    def invAll(cls, *ms):
        return [ cls.inv(m) for m in ms ]
    inv2 = inv5 = classmethod(lambda cls, *ms: tuple( cls.inv(m) for m in ms ))

def main(nreps):
    import libGooPy.consts as lc
    import libGooPy.test_util as tu

    # test on random RSA modulus
    (p, q) = lutil.rand.sample(tu.primes_1024, 2)
    n = p * q
    Grandom = lc.gen_RSA_group_obj(n, 5, 7)

    t1 = RSAGroupOps(lc.Grsa2048, 2048)
    t2 = RSAGroupOps(Grandom, 2048)

    def test_pow2():
        "pow2_wnaf,RSA_chal,RSA_rand"

        (b1, b2, e1, e2) = ( lutil.rand.getrandbits(2048) for _ in range(0, 4) )
        (b1Inv, b2Inv)= t1.inv2(b1, b2)
        out1 = (pow(b1, e1, t1.n) * pow(b2, e2, t1.n)) % t1.n
        t1o = t1.pow2(b1, b1Inv, e1, b2, b2Inv, e2)

        (b1Inv, b2Inv)= t2.inv2(b1, b2)
        out2 = (pow(b1, e1, t2.n) * pow(b2, e2, t2.n)) % t2.n
        t2o = t2.pow2(b1, b1Inv, e1, b2, b2Inv, e2)

        return (out1 == t1o, out2 == t2o)

    def test_powgh():
        "powgh,RSA_chal,RSA_rand"

        (e1, e2) = ( lutil.rand.getrandbits(2 * 2048 + Defs.chalbits + 2) for _ in range(0, 2) )

        out1 = (pow(2, e1, t1.n) * pow(3, e2, t1.n)) % t1.n
        t1o = t1.powgh(e1, e2)

        (e1_s, e2_s) = ( x >> (2048 + Defs.chalbits) for x in (e1, e2) )
        out2 = (pow(5, e1_s, t2.n) * pow(7, e2_s, t2.n)) % t2.n
        t2o = t2.powgh(e1_s, e2_s)

        return (out1 == t1o, out2 == t2o)

    def test_inv2():
        "inv2,RSA_chal,RSA_rand"

        (e1, e2) = ( lutil.rand.getrandbits(2048) for _ in range(0, 2) )
        (e1Inv, e2Inv) = t1.inv2(e1, e2)
        t1pass = t1.reduce((e1 * e1Inv) % t1.n) == 1 and t1.reduce((e2 * e2Inv) % t1.n) == 1

        (e1_s, e2_s) = ( x >> 1536 for x in (e1, e2) )
        (e1_sInv, e2_sInv) = t2.inv2(e1_s, e2_s)
        t2pass = t2.reduce((e1_s * e1_sInv) % t2.n) == 1 and t2.reduce((e2_s * e2_sInv) % t2.n) == 1

        return (t1pass, t2pass)

    def test_inv5():
        "inv5,RSA_chal,RSA_rand"

        eVals = tuple( lutil.rand.getrandbits(2048) for _ in range(0, 5) )
        eInvs = t1.inv5(*eVals)
        t1pass = all( t1.reduce((e * eInv) % t1.n) == 1 for (e, eInv) in zip(eVals, eInvs) )

        eInvs = t2.inv5(*eVals)
        t2pass = all( t2.reduce((e * eInv) % t2.n) == 1 for (e, eInv) in zip(eVals, eInvs) )

        return (t1pass, t2pass)

    tu.run_all_tests(nreps, "group_ops", test_pow2, test_powgh, test_inv2, test_inv5)

if __name__ == "__main__":
    try:
        nr = int(sys.argv[1])
    except:
        nr = 32
    main(nr)
