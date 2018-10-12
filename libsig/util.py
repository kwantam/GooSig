#!/usr/bin/python
#
# (C) 2018 Riad S. Wahby <rsw@cs.stanford.edu>

import math
import random
import sys

# python 2/3 hack
if sys.version_info[0] == 2:
    range = xrange      # pylint: disable=redefined-builtin,undefined-variable

rand = random.SystemRandom()

# ceiling of log2
def clog2(val):
    if isinstance(val, float):
        val = int(math.ceil(val))
    return (val-1).bit_length()

def invert_modp(val, prime):
    if val % prime == 0:
        return 0
    (inv, _) = ext_euclid(val % prime, prime)
    if (inv * val - 1) % prime != 0:
        assert False
    return inv % prime

def ext_euclid(a, b):
    s  = t_ = 0
    s_ = t  = 1
    r  = a
    r_ = b

    while r != 0:
        q = r_ // r
        (r_, r) = (r, r_ - q * r)
        (s_, s) = (s, s_ - q * s)
        (t_, t) = (t, t_ - q * t)

    return (t_, s_)

# compute Jacobi symbol, n prime or composite
def jacobi(a, n):
    if n <= 0 or n % 2 == 0:
        return 0

    negate = False
    a = a % n
    while a != 0:
        while a % 2 == 0:
            a = a // 2
            if n % 8 == 3 or n % 8 == 5:
                negate = not negate

        if a % 4 == 3 and n % 4 == 3:
            negate = not negate

        (n, a) = (a, n)

        a = a % n

    if n == 1:
        return -1 if negate else 1

    return 0

def isqrt(n):
    if n < 0:
        return None
    if n < 2:
        return n

    smallCand = isqrt(n >> 2) << 1
    largeCand = smallCand + 1
    if largeCand * largeCand > n:
        return smallCand
    return largeCand

def is_square(n):
    isqn = isqrt(n)
    if isqn * isqn == n:
        return True
    return False

def num_to_bits(n):
    nbits = int(n).bit_length()
    res = [None] * nbits
    for idx in range(0, nbits):
        res[nbits - 1 - idx] = True if n & 1 else False
        n = n >> 1
    return res

def factor_twos(n):
    d = n
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1
    return (d, s)

def is_prime_lucas(n, nreps):
    if is_square(n):
        return False

    half = invert_modp(2, n)
    if (2 * half) % n != 1:
        return False

    (d, s) = factor_twos(n + 1)
    dbits = num_to_bits(d)
    assert dbits[0]
    dbits = dbits[1:]

    def lucas_double(u, v, Qk):
        u = (u * v) % n
        v = (v ** 2 - 2 * Qk) % n
        Qk = (Qk * Qk) % n
        return (u, v, Qk)

    def lucas_add1(u, v, Qk, D):
        (u, v) = (((u + v) * half) % n, ((D * u + v) * half) % n)
        Qk = (Qk * Q) % n
        return (u, v, Qk)

    i = 0
    for _ in range(0, nreps):
        while True:
            i += 1
            D = pow(-1, i + 1) * (3 + 2 * i)
            if jacobi(D, n) == -1:
                break
        Q = (1 - D) // 4

        (u, v, Qk) = (1, 1, Q)
        for db in dbits:
            (u, v, Qk) = lucas_double(u, v, Qk)
            if db:
                (u, v, Qk) = lucas_add1(u, v, Qk, D)

        # now we have Ud and Vd
        if u % n == 0:
            continue

        # check V_{d*2^r}, 0 <= r < s
        cont = False
        for _ in range(0, s):
            if v % n == 0:
                cont = True
                break
            (u, v, Qk) = lucas_double(u, v, Qk)

        if cont:
            continue

        return False

    return True

# Rabin-Miller
def is_prime_rm(n, nreps):
    if n < 7:
        if n in (3, 5):
            return True
        return False

    (d, r) = factor_twos(n - 1)

    for _ in range(0, nreps):
        a = int(rand.randint(2, n - 2))
        x = pow(a, d, n)

        if x in (1, n-1):
            continue

        cont = False
        for _ in range(0, r - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                cont = True
                break

        if cont:
            continue

        return False

    return True

# Baillie-PSW primality test (default #reps is massive overkill)
def is_prime(n, nreps=8):
    return is_prime_rm(n, 16 * nreps) and is_prime_lucas(n, nreps)

# Tonelli-Shanks
def sqrt_modp(n, p):
    n = n % p
    if n == 0:
        return 0

    if jacobi(n, p) == -1:
        return None

    if p % 4 == 3:
        return pow(n, (p + 1) // 4, p)

    # factor out 2^s from p - 1
    (Q, s) = factor_twos(p - 1)

    # find a non-residue mod p
    w = 2
    while jacobi(w, p) != -1:
        w += 1

    w = pow(w, Q, p)
    y = pow(n, Q, p)
    q = pow(n, (Q + 1) // 2, p)

    while True:
        y_save = y
        i = 0
        while i < s and y != 1:
            y = pow(y, 2, p)
            i += 1
        if i == 0:
            break
        if i == s:
            return None

        w = pow(w, 1 << (s - i - 1), p)
        s = i
        q = (q * w) % p
        w = pow(w, 2, p)
        y = (y_save * w) % p

    if q > (p // 2):
        q = p - q

    assert n == (q * q) % p
    return q

# sqrt mod a product of primes
def sqrt_modn(x, p, q):
    sqrtP = sqrt_modp(x, p)
    sqrtQ = sqrt_modp(x, q)

    if sqrtP is None or sqrtQ is None:
        return None

    (mP, mQ) = ext_euclid(p, q)

    return (sqrtQ * mP * p + sqrtP * mQ * q) % (p * q)

def random_prime(nbits, rng=None):
    p = 1
    while p.bit_length() != nbits or not is_prime(p):
        p = rng.getrandbits(nbits) | 1
        while p.bit_length() == nbits and not is_prime(p):
            p += 2
    return p
