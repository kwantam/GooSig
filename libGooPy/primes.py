#!/usr/bin/python
#
# (C) 2018 Dan Boneh, Riad S. Wahby <rsw@cs.stanford.edu>

import heapq
from itertools import cycle

import libGooPy.util as lutil

class PrimeDefs(object):
    wheel_incs = [2,4,2,4,6,2,6,4,2,4,6,6,2,6,4,2,6,4,6,8,4,2,4,2,4,8,6,4,6,2,4,6,2,6,6,4,2,4,6,2,6,4,2,4,2,10,2,10]
    wheel_ps = [2, 3, 5, 7]
    test_primes = ()        # gets set below (after primes() is defined)

# a "prime wheel" of circumference 210. Skips all multiples of 2, 3, 5, 7
def _wheel():
    incs = cycle(PrimeDefs.wheel_incs)
    nval = 11
    while True:
        ret = nval
        inc = next(incs, None)
        nval += inc
        yield ret

# A lazy Sieve of Eratosthenes based on
#     O'Neill, M. E. "The Genuine Sieve of Eratosthenes." J Functional Programming,
#     Vol 19, Iss 1, 2009, pp. 95--106.
def primes():
    w = _wheel()
    f = []

    # first, the initial primes of the wheel
    for p in PrimeDefs.wheel_ps:
        yield p

    # then, a prime heap--based iterator
    while True:
        n = next(w, None)
        prime = True
        while f:
            (nx, inc) = f[0]
            if nx > n:
                break
            if nx == n:
                prime = False
            heapq.heapreplace(f, (nx + inc, inc))
        if prime:
            heapq.heappush(f, (n * n, n))
            yield n
PrimeDefs.test_primes = [ next(p__) for p__ in (primes(),) for _ in range(0, 1000) ]

def primes_skip(nskip):
    p = primes()
    for _ in range(0, nskip):
        next(p)
    return p

def is_square(n):
    isqn = lutil.isqrt(n)
    if isqn * isqn == n:
        return True
    return False

def is_prime_lucas(n, nreps):
    half = (n + 1) // 2
    if (2 * half) % n != 1:
        return False

    (d, s) = lutil.factor_twos(n + 1)
    dbits = lutil.num_to_bits(d)
    msbit = next(dbits)
    assert msbit
    dbits = list(dbits)

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
    ilim = 20
    for _ in range(0, nreps):
        while True:
            if i == ilim:
                if is_square(n):
                    return False
            i += 1
            D = pow(-1, i + 1) * (3 + 2 * i)
            if lutil.jacobi(D, n) == -1:
                break
        ilim = -1
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

    (d, r) = lutil.factor_twos(n - 1)

    for _ in range(0, nreps):
        a = int(lutil.rand.randint(2, n - 2))
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

def is_prime_div(n):
    for p in PrimeDefs.test_primes:
        if n % p == 0:
            return False
    return True

# Baillie-PSW primality test (default #reps is overkill)
def is_prime(n, nreps=2):
    return is_prime_div(n) and is_prime_rm(n, 16 * nreps) and is_prime_lucas(n, nreps)

def next_prime(p):
    p |= 1
    while not is_prime(p):
        p += 2
    return p

def primeprod_and_carmichael(nbits):
    p = primes_skip(1)
    prod = 1
    carm = 1
    while True:
        (prod_, carm_) = (prod, carm)
        np = next(p)
        prod *= np
        carm = (carm * (np - 1)) // lutil.gcd(carm, np - 1)
        if prod.bit_length() > nbits:
            return (prod_, carm_)

def find_mindelta(m, maxmult):
    mindelta = 1
    iii = 1
    for i in range(1, maxmult):
        mm = m * i
        delta = ((1 << mm.bit_length()) - mm + 0.0) / (1 << mm.bit_length())
        if delta < mindelta:
            iii = i
            mindelta = delta
    return iii

def gen_ft_prime_opts(nbits, nfix):
    (m, lamm) = primeprod_and_carmichael(nbits - nfix)
    m_multiplier = find_mindelta(m, 1024)
    amax = (1 << nbits) // m
    a_multiplier = find_mindelta(amax, 1024)
    return (m, m_multiplier, lamm, amax, a_multiplier)

# from Fouque and Tibouchi, "Close to uniform prime number generation with fewer random bits."
#      https://eprint.iacr.org/2011/418
def fouque_tibouchi_primegen(opts, rng):
    (m, m_multiplier, lamm, amax, a_multiplier) = opts
    mlimit = m * m_multiplier
    alimit = amax * a_multiplier
    while True:
        u = 1
        b = 0
        while u != 0:
            r = (rng.randrange(mlimit) * u // m_multiplier) % m
            b = (b + r) % m
            u = 1 - pow(b, lamm, m)
            u = u + m if u < 0 else u

        p = 2
        i = 0
        cont = False
        while not is_prime(p):
            if i > amax // 10:
                # did we choose a "bad" b?
                cont = True
                break
            i += 1
            a = rng.randrange(alimit) // a_multiplier
            p = (a * m + b) | 1
        if not cont:
            break
    return p
