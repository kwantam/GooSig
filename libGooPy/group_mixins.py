#!/usr/bin/python
#
# (C) 2018 Riad S. Wahby <rsw@cs.stanford.edu>

import itertools
import sys

from libGooPy.defs import Defs
import libGooPy.util as lutil

# python 2/3 hack
if sys.version_info[0] == 2:
    zip = itertools.izip    # pylint: disable=redefined-builtin,no-member
    range = xrange          # pylint: disable=redefined-builtin,undefined-variable

class _CombPrecomp(list):
    def __init__(self, g, combspec, gops):
        (ppa, aps, nshifts, bpw, _, size) = combspec
        self.gops = gops
        self.points_per_add = ppa
        self.adds_per_shift = aps
        self.nshifts = nshifts
        self.bits_per_window = bpw
        self.nbits = bpw * ppa
        nskip = 2 ** ppa - 1
        self.points_per_subcomb = nskip

        # allocate space
        list.__init__(self, [None] * size)

        # compute bottom comb
        self[0] = g
        for i in range(1, ppa):
            oval = 1 << i
            ival = oval // 2
            self[oval - 1] = gops.pow(self[ival - 1], None, 1 << bpw)
            for j in range(oval + 1, 2 * oval):
                self[j - 1] = gops.mul(self[j - oval - 1], self[oval - 1])

        # compute upper combs from bottom
        powval = 1 << nshifts
        for i in range(1, aps):
            for j in range(0, nskip):
                self[i * nskip + j] = gops.pow(self[i * nskip + j - nskip], None, powval)

    @staticmethod
    def _from_win(vs):
        ret = 0
        for v in vs:
            ret <<= 1
            ret += 1 if v else 0
        return ret

    def to_comb_exp(self, e):
        ebits = list(lutil.num_to_bits(e, self.nbits))
        nsh = self.nshifts
        nsh_tot = self.adds_per_shift * self.points_per_add
        ebits_split = [ ebits[j*nsh:(j+1)*nsh] for j in range(0, nsh_tot) ]
        iters = []
        for i in reversed(range(0, self.adds_per_shift)):
            iters.append( self._from_win(x) for x in zip(*( ebits_split[i + j*self.adds_per_shift] for j in range(0, self.points_per_add) )) )
        return zip(*iters)

    @staticmethod
    def gen_opt_combs(nbits, maxsize=None):
        # an "optimal" comb for a given #ops is the one that uses the least storage
        #                   for a given storage size is the one that uses the least ops
        opt_combs = {}

        def _gen_comb_result(nshifts, adds_per_shift, points_per_add, bits_per_window):
            # NOTE: this assumes add/mull and double/square have the same cost;
            #       you might adjust this to get a better optimzation result!
            nops = nshifts * (adds_per_shift + 1) - 1
            size = (2 ** points_per_add - 1) * adds_per_shift
            result = (points_per_add, adds_per_shift, nshifts, bits_per_window, nops, size)

            best_so_far = opt_combs.get(nops, None)
            if best_so_far is None or best_so_far[5] > size:
                opt_combs[nops] = result

        for points_per_add in range(2, 18):
            bits_per_window = (nbits + points_per_add - 1) // points_per_add
            for adds_per_shift in range(1, lutil.isqrt(bits_per_window) + 2):
                if bits_per_window % adds_per_shift != 0:
                    # only factorizations of bits_per_window are useful
                    continue
                nshifts = bits_per_window // adds_per_shift

                _gen_comb_result(nshifts, adds_per_shift, points_per_add, bits_per_window)
                _gen_comb_result(adds_per_shift, nshifts, points_per_add, bits_per_window)

        ret_all = []
        ret = None
        sm = None
        for nops in sorted(opt_combs.keys()):
            opt_comb_val = opt_combs[nops]

            if sm is not None and sm <= opt_comb_val[5]:
                continue

            sm = opt_comb_val[5]
            ret_all.append(opt_comb_val)

            if ret is None and maxsize is not None and opt_comb_val[5] <= maxsize:
                ret = opt_comb_val
                break

        if maxsize is None:
            return ret_all
        return ret

class _CombMixin(object):
    def __init__(self, max_comb_size, modbits):
        # combs for g and h
        # NOTE: you really want to store all combs on disk!
        #       I'd recommend having a telescope of combs supporting up to (say)
        #       8192-bit RSA keys (i.e., a ~(2 * 8192 + chalbits + 1) sized comb)
        #       sized for the group of unknown order (see big_prod_size below)
        #
        # P needs two comb sizes, V needs one
        #
        # figure out comb sizes
        if modbits is not None:
            # largest exponent P handles is the greater of
            #       chalbits + 2 * log2(P's RSA modulus) + 1
            #       chalbits + log2(P's RSA modulus) + log2(n) + 1
            big_nbits = max(2 * modbits, modbits + self.nbits_rand) + Defs.chalbits + 1
            big_combspec = _CombPrecomp.gen_opt_combs(big_nbits, max_comb_size)
            small_nbits = self.nbits_rand
            small_combspec = _CombPrecomp.gen_opt_combs(small_nbits, max_comb_size)
            self.combs = [(_CombPrecomp(self.g, small_combspec, self), _CombPrecomp(self.h, small_combspec, self))
                         ,(_CombPrecomp(self.g, big_combspec, self), _CombPrecomp(self.h, big_combspec, self))
                         ]
        else:
            tiny_nbits = Defs.chalbits
            tiny_combspec = _CombPrecomp.gen_opt_combs(tiny_nbits, max_comb_size)
            self.combs = [(_CombPrecomp(self.g, tiny_combspec, self), _CombPrecomp(self.h, tiny_combspec, self))]

    def powgh(self, e1, e2):
        loge = max(int(e1).bit_length(), int(e2).bit_length())
        gcomb = hcomb = None
        for c in self.combs:
            if loge <= c[0].nbits:
                (gcomb, hcomb) = c
                break
        if gcomb is None or hcomb is None:
            raise ValueError("got unexpectedly large exponent in powgh")

        ret = self.id
        for (e1vs, e2vs) in zip(gcomb.to_comb_exp(e1), hcomb.to_comb_exp(e2)):
            if ret != self.id:
                ret = self.sqr(ret)

            for (idx, (e1v, e2v)) in enumerate(zip(e1vs, e2vs)):
                if e1v != 0:
                    ret = self.mul(ret, gcomb[idx * gcomb.points_per_subcomb + e1v - 1])
                if e2v != 0:
                    ret = self.mul(ret, hcomb[idx * hcomb.points_per_subcomb + e2v - 1])

        return ret

class _WnafMixin(object):
    def __init__(self, cheap_inv):
        if cheap_inv:
            self._one_mul = self._one_mul_cheapinv
            self._precomp_wnaf = self._precomp_wnaf_cheapinv
        else:
            self._one_mul = self._one_mul_expinv
            self._precomp_wnaf = self._precomp_wnaf_expinv

    def _wnaf_pc_help(self, b, winsize):
        tablen = 1 << (winsize - 2)
        pctab = [None] * tablen
        bSq = self.sqr(b)
        pctab[0] = b
        for i in range(1, tablen):
            pctab[i] = self.mul(pctab[i-1], bSq)
        return pctab

    def _one_mul_cheapinv(self, ret, w, pctabP, _):
        if w > 0:
            ret = self.mul(ret, pctabP[(w-1)//2])
        elif w < 0:
            ret = self.mul(ret, self.inv(pctabP[(-1-w)//2]))
        return ret

    def _one_mul_expinv(self, ret, w, pctabP, pctabN):
        if w > 0:
            ret = self.mul(ret, pctabP[(w-1)//2])
        elif w < 0:
            ret = self.mul(ret, pctabN[(-1-w)//2])
        return ret

    def _precomp_wnaf_cheapinv(self, b, _, winsize):
        return (self._wnaf_pc_help(b, winsize), None)

    def _precomp_wnaf_expinv(self, b, bInv, winsize):
        return (self._wnaf_pc_help(b, winsize), self._wnaf_pc_help(bInv, winsize))

    @staticmethod
    def _wnaf(r, w, bitlen=None):
        if bitlen is None:
            bitlen = r.bit_length() + 1
        out = [None] * bitlen
        for i in reversed(range(0, bitlen)):
            val = 0
            if r % 2:
                val = r & ((1<<w)-1)
                if val & (1<<(w-1)):
                    val -= (1<<w)
                r -= val
            out[i] = val
            r = r >> 1
        assert r == 0
        return out

    def pow(self, b, bInv, e):
        (pctabP, pctabN) = self._precomp_wnaf(b, bInv, Defs.winsize)
        ebits = self._wnaf(e, Defs.winsize)

        ret = self.id
        for w in ebits:
            if ret != self.id:
                ret = self.sqr(ret)
            ret = self._one_mul(ret, w, pctabP, pctabN)

        return ret

    def pow2(self, b1, b1Inv, e1, b2, b2Inv, e2):
        (pctabP1, pctabN1) = self._precomp_wnaf(b1, b1Inv, Defs.winsize)
        (pctabP2, pctabN2) = self._precomp_wnaf(b2, b2Inv, Defs.winsize)

        totlen = max(e1.bit_length(), e2.bit_length()) + 1
        e1bits = self._wnaf(e1, Defs.winsize, totlen)
        e2bits = self._wnaf(e2, Defs.winsize, totlen)

        ret = self.id
        for (w1, w2) in zip(e1bits, e2bits):
            if ret != self.id:
                ret = self.sqr(ret)
            ret = self._one_mul(ret, w1, pctabP1, pctabN1)
            ret = self._one_mul(ret, w2, pctabP2, pctabN2)

        return ret

class _RandMixin(object):
    def __init__(self, nbits, prng):
        self.nbits_rand = nbits
        if prng is None:
            self.prng = lutil.rand
        else:
            self.prng = prng

    def rand_scalar(self):
        return self.prng.getrandbits(self.nbits_rand)
