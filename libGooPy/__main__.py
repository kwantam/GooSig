#
#
# (C) 2018 Riad S. Wahby <rsw@cs.stanford.edu>

from __future__ import print_function
import sys
import time

try:
    from libGooPy.defs import Defs
    import libGooPy.group_ops as lg
    from libGooPy.sign import GooSigSigner
    import libGooPy.test_util as tu
    import libGooPy.util as lu
    from libGooPy.verify import GooSigVerifier
except:
    print("ERROR: Could not import libGooPy. Try invoking as `python -m libGooPy`.")
    sys.exit(1)

def main(run_submodules, nreps):
    if run_submodules:
        lu.main(nreps)
        lg.main(nreps)

    # reuse Gops throughout. Notice that you can reuse gops for different
    # Signer modulus as long as the *size* of the Signer's modulus stays the same.
    #
    # group ops to support 2048-bit Signer modulus and 4096-bit AOL root cert group
    gops_2048_p = lg.RSAGroupOps(Defs.Gaol, 2048)
    gops_2048_v = lg.RSAGroupOps(Defs.Gaol, None)
    # group ops to support 4096-bit Signer modulus and 2048-bit RSA challenge group
    gops_4096_p = lg.RSAGroupOps(Defs.Grsa2048, 4096)
    gops_4096_v = lg.RSAGroupOps(Defs.Grsa2048, None)

    # keep rough time measuremnts
    pv_times = [([], []), ([], [])]

    def test_sign_verify():
        "sign_and_verify,AOLx2048,RSAx4096"

        res = [None, None]
        for (idx, (plist, msg, gops_p, gops_v)) in \
                enumerate(( (Defs.primes_1024, "AOL group, 2048-bit RSA key. Test #1.", gops_2048_p, gops_2048_v) \
                          , (Defs.primes_2048, "RSA group, 4096-bit RSA key. Test #2.", gops_4096_p, gops_4096_v) \
                          )):
            (p, q) = lu.rand.sample(plist, 2)
            prv = GooSigSigner(p, q, gops_p)
            s = prv.gops.rand_scalar()
            C1 = prv.gops.reduce(prv.gops.powgh(p * q, s))
            start_time = time.time()
            (C2, t, sigma) = prv.sign(C1, s, msg)
            stop_time = time.time()
            pv_times[idx][0].append(stop_time - start_time)

            ver = GooSigVerifier(gops_v)
            start_time = time.time()
            res[idx] = ver.verify((C1, C2, t), msg, sigma)
            stop_time = time.time()
            pv_times[idx][1].append(stop_time - start_time)

        return res

    tu.run_all_tests(nreps, "end-to-end", test_sign_verify)
    tu.show_timing_pair("4096-bit GOO, 2048-bit PK", pv_times[0])
    tu.show_timing_pair("2048-bit GOO, 4096-bit PK", pv_times[1])

if __name__ == "__main__":
    run_all = False
    nr = 16
    for i in range(1, len(sys.argv)):
        if sys.argv[i] == "-a":
            run_all = True
        else:
            try:
                nr = int(sys.argv[i])
            except:
                nr = 16
    main(run_all, nr)
