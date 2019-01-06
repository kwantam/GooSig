#!/usr/bin/python
#
# (C) 2018 Dan Boneh, Riad S. Wahby <rsw@cs.stanford.edu>

from __future__ import print_function
import sys
import time

try:
    import libGooPy.consts as lc
    import libGooPy.group_ops as lg
    import libGooPy.rsa as lrsa
    from libGooPy.sign import GooSigSigner
    import libGooPy.test_util as tu
    from libGooPy.tokengen import GooSigTokGen
    import libGooPy.util as lutil
    from libGooPy.verify import GooSigVerifier
except Exception as e:  # pylint: disable=broad-except
    print("ERROR: Could not import libGooPy. Try invoking as `python -m libGooPy`.")
    print(str(e))
    sys.exit(1)

def main(run_submodules, nreps):
    if run_submodules:
        lutil.main(nreps)
        lg.main(nreps)
        lrsa.main(nreps)

    # reuse Gops throughout. Notice that you can reuse gops for different
    # Signer modulus as long as the *size* of the Signer's modulus is no
    # larger than the one the gops object was built for.
    #
    # 4096-bit GoUO
    gops_4_2_p = lg.RSAGroupOps(lc.Gaol, 2048)          # 4096-bit RSA GoUO, 2048-bit Signer key
    gops_4_4_p = lg.RSAGroupOps(lc.Gaol, 4096)          # 4096-bit RSA GoUO, 4096-bit Signer key
    gops_4_v = lg.RSAGroupOps(lc.Gaol, None)            # 4096-bit RSA GoUO (verification)
    # 2048-bit GoUO
    gops_2_2_p = lg.RSAGroupOps(lc.Grsa2048, 2048)      # 2048-bit RSA GoUO, 2048-bit Signer key
    gops_2_4_p = lg.RSAGroupOps(lc.Grsa2048, 4096)      # 2048-bit RSA GoUO, 4096-bit Signer key
    gops_2_v = lg.RSAGroupOps(lc.Grsa2048, None)        # 2048-bit RSA GoUO (verification)
    # 2048-bit BQF discriminant
    gops_c2_2_p = lg.ClassGroupOps(lc.Ggoo2048, 2048)   # 2048-bit BQF GoUO, 2048-bit Signer key
    gops_c2_4_p = lg.ClassGroupOps(lc.Ggoo2048, 4096)   # 2048-bit BQF GoUO, 4096-bit Signer key
    gops_c2_v = lg.ClassGroupOps(lc.Ggoo2048, None)     # 2048-bit BQF GoUO (verification)
    # 1024-bit BQF discriminant
    gops_c1_2_p = lg.ClassGroupOps(lc.Ggoo1024, 2048)   # 1024-bit BQF GoUO, 2048-bit Signer key
    gops_c1_4_p = lg.ClassGroupOps(lc.Ggoo1024, 4096)   # 1024-bit BQF GoUO, 2048-bit Signer key
    gops_c1_v = lg.ClassGroupOps(lc.Ggoo1024, None)     # 1024-bit BQF GoUO, 4096-bit Signer key

    # measure times
    pv_expts = [ ("4096-bit RSA GoUO, 2048-bit Signer PK", gops_4_2_p, gops_4_v)
               , ("4096-bit RSA GoUO, 4096-bit Signer PK", gops_4_4_p, gops_4_v)
               , ("2048-bit RSA GoUO, 2048-bit Signer PK", gops_2_2_p, gops_2_v)
               , ("2048-bit RSA GoUO, 4096-bit Signer PK", gops_2_4_p, gops_2_v)
               , ("2048-bit BQF GoUO, 2048-bit Signer PK", gops_c2_2_p, gops_c2_v)
               , ("2048-bit BQF GoUO, 4096-bit Signer PK", gops_c2_4_p, gops_c2_v)
               , ("1024-bit BQF GoUO, 2048-bit Signer PK", gops_c1_2_p, gops_c1_v)
               , ("1024-bit BQF GoUO, 4096-bit Signer PK", gops_c1_4_p, gops_c1_v)
               ]
    pv_times = [ ([], [], []) for _ in range(0, len(pv_expts)) ]
    pv_plsts = [tu.primes_1024, tu.primes_2048]

    def test_sign_verify():
        "sign_and_verify,4x2,4x4,2x2,2x4,c2x2,c2x4,c1x2,c1x4"

        res = [None] * len(pv_times)
        for (idx, (msg, gops_p, gops_v)) in enumerate(pv_expts):
            ## create a GooSig for a random Signer key
            (p, q) = lutil.rand.sample(pv_plsts[idx % 2], 2)
            rsakey = lrsa.RSAKey(p, q)              # the Signer's private key
            rsapubkey = rsakey.get_public_key()     # the Signer's public key
            gen = GooSigTokGen(gops_p)
            prv = GooSigSigner(rsakey, gops_p)
            ver = GooSigVerifier(gops_v)

            # generate the challenge token
            start_time = time.time()
            (C0, C1) = gen.send_tokens(rsapubkey)
            stop_time = time.time()
            pv_times[idx][0].append(stop_time - start_time)

            # generate the signature
            start_time = time.time()
            ## UPDATE 2019 Jan 06: unpack C3 from sign()'s return value; sigma also contains one new element
            (C2, C3, t, sigma) = prv.sign(C0, C1, msg)
            stop_time = time.time()
            pv_times[idx][1].append(stop_time - start_time)

            # verify the signature
            start_time = time.time()
            res[idx] = ver.verify((C1, C2, t), msg, sigma)
            stop_time = time.time()
            pv_times[idx][2].append(stop_time - start_time)

        return res

    tu.run_all_tests(nreps, "end-to-end", test_sign_verify)
    if sys.flags.optimize == 0:
        tu.show_warning("you should call Python with the -O flag to get reasonable timing numbers.")
    for (idx, (n, _, _)) in enumerate(pv_expts):
        tu.show_timing_triple(n, pv_times[idx])

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
