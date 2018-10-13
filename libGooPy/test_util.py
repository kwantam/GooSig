#!/usr/bin/python
#
# (C) 2018 Riad S. Wahby <rsw@cs.stanford.edu>

import sys

def show_one_result(name, fails, nreps):
    istr = { False: u" \u2717 ", True: u" \u2713 " }
    if fails is None:
        sys.stdout.write(u"\033[92m%sPASS\033[0m: %s\n" % (istr[True], name))
    elif fails > 0:
        sys.stdout.write(u"\033[91m%sFAIL\033[0m: %s (%d/%d failed)\n" % (istr[False], name, fails, nreps))
    else:
        sys.stdout.write(u"\033[92m%sPASS\033[0m: %s\n" % (istr[True], name))
    sys.stdout.flush()

def show_test(name, just=32):
    sys.stdout.write(("\033[38;5;33m%s\033[0m: " % name).ljust(just))
    sys.stdout.flush()

def show_progress(failed):
    if failed:
        sys.stdout.write('\033[91m.\033[0m')
    else:
        sys.stdout.write('\033[92m.\033[0m')
    sys.stdout.flush()

def run_test(f, nreps):
    ndisp = max(1, nreps >> 6)

    fail_names = f.__doc__.split(",")
    (test_name, fail_names) = (fail_names[0], fail_names[1:])

    show_test(test_name)
    fails = [0] * len(fail_names)
    failed = False
    for idx in range(0, nreps):
        checks = f()
        cidx = 0
        for (cidx, c) in enumerate(checks):
            if not c:
                failed = True
                fails[cidx] += 1
        assert cidx + 1 == len(fails)

        if idx % ndisp == ndisp - 1:
            show_progress(failed)
            failed = False

    sys.stdout.write("\n")
    if all( x == 0 for x in fails ):
        show_one_result("all %s subtests passed (%d)" % (test_name, len(fails)), None, None)
    else:
        for (nf, nn) in zip(fails, fail_names):
            show_one_result("%s_%s" % (test_name, nn), nf, nreps)

    return (sum( 1 for x in fails if x > 0 ), len(fails))

def run_all_tests(nreps, modname, *tests):
    show_test("%s tests" % modname, 0)
    sys.stdout.write("\n")

    fails = subtests = 0
    for test in tests:
        (f, s) = run_test(test, nreps)
        fails += f
        subtests += s

    show_test("Summary", 0)
    sys.stdout.write("\n")
    if fails == 0:
        show_one_result("all %d subtests passed" % subtests, None, None)
    else:
        show_one_result("some subtests did not pass", fails, subtests)
    sys.stdout.write("\n")
