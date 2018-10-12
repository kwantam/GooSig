#!/usr/bin/python

import sys

def main(nreps):
    import libsig.util as lu
    import libsig.group_ops as lg

    lu.main(16*nreps)
    lg.main(nreps)

if __name__ == "__main__":
    try:
        nr = int(sys.argv[1])
    except:
        nr = 16
    main(nr)
