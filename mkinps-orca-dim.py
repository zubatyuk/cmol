#!/usr/bin/env python

import os
import pybel
import sys

ROUTE = '!blyp def2-tzvp def2-tzvp/j ecp{def2-tzvp,def2-tzvp/j} vdw10 nososcf\n'

def mkinp_orca(mols, bq=()):
    etab = pybel.ob.OBElementTable()
    res = list()
    res.append(ROUTE)
    res.append('*xyz 0 1\n')
    for i in range(len(mols)):
        if i in bq: bqflag = ':'
        else: bqflag = ' '
        for j in range(len(mols[i].atoms)):
            atom = mols[i].atoms[j].OBAtom
            res.append('%2s %s % 10.6f % 10.6f % 10.6f\n' % (etab.GetSymbol(atom.GetAtomicNum()),
                                    bqflag, atom.GetX(), atom.GetY(), atom.GetZ()))
    res.append('*\n')
    return res

if __name__ == '__main__':
    xyzfile = sys.argv[1]  # should contain separated mols
    basename = os.path.splitext(xyzfile)[0]

    mols = [m for m in pybel.readfile('xyz', xyzfile)]

    # dimers
    f = file(basename + '_m0.inp', 'w')
    f.writelines(mkinp_orca(mols, (1,)))
    f.close()
    f = file(basename + '_m1.inp', 'w')
    f.writelines(mkinp_orca(mols, (0,)))
    f.close()
    f = file(basename + '_m0-1.inp', 'w')
    f.writelines(mkinp_orca(mols))
    f.close()
