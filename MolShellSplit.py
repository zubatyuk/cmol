#!/usr/bin/env python

from cMol import *
import os
import sys


def main(options,args):
    basename=args[0]

    cmol=read_cmol(basename)
    if not cmol: sys.exit()

    convs=ob.OBConversion()
    convs.SetOutFormat('xyz')
    conv=ob.OBConversion()
    conv.SetOutFormat('xyz')

    with open(basename+'.symm'+) as symm:
        id=0
        for id0 in xrange(len(cmol.mol_map)):
            m0=ob.OBMol()
            for i in cmol.mol_map[id0]:
                m0.AddAtom(cmol.atoms[i])
            ms=list()
            for s,id1 in cmol.iter_close(id0):
                symm.write('%s %s %s\n' %(s,id0,id1 ))
                id+=1
                m1=ob.OBMol()
                atoms1=[cmol.atoms[i] for i in cmol.mol_map[id1]]
                for a in atoms1:
                    na=ob.OBAtom()
                    na.SetVector(cmol.f2c(s.apply(cmol.c2f(a.GetVector()))))
                    na.SetAtomicNum(a.GetAtomicNum())
                    m1.AddAtom(na)
                conv.WriteFile(m0, '%s-d%03i.xyz' %(basename,id))
                conv.Write(m1)
                conv.CloseOutFile()
                ms.append(m1)
            convs.WriteFile(m0,'%s-s%02i.xyz' %(basename,id0))
            for m in ms:
                convs.Write(m)
            convs.CloseOutFile()

if __name__ == '__main__':
    from optparse import OptionParser

    usage = "%prog [-h|--help] [-w|--wdvinc 1.0]  <name>"
    description = "Script for constructions of molecular shells from CIF. Reads files <name>.cif. Writes symm and xyz files."

    parser = OptionParser(usage=usage, description=description)

    # parser.add_option
    parser.add_option('-w', '--wdvinc', dest='wdvinc', help='Find molecules with close contact of WdV radii sum + vdwinc', type='float')

    (options, args) = parser.parse_args()

    # arg checks
    if not len(args) == 1:
        parser.error('Single argument <name> required.\nGet more help with -h option.')

    # go..
    main(options, args)





