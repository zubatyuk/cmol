#!/usr/bin/env python

import os
import pybel
import sys

ROUTE = '!b97-d3 def2-tzvp def2-tzvp/j ecp(def2-tzvp,def2-tzvp/j) slowconv nososcf'

def mkinp_orca(mols, bq=(), route=""):
    etab = pybel.ob.OBElementTable()
    res = list()
    res.append(route+"\n")
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
    from optparse import OptionParser
    from glob import glob
    
    usage = "%prog [-h|--help] [-b|--bsse]  <name>"
    
    description = "Script for construction of ORCA inputs from XYZ files. Optionally creates input files for BSSE correction. Reads <name>-d???.xyz and <name>-m??.xyz files."

    parser = OptionParser(usage=usage, description=description)
    
    # parser.add_option
    parser.add_option('-b', '--bsse', action="store_true", dest='bsse', help='Flag to create input files for BSSE correction.', default=False)
    parser.add_option('-r', '--route', dest='route', help='ORCA instructions.', type='str', default=ROUTE)
    
    (options, args) = parser.parse_args()

    # arg checks
    if not len(args) == 1:
        parser.error('Single argument <name> required.\nGet more help with -h option.')
    name=args[0]
        
    #..go
    xyz_d=glob('%s-d???.xyz' %name)
    assert xyz_d, "Could not find files matching mask %s-d???.xyz" %name
    if not options.bsse:
        xyz_m=glob('%s-m??.xyz' %name)
        assert xyz_m, "Could not find files matching mask %s-m??.xyz" %name
    
    if options.bsse:
        #create name-d???-m0-1.inp, name-d???-m0.inp and name-d???-m1.inp
        for xyzfile in xyz_d:
            basename=os.path.splitext(xyzfile)[0]
            mols = [m for m in pybel.readfile('xyz', xyzfile)]
            assert len(mols) == 2, "Expected 2 molecules in file %s, found %i" %(xyzfile, len(mols))
            with open(basename + '-m0-1.inp', 'w') as f:
                f.writelines(mkinp_orca(mols,route=options.route))
            with open(basename + '-m0.inp', 'w') as f:
                f.writelines(mkinp_orca(mols,(1,),route=options.route))
            with open(basename + '-m1.inp', 'w') as f:
                f.writelines(mkinp_orca(mols,(0,),route=options.route))
    if not options.bsse:
        #create name-d???.inp amd name-m??.inp
        for xyzfile in xyz_d+xyz_m:
            basename=os.path.splitext(xyzfile)[0]
            mols = [m for m in pybel.readfile('xyz', xyzfile)]
            with open(basename + '.inp', 'w') as f:
                f.writelines(mkinp_orca(mols,route=options.route))
