#!/usr/bin/env python

import sys, os

def get_sp_en_orca(fil):
    if not os.path.isfile(fil):
        raise IOError("Could not open file %s" %fil)
    f=file(fil)
    with open(fil) as f:
        for l in f:
            if l.startswith('FINAL SINGLE POINT ENERGY'):
                s=l.split()
                return float(s[4])
    raise IOError('Failed to read energy from file %s' %fil)

def main():
    pass

if __name__ == '__main__':
    from optparse import OptionParser
    from glob import glob

    usage = "%prog [-h|--help] [-b|--bsse]  <name>"
    description = '''Script to get interaction energies rom ORCA output files. Reads <name>.symm and outpus files, writes <name>.ene'''

    parser = OptionParser(usage=usage, description=description)

    # parser.add_option
    parser.add_option('-b', '--bsse', action="store_true", dest='bsse', help='Flag to create input files for BSSE correction.', default=False)

    (options, args) = parser.parse_args()

    # arg checks
    if not len(args) == 1:
        parser.error('Single argument <name> required.\nGet more help with -h option.')
    name=args[0]

    #read .symm file
    symm=list()
    with open("%s.symm" %name) as f:
        for line in f:
            if not line.split()>=3:
                break

            symm.append(line.split())

    #read energies
    en_list=list()
    if options.bsse:
        for i in xrange(len(symm)):
            e=[]
            for suff in ('-m0-1.out', '-m0.out', '-m1.out'):
                e.append(get_sp_en_orca(name+'-d%03i'%(i+1)+suff))
            energy=(e[0]-e[1]-e[2])*627.5095
            en_list.append(energy)

    if not options.bsse:
        ene_m={}
        for i in xrange(len(symm)):
            for m in symm[i][1], symm[i][2]:
                if not ene_m.has_key(m):
                    ene_m[m]=get_sp_en_orca(name+'-m%02i.out' %int(m))
            e0=get_sp_en_orca(name+'-d%03i.out'%(i+1))
            energy=(e0-ene_m[symm[i][1]]-ene_m[symm[i][2]])*627.5095
            en_list.append(energy)

    with open (name+'.ene','w') as ene:
        for i in en_list:
            ene.write(str(i)+'\n')
