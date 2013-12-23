#!/usr/bin/env python

from cMol import *
import pybel, types

def calc_center(vecs):
    center=pybel.ob.vector3(0,0,0)
    l=0
    for v in vecs:
        center+=v
        l+=1
    center/=l
    return center

def read_symm(fil):
    assert os.path.isfile(fil), "File not found: %s" %fil
    with open(fil,'r') as f:
        p=filter(None,(line.rstrip() for line in f))
    symm=list()
    for i in p:
        ii=i.split()
        assert len(ii)>=3, 'Each line of symm file should contain three records reparated by spaces.'
        ii[0]=SymOp(ii[0])
        try:
            ii[1]=int(ii[1])
            ii[2]=int(ii[2])
        except:
            raise ValueError, "Failed to parse line in symm file:\n%s\n" %i
        symm.append(ii)
    return symm

def read_ene(fil):
    assert os.path.isfile(fil), "File not found: %s" %fil
    with open(fil,'r') as f:
        p=filter(None,(line.rstrip() for line in f))
    ene=list()
    for i in p:
        try:
            i=float(i)
        except:
            raise ValueError, "Failed to parse line in ene file:\n%s\n" %i
        #ignore zero and positive values, just replace with small negative
        if i>=0:
            i=-0.05
        ene.append(i)
    return ene

def resfile_atom(name, ntype, x, y, z):
    return '%-4s %1i %7.4f %7.4f %7.4f 1.0 0.05\n' %(name,ntype,x,y,z)

def main(options,args):
    basename=args[0]
    symm=read_symm(basename+'.symm')
    ene=read_ene(basename+'.ene')
    assert len(symm)==len(ene), "Different number of records in symm and ene files!"
    cmol=read_cmol(basename)

    centers=list()
    for m in cmol.mol_map:
        centers.append(calc_center([cmol.c2f(cmol.atoms[a].GetVector()) for a in m]))

    #calculate vectors
    vectors=list()
    en_scale=1/float(2*min(ene))
    for i in xrange(len(symm)):
        s=symm[i][0]
        e=ene[i]
        v0=centers[symm[i][1]]
        v=s.apply(centers[symm[i][2]])
        v-=v0
        v*=(en_scale*e)
        v+=v0
        vectors.append(v)

    #write mol2
    if options.mol2:
        mol=pybel.Molecule(pybel.ob.OBMol())
        for c in centers:
            atom=pybel.ob.OBAtom()
            atom.SetAtomicNum(6)
            atom.SetVector(cmol.f2c(c))
            mol.OBMol.AddAtom(atom)
        for i in xrange(len(vectors)):
            atom=pybel.ob.OBAtom()
            atom.SetAtomicNum(1)
            atom.SetVector(cmol.f2c(vectors[i]))
            mol.OBMol.AddAtom(atom)
            bond=pybel.ob.OBBond()
            #bond.SetBondOrder=1
            bond.SetBegin(mol.atoms[i+len(centers)].OBAtom)
            bond.SetEnd=atom
            mol.OBMol.AddBond(bond)
        mol.write('mol2', basename+'.mol2', overwrite=True)

    #write res
    if options.res:
        with open(basename+'.res','w') as resfile:
            for i in xrange(len(centers)):
                (x,y,z)=(centers[i].GetX(), centers[i].GetY(), centers[i].GetZ())
                idx=i
                resfile.write(resfile_atom('M%i' %idx, 1, x, y, z))
            for i in xrange(len(vectors)):
                (x,y,z)=(vectors[i].GetX(), vectors[i].GetY(), vectors[i].GetZ())
                idx=i+len(centers)
                resfile.write(resfile_atom('M%i' %idx, 2, x, y, z))
            resfile.write('END\n')
            resfile.write('FMOL\n')
            resfile.write('UNDO\n')
            for i in xrange(len(symm)):
                resfile.write('link 1 m%i m%i\n' %(symm[i][1], len(centers)+i))
            resfile.write('END\n')

    sys.exit(0)

if __name__=='__main__':
    import os,sys
    from optparse import OptionParser

    usage = "%prog [-h|--help] [-m|--mol2] [-r|--res] <name>"
    description = "Script for constructions of EVD. Reads files <name>.cif, <name>.symm and <name>.ene. Writes MOL2 or SHELX RES files."

    parser=OptionParser(usage=usage,description=description)

    #parser.add_option
    parser.add_option('-m','--mol2',dest='mol2',help='Write MOL2 file',action='store_true')
    parser.add_option('-r','--res',dest='res',help='Write SHELX RES file',action='store_true')

    (options, args) = parser.parse_args()

    #arg checks
    if not len(args)==1:
        parser.error('Single argument <name> required.\nGet more help with -h option.')
    if not (options.mol2 or options.res):
        parser.error('Nothing to do. Specify either -mol2 or -res')

    for ext in ('cif','symm','ene'):
        f=args[0]+'.'+ext
        if not os.path.isfile(f):
            parser.error('File not found: %s' %f)

    #go..
    main(options,args)


