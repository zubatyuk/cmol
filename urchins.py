#!/usr/bin/env python

from cMol import *
import os
import openbabel as ob
import sys
import types


def read_symm(fil):
    assert os.path.isfile(fil), "File not found: %s" % fil
    with open(fil, 'r') as f:
        p = filter(None, (line.rstrip() for line in f))
    symm = list()
    for i in p:
        ii = i.split()
        assert len(ii) >= 3, 'Each line of symm file should contain three records reparated by spaces.'
        ii[0] = SymOp(ii[0])
        try:
            ii[1] = int(ii[1])
            ii[2] = int(ii[2])
        except:
            raise ValueError, "Failed to parse line in symm file:\n%s\n" % i
        symm.append(ii)
    return symm

def read_ene(fil):
    assert os.path.isfile(fil), "File not found: %s" % fil
    with open(fil, 'r') as f:
        p = filter(None, (line.rstrip() for line in f))
    ene = list()
    for i in p:
        try:
            i = float(i)
        except:
            raise ValueError, "Failed to parse line in ene file:\n%s\n" % i
        # ignore zero and positive values, just replace with small negative
        if i >= 0:
            i = -0.05
        ene.append(i)
    return ene

def calc_center(vecs):
    center = pybel.ob.vector3(0, 0, 0)
    l = 0
    for v in vecs:
        center += v
        l += 1
    center /= l
    return center

def get_evd_centers(cmol):
    centers = list()
    for m in cmol.mol_map:
        centers.append(calc_center([cmol.c2f(cmol.atoms[a].GetVector()) for a in m]))
    return centers

def get_evd_vectors(centers,symm,ene):
    vectors = list()
    en_scale = 1 / float(2 * min(ene))
    for i in xrange(len(symm)):
        s = symm[i][0]
        e = ene[i]
        v0 = centers[symm[i][1]]
        v = s.apply(centers[symm[i][2]])
        v -= v0
        v *= (en_scale * e)
        v += v0
        vectors.append(v)
    return vectors

def make_evd_struct(cmol,centers,vectors,symm,evd=True,mol=True):
    evdmol = ob.OBMol()
    evdmol.CloneData(cmol.OBMol.GetData(ob.UnitCell))
    if mol:
        for m in cmol.mol_map:
            for i in m:
                evdmol.AddAtom(cmol.OBMol.GetAtom(i+1))
        evdmol.ConnectTheDots()
    if evd:
        nat=evdmol.NumAtoms()
        for c in centers:
            atom = ob.OBAtom()
            atom.SetAtomicNum(0)
            atom.SetVector(cmol.f2c(c))
            evdmol.AddAtom(atom)
        for i in xrange(len(vectors)):
            atom = ob.OBAtom()
            atom.SetAtomicNum(0)
            atom.SetVector(cmol.f2c(vectors[i]))
            evdmol.AddAtom(atom)
        for i in xrange(len(vectors)):
            bond = ob.OBBond()
            bond.SetBegin(evdmol.GetAtom(symm[i][1]+nat+1))
            bond.SetEnd(evdmol.GetAtom(i+len(centers)+nat+1))
            bond.SetBondOrder=1
            evdmol.AddBond(bond)
    return evdmol

def writepdb(mol,filename):
    obconversion = ob.OBConversion()
    formatok = obconversion.SetOutFormat('pdb')
    obconversion.WriteFile(mol,filename)
    obconversion.CloseOutFile()

def main(options, args):
    basename = args[0]
    symm = read_symm(basename + '.symm')
    ene = read_ene(basename + '.ene')
    assert len(symm) == len(ene), "Different number of records in symm and ene files!"
    cmol = read_cmol(basename)

    centers = get_evd_centers(cmol)
    vectors = get_evd_vectors(centers,symm,ene)

    #molecule
    mol=make_evd_struct(cmol,centers,vectors,symm,evd=False,mol=True)
    writepdb(mol,basename+'_m.pdb')
    #molrecule+evd
    mol=make_evd_struct(cmol,centers,vectors,symm,evd=True,mol=True)
    writepdb(mol,basename+'_c.pdb')
    #evd
    mol=make_evd_struct(cmol,centers,vectors,symm,evd=True,mol=False)
    writepdb(mol,basename+'.pdb')

    sys.exit(0)

if __name__ == '__main__':
    from optparse import OptionParser

    usage = "%prog [-h|--help] <name>"
    description = "Script for constructions of EVD. Reads files <name>.cif, <name>.symm and <name>.ene. Writes PDB files."

    parser = OptionParser(usage=usage, description=description)

    (options, args) = parser.parse_args()

    # arg checks
    if not len(args) == 1:
        parser.error('Single argument <name> required.\nGet more help with -h option.')

    for ext in ('cif', 'symm', 'ene'):
        f = args[0] + '.' + ext
        if not os.path.isfile(f):
            parser.error('File not found: %s' % f)

    # go..
    main(options, args)



