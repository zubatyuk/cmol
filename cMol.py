#!/usr/bin/env python

from __future__ import division

from fractions import Fraction
from libcif import parsecif, readcif
import openbabel as ob
import os
import pybel
import sys
import types


def get_symops(cifname):
    symops = list()
    cif = parsecif(readcif(cifname))
    for l in cif['loop_']:
        for t in '_space_group_symop_operation_xyz', '_symmetry_equiv_pos_as_xyz':
            if t in l[0]:
                for s in l[1]:
                    symops.append(SymOp(s[t]))
            continue
    return symops

def read_cmol(basename):
    if not os.path.isfile(basename + '.cif'):
        sys.stderr.write('No CIF to work with!\n')
        return None
    if os.path.isfile(basename + '.map'):
        '''
        file should look like:
        2
        (SymOp('1/2-x,1-y,1/2+z'),)
        ((0,1),)
        '''
        separate_args = [eval(l.rstrip().rstrip()) for l in file(basename + '.map').readlines() if l]
    else:
        separate_args = ()
    symops = get_symops(basename + '.cif')
    mol = pybel.readfile('cif', basename + '.cif').next()
    return cMol(mol.OBMol, mol.unitcell, symops, separate_args=separate_args)

def printvec(v):
    return str((v.GetX(), v.GetY(), v.GetZ()))

class SymOp:
    def __init__(self,r=ob.transform3d(1),t=ob.vector3(0,0,0)):
        self.T=ob.vector3(0,0,0)
        if type(r) is types.StringType:
            self.set_from_str(r)
        else:
            self.R=r
        self.add_trans(t)

    def __eq__(self,other):
        assert type(other)==type(self), '%s is not %s !' %(type(other),type(self))
        return self.R.DescibeAsString() == other.R.DescibeAsString() and \
               self.T.IsApprox(other.T,0.01)
            
    def __repr__(self):
        #describe as 4x4 matrix
        v=[float(i) for i in self.R.DescribeAsValues().split()]
        #construct representation of symmetry operation without translations
        s=ob.transform3d(ob.vector3(v[0],v[1],v[2]),ob.vector3(v[4],v[5],v[6]), 
               ob.vector3(v[8],v[9],v[10]),ob.vector3(0)).DescribeAsString().split(',')
        #add translational part from 4x4 matrix to true translation from
        #symmetry operation
        t=[v[3]+self.T.GetX(),v[7]+self.T.GetY(),v[11]+self.T.GetZ()]
        #construct scring
        for i in xrange(3):
            t[i]=str(Fraction.from_float(t[i]).limit_denominator(6))
            t[i]+='+'+s[i]
            t[i]=t[i].replace('+-','-')
            if t[i].startswith('0'):
                t[i]=t[i][1:]
            if t[i].startswith('+'):
                t[i]=t[i][1:]
        return ','.join(t)
        
    def set_from_str(self,s):
        assert type(s) is types.StringType, 's is not string: %s' %s
        assert 'x' in s and 'y' in s and 'z' in s, 's sould containt x,y,z: %s' %s
        (x,y,z)=(0,0,0)
        s0=eval(s)
        vt2=[int(i) for i in s0]
        vt1=[s0[i]-vt2[i] for i in range(3)]
        r=list()
        for (x,y,z) in ((1,0,0),(0,1,0),(0,0,1)):
            s1=eval(s)
            r.append([s1[i]-s0[i] for i in range(3)])
        self.R=ob.transform3d(ob.vector3(r[0][0],r[1][0],r[2][0]),
                              ob.vector3(r[0][1],r[1][1],r[2][1]),
                              ob.vector3(r[0][2],r[1][2],r[2][2]),
                              ob.vector3(vt1[0],vt1[1],vt1[2]))
        r1=[float(i) for i in self.R.DescribeAsValues().split()]
        self.T=ob.vector3(s0[0]-r1[3],s0[1]-r1[7],s0[2]-r1[11])
        #print map(list,zip(*r)),vt1,vt2,s, str(self)
        #assert str(self) == s.replace(' ','',99)

    def is_identity(self):
        v0=ob.vector3(.1,.2,.3)
        v=self.apply(v0)
        if v.IsApprox(v0,0.01):
            return True
        return False
        
    def add_trans(self,t):
        self.T+=t
    
    def set_trans(self,t):
        self.T=t
        
    def apply(self,v):
        res=self.R*v
        res+=self.T
        return res

class cMol(object):
    def __init__(self, OBMol, unitcell, symops=list(), separate_args=(1,)):
        self.OBMol = OBMol
        self.unitcell = unitcell
        self.f2c = self.unitcell.FractionalToCartesian
        self.c2f = self.unitcell.CartesianToFractional
        self.symops = symops

        self.norm_h()
        self.separate(separate_args)
        #self.move_to_center()

        self.vc = [a.GetVector() for a in self.atoms]
        self.vf = [self.c2f(v) for v in self.vc]
        obet = ob.OBElementTable()
        self.vdw = [obet.GetVdwRad(a.GetAtomicNum()) for a in self.atoms]

    def iter_trans2(self):
        d = (0, 1, -1, 2, -2, 3, -3)
        for i in d:
            for j in d:
                for k in d:
                    yield ob.vector3(i, j, k)

    def iter_symop_t2(self):
        for s in self.symops:
            for t in self.iter_trans2():
                yield SymOp(s.R, t)

    @property
    def center_c(self):
        c = ob.vector3(0, 0, 0)
        n = 0
        for atom in self.atoms:
            c += atom.GetVector()
            n += 1
        return ob.vector3(c.GetX() / n, c.GetY() / n, c.GetZ() / n)

    @property
    def center_f(self):
        return self.c2f(self.center_c)

    def move_to_center(self):
       c = self.center_f
       vc = list()
       for i in (c.GetX(), c.GetY(), c.GetZ()):
           if i < 0: i -= 1
           vc.append(int(i))
       t = ob.vector3(vc[0], vc[1], vc[2])
       if not t.IsApprox(ob.vector3(0, 0, 0), 0.001):
            print "WARNING: moving by [ %4.1f %4.1f %4.1f ]" %(vc[0], vc[1], vc[2])
            tc=self.f2c(t)
            for i in xrange(self.OBMol.NumAtoms()):
                atom=self.OBMol.GetAtom(i+1)
                v = ob.vector3(atom.GetX(), atom.GetY(), atom.GetZ()) #wtf! GetVector fails with segfault
                v -= tc
                atom.SetVector(v)
    @property
    def atoms(self):
        return [self.OBMol.GetAtom(i + 1) for i in range(self.OBMol.NumAtoms())]

    def norm_h(self):
        _bond_data = {6: 1.089, 7: 1.015, 8: 0.993}
        for bond in ob.OBMolBondIter(self.OBMol):
            atoms = (bond.GetBeginAtom(), bond.GetEndAtom())
            anums = (atoms[0].GetAtomicNum(), atoms[1].GetAtomicNum())
            if anums[0] == 1: bond.SetLength(atoms[1], _bond_data[anums[1]])
            if anums[1] == 1: bond.SetLength(atoms[0], _bond_data[anums[0]])


    def set_mon_map(self, mons=[]):
        self.mol_map = [list() for i in mons]
        for mon in mons:
            for i in range(mon.NumAtoms()):
                v1 = mon.GetAtom(i + 1).GetVector()
                for j in range(self.OBMol.NumAtoms()):
                    v2 = self.OBMol.GetAtom(j + 1).GetVector()
                    if v1.IsApprox(v2, 0.001):
                        self.mol_map[mons.index(mon)].append(j)
                        continue

    def separate(self, args):
        try:        mode = args[0]
        except:     mode = 1
        try:        symops = args[1]
        except:     symops = list()
        try:        molmap = args[2]
        except:     molmap = list()


        if mode == 0:
            # no separation
            self.set_mon_map((self.OBMol,))
        if mode == 1:
            # OBMol.Separate
            self.set_mon_map(self.OBMol.Separate())
        if mode == 2:
            # user mode
            for s in symops:
                for a in self.atoms:
                    if a.GetVector() == s.apply(a.GetVector()):
                        continue
                    atom = ob.OBAtom()
                    atom.SetVector(self.f2c(s.apply(self.c2f(a.GetVector()))))
                    if a.GetVector().IsApprox(atom.GetVector(),0.01):
                        continue
                    atom.SetAtomicNum(a.GetAtomicNum())
                    self.OBMol.AddAtom(atom)
                self.OBMol.ConnectTheDots()
            monsx = self.OBMol.Separate()
            if molmap:
                res = list()
                for i in molmap:
                    m = ob.OBMol()
                    for j in i:
                        for k in xrange(monsx[j].NumAtoms()):
                            m.AddAtom(monsx[j].GetAtom(k + 1))
                    m.ConnectTheDots()
                    res.append(m)
                self.set_mon_map(res)
            else:
                self.set_mon_map(monsx)


    def iter_eqiv(self, v):
        for symop in self.symops:
            vs = symop.apply(v)
            yield vs, symop

    def iter_close(self, molid=0, vdwinc=1.0):
        s0 = ob.vector3()
        for i in self.mol_map[molid]: s0 += self.vf[i]
        slist = list()
        slist.append(s0)
        for s in self.iter_symop_t2():
            _found = []
            I = s.is_identity()
            for a0idx in self.mol_map[molid]:
                for molid1 in xrange(len(self.mol_map)):
                    if molid1 in _found: continue
                    if I and molid == molid1: continue
                    for a1idx in self.mol_map[molid1]:
                        vt = self.f2c(s.apply(self.vf[a1idx]))
                        svdw2 = self.vdw[a0idx] + self.vdw[a1idx] + vdwinc
                        svdw2 *= svdw2
                        if self.vc[a0idx].distSq(vt) < svdw2:
                            s1 = ob.vector3(0, 0, 0)
                            for i in self.mol_map[molid1]: s1 += s.apply(self.vf[i])
                            _found1 = False
                            for i in slist:
                                if i.distSq(s1) < 0.001:
                                    _found1 = True
                                    break
                            if _found1: break
                            _found.append(molid1)
                            slist.append(s1)
                            yield s, molid1
                            break

