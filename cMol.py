#!/usr/bin/env python

from __future__ import division
import openbabel as ob
import sympy, types, re, fractions,os, sys
import libcif
import pybel

def get_symops(cifname):
    symops=list()
    cif=libcif.parsecif(libcif.readcif(cifname))
    for l in cif['loop_']:
        for t in '_space_group_symop_operation_xyz','_symmetry_equiv_pos_as_xyz':
            if t in l[0]:
                for s in l[1]:
                    symops.append(SymOp(s[t]))
            continue
    return symops

def read_cmol(basename):
    if not os.path.isfile(basename+'.cif'):
        sys.stderr.write('No CIF to work with!\n')
        return None
    if os.path.isfile(basename+'.map'):
        '''
        file should look like:
        2
        (SymOp('1/2-x,1-y,1/2+z'),)
        ((0,1),)
        '''
        separate_args=[eval(l.rstrip().rstrip()) for l in file(basename+'.map').readlines() if l]
    else:
        separate_args=()
    symops=get_symops(basename+'.cif')
    mol=pybel.readfile('cif',basename+'.cif').next()
    return cMol(mol.OBMol,mol.unitcell,symops,separate_args=separate_args)

def printvec(v):
    return str((v.GetX(),v.GetY(),v.GetZ()))

class SymOp:
    def __init__(self,r=ob.transform3d(1),t=ob.vector3(0,0,0)):
        self.T=ob.vector3(0,0,0)
        if type(r) is types.StringType:
            self.set_from_str(r)
        else:
            self.R=r
        self.add_trans(t)
        
    def __eq__(self,other):
        assert type(other)==type(self)
        v1=self.apply(ob.vector3(.1,.2,.3))
        v2=other.apply(ob.vector3(.1,.2,.3))
        if v1.IsApprox(v2,0.01):
            return True
        else:
            return False
        
    def __repr__(self):
        re1=re.compile('([+-]?[xyz])+')
        re2=re.compile('([-+]?(\d+(\.\d*)?|\.\d+))')
        s=self.R.DescribeAsString().split(',')
        t=(self.T.GetX(),self.T.GetY(),self.T.GetZ())
        (x,y,z)=t
        st=[eval(i) for i in s]
        r=list()
        for ii in [str(t[i]+sympy.expand(s[i]+'+1.0-1.0')) for i in (0,1,2)]:
            #print '>',ii
            ii=ii.replace(' ','')
            m2=re2.search(ii)
            #print '>>', re1.search(ii).group(0)
            if m2:
                f=fractions.Fraction.from_float(float(m2.group(0))).limit_denominator(6)
                if f.denominator==1:
                    ii=str(f.numerator) +re1.search(ii).group(0)
                else:
                    ii='%i/%i' %(f.numerator,f.denominator) +re1.search(ii).group(0)
            ii=re.sub('(\d)([xyz])',lambda x: x.group(1)+'+'+x.group(2),ii)
            ii=re.sub('^\+','',ii)
            r.append(ii)
        return ','.join(r)
            
    def is_identity(self):
        v0=ob.vector3(.1,.2,.3)
        v=self.apply(v0)
        if v.IsApprox(v0,0.01):
            return True
        return False

    def set_from_str(self,s):
        assert type(s) is types.StringType, 's is not string: %s' %s
        assert 'x' in s and 'y' in s and 'z' in s, 's sould containt x,y,z: %s' %s
        c1=((1,0,0),(0,1,0),(0,0,1))
        v1=[eval(s) for (x,y,z) in c1]
        (x,y,z)=(0,0,0)
        v0=eval(s)
        v=[(v1[0][i]-v0[i],v1[1][i]-v0[i],v1[2][i]-v0[i]) for i in xrange(3)]
        vi=[int(i) for i in v0]
        vt=[i-int(i) for i in v0]
        vec=[ob.vector3(i,j,k) for (i,j,k) in v]
        self.R=ob.transform3d(vec[0],vec[1],vec[2],ob.vector3(vt[0],vt[1],vt[2]))
        self.T=ob.vector3(vi[0],vi[1],vi[2])
        
    def add_trans(self,t):
        '''TODO: sanity check
        '''
        self.T+=t
    
    def set_trans(self,t):
        self.T=t
        
    def apply(self,v):
        res=self.R*v
        res+=self.T
        return res

class cMol(object):
    def __init__(self,OBMol,unitcell,symops=list(),separate_args=(1,)):
        self.OBMol=OBMol
        self.unitcell=unitcell
        self.f2c=self.unitcell.FractionalToCartesian
        self.c2f=self.unitcell.CartesianToFractional
        self.symops=symops
        
        self.norm_h()
        self.separate(separate_args)
        #self.move_to_center()
        
        self.vc=[a.GetVector() for a in self.atoms]
        self.vf=[self.c2f(v) for v in self.vc]
        obet=ob.OBElementTable()
        self.vdw=[obet.GetVdwRad(a.GetAtomicNum()) for a in self.atoms]
        
    def iter_trans2(self):
        d=(0,1,-1,2,-2,3,-3)
        for i in d:
            for j in d:
                for k in d:
                    yield ob.vector3(i,j,k)
                
    def iter_symop_t2(self):
        for s in self.symops:
            for t in self.iter_trans2():
                yield SymOp(s.R,t) 
        
    @property    
    def center_c(self):
        c=ob.vector3(0,0,0)
        n=0
        for atom in self.atoms:
            c+=atom.GetVector()
            n+=1
        return ob.vector3(c.GetX()/n,c.GetY()/n,c.GetZ()/n)
    
    @property
    def center_f(self):
        return self.c2f(self.center_c)
    
    def move_to_center(self):
       c=self.center_f
       vc=list()
       for i in (c.GetX(), c.GetY(), c.GetZ()):
           if i < 0: i-=1
           vc.append(int(i))
       t=ob.vector3(vc[0],vc[1],vc[2])
       if not t.IsApprox(ob.vector3(0,0,0),0.001):
            print "WARNING: moving!"
            for atom in self.atoms:
                v=atom.GetVector()
                v-=t
                atom.SetVector(v)
    @property
    def atoms(self):
        return [self.OBMol.GetAtom(i+1) for i in range(self.OBMol.NumAtoms())]
    
    def norm_h(self):
        _bond_data={6: 1.089, 7: 1.015, 8: 0.993}
        for bond in ob.OBMolBondIter(self.OBMol):
            atoms = (bond.GetBeginAtom(), bond.GetEndAtom())
            anums = (atoms[0].GetAtomicNum(), atoms[1].GetAtomicNum())
            if anums[0] == 1: bond.SetLength(atoms[1], _bond_data[anums[1]])
            if anums[1] == 1: bond.SetLength(atoms[0], _bond_data[anums[0]])

        
    def set_mon_map(self,mons=[]):
        self.mol_map=[list() for i in mons]
        for mon in mons:
            for i in range(mon.NumAtoms()):
                v1=mon.GetAtom(i+1).GetVector()
                for j in range(self.OBMol.NumAtoms()):
                    v2=self.OBMol.GetAtom(j+1).GetVector()
                    if v1.IsApprox(v2,0.001):
                        self.mol_map[mons.index(mon)].append(j)
                        continue
                    
    def separate(self,args):
        try:        mode=args[0]
        except:     mode=1
        try:        symops=args[1]
        except:     symops=list()
        try:        molmap=args[2]
        except:     molmap=list()

        
        if mode==0:
            #no separation
            self.set_mon_map((self.OBMol,))
        if mode==1:
            #OBMol.Separate
            self.set_mon_map(self.OBMol.Separate())
        if mode==2:
            #user mode
            for s in symops:
                for a in self.atoms:
                    atom=ob.OBAtom()
                    atom.SetVector(self.f2c(s.apply(self.c2f(a.GetVector()))))
                    atom.SetAtomicNum(a.GetAtomicNum())
                    self.OBMol.AddAtom(atom)
                self.OBMol.ConnectTheDots()
            monsx=self.OBMol.Separate()
            if molmap:
                res=list()
                for i in molmap:
                    m=ob.OBMol()
                    for j in i:
                        for k in xrange(monsx[j].NumAtoms()):
                            m.AddAtom(monsx[j].GetAtom(k+1))
                    m.ConnectTheDots()
                    res.append(m)
                self.set_mon_map(res)
            else:
                self.set_mon_map(monsx)
                
                    
    def iter_eqiv(self,v):
        for symop in self.symops:
            vs=symop.apply(v)
            yield vs,symop    

    def iter_close(self,molid=0,vdwinc=1.0):
        s0=ob.vector3()
        for i in self.mol_map[molid]: s0+=self.vf[i]
        slist=list()
        slist.append(s0)
        for s in self.iter_symop_t2():
            _found=[]
            I=s.is_identity()
            for a0idx in self.mol_map[molid]:
                for molid1 in xrange(len(self.mol_map)):
                    if molid1 in _found: continue
                    if I and molid==molid1: continue
                    for a1idx in self.mol_map[molid1]:
                        vt=self.f2c(s.apply(self.vf[a1idx]))
                        svdw2=self.vdw[a0idx]+self.vdw[a1idx]+vdwinc
                        svdw2*=svdw2
                        if self.vc[a0idx].distSq(vt) < svdw2:
                            s1=ob.vector3(0,0,0)
                            for i in self.mol_map[molid1]: s1+=s.apply(self.vf[i])
                            _found1=False
                            for i in slist:
                                if i.distSq(s1)< 0.001:
                                    _found1=True
                                    break
                            if _found1: break    
                            _found.append(molid1)
                            slist.append(s1)
                            yield s,molid1
                            break
                         
    
if __name__ == '__main__':
    import sys,os
    import pybel
    import libcif
    struct=sys.argv[1]
    basename=os.path.splitext(struct)[0]

    #separate_args=()
    separate_args=(2, (SymOp('1/2-x,1-y,1/2+z'),), ((0,1),) )

    #parse cif
    symops=list()
    cif=libcif.parsecif(libcif.readcif(struct))
    for l in cif['loop_']:
        for t in '_space_group_symop_operation_xyz','_symmetry_equiv_pos_as_xyz':
            if t in l[0]:
                for s in l[1]:
                    symops.append(s[t])
            continue 

    mol=pybel.readfile(os.path.splitext(struct)[1][1:],struct).next()
    cmol=cMol(mol.OBMol,mol.unitcell,symops,separate_args=separate_args)
    
    convs=ob.OBConversion()
    convs.SetOutFormat('xyz')
    conv=ob.OBConversion()
    conv.SetOutFormat('xyz')

    id=0
    for id0 in xrange(len(cmol.mol_map)):
        m0=ob.OBMol()
        for i in cmol.mol_map[id0]:
            m0.AddAtom(cmol.atoms[i])
        ms=list()
        for s,id1 in cmol.iter_close(id0):
            print s,id0,id1 
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
            
