#!/usr/bin/env python

from __future__ import division
from cMol import *

def calc_c(vecs):
    c=pybel.ob.vector3(0,0,0)
    l=len(vecs)
    for v in vecs:
        c+=v
    return pybel.ob.vector3(c.GetX()/l,c.GetY()/l,c.GetZ()/l)

def resfile_atom(name, ntype, x, y, z):
    return '%-4s %1i %7.4f %7.4f %7.4f 1.0 0.05\n' %(name,ntype,x,y,z)


if __name__=='__main__':
    import sys,os
    import pybel

    basename=sys.argv[1]

    if not os.path.isfile(basename+'.symm'):
        sys.stdout.write('No symm file!')
        sys.exit()
    if not os.path.isfile(basename+'.ene'):
        sys.stdout.write('No ene file!')
        sys.exit()

    symm=[i.split() for i in file(basename+'.symm').readlines()]
    ene=[float(i) for i in file(basename+'.ene').readlines()]
    for i in range(len(ene)): 
      if ene[i]>0: ene[i]=-0.05

    cmol=read_cmol(basename)
    if not cmol: sys.exit()

    newmol=pybel.Molecule(pybel.ob.OBMol())
    centers=list()
    for m in cmol.mol_map:
        centers.append(calc_c([cmol.c2f(cmol.atoms[a].GetVector()) for a in m]))
    en_sorted = [i for i in ene]
    en_sorted.sort()
    en_scale=1/(2*en_sorted[0])

    f=file(basename+'.res','w')
    #a) coordinates
    iat=0
    for i in centers:
        f.write(resfile_atom('M%i' %iat, 2, i.GetX(), i.GetY(), i.GetZ()))
        at=pybel.ob.OBAtom()
        at.SetAtomicNum(6)
        at.SetVector(pybel.ob.vector3(cmol.f2c(i)))
        newmol.OBMol.AddAtom(at)
        iat+=1
    for i in symm:
        i[1]=int(i[1])
        i[2]=int(i[2])
        x=centers[i[2]].GetX()
        y=centers[i[2]].GetY()
        z=centers[i[2]].GetZ()
        c=eval(i[0])
        #print i,(x,y,z),c
        #print SymOp(i[0])
        #v=SymOp(i[0]).apply(centers[i[2]])
        v=ob.vector3(c[0],c[1],c[2])
        #scale
        s=en_scale*ene[symm.index(i)]
        v-=centers[i[1]]
        v*=s
        v+=centers[i[1]]
        f.write(resfile_atom('M%i' %iat, 2, v.GetX(), v.GetY(), v.GetZ()))
        at=pybel.ob.OBAtom()
        at.SetAtomicNum(1)
        at.SetVector(pybel.ob.vector3(cmol.f2c(v)))
        newmol.OBMol.AddAtom(at)
        iat+=1
    
    f.write('END\n')
    
    #xp commands
    f.write('FMOL\n')
    f.write('UNDO\n')
    
    for i in symm:
        f.write('link 1 m%i m%i\n' %(int(i[1]), symm.index(i)+len(cmol.mol_map)))
        bo=pybel.ob.OBBond()
        bo.Set(0,newmol.OBMol.GetAtom(i[1]+1),newmol.OBMol.GetAtom(symm.index(i)+len(cmol.mol_map)+1),1,0)
        newmol.OBMol.AddBond(bo)
        
    f.write('END\n')
    
    output = pybel.Outputfile("mol2", basename+".mol2", overwrite=True)
    output.write(newmol)
    output.close()

