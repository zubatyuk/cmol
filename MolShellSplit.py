#!/usr/bin/env python

from cMol import *

if __name__ == '__main__':
    import sys,os
    import pybel
    basename=sys.argv[1]

    cmol=read_cmol(basename)
    if not cmol: sys.exit()
    
    convs=ob.OBConversion()
    convs.SetOutFormat('xyz')
    conv=ob.OBConversion()
    conv.SetOutFormat('xyz')

    symm=file(basename+'.symm','w')

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
    symm.close()
            
