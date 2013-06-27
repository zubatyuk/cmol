#!/usr/bin/env python


'''
Created on May 6, 2011

@author: Roman
'''

import pybel, sys, os, copy

def mkinp_orca(mols, bq=()):
    etab = pybel.ob.OBElementTable()
    res=list()
    res.append('!blyp def2-tzvp def2-tzvp/j ecp{def2-tzvp,def2-tzvp/j} vdw10 nososcf\n')
    res.append('*xyz 0 1\n')
    for i in range(len(mols)):
        if i in bq: bqflag=':'
        else: bqflag=' '
        for j in range(len(mols[i].atoms)):
            atom=mols[i].atoms[j].OBAtom
            res.append('%2s %s % 10.6f % 10.6f % 10.6f\n' %(etab.GetSymbol(atom.GetAtomicNum()),
                                    bqflag, atom.GetX(), atom.GetY(), atom.GetZ()))
    res.append('*\n')
    return res

def mkxyz(mols):
    etab = pybel.ob.OBElementTable()
    res=list()
    for i in range(len(mols)):
        for j in range(len(mols[i].atoms)):
            atom=mols[i].atoms[j].OBAtom
            res.append('%2s % 10.6f % 10.6f % 10.6f\n' %(etab.GetSymbol(atom.GetAtomicNum()),
                                    atom.GetX(), atom.GetY(), atom.GetZ()))
    l=len(res)
    res.insert(0,'%i\n\n' %l)
    return res
            
if __name__ == '__main__':
    xyzfile=sys.argv[1] #should contain separated mols
    basename=os.path.splitext(xyzfile)[0]
    
    mols=[m for m in pybel.readfile('xyz', xyzfile)]


    #dimers
    f=file(basename+'_m0.inp', 'w')
    f.writelines(mkinp_orca(mols,(1,)))
    f.close()
    f=file(basename+'_m1.inp', 'w')
    f.writelines(mkinp_orca(mols,(0,)))
    f.close()
    f=file(basename+'_m0-1.inp', 'w')
    f.writelines(mkinp_orca(mols))
    f.close()

    sys.exit()
    #complex eint
    f=file('compl.inp', 'w')
    f.writelines(mkinp_orca(mols))
    f.close()
    for i in range(len(mols)):
        mm=copy.copy(mols)
        mm.pop(i)
        f=file('compl-%i.inp' %i, 'w')
        f.writelines(mkinp_orca(mm, (i,)))
        f.close()
        ii=range(len(mols))
        ii.pop(i)
        f=file('compl-%i_.inp' %i, 'w')
        f.writelines(mkinp_orca(mm, ii))
        f.close()
        #
        if i>0:continue
        #
        for j in range(i+1,len(mols)):
            mmm=copy.copy(mm)
            mmm.pop(j-1)
            f=file('compl-%i-%i.inp' %(i,j), 'w')
            f.writelines(mkinp_orca(mmm, (i,j)))
            f.close()
            iii=copy.copy(ii)
            iii.pop(j-1)
            f=file('compl-%i-%i_.inp' %(i,j), 'w')
            f.writelines(mkinp_orca(mm, iii))
            f.close()

