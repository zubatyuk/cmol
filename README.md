#CMOL
##### *Collection of Python scripts for Energy Vector Digram analysis of crystal structures.*

###Introduction

Energy Vector Diagram method provides _ab initio_ methodology for analysis of inter-molecular interactions in molecular crystals. 

For details see:

  [Zubatyuk, R. I.; Sinelshchikova, A. A.; Enakieva, Y. Y.; Gorbunova, Y. G.; Tsivadze, A. Y.; Nefedov, S. E.; Bessmertnykh-Lemeune, A.; Guilard, R.; Shishkin, O. V. Insights into the crystal packing of phosphorylporphyrins based on the topology of their intermolecular interaction energies. CrystEngComm 2014, 16, 10428–10438.](http://xlink.rsc.org/?DOI=C4CE01623H)

  [Shishkin, O. V; Zubatyuk, R. I.; Shishkina, S. V; Dyakonenko, V. V; Medviediev, V. V. Role of supramolecular synthons in the formation of the supramolecular architecture of molecular crystals revisited from an energetic viewpoint. Phys. Chem. Chem. Phys. 2014, 16, 6773–6786.](http://dx.doi.org/10.1039/c3cp55390f)

  [Yufit, D. S.; Shishkin, O. V.; Zubatyuk, R. I.; Howard, J. A. K. Trimethyltrioxane (Paraldehyde) and Its Halomethanes Complexes: Crystallization, Structures, and Analysis of Packing Motifs. Cryst. Growth Des. 2014, 14, 4303–4309.](http://dx.doi.org/10.1021/cg500354t)

  [Shishkin, O. V.; Medvediev, V. V.; Zubatyuk, R. I. Supramolecular architecture of molecular crystals possessing shearing mechanical properties: columns versus layers. CrystEngComm 2013, 15, 160–167.](http://xlink.rsc.org/?DOI=c2ce26126j)

  [Yufit, D. S.; Zubatyuk, R.; Shishkin, O. V.; Howard, J. a. K. Low-melting molecular complexes. Halogen bonds in molecular complexes of bromoform. CrystEngComm 2012, 14, 8222.](http://xlink.rsc.org/?DOI=c2ce26191j)

  [Shishkin, O. V.; Medvediev, V. V.; Zubatyuk, R. I.; Shyshkina, O. O.; Kovalenko, N. V.; Volovenko, J. M. Role of different molecular fragments in formation of the supramolecular architecture of the crystal of 1,1-dioxo-tetrahydro-1λ6-thiopyran-3-one. CrystEngComm 2012, 14, 8698.](http://xlink.rsc.org/?DOI=c2ce26332g)

EVD of crystal structure could be constructed starting from proper  [CIF](http://en.wikipedia.org/wiki/Crystallographic_Information_File) file with now disorder or similar shortcomingsin following steps:

  1. Normalize of positions of hydrogen atoms

  2. For each symmetry unique molecule, build first coordination sphere

  3. For each molecule from first coordination sphere, calculate interaction energy with base molecule
  
  4. In the crystal structure, build set of vectors starting from geometrical center of base molecule, pointing to centers of neighboring molecules, with lenght Rij = Eij * Dij / 2Emax
  
  5. Visualize resulted diagram.

###Dependencies

 _Download links for Windows_

  Python : [http://python.org/ftp/python/2.7.5/python-2.7.5.msi]
  
  Openbabel : [http://sourceforge.net/projects/openbabel/files/openbabel/2.3.2/OpenBabel2.3.2a_Windows_Installer.exe/download]
  
  Openbabel Python module : [http://sourceforge.net/projects/openbabel/files/openbabel-python/1.8/openbabel-python-1.8.py27.exe/download]
  
###Typical usage

Start with `MyStructure.cif`

####Calculate coordination spheres

  `MolShellSplit.py MyStructure`

  The script will read CIF file, normalize positions of hydrogen atoms (C-H: 1.089, N-H: 1.015, O-H: 0.993), find and for each symmetry unique molecule will calculate its first coordination sphere. Will output

  `MyStructure.symm` -  each line contains symmemetry operation and IDs of base and symmetry related molecule for all molecules from first coordination sphere (defined as all molecule which have contacs with base molecule within sum of VdW radii + 1A).

  `MyStructure-s##.xyz` - geometry of first coordination sphere for each symmetry unique molecule

  `MyStructure-d###.xyz` - geometry of dimers (base molecule + neighbiring molecule)

####Perform quantum-chemical calculations of interaction energies

  You can prepare input files for DFT+BSSE calculations with (ORCA)[http://cec.mpg.de/forum/] using mkinps-orca-dim.py script. 

    `mkinps-orca-dim.py -b -r '!b97-d3 tzvp tzvp/j nososcf' MyStructure`

Run calculations. You can extract interaction energies with geten_orca.py bash script. This will generate list on interaction energies in MyStructure.ene file.

  `geten_orca.py -b MyStructure`
   
####Build EVDs

  `urchins.py MyStructure` will read MySTructure.cif MyStructure.symm and MyStructure.ene files and output PDB file suitable for viewing with Mercury. 
  
  
