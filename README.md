#CMOL
##### *Collection of Python scripts for Energy Vector Digram analysis of crystal structures.*

###Introduction

Energy Vector Diagram method provides _ab initio_ methodology for analysis of inter-molecular interactions in molecular crystals. 

For details see:

  [Dyakonenko, Viktoriya V., Andrey V. Maleev, Alexander I. Zbruyev, Valentin A. Chebanov, Sergey M. Desenko, and Oleg V. Shishkin. 2010. “Layered Crystal Structure of Bicyclic Aziridines as Revealed by Analysis of Intermolecular Interactions Energy.” _CrystEngComm_ 12 (6): 1816. doi:10.1039/b922131j.](http://xlink.rsc.org/?DOI=b922131j)

  [Shishkin, Oleg V., Volodymyr V. Medvediev, Roman I. Zubatyuk, Olena O. Shyshkina, Nataliya V. Kovalenko, and Julian M. Volovenko. 2012. “Role of Different Molecular Fragments in Formation of the Supramolecular Architecture of the Crystal of 1,1-Dioxo-Tetrahydro-1λ6-Thiopyran-3-One.” _CrystEngComm_ 14 (24): 8698. doi:10.1039/c2ce26332g.](http://xlink.rsc.org/?DOI=c2ce26332g)

  [Yufit, Dmitry S., Roman Zubatyuk, Oleg V. Shishkin, and Judith a. K. Howard. 2012. “Low-Melting Molecular Complexes. Halogen Bonds in Molecular Complexes of Bromoform.” _CrystEngComm_ 14 (23): 8222. doi:10.1039/c2ce26191j.](http://xlink.rsc.org/?DOI=c2ce26191j)

  [Shishkin, Oleg V., Volodymyr V. Medvediev, and Roman I. Zubatyuk. 2013. “Supramolecular Architecture of Molecular Crystals Possessing Shearing Mechanical Properties: Columns Versus Layers.” _CrystEngComm_ 15 (1): 160. doi:10.1039/c2ce26126j.](http://xlink.rsc.org/?DOI=c2ce26126j)

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
  
  SymPy library : [https://code.google.com/p/sympy/downloads/detail?name=sympy-0.7.2.win32.exe]

###Typical usage

Start with `MyStructure.cif`

####Calculate coordination spheres

  `MolShellSplit.py MyStructure`

  The script will read CIF file, normalize positions of hydrogen atoms (C-H: 1.089, N-H: 1.015, O-H: 0.993), find and for each symmetry unique molecule will calculate its first coordination sphere. Will output

  `MyStructure.symm` -  each line contains symmemetry operation and IDs of base and symmetry related molecule for all molecules from first coordination sphere (defined as all molecule which have contacs with base molecule within sum of VdW radii + 1A).

  `MyStructure-s##.xyz` - geometry of first coordination sphere for each symmetry unique molecule

  `MyStructure-d###.xyz` - geometry of dimers (base molecule + neighbiring molecule)

####Perform quantum-chemical calculations of interaction energies

  You can prepare input files for DFT+BSSE calculations with (ORCA)[http://cec.mpg.de/forum/] using mkinps-orca-dim.py script (look at the script if you want to change level of theory). 

  For each `MyStructure-d###.xyz` run:

  `mkinps-orca-dim.py MyStructure-d001.xyz`

  More conveniently, using batch command in Windows:
  
  `for %f in (MyStructure-d???.xyz); do mkinps-orca-dim.py %f`

  or in bash:

  `for f in MyStructure-d???.xyz; do mkinps-orca-dim.py $f; done`

####Run calculations. You can extract interaction energies with geten-bsse.sh bash script 

  `geten-bsse.sh MyStructure |awk '{print $2}' > MyStructure.ene`
   
####Build EVDs

  `urchins.py -m MyStructure` will read MySTructure.cif MyStructure.symm and MyStructure.ene files and output MOL2 file suitable for viewing with Mercury. 
  
  Did not found any better way to put unit cell and symmetry operations in MOL2 file except open CIF in Mercury, save to a temporary MOL2 file and just copy two last lines from that file.


  
