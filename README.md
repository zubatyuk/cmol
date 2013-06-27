# CMol


Collection of Python scripts and libraries for performing packing analisys of molecular clsystals.

#### Dependencies 

Download links are for Windows. Install in listed order.

* Python 2.7 http://python.org/ftp/python/2.7.5/python-2.7.5.msi
* Openbabel http://sourceforge.net/projects/openbabel/files/openbabel/2.3.2/OpenBabel2.3.2a_Windows_Installer.exe/download
* Openbabel Python module http://sourceforge.net/projects/openbabel/files/openbabel-python/1.8/openbabel-python-1.8.py27.exe/download
* SymPy library https://code.google.com/p/sympy/downloads/detail?name=sympy-0.7.2.win32.exe

#### How to run the code

* Clone GitHub repository from https://github.com/zubatyuk/cmol
* Put MyStructure.cif in this directory
* From command line execute
<pre><code>MolShellSplit.py MyStructure
</code></pre>
This will produce some files:

  * MyStructure.symm - each line contains symmemetry operation and IDs of base and symmetry related molecule for all molecules from first coordination sphere (defined as all molecule which have contacs with base molecule within sum of VdW radii + 1A). 
  * MyStructure-s??.xyz - XYZ file(s) contatining base molecule and its coordination sphere.
  * MyStructure-d???.xyz - XYZ files for dimers containing base and symmetry related molecules.
  
* For each MyStructure-d???.xyz run:
<pre><code>mkinps-orca-dim.py MyStructure-d001.xyz
</code></pre>
Look at line 44 of the script to change computational method.

The script will produce impuut files for BSSE corrected ORCA calculations. 

* Collect computed interaction energies in file MyStructure.ene
* Build energy vector diagrams of the structure
<pre><code>urchins.py MyStructure
</code></pre>
This will produce MyStructure.res and MyStructure.mol2 files (SHELX and MOL2 formats) almost ready for visualization. You need to add unit cell and symmetry info to these files. Simpliest way: use Mercury to produce res and mol2 files from your cif.

Use XP to visualize res file and Mercury for mol2.






