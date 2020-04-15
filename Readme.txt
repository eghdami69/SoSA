SoSA is developed to produce, for each protein,  all possible conformation instances corresponding to each candidate topology and, from among them, find and report the best one. In other words, it reports the best conformation instance for each candidate topology.
=================================================================
Knowledge Engineering Research Group (KERG), Department of Computer Engineering,
Faculty of Engineering, Ferdowsi University of Mashhad, Iran. 
Home page: https://kerg.um.ac.ir/  
Contact: Dehgani.toktam@mail.um.ac.ir; naghibzadeh@um.ac.ir 
=================================================================

Input Files Format
=============
1- Proteins ' information (Proteins Folder)
First line: pdb-id
Second line: target protein sequence without gaps
Third line: Secondary structure of the target protein

2- Contact Scores (ultra_deep_learning_contactmap Folder, Raptor-X output files)
residue contact scores matrix

Output File Format (Predictions Folder)
==============
For each selected conformation instances corresponding to each candidate topology, two files are created:

1- Conformation File
For each beta-strand pairing, two lines are created:
First line: strand index of partner 1(starting from 0) strand index of partner 2(starting from 0)
Second line: pairing direction => P = parallel A = anti-parallel B = isolated beta-bridge

2- Contact-map File
b-residue contacts are reported in a matrix