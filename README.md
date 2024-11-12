PRA-MutPred

This server predicts protein-RNA binding affinity change upon mutation using sequence and structure-based features of protein-RNA complexes. We found that number of atomic contacts between the RNA and protein,contact potentials and accessibility of the residues are important to understand the binding affinity. PRA-MutPred shows a correlation of 0.85 and a mean absolute error(MAE) of 0.59 in Self-consistency, correlation of 0.75 and a MAE of 0.84 kcal/mol in jack-knife test.

The reuirements for the instalation
Install packages such sklearn, 	Bio, and all the neccessary packages
Install other softwares such as 
Foldx software. (Source: https://foldxsuite.crg.eu/) Execute x3dna_setup and add the path to environment variable (http://forum.x3dna.org/downloads/3dna-download/). Request and download the 3DNA. 
Stride(https://webclu.bio.wzw.tum.de/stride/)
PyMol(https://www.pymol.org/). Make sure "pymol" executable is inside Pymol/pymol/. 
PsiBlast, UniProtDB (https://www.ncbi.nlm.nih.gov/books/NBK52640/)
Run "python3 pdpredict.py" Enter option for PDB-ID and PDB-file Provide the PDB ID or PDB-file Please the mutations in the following format 'chain:wt-position-mutation'.(For eg. A:F51A For mutations has to be made in multiple chains of homomer, use as A-B:H49D)
The program will be executed, and the result will be provided in result.txt In the Form of PDB ID (Mutation), Predicted affinity changes(kcal/mol), effect of mutation. You can also access it through the webserver https://web.iitm.ac.in/bioinfo2/pramutpred/index.html.
