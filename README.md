# alterPDBbfactor

This R package attempts to obtain all or a specific PDB structure associated to a UniProt AC of your interest and replaces their B-factor column with scores of your preference in the corresponding PDB chain

## Installation

### Dependencies

Before installing this package, ensure you have the following: 
  - By now this package has only been tested in Unix based OS
  - [R](https://cran.r-project.org/) version >= 4.0.X. Full list of R packages is in the DESCRIPTION file ('Imports') 

### How to install

From R, run any of the following sets of commands:
	
a) If installing from GitHub:

 ``install.packages("devtools")``


 ``devtools::install_github("marciaah/alterPDBbfactor")``


b) If installing from a gzipped copy of the package:

``install.packages("../alterPDBbfactor_1.0.0.tar.gz", repos = NULL)``

Then just load the package and it will be ready to use:

``library(alterPDBbfactor)``

## How to use

To run an example, use this command:

``alterPDBbfactor(pdbid = "all",uniprotAcc = "P51587",seq_scores = "./test_input/P51587_random_scores.tsv",outfolder = "./out/")``

Where:

``pdbid`` = a PDB identifier for only one PDB or the word ``"all"`` for all the available PDBs that correspond to the UniProtAC

``uniprotAcc`` = the uniProtKB accession code of your protein 

``seq_scores`` = Path to a file where you have a table with two columns: 


   ``aacPos`` : residue position number  

   ``score`` : by-residue scores that you want to add in the B-factor field of the PDB

   Use this header names otherwise it won't work! Check the ``./test_input/P51587_random_scores.tsv`` file for reference.


``outfolder`` = Path to a folder where the altered PDB structures will be saved. 



All parameters are mandatory!

