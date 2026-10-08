# RBFE-Hybrid-Topology
Created to run a RBFE simulations with OpenFE 
## Pre-requisites 
- 1. Have an installation of micromamba and run:
``` bash
micromamba env create -n -f openfe.yml 
```
- 2. Make sure you have an SDF file with all your ligands and named.
- 3. Have a PDB file, that must be prepared.                  
## Run
- Open 01\_Create\_plan.py file and change the name of the files (sdf, pdb), set temperature, etc.
- Run as:
``` bash   
python 01_Create_plan.py
```

