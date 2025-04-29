# DockInspect Demo
This file demonstrates how to use the DockInspect app with sample data from `demo/test_files`.

Run the following to see available commands:

```bash
dockinspect
```
Expected output:
```
Usage: dockinspect [OPTIONS] COMMAND [ARGS]...

  Docking analysis tool

Options:
  --help  Show this message and exit.

Commands:
  ligand_info  Shows properties of a ligand from a SMILES string.
  pocket_info  Shows SASA, GRAVY and charge for predicted pockets.
  poses_info   Shows ligand–pocket interaction summary per docking pose.
  shell        Initializes a docking analysis session and launches the...
  visualize    Launches PyMOL visualization of a selected docking pose.
```

## ligand_info

```bash
dockinspect ligand_info "NC(=O)N"
```

```
Ligand Properties for SMILES: NC(=O)N
----------------------------------------
logP:        -0.98
SASA (A^2):  95.97
TPSA (A^2):  69.11
Volume (A^3): 173.03
Charge:      0
```

## pocket_info

```bash
dockinspect pocket_info \
demo/test_files/autodock_vina_results_2SRC_urea/structure.pdbqt \
demo/test_files/prankweb-2SRC/structure.cif_predictions.csv \
demo/test_files/prankweb-2SRC/structure.cif_residues.csv 
```

```
Pocket ID  SASA (A^2)      GRAVY     Charge
-------------------------------------------
pocket1       1907.37       0.59          0
pocket2        699.67      -2.51          0
pocket3         50.27      -2.75          2
pocket4        443.51      -1.15         -2
pocket5        239.40      -0.10          0
pocket6        757.28      -0.80          0
```

## poses_info
```bash
dockinspect poses_info \
"NC(=O)N" \
"2src" \
demo/test_files/autodock_vina_results_2SRC_urea/out_vina.pdbqt \
demo/test_files/autodock_vina_results_2SRC_urea/structure.pdbqt \
demo/test_files/prankweb-2SRC/structure.cif_predictions.csv \
demo/test_files/prankweb-2SRC/structure.cif_residues.csv 
```

```
Pose   Pocket     HBonds   GRAVY(P)/LogP(L) SASA(P/L/R)               Charge(P/L)    
-------------------------------------------------------------------------------------
1      pocket6    2        -0.80/-0.98     757.28/95.97/0.13         0/0            
2      pocket6    1        -0.80/-0.98     757.28/95.97/0.13         0/0            
3      pocket6    1        -0.80/-0.98     757.28/95.97/0.13         0/0            
4      pocket1    3        0.59/-0.98      1907.37/95.97/0.05        0/0            
5      pocket1    0        0.59/-0.98      1907.37/95.97/0.05        0/0            
6      pocket1    0        0.59/-0.98      1907.37/95.97/0.05        0/0            
7      pocket1    1        0.59/-0.98      1907.37/95.97/0.05        0/0            
8      pocket1    3        0.59/-0.98      1907.37/95.97/0.05        0/0            
9      pocket1    2        0.59/-0.98      1907.37/95.97/0.05        0/0  
```

## shell
Note: It is necessary to provide absolute path to the out_vina.pdbqt file.
```bash
dockinspect shell \
"NC(=O)N" \
"2src" \
<absolute path to dockinspect>/dockinspect/demo/test_files/autodock_vina_results_2SRC_urea/out_vina.pdbqt \
demo/test_files/autodock_vina_results_2SRC_urea/structure.pdbqt \
demo/test_files/prankweb-2SRC/structure.cif_predictions.csv \
demo/test_files/prankweb-2SRC/structure.cif_residues.csv 
```

```
Session initialized. Starting interactive mode...

Type 'help' or '?' to list commands.

> 
```

Demo workflow inside the shell:
```
> ligand_info
Ligand Properties for SMILES: NC(=O)N
----------------------------------------
logP:        -0.98
SASA (A^2):  95.81
TPSA (A^2):  69.11
Volume (A^3): 172.67
Charge:      0
> pocket_info --pocket_id pocket1
Pocket ID  SASA (A^2)      GRAVY     Charge
-------------------------------------------
pocket1       1907.37       0.59          0
> poses_info 1
Pose   Pocket     HBonds   GRAVY/LogP      SASA(P/L/R)               Charge(P/L)    
------------------------------------------------------------------------------------
1      pocket6    2        -0.80/-0.98     757.28/95.81/0.13         0/0            
> visualize
Launching PyMOL with mode 'hbonds'...
> visualize --save hbonds_session.pml
Saving PyMOL script...
PML script saved to hbonds_session.pml.
```