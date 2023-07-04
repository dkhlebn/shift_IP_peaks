# shift_IP_peaks
This is the repo for scripts used to estimate significance of found RNA-protein-DNA triple interactions 
in [insert citation].


## SHIFT_AUX_FILES directory
Directory with auxilliary files for shifting RNA-Protein and DNA-protein peaks.

  * `ChromHMM_*` - simplified ChromHMM annotation files for K562 [1] and mESC [2] cell lines respectively;

  * `K562_SPIN_25kb_hg38.bed` - SPIN annotation files for K562 cell line [3];

  * `*.compartments.bed` - files with A/B-compartments annotation for K562 and mESC;

  * `close_proteins.txt` - text file containing dictionary of close protein for each protein studied. Here 
     protein closeness is defined as size of RNA set that have a significant peak with both proteins for a 
     pair of them;

  * `func_proteins.txt` - text file containing dictionary of functionally close protein for each protein 
     studied. Here protein closeness is defined by their functional group inferred from various sources of 
     biological data;

  * `*_crd_map.txt` - files with absolute to relative mapping in terms of 
     A/B-compartment structure for K562 and mESC.

## shiftlib directory
Directory with main driver scripts for shifting and intersecting data.

### `mapping.py` 
Defines a function `GENERATE_ABSREL_MAPPING`, taking **cell line** 
as an input and writing dict as string into file: `SHIFT_AUX_FILES/*_crd_map.txt`;
Generates a utility file to map set of genomic coordinates to a set of coordinates
relative to A/B compartments. 

### `annotTAD.py`
Defines a function `INTERSECT_W_TAD`, taking **ChIP_PATH** path to ChIP-Seq peaks file,
**COMPTS_PATH** path to compartments file (`SHIFT_AUX_FILES/*.compartments.bed`) and 
outputs ChIP-Seq peaks with A/B-compartments annotation into **OUTPATH** file;
Is used to obtain A/B-compartment annotated ChIP-Seq peaks.

### `calculate_close.py`
Defines two functions: CALC_CLOSENESS (takes **RNA_LIST_DIR** path to directory with lists of RNA 
for each protein and an optional **threshold** set to 2/3 by default) that generates a 
`close_proteins` file in aux files directory; CALC_FUNCTIONAL that manually loads `func_proteins` 
into the same directory.

### `ChIPShift.py` 
Defines a function `SHIFT_DNA_Protein`, that shifts ChIP-Seq peaks according to 
A/B-compartmental structure. Arguments: **PEAKFILE_PATH** path to ChIP-Seq peaks file with A/B 
annotation (as shown below), **OUTFILE_PATH** path to shifted ChIP-Seq peaks file with A/B annotation, 
**SHIFTS_DICT** mapping from chromosome karyotype to shift size for a given organism (accounted 
for in `main.py`), **CMP_MAP_PATH** path to file of absolute-relative coordinate mapping in 
A/B-compartments in a given organism;
Is used in simulation to get shifted ChIP-Seq peaks.

###  `RIPShift.py`
Defines a function `PERMUTE_RNA_Protein`, that permutes the set of RNA-protein 
interactions peaks as follows. First, a pool of RNA-protein peaks is created from RNA-protein peaks 
of currently simulated protein and its close ploteins as defined by either `func_proteins` or 
`close_proteins` file; another pool of peaks is created from RNA-protein interactions of unrelated
proteins. Then, a pool of RNA-protein interactions is drawn from two pools with ratio 4:1 as 
for protein and close proteins to unrelated proteins interaction sets, while preserving the 
RNA biotype content of initial RNA-protein interaction dataset. Arguments: **PEAKFILE_PATH** path to 
fRIP/RIP-Seq or eCLIP peaks file with peak annotation (as shown below), **OUTFILE_PATH** path to 
permuted RNA-protein peaks file, **ORGANISM** - K562/mESC for different processing of mouse and human
RNA-protein peaks, **FUNC_FLAG** - a flag signifying the use of `func_proteins` file instead of 
defaulting to `close_proteins` usage.

###  `construct_triads.py`
Defines a function that runs a bash script `` that constructs the RNA-protein-DNA triad from simulated
peaks files and real RNA-DNA data, writing the output to specified inner location. Arguments: **RD_DATA** 
path to RNA-DNA interactions data file, **RIP_TYPE** string representing the type of RNA-protein interactions
capture protocol, **OUTDIR** a dictionary with locations of simulated RIP and ChIP peak files, **PROT** name 
of the protein for which the triads are constructed, **ORG** string specifying the cell line, **TMP_dir** working
dir for the script, **fname** string specifying the filename of the permuted RNA-protein peaks file. Returns 
the filename where the simulated RNA-protein-DNA interactions are stored as well as the number of simulated interactions. 

###  `annot_simtriads.py`
Defines a function that obtains the ChromHMM (for K562 and mESC) and SPIN (K562 only) annotation of DNA parts
of simulated interaction triads. Arguments: **intersectBed** path to bedtools executable,  **out_tri_fname** - file where 
simulated RNA-protein-DNA interactions are stored, **ORG** string specifying the cell line, **TMP_dir** working dir for
the script. Returns somewhat parsable wtring for annotation state counts of the DNA parts of simulated triads.

### `intersect_script.sh`
A bash script to construct triads.

###  `shift_intersect.py`
Defines a driver function `SHIFT_THEN_INTERSECT`, that shifts and 
intersects RIP- and ChIP-Seq peaks for a given protein in a given cell_line. Arguments:
**PROT_ORG** tuple (protein, cell_line), **INIT_LOCS** initial locations for RIP- and ChIP-Seq 
peak files, **OUTDIR** output directory for shifted RIP- and ChIP-Seq files to be put into, 
**SHIFTS_DICT** - dictionary for shifting RIP- and ChIP-Seq data for a given organism,
**RD_DATA** - path to RNA-DNA interactions data, **AUX_DIR** directory where auxilliary files 
are situated), **ANNOT_FLAG** flag used to specify whether the assessment of annotation 
states is needed, **FUNCTIONAL_FLAG** flag used to specify the usage of `func_proteins` file 
instead of defaulting to `close_proteins` for RNA-protein peaks permutations. 
Returns a parsable string of triad file, triad number and, optionally, annotation stated for
DNA parts.

###  `main.py`
Main file, runs simulations. Takes arguments in argv style: **ChIP_DIR**, **RIP_DIR**, 
**RD_DIR**, **TMP_DIR**, **ANNOT_FLAG**, **FUNCTIONAL_FLAG**.

# Examples of data formats to be used:

Examples of some files are provided in SHIFT_AUX_FILES directory. Others are:
## 1. Un-annotated ChIP-Seq files:
|chrN|start|end|name|score|strand|something|pvalue|qvalue|something|
|----|-----|---|----|-----|------|---------|------|------|---------|
|chr1|42959695|42960055|.|538|.|12.12817|-1.00000|0.47207|180|
|chr1|22026852|22027212|.|597|.|12.15709|-1.00000|0.47387|180|
## 2. Annotated ChIP-Seq files:
```
|chrN|start|end|name|score|strand|something|pvalue|qvalue|something|cmpt|
|----|-----|---|----|-----|------|---------|------|------|---------|----|
|chr1|42959695|42960055|.|538|.|12.12817|-1.00000|0.47207|180|A|
|chr1|22026852|22027212|.|597|.|12.15709|-1.00000|0.47387|180|A|
```
## 3. RNA-protein peaks files:
```
|chrN|start|end|strand|gene_type|gene_name|pvalue|qvalue|
|----|-----|---|------|---------|---------|------|------|
|chr10|992400|992700|+|protein_coding|GTPBP4|0.00247773|0.0034409135027472526|
|chr10|995700|996300|+|protein_coding|GTPBP4|0.00421554|0.005076761488981537|
```
## 4. RNA list file:
```
1016_H1-Esc
1031_H1-Esc
1045_H1-Esc
104_GM12878
1086_HelaS3
...
```
## 5. RNA-DNA contacts
```
dna_chr	dna_start	dna_end		rna_chr	rna_start	rna_end		gene_type	gene_name	scaling_value	pvalue	qvalue
chr10	53865689	53865713	chr1	167256826	167256851	-	POU2F1	protein_coding	1.09306870787035	0.13471874050772295	0.2448020706262164
chr10	64710524	64710548	chr1	145993344	145993388	+	TXNIP	protein_coding	0.901432832277764	0.03145473554779507	0.20506626292109298

```

# Example usage (as of source article):
  You could notice the number of iterations in main script is set to 20. This is intentional since 
  the script was run on 50 instances simultaneously to obtain results faster.
```{python}
# generate A/B compartment mapping
python3 shiftlib/mapping.py SHIFT_AUX_FILES;
# annotate ChIP-Seq with A/B-compartments
python3 shiftlib/annotTAD.py ChIP-Seq_Peaks BG_WORKING DIR;
# calculate close proteins from the RNA they interact with
python3 shiftlib/calculate_close.py PROT_RNA_LISTS;

#run the simulation
python3 shiftlib/main.py ChIP_w_cmpts/WINDOW_0 RIP_DATA/annot_WINDOW_0 RNA_DNA_DATA SHIFT_WORKING_DIR_1 1 0
```
