# shift_IP_peaks
A pack of utility scripts to shift DNA-Protein and RNA-Protein peaks along the genome


## SHIFT_AUX_FILES directory
Directory with auxilliary files for shifting RNA-Protein and DNA-protein peaks.

  * `*.compartments.bed` - files with A/B-compartments annotation for K562 and mESC;

  * `*_crd_map.txt` - files with absolute to relative mapping in terms of 
    A/B-compartment structure for K562 and mESC.

## shiftlib directory
Directory with main driver scripts for shifting and intersecting data.

  * `mapping.py` - defines a function `GENERATE_ABSREL_MAPPING`, taking **cell line** 
    as an input and writing dict as string into file: `SHIFT_AUX_FILES/*_crd_map.txt`;

  * `annotTAD.py` - defines a function `INTERSECT_W_TAD`, taking **ChIP_PATH** path to ChIP-Seq peaks file,
    **COMPTS_PATH** path to compartments file (`SHIFT_AUX_FILES/*.compartments.bed`) and outputs ChIP-Seq 
    peaks with A/B-compartments annotation into **OUTPATH** file;

  * `ChIPShift.py` - defines a function `SHIFT_DNA_Protein`, that shifts ChIP-Seq peaks according to 
    A/B-compartmental structure. Arguments: **PEAKFILE_PATH** path to ChIP-Seq peaks file with A/B 
    annotation, **OUTFILE_PATH** path to shifted ChIP-Seq peaks file with A/B annotation, 
    **SHIFTS_DICT** mapping from chromosome karyotype to shift size for a given organism (accounted 
    for in `main.py`), **CMP_MAP_PATH** path to file of absolute-relative coordinate mapping in 
    A/B-compartments in a given organism;

  * RIPShift.py - defines a function `SHIFT_RNA_Protein`, that shifts RIP-Seq peaks file. Arguments:
    **PEAKFILE_PATH** path to RIP-Seq peaks file, **OUTFILE_PATH** path to shifted RIP-Seq peaks file, 
    **SHIFT_DICT** mapping from chromosome karyotype to shift size for a given organism (accounted
    for in `main.py`), **CHR_LIST** karyotype for a cell_line;

  * shift_intersect.py - defines a driver function `SHIFT_THEN_INTERSECT`, that shifts (sic!) and 
    intersects RIP- and ChIP-Seq peaks for a given protein in a given cell_line. Arguments:
    **PROT_ORG** tuple (protein, cell_line), **INIT_LOCS** initial locations for RIP- and ChIP-Seq 
    peak files, **OUTDIR** output directory for shifted RIP- and ChIP-Seq files to be put into, 
    **SHIFTS_DICT** - dictionary for shifting RIP- and ChIP-Seq data for a given organism, 
    **CHR_LIST** keryotype for a cell_line, **RD_DATA** path to RNA-DNA interactions file, 
    **AUX_DIR** directory where auxilliary files are situated);

  * main.py
