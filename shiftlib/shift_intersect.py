import os
import sys
import math
import glob
import random
import logging
import pandas as pd
import shlex as sh
import subprocess as sp

from RIPShift import SHIFT_RNA_Protein
from ChIPShift import SHIFT_DNA_Protein
intersectBed="/home/tools/bedtools/bedtools-v2.16.2/bin/intersectBed"

def SHIFT_THEN_INTERSECT(PROT_ORG, INIT_LOCS, OUTDIR, SHIFTS_DICT, CHR_LIST, RD_DATA, AUX_DIR):
  ''' 
       PROT_ORG : (protein, cell_line), 
       INIT_LOCS : dict( RIP | ChIP => initial peak location),
       OUTDIR    : dict( RIP | ChIP => TMP_dir / RP | DP) - directory, NOT file,
       SHIFTS_DICT : dict( RIP | ChIP => {chr => shift}),
       CHR_LIST : karyotype passed to not to redefine it again,
       RD_DATA : path to RNA-DNA experiment file,
       AUX_DIR : path to auxilliary file (ABS-REL coord mapping) - directory, NOT file.
  '''
  PROT, ORG = PROT_ORG[0], PROT_ORG[1]
  chip_TADmap = f"{AUX_DIR}/{ORG}_crd_map.txt"
  fname = INIT_LOCS["RIP"].split('/')[-1]
  
  # DEFINE ANNOT LOCATION
  if ORG == "mESC":
    chromhmm_loc = "SHIFT_AUX_FILES/ChromHMM_mm10.mESC.bed"
  if ORG == "K562":
    chromhmm_loc = "SHIFT_AUX_FILES/ChromHMM_hg38.K562.bed"
    spinanno_loc = "SHIFT_AUX_FILES/K562_SPIN_25kb_hg38.bed"
    spin_annot_cmd = f"{intersectBed} -a stdin -b {spinanno_loc} -wb"

  # SHIFT RIP PEAKS
  SHIFT_RNA_Protein(INIT_LOCS["RIP"], f"{OUTDIR['RIP']}/{fname}", SHIFTS_DICT["RIP"], CHR_LIST)
  SHIFT_DNA_Protein(INIT_LOCS["ChIP"], f"{OUTDIR['ChIP']}/{PROT}.ChIP.{ORG}.bed", SHIFTS_DICT["ChIP"], chip_TADmap)

  # INTERSECT SHIFTED
  intsec_DP_cmd = f"{intersectBed} -a {RD_DATA} -b {OUTDIR['ChIP']}/{PROT}.ChIP.{ORG}.bed -wb"
  rearr_cmd = "awk \'$9 != \"protein_coding\" {print $4\"\t\"$5\"\t\"$6\"\t\"$7\"\t\"$1\"\t\"$2\"\t\"$3\"\t\"$8\"\t\"$9\"\t\"$10\"\t\"$11\"\t\"$19\"\t\"$20}\'"
  intsec_RP_cmd = f"{intersectBed} -a stdin -b {OUTDIR['RIP']}/{fname} -wb"
  awk_filt_cmd = "awk \'$4 == $17 {print $5\"\t\"$6\"\t\"$7}\'"
  wc_cmd = "wc -l"
  chromHmm_annot_cmd =  f"{intersectBed} -a stdin -b {chromhmm_loc} -wb"
  cut_sort_cmd = "cut -f 7 | sort -u | uniq -c"

  intersect_ChIP = sp.Popen(sh.split(intsec_DP_cmd),                                stdout = sp.PIPE, text = True)
  rearrange      = sp.Popen(sh.split(rearr_cmd),     stdin = intersect_ChIP.stdout, stdout = sp.PIPE, text = True)
  intersect_ChIP.stdout.close()
  intersect_RIP  = sp.Popen(sh.split(intsec_RP_cmd), stdin = rearrange.stdout,      stdout = sp.PIPE, text = True)
  rearrange.stdout.close()
  awk_filt       = sp.Popen(sh.split(awk_filt_cmd),  stdin = intersect_RIP.stdout,  stdout = sp.PIPE, text = True)
  intersect_RIP.stdout.close()

  ## ---------------------------------------------------------------------
  ## NEXT LINES ARE FOR ANNOTATION SAMPLING - (DE)COMMENT if needed
#  chrom_hmm_pro  = sp.Popen(sh.split(chromHmm_annot_cmd), stdin = awk_filt.stdout,      stdout = sp.PIPE, text = True)
#  awk_filt.stdout.close()
#  cut_sort_proc  = sp.Popen(cut_sort_cmd, shell=True,     stdin = chrom_hmm_pro.stdout, stdout = sp.PIPE, text = True)
#  chrom_hmm_pro.stdout.close()
#  state_count, _ = cut_sort_proc.communicate()

  ## ---------------------------------------------------------------------
  ## NEXT LINES ARE FOR ANNOTATION SPIN SAMPLING - (DE)COMMENT if needed
  spin_ann_pro = sp.Popen(sh.split(spin_annot_cmd),     stdin = awk_filt.stdout,     stdout = sp.PIPE, text = True)
  awk_filt.stdout.close()
  cut_sort_proc  = sp.Popen(cut_sort_cmd, shell=True,   stdin = spin_ann_pro.stdout, stdout = sp.PIPE, text = True)
  spin_ann_pro.stdout.close()
  state_count, _ = cut_sort_proc.communicate()

  return (ORG, state_count.replace('\n', '\t'))

  ## ---------------------------------------------------------------------
  ## NEXT LINES ARE FOR SIGNIFICANCE ESTIMATION - (DE)COMMENT if needed
#  wc_results     = sp.Popen(sh.split(wc_cmd),        stdin = awk_filt.stdout,  stdout = sp.PIPE, text = True)
#  awk_filt.stdout.close()

#  triad_count, _ = wc_results.communicate()
#  return (f"{PROT}_{ORG}", triad_count.strip())

  
