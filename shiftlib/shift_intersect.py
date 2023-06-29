import os
import sys
import math
import glob
import random
import logging
import pandas as pd
import shlex as sh
import subprocess as sp

from RIPShift import PERMUTE_RNA_Protein
from ChIPShift import SHIFT_DNA_Protein
intersectBed="/home/tools/bedtools/bedtools-v2.16.2/bin/intersectBed"

def SHIFT_THEN_INTERSECT(PROT_ORG, INIT_LOCS, OUTDIR, SHIFTS_DICT,
                         CHR_LIST, RD_DATA, AUX_DIR, ANNOT_FLAG,
                         FUNCTIONAL):
  '''
       PROT_ORG : (protein, RNA-protein experiment, cell_line),
       INIT_LOCS : dict( RIP | ChIP => initial peak location),
       OUTDIR    : dict( RIP | ChIP => TMP_dir / RP | DP) - directory, NOT file,
       SHIFTS_DICT : dict( RIP | ChIP => {chr => shift}),
       CHR_LIST : karyotype passed to not to redefine it again,
       RD_DATA : path to RNA-DNA experiment file,
       AUX_DIR : path to auxilliary file (ABS-REL coord mapping) - directory, NOT file,
       ANNOT_FLAG : True/False indicating whether the annotation is needed,
       FUNCTIONAL : whether we use functional groups or not.
  '''
  PROT, RIP_EXP, ORG = PROT_ORG[0], PROT_ORG[1], PROT_ORG[2]
  chip_TADmap = f"{AUX_DIR}/{ORG}_crd_map.txt"
  fname = INIT_LOCS["RIP"].split('/')[-1]
  TMP_dir = OUTDIR["RIP"].split('/')[0]

  # DEFINE ANNOT LOCATION
  if ORG == "mESC":
    chromhmm_loc = "SHIFT_AUX_FILES/ChromHMM_mm10.mESC.bed"
  if ORG == "K562":
    chromhmm_loc = "SHIFT_AUX_FILES/ChromHMM_hg38.K562.bed"
    spinanno_loc = "SHIFT_AUX_FILES/K562_SPIN_25kb_hg38.bed"
    spin_annot_cmd = f"{intersectBed} -a {TMP_dir}/SIM_TRIADS/{fname}.tri.bed -b {spinanno_loc} -wb"

  # SHIFT RIP PEAKS
  PERMUTE_RNA_Protein(INIT_LOCS["RIP"], f"{OUTDIR['RIP']}/{fname}", ORG, FUNCTIONAL)
  SHIFT_DNA_Protein(INIT_LOCS["ChIP"], f"{OUTDIR['ChIP']}/{PROT}.ChIP.{ORG}.bed", SHIFTS_DICT["ChIP"], chip_TADmap)

  # DEFINE SUBPROCESS CMDS TO RUN
  intsec_DP_cmd = f"{intersectBed} -a {RD_DATA} -b {OUTDIR['ChIP']}/{PROT}.ChIP.{ORG}.bed -wao"
  rearr_cmd = "awk \'$23 != 0 {print $4\"\t\"$5\"\t\"$6\"\t\"$7\"\t\"$1\"\t\"$2\"\t\"$3\"\t\"$8\"\t\"$9\"\t\"$10\"\t\"$11\"\t\"$12\"\t\"$20\"\t\"$21\"\t\"$23}\'"
  intsec_RP_cmd = f"{intersectBed} -a stdin -b {OUTDIR['RIP']}/{fname} -wao"
  awk_filt_cmd = "awk \'($22 != 0) && ($4 == $19) {print $1\"\t\"$2\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6\"\t\"$7\"\t\"$8\"\t\"$9\"\t\"$11\"\t\"$12\"\t\"$13\"\t\"$14\"\t\"$20\"\t\"$21\"\t\"$15\"\t\"$22\"\t\"$10}\'"
  uniq_filt_cmd = "sort -u"
  wc_cmd = f"wc -l {TMP_dir}/SIM_TRIADS/{fname}.tri.bed"
  chromHmm_annot_cmd =  f"{intersectBed} -a {TMP_dir}/SIM_TRIADS/{fname}.tri.bed -b {chromhmm_loc} -wb"
  cut_sort_cmd = "cut -f 22 | sort | uniq -c"

  # INTERSECT SHIFTED AND PERMUTED DATA
  intersect_ChIP = sp.Popen(sh.split(intsec_DP_cmd),                                stdout = sp.PIPE, text = True)
  rearrange      = sp.Popen(sh.split(rearr_cmd),     stdin = intersect_ChIP.stdout, stdout = sp.PIPE, text = True)
  intersect_ChIP.stdout.close()
  intersect_RIP  = sp.Popen(sh.split(intsec_RP_cmd), stdin = rearrange.stdout,      stdout = sp.PIPE, text = True)
  rearrange.stdout.close()
  awk_filt       = sp.Popen(sh.split(awk_filt_cmd),  stdin = intersect_RIP.stdout,  stdout = sp.PIPE, text = True)
  intersect_RIP.stdout.close()
  with open(f"{TMP_dir}/SIM_TRIADS/{fname}.tri.bed", 'w') as triad_fn:
    uniq_filt    = sp.Popen(sh.split(uniq_filt_cmd), stdin = awk_filt.stdout,       stdout = triad_fn)
  awk_filt.stdout.close()

  # CHECK ANNOTATION IF NEEDED
  if ANNOT_FLAG:
    chrom_hmm_pro  = sp.Popen(sh.split(chromHmm_annot_cmd),                               stdout = sp.PIPE, text = True)
    cut_sort_proc  = sp.Popen(cut_sort_cmd, shell=True,     stdin = chrom_hmm_pro.stdout, stdout = sp.PIPE, text = True)
    chrom_hmm_pro.stdout.close()
    state_count_chromhmm, _ = cut_sort_proc.communicate()

    if ORG == "K562":
      spin_ann_pro = sp.Popen(sh.split(spin_annot_cmd),                                  stdout = sp.PIPE, text = True)
      awk_filt.stdout.close()
      cut_sort_proc  = sp.Popen(cut_sort_cmd, shell=True,   stdin = spin_ann_pro.stdout, stdout = sp.PIPE, text = True)
      spin_ann_pro.stdout.close()
      state_count_spin, _ = cut_sort_proc.communicate()

      ret_str = (state_count_chromhmm + '\n' + state_count_spin).replace('\n', '\t')
      return (f"{PROT}_{ORG}_{RIP_EXP}", ret_str)

    return (f"{PROT}_{ORG}_{RIP_EXP}", state_count_chromhmm.replace('\n', '\t'))

  # COUNT TRIADS FOR SIGNIFICANCE ESTIMATION
  wc_results     = sp.Popen(sh.split(wc_cmd),  stdout = sp.PIPE, text = True)
  triad_count, _ = wc_results.communicate()
  return (f"{PROT}_{ORG}_{RIP_EXP}", triad_count.strip())
