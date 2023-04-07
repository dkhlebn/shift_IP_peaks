import os
import math
import glob
import random
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

  # SHIFT RIP PEAKS
  SHIFT_RNA_Protein(INIT_LOCS["RIP"], f"{OUTDIR['RIP']}/fname", SHIFTS_DICT["RIP"], CHR_LIST)
  SHIFT_DNA_Protein(INIT_LOCS["ChIP"], f"{OUTDIR['ChIP']}/{PROT}.ChIP.{ORG}.bed", SHIFTS_DICT["ChIP"], chip_TADmap)

  # INTERSECT SHIFTED
  intsec_DP_cmd = f"{intersectBed} -a {RD_DATA} -b {OUTDIR['ChIP']}/{PROT}.ChIP.{ORG}.bed -wb"
  rearr_cmd = "awk '{print $4"\t"$5"\t"$6"\t"$7"\t"$1"\t"$2"\t"$3"\t"$8"\t"$9"\t"$10"\t"$11"\t"$19"\t"$20}'"
  intsec_RP_cmd = f"{intersectBed} -a stdin -b {OUTDIR['RIP']}/fname -wb"
  rear_filt_cmd = "awk '$4 == $17 {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$18"\t"$19}'"
  wc_cmd = "wc -l"

  intersect_ChIP = sp.Popen(sh.split(intsec_DP_cmd),                                stdout = sp.PIPE)
  rearr_1        = sp.Popen(sh.split(rearr_cmd),     stdin = intersect_ChIP.stdout, stdout = sp.PIPE)
  intersect_RIP  = sp.Popen(sh.split(intsec_RP_cmd), stdin = rearr_1.stdout,        stdout = sp.PIPE)
  rear_filt      = sp.Popen(sh.split(rear_filt_cmd), stdin = intersect_RIP.stdout,  stdout = sp.PIPE)
  wc_results     = sp.Popen(sh.split(wc_cmd),        stdin = rear_filt.stdout,      stdout = sp.PIPE, text = True)

  triad_count, _ = wc_results.communicate()
  return (f"{PROT}_{ORG}", triad_count)

  
