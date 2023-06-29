import os
import math
import glob
import random
import logging
import pandas as pd


def PERMUTE_RNA_Protein(PEAKFILE_PATH, OUTFILE_PATH, ORGANISM, FUNC_FLAG):
  """
  Driver function branching by organism
  """
  if ORGANISM == "K562":
    return _permute_K562(PEAKFILE_PATH, OUTFILE_PATH, FUNC_FLAG)
  elif ORGANISM == "mESC":
    return _permute_mESC(PEAKFILE_PATH, OUTFILE_PATH)


def _permute_K562(PEAKFILE_PATH, OUTFILE_PATH, FUNC_FLAG):
  """
  Permute RNA-protein peaks to generate background data for triads construction
  in K562 cell line
  """
  PEAKDIR = PEAKFILE_PATH.split('/')[-2]
  fname = PEAKFILE_PATH.split('/')[-1]
  all_proteins = set(['LARP7', 'TARDBP', 'MORC2', 'SAFB2', 'FUS', 'CTCF', 'TAF15',
                 'EWSR1', 'ZC3H8', 'ILF3', 'U2AF2', 'ZC3H11A', 'NUP98',
                 'SUZ12', 'GTF2F1', 'U2AF1', 'WDR5', 'HNRNPH', 'HNRNPL',
                 'PHF8', 'NONO', 'YBX3', 'CHD1', 'SAFB', 'HNRNPUL1', 'EZH2',
                 'SRSF1', 'LSD1', 'RBM22', 'SRSF7', 'PCBP1', 'PTBP1', 'HNRNPK',
                 'RBFOX2', 'DNMT1', 'SRSF9', 'HDAC1', 'PCAF', 'CBP', 'RBM15',
                 'RBBP5', 'HNRNPC', 'CHD7', 'CHD4', 'HLTF', 'KHSRP', 'CBX3'])
  if FUNC_FLAG:
    rel_file = "SHIFT_AUX_FILES/func_proteins.txt"
  else:
    rel_file = "SHIFT_AUX_FILES/close_proteins.txt"

  with open(rel_file, 'r') as fn:
    prot_closeness = eval(fn.readline())

  # construct related and unrelated peak datasets
  rellist, unrlist = [], []
  cols = ["chr", "start", "end", "strand", "rna_type", "name", "pval", "qval"]
  dt_cols = {k : v for k, v in enumerate(cols)}
  for relprot in prot_closeness[protein]:
      fname = glob.glob(f"{PEAKDIR}/{relprot}*K562*")[0]
      df = pd.read_table(fname, header=None).rename(columns=dt_cols)
      if relprot == protein:
        protein_peaks = df
      rellist.append(df)
  for unrprot in all_proteins.difference(prot_closeness[protein]):
    fname = glob.glob(f"{PEAKDIR}/{unrprot}*K562*")[0]
    df = pd.read_table(fname, header=None).rename(columns=dt_cols)
    unrlist.append(df)
  related = pd.concat(rellist, axis=0, ignore_index=True)
  unrelated = pd.concat(unrlist, axis=0, ignore_index=True)

  # draw peaks from RNA types
  draft_list = protein_peaks.rna_type.value_counts().to_dict()

  simulated_peaks = []
  for gene_type, peak_cnt in draft_list.items():
    rel_num = round(0.8 * peak_cnt)
    unr_num = peak_cnt - rel_num

    from_rel = (related
                .query(f"rna_type == '{gene_type}'")
                .sample(n=rel_num, replace=True))
    from_unr = (unrelated
                .query(f"rna_type == '{gene_type}'")
                .sample(n=unr_num, replace=True))
    simulated_peaks.extend([from_rel, from_unr])

  sim_peaks = (pd.concat(simulated_peaks, axis=0, ignore_index=True)
               .sort_values(by=["chr", "start", "end"])
               .iloc[:, [0,1,2,3,4,6,7]])
  sim_peaks.to_csv(f"{OUTFILE_PATH}/{fname}", index=False, header=False, sep="\t")
  return 0
