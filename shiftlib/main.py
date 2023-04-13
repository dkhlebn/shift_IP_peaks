import glob
import random
import subprocess
import shlex as sh
import concurrent.futures as cf
import pandas as pd
from sys import argv

from annotTAD import INTERSECT_W_TAD
from shift_intersect import SHIFT_THEN_INTERSECT

ChIP_dir = argv[1]
RIP_dir = argv[2]
RD_dir = argv[3]
TMP_dir = argv[4]

# COMMON DEFS
aux_files_dir = "/home/dkhlebnikov/TRIADS/SHIFT_AUX_FILES"
chr_org = {"K562" : ["chr" + str(val) for val in [*range(1,23), "X", "Y"]],
           "mESC" : ["chr" + str(val) for val in [*range(1,20), "X", "Y"]]}
K562_path = f"{aux_files_dir}/K562.compartments.bed"
mESC_path = f"{aux_files_dir}/mESC.compartments.bed"

# MAKE WORKING DIRS
subprocess.run(sh.split(f"mkdir {TMP_dir}"))
subprocess.run(sh.split(f"mkdir {TMP_dir}/DP_tads"))
subprocess.run(sh.split(f"mkdir {TMP_dir}/RP_peaks"))
subprocess.run(sh.split(f"mkdir {TMP_dir}/DP_peaks"))

# GET LIST OF PROTEINS AND ORGANISMS
PROT_ORG = []
for el in glob.glob(f"{ChIP_dir}/*"):
  fname = el.split('.')[2]
  prot =  fname.split('.')[0]
  org =   fname.split('.')[2]
  PROT_ORG.append((prot, org))

  # RUN ChIP-SEQ INTERSECT WITH TADs
  if org == "K562":
    INTERSECT_W_TAD(el, K562_path, f"{TMP_dir}/DP_tads/{fname}")
  elif org == "mESC":
    INTERSECT_W_TAD(el, mESC_path, f"{TMP_dir}/DP_tads/{fname}")

# DRIVER CODE
sim_res = []
for i in tqdm(range(10_000)):
  
  # GENERATE SHIFT DICT FOR BOTH ORGANISMS
  random.seed(SEED)
  shifts = {}
  shifts["K562"] = { IP : {ch : 0 for ch in chr_org["K562"]} for IP in ("ChIP", "RIP")}
  shifts["mESC"] = { IP : {ch : 0 for ch in chr_org["mESC"]} for IP in ("ChIP", "RIP")}

  for ch in chr_org["K562"]:
    shifts["K562"]["ChIP"][ch] = random.choice([-1, 1]) * random.choice([i * 10e6 for i in (1, 3, 5, 7, 10)])
    shifts["K562"]["RIP"][ch] =  random.choice([-1, 1]) * random.choice([i * 10e6 for i in (1, 3, 5, 7, 10)])

  for ch in chr_org["mESC"]:
    shifts["mESC"]["ChIP"][ch] = random.choice([-1, 1]) * random.choice([i * 10e6 for i in (1, 3, 5, 7, 10)])
    shifts["mESC"]["RIP"][ch] =  random.choice([-1, 1]) * random.choice([i * 10e6 for i in (1, 3, 5, 7, 10)])


  # GENERATE ARGUMENTS FOR SHIFTS AND INTERSECTS
  exec_arguments = []
  for el in PROT_ORG:
    prot, org = PROT_ORG[0], PROT_ORG[1]
    init_locs = {"RIP" : glob.glob(f"{RIP_dir}/{prot}*{org}*bed")[0],
                 "ChIP" : glob.glob(f"{ChIP_dir}/{prot}.ChIP.{org}.bed")[0]}
    out_locs = {"RIP" : f"{TMP_dir}/RP_peaks",
                "ChIP" : f"{TMP_dir}/DP_peaks"}
    rd_path = f"{RD_dir}/grid_mm10_mESC.bed" if org == "mESC" else f"{RD_dir}/redc_hg38_K562.bed"
    exec_arguments.append((el, init_locs, out_locs, shifts[org], chr_org[org], rd_path, aux_files_dir))

  # SHIFT PEAKS AND INTERSECT
  with cf.ThreadPoolExecutor() as executor:
    futures = {executor.submit(SHIFT_INTERSECT, *argts) for argts in exec_arguments}
    
    for futr in cf.as_completed(futures):
      triad_cnt = futr.result()
      sim_res.append((i, triad_cnt[0], triad_cnt[1]))

pd.DataFrame(sim_res).to_csv("sim_res.tsv", sep='\t', header=False, index=False)

