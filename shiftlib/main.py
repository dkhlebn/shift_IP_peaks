import sys
import glob
import time
import random
import logging
import shlex as sh
import pandas as pd
import subprocess as sp
import concurrent.futures as cf

from shift_intersect import SHIFT_THEN_INTERSECT


# READ CMD ARGUMENTS
ChIP_dir = sys.argv[1]
RIP_dir = sys.argv[2]
RD_dir = sys.argv[3]
TMP_dir = sys.argv[4]
DPt_dir = sys.argv[5]
idx = TMP_dir.split('_')[-1]

# SET UP LOGS

root = logging.getLogger()
root.setLevel(logging.DEBUG)
handler = logging.FileHandler(f"LOGS_SHIFT_{idx}.sh.txt")
handler.setLevel(logging.DEBUG)
root.addHandler(handler)
logging.info(f"ARGUMENTS passsed: {ChIP_dir, RIP_dir, RD_dir, TMP_dir, DPt_dir}...")

# COMMON DEFS
aux_files_dir = "/home/dkhlebnikov/TRIADS/SHIFT_AUX_FILES"
chr_org = {"K562" : ["chr" + str(val) for val in [*range(1,23), "X", "Y"]],
           "mESC" : ["chr" + str(val) for val in [*range(1,20), "X", "Y"]]}
K562_path = f"{aux_files_dir}/K562.compartments.bed"
mESC_path = f"{aux_files_dir}/mESC.compartments.bed"

# MAKE WORKING DIRS
sp.run(sh.split(f"mkdir {TMP_dir}"))
sp.run(sh.split(f"mkdir {TMP_dir}/RP_peaks"))
sp.run(sh.split(f"mkdir {TMP_dir}/DP_peaks"))
logging.info("DIRS created, PREREQUISITES defined...")

# GET LIST OF PROTEINS AND ORGANISMS
PROT_ORG = []
for el in glob.glob(f"{ChIP_dir}/*"):
  fname = el.split('/')[1]
  prot =  fname.split('.')[0]
  org =   fname.split('.')[2]
  PROT_ORG.append((prot, org))

# DRIVER CODE
sim_res = []
for i in range(20):
  logging.info(f"SHIFTING started {i+1}-th time. Shifting & Intersecting...")
  init_time = time.time()

  # GENERATE SHIFT DICT FOR BOTH ORGANISMS
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
    prot, org = el[0], el[1]
    init_locs = {"RIP" : glob.glob(f"{RIP_dir}/{prot}*{org}*bed")[0],
                 "ChIP" : glob.glob(f"{DPt_dir}/{prot}.ChIP.{org}.bed")[0]}
    out_locs = {"RIP" : f"{TMP_dir}/RP_peaks",
                "ChIP" : f"{TMP_dir}/DP_peaks"}
    rd_path = f"{RD_dir}/grid_mm10_mESC.bed" if org == "mESC" else f"{RD_dir}/redc_hg38_K562.bed"
    exec_arguments.append((el, init_locs, out_locs, shifts[org], chr_org[org], rd_path, aux_files_dir))

  # SHIFT PEAKS AND INTERSECT
  with cf.ProcessPoolExecutor(max_workers=72) as executor:
    futures = {executor.submit(SHIFT_THEN_INTERSECT, *argts) for argts in exec_arguments}
  
  time_taken = time.time() - init_time
  logging.info(f"Shifts and intersects done (Time: {time_taken}). Writing!")    
  for futr in cf.as_completed(futures):
    triad_cnt = futr.result()
    sim_res.append((i, triad_cnt[0], triad_cnt[1]))

  removal_rip =  f"rm {TMP_dir}/*P_peaks/* -rf"
  remove_sp = sp.Popen(sh.split(removal_rip))
  _, _ = remove_sp.communicate()

# CLEAN UP
clean_up =  f"rm {TMP_dir} -rf"
clean_up_proc = sp.Popen(sh.split(clean_up))
_, _ = clean_up_proc.communicate()

# SAVE RESULTS
pd.DataFrame(sim_res).to_csv(f"sim_res_{idx}.tsv", sep='\t', header=False, index=False)

