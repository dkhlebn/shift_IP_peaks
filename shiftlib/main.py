import sys
import glob
import time
import random
import logging
import shlex as sh
import subprocess as sp
import concurrent.futures as cf
import pandas as pd

from shift_intersect import SHIFT_THEN_INTERSECT


# READ CMD ARGUMENTS
ChIP_dir = sys.argv[1]
RIP_dir = sys.argv[2]
RD_dir = sys.argv[3]
TMP_dir = sys.argv[4]
ANNOT_FLAG = bool(sys.argv[5])
FUNCTIONAL = bool(sys.argv[6])
idx = TMP_dir.rsplit("_", maxsplit=1)[-1]

# SET UP LOGS

root = logging.getLogger()
root.setLevel(logging.DEBUG)
handler = logging.FileHandler(f"LOGS_SHIFT_{idx}.sh.txt")
handler.setLevel(logging.DEBUG)
root.addHandler(handler)
logging.info(f"ARGUMENTS passsed: {', '.join(sys.argv[1:])}...")

# COMMON DEFS
aux_files_dir = f"shift_IP_peaks/SHIFT_AUX_FILES"
chr_org = {
    "K562": ["chr" + str(val) for val in [*range(1, 23), "X", "Y"]],
    "mESC": ["chr" + str(val) for val in [*range(1, 20), "X", "Y"]],
}
K562_path = f"{aux_files_dir}/K562.compartments.bed"
mESC_path = f"{aux_files_dir}/mESC.compartments.bed"

# MAKE WORKING DIRS
sp.run(sh.split(f"mkdir {TMP_dir}"))
sp.run(sh.split(f"mkdir {TMP_dir}/RP_peaks"))
sp.run(sh.split(f"mkdir {TMP_dir}/DP_peaks"))
sp.run(sh.split(f"mkdir {TMP_dir}/SIM_TRIADS"))
logging.info("DIRS created, PREREQUISITES defined...")

# GET LIST OF PROTEINS AND ORGANISMS
PROT_ORG = []
for el in glob.glob(f"{RIP_dir}/*"):
    fname = el.split("/")[-1]
    prot = fname.split(".")[0]
    exp = fname.split(".")[1]
    org = fname.split(".")[2]
    if glob.glob(f"{ChIP_dir}/{prot}.*"):
        PROT_ORG.append((prot, exp, org))

mESC_RD_exps = [
    "grid_mm10_mESC.bed",
    "radicl2FA_mm10_mESC.bed",
    "radiclNPM_mm10_mESC.bed",
]
K562_RD_exps = ["redc_hg38_K562.bed"]

# DRIVER CODE
sim_res = []
for i in range(1):
    logging.info(
        f"GENERATION OF BG DATA started {i+1}-th time. Shifting, Permuting & Intersecting..."
    )
    init_time = time.time()

    # GENERATE SHIFT DICT FOR BOTH ORGANISMS
    shifts = {}
    shifts["K562"] = {ch: 0 for ch in chr_org["K562"]}
    shifts["mESC"] = {ch: 0 for ch in chr_org["mESC"]}
    for ch in chr_org["K562"]:
        shifts["K562"][ch] = random.choice([-1, 1]) * random.choice(
            [i * 10e6 for i in (1, 3, 5, 7, 10)]
        )
    for ch in chr_org["mESC"]:
        shifts["mESC"][ch] = random.choice([-1, 1]) * random.choice(
            [i * 10e6 for i in (1, 3, 5, 7, 10)]
        )

    # GENERATE ARGUMENTS FOR SHIFTS AND INTERSECTS
    exec_arguments = []
    for el in PROT_ORG:
        prot, exp, org = el[0], el[1], el[2]
        init_locs = {
            "RIP": glob.glob(f"{RIP_dir}/{prot}.{exp}.{org}.bed")[0],
            "ChIP": glob.glob(f"{ChIP_dir}/{prot}.ChIP.{org}.bed")[0],
        }
        out_locs = {"RIP": f"{TMP_dir}/RP_peaks", "ChIP": f"{TMP_dir}/DP_peaks"}
        rd_exps = mESC_RD_exps if org == "mESC" else K562_RD_exps
        rd_paths = [f"{RD_dir}/{exp}" for exp in rd_exps]
        exec_arguments.append(
            (
                el,
                init_locs,
                out_locs,
                shifts[org],
                rd_paths,
                aux_files_dir,
                ANNOT_FLAG,
                FUNCTIONAL,
            )
        )

    # SHIFT PEAKS AND INTERSECT
    with cf.ProcessPoolExecutor(max_workers=len(exec_arguments)) as executor:
        futures = {
            executor.submit(SHIFT_THEN_INTERSECT, *argts) for argts in exec_arguments
        }
    time_taken = time.time() - init_time
    logging.info(f"Shifts and intersects done (Time: {time_taken}). Writing!")
    for futr in cf.as_completed(futures):
        res_list = futr.result()
        for el in res_list:
            sim_res.append((str(i), *el))

    # CLEAN TMP DIRS
    removal_rip = f"rm {TMP_dir}/*P_peaks/* {TMP_dir}/SIM_TRIADS/* -rf"
    remove_sp = sp.Popen(sh.split(removal_rip))
    _, _ = remove_sp.communicate()

# CLEAN UP
clean_up = f"rm {TMP_dir} -rf"
clean_up_proc = sp.Popen(sh.split(clean_up))
_, _ = clean_up_proc.communicate()

# SAVE RESULTS
with open(f"sim_annot_{idx}.tsv", "w") as fout_handle:
    for el in sim_res:
        fout_handle.write("\t".join(el) + "\n")
