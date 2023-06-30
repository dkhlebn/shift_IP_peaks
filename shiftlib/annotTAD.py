import glob
import sys
import logging
import shlex as sh
import subprocess as sp


def INTERSECT_W_TAD(ChIP_PATH, COMPTS_PATH, OUTPATH):
    intersectBed = "/home/tools/bedtools/bedtools-v2.16.2/bin/intersectBed"

    cmd_inter = f"{intersectBed} -a {ChIP_PATH} -b {COMPTS_PATH} -wb"
    cmd_cut = "cut -f 1,2,3,4,5,6,7,8,9,10,14"

    with open(f"{OUTPATH}", "w") as stdout_direction:
        intersection = sp.Popen(sh.split(cmd_inter), stdout=sp.PIPE)
        cut = sp.Popen(
            sh.split(cmd_cut), stdin=intersection.stdout, stdout=stdout_direction
        )
        return cut.wait()


ChIP_dir = sys.argv[1]
TARGET_d = sys.argv[2]

aux_files_dir = "/home/dkhlebnikov/TRIADS/SHIFT_AUX_FILES"
K562_path = f"{aux_files_dir}/K562.compartments.bed"
mESC_path = f"{aux_files_dir}/mESC.compartments.bed"

for el in glob.glob(f"{ChIP_dir}/*"):
    fname = el.split("/")[-1]
    prot = fname.split(".")[0]
    org = fname.split(".")[2]
    logging.info(f"{prot, org}")

    # RUN ChIP-SEQ INTERSECT WITH TADs
    if org == "K562":
        INTERSECT_W_TAD(el, K562_path, f"{TARGET_d}/DP_tads/{fname}")
    elif org == "mESC":
        INTERSECT_W_TAD(el, mESC_path, f"{TARGET_d}/DP_tads/{fname}")
logging.info("A/BS and ChIP-Seq intersected...")
