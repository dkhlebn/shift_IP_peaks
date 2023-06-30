import shlex as sh
import subprocess as sp


def _ANNOT_SPIN_K562(intersectBed, out_tri_fname, TMP_dir, cut_sort_cmd, states):
    spinanno_loc = "shift_IP_peaks/SHIFT_AUX_FILES/K562_SPIN_25kb_hg38.bed"
    spin_annot_cmd = (
        f"{intersectBed} -a {TMP_dir}/SIM_TRIADS/{out_tri_fname} -b {spinanno_loc} -wb"
    )

    spin_ann_proc = sp.Popen(sh.split(spin_annot_cmd), stdout=sp.PIPE, text=True)
    cut_sort_proc = sp.Popen(
        cut_sort_cmd,
        shell=True,
        stdin=spin_ann_proc.stdout,
        stdout=sp.PIPE,
        text=True,
    )
    spin_ann_proc.stdout.close()
    state_count_spin, _ = cut_sort_proc.communicate()
    if state_count_spin:
        ret_str = state_count_spin.replace("\n", "\t")
    else:
        ret_str = "\t".join([f"0 {state}" for state in states])
    return ret_str


def ANNOT_SIM_TRIADS(
    intersectBed,
    out_tri_fname,
    ORG,
    TMP_dir,
):
    """
    Aux function to have cleaner code, intersects
    simulated triads with ChromHMM annotation, and,
    if built for K562, with SPIN annotation.
    """

    ChromHMM_STATES = [
        "Enhancer",
        "HetChrom_CNV",
        "Insulator",
        "PRC_Repressed",
        "Promoter",
        "Transcribed",
    ]
    SPIN_STATES = [
        "A_compartment",
        "B_compartment",
        "Lamina",
        "Lamina_Like",
        "Near_Lm",
        "Speckle",
    ]

    if ORG == "mESC":
        chromhmm_loc = "shift_IP_peaks/SHIFT_AUX_FILES/ChromHMM_mm10.mESC.bed"
    if ORG == "K562":
        chromhmm_loc = "shift_IP_peaks/SHIFT_AUX_FILES/ChromHMM_hg38.K562.bed"

    chromHmm_annot_cmd = (
        f"{intersectBed} -a {TMP_dir}/SIM_TRIADS/{out_tri_fname} -b {chromhmm_loc} -wb"
    )
    cut_sort_cmd = "cut -f 22 | sort | uniq -c"

    # CHECK ANNOTATION IF NEEDED
    chrom_hmm_proc = sp.Popen(sh.split(chromHmm_annot_cmd), stdout=sp.PIPE, text=True)
    cut_sort_proc = sp.Popen(
        cut_sort_cmd,
        shell=True,
        stdin=chrom_hmm_proc.stdout,
        stdout=sp.PIPE,
        text=True,
    )
    chrom_hmm_proc.stdout.close()
    state_count_chromhmm, _ = cut_sort_proc.communicate()

    if state_count_chromhmm:
        chromhmm_annot = state_count_chromhmm.replace("\n", "\t")
    else:
        chromhmm_annot = "\t".join([f"0 {state}" for state in ChromHMM_STATES])

    if ORG == "K562":
        spin_annot = _ANNOT_SPIN_K562(
            intersectBed, out_tri_fname, TMP_dir, cut_sort_cmd, SPIN_STATES
        )
    else:
        spin_annot = "\t".join([f"NA {state}" for state in SPIN_STATES])

    ret_str = "\t".join([chromhmm_annot, spin_annot])
    return ret_str
