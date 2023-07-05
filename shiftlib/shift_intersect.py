from RIPShift import PERMUTE_RNA_Protein
from ChIPShift import SHIFT_DNA_Protein
from construct_triads import CONSTRUCT_TRIADS
from annot_simtriads import ANNOT_SIM_TRIADS

intersectBed = "/home/tools/bedtools/bedtools-v2.16.2/bin/intersectBed"


def SHIFT_THEN_INTERSECT(
    PROT_ORG,
    INIT_LOCS,
    OUTDIR,
    SHIFTS_DICT,
    RD_DATA,
    AUX_DIR,
    ANNOT_FLAG,
    FUNCTIONAL,
):
    """
    PROT_ORG : (protein, RNA-protein experiment, cell_line),
    INIT_LOCS : dict( RIP | ChIP => initial peak location),
    OUTDIR    : dict( RIP | ChIP => TMP_dir / RP | DP) - directory, NOT file,
    SHIFTS_DICT : dict( RIP | ChIP => {chr => shift}),
    RD_DATA : paths to RNA-DNA experiment file,
    AUX_DIR : path to auxilliary file (ABS-REL coord mapping) - directory, NOT file,
    ANNOT_FLAG : True/False indicating whether the annotation is needed,
    FUNCTIONAL : whether we use functional groups or not.
    """
    PROT, RIP_EXP, ORG = PROT_ORG[0], PROT_ORG[1], PROT_ORG[2]
    chip_TADmap = f"{AUX_DIR}/{ORG}_crd_map.txt"
    fname = INIT_LOCS["RIP"].split("/")[-1]
    TMP_dir = OUTDIR["RIP"].split("/")[0]

    # SHIFT ChIP PEAKS, PERMUTE RIP PEAKS
    PERMUTE_RNA_Protein(INIT_LOCS["RIP"], f"{OUTDIR['RIP']}", ORG, FUNCTIONAL)
    SHIFT_DNA_Protein(
        INIT_LOCS["ChIP"],
        f"{OUTDIR['ChIP']}",
        SHIFTS_DICT,
        chip_TADmap,
    )
    results = []
    for RD_experiment in RD_DATA:
        # CONSTRUCT SIMULATED TRIADS AND ASSESS QUANTITY
        out_tri_fname, tri_cnt = CONSTRUCT_TRIADS(
            RD_experiment, RIP_EXP, OUTDIR, PROT, ORG, TMP_dir, fname
        )
        # CHECK ANNOTATION IF NEEDED
        if ANNOT_FLAG:
            annot_cnts = ANNOT_SIM_TRIADS(intersectBed, out_tri_fname, ORG, TMP_dir)
        results.append([out_tri_fname, tri_cnt, annot_cnts])
    return results
