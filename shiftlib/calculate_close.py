import glob
import pandas as pd


def CALC_CLOSENESS(RNA_LIST_DIR, threshold=2 / 3):
    """
    Calculates distance matrices by intersection sizes,
    define lists of close proteins for each protein in the set
    -----
    Argts: path to dir and optional threshold
    """

    def read_in(file):
        rnas = []
        fh = open(file, "r")
        for line in fh:
            rnas.append(line.strip())
        fh.close()
        return set(rnas)

    # READ RNA LISTS
    dt = {}
    for file in glob.glob(f"{RNA_LIST_DIR}/*"):
        prot = file.split("/")[1].split(".")[0]
        dt[prot] = read_in(file)

    # CALCULATE INTERSECTION SIZES
    sim_dt = {prot: {} for prot in dt.keys()}
    for pr1 in dt.keys():
        for pr2 in dt.keys():
            common1 = len(dt[pr1].intersection(dt[pr2]))
            sim_dt[pr1][pr2] = common1
    df = pd.DataFrame.from_dict(sim_dt)
    for col in df.columns:
        df[col] = df[col] / len(dt[col])

    # DEFINE CLOSE PROTEINS
    dt_wrt = {}
    for col in df.columns:
        ser = df[col]
        thr = ser.quantile(threshold)
        dt_wrt[col] = sorted(ser[ser >= thr].index.to_list())
    with open("SHIFT_AUX_FILES/close_proteins.txt", "w") as fh:
        fh.write(str(dt_wrt))

    return 0


def CALC_FUNCTIONAL():
    """
    Aux function to manually load functionally related proteins
    """
    remodellers = [
        "CBP",
        "CBX3",
        "CHD1",
        "CHD4",
        "CHD7",
        "CTCF",
        "DNMT1",
        "EZH2",
        "HDAC1",
        "HLTF",
        "LSD1",
        "MORC2",
        "PCAF",
        "PHF8",
        "RBBP5",
        "SAFB",
        "SAFB2",
        "SUZ12",
        "WDR5",
    ]
    transcribers = [
        "CHD7",
        "EWSR1",
        "FUS",
        "GTF2F1",
        "HLTF",
        "HNRNPK",
        "HNRNPUL1",
        "LARP7",
        "PCBP1",
        "RBM15",
        "TAF15",
        "U2AF2",
        "YBX3",
        "ZC3H8",
    ]
    processors = [
        "FUS",
        "HNRNPC",
        "HNRNPH",
        "HNRNPK",
        "HNRNPL",
        "HNRNPUL1",
        "ILF3",
        "KHSRP",
        "LARP7",
        "NONO",
        "NUP98",
        "PTBP1",
        "RBFOX2",
        "RBM15",
        "RBM22",
        "SRSF1",
        "SRSF7",
        "SRSF9",
        "TARDBP",
        "U2AF1",
        "U2AF2",
        "ZC3H11A",
    ]
    all_proteins = list(
        set(remodellers).union(set(transcribers)).union(set(processors))
    )
    lsts = [remodellers, transcribers, processors]

    dt = {}
    for protein in all_proteins:
        functionally_close = []
        for lst in [lst for lst in lsts if protein in lst]:
            functionally_close.extend([elem for elem in lst])
        dt[protein] = list(set(functionally_close))
    with open("SHIFT_AUX_FILES/func_proteins.txt", "w") as fh:
        fh.write(str(dt))
    return 0


RNA_list_for_proteins = sys.argv[1]
CALC_CLOSENESS(RNA_list_for_proteins)
CALC_FUNCTIONAL()
