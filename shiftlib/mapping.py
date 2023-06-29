import os
import math
import glob
import random
import pandas as pd


def GENERATE_ABSREL_MAPPING(ORGANISM, CMPTDIR):
    """
    Generate absolute to relative coordinates mapping for ChIP-Seq shifting
    """

    def make_coords(tad_df, chrom, cm_type):
        """
        Returns a dict mapping real A/B compartment start to relative A/B compartment
        start position. Keys are absolute start positions, values are relative start
        positions.
        """
        tmp = tad_df.query(f"chr == '{chrom}' and compt == '{cm_type}'")
        ranges = [(a, b) for a, b in zip(tmp.start, tmp.end)]
        shift = ranges[0][0]
        ranges_new = [(ranges[0][0] - shift, ranges[0][1] - shift)]
        for i, val in enumerate(ranges[1:], 1):
            shift = shift + (val[0] - ranges[i - 1][1])
            ranges_new.append((val[0] - shift, val[1] - shift))
        return dict(zip(ranges_new, map(lambda x: x, ranges)))

    # LOAD A/B compartments
    TADs = pd.read_csv(
        f"{CMPTDIR}/{ORGANISM}_compartments.bed", sep="\t", header=None
    ).rename(columns={0: "chr", 1: "start", 2: "end", 3: "compt"})

    # generate dict prerequisites
    chrlist = ["chr" + str(i) for i in [*range(1, 20), "X", "Y"]]
    if ORGANISM == "K562":
        chrlist.extend(["chr20", "chr21", "chr22"])

    chr2cmp2crds = {chrom: {"A": {}, "B": {}} for chrom in chrlist}

    for chrom in chrlist:
        for cmp_type in ["A", "B"]:
            chr2cmp2crds[chrom][cmp_type] = make_coords(TADs, chrom, cmp_type)

    with open(f"{CMPTDIR}/{ORGANISM}_crd_map.txt", "w") as f:
        f.write(str(chr2cmp2crds))
    return 0


cmp_dir = sys.argv[1]

for org in ["K562", "mESC"]:
    GENERATE_ABSREL_MAPPING(org, cmp_dir)
