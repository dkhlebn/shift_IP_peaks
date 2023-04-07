import os
import math
import glob
import random
import pandas as pd


def SHIFT_DNA_Protein(PEAKFILE_PATH, OUTFILE_PATH, SHIFTS_DICT, CMP_MAP_PATH):
    '''
    Shifts given ChIP-Seq peakfile, writes into outfile, uses passed seed.
    '''
    
    # MODIFIED BINSEARCH FOR TAD BORDER
    def lowerBorder(arr, threshold):
        low = 0
        high = len(arr)

        while low < high:
            mid = math.floor((low + high) / 2)
            if arr[mid][0] == threshold:
                return mid
            elif arr[mid][0] < threshold and mid != low:
                low = mid
            elif arr[mid][0] > threshold and mid != high:
                high = mid
            else:
                high = low = low + 1
        return low - 1

    # DEFINE MAIN SHIFT FUNCTION
    def shift_peak(shift_dict, peak_row, chr_comp_dct):
        '''
        Shift given peak according to chrN and A/B-compartment characteristics.
        '''

        # retrieve info and dict
        chrom, tad = peak_row[1]["chr"], peak_row[1]["tad"]
        start_init, end_init = peak_row[1]["start"], peak_row[1]["end"]
        dct = chr_comp_dct[chrom][tad]
        rel_crds = list(dct.keys())
        abs_crds = list(dct.values())
        rel_length = rel_crds[-1][1]
        shift = shift_dict[chrom]

        # find shifted relative coordinates; *_shift denotes the shifted coordinate
        init_tad = lowerBorder(abs_crds, start_init)
        start_shift = rel_crds[init_tad][0] + (start_init - abs_crds[init_tad][0]) + shift
        end_shift = rel_crds[init_tad][0] + (end_init - abs_crds[init_tad][0]) + shift

        # check if st/end_shift > total relative length (i.e. cycle over chr)
        start_shift = start_shift - rel_length * (start_shift > rel_length) + rel_length * (start_shift < 0)
        end_shift = end_shift - rel_length * (end_shift > rel_length) + rel_length * (end_shift < 0)

        # find compartment for *_shift. If newPeak is split between 2, create two peaks
        new_tad_st = lowerBorder(rel_crds, start_shift)
        new_tad_en = lowerBorder(rel_crds, end_shift)

        if new_tad_st == new_tad_en:
            start_new = abs_crds[new_tad_st][0] + (start_shift - rel_crds[new_tad_st][0])
            end_new = abs_crds[new_tad_en][0] + (end_shift - rel_crds[new_tad_en][0])
            peak_row[1]["start"] = int(start_new)
            peak_row[1]["end"] = int(end_new)
            return [peak_row[1]]
        else:
            peak_fromStart = (0, peak_row[1].copy())
            peak_fromStart[1]["start"] = int(abs_crds[new_tad_st][0] + (start_shift - rel_crds[new_tad_st][0]))
            peak_fromStart[1]["end"] = int(abs_crds[new_tad_st][1])

            peak_fromEnd = (0, peak_row[1].copy())
            peak_fromEnd[1]["start"] = int(abs_crds[new_tad_en][0])
            peak_fromEnd[1]["end"] = int(abs_crds[new_tad_en][0] + (end_shift - rel_crds[new_tad_en][0]))
            return [peak_fromStart[1], peak_fromEnd[1]]

    
    # LOAD TAD ABS/REL MAPPING
    chr2cmp2crds = ''
    with open(f"{CMP_MAP_PATH}", 'r') as f:
        for i in f.readlines():
            chr2cmp2crds = i
    chr2cmp2crds = eval(chr2cmp2crds)

    # PREPARE PEAKFILE - WITH TAD ANNOTATION
    columns = ["chr", "start", "end", "name", "score", "strand", "signalValue",
              "pval", "qval", "peak", "tad"]
    coldct = {i : col for i, col in enumerate(columns)}
    peaks = pd.read_csv(f"{PEAKFILE_PATH}", sep='\t', header = None).rename(columns = coldct)
    
    # SHIFT PEAKS
    res_l = []
    for row in peaks.iterrows():
        shifted = shift_peak(SHIFTS_DICT, row, chr2cmp2crds)
        for el in shifted:
            res_l.append(el)

    pd.DataFrame.from_records(res_l).to_csv(f"{OUTFILE_PATH}", header = False, index = False, sep='\t')
    return 0
