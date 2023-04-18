import os
import math
import glob
import random
import logging
import pandas as pd


def SHIFT_RNA_Protein(PEAKFILE_PATH, OUTFILE_PATH, SHIFT_DICT, CHR_LIST):
    '''
    Generates a shift for RNA-protein peak file, uses passed seed
    ''' 
    #DEFINE MAIN SHIFT FUNCTION
    def shift_rip_peak(shifts_dict, peak, limits_dct):
        dw_limit = limits_dct[peak[1].chr][0]
        up_limit = limits_dct[peak[1].chr][1]
        ft_len = up_limit - dw_limit

        shift = shifts_dict[peak[1].chr]
        start_init, end_init = peak[1].start, peak[1].end

        start_shift = start_init + shift
        end_shift = end_init + shift

        # Did not invent anything better. Accounting for outbounds
        while start_shift < dw_limit:
          start_shift += ft_len
        
        while start_shift > up_limit:
          start_shift -= ft_len

        while end_shift < dw_limit:
          end_shift += ft_len

        while end_shift > up_limit:
          end_shift -= ft_len

        if start_shift > end_shift:
            peak_fromStart = (0, peak[1].copy())
            peak_fromStart[1]["start"] = int(start_shift)
            peak_fromStart[1]["end"] = up_limit

            peak_fromEnd = (0, peak[1].copy())
            peak_fromEnd[1]["start"] = dw_limit
            peak_fromEnd[1]["end"] = int(end_shift)
            return [peak_fromStart[1], peak_fromEnd[1]]
        else:
            peak[1]["start"] = int(start_shift)
            peak[1]["end"] = int(end_shift)
            return [peak[1]]

    
    # PREPARE PEAKFILE
    columns = ["chr", "start", "end", "strand", "pval", "qval"]
    coldct = {i : col for i, col in enumerate(columns)}
    peaks = pd.read_csv(f"{PEAKFILE_PATH}", sep='\t', header = None).rename(columns = coldct)

    # GENERATE PREREQUISITES
    minmax_chr = {ch : [0, 0] for ch in CHR_LIST}
    for ch in CHR_LIST:
        minmax_chr[ch][0] = peaks.query(f"chr == '{ch}'").start.min()
        minmax_chr[ch][1] = peaks.query(f"chr == '{ch}'").end.max()

    # SHIFT PEAKS AND WRITE IN OUTFILE
    data = []
    for row in peaks.iterrows():
        shifted = shift_rip_peak(SHIFT_DICT, row, minmax_chr)
        for el in shifted:
            data.append(el)
    pd.DataFrame.from_records(data).to_csv(f"{OUTFILE_PATH}", header = False, index = False, sep='\t')
    return 0
