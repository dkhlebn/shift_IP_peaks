import glob
import shlex as sh
import subprocess as sp

def INTERSECT_W_TAD(ChIP_PATH, COMPTS_PATH, OUTPATH):

  intersectBed='/home/tools/bedtools/bedtools-v2.16.2/bin/intersectBed'

  cmd_inter = f"{intersectBed} -a {ChIP_PATH} -b {COMPTS_PATH} -wb"
  cmd_cut = f"cut -f 1,2,3,4,5,6,7,8,9,10,14"

  with open(f"{OUTPATH}", 'w') as stdout_direction:
    intersection = sp.Popen(sh.split(cmd_inter), stdout = sp.PIPE)
    cut = sp.Popen(sh.split(cmd_cut), stdin = intersection.stdout, stdout = stdout_direction)
    return cut.wait()
