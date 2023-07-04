import shlex as sh
import subprocess as sp


def CONSTRUCT_TRIADS(
    RD_DATA, RIP_TYPE, OUTDIR, PROT, ORG, TMP_dir, fname
):
    """
    Aux script to construct simulated triads file. Has technical arguments
    and shouldn't be used not in the shift_intersect.py script.
    """
    RD_TYPE = RD_DATA.split("/")[-1]
    out_tri_fname = f"{PROT}.{RD_TYPE}.{RIP_TYPE}.{ORG}.bed"
    wc_cmd = "wc -l {TMP_dir}/SIM_TRIADS/{out_tri_fname}"
    construct_proc = sp.Popen(
        sh.split(
            f"bash shift_IP_peaks/shiftlib/tmp_script.sh {RD_DATA} {OUTDIR['ChIP']}/{PROT}.ChIP.{ORG}.bed {OUTDIR['RIP']}/{fname}"
        ),
        stdout=sp.PIPE,
        text=True,
    )
    result, _ = construct_proc.communicate()
    with open(f"{TMP_dir}/SIM_TRIADS/{out_tri_fname}", "w") as triad_fn:
        print(result, file=triad_fn)
    wc_l = sp.Popen(sh.split(wc_cmd), stdout=sp.PIPE)
    triad_count, _ = wc_l.communicate()

    return (out_tri_fname, triad_count.decode("ascii"))
