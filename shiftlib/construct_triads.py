import shlex as sh
import subprocess as sp


def CONSTRUCT_TRIADS(
    intersectBed, RD_DATA, RIP_TYPE, OUTDIR, PROT, ORG, TMP_dir, fname
):
    """
    Aux script to construct simulated triads file. Has technical arguments
    and shouldn't be used not in the shift_intersect.py script.
    """
    RD_TYPE = RD_DATA.split("/")[-1]
    out_tri_fname = f"{PROT}.{RD_TYPE}.{RIP_TYPE}.{ORG}.bed"
    intsec_DP_cmd = (
        f"{intersectBed} -a {RD_DATA} -b {OUTDIR['ChIP']}/{PROT}.ChIP.{ORG}.bed -wao"
    )
    rearr_cmd = 'awk \'$23 != 0 {print $4"\t"$5"\t"$6"\t"$7"\t"$1"\t"$2"\t"$3"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$20"\t"$21"\t"$23}\''
    intsec_RP_cmd = f"{intersectBed} -a stdin -b {OUTDIR['RIP']}/{fname} -wao"
    awk_filt_cmd = 'awk \'($22 != 0) && ($4 == $19) {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$11"\t"$12"\t"$13"\t"$14"\t"$20"\t"$21"\t"$15"\t"$22"\t"$10}\''
    uniq_filt_cmd = "sort -u"
    wc_cmd = f"wc -l {TMP_dir}/SIM_TRIADS/{out_tri_fname}"

    # INTERSECT SHIFTED AND PERMUTED DATA
    intersect_ChIP = sp.Popen(sh.split(intsec_DP_cmd), stdout=sp.PIPE, text=True)
    rearrange = sp.Popen(
        sh.split(rearr_cmd), stdin=intersect_ChIP.stdout, stdout=sp.PIPE, text=True
    )
    intersect_ChIP.stdout.close()
    intersect_RIP = sp.Popen(
        sh.split(intsec_RP_cmd), stdin=rearrange.stdout, stdout=sp.PIPE, text=True
    )
    rearrange.stdout.close()
    awk_filt = sp.Popen(
        sh.split(awk_filt_cmd), stdin=intersect_RIP.stdout, stdout=sp.PIPE, text=True
    )
    intersect_RIP.stdout.close()
    print("Intersection done")
    with open(f"{TMP_dir}/SIM_TRIADS/{out_tri_fname}", "w") as triad_fn:
        sp.Popen(sh.split(uniq_filt_cmd), stdin=awk_filt.stdout, stdout=triad_fn)
    awk_filt.stdout.close()

    wc_l = sp.Popen(sh.split(wc_cmd), stdout=sp.PIPE)
    triad_count, _ = wc_l.communicate()

    return (out_tri_fname, triad_count)
