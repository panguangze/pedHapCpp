from pyfaidx import Fasta,Faidx
import sys
import os
import random
import argparse
import logging
import re
logging.basicConfig(level=logging.DEBUG)

PYTHON = "~/miniconda3/envs/py3/bin/python"
LOCALHAP = "/home/gzpan2/app/localhaptgs/debug/localHap"
SAMTOOLS = "~/app/samtools/bin/samtools"
CBC = "~/miniconda3/envs/py2/bin/cbc"
pbsim = "~/app/pbsim2/src/pbsim"
pbmodel = "~/app/pbsim2/data/P6C4.model"
nanomodel = "~/app/pbsim2/data/R95.model"
hpvpip_root = "~/app/hpvpipe"
faToTwoBit = "~/app/faToTwoBit"
computeGCBias = "~/miniconda3/envs/py3.6/bin/computeGCBias"
correctGCBias = "~/miniconda3/envs/py3.6/bin/correctGCBias"
samtools = "~/app/samtools/bin/samtools"
sim3c = "/home/gzpan2/.local/bin/sim3C"
ISPB=1
DRY = True
# ~/app/pbsim2/src/pbsim --depth 20 --prefix tgs --hmm_model ~/app/pbsim2/data/P6C4.model test.out.fa
def execmd(cmd):
    # print("Exec: {}".format(cmd))
    if DRY:
        logging.info(cmd)
    else:
        logging.info(cmd)
        os.system(cmd)
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--host_ref',
                        required=True,
                        help='host_ref')
    parser.add_argument('--out',
                        required=True,
                        help='out')
    parser.add_argument('--script_root',
                        required=True,
                        help='out')
    parser.add_argument('--simple_par',
                        required=True,
                        help='simple_par')
    parser.add_argument('--depth',
                        required=True,
                        help='simple_par')
    parser.add_argument('--pb',
                        required=True,
                        type=int,
                        default=1,
                        help='simple_par')
    args = parser.parse_args()
    if not os.path.exists(args.out):
        os.mkdir(args.out)
    # choose host chr
    # mk_fa(args.host_ref, host_chrs, args.v_ref, args.v_chr, args.out)
    # generate_var(host_chrs, args.v_chr, v_len,args.out,args.v_ref)
    simulate(args.host_ref,args.simple_par, args.out, args.script_root, args.depth, args.pb)
    ISPB = args.pb


def simulate(ref, par, script_root, out_dir, depth, type):
    ref_fa = Fasta(ref)
    sim_fa = ref_fa["chr1"]
    fa_len = len(sim_fa)
    if type == "NGS":
        reads_num = depth * fa_len / (150 * 2)
        cmd = "bash {}/ngs.sh {} {} {} {}".format(script_root, ref, out_dir, reads_num, order_file)
        execmd(cmd)
def check_dir(dir):
    if not os.path.exists(dir):
        os.makedirs(dir,exist_ok=True)



if __name__ == "__main__":
    main()
