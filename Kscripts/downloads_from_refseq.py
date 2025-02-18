#!/media/gsadmin/vd3/mengxiangfu/software/miniconda3/envs/workenv/bin/python
# -*- coding:utf-8 -*- 

from pprint import pprint
import subprocess
from multiprocessing import Pool


def download_seq(args):
    accession, out_path = args
    cmd_line = r"/media/gsadmin/vd3/mengxiangfu/software/EDirect/edirect/esearch -db nuccore -query '{0}' | \
        /media/gsadmin/vd3/mengxiangfu/software/EDirect/edirect/efetch -db nuccore -format fasta > {1}/{0}.fasta".format(accession, out_path)
    try:
        subprocess.run(cmd_line, shell=True, timeout=60, check=True)
    except subprocess.TimeoutExpired:
        print("TimeoutExpired: " + accession)
    except subprocess.CalledProcessError:
        print("CalledProcessError: " + accession)

if __name__ == '__main__':
    out_path = r'ref_seq/sequences/'
    TASKS = list()

    ## first
    with open(r'Organism_GBID.txt', 'rt') as f:
        for line in f:
            accession = line.strip().split('\t')[1]
            TASKS.append((accession, out_path))

    ## second
    # accession_list = [
    #     "NC_001803.1",
    #     "NC_006213.1",
    #     "NC_006577.2",
    #     "NZ_LN831051.1",
    #     "NC_002204.1",
    #     "NC_038311.1",
    # ]
    # for a in accession_list:
    #     TASKS.append((a, out_path))

    # pprint(TASKS)
    with Pool(4) as p:
        p.map(download_seq, TASKS)

