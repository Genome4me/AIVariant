import os

cmdl = 'python main.py --tumor data/tumor.bam --normal data/normal.bam --reference data/Homo_sapiens_assembly19.fasta --hg hg19 --out_dir run'    # output run directory would be like run_example
os.system(cmdl)
