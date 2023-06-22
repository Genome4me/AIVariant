###
# Written by Junhak Ahn
# email : eyeahn11@snu.ac.kr
###

import os
from setting import samtools
import pysam


def filter_etc_bam(in_file_bam, out_file_bam):
    in_bam = pysam.AlignmentFile(in_file_bam)
    out_bam = pysam.AlignmentFile(out_file_bam, 'wb', template=in_bam)

    for read in in_bam.fetch():
        if read.is_unmapped or read.is_secondary or read.is_qcfail or read.is_duplicate:
            continue

        mq = read.mapping_quality
        if mq == 0:
            continue

        out_bam.write(read)

    in_bam.close()
    out_bam.close()

    bam_index(out_file_bam)

def bam_index(in_file_bam):
    cmdl = f'{samtools} index {in_file_bam}'
    os.system(cmdl)

