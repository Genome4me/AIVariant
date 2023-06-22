###
# Written by Junhak Ahn
# email : eyeahn11@snu.ac.kr
###

import time
import os

def get_contig_lengths_from_fai(fa_file_):
    contig_lengths = {}

    fa_file = open(fa_file_, 'r')
    for line in fa_file:
        cols = line.strip('\n').split('\t')
        contig = cols[0]
        contig_len = int(cols[1])
        contig_lengths[contig] = contig_len
    fa_file.close()

    return contig_lengths

def get_reference_info(fa_file_, default_chromosomes_for_AIVariant):
    whole_contigs = []
    contigs_for_AIVariant = []

    is_with_chr = False

    fa_file = open(fa_file_, 'r')
    for line in fa_file:
        cols = line.strip('\n').split('\t')
        contig = cols[0]
        whole_contigs.append(contig)
        if 'chr' in contig:
            is_with_chr = True
    fa_file.close()

    if is_with_chr:
        reference_type = 'with_chr'
        for chrom in default_chromosomes_for_AIVariant:
            chrom_with_chr = 'chr%s' % chrom
            if chrom_with_chr in whole_contigs:
                contigs_for_AIVariant.append(chrom_with_chr)
    else:
        reference_type = 'without_chr'
        for chrom in default_chromosomes_for_AIVariant:
            if chrom in whole_contigs:
                contigs_for_AIVariant.append(chrom)

    return contigs_for_AIVariant

def convert_to_vcf(aiv_call_file, out_vcf_path, reference, reference_index,
                   tumor_bam_norm, tumor_bam_unnorm, normal_bam_norm, normal_bam_unnorm):
    vcf_fixed_header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n"
    contig_lengths = get_contig_lengths_from_fai(reference_index)

    out_vcf_handle = open(out_vcf_path, 'w')
    out_vcf_handle.write('##fileformat=VCFv4.1\n')
    out_vcf_handle.write(f'##fileDate={time.ctime}\n')
    out_vcf_handle.write(f"##reference=file:{os.path.abspath(reference)}\n")
    out_vcf_handle.write(f"##TumorBAM_normalized={os.path.abspath(tumor_bam_norm)}\n")
    out_vcf_handle.write(f"##TumorBAM_unnormalized={os.path.abspath(tumor_bam_unnorm)}\n")
    out_vcf_handle.write(f"##NomralBAM_normalized={os.path.abspath(normal_bam_norm)}\n")
    out_vcf_handle.write(f"##NormalBAM_unnormalized={os.path.abspath(normal_bam_unnorm)}\n")

    for contig in contig_lengths:
        contig_len = contig_lengths[contig]
        out_vcf_handle.write(f"##contig=<ID={contig},length={contig_len}>\n")

    out_vcf_handle.write('##INFO=<ID=score,Number=1,Type=Float,Description="AIVariant score">\n')
    out_vcf_handle.write("##FILTER=<ID=.,Number=0,Type=String>\n")
    out_vcf_handle.write("##FORMAT=<ID=.,Number=0,Type=String>\n")
    out_vcf_handle.write(vcf_fixed_header)

    for line in open(aiv_call_file):
        if line.startswith('/'): continue
        cols = line.strip('\n').split('\t')  # 0-based
        contig = cols[0][3:]
        pos = cols[1]
        ref = cols[2]
        alt = cols[3]
        score = cols[5]
        pos_1based = int(pos) + 1
        info = f'score={score}'
        vcf_line = f'{contig}\t{pos_1based}\t.\t{ref}\t{alt}\t.\t.\t{info}\t.\n'
        out_vcf_handle.write(vcf_line)

    out_vcf_handle.close()