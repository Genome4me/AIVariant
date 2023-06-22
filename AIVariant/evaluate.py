###
# Written by Junhak Ahn
# email : eyeahn11@snu.ac.kr
###

import sys
import os
import argparse
import time
sys.path.append('./bin')
# from low_qual_filter import filter_etc_bam
# import bq_normalize
# import image_w_unnorm
from vcf import convert_to_vcf, get_reference_info
from setting import max_process_num, chromosomes, ref_bq_dist, \
    window_size, candsearch_cpp, epi_dir, aivariant_call, gpu_id


def main():
    # 0. Input
    parser = argparse.ArgumentParser(prog='AIVariant program')
    parser.add_argument('--tumor', '-t', help='tumor bam file path',
                        required=True)
    parser.add_argument('--normal', '-n', help='normal bam file path',
                        required=True)
    parser.add_argument('--reference', '-r', help='reference genome path',
                        required=True)
    parser.add_argument('--out_dir', '-o',
                        help='output variant call VCF file path',
                        required=True)
    parser.add_argument('--hg_version', '-g',
                        help='human genome version : hg19 / hg38',
                        required=True)
    parser.add_argument('--depth', '-d', help='bam file average depth',
                        default=30)
    args = parser.parse_args()

    tumor_bam = args.tumor
    normal_bam = args.normal
    reference = args.reference
    depth = args.depth
    hg_version = args.hg_version
    out_dir = args.out_dir

    reference_index = reference + '.fai'
    tumor_bam_index = tumor_bam + '.bai'
    normal_bam_index = normal_bam + '.bai'

    os.makedirs(out_dir, exist_ok=True)
    if tumor_bam[-4:] != '.bam' or normal_bam[-4:] != '.bam':
        sys.exit(
            f"Input bam file {tumor_bam} or {normal_bam} does not ends with '.bam'")
    if not os.path.isfile(tumor_bam_index):
        sys.exit(f"{tumor_bam_index} file does not exists")
    if not os.path.isfile(normal_bam_index):
        sys.exit(f"{normal_bam_index} file does not exists")
    if not os.path.isfile(reference_index):
        sys.exit(f"{reference_index} file does not exists")
    if hg_version != 'hg19' and hg_version != 'hg38':
        sys.exit(f"{hg_version} not supported : hg19 or hg38 required")

    tumor_bam_name = tumor_bam.split('/')[-1][:-4]
    normal_bam_name = normal_bam.split('/')[-1][:-4]

    contigs_for_AIVariant = get_reference_info(reference_index, chromosomes)

    # # 1. Low quality filter
    # print(f"[{time.ctime()}] AIVariant program starts")
    # print(f"[{time.ctime()}] Low qulity filter starts")

    # '''1)Tumor'''
    tumor_bam_filtered = os.path.join(out_dir,
                                      tumor_bam_name + '_filtered.bam')
    # filter_etc_bam(tumor_bam, tumor_bam_filtered)

    # '''2)Normal'''
    normal_bam_filtered = os.path.join(out_dir,
                                       normal_bam_name + '_filtered.bam')
    # filter_etc_bam(normal_bam, normal_bam_filtered)

    # print(f"[{time.ctime()}] Low qulity filter ends")

    # # 2. base-quality normalize
    # print(f"[{time.ctime()}] BQnorm starts")
    # bqnorm_dir = os.path.join(out_dir, 'bqnorm')
    # os.makedirs(bqnorm_dir, exist_ok=True)

    # '''1)Tumor'''
    # tumor_bq_dist = os.path.join(bqnorm_dir, f'{tumor_bam_name}_dist.pkl')
    # tumor_lookup_table = os.path.join(bqnorm_dir,
    #                                   f'{tumor_bam_name}_lookup.pkl')
    # bq_normalize.BamBqScanner(tumor_bam_filtered, reference,
    #                           contigs_for_AIVariant, tumor_bam_name,
    #                           bqnorm_dir, max_process_num).run()
    # bq_normalize.LookupTableCreator(ref_bq_dist, tumor_bq_dist, tumor_bam_name,
    #                                 bqnorm_dir).run()
    # bq_normalize.BamBqNormalizer(tumor_bam_filtered, tumor_lookup_table,
    #                              max_process_num, out_dir,
    #                              tumor_bam_name + '_normalized',
    #                              chromosome=None).run()

    # '''2)Normal'''
    # normal_bq_dist = os.path.join(bqnorm_dir, f'{normal_bam_name}_dist.pkl')
    # normal_lookup_table = os.path.join(bqnorm_dir,
    #                                    f'{normal_bam_name}_lookup.pkl')
    # bq_normalize.BamBqScanner(normal_bam_filtered, reference,
    #                           contigs_for_AIVariant, normal_bam_name,
    #                           bqnorm_dir, max_process_num).run()
    # bq_normalize.LookupTableCreator(ref_bq_dist, normal_bq_dist,
    #                                 normal_bam_name, bqnorm_dir).run()
    # bq_normalize.BamBqNormalizer(normal_bam_filtered, normal_lookup_table,
    #                              max_process_num, out_dir,
    #                              normal_bam_name + '_normalized',
    #                              chromosome=None).run()

    # print(f"[{time.ctime()}] BQnorm ends")

    # # 3. candidate search
    # print(f"[{time.ctime()}] candidate search starts")
    tumor_bam_normalized = os.path.join(out_dir,
                                        tumor_bam_name + '_normalized.bam')
    normal_bam_normalized = os.path.join(out_dir,
                                         normal_bam_name + '_normalized.bam')

    # candidates = []
    # for chrom in contigs_for_AIVariant:
    #     candidate_dir = os.path.join(out_dir, 'candSearch', chrom)
    #     os.makedirs(candidate_dir, exist_ok=True)
    #     candidate_file = os.path.join(candidate_dir, 'candidates.txt')
    #     candidates.append(candidate_file)

    #     cmdl = f'{candsearch_cpp} -t {tumor_bam_normalized} -n {normal_bam_normalized}' \
    #            f' -x {tumor_bam_filtered} -u {normal_bam_filtered} ' \
    #            f'-r {reference} -c chr{chrom} -o {candidate_dir} -j {max_process_num}' \
    #            f' -w {window_size} -d {depth} -g train -p forward'
    #     os.system(cmdl)

    # print(f"[{time.ctime()}] candidate search ends")

    # 4. make Image
    # print(f"[{time.ctime()}] makeImage starts")

    image_base_dir = os.path.join(out_dir, 'image')
    os.makedirs(image_base_dir, exist_ok=True)
    eval_image_path_file = os.path.join(image_base_dir, 'image_path4eval.txt')

    eval_image_path_handle = open(eval_image_path_file, 'w')
    for chrom in contigs_for_AIVariant:
        image_dir = os.path.join(image_base_dir, chrom)
        os.makedirs(image_dir, exist_ok=True)
        candidate_file = os.path.join(out_dir, 'candSearch', chrom,
                                      'candidates.txt')
        candidate_file_handle = open(candidate_file, 'r')
        candidate_file_handle.readline()  # header
        for line in candidate_file_handle:
            contig, pos, ref, alt = line.strip('\n').split('\t')
            if 'I' in alt or 'D' in alt: continue
            image_path = os.path.join(os.path.abspath(image_dir),
                                      f'{contig}:{pos}:{ref}:{alt}.pkl')
    #         image_generator = image_w_unnorm.ImageGenerator(
    #             tumor_bam_normalized, normal_bam_normalized, reference,
    #             contig, int(pos), ref, alt,
    #             tumor_bam_filtered, normal_bam_filtered, image_dir,
    #             hg_version=hg_version,
    #             epi_path=epi_dir)
    #         image_generator.run()
            eval_image_path_handle.write(image_path + '\n')

    #     candidate_file_handle.close()
    eval_image_path_handle.close()

    # print(f"[{time.ctime()}] makeImage ends")

    # 5. variant call
    print(f"[{time.ctime()}] variant call starts")
    call_dir = os.path.join(out_dir, 'call')
    os.makedirs(call_dir, exist_ok=True)
    aivariant_call_file = os.path.join(call_dir, 'preds.tsv')
    aivariant_vcf_file = os.path.join(call_dir, 'preds.vcf')
    cmdl = f'python {aivariant_call} --gpu_id {gpu_id} --test_data_file {eval_image_path_file} --prediction_file {aivariant_call_file}'
    os.system(cmdl)
    convert_to_vcf(aivariant_call_file, aivariant_vcf_file, reference,
                   reference_index,
                   tumor_bam_normalized, tumor_bam_filtered,
                   normal_bam_normalized, normal_bam_filtered)

    print(f"[{time.ctime()}] variant call ends")
    print(f"[{time.ctime()}] AIVariant program ends")


if __name__ == '__main__':
    main()

