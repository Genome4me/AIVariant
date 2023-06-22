import os

bin_dir = './bin'
ref_bq_dist = os.path.join(bin_dir, 'ref_dist_v2.pkl')  # 'reference base quality distribution'
candsearch_cpp = os.path.join(bin_dir, 'candidate_search')
# candsearch_cpp = os.path.join(bin_dir, 'candidate_search_221211')
window_size = 500000
epi_dir = os.path.join(bin_dir, 'encode')
aivariant_call = os.path.join(bin_dir, 'evaluation.py')


"""
Path settings
#TODO: Modify below to fit your environment
"""
# samtools = '~/baeklab/Junhak/anaconda3/envs/AIVariant/bin/samtools'   # v1.9
samtools = 'samtools'
max_process_num = 16
chromosomes = [str(i) for i in range(1,23)] + ['X', 'Y']
gpu_id = 0
