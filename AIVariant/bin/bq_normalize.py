###
# Written by Junhak Ahn
# email : eyeahn11@snu.ac.kr
###

from typing import *
import time
import pickle
import os
from multiprocessing import Pool
import random
import shutil
import math
import operator
import pysam
from setting import samtools

FLOATING_POINT_DIGIT_LIMIT = 16

class ChromBqScanner():
    """Scan base qualities of reads aligned to bam file for a given chromosome
    to generate base quality distribution.
    """

    def __init__(self, input_bam_file: str, chromosome: str, out_file: str,
                 total_read_count: int, base_num_to_scan: int,
                 read_count_per_batch: int):
        """Create Scanner to scan base qualities of an input bam file
        for a given chromosome

        param input_bam_file: path of a bam file for base quality scanning
        param chromosome: a target chromosome for base quality scanning
        param out_file: path of output base quality count file
        param total_read_count: total read count of bam file for a given
                                chromosome
        param base_num_to_scan: total base number (base quality number) to scan
        param read_count_per_batch: number of reads for each scanning batch
        """

        self.input_bam_file = input_bam_file
        self.chromosome = chromosome
        self.out_file = out_file

        self.total_read_count = total_read_count
        self.base_num_to_scan = base_num_to_scan
        self.read_count_per_batch = read_count_per_batch

        self.batch_read_count = None
        self.read_to_scan = {}
        self.scanned_read_count = 0

        self.bqs = []
        self.bq_counts = {}

    def add_bq(self, read):
        """Append base qualities of given read without 'None' values

        param read: a given read to scan base qualities
        """

        for bq in read.query_qualities:
            if bq is not None:
                self.bqs.append(bq)

    def generate_read_batch(self):
        """Generate a read batch for next base quality scanning
        """

        not_checked_reads = [read_idx for read_idx in self.read_to_scan
                             if self.read_to_scan[read_idx] is None]

        if len(not_checked_reads) < self.read_count_per_batch:
            self.batch_read_count = len(not_checked_reads)
        else:
            self.batch_read_count = self.read_count_per_batch

        read_batch = random.sample(not_checked_reads, self.batch_read_count)

        for read_idx in read_batch:
            self.read_to_scan[read_idx] = True

    def init_scanner(self):
        """Initialize scanning states for total reads to None
        """

        for read_idx in range(self.total_read_count):
            self.read_to_scan[read_idx] = None

        self.batch_read_count = self.read_count_per_batch

    def scan(self):
        """Scan base qualities of reads in the given read batch
        Scanning state : True (need to scan), False (already scanned),
                         None (no need to considered, just pass)
        """

        input_bam_handle = pysam.AlignmentFile(self.input_bam_file, "rb")

        while len(self.bqs) < self.base_num_to_scan and self.batch_read_count == self.read_count_per_batch:
            self.generate_read_batch()
            for read_idx, read in enumerate(input_bam_handle.fetch(self.chromosome)):
                if self.read_to_scan[read_idx]:
                    self.add_bq(read)
                    self.scanned_read_count += 1
                    self.read_to_scan[read_idx] = False

        input_bam_handle.close()

        del self.read_to_scan

    def postprocess(self):
        """Over-sample or remove scanned base qualities to fit given
        'base_num_to_scan' and count the number of each base quality score
        """

        scanned_base_num = len(self.bqs)

        if scanned_base_num < self.base_num_to_scan:
            print(f"{self.chromosome} "
                  f": Not enough read count to scan, oversampling proceeds")
            self.bqs = random.choices(self.bqs, k=self.base_num_to_scan)

        elif scanned_base_num > self.base_num_to_scan:
            for i in range(scanned_base_num - self.base_num_to_scan):
                self.bqs.pop()

        for bq in self.bqs:
            try:
                self.bq_counts[bq] += 1
            except KeyError:
                self.bq_counts[bq] = 1

        del self.bqs

    def dump(self):
        """Dump a pickle file containing base quality score count
        """

        with open(self.out_file, "wb") as f:
            pickle.dump(self.bq_counts, f)

    def run(self):
        """Run this ChromBqScanner to scan base qualities of given chromosomes
        and dump base quality score count
        """

        self.init_scanner()
        self.scan()
        self.postprocess()
        self.dump()


class BamBqScanner():
    """Scan base qualities of reads aligned to bam file to generate
    base quality distribution. Generated base quality distributions
    are used for creating lookup table.
    """

    buffer_ratio = 0.05
    total_base_num_to_scan = 10 ** 7

    def __init__(self, input_bam_file: str, ref_file: str, chromosomes: List[str],
                 sample_name: str, output_dir: str, process_num: int):
        """Create Scanner to scan base qualities of an input bam file

        param input_bam_file: path of BAM file for base quality scanning
        param ref_file: reference FASTA file
        param chromosomes: target chromosomes to restrict scanning
        garam sample_name: the name which will be used for output file prefix
        param output_dir: directory path which holds output file
        param process_num: # of process number used for the scanning.
                           However, process number for the actual scanning
                           would be number of target chromosomes at most.
        """

        self.input_bam_file = input_bam_file
        self.ref_file = ref_file
        self.chromosomes = chromosomes
        self.sample_name = sample_name
        self.output_dir = output_dir

        self.process_num = min(process_num, len(self.chromosomes))

        self.chrom_length = {}
        self.chrom_base_num = {}

        self.intermediate_files = {}

        self.bq_counts = {}
        self.bq_dist = {}

    def correct_jobs(self):
        """Fill up scanning base number to 'total_base_num_to_scan',
        from the longest chromosome
        """

        base_num_remained = self.total_base_num_to_scan \
                            - sum(self.chrom_base_num.values())
        if base_num_remained == 0: return
        for chromosome in self.chromosomes:
            self.chrom_base_num[chromosome] += 1
            base_num_remained -= 1
            if base_num_remained == 0: break

    def get_chromosome_len(self):
        """Get chromosome length from bam header and check whether every input
        chromosome length is normally got from the header
        """

        input_bam_handle = pysam.AlignmentFile(self.input_bam_file, "rb")
        bam_chrom_lengths = input_bam_handle.lengths

        for chrom_idx in range(0, len(bam_chrom_lengths)):
            chrom = input_bam_handle.get_reference_name(chrom_idx)
            if chrom in self.chromosomes:
                self.chrom_length[chrom] = bam_chrom_lengths[chrom_idx]

        self.chromosomes = [c[0] for c in sorted(self.chrom_length.items(),
                            key=operator.itemgetter(1), reverse=True)]

        assert len(self.chromosomes) == len(self.chrom_length)

        input_bam_handle.close()

    def get_bam_stats(self, chromosome):
        """Calculate total read count and average read length of an input bam
        file for an input chromosome

        param chromosome: target chromosome for statistics
        """

        input_bam_handle = pysam.AlignmentFile(self.input_bam_file, "rb")

        read_count = 0
        sum_read_length = 0
        for read in input_bam_handle.fetch(chromosome):
            read_count += 1
            try:
                sum_read_length += read.reference_length
            except TypeError:
                continue

        #assert (read_count > 0) #TODO - 22 Nov removed

        input_bam_handle.close()

        if read_count == 0: return 0, 0 # TODO - 22 Nov added
        return read_count, sum_read_length // read_count

    def split_jobs_by_chromosome(self):
        """Distribute base numbers to scan for each chromosome depending on
        chromosome length and set intermediate file names for chromosomes
        """

        self.get_chromosome_len()

        sum_chrom_length = sum(self.chrom_length.values())
        for chrom in self.chromosomes:
            self.intermediate_files[chrom] \
                = os.path.join(self.output_dir,
                               f'{self.sample_name}_bqcount.{chrom}.pkl')
            self.chrom_base_num[chrom] \
                = math.floor(self.total_base_num_to_scan *
                             self.chrom_length[chrom]/float(sum_chrom_length))

        self.correct_jobs()

    def scan_by_chromosome(self, chromosome):
        """Scan a chromosome to generate base quality distribution

        param chromosome: target chromosome for scanning
        """

        total_read_count, avg_read_length = self.get_bam_stats(chromosome)
        if total_read_count == 0: return # TODO : 22 Nov added
        base_num_to_scan = self.chrom_base_num[chromosome]
        read_count_per_batch = round(base_num_to_scan / avg_read_length
                                     * (1+self.buffer_ratio))
        ChromBqScanner(self.input_bam_file, chromosome,
                       self.intermediate_files[chromosome],
                       total_read_count, base_num_to_scan,
                       read_count_per_batch).run()

    def merge_intermediates(self):
        """Merge intermediate files results and calculate each base quality
        score proportion
        """

        for chromosome in self.chromosomes:
            if not os.path.isfile(self.intermediate_files[chromosome]): continue    #TODO - 22 Nov added

            with open(self.intermediate_files[chromosome], "rb") as f:
                chrom_bq_counts = pickle.load(f)

            for bq in chrom_bq_counts:
                try:
                    self.bq_counts[bq] += chrom_bq_counts[bq]
                except KeyError:
                    self.bq_counts[bq] = chrom_bq_counts[bq]

        #assert sum(self.bq_counts.values()) == self.total_base_num_to_scan #TODO - 22 Nov modified

        sum_bq_counts = sum(self.bq_counts.values())
        for bq in self.bq_counts:
            self.bq_dist[bq] = self.bq_counts[bq] / float(sum_bq_counts)

    def postprocess(self):
        """Fill up the missing base quality within range from the lowest
        base quality to the highest base quality from scanning
        """

        bqs = self.bq_dist.keys()
        for bq in range(min(bqs), max(bqs)):
            if bq not in bqs:
                self.bq_dist[bq] = 0.0

    def dump(self):
        """Dump a pickle file containing base quality distribution
        """

        bq_dist_path = os.path.join(self.output_dir,
                                    f'{self.sample_name}_dist.pkl')
        with open(bq_dist_path, "wb") as f:
            pickle.dump(self.bq_dist, f)

    def remove_intermediates(self):
        """Remove intermediate files
        """

        for intermediate_file in self.intermediate_files.values():
            if not os.path.isfile(intermediate_file): continue    #TODO - 22 Nov added
            os.remove(intermediate_file)

    def scan(self):
        """Scan bam files by each chromosome to generate base quality
        distribution, multiprocessing is used if process_num > 1
        """

        self.split_jobs_by_chromosome()
        os.makedirs(self.output_dir, exist_ok=True)
        with Pool(self.process_num) as pool:
            pool.map(self.scan_by_chromosome, self.chromosomes)

    def generate_bq_distribution(self):
        """Generate base quality distribution by merging individual chromosome
        scanning results (base quality score count files)
        """

        self.merge_intermediates()
        self.postprocess()
        self.dump()
        self.remove_intermediates()

    def run(self):
        """Run this BamBqScanner to scan input bam file and generate
        base quality distribution of it
        """

        self.scan()
        self.generate_bq_distribution()


class LookupTableCreator():
    """Create a lookup table for base qualities of input base quality
    distribution based on reference base quality distribution.
    Output values for input base qualities are matched reference base qualities
    and probability score. The created lookup table is used in BamBqNormalizer.
    """

    def __init__(self, ref_bq_dist_path: str, input_bq_dist_path: str,
                 sample_name: str, output_dir: str):

        """Lookup table creator based on reference and input base quality
        distribution

        param ref_bq_dist_path: path of reference base quality distribution,
                                base quality distribution should be from
                                BamBqScanner
                                {RefBq1 : Prop1, RefBq2 : Prop2, ...}
        param input_bq_dist_path: path of input base quality distribution
                                  {InBq1 : Prop1, InBq2 : Prop2, ...}
        param sample_name: the name which will be used for output file prefix
        param output_dir: path of directory where lookup table pickled object
                          file will be dumped
        """

        self.ref_bq_dist_path = ref_bq_dist_path
        self.input_bq_dist_path = input_bq_dist_path
        self.sample_name = sample_name
        self.output_dir = output_dir

        self.ref_bq_dist = {}
        self.input_bq_dist = {}
        self.ref_bqs = []
        self.input_bqs = []

        self.bq_correspondence = {}
        self.lookup_table = {}

    def check_bq_distribution(self):
        """Check whether base quality distribution files are truncated
        """

        ref_bq_proportions = self.ref_bq_dist.values()
        assert min(ref_bq_proportions) >= 0
        assert round(sum(ref_bq_proportions), FLOATING_POINT_DIGIT_LIMIT-2) == 1.0
        assert len(self.ref_bqs) == max(self.ref_bqs)-min(self.ref_bqs)+1

        input_bq_proportions = self.input_bq_dist.values()
        assert min(input_bq_proportions) >= 0
        assert round(sum(input_bq_proportions), FLOATING_POINT_DIGIT_LIMIT-2) == 1.0
        assert len(self.input_bqs) == max(self.input_bqs)-min(self.input_bqs)+1

    def load_bq_distribution(self):
        """Load reference and input base quality distribution
        """

        with open(self.ref_bq_dist_path, "rb") as f1:
            self.ref_bq_dist = pickle.load(f1)
        with open(self.input_bq_dist_path, "rb") as f2:
            self.input_bq_dist = pickle.load(f2)

        self.ref_bqs = sorted(self.ref_bq_dist.keys(), reverse=True)
        self.input_bqs = sorted(self.input_bq_dist.keys(), reverse=True)

        self.check_bq_distribution()

    def correspond_to_ref_bq_distribution(self, input_bq_proportion_start,
                                          input_bq_proportion_end):
        """Correspond reference base qualities to a input base quality
        """

        ref_bq_correspondence = {}

        ref_bq_proportion_end = 0.0
        for ref_bq in self.ref_bqs:
            ref_bq_proportion_start = ref_bq_proportion_end
            ref_bq_proportion_end = ref_bq_proportion_start + self.ref_bq_dist[ref_bq]

            if ref_bq_proportion_end <= input_bq_proportion_start: continue
            if input_bq_proportion_end <= ref_bq_proportion_start: break

            if input_bq_proportion_end < ref_bq_proportion_end:
                if ref_bq_proportion_start < input_bq_proportion_start:
                    proportion = input_bq_proportion_end - input_bq_proportion_start
                else:
                    proportion = input_bq_proportion_end - ref_bq_proportion_start
            else:
                if ref_bq_proportion_start < input_bq_proportion_start:
                    proportion = ref_bq_proportion_end - input_bq_proportion_start
                elif input_bq_proportion_start <= ref_bq_proportion_start:
                    proportion = ref_bq_proportion_end - ref_bq_proportion_start

            ref_bq_correspondence[ref_bq] = proportion

        return ref_bq_correspondence

    def correspond_bq_distribution(self):
        """Correspond reference and input base qualities based on proportion
        """

        input_bq_proportion_end = 0.0
        for input_bq in self.input_bqs:
            input_bq_proportion_start = input_bq_proportion_end
            input_bq_proportion_end = input_bq_proportion_start \
                                      + self.input_bq_dist[input_bq]
            self.bq_correspondence[input_bq] = \
                self.correspond_to_ref_bq_distribution(input_bq_proportion_start,
                                                       input_bq_proportion_end)

    def set_probability(self):
        """Set the probability of matched (reference) base quality for each
        input base quality based on weight of proportion
        """

        for input_bq, ref_bq_correspondence in self.bq_correspondence.items():
            proportion_sum = sum(ref_bq_correspondence.values())
            for ref_bq, proportion in ref_bq_correspondence.items():
                try:
                    prob = proportion / proportion_sum
                except ZeroDivisionError:
                    assert (len(ref_bq_correspondence) == 1)
                    prob = 1.0

                try:
                    self.lookup_table[input_bq][ref_bq] = prob
                except KeyError:
                    self.lookup_table[input_bq] = {ref_bq : prob}

    def dump_lookup_table(self):
        """Dump a pickle file containing lookup table
        """

        lookup_table_path = \
            os.path.join(self.output_dir, f'{self.sample_name}_lookup.pkl')
        with open(lookup_table_path, "wb") as f:
            pickle.dump(self.lookup_table, f)

    def run(self):
        """Run this LookupTableCreator to load base quality distribution files
        and create a lookup table
        """

        self.load_bq_distribution()
        self.correspond_bq_distribution()
        self.set_probability()
        self.dump_lookup_table()


class BamBqNormalizer():
    """bq-normalize the BAM file : normalize base qualities of BAM file reads,
     needed for BAM preprocessing for AIVariant. bq-normalization is likely
     quantile normalization method. Difference from quantile normalization is
     how to handling ties. For handling ties, random selection base on
     probability (or weight) is used, rather than average.
    """

    bq_lower_cap = 0
    bq_upper_cap = 93
    linux_ulimit = 1024

    def __init__(self, input_bam_file: str, lookup_table_path: str,
                 process_num: int, out_bam_dir: str, out_bam_prefix: str,
                 chromosome=None):
        """bq-normalizer of bam file based on lookup table

        param input_bam_file: path of bam file for bq-normalization
        param lookup_table_path_: file path of lookup table made from LookupTableCreator.
                                 Key : base quality of input_bam_file,
                                 Value : sub dictionaries containing probability
                                 information for base quality to be converted.
                                 Sum of probabilities of each sub dictionary should be 1.0.
                                 {InBq1 : {ConvBq1 : Prob1, ConvBq2 : Prob2, ConvBq3 : Prob3},
                                  InBq2 : {ConvBq3 : Prob4, ConvBq4 : Prob5}, ... }
        param process_num: # of process to for bq conversion.
        param out_bam_dir: directory path which holds output bq converted bam file
        param out_bam_prefix: prefix of output bq converted BAM file,
        param chromosome: chromosome [str] (set None automatically if not specified)
        """

        self.input_bam_file = input_bam_file
        self.lookup_table_path = lookup_table_path
        self.process_num = process_num
        self.out_bam_dir = out_bam_dir
        self.out_bam_prefix = out_bam_prefix
        self.chromosome = chromosome

        self.lookup_table = {}
        self.read_count = None
        self.read_batch_size = 10 ** 6
        self.read_batches = []
        self.sub_bam_num = None
        self.sub_bam_dir = os.path.join(self.out_bam_dir,
                                        f'tmp_{self.out_bam_prefix}')
        self.out_bam_file = os.path.join(self.out_bam_dir,
                                         f'{self.out_bam_prefix}.bam')

    def count_reads(self):
        """Count the read number of an input bam file
        """

        input_bam_handle = pysam.AlignmentFile(self.input_bam_file, "rb")
        self.read_count = input_bam_handle.count(self.chromosome)
        input_bam_handle.close()

    def generate_read_batches(self):
        """Split reads of input bam file to read batches based on size of a
        read batch
        """

        sub_bam_num = self.read_count//self.read_batch_size + 1
        if sub_bam_num > self.linux_ulimit:
            self.read_batch_size = self.read_count//self.linux_ulimit + 1
        self.sub_bam_num = self.read_count//self.read_batch_size + 1

        for read_start_idx in range(0, self.read_count, self.read_batch_size):
            read_end_idx = read_start_idx + self.read_batch_size
            if read_end_idx > self.read_count:
                read_end_idx = self.read_count
            self.read_batches.append((read_start_idx, read_end_idx))

    def count_bq(self, read_idx_start, read_idx_end):
        """Count the number of base qualities for input read index range

        param read_idx_start: 0-based, open
        param read_idx_end: 0-based, closed
        """

        input_bam_handle = pysam.AlignmentFile(self.input_bam_file, "rb")
        bam_itr = input_bam_handle.fetch(self.chromosome, multiple_iterators=True)

        reads = []
        bq_counts = {}

        for read_idx, read in enumerate(bam_itr):
            if read_idx < read_idx_start: continue
            elif read_idx >= read_idx_end: break
            else:
                for bq in read.query_qualities:
                    if bq is not None:
                        try:
                            bq_counts[bq] += 1
                        except KeyError:
                            bq_counts[bq] = 1
                reads.append(read)

        del bam_itr

        input_bam_handle.close()

        return reads, bq_counts

    def random_convert_reads_bq(self, read_idx_start, read_idx_end):
        """Convert original base qualities based on lookup table, if there are
        more than one matched base qualities, randomly choose a base quality
        based on probability written in lookup table
        """
        reads, bq_counts = self.count_bq(read_idx_start, read_idx_end)

        converted_bqs_4_input_bam_bq = {}
        for bq in bq_counts:
            try:
                matched_bq_probs = self.lookup_table[bq]
            except KeyError:
                assert self.bq_lower_cap <= bq <= self.bq_upper_cap
                raise Exception (f"[Error] bq {bq} does not exists in the lookup table, "
                                 f"means bq {bq} is within unexpected bq range.")

            matched_bqs = list(matched_bq_probs.keys())
            probs = list(matched_bq_probs.values())
            bq_count = bq_counts[bq]

            converted_bqs_4_input_bam_bq[bq] \
                = random.choices(matched_bqs, k=bq_count, weights=probs)

        bq_converted_reads = []
        for read in reads:
            read_bqs = read.query_qualities

            for bq_idx, prior_bq in enumerate(read_bqs):
                if prior_bq is not None:
                    read_bqs[bq_idx] = converted_bqs_4_input_bam_bq[prior_bq].pop()
            read.query_qualities = read_bqs

            bq_converted_reads.append(read)

        del reads
        del bq_counts
        del converted_bqs_4_input_bam_bq

        return bq_converted_reads

    def write_sub_bam(self, reads, bam_suffix):
        """Write bq-normalized reads to a sub bam file
        """

        input_bam_handle = pysam.AlignmentFile(self.input_bam_file, "rb")
        sub_bam_file = os.path.join(self.sub_bam_dir,
                                    f'{self.out_bam_prefix}.{bam_suffix}.bam')
        out_bam_handle = pysam.AlignmentFile(sub_bam_file, "wb",
                                             template = input_bam_handle)

        for read in reads: out_bam_handle.write(read)

        input_bam_handle.close()
        out_bam_handle.close()

    def convert_bq_by_batch(self, read_batch, sub_bam_idx):
        """Convert base qualities of read batch based on lookup table

        param read_batch: range of read batch, 0-based and half open,
                          [read_start_idx, read_end_idx) tuple
        param sub_bam_idx: index of a sub bam file
        """

        bq_converted_reads = self.random_convert_reads_bq(read_batch[0], read_batch[1])
        self.write_sub_bam(bq_converted_reads, sub_bam_idx)

    def preprocess_lookup_table(self):
        """Check whether the input lookup table is truncated.
        Next, add additional(and possible) base qualities which could not be
        detected in creating lookup table step due to subsampling.
        """

        lookup_table_handle = open(self.lookup_table_path, "rb")
        self.lookup_table = pickle.load(lookup_table_handle)
        lookup_table_handle.close()

        input_bam_bqs = self.lookup_table.keys()
        max_input_bam_bq = max(input_bam_bqs)
        min_input_bam_bq = min(input_bam_bqs)
        assert (max_input_bam_bq <= self.bq_upper_cap
                and min_input_bam_bq >= self.bq_lower_cap)

        for bq in range(min_input_bam_bq, max_input_bam_bq+1):
            try:
                matched_bq_probs = self.lookup_table[bq]
                assert round(sum(matched_bq_probs.values()), FLOATING_POINT_DIGIT_LIMIT-2) == 1.0
            except KeyError:
                raise Exception(f"[Error] bq {bq} does not exist in the lookup table, "
                                f"normally means the lookup table is truncated.")

        max_out_bam_bq = \
            max(max(matched_bq_probs.keys()) for matched_bq_probs in self.lookup_table.values())
        for bq in range(max_input_bam_bq+1, self.bq_upper_cap+1):
            self.lookup_table[bq] = {max_out_bam_bq : 1.0}
        min_out_bam_bq = \
            min(min(matched_bq_probs.keys()) for matched_bq_probs in self.lookup_table.values())
        for bq in range(self.bq_lower_cap, min_input_bam_bq):
            self.lookup_table[bq] = {min_out_bam_bq : 1.0}

    def convert_bq_sub_bam(self):
        """Convert base qualities of an input bam file based on lookup table,
        multiprocessing is used if process_num > 1
        """

        self.count_reads()
        self.generate_read_batches()
        os.makedirs(self.sub_bam_dir, exist_ok = True)
        with Pool(self.process_num) as pool:
            pool.starmap(self.convert_bq_by_batch,
                         [(batch, idx) for idx, batch in enumerate(self.read_batches)])

    def merge_bam(self):
        """Merge bq-normalized sub bam files to a merged bam file
        """

        sub_bam_paths = [os.path.join(self.sub_bam_dir,
                                      f'{self.out_bam_prefix}.{out_bam_suffix}.bam')
                         for out_bam_suffix in range(self.sub_bam_num)]

        if len(sub_bam_paths) == 1:
            shutil.move(sub_bam_paths[0], self.out_bam_file)
            BamMerger([], self.process_num , self.out_bam_file).index()
        else:
            BamMerger(sub_bam_paths, self.process_num, self.out_bam_file).run()

        shutil.rmtree(self.sub_bam_dir)

    def run(self):
        """Run this BamBqNormalizer to bq-normalize an input bam files
        """

        self.preprocess_lookup_table()
        self.convert_bq_sub_bam()
        self.merge_bam()


class BamMerger():
    """Merge sub bam files to a bam file
    """

    def __init__(self, sub_bam_files: List[str], thread_num: int,
                 out_bam_file: str):
        """Merge bam files using 'Samtools'

        param sub_bam_files: list of sub bam file paths to merge
        param thread_n: # of thread for 'Samtools merge'
        param out_bam_file: output bam file path (a merged bam)
        """

        self.sub_bam_files = sub_bam_files
        self.out_bam_file = out_bam_file
        self.thread_num = thread_num

    def index(self):
        """Index a merged bam file
        """

        cmd = f'{samtools} index -@ {self.thread_num} {self.out_bam_file}'
        os.system(cmd)

    def merge(self):
        """Merge bam files
        """

        input_bam_files = ' '.join(self.sub_bam_files)
        cmd = f'{samtools} merge -c -p --threads {self.thread_num} ' \
              f'{self.out_bam_file} {input_bam_files}'
        os.system(cmd)

    def run(self):
        """Run this BamMerger to merge input bam files to a merged bam file
        """

        self.merge()
        self.index()