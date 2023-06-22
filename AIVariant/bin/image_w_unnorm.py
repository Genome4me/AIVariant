###
# Written by Hyeonseong Jeon
# email: jun0605@snu.ac.kr
###

from typing import *
import os
import sys
import pickle
from pathlib import Path
import pysam
import numpy as np
import pyBigWig

from utils import Utils
from image_utils_w_unnorm import ReadFilter, ImageChannel, ImageNormalizer, \
    StatImageChannel, ReadsStat


class ImageGenerator():
    """Generate images from bam file
    """

    def __init__(self, tumor_bam_file: str, normal_bam_file: str,
                 ref_file: str,
                 chromosome: str, position: int,
                 reference: str, alternate: str,
                 tumor_unnorm_bam_file: str, normal_unnorm_bam_file: str,
                 output_path: str,
                 hg_version: str,
                 epi_features=None,
                 epi_path: str = None,
                 target_depth: int = 30,
                 img_height: int = 100):
        """Construct image generator

        param tumor_bam_file: BAM file from tumor sample
        param normal_bam_file: BAM file from normal sample
        param ref_file: reference FASTA file

        param chromosome: chromosome
        param position: position on the genome, 0-based
        param reference: reference base of a given position on a given chromosome
        param alternate: alternate base of a given position on a given chromosome
        param epi_features: epigenetic features to use
        param epi_path: path where it holds the files for epigenetic features

        param output_path: the directory in which image files will be stored

        param tumor_unnorm_bam_file: BAM file from tumor sample with unnormalized bqs
        param normal_unnorm_bam_file: BAM file from normal sample with unnormalized bqs
        """

        self.IMG_HEIGHT = img_height
        self.target_depth = target_depth

        self.unnorm_bam_files = [tumor_unnorm_bam_file, normal_unnorm_bam_file]

        self.unnorm_hash = dict()

        self.tumor_bam_file = tumor_bam_file
        self.normal_bam_file = normal_bam_file
        self.ref_file = ref_file

        self.chromosome = chromosome
        if self.chromosome.startswith('chrchr'):
            self.contig = self.chromosome[3:]
        else:
            self.contig = self.chromosome

        self.reference = reference
        self.position = position

        if hg_version == 'hg19':
            self.hg_version = 'hg19'
            self.hg19_position = position
        elif hg_version == 'hg38':
            self.hg_version = 'hg38'
            try:
                self.hg19_position = get_lifter('hg38', 'hg19')[self.contig][self.position][0][1]
            except IndexError:
                self.hg19_position = None
        else:
            sys.exit(f"{hg_version} not supported : hg19 or hg38 required")

        self.alternate = alternate

        self.tumor_raw_reads = []
        self.tumor_processed_reads = []
        self.normal_raw_reads = []
        self.normal_processed_reads = []

        self.tumor_unnorm_reads = []
        self.tumor_processed_unnorm_reads = []
        self.normal_unnorm_reads = []
        self.normal_processed_unnorm_reads = []

        # [(col_idx, read_pos, cigar_char)]
        self.tumor_raw_reads_coordinate = []
        self.normal_raw_reads_coordinate = []
        self.tumor_processed_reads_coordinate = []
        self.normal_processed_reads_coordinate = []

        self.coordinate = dict()  # (reference_position, indel_idx) => matrix column index
        self.insertion_size = dict()  # reference_position => maximum insertion size after the reference position
        self.img_shape = (self.IMG_HEIGHT, 101, 16)

        self.stat_img_shape = (5, 101, 9)

        self.tumor_raw_img = None
        self.tumor_processed_img = None
        self.normal_raw_img = None
        self.normal_processed_img = None

        self.tumor_processed_stat_img = None
        self.normal_processed_stat_img = None

        self.img_dict = {
            "tumor_raw_img": None,
            "tumor_processed_img": None,
            "normal_raw_img": None,
            "normal_processed_img": None,
            "tumor_processed_stat_img": None,
            "normal_processed_stat_img": None,
            "lodt": None,
            "lodn": None,
            "unnorm_lodt": None,
            "unnorm_lodn": None,
            "label": None
        }

        self.image_keys = ['tumor_raw_img', 'tumor_processed_img',
                           'normal_raw_img', 'normal_processed_img']

        if epi_features is None:
            self.epi_features = ['H3K9me3']
        else:
            self.epi_features = epi_features

        self.epi_path = epi_path
        self.epi_pval_path = Path(os.path.join(self.epi_path, "pvals"))
        self.epi_fold_path = Path(os.path.join(self.epi_path, "folds"))

        self.output_path = output_path
        self.output_file = os.path.join(output_path,
                                        f"{chromosome}:{position}:{reference}:{alternate}.pkl")

        self.depth_limit = 0

    def gen_unnorm_hash(self):
        """Helper function which fetch reads of unnormalized base qualities from bam files and hashing them
        """

        for bam_file in self.unnorm_bam_files:
            with pysam.AlignmentFile(bam_file, "rb") as bam_handle:
                for pileupcolumn in bam_handle.pileup(self.chromosome[3:],
                                                      self.position,
                                                      self.position + 1,
                                                      truncate=True,
                                                      stepper="all",
                                                      max_depth=self.depth_limit):
                    for read in pileupcolumn.pileups:
                        key = (
                            read.alignment.query_name, read.alignment.cigarstring,
                            read.alignment.query_sequence)
                        # assert key not in self.unnorm_hash
                        self.unnorm_hash[key] = read

    def order(self, read):
        """Order of reads on the image,
        Firstly sorted in order with max alternate base, alternate base, reference base, other base(ex. deletion)
        Secondly sorted in order with mapping qualiy

        param read: aligned read (pysam pileupread)
        """

        return Utils.order(self.position, self.reference, self.alternate, read)

    def fetch_raw_reads(self, normal=False):
        """Fetch raw (unprocessed except preprocess) reads that covers the given position from BAM file

        param nomral: indicate whether the BAM file is from tumor sample or normal sample
                      (True, False)
        """

        if not normal:
            bam_handle = pysam.AlignmentFile(self.tumor_bam_file, "rb")
        else:
            bam_handle = pysam.AlignmentFile(self.normal_bam_file, "rb")

        for pileupcolumn in bam_handle.pileup(self.chromosome[3:],
                                              self.position, self.position + 1,
                                              truncate=True, stepper="all",
                                              max_depth=self.depth_limit):
            reads = Utils.fetch_reads_from_pileupcolumn(pileupcolumn)

        reads.sort(key=self.order, reverse=True)

        if not normal:
            self.tumor_raw_reads = reads[:]
        else:
            self.normal_raw_reads = reads[:]
        del reads

        bam_handle.close()

    def fetch_processed_reads(self, normal=False):
        """Fetch processed reads that cover the given position from BAM file

        param normal: indicate whether fetch from tumor sample BAM file or normal sample BAM file
                      (True, False)
        """

        if not normal:
            raw_reads = self.tumor_raw_reads
            processed_reads = self.tumor_processed_reads
            filters = [ReadFilter.delete_filter,
                       ReadFilter.base_quality_filter,
                       ReadFilter.mapping_quality_filter,
                       ReadFilter.softclip_filter,
                       lambda read: ReadFilter.mismatch_filter(self.ref_file,
                                                               self.chromosome,
                                                               read)]
        else:
            raw_reads = self.normal_raw_reads
            processed_reads = self.normal_processed_reads
            filters = [ReadFilter.delete_filter,
                       ReadFilter.base_quality_filter]

        for read in raw_reads:
            if all([_filter(read) for _filter in
                    filters]): processed_reads.append(read)


    def fetch_unnorm_reads(self, normal=False, processed=False):
        """Fetch reads with unnormalized base qualities from unnormalized BAM file

        param normal: indicate whether the BAM file is from tumor sample or normal sample
                      (True, False)
        """

        if not normal:
            if not processed:
                fetched_reads = self.tumor_raw_reads
            else:
                fetched_reads = self.tumor_processed_reads
        else:
            if not processed:
                fetched_reads = self.normal_raw_reads
            else:
                fetched_reads = self.normal_processed_reads

        matched_reads = []
        for fr in fetched_reads:
            key = (fr.alignment.query_name, fr.alignment.cigarstring,
                   fr.alignment.query_sequence)
            matched_reads.append(self.unnorm_hash[key])

        if not normal:
            if not processed:
                self.tumor_unnorm_reads = matched_reads
            else:
                self.tumor_processed_unnorm_reads = matched_reads
        else:
            if not processed:
                self.normal_unnorm_reads = matched_reads
            else:
                self.normal_processed_unnorm_reads = matched_reads


    def _get_insertion_size(self, reads):
        for read in reads:
            cigar_operations = [pair[0] for pair in read.alignment.cigartuples for
                                _ in range(pair[1])]
            aligned_pairs = read.alignment.get_aligned_pairs()

            assert len(cigar_operations) == len(aligned_pairs)

            for i, p in enumerate(aligned_pairs):
                read_pos, ref_pos = p
                if ref_pos is None:
                    continue
                if ref_pos not in self.insertion_size:
                    self.insertion_size[ref_pos] = 0

                I, I_i = 0, i + 1
                # D, D_i = 0, i+1
                while I_i < len(aligned_pairs) and aligned_pairs[I_i][
                    1] is None and cigar_operations[I_i] == 1:
                    I_i += 1
                    I += 1
                # while D_i < len(aligned_pairs) and aligned_pairs[D_i][0] is None and cigar_operations[D_i] == 2:
                #     D_i += 1
                #     D += 1

                self.insertion_size[ref_pos] = max(self.insertion_size[ref_pos], I)


    def get_coordinate(self):
        self._get_insertion_size(self.tumor_raw_reads)
        self._get_insertion_size(self.normal_raw_reads)
        self._get_insertion_size(self.tumor_unnorm_reads)
        self._get_insertion_size(self.normal_unnorm_reads)

        positions = []
        for ref_pos in range(self.position - 50, self.position + 51):
            if ref_pos not in self.insertion_size:
                positions.append((ref_pos, 0))
            else:
                for i in range(self.insertion_size[ref_pos] + 1):
                    positions.append((ref_pos, i))

        pos_idx = positions.index((self.position, 0))
        positions = positions[pos_idx - 50:pos_idx + 51]

        for i, position in enumerate(positions):
            self.coordinate[position] = i


    def get_read_coordinate(self, read_coord, reads):
        for read in reads:
            read_coord.append([])
            aligned_pairs = read.alignment.get_aligned_pairs()
            cigar_chars = Utils.get_unrolled_cigar_chars(read)

            assert len(aligned_pairs) == len(cigar_chars)

            insertion_idx = 0
            prev_ref_pos = None
            for idx, pair in enumerate(aligned_pairs):
                read_pos, ref_pos = pair

                if read_pos is not None and ref_pos is not None:
                    insertion_idx = 0
                    prev_ref_pos = ref_pos
                    if (ref_pos, insertion_idx) in self.coordinate:
                        col_idx = self.coordinate[(ref_pos, insertion_idx)]
                        read_coord[-1].append(
                            (col_idx, read_pos, cigar_chars[idx]))
                elif read_pos is not None and ref_pos is None:
                    insertion_idx += 1
                    if (prev_ref_pos, insertion_idx) in self.coordinate:
                        col_idx = self.coordinate[(prev_ref_pos, insertion_idx)]
                        read_coord[-1].append(
                            (col_idx, read_pos, cigar_chars[idx]))
                elif read_pos is None and ref_pos is not None:
                    insertion_idx = 0
                    prev_ref_pos = ref_pos
                    if (ref_pos, insertion_idx) in self.coordinate:
                        col_idx = self.coordinate[(ref_pos, insertion_idx)]
                        read_coord[-1].append(
                            (col_idx, read_pos, cigar_chars[idx]))


    def fetch_reference_channel(self):
        """Fetch referece channel, reference channel index is [0:4)
        """

        ref_handle = pysam.FastaFile(self.ref_file)
        offset = (self.img_shape[1] - 1) // 2
        ref_seq = ref_handle.fetch(self.chromosome[3:], self.position - offset,
                                   self.position + offset + 1).upper()  # 220802 fixed
        ref_handle.close()

        start_ref_pos = self.position - offset

        ImageChannel.fetch_reference_channel(
            min(self.IMG_HEIGHT, len(self.tumor_raw_reads)), self.tumor_raw_img,
            ref_seq, start_ref_pos, self.coordinate)
        ImageChannel.fetch_reference_channel(
            min(self.IMG_HEIGHT, len(self.tumor_processed_reads)),
            self.tumor_processed_img, ref_seq, start_ref_pos, self.coordinate)
        ImageChannel.fetch_reference_channel(
            min(self.IMG_HEIGHT, len(self.normal_raw_reads)), self.normal_raw_img,
            ref_seq, start_ref_pos, self.coordinate)
        ImageChannel.fetch_reference_channel(
            min(self.IMG_HEIGHT, len(self.normal_processed_reads)),
            self.normal_processed_img, ref_seq, start_ref_pos, self.coordinate)


    def fetch_stat_reference_channel(self):
        """Fetch reference channel for stat img, channel index is 0
        """

        channel_idx = 0

        ref_handle = pysam.FastaFile(self.ref_file)
        offset = (self.img_shape[1] - 1) // 2
        ref_seq = ref_handle.fetch(self.chromosome[3:], self.position - offset,
                                   self.position + offset + 1)
        ref_handle.close()

        start_ref_pos = self.position - offset
        StatImageChannel.fetch_reference_channel(
            self.tumor_processed_stat_img[:, :, channel_idx], ref_seq,
            start_ref_pos, self.coordinate)
        StatImageChannel.fetch_reference_channel(
            self.normal_processed_stat_img[:, :, channel_idx], ref_seq,
            start_ref_pos, self.coordinate)


    def fetch_read_channels(self, img, reads, unnorm_reads, read_coord):
        """Fetch read channels for img
        """

        # nucleotide channel
        ImageChannel.fetch_nucleotide_channel(img[:, :, 4:8], reads, read_coord)

        # mapping quality channel
        ImageChannel.fetch_mq_channel(img[:, :, 8], reads, read_coord)

        # strand channel
        ImageChannel.fetch_strand_channel(img[:, :, 9], reads, read_coord)

        # base quality channel
        ImageChannel.fetch_bq_channel(img[:, :, 10], reads, read_coord)

        # soft clip channel
        ImageChannel.fetch_softclip_channel(img[:, :, 11], reads, read_coord)

        # distance channel
        ImageChannel.fetch_distance_channel(img[:, :, 12], reads, read_coord)

        # unnorm base quality channel
        ImageChannel.fetch_bq_channel(img[:, :, 13], unnorm_reads, read_coord)


    def fetch_stat_read_channels(self):
        """Fetch read channels for stat img
        """

        StatImageChannel.fetch_frequency_channel(
            self.tumor_processed_stat_img[:, :, 1], self.tumor_processed_reads,
            self.tumor_processed_reads_coordinate)
        StatImageChannel.fetch_frequency_channel(
            self.normal_processed_stat_img[:, :, 1], self.normal_processed_reads,
            self.normal_processed_reads_coordinate)

        StatImageChannel.fetch_mq_channel(self.tumor_processed_stat_img[:, :, 2],
                                          self.tumor_processed_reads,
                                          self.tumor_processed_reads_coordinate)
        StatImageChannel.fetch_mq_channel(self.normal_processed_stat_img[:, :, 2],
                                          self.normal_processed_reads,
                                          self.normal_processed_reads_coordinate)

        StatImageChannel.fetch_strand_channel(
            self.tumor_processed_stat_img[:, :, 3], self.tumor_processed_reads,
            self.tumor_processed_reads_coordinate)
        StatImageChannel.fetch_strand_channel(
            self.normal_processed_stat_img[:, :, 3], self.normal_processed_reads,
            self.normal_processed_reads_coordinate)

        StatImageChannel.fetch_bq_channel(self.tumor_processed_stat_img[:, :, 4],
                                          self.tumor_processed_reads,
                                          self.tumor_processed_reads_coordinate)
        StatImageChannel.fetch_bq_channel(self.normal_processed_stat_img[:, :, 4],
                                          self.normal_processed_reads,
                                          self.normal_processed_reads_coordinate)

        StatImageChannel.fetch_distance_channel(
            self.tumor_processed_stat_img[:, :, 5], self.tumor_processed_reads,
            self.tumor_processed_reads_coordinate)
        StatImageChannel.fetch_distance_channel(
            self.normal_processed_stat_img[:, :, 5], self.normal_processed_reads,
            self.normal_processed_reads_coordinate)

        StatImageChannel.fetch_bq_channel(self.tumor_processed_stat_img[:, :, 6],
                                          self.tumor_processed_unnorm_reads,
                                          self.tumor_processed_reads_coordinate)
        StatImageChannel.fetch_bq_channel(self.normal_processed_stat_img[:, :, 6],
                                          self.normal_processed_unnorm_reads,
                                          self.normal_processed_reads_coordinate)


    def fetch_img_channels(self):
        """Fetch channels for each image
        """

        self.fetch_reference_channel()

        self.fetch_read_channels(self.tumor_raw_img,
                                 self.tumor_raw_reads[:self.IMG_HEIGHT],
                                 self.tumor_unnorm_reads[:self.IMG_HEIGHT],
                                 self.tumor_raw_reads_coordinate[:self.IMG_HEIGHT])
        self.fetch_read_channels(self.tumor_processed_img,
                                 self.tumor_processed_reads[:self.IMG_HEIGHT],
                                 self.tumor_processed_unnorm_reads[
                                 :self.IMG_HEIGHT],
                                 self.tumor_processed_reads_coordinate[
                                 :self.IMG_HEIGHT])
        self.fetch_read_channels(self.normal_raw_img,
                                 self.normal_raw_reads[:self.IMG_HEIGHT],
                                 self.normal_unnorm_reads[:self.IMG_HEIGHT],
                                 self.normal_raw_reads_coordinate[
                                 :self.IMG_HEIGHT])
        self.fetch_read_channels(self.normal_processed_img,
                                 self.normal_processed_reads[:self.IMG_HEIGHT],
                                 self.normal_processed_unnorm_reads[
                                 :self.IMG_HEIGHT],
                                 self.normal_processed_reads_coordinate[
                                 :self.IMG_HEIGHT])


    def fetch_stat_img_channels(self):
        """Fetch channels for each stat image
        """

        self.fetch_stat_reference_channel()
        self.fetch_stat_read_channels()


    def init_reads(self):
        """Initialize each category of reads
        """

        self.tumor_raw_reads = []
        self.tumor_processed_reads = []
        self.normal_raw_reads = []
        self.normal_processed_reads = []


    def init_imgs(self):
        """Initialize each category of image
        """

        self.get_coordinate()
        self.get_read_coordinate(self.tumor_raw_reads_coordinate,
                                 self.tumor_raw_reads)
        self.get_read_coordinate(self.normal_raw_reads_coordinate,
                                 self.normal_raw_reads)
        self.get_read_coordinate(self.tumor_processed_reads_coordinate,
                                 self.tumor_processed_reads)
        self.get_read_coordinate(self.normal_processed_reads_coordinate,
                                 self.normal_processed_reads)

        self.tumor_raw_img = np.zeros(self.img_shape)
        self.tumor_processed_img = np.zeros(self.img_shape)
        self.normal_raw_img = np.zeros(self.img_shape)
        self.normal_processed_img = np.zeros(self.img_shape)
        self.tumor_processed_stat_img = np.zeros(self.stat_img_shape)
        self.normal_processed_stat_img = np.zeros(self.stat_img_shape)


    def generate_reads(self):
        """Generate each category of reads by initializing and fetching
        """

        self.init_reads()
        self.fetch_raw_reads()
        self.fetch_raw_reads(True)
        self.fetch_processed_reads()
        self.fetch_processed_reads(True)

        self.fetch_unnorm_reads(normal=False, processed=False)
        self.fetch_unnorm_reads(normal=False, processed=True)
        self.fetch_unnorm_reads(normal=True, processed=False)
        self.fetch_unnorm_reads(normal=True, processed=True)


    def generate_imgs(self):
        """Generate each category of image by initializing, fetching and normalizing
        """

        self.init_imgs()

        self.fetch_img_channels()
        ImageNormalizer.normalize_img(self.tumor_raw_img)
        ImageNormalizer.normalize_img(self.tumor_processed_img)
        ImageNormalizer.normalize_img(self.normal_raw_img)
        ImageNormalizer.normalize_img(self.normal_processed_img)
        self.img_dict["tumor_raw_img"] = self.tumor_raw_img
        self.img_dict["tumor_processed_img"] = self.tumor_processed_img
        self.img_dict["normal_raw_img"] = self.normal_raw_img
        self.img_dict["normal_processed_img"] = self.normal_processed_img

        self.fetch_stat_img_channels()
        self.img_dict["tumor_processed_stat_img"] = self.tumor_processed_stat_img
        self.img_dict["normal_processed_stat_img"] = self.normal_processed_stat_img


    def generate_lodt(self):
        """Generate LODT score for the given position
        """

        lodt = ReadsStat.calc_lodt(self.tumor_processed_reads, self.reference,
                                   self.alternate)
        lodt = ReadsStat.normalize_lodt(lodt)
        self.img_dict["lodt"] = lodt


    def generate_lodn(self):
        """Generator LODN score for the given position
        """

        lodn = ReadsStat.calc_lodn(self.normal_processed_reads, self.reference,
                                   self.alternate)
        lodn = ReadsStat.normalize_lodn(lodn)
        self.img_dict["lodn"] = lodn


    def generate_unnorm_lodt(self):
        """Generate LODT score for the given position with unnormalized bqs
        """

        unnorm_lodt = ReadsStat.calc_lodt(self.tumor_processed_unnorm_reads,
                                          self.reference, self.alternate)
        unnorm_lodt = ReadsStat.normalize_lodt(unnorm_lodt)
        self.img_dict["unnorm_lodt"] = unnorm_lodt


    def generate_unnorm_lodn(self):
        """Generate LODN score for the given position with unnormalized bqs
        """

        unnorm_lodn = ReadsStat.calc_lodn(self.normal_processed_unnorm_reads,
                                          self.reference, self.alternate)
        unnorm_lodn = ReadsStat.normalize_lodn(unnorm_lodn)
        self.img_dict["unnorm_lodn"] = unnorm_lodn


    def _get_epigenetic_features(self):
        bw_pval_handles = []
        bw_fold_handles = []
        for pval_file in self.epi_pval_path.iterdir():
            if any([filter_str in str(pval_file) for filter_str in
                    self.epi_features]):
                bw_pval_handles.append(pyBigWig.open(str(pval_file)))
        for fold_file in self.epi_fold_path.iterdir():
            if any([filter_str in str(fold_file) for filter_str in
                    self.epi_features]):
                bw_fold_handles.append(pyBigWig.open(str(fold_file)))

        pvals = []
        for bw_pval in bw_pval_handles:
            if self.hg19_position == None:  #221216 modified
                pval = np.empty(101)
                pval[:] = np.nan
            else:
                try:
                    pval = np.array(bw_pval.values(self.contig, self.hg19_position - 50, self.hg19_position + 51))
                except RuntimeError:
                    pval = np.empty(101)
                    pval[:] = np.nan

            pval = 10 ** (-1. * pval)
            for i, p in enumerate(pval):
                if np.isnan(p) or np.isinf(p):
                    # print(f"pval {pval[i]} results in {p}")
                    pval[i] = 0
                    # print(f"{p} set to {pval[i]}")
            # pval = pval.reshape(1, -1, 1).repeat(self.img_shape[0], axis=0)
            pvals.append(pval)
        folds = []
        for bw_fold in bw_fold_handles:
            if self.hg19_position == None:  #221216 modified
                fold_ = np.empty(101)
                fold_[:] = np.nan
            else:
                try:
                    fold_ = np.array(bw_fold.values(self.contig, self.hg19_position - 50, self.hg19_position + 51))
                except RuntimeError:
                    fold_ = np.empty(101)
                    fold_[:] = np.nan

            fold = np.log2(fold_ + 1) / 3  # should check for inf or NaN
            for i, f in enumerate(fold):
                if np.isnan(f) or np.isinf(f):
                    # print(f"fold {fold_[i]} results in {f}")
                    fold[i] = 0
                    # print(f"{f} set to {fold[i]}")
            fold[np.where(fold > 1)] = 1
            # fold = fold.reshape(1, -1, 1).repeat(self.img_shape[0], axis=0)
            folds.append(fold)

        # for image_key in self.image_keys:
        #     image = self.img_dict[image_key]
        #     for pval in pvals:
        #         image = np.concatenate([image, pval], axis=-1)
        #     for fold in folds:
        #         image = np.concatenate([image, fold], axis=-1)

        #     self.img_dict[image_key] = image

        return pvals, folds


    def fetch_epigenetic_channel(self):
        """Fetch epigenetic channel, epigenetic channel index is [14:16)
        """

        offset = (self.img_shape[1] - 1) // 2
        pvals, folds = self._get_epigenetic_features()

        start_ref_pos = self.position - offset

        ImageChannel.fetch_epigenetic_channel(
            min(self.IMG_HEIGHT, len(self.tumor_raw_reads)), self.tumor_raw_img,
            pvals, folds, start_ref_pos, self.coordinate)
        ImageChannel.fetch_epigenetic_channel(
            min(self.IMG_HEIGHT, len(self.tumor_processed_reads)),
            self.tumor_processed_img, pvals, folds, start_ref_pos, self.coordinate)
        ImageChannel.fetch_epigenetic_channel(
            min(self.IMG_HEIGHT, len(self.normal_raw_reads)), self.normal_raw_img,
            pvals, folds, start_ref_pos, self.coordinate)
        ImageChannel.fetch_epigenetic_channel(
            min(self.IMG_HEIGHT, len(self.normal_processed_reads)),
            self.normal_processed_img, pvals, folds, start_ref_pos,
            self.coordinate)


    def fetch_stat_epigenetic_channel(self):
        offset = (self.stat_img_shape[1] - 1) // 2
        pvals, folds = self._get_epigenetic_features()

        start_ref_pos = self.position - offset

        for pval in pvals:
            StatImageChannel.fetch_epigenetic_channel(
                self.tumor_processed_stat_img[:, :, 7], pval, start_ref_pos,
                self.coordinate)
            StatImageChannel.fetch_epigenetic_channel(
                self.normal_processed_stat_img[:, :, 7], pval, start_ref_pos,
                self.coordinate)
        for fold in folds:
            StatImageChannel.fetch_epigenetic_channel(
                self.tumor_processed_stat_img[:, :, 8], fold, start_ref_pos,
                self.coordinate)
            StatImageChannel.fetch_epigenetic_channel(
                self.normal_processed_stat_img[:, :, 8], fold, start_ref_pos,
                self.coordinate)


    def generate_info(self):
        """Generate image information
        """

        self.gen_unnorm_hash()

        self.generate_reads()
        self.generate_imgs()
        self.generate_lodt()
        self.generate_lodn()
        self.generate_unnorm_lodt()
        self.generate_unnorm_lodn()
        self.fetch_epigenetic_channel()
        self.fetch_stat_epigenetic_channel()


    def write_info(self):
        """Write image information
        """

        os.makedirs(self.output_path, exist_ok=True)
        with open(self.output_file, "wb") as f:
            pickle.dump(self.img_dict, f, pickle.HIGHEST_PROTOCOL)


    def run(self):
        """Run image generator to generate and write image information
        """

        self.generate_info()
        self.write_info()

        del self.insertion_size
        del self.coordinate

        del self.tumor_raw_reads_coordinate
        del self.normal_raw_reads_coordinate
        del self.tumor_processed_reads_coordinate
        del self.normal_processed_reads_coordinate

        del self.img_dict

        del self.tumor_raw_img
        del self.tumor_processed_img
        del self.normal_raw_img
        del self.normal_processed_img

        del self.tumor_raw_reads
        del self.tumor_processed_reads
        del self.normal_raw_reads
        del self.normal_processed_reads

        del self.tumor_unnorm_reads
        del self.tumor_processed_unnorm_reads
        del self.normal_unnorm_reads
        del self.normal_processed_unnorm_reads

        del self.unnorm_hash
        del self.unnorm_bam_files
