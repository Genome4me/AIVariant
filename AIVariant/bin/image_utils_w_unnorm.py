###
# Written by Hyeonseong Jeon
# email: jun0605@snu.ac.kr
###

import pysam
import numpy as np

class ReadFilter():
    """Collection of read filters that can be applied to get processed reads
    """

    base_quality_threshold = 5
    mapping_quality_threshold = 0
    softclip_threshold = 0.3
    mismatch_threshold = 100

    @staticmethod
    def delete_filter(read):
        """Check query position on the read is deleted or not
        """

        condition = read.query_position is not None
        return condition

    @staticmethod
    def base_quality_filter(read):
        """Check base quality at the query position is greater than or equal to base quality threshold
        """

        query_pos = read.query_position
        condition = (query_pos is not None) and (
                read.alignment.query_qualities[
                    query_pos] >= ReadFilter.base_quality_threshold)
        return condition

    @staticmethod
    def mapping_quality_filter(read):
        """Check mapping qulaity is greater than mapping qualtiy threshold
        """

        condition = read.alignment.mapping_quality > ReadFilter.mapping_quality_threshold
        return condition

    @staticmethod
    def softclip_filter(read):
        """Check softclip ration of read is less than softclip threshold
        """

        confirm_cigar = read.alignment.infer_read_length()
        if confirm_cigar is None: return False

        num_softclips = sum(
            [k[1] for k in read.alignment.cigartuples if k[0] == 4])

        read_sequence = read.alignment.query_sequence
        if read_sequence is None: return False
        read_length = len(read_sequence)

        condition = num_softclips / read_length < ReadFilter.softclip_threshold
        return condition

    @staticmethod
    def mismatch_filter(ref_file, chromosome, read):
        """Check mismatches of the read is not greater than mismatch threshold
        """

        read_sequence = read.alignment.query_sequence
        read_qualities = read.alignment.query_qualities
        if read_sequence is None or read_qualities is None: return False

        ref_data = pysam.FastaFile(ref_file)
        ref_pos = read.alignment.get_reference_positions(full_length=True)
        ref_seq = ""
        for pos in ref_pos:
            if pos is None:
                ref_seq += "x"
            else:
                ref_seq += ref_data.fetch(chromosome[3:], pos, pos + 1).upper()
        ref_data.close()

        num_mismatches = 0
        for idx, ref in enumerate(ref_seq):
            if ref != "x" and ref != (
                    read_sequence[idx].upper()): num_mismatches += read_qualities[idx]

        condition = num_mismatches <= ReadFilter.mismatch_threshold
        return condition


class ImageChannel():
    """Collection of image channel fetchers
    """

    @staticmethod
    def fetch_reference_channel(rows, img, ref_seq, start_ref_pos, coordinate):
        one_hot_idx = {"A": 0, "C": 1, "G": 2, "T": 3}
        for idx, ref in enumerate(ref_seq):
            ref_pos = start_ref_pos + idx
            if ref in one_hot_idx and (ref_pos, 0) in coordinate:
                col_idx = coordinate[(ref_pos, 0)]
                img[:rows, col_idx, one_hot_idx[ref]] = 1

    @staticmethod
    def fetch_epigenetic_channel(rows, img, pvals, folds, start_ref_pos,
                                 coordinate):
        for pval in pvals:
            for idx, pv in enumerate(pval):
                ref_pos = start_ref_pos + idx
                if (ref_pos, 0) in coordinate:
                    col_idx = coordinate[(ref_pos, 0)]
                    img[:rows, col_idx, 14] = pv
        for fold in folds:
            for idx, fl in enumerate(fold):
                ref_pos = start_ref_pos + idx
                if (ref_pos, 0) in coordinate:
                    col_idx = coordinate[(ref_pos, 0)]
                    img[:rows, col_idx, 15] = fl

    @staticmethod
    def fetch_nucleotide_channel(img, reads, read_coord):
        one_hot_idx = {"A": 0, "C": 1, "G": 2, "T": 3}

        for idx, read in enumerate(reads):
            query_seq = read.alignment.query_sequence
            for col_idx, read_pos, cigar_char in read_coord[idx]:
                if cigar_char in "IM" and query_seq[read_pos] in one_hot_idx:
                    img[idx, col_idx, one_hot_idx[query_seq[read_pos]]] = 1

    @staticmethod
    def fetch_indel_channel(indel_channel, arr_idx, cigarstring):
        window_size = indel_channel.shape[0]
        del_indexes = []
        ins_indexes = []
        pos_in_read = 0
        for op in cigarstring:
            if op == "D":
                del_indexes.append(pos_in_read)
            elif op == "I":
                ins_indexes.append((pos_in_read, pos_in_read + 1))
            if op == "D" or op == "M":
                pos_in_read += 1

        for del_index in del_indexes:
            del_index_in_arr = arr_idx + del_index
            if del_index_in_arr >= 0 and del_index_in_arr < window_size:
                indel_channel[del_index_in_arr] -= 1

        for ins_index in ins_indexes:
            ins_index_in_arr = (arr_idx + ins_index[0], arr_idx + ins_index[1])
            if ins_index_in_arr[0] >= 0 and ins_index_in_arr[0] < window_size:
                indel_channel[ins_index_in_arr[0]] += 0.1
            if ins_index_in_arr[1] >= 0 and ins_index_in_arr[1] < window_size:
                indel_channel[ins_index_in_arr[1]] += 0.1

    @staticmethod
    def fetch_mq_channel(img, reads, read_coord):

        for idx, read in enumerate(reads):
            mq = read.alignment.mapping_quality
            for col_idx, _, cigar_char in read_coord[idx]:
                if cigar_char in "IDM":
                    img[idx, col_idx] = mq

    @staticmethod
    def fetch_strand_channel(img, reads, read_coord):

        for idx, read in enumerate(reads):
            strand = -1 if read.alignment.is_reverse else 1
            for col_idx, _, cigar_char in read_coord[idx]:
                if cigar_char in "IDM":
                    img[idx, col_idx] = strand

    @staticmethod
    def fetch_bq_channel(img, reads, read_coord):

        for idx, read in enumerate(reads):
            bqs = read.alignment.query_qualities
            for col_idx, read_pos, cigar_char in read_coord[idx]:
                if cigar_char in "IM":
                    img[idx, col_idx] = bqs[read_pos]

    @staticmethod
    def fetch_softclip_channel(img, reads, read_coord):

        for idx, read in enumerate(reads):
            cigartuples = read.alignment.cigartuples

            left_softclip_len = [0, cigartuples[0][1]][cigartuples[0][0] == 4]
            right_softclip_len = [0, cigartuples[-1][1]][
                cigartuples[-1][0] == 4]

            start_col_idx = read_coord[idx][0][0]
            end_col_idx = read_coord[idx][-1][0]

            for offset in range(1, left_softclip_len + 1):
                softclip_idx = start_col_idx - offset
                if softclip_idx >= 0:
                    img[idx, softclip_idx] = 1

            for offset in range(1, right_softclip_len + 1):
                softclip_idx = end_col_idx + offset
                if softclip_idx < len(img[0]):
                    img[idx, softclip_idx] = 1

    @staticmethod
    def fetch_distance_channel(img, reads, read_coord):

        for idx, read in enumerate(reads):
            read_len = read.alignment.query_length
            for col_idx, read_pos, cigar_char in read_coord[idx]:
                if cigar_char in "IM":
                    distance = min(read_pos + 1, read_len - read_pos) / float(
                        read_len // 2)
                    img[idx, col_idx] = distance


class ImageNormalizer():
    """Collection of image channel normalizers
    """

    mq_channel_idx = 8
    bq_channel_idx = 10
    # indel_channel_idx = 8
    unnorm_bq_channel_idx = 13

    max_mq = 60
    max_bq = 41

    @staticmethod
    def normalize_img(img):
        ImageNormalizer.normalize_mq(img)
        ImageNormalizer.normalize_bq(img)
        # ImageNormalizer.normalize_indel(img)
        ImageNormalizer.normalize_unnorm_bq(img)

    @classmethod
    def normalize_mq(cls, img):
        mq_img = img[:, :, cls.mq_channel_idx]
        mq_img = np.round(2. * mq_img / cls.max_mq - 1., 2)
        mq_img[np.where(mq_img > 1.0)] = 1.0
        img[:, :, cls.mq_channel_idx] = mq_img

    @classmethod
    def normalize_bq(cls, img):
        bq_img = img[:, :, cls.bq_channel_idx]
        bq_img = np.round(2. * bq_img / cls.max_bq - 1., 2)
        bq_img[np.where(bq_img > 1.0)] = 1.0
        img[:, :, cls.bq_channel_idx] = bq_img

    @classmethod
    def normalize_indel(cls, img):
        indel_img = img[:, :, cls.indel_channel_idx]
        indel_img[np.where(indel_img > 1.0)] = 1.0
        img[:, :, cls.indel_channel_idx] = indel_img

    @classmethod
    def normalize_unnorm_bq(cls, img):
        unnorm_bq_img = img[:, :, cls.unnorm_bq_channel_idx]
        unnorm_bq_img = np.round(2. * unnorm_bq_img / cls.max_bq - 1., 2)
        unnorm_bq_img[np.where(unnorm_bq_img > 1.0)] = 1.0
        img[:, :, cls.unnorm_bq_channel_idx] = unnorm_bq_img


class StatImageChannel():
    """Collection of stat image channel fetchers
    """

    max_mq = 60
    max_bq = 41

    def aggregate(img, reads, read_coord, read_parser, fill_gap=False):
        cnt = np.zeros(img.shape)
        row_idx = {"A": 1, "C": 2, "G": 3, "T": 4}
        for idx, read in enumerate(reads):
            query_seq = read.alignment.query_sequence
            val = read_parser(read)
            marks = [0 for _ in range(len(img[0]))]
            mark_start = None
            mark_end = None
            for col_idx, read_pos, cigar_char in read_coord[idx]:
                if cigar_char in "IM" and query_seq[read_pos] in row_idx:
                    if not hasattr(val, "__iter__"):
                        img[row_idx[query_seq[read_pos]], col_idx] += val
                        cnt[row_idx[query_seq[read_pos]], col_idx] += 1
                        marks[col_idx] = 1
                        if mark_start is None:
                            mark_start = col_idx
                        mark_end = col_idx
                    else:
                        img[row_idx[query_seq[read_pos]], col_idx] += val[
                            read_pos]
                        cnt[row_idx[query_seq[read_pos]], col_idx] += 1
                        marks[col_idx] = 1
                        if mark_start is None:
                            mark_start = col_idx
                        mark_end = col_idx
                if cigar_char == 'D':
                    if not hasattr(val, "__iter__"):
                        img[0, col_idx] += val
                        cnt[0, col_idx] += 1
                        marks[col_idx] = 1
                        if mark_start is None:
                            mark_start = col_idx
                        mark_end = col_idx

            if fill_gap:
                for i in range(mark_start, mark_end):
                    if marks[i] == 0:
                        img[0, i] += val
                        cnt[0, i] += 1

        return cnt

    def normalize_col(img, cnt):
        for j in range(len(cnt[0])):
            divisor = float(sum(cnt[:, j]))
            if divisor != 0:
                img[:, j] /= divisor

    def normalize_elem(img, cnt):
        for i in range(len(cnt)):
            for j in range(len(cnt[0])):
                divisor = float(cnt[i, j])
                if divisor != 0:
                    img[i, j] /= divisor

    @staticmethod
    def fetch_reference_channel(stat_img_channel, ref_seq, start_ref_pos,
                                coordinate):
        row_idx = {"A": 1, "C": 2, "G": 3, "T": 4}
        for idx, ref in enumerate(ref_seq):
            ref_pos = start_ref_pos + idx
            if ref in row_idx and (ref_pos, 0) in coordinate:
                col_idx = coordinate[(ref_pos, 0)]
                stat_img_channel[row_idx[ref], col_idx] = 1

            insert_idx = 1
            while (ref_pos, insert_idx) in coordinate:
                col_idx = coordinate[(ref_pos, insert_idx)]
                stat_img_channel[0, col_idx] += 1
                insert_idx += 1

    @staticmethod
    def fetch_epigenetic_channel(stat_img_channel, epi_seq, start_ref_pos,
                                 coordinate):
        for idx, epi in enumerate(epi_seq):
            ref_pos = start_ref_pos + idx
            if (ref_pos, 0) in coordinate:
                col_idx = coordinate[(ref_pos, 0)]
                stat_img_channel[:, col_idx] = epi

    @classmethod
    def fetch_frequency_channel(cls, img, reads, read_coord):
        cnt = cls.aggregate(img, reads, read_coord, lambda read: 1,
                            fill_gap=True)
        cls.normalize_col(img, cnt)

    @classmethod
    def fetch_mq_channel(cls, img, reads, read_coord):
        def read_parser(read):
            return np.round(
                2. * read.alignment.mapping_quality / cls.max_mq - 1., 2)

        cnt = cls.aggregate(img, reads, read_coord, read_parser, fill_gap=True)
        cls.normalize_elem(img, cnt)

    @classmethod
    def fetch_strand_channel(cls, img, reads, read_coord):
        def read_parser(read):
            return -1 if read.alignment.is_reverse else 1

        cnt = cls.aggregate(img, reads, read_coord, read_parser, fill_gap=True)
        cls.normalize_col(img, cnt)

    @classmethod
    def fetch_bq_channel(cls, img, reads, read_coord):
        def read_parser(read):
            return np.round(2. * np.array(
                read.alignment.query_qualities) / cls.max_bq - 1., 2)

        cnt = cls.aggregate(img, reads, read_coord, read_parser)
        cls.normalize_elem(img, cnt)

    @classmethod
    def fetch_distance_channel(cls, img, reads, read_coord):
        def read_parser(read):
            read_len = read.alignment.query_length
            return [
                min(read_pos + 1, read_len - read_pos) / float(read_len // 2)
                for read_pos in range(read_len)]

        cnt = cls.aggregate(img, reads, read_coord, read_parser)
        cls.normalize_elem(img, cnt)


class ReadsStat():
    """Collection of calculators for statistics of reads
    """

    max_lodt = 100
    max_lodn = 50

    @staticmethod
    def count_reads(reads, ref, alt):
        bases = []
        quals = []
        ref_count = 0
        alt_count = 0
        for read in reads:
            query_pos = read.query_position
            if query_pos is None: continue
            base = read.alignment.query_sequence[query_pos].upper()
            qual = read.alignment.query_qualities[query_pos]
            if base == ref:
                ref_count += 1
            elif base == alt:
                alt_count += 1
            bases.append(base)
            quals.append(qual)
        return bases, quals, ref_count, alt_count

    @staticmethod
    def calc_likelihood(base, ref, alt, frac, err):
        if base == ref:
            return frac * (err / 3) + (1 - frac) * (1 - err)
        elif base == alt:
            return frac * (1 - err) + (1 - frac) * (err / 3)
        else:
            return err / 3

    @staticmethod
    def calc_lodt(reads, ref, alt):
        bases, quals, ref_count, alt_count = ReadsStat.count_reads(reads, ref,
                                                                   alt)

        if alt_count == 0:
            lodt = 0.0
        else:
            likelihood_var = 0.0;
            likelihood_no_var = 0.0;
            frac = alt_count / (ref_count + alt_count)
            for idx, base in enumerate(bases):
                err = 10 ** (-(quals[idx]) / 10)
                likelihood_var += np.log10(
                    ReadsStat.calc_likelihood(base, ref, alt, frac, err))
                likelihood_no_var += np.log10(
                    ReadsStat.calc_likelihood(base, ref, alt, 0.0, err))
            lodt = likelihood_var - likelihood_no_var

        return lodt

    @staticmethod
    def calc_lodn(reads, ref, alt):
        bases, quals, _, _ = ReadsStat.count_reads(reads, ref, alt)

        likelihood_var = 0.0;
        likelihood_no_var = 0.0;
        for idx, base in enumerate(bases):
            err = 10 ** (-(quals[idx]) / 10)
            likelihood_var += np.log10(
                ReadsStat.calc_likelihood(base, ref, alt, 0.5, err))
            likelihood_no_var += np.log10(
                ReadsStat.calc_likelihood(base, ref, alt, 0.0, err))

        lodn = likelihood_no_var - likelihood_var

        return lodn

    @staticmethod
    def normalize_lodt(lodt):
        if lodt > ReadsStat.max_lodt:
            return 1.0
        else:
            return lodt / ReadsStat.max_lodt

    @staticmethod
    def normalize_lodn(lodn):
        if lodn > ReadsStat.max_lodn:
            return 1.0
        else:
            return lodn / ReadsStat.max_lodn
