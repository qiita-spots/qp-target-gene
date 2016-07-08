# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from unittest import TestCase, main
from os.path import isdir, exists, join, dirname
from os import remove, environ, close
from shutil import rmtree, copyfile
from tempfile import mkdtemp, mkstemp
from json import dumps
from functools import partial

from qiita_client import QiitaClient, ArtifactInfo

from qp_target_gene.split_libraries.split_libraries import (
    generate_parameters_string, generate_process_sff_commands,
    generate_split_libraries_cmd, split_libraries)

CLIENT_ID = '19ndkO3oMKsoChjVVWluF7QkxHRfYhTKSFbAVt8IhK7gZgDaO4'
CLIENT_SECRET = ('J7FfQ7CQdOxuKhQAf1eoGgBAE81Ns8Gu3EKaWFm3IO2JKh'
                 'AmmCWZuabe0O5Mp28s1')


class SplitLibrariesTests(TestCase):
    @classmethod
    def setUpClass(cls):
        server_cert = environ.get('QIITA_SERVER_CERT', None)
        cls.qclient = QiitaClient("https://localhost:21174", CLIENT_ID,
                                  CLIENT_SECRET, server_cert=server_cert)

    @classmethod
    def tearDownClass(cls):
        cls.qclient.post('/apitest/reset/')

    def setUp(self):
        self._clean_up_files = []

    def tearDown(self):
        for fp in self._clean_up_files:
            if exists(fp):
                if isdir(fp):
                    rmtree(fp)
                else:
                    remove(fp)

    def test_generate_parameters_string(self):
        parameters = {
            "min_seq_len": 200, "max_seq_len": 1000, "trim_seq_length": False,
            "min_qual_score": 25, "max_ambig": 6, "max_homopolymer": 6,
            "max_primer_mismatch": 0, "barcode_type": "golay_12",
            "max_barcode_errors": 1.5, "disable_bc_correction": False,
            "qual_score_window": 0, "disable_primers": False,
            "reverse_primers": "disable", "reverse_primer_mismatches": 0,
            "truncate_ambi_bases": False, "input_data": 1}
        obs = generate_parameters_string(parameters)
        exp = ("--min_seq_len 200 --max_seq_len 1000 --min_qual_score 25 "
               "--max_ambig 6 --max_homopolymer 6 --max_primer_mismatch 0 "
               "--barcode_type golay_12 --max_barcode_errors 1.5 "
               "--qual_score_window 0 --reverse_primer_mismatches 0 "
               "--reverse_primers disable")
        self.assertEqual(obs, exp)

    def test_generate_parameters_string_error(self):
        parameters = {
            "min_seq_len": 200, "max_seq_len": 1000, "trim_seq_length": False,
            "min_qual_score": 25, "max_ambig": 6, "max_homopolymer": 6,
            "max_primer_mismatch": 0, "barcode_type": "golay_12",
            "max_barcode_errors": 1.5, "disable_bc_correction": False,
            "qual_score_window": 0, "disable_primers": False,
            "reverse_primers": "whops", "reverse_primer_mismatches": 0,
            "truncate_ambi_bases": False, "input_data": 1}
        with self.assertRaises(ValueError):
            generate_parameters_string(parameters)

    def test_generate_process_sff_commands(self):
        out_dir = "/directory/output/"
        sff_fps = ["/directory/file1.sff", "/directory/file2.sff.gz"]
        obs_cmds, obs_seqs, obs_quals = generate_process_sff_commands(
            sff_fps, out_dir)

        exp_cmds = [
            "process_sff.py -i /directory/file1.sff -o /directory/output/",
            "process_sff.py -i /directory/file2.sff.gz -o /directory/output/"]
        exp_seqs = ["/directory/output/file1.fna",
                    "/directory/output/file2.fna"]
        exp_quals = ["/directory/output/file1.qual",
                     "/directory/output/file2.qual"]

        self.assertEqual(obs_cmds, exp_cmds)
        self.assertEqual(obs_seqs, exp_seqs)
        self.assertEqual(obs_quals, exp_quals)

    def test_generate_split_libraries_cmd_single(self):
        test_dir = mkdtemp()
        self._clean_up_files.append(test_dir)
        seqs = [join(test_dir, "seqs.fna")]
        quals = [join(test_dir, "seqs.qual")]
        mapping_file = join(test_dir, "mapping_file.txt")
        with open(mapping_file, 'w') as f:
            f.write(MAPPING_FILE_SINGLE)
        out_dir = join(test_dir, 'sl_out')
        params = {
            "min_seq_len": 200, "max_seq_len": 1000, "trim_seq_length": False,
            "min_qual_score": 25, "max_ambig": 6, "max_homopolymer": 6,
            "max_primer_mismatch": 0, "barcode_type": "golay_12",
            "max_barcode_errors": 1.5, "disable_bc_correction": False,
            "qual_score_window": 0, "disable_primers": False,
            "reverse_primers": "disable", "reverse_primer_mismatches": 0,
            "truncate_ambi_bases": False, "input_data": 1}
        obs_cmd, obs_outdir = generate_split_libraries_cmd(
            seqs, quals, mapping_file, out_dir, params)
        exp_cmd = [
            "split_libraries.py -f {0}/seqs.fna -m {0}/mapping_file.txt "
            "-q {0}/seqs.qual -d -o {0}/sl_out --min_seq_len 200 "
            "--max_seq_len 1000 --min_qual_score 25 --max_ambig 6 "
            "--max_homopolymer 6 --max_primer_mismatch 0 "
            "--barcode_type golay_12 --max_barcode_errors 1.5 "
            "--qual_score_window 0 --reverse_primer_mismatches 0 "
            "--reverse_primers disable".format(test_dir)]
        self.assertEqual(obs_cmd, exp_cmd)
        self.assertEqual(obs_outdir, [out_dir])

    def test_generate_split_libraries_cmd_mutliple(self):
        test_dir = mkdtemp()
        self._clean_up_files.append(test_dir)
        seqs = [join(test_dir, "prefix_1_seqs.fna"),
                join(test_dir, "prefix_2_seqs.fna")]
        quals = [join(test_dir, "prefix_1_seqs.qual"),
                 join(test_dir, "prefix_2_seqs.qual")]
        mapping_file = join(test_dir, "mapping_file.txt")
        with open(mapping_file, 'w') as f:
            f.write(MAPPING_FILE_MULT)
        out_dir = join(test_dir, 'sl_out')
        params = {
            "min_seq_len": 200, "max_seq_len": 1000, "trim_seq_length": False,
            "min_qual_score": 25, "max_ambig": 6, "max_homopolymer": 6,
            "max_primer_mismatch": 0, "barcode_type": "golay_12",
            "max_barcode_errors": 1.5, "disable_bc_correction": False,
            "qual_score_window": 0, "disable_primers": False,
            "reverse_primers": "disable", "reverse_primer_mismatches": 0,
            "truncate_ambi_bases": False, "input_data": 1}
        obs_cmd, obs_outdir = generate_split_libraries_cmd(
            seqs, quals, mapping_file, out_dir, params)
        exp_cmd = [
            "split_libraries.py -f {0}/prefix_1_seqs.fna -m "
            "{0}/sl_out/mappings/prefix_1_mapping_file.txt -q "
            "{0}/prefix_1_seqs.qual -d -o {0}/sl_out/prefix_1_mapping_file "
            "-n 1 --min_seq_len 200 --max_seq_len 1000 --min_qual_score 25 "
            "--max_ambig 6 --max_homopolymer 6 --max_primer_mismatch 0 "
            "--barcode_type golay_12 --max_barcode_errors 1.5 "
            "--qual_score_window 0 --reverse_primer_mismatches 0 "
            "--reverse_primers disable".format(test_dir),
            "split_libraries.py -f {0}/prefix_2_seqs.fna -m "
            "{0}/sl_out/mappings/prefix_2_mapping_file.txt -q "
            "{0}/prefix_2_seqs.qual -d -o {0}/sl_out/prefix_2_mapping_file "
            "-n 800000 --min_seq_len 200 --max_seq_len 1000 --min_qual_score "
            "25 --max_ambig 6 --max_homopolymer 6 --max_primer_mismatch 0 "
            "--barcode_type golay_12 --max_barcode_errors 1.5 "
            "--qual_score_window 0 --reverse_primer_mismatches 0 "
            "--reverse_primers disable".format(test_dir)]
        self.assertEqual(obs_cmd, exp_cmd)
        exp_outdir = [join(out_dir, 'prefix_1_mapping_file'),
                      join(out_dir, 'prefix_2_mapping_file')]
        self.assertEqual(obs_outdir, exp_outdir)

    def _create_qiita_bits(self, files):
        # Create a new prep template
        prep_info = {
            'SKB2.640194': {'barcode': 'AGCACGAGCCTA',
                            'primer': 'YATGCTGCCTCCCGTAGGAGT'},
            'SKM4.640180': {'barcode': 'AACTCGTCGATG',
                            'primer': 'YATGCTGCCTCCCGTAGGAGT'},
            'SKB3.640195': {'barcode': 'ACAGACCACTCA',
                            'primer': 'YATGCTGCCTCCCGTAGGAGT'},
            'SKB6.640176': {'barcode': 'ACCAGCGACTAG',
                            'primer': 'YATGCTGCCTCCCGTAGGAGT'},
            'SKD6.640190': {'barcode': 'AGCAGCACTTGT',
                            'primer': 'YATGCTGCCTCCCGTAGGAGT'},
            'SKM6.640187': {'barcode': 'AACTGTGCGTAC',
                            'primer': 'YATGCTGCCTCCCGTAGGAGT'},
            'SKD9.640182': {'barcode': 'ACAGAGTCGGCT',
                            'primer': 'YATGCTGCCTCCCGTAGGAGT'},
            'SKM8.640201': {'barcode': 'ACCGCAGAGTCA',
                            'primer': 'YATGCTGCCTCCCGTAGGAGT'},
            'SKM2.640199': {'barcode': 'ACGGTGAGTGTC',
                            'primer': 'YATGCTGCCTCCCGTAGGAGT'}}
        data = {'prep_info': dumps(prep_info),
                'study': 1,
                'data_type': '16S'}
        template = self.qclient.post(
            '/apitest/prep_template/', data=data)['prep']

        # Create a new artifact attached to the previous template
        data = {'filepaths': dumps(files),
                'type': 'FASTA',
                'name': 'Test artitact for split libraries',
                'prep': template}
        artifact = self.qclient.post(
            '/apitest/artifact/', data=data)['artifact']

        # Create a new job
        parameters = {"min_seq_len": 200,
                      "max_seq_len": 1000,
                      "trim_seq_length": False,
                      "min_qual_score": 25,
                      "max_ambig": 6,
                      "max_homopolymer": 6,
                      "max_primer_mismatch": 0,
                      "barcode_type": "golay_12",
                      "max_barcode_errors": 1.5,
                      "disable_bc_correction": False,
                      "qual_score_window": 0,
                      "disable_primers": False,
                      "reverse_primers": "disable",
                      "reverse_primer_mismatches": 0,
                      "truncate_ambi_bases": False,
                      "input_data": artifact}
        data = {'user': 'demo@microbio.me',
                'command': 2,
                'status': 'running',
                'parameters': dumps(parameters)}
        job_id = self.qclient.post(
            '/apitest/processing_job/', data=data)['job']

        return job_id, parameters

    def test_split_libraries(self):
        fd, reads_fp = mkstemp(suffix='.fna.gz')
        close(fd)
        self._clean_up_files.append(reads_fp)
        fd, qual_fp = mkstemp(suffix='.qual.gz')
        close(fd)
        self._clean_up_files.append(qual_fp)

        path_builder = partial(join, dirname(__file__), 'test_data')
        copyfile(path_builder('example_seqs.fna.gz'), reads_fp)
        copyfile(path_builder('example_quals.fna.gz'), qual_fp)
        files = [(reads_fp, 'raw_fasta'), (qual_fp, 'raw_qual')]
        job_id, parameters = self._create_qiita_bits(files)

        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        obs_success, obs_ainfo, obs_msg = split_libraries(
            self.qclient, job_id, parameters, out_dir)

        self.assertTrue(obs_success)
        path_builder = partial(join, out_dir, 'sl_out')
        fps = [(path_builder('seqs.fna'), 'preprocessed_fasta'),
               (path_builder('seqs.fastq'), 'preprocessed_fastq'),
               (path_builder('seqs.demux'), 'preprocessed_demux'),
               (path_builder('split_library_log.txt'), 'log')]
        exp_ainfo = [ArtifactInfo('demultiplexed', 'Demultiplexed', fps)]
        self.assertEqual(obs_ainfo, exp_ainfo)
        self.assertEqual(obs_msg, "")

    def test_split_libraries_sff(self):
        fd, sff_fp = mkstemp(suffix='.sff.gz')
        close(fd)
        self._clean_up_files.append(sff_fp)
        copyfile(join(dirname(__file__), 'test_data', 'example_sff.sff.gz'),
                 sff_fp)
        job_id, parameters = self._create_qiita_bits([(sff_fp, 'raw_sff')])

        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        obs_success, obs_ainfo, obs_msg = split_libraries(
            self.qclient, job_id, parameters, out_dir)

        self.assertTrue(obs_success)
        path_builder = partial(join, out_dir, 'sl_out')
        fps = [(path_builder('seqs.fna'), 'preprocessed_fasta'),
               (path_builder('seqs.fastq'), 'preprocessed_fastq'),
               (path_builder('seqs.demux'), 'preprocessed_demux'),
               (path_builder('split_library_log.txt'), 'log')]
        exp_ainfo = [ArtifactInfo('demultiplexed', 'Demultiplexed', fps)]
        self.assertEqual(obs_ainfo, exp_ainfo)
        self.assertEqual(obs_msg, "")


MAPPING_FILE_SINGLE = (
    "#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tDescription\n"
    "Sample1\tGTCCGCAAGTTA\tGTGCCAGCMGCCGCGGTAA\tTGP test\n"
    "Sample2\tCGTAGAGCTCTC\tGTGCCAGCMGCCGCGGTAA\tTGP test\n"
)


MAPPING_FILE_MULT = (
    "#SampleID\tBarcodeSequence\tLinkerPrimerSequence\trun_prefix\t"
    "Description\n"
    "Sample1\tGTCCGCAAGTTA\tGTGCCAGCMGCCGCGGTAA\tprefix_1\tTGP test\n"
    "Sample2\tCGTAGAGCTCTC\tGTGCCAGCMGCCGCGGTAA\tprefix_2\tTGP test\n"
)


if __name__ == '__main__':
    main()
