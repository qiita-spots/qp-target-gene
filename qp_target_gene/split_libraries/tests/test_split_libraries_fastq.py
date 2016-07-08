# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from unittest import TestCase, main
from os.path import isdir, exists, join
from os import remove, close, environ
from shutil import rmtree
from tempfile import mkstemp, mkdtemp
from json import dumps
from gzip import GzipFile
from functools import partial

from qiita_client import QiitaClient, ArtifactInfo

from qp_target_gene.split_libraries.split_libraries_fastq import (
    generate_parameters_string, get_sample_names_by_run_prefix,
    generate_per_sample_fastq_command, generate_split_libraries_fastq_cmd,
    split_libraries_fastq)


CLIENT_ID = '19ndkO3oMKsoChjVVWluF7QkxHRfYhTKSFbAVt8IhK7gZgDaO4'
CLIENT_SECRET = ('J7FfQ7CQdOxuKhQAf1eoGgBAE81Ns8Gu3EKaWFm3IO2JKh'
                 'AmmCWZuabe0O5Mp28s1')


class SplitLibrariesFastqTests(TestCase):
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
            "max_bad_run_length": 3, "min_per_read_length_fraction": 0.75,
            "sequence_max_n": 0, "rev_comp_barcode": False,
            "rev_comp_mapping_barcodes": True, "rev_comp": False,
            "phred_quality_threshold": 3, "barcode_type": "golay_12",
            "max_barcode_errors": 1.5, "input_data": 1, "phred_offset": ""}

        obs = generate_parameters_string(parameters)
        exp = ("--max_bad_run_length 3 --min_per_read_length_fraction 0.75 "
               "--sequence_max_n 0 --phred_quality_threshold 3 "
               "--barcode_type golay_12 --max_barcode_errors 1.5 "
               "--rev_comp_mapping_barcodes")
        self.assertEqual(obs, exp)

    def test_get_sample_names_by_run_prefix(self):
        fd, fp = mkstemp()
        close(fd)
        with open(fp, 'w') as f:
            f.write(MAPPING_FILE)
        self._clean_up_files.append(fp)

        obs = get_sample_names_by_run_prefix(fp)
        exp = {'s3': 'SKB7.640196', 's2': 'SKD8.640184', 's1': 'SKB8.640193'}
        self.assertEqual(obs, exp)

    def test_get_sample_names_by_run_prefix_error(self):
        fd, fp = mkstemp()
        close(fd)
        with open(fp, 'w') as f:
            f.write(MAPPING_FILE_2)
        self._clean_up_files.append(fp)

        with self.assertRaises(ValueError):
            get_sample_names_by_run_prefix(fp)

    def test_generate_per_sample_fastq_command(self):
        fd, fp = mkstemp()
        close(fd)
        with open(fp, 'w') as f:
            f.write(MAPPING_FILE)
        self._clean_up_files.append(fp)
        forward_seqs = ["s1.fastq.gz", "s2.fastq.gz", "s3.fastq.gz"]
        reverse_seqs = ["s1_rev.fastq.gz", "s2_rev.fastq.gz",
                        "s3_rev.fastq.gz"]
        barcode_fps = []
        mapping_file = fp
        output_dir = "/output/dir"
        params_str = (
            "--max_bad_run_length 3 --min_per_read_length_fraction 0.75 "
            "--sequence_max_n 0 --phred_quality_threshold 3 "
            "--barcode_type golay_12 --max_barcode_errors 1.5 "
            "--rev_comp_mapping_barcodes")
        obs = generate_per_sample_fastq_command(
            forward_seqs, reverse_seqs, barcode_fps,
            mapping_file, output_dir, params_str)
        exp = ("split_libraries_fastq.py --store_demultiplexed_fastq -i "
               "s1.fastq.gz,s2.fastq.gz,s3.fastq.gz --sample_ids "
               "SKB8.640193,SKD8.640184,SKB7.640196 -o /output/dir "
               "--max_bad_run_length 3 --min_per_read_length_fraction 0.75 "
               "--sequence_max_n 0 --phred_quality_threshold 3 "
               "--barcode_type golay_12 --max_barcode_errors 1.5 "
               "--rev_comp_mapping_barcodes")
        self.assertEqual(obs, exp)

    def test_generate_per_sample_fastq_command_regex(self):
        fd, fp = mkstemp()
        close(fd)
        with open(fp, 'w') as f:
            f.write(MAPPING_FILE)
        self._clean_up_files.append(fp)
        forward_seqs = ["1_s1.fastq.gz", "1_s2.fastq.gz", "1_s3.fastq.gz"]
        reverse_seqs = ["1_s1_rev.fastq.gz", "1_s2_rev.fastq.gz",
                        "1_s3_rev.fastq.gz"]
        barcode_fps = []
        mapping_file = fp
        output_dir = "/output/dir"
        params_str = (
            "--max_bad_run_length 3 --min_per_read_length_fraction 0.75 "
            "--sequence_max_n 0 --phred_quality_threshold 3 "
            "--barcode_type golay_12 --max_barcode_errors 1.5 "
            "--rev_comp_mapping_barcodes")
        obs = generate_per_sample_fastq_command(
            forward_seqs, reverse_seqs, barcode_fps,
            mapping_file, output_dir, params_str)
        exp = ("split_libraries_fastq.py --store_demultiplexed_fastq -i "
               "1_s1.fastq.gz,1_s2.fastq.gz,1_s3.fastq.gz --sample_ids "
               "SKB8.640193,SKD8.640184,SKB7.640196 -o /output/dir "
               "--max_bad_run_length 3 --min_per_read_length_fraction 0.75 "
               "--sequence_max_n 0 --phred_quality_threshold 3 "
               "--barcode_type golay_12 --max_barcode_errors 1.5 "
               "--rev_comp_mapping_barcodes")
        self.assertEqual(obs, exp)

    def test_generate_per_sample_fastq_command_error_barcodes(self):
        fd, fp = mkstemp()
        close(fd)
        with open(fp, 'w') as f:
            f.write(MAPPING_FILE)
        self._clean_up_files.append(fp)
        forward_seqs = ["s1.fastq.gz", "s2.fastq.gz", "s3.fastq.gz"]
        reverse_seqs = ["s1_rev.fastq.gz", "s2_rev.fastq.gz",
                        "s3_rev.fastq.gz"]
        barcode_fps = ["s1_barcodes.fastq.gz", "s2_barcodes.fastq.gz",
                       "s3_barcodes.fastq.gz"]
        mapping_file = fp
        output_dir = "/output/dir"
        params_str = (
            "--max_bad_run_length 3 --min_per_read_length_fraction 0.75 "
            "--sequence_max_n 0 --phred_quality_threshold 3 "
            "--barcode_type golay_12 --max_barcode_errors 1.5 "
            "--rev_comp_mapping_barcodes")
        with self.assertRaises(ValueError):
            generate_per_sample_fastq_command(
                forward_seqs, reverse_seqs, barcode_fps,
                mapping_file, output_dir, params_str)

    def test_generate_per_sample_fastq_command_error_prefixes(self):
        fd, fp = mkstemp()
        close(fd)
        with open(fp, 'w') as f:
            f.write(MAPPING_FILE)
        self._clean_up_files.append(fp)
        forward_seqs = ["s1.fastq.gz", "s2.fastq.gz", "sX.fastq.gz"]
        reverse_seqs = ["s1_rev.fastq.gz", "s2_rev.fastq.gz",
                        "sX_rev.fastq.gz"]
        barcode_fps = []
        mapping_file = fp
        output_dir = "/output/dir"
        params_str = (
            "--max_bad_run_length 3 --min_per_read_length_fraction 0.75 "
            "--sequence_max_n 0 --phred_quality_threshold 3 "
            "--barcode_type golay_12 --max_barcode_errors 1.5 "
            "--rev_comp_mapping_barcodes")
        with self.assertRaises(ValueError):
            generate_per_sample_fastq_command(
                forward_seqs, reverse_seqs, barcode_fps,
                mapping_file, output_dir, params_str)

    def test_generate_split_libraries_fastq_cmd_per_sample_FASTQ(self):
        fps = {
            "raw_forward_seqs": ["s1.fastq.gz", "s2.fastq.gz", "s3.fastq.gz"],
            "raw_reverse_seqs": ["s1_rev.fastq.gz", "s2_rev.fastq.gz",
                                 "s3_rev.fastq.gz"]}
        fd, fp = mkstemp()
        with open(fp, 'w') as f:
            f.write(MAPPING_FILE)
        self._clean_up_files.append(fp)
        mapping_file = fp
        atype = "per_sample_FASTQ"
        out_dir = "/output/dir"
        parameters = {
            "max_bad_run_length": 3, "min_per_read_length_fraction": 0.75,
            "sequence_max_n": 0, "rev_comp_barcode": False,
            "rev_comp_mapping_barcodes": True, "rev_comp": False,
            "phred_quality_threshold": 3, "barcode_type": "golay_12",
            "max_barcode_errors": 1.5, "input_data": 1, "phred_offset": ""}
        obs_cmd, obs_outdir = generate_split_libraries_fastq_cmd(
            fps, mapping_file, atype, out_dir, parameters)

        exp_cmd = (
            "split_libraries_fastq.py --store_demultiplexed_fastq -i "
            "s1.fastq.gz,s2.fastq.gz,s3.fastq.gz --sample_ids "
            "SKB8.640193,SKD8.640184,SKB7.640196 -o /output/dir/sl_out "
            "--max_bad_run_length 3 --min_per_read_length_fraction 0.75 "
            "--sequence_max_n 0 --phred_quality_threshold 3 "
            "--barcode_type golay_12 --max_barcode_errors 1.5 "
            "--rev_comp_mapping_barcodes")
        self.assertEqual(obs_cmd, exp_cmd)
        self.assertEqual(obs_outdir, "/output/dir/sl_out")

    def test_generate_split_libraries_fastq_cmd(self):
        out_dir = mkdtemp()
        fps = {
            "raw_forward_seqs": ["s1.fastq.gz", "s2.fastq.gz", "s3.fastq.gz"],
            "raw_reverse_seqs": ["s1_rev.fastq.gz", "s2_rev.fastq.gz",
                                 "s3_rev.fastq.gz"],
            "raw_barcodes": ["s1_barcodes.fastq.gz", "s2_barcodes.fastq.gz",
                             "s3_barcodes.fastq.gz"],
            "html_summary": ["artifact_summary.html"]}
        self._clean_up_files.append(out_dir)
        fd, fp = mkstemp()
        with open(fp, 'w') as f:
            f.write(MAPPING_FILE)
        self._clean_up_files.append(fp)
        mapping_file = fp
        atype = "FASTQ"
        parameters = {
            "max_bad_run_length": 3, "min_per_read_length_fraction": 0.75,
            "sequence_max_n": 0, "rev_comp_barcode": False,
            "rev_comp_mapping_barcodes": True, "rev_comp": False,
            "phred_quality_threshold": 3, "barcode_type": "golay_12",
            "max_barcode_errors": 1.5, "input_data": 1, "phred_offset": ""}
        obs_cmd, obs_outdir = generate_split_libraries_fastq_cmd(
            fps, mapping_file, atype, out_dir, parameters)
        exp_cmd = (
            "split_libraries_fastq.py --store_demultiplexed_fastq -i "
            "s1.fastq.gz,s2.fastq.gz,s3.fastq.gz -b "
            "s1_barcodes.fastq.gz,s2_barcodes.fastq.gz,s3_barcodes.fastq.gz "
            "-m {0}/mappings/s1_mapping_file.txt,"
            "{0}/mappings/s2_mapping_file.txt,"
            "{0}/mappings/s3_mapping_file.txt "
            "-o {0}/sl_out --max_bad_run_length 3 "
            "--min_per_read_length_fraction 0.75 --sequence_max_n 0 "
            "--phred_quality_threshold 3 --barcode_type golay_12 "
            "--max_barcode_errors 1.5 "
            "--rev_comp_mapping_barcodes".format(out_dir))
        self.assertEqual(obs_cmd, exp_cmd)
        self.assertEqual(obs_outdir, join(out_dir, "sl_out"))

    def test_generate_split_libraries_fastq_cmd_valueerror(self):
        out_dir = "/output/dir"
        fps = {
            "raw_forward_seqs": ["s1.fastq.gz", "s2.fastq.gz", "s3.fastq.gz"],
            "raw_reverse_seqs": ["s1_rev.fastq.gz", "s2_rev.fastq.gz",
                                 "s3_rev.fastq.gz"],
            "raw_barcodes": ["s1_barcodes.fastq.gz", "s2_barcodes.fastq.gz"]}
        mapping_file = "mapping_file.txt"
        atype = "FASTQ"
        parameters = {
            "max_bad_run_length": 3, "min_per_read_length_fraction": 0.75,
            "sequence_max_n": 0, "rev_comp_barcode": False,
            "rev_comp_mapping_barcodes": True, "rev_comp": False,
            "phred_quality_threshold": 3, "barcode_type": "golay_12",
            "max_barcode_errors": 1.5, "input_data": 1, "phred_offset": ""}
        with self.assertRaises(ValueError):
            generate_split_libraries_fastq_cmd(
                fps, mapping_file, atype, out_dir, parameters)

    def test_split_libraries_fastq(self):
        # Create a new job
        parameters = {"max_bad_run_length": 3,
                      "min_per_read_length_fraction": 0.75,
                      "sequence_max_n": 0,
                      "rev_comp_barcode": False,
                      "rev_comp_mapping_barcodes": True,
                      "rev_comp": False,
                      "phred_quality_threshold": 3,
                      "barcode_type": "golay_12",
                      "max_barcode_errors": 1.5,
                      "phred_offset": "",
                      "input_data": 1}
        data = {'user': 'demo@microbio.me',
                'command': 1,
                'status': 'running',
                'parameters': dumps(parameters)}
        job_id = self.qclient.post(
            '/apitest/processing_job/', data=data)['job']

        # Create the files that Qiita returns (they don't exist in the test)
        files = self.qclient.get('/qiita_db/artifacts/1/')['files']
        bcds_fp = files['raw_barcodes'][0]
        self._clean_up_files.append(bcds_fp)
        with GzipFile(bcds_fp, mode='w') as fh:
            fh.write(BARCODES)
        fwd_fp = files['raw_forward_seqs'][0]
        self._clean_up_files.append(fwd_fp)
        with GzipFile(fwd_fp, mode='w') as fh:
            fh.write(READS)

        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        obs_success, obs_ainfo, obs_msg = split_libraries_fastq(
            self.qclient, job_id, parameters, out_dir)
        self.assertTrue(obs_success)
        path_builder = partial(join, out_dir, 'sl_out')
        filepaths = [
            (path_builder('seqs.fna'), 'preprocessed_fasta'),
            (path_builder('seqs.fastq'), 'preprocessed_fastq'),
            (path_builder('seqs.demux'), 'preprocessed_demux'),
            (path_builder('split_library_log.txt'), 'log')]
        exp_ainfo = [ArtifactInfo('demultiplexed', 'Demultiplexed', filepaths)]
        self.assertEqual(obs_ainfo, exp_ainfo)
        self.assertEqual(obs_msg, "")


MAPPING_FILE = (
    "#SampleID\tplatform\tbarcode\texperiment_design_description\t"
    "library_construction_protocol\tcenter_name\tprimer\trun_prefix\t"
    "instrument_model\tDescription\n"
    "SKB7.640196\tILLUMINA\tA\tA\tA\tANL\tA\ts3\tIllumina MiSeq\tdesc1\n"
    "SKB8.640193\tILLUMINA\tA\tA\tA\tANL\tA\ts1\tIllumina MiSeq\tdesc2\n"
    "SKD8.640184\tILLUMINA\tA\tA\tA\tANL\tA\ts2\tIllumina MiSeq\tdesc3\n"
)

MAPPING_FILE_2 = (
    "#SampleID\tplatform\tbarcode\texperiment_design_description\t"
    "library_construction_protocol\tcenter_name\tprimer\t"
    "run_prefix\tinstrument_model\tDescription\n"
    "SKB7.640196\tILLUMINA\tA\tA\tA\tANL\tA\ts3\tIllumina MiSeq\tdesc1\n"
    "SKB8.640193\tILLUMINA\tA\tA\tA\tANL\tA\ts1\tIllumina MiSeq\tdesc2\n"
    "SKD8.640184\tILLUMINA\tA\tA\tA\tANL\tA\ts1\tIllumina MiSeq\tdesc3\n"
)

READS = """@M00176:18:000000000-A0DK4:1:1:15579:1518 1:N:0:0
GACAGAGGGTGCAAACGTTGCTCGGAATCACTGGGCGTAAAGGGCGTGTAGGCGGGAAGGATAGTCAGATGTGAAATCCCTGGGCTCAACCCAGGAACTGCATTTGAAACTCCCTGTCTTGAGTGTCGGAGAGGGTAGCGGTATTCCTGGT
+
==?+ADDDD2ADDEAEEEEIAEFFIIIIEEIDEIIDI@DD>CCDIA-:=CDDDDA8&)88(+:>+:>A>AADAEEAAAAA:A?>>9<4>?<;(+8<88:>A3>AA>:>>AA:>8((+:>A>>A>:>(5939>??.+><())))+4(++:3(
@M00176:18:000000000-A0DK4:1:1:14647:1519 1:N:0:0
TACGTAGGGGCCGAGCGTTGTCCGGAATTATTGGGCGTAAAGAGCTCGTAGGCGGGTTGTCACGTCTGCGGGGAAATCCCGAGGCTCAACCGCGGATCTGCAGCTGTGCCGGTCTGACTGGCGTTCGGGAGTGTAGGAGTGAATTCCTAGT
+
=?@=ABD?+<<C)<CF1:?DGGFFF@@DD=4CF>=FF4AC=A3A@DB6,(51;;')507&>::52830>&)&)05?4+4><-9-5+(44?()&0&&)04::(((((((++)&)9&3(+:3((0&&()))5<((+((+(((+(4((+(((++
@M00176:18:000000000-A0DK4:1:1:17271:1520 1:N:0:0
TACAGAGGTCCCAAGCGTTGTTCGGATTTACTGGGCGTAAAGGGTGCGTAGGCGGTCGGTTAAGTCTGTTGTGAAATCTCCCGGCTCAACCGTGAAACTGCATAGGATACTATTTAGCTCGGGGACTGAAGTGGGGACTTGTATACTCGGT
+
??@=?DDDD?<C<AFEFII12AHFE7FGIED:FFDGAFFBBF>AFE1AEFF<)99?/'/'35::@BAA+(49?BBB>ABA(+09999@:9<9(+<<@(:(44@++(((+:@:>3>(:3++&&&)&&(+((4+((2&(&((+(+44:((+2<
@M00176:18:000000000-A0DK4:1:1:14752:1521 1:N:0:0
TACGGGAGTGGCAAGCGTTATCCGGAATTATTGGGCGTAAAGCGCGCGTAGGCGGTTTGTTAAGTTGGATGTGAAAGCCCTGGGCTCAACCTAGGAACTGCATCCGAAACTAACTCACTAGAGTACGATAGAGGGAGGTGGAATTCCTAGT
+
=?@D+=8@:AFFFHGFGFHIEHIGHC:DHG4DGHCHFDEEICC@8/9A'9?B(8805&++9ADCE3:A3:>C@>@C4>(23<??BC1:CC294>(+>8@C@+++(&2(05A((+((+3:(::>>(+8(+(53:?9595?(28:@CC(4(+:
@M00176:18:000000000-A0DK4:1:1:16117:1522 1:N:0:0
AACGCCCTCTTAAGGATATTCGCGATGCGTATAATTACCCCAAAAAGAAAGGTATTAAGGATGAGTGTTCAAGATTGCTGTAGGCCTCCACTATGAAATCGCGTAGAGGCTTTGCTATTCAGCGTTTGATGAATGCAATGCGACAGGCTCA
+
::?=A)=AFFADF4:::+3AADFD6?F10?60??BB?DEDFFFA@FFI()7=7;)7?7=;);???@@;>(6;5(;:>BBA(,,::>>?@B943>BB(+:4(+00<.<A18?@(44:@>@>+445820?(4>:(+:(+(++3&(0098@?B4
"""

BARCODES = """@M00176:18:000000000-A0DK4:1:1:15579:1518 1:N:0:0
TAGTCAGGCCAT
+
=:8BBD;DAH:<
@M00176:18:000000000-A0DK4:1:1:14647:1519 1:N:0:0
TAGTCAGGCCAT
+
=:8BBD;DAH:<
@M00176:18:000000000-A0DK4:1:1:17271:1520 1:N:0:0
CGTAGAGCTCTC
+
=:8BBD;DAH:<
@M00176:18:000000000-A0DK4:1:1:14752:1521 1:N:0:0
CGTAGAGCTCTC
+
=:8BBD;DAH:<
@M00176:18:000000000-A0DK4:1:1:16117:1522 1:N:0:0
CCTCTGAGAGCT
+
=:8BBD;DAH:<
"""

if __name__ == '__main__':
    main()
