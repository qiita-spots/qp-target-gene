# -*- coding: utf-8 -*-

# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from unittest import TestCase, main
from os.path import isdir, exists, join, basename
from shutil import rmtree
from os import remove, close, environ
from tempfile import mkdtemp, mkstemp

from qiita_client import QiitaClient, ArtifactInfo

from qp_target_gene.split_libraries.util import (
    get_artifact_information, split_mapping_file, generate_demux_file,
    generate_artifact_info)


CLIENT_ID = '19ndkO3oMKsoChjVVWluF7QkxHRfYhTKSFbAVt8IhK7gZgDaO4'
CLIENT_SECRET = ('J7FfQ7CQdOxuKhQAf1eoGgBAE81Ns8Gu3EKaWFm3IO2JKh'
                 'AmmCWZuabe0O5Mp28s1')


class UtilTests(TestCase):
    @classmethod
    def setUpClass(cls):
        server_cert = environ.get('QIITA_SERVER_CERT', None)
        cls.qclient = QiitaClient("https://localhost:21174", CLIENT_ID,
                                  CLIENT_SECRET, server_cert=server_cert)

    @classmethod
    def tearDownClass(cls):
        cls.qclient.post("/apitest/reset/")

    def setUp(self):
        self._clean_up_files = []

    def tearDown(self):
        for fp in self._clean_up_files:
            if exists(fp):
                if isdir(fp):
                    rmtree(fp)
                else:
                    remove(fp)

    def test_get_artifact_information(self):
        obs_fps, obs_map, obs_at = get_artifact_information(self.qclient, 1)

        for k in obs_fps:
            obs_fps[k] = [basename(v) for v in obs_fps[k]]
        exp_fps = {
            'raw_barcodes': ['1_s_G1_L001_sequences_barcodes.fastq.gz'],
            'raw_forward_seqs': ['1_s_G1_L001_sequences.fastq.gz']}
        self.assertEqual(obs_fps, exp_fps)
        self.assertEqual(obs_at, "FASTQ")
        self.assertEqual(basename(obs_map),
                         '1_prep_1_qiime_19700101-000000.txt')

    def test_split_mapping_file_single(self):
        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        fd, fp = mkstemp(suffix='_map.txt')
        close(fd)
        self._clean_up_files.append(fp)

        with open(fp, 'w') as f:
            f.write(MAPPING_FILE_SINGLE)

        obs = split_mapping_file(fp, out_dir)
        self.assertEqual(obs, [fp])

    def test_split_mapping_file_multiple(self):
        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        fd, fp = mkstemp(suffix='_map.txt')
        close(fd)
        self._clean_up_files.append(fp)

        with open(fp, 'w') as f:
            f.write(MAPPING_FILE_MULT)

        obs = split_mapping_file(fp, out_dir)
        exp = [join(out_dir, 'prefix_1_mapping_file.txt'),
               join(out_dir, 'prefix_2_mapping_file.txt')]
        self.assertItemsEqual(obs, exp)

        obs = sorted(obs)
        with open(obs[0], "U") as f:
            self.assertEqual(f.read(), EXP_MAPPING_FILE_1)

        with open(obs[1], "U") as f:
            self.assertEqual(f.read(), EXP_MAPPING_FILE_2)

    def test_generate_demux_file(self):
        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)
        with open(join(out_dir, 'seqs.fastq'), "w") as f:
            f.write(DEMUX_SEQS)

        obs_fp = generate_demux_file(out_dir)

        exp_fp = join(out_dir, 'seqs.demux')
        self.assertEqual(obs_fp, exp_fp)
        self.assertTrue(exists(exp_fp))

    def test_generate_demux_file_error(self):
        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)
        with self.assertRaises(ValueError):
            generate_demux_file(out_dir)

    def test_generate_demux_file_empty(self):
        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)
        with open(join(out_dir, 'seqs.fastq'), "w") as f:
            f.write('')

        with self.assertRaises(ValueError):
            generate_demux_file(out_dir)

    def test_generate_artifact_info(self):
        obs = generate_artifact_info("/sl/output/")
        fps = [("/sl/output/seqs.fna", "preprocessed_fasta"),
               ("/sl/output/seqs.fastq", "preprocessed_fastq"),
               ("/sl/output/seqs.demux", "preprocessed_demux"),
               ("/sl/output/split_library_log.txt", "log")]
        exp = [ArtifactInfo('demultiplexed', 'Demultiplexed', fps)]
        self.assertEqual(obs, exp)


DEMUX_SEQS = """@a_1 orig_bc=abc new_bc=abc bc_diffs=0
xyz
+
ABC
@b_1 orig_bc=abw new_bc=wbc bc_diffs=4
qwe
+
DFG
@b_2 orig_bc=abw new_bc=wbc bc_diffs=4
qwe
+
DEF
"""

MAPPING_FILE_SINGLE = (
    "#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tDescription\n"
    "Sample1\tGTCCGCAAGTTA\tGTGCCAGCMGCCGCGGTAA\tTGP test\n"
    "Sample2\tCGTAGAGCTCTC\tGTGCCAGCMGCCGCGGTAA\tTGP test\n"
)

MAPPING_FILE_MULT = (
    "#SampleID\tBarcodeSequence\tLinkerPrimerSequence\trun_prefix\t"
    "Description\n"
    "Sample1\tGTCCGCAAGTTA\tGTGCCAGCMGCCGCGGTAA\tprefix_1\tTGP øtest\n"
    "Sample2\tCGTAGAGCTCTC\tGTGCCAGCMGCCGCGGTAA\tprefix_2\tTGP øtest\n"
)

EXP_MAPPING_FILE_1 = (
    "#SampleID\tBarcodeSequence\tLinkerPrimerSequence\trun_prefix\t"
    "Description\n"
    "Sample1\tGTCCGCAAGTTA\tGTGCCAGCMGCCGCGGTAA\tprefix_1\tTGP øtest\n"
)

EXP_MAPPING_FILE_2 = (
    "#SampleID\tBarcodeSequence\tLinkerPrimerSequence\trun_prefix\t"
    "Description\n"
    "Sample2\tCGTAGAGCTCTC\tGTGCCAGCMGCCGCGGTAA\tprefix_2\tTGP øtest\n"
)

if __name__ == '__main__':
    main()
