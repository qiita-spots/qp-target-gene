# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from unittest import main
from os.path import isdir, exists, join
from os import remove, close, mkdir
from shutil import rmtree
from tempfile import mkstemp, mkdtemp
from json import dumps
from functools import partial
from glob import glob

from qiita_client import ArtifactInfo
from qiita_client.testing import PluginTestCase

from qp_target_gene.pick_otus import (
    write_parameters_file, generate_artifact_info,
    generate_pick_closed_reference_otus_cmd, generate_sortmerna_tgz,
    pick_closed_reference_otus)

CLIENT_ID = '19ndkO3oMKsoChjVVWluF7QkxHRfYhTKSFbAVt8IhK7gZgDaO4'
CLIENT_SECRET = ('J7FfQ7CQdOxuKhQAf1eoGgBAE81Ns8Gu3EKaWFm3IO2JKh'
                 'AmmCWZuabe0O5Mp28s1')


class PickOTUsTests(PluginTestCase):
    def setUp(self):
        self._clean_up_files = []
        self.parameters = {
            'reference-seq': '/databases/gg/13_8/rep_set/97_otus.fasta',
            'reference-tax': '/databases/gg/13_8/taxonomy/97_otu_taxonomy.txt',
            "sortmerna_e_value": 1, "sortmerna_max_pos": 10000,
            "similarity": 0.97, "sortmerna_coverage": 0.97, "threads": 1,
            "input_data": 2}

    def tearDown(self):
        for fp in self._clean_up_files:
            if exists(fp):
                if isdir(fp):
                    rmtree(fp)
                else:
                    remove(fp)

    def test_write_parameters_file(self):
        fd, fp = mkstemp()
        close(fd)
        self._clean_up_files.append(fp)

        write_parameters_file(fp, self.parameters)

        with open(fp, 'U') as f:
            obs = f.read()
        exp = EXP_PARAMS
        self.assertEqual(obs, exp)

    def test_generate_pick_closed_reference_otus_cmd(self):
        output_dir = mkdtemp()
        self._clean_up_files.append(output_dir)
        filepaths = {'preprocessed_fasta': ['/directory/seqs.fna'],
                     'preprocessed_demux': ['/directory/seqs.demux']}

        obs, obs_dir = generate_pick_closed_reference_otus_cmd(
            filepaths, output_dir, self.parameters, True)
        exp = ("pick_closed_reference_otus.py -i /directory/seqs.fna "
               "-r /databases/gg/13_8/rep_set/97_otus.fasta -o {0}/cr_otus "
               "-p {0}/cr_params.txt -t "
               "/databases/gg/13_8/taxonomy/97_otu_taxonomy.txt".format(
                  output_dir))

        self.assertEqual(obs, exp)
        self.assertEqual(obs_dir, join(output_dir, 'cr_otus'))

    def test_generate_sortmerna_tgz(self):
        outdir = mkdtemp()
        self._clean_up_files.append(outdir)
        mkdir(join(outdir, 'sortmerna_picked_otus'))
        self.assertIsNone(generate_sortmerna_tgz(outdir))

    def test_generate_artifact_info(self):
        outdir = mkdtemp()
        self._clean_up_files.append(outdir)
        log_fp = join(outdir, "log_20151204223007.txt")
        with open(log_fp, 'w') as f:
            f.write("\n")
        self._clean_up_files.append(log_fp)

        obs = generate_artifact_info(outdir)
        fps = [(join(outdir, "otu_table.biom"), "biom"),
               (join(outdir, "sortmerna_picked_otus"), "directory"),
               (join(outdir, "sortmerna_picked_otus.tgz"), "tgz"),
               (log_fp, "log")]
        exp = [ArtifactInfo('OTU table', 'BIOM', fps)]
        self.assertEqual(obs, exp)

    def test_pick_closed_reference_otus(self):
        # Create a new job
        data = {'user': 'demo@microbio.me',
                'command': dumps(['QIIMEq2', '1.9.1',
                                  'Pick closed-reference OTUs']),
                'status': 'running',
                'parameters': dumps(self.parameters)}
        job_id = self.qclient.post(
            '/apitest/processing_job/', data=data)['job']

        # These filepaths do not exist in Qiita - create them
        fps = self.qclient.get('/qiita_db/artifacts/2/')['files']
        fasta_fp = fps['preprocessed_fasta'][0]['filepath']
        self.parameters['reference-seq'] = '/tmp/seq.fna'
        self.parameters['reference-tax'] = '/tmp/tax.txt'
        with open(fasta_fp, 'w') as f:
            f.write(READS)
        # self._clean_up_files.append(fasta_fp)
        with open(self.parameters['reference-seq'], 'w') as f:
            f.write(REF_SEQ)
        # self._clean_up_files.append(self.parameters['reference-seq'])
        with open(self.parameters['reference-tax'], 'w') as f:
            f.write(REF_TAX)
        # self._clean_up_files.append(self.parameters['reference-tax'])

        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        obs_success, obs_ainfo, obs_msg = pick_closed_reference_otus(
            self.qclient, job_id, self.parameters, out_dir)
        self.assertEqual(obs_msg, "")
        self.assertTrue(obs_success)
        path_builder = partial(join, out_dir, 'cr_otus')
        log_fp = glob(path_builder("log_*.txt"))[0]
        fps = [(path_builder("otu_table.biom"), "biom"),
               (path_builder("sortmerna_picked_otus"), "directory"),
               (path_builder("sortmerna_picked_otus.tgz"), "tgz"),
               (log_fp, "log")]
        exp_ainfo = [ArtifactInfo('OTU table', 'BIOM', fps)]
        self.assertEqual(obs_ainfo, exp_ainfo)


EXP_PARAMS = """pick_otus:otu_picking_method\tsortmerna
pick_otus:sortmerna_max_pos\t10000
pick_otus:similarity\t0.97
pick_otus:sortmerna_coverage\t0.97
pick_otus:threads\t1
"""

READS = """>1001.SKB1_0 orig_bc=TAACTTGCGGAC new_bc=TAACTTGCGGAC bc_diffs=0
TACGTAGGTGGCAAGCGTTGTCCGGAATTATTGGGCGTAAAGCGCGCGCAGGCGGTCCTTTAAGTCTGATGTGAAAGC\
CCACGGCTTAACCGTGGAGGGTCATTGGAAACTGGAGGACTTGAGTACAGAAGAGGAGAGAGGAATTCCACGT
>1001.SKB1_1 orig_bc=TAACTTGCGGAC new_bc=TAACTTGCGGAC bc_diffs=0
TACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGTTGTATAAGTCAGTGCTGAAATA\
TCCCGGCTTAACCGGGAGGGTGGCATTGATACTGCGGGGCTTGAGAACGGGTGAGGTAGGCGGAATTGACGGT
>1001.SKB1_2 orig_bc=TAACTTGCGGAC new_bc=TAACTTGCGGAC bc_diffs=0
TACGTAGGGGGCAAGCGTTGTCCGGAATTATTGGGCGTAAAGCGCGCGCAGGCGGTTCTTTAAGTCTGGTGTTTAAAC\
CCGAGGCTCAACTTCGGGTCGCACTGGAAACTGGTGAACTTGCGGGCAGATGAGGAAAGCGGAATTCCACGTG
>1001.SKB1_3 orig_bc=TAACTAGCGGAC new_bc=TAACTTGCGGAC bc_diffs=1
TACGAAGGGTGCAAGCGTTACTCGGAATTACTGGGCGTAAAGCGTGCGTAGGCGGTTTGTTAAGTCTGATGTGAAAGC\
CCTGGGCTCAACCTGAGAATGGCATTGGATACTGTCAGTCTAGAGTGCGGTATAGGCAAGCGGAATTCCTGGT
>1001.SKB1_4 orig_bc=TAACTTGCGAAC new_bc=TAACTTGCGGAC bc_diffs=1
GACGTAGGGGCCGAGCGTTGTCCGGAGTTACTGGGCGTAAAGCGCGCGCAGGCGGATCAGCGCATCGTCGGTGAAAGC\
CCCCCGCTCAACGGGGGAGGGTCCGGCGAGACGGCTGGGCTGGAGGCAGGCAGAGGCGAGTGGTATTCCAGGT
>1001.SKB1_5 orig_bc=TAACTTGCGGAC new_bc=TAACTTGCGGAC bc_diffs=0
TACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGTCCTTTAAGTCAGTGCTGAAATA\
CTCCAGCTTAACTGGAGGGGTGGCATTGATACTGGGGGACTTGAATGAAGTCGAGGTAGGCGGACTTGACGGG
>1001.SKB1_6 orig_bc=TAACTTGCGGAC new_bc=TAACTTGCGGAC bc_diffs=0
TACGGAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGCGGTTATCTAAGCTAGATGTGAAATC\
CCCGGGCTTAACCTGGGAATTGCATTTAGACCTGGATGGCTAGAGTATGGGAGAGGAGTGTGGCATTTCAGGT
>1001.SKB1_7 orig_bc=TAACTTGCGGAC new_bc=TAACTTGCGGAC bc_diffs=0
TACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCTTAAAGCGTGCGTAGTCGGTTATTCAAGTCGGGGGTGAAAGC\
CCCGGGCTCAACCTGGGAATTGCATTCGATACTGTTTAGCTAGAGTTCGGCAGAGGGAAGTGGAATTTCCGGT
>1001.SKB1_8 orig_bc=TAACTTGCGGAC new_bc=TAACTTGCGGAC bc_diffs=0
TACGTAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGTGCGCAGGCGGTTTTGTAAGACCGATGTGAAATC\
CCCGGGCTTAACCTGGGAACTGCATTGGTGCCTGCAAGGCTTGAGTGTGTCAGAGGGAGGTGGAATTCCGCGT
>1001.SKB1_9 orig_bc=TAACTTGCGGAC new_bc=TAACTTGCGGAC bc_diffs=0
TACGAAGGGGACTAGCGTTGTTCGGAATCACTGGGCGTAAAGCGCACGTAGGCGGATATGTCAGTCAGGGGTGAAATC\
CCGGGGCTCAACCTCGGAACTGCCTTTGATACAGCGTCTCTTGAGTCCGATAGAGGCGGGTGGCATTCCTAGT
"""


REF_TAX = """367523\tk__Bacteria; p__Bacteroidetes; c__Flavobacteriia; \
o__; f__; g__; s__
187144\tk__Bacteria; p__Firmicutes; c__Clostridia; o__; f__; g__; s__
836974\tk__Bacteria; p__Cyanobacteria; c__Chloroplast; o__; f__; g__; s__
310669\tk__Bacteria; p__Firmicutes; c__Clostridia; o__; f__; g__; s__"""


REF_SEQ = """>367523
TACGTAGGTGGCAAGCGTTGTCCGGAATTATTGGGCGTAAAGCGCGCGCAGGCGGTCCTTTAAGTCTGATGTGAAAGC\
CCACGGCTTAACCGTGGAGGGTCATTGGAAACTGGAGGACTTGAGTACAGAAGAGGAGAGAGGAATTCCACGT
>187144
TACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGTTGTATAAGTCAGTGCTGAAATA\
TCCCGGCTTAACCGGGAGGGTGGCATTGATACTGCGGGGCTTGAGAACGGGTGAGGTAGGCGGAATTGACGGT
>836974
TACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCTTAAAGCGTGCGTAGTCGGTTATTCAAGTCGGGGGTGAAAGC\
CCCGGGCTCAACCTGGGAATTGCATTCGATACTGTTTAGCTAGAGTTCGGCAGAGGGAAGTGGAATTTCCGGT
>310669
TACGAAGGGGACTAGCGTTGTTCGGAATCACTGGGCGTAAAGCGCACGTAGGCGGATATGTCAGTCAGGGGTGAAATC\
CCGGGGCTCAACCTCGGAACTGCCTTTGATACAGCGTCTCTTGAGTCCGATAGAGGCGGGTGGCATTCCTAGT
"""


if __name__ == '__main__':
    main()
