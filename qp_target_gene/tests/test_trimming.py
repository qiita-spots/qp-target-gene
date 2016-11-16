# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from unittest import main
from os.path import isdir, exists, join
from os import remove, close
from shutil import rmtree, copyfile
from tempfile import mkstemp, mkdtemp
from json import dumps
from functools import partial

from qiita_client import ArtifactInfo
from qiita_client.testing import PluginTestCase

from qp_target_gene.trimming import (trimming)
from qp_target_gene import plugin


CLIENT_ID = '19ndkO3oMKsoChjVVWluF7QkxHRfYhTKSFbAVt8IhK7gZgDaO4'
CLIENT_SECRET = ('J7FfQ7CQdOxuKhQAf1eoGgBAE81Ns8Gu3EKaWFm3IO2JKh'
                 'AmmCWZuabe0O5Mp28s1')


class TrimmingTest(PluginTestCase):
    def setUp(self):
        plugin("https://localhost:21174", 'register', 'ignored')
        self._clean_up_files = []

    def tearDown(self):
        for fp in self._clean_up_files:
            if exists(fp):
                if isdir(fp):
                    rmtree(fp)
                else:
                    remove(fp)

    def test_trimming(self):
        # generating filepaths
        fd, fp = mkstemp(suffix='_seqs.demux')
        close(fd)
        self._clean_up_files.append(fp)
        copyfile('support_files/filtered_5_seqs.demux', fp)

        # inserting new prep template
        prep_info_dict = {
            'SKB7.640196': {'description': 'SKB7'},
            'SKB8.640193': {'description': 'SKB8'}
        }
        data = {'prep_info': dumps(prep_info_dict),
                # magic #1 = testing study
                'study': 1,
                'data_type': '16S'}
        pid = self.qclient.post('/apitest/prep_template/', data=data)['prep']

        # inserting artifacts
        data = {
            'filepaths': dumps([(fp, 'preprocessed_demux')]),
            'type': "Demultiplexed",
            'name': "New demultiplexed artifact",
            'prep': pid}
        aid = self.qclient.post('/apitest/artifact/', data=data)['artifact']

        params = {'input_data': aid, 'length': 50}
        data = {'user': 'demo@microbio.me',
                'command': dumps(['QIIME', '1.9.1', 'Trimming']),
                'status': 'running', 'parameters': dumps(params)}
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']

        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        success, ainfo, msg = trimming(self.qclient, jid, params, out_dir)
        self.assertTrue(success)
        pb = partial(join, out_dir)
        exp_ainfo = [
            ArtifactInfo(
                'Trimmed Demultiplexed', 'Demultiplexed',
                [(pb('seqs.fna'), 'preprocessed_fasta'),
                 (pb('seqs.fastq'), 'preprocessed_fastq'),
                 (pb('seqs.demux'), 'preprocessed_demux')])]
        self.assertEqual(ainfo, exp_ainfo)
        self.assertEqual(msg, "")


if __name__ == '__main__':
    main()
