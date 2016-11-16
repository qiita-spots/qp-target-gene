# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from os.path import join
from functools import partial
from h5py import File

from qiita_client import ArtifactInfo

from qp_target_gene.split_libraries.util import generate_demux_file
from qiita_files.demux import fetch
from qiita_files.format.fasta import format_fasta_record
from qiita_files.format.fastq import format_fastq_record


def generate_trimming(filepaths, out_dir, parameters):
    """Generate the trimming of the filepaths

    Parameters
    ----------
    filepaths : list of str
        The demux filepaths
    out_dir : str
        The job output directory
    parameters : dict
        The command's parameters, keyed by parameter name

    Returns
    -------
    str
        The pick_closed_reference_otus.py command
        The output directory
    """
    length = parameters['length']

    id_fmt = (b"%(sample)s_%(idx)d orig_bc=%(bc_ori)s new_bc=%(bc_cor)s "
              b"bc_diffs=%(bc_diff)d")
    pd = partial(join, out_dir)
    ffp = pd('seqs.fna')
    qfp = pd('seqs.fastq')
    for f in filepaths:
        with open(ffp, 'w') as ffh, open(qfp, 'w') as qfh, File(f) as fh:
            for samp, idx, seq, qual, bc_ori, bc_cor, bc_err in fetch(fh):
                seq_id = id_fmt % {b'sample': samp, b'idx': idx,
                                   b'bc_ori': bc_ori, b'bc_cor': bc_cor,
                                   b'bc_diff': bc_err}
                ffh.write(
                    format_fasta_record(seq_id, seq[:length], qual[:length]))
                qfh.write(
                    format_fastq_record(seq_id, seq[:length], qual[:length]))


def trimming(qclient, job_id, parameters, out_dir):
    """Run trimming over the given parameters

    Parameters
    ----------
    qclient : tgp.qiita_client.QiitaClient
        The Qiita server client
    job_id : str
        The job id
    parameters : dict
        The parameter values to run split libraries
    out_dir : str
        Yhe path to the job's output directory

    Returns
    -------
    bool, list, str
        The results of the job
    """
    qclient.update_job_step(job_id, "Step 1 of 3: Collecting information")
    artifact_id = parameters['input_data']
    a_info = qclient.get("/qiita_db/artifacts/%s/" % artifact_id)
    fps = a_info['files']
    if 'preprocessed_demux' not in fps:
        error_msg = "Artifact doesn't contain a preprocessed demux"
        return False, None, error_msg

    qclient.update_job_step(job_id, "Step 2 of 3: Executing Trimming")
    generate_trimming(fps['preprocessed_demux'], out_dir, parameters)

    qclient.update_job_step(job_id, "Step 3 of 3: Generating new Demuxed")
    generate_demux_file(out_dir)

    pb = partial(join, out_dir)
    ainfo = [
        ArtifactInfo(
            'Trimmed Demultiplexed', 'Demultiplexed',
            [(pb('seqs.fna'), 'preprocessed_fasta'),
             (pb('seqs.fastq'), 'preprocessed_fastq'),
             (pb('seqs.demux'), 'preprocessed_demux')])]

    return True, ainfo, ""
