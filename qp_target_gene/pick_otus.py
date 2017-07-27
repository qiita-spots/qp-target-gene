# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from os.path import join, basename
from functools import partial
from glob import glob
from tarfile import open as taropen

from qiita_client import ArtifactInfo

from qp_target_gene.util import system_call


def write_parameters_file(fp, parameters):
    """Write the QIIME parameters file

    Parameters
    ----------
    fp : str
        The paramters file filepath
    parameters : dict
        The commands parameters, keyed by parameter name
    """
    params = ['sortmerna_max_pos', 'similarity', 'sortmerna_coverage',
              'threads']
    with open(fp, 'w') as f:
        f.write("pick_otus:otu_picking_method\tsortmerna\n")
        for p in params:
            f.write("pick_otus:%s\t%s\n" % (p, parameters[p]))


def generate_pick_closed_reference_otus_cmd(filepaths, out_dir, parameters):
    """Generates the pick_closed_reference_otus.py command

    Parameters
    ----------
    filepaths : list of (str, str)
        The artifact's filepaths and their types
    out_dir : str
        The job output directory
    parameters : dict
        The command's parameters, keyed by parameter name

    Returns
    -------
    str, str
        The pick_closed_reference_otus.py command
        The output directory
    """
    # It should be only a single preprocessed fasta file
    seqs_fp = filepaths['preprocessed_fasta'][0]

    output_dir = join(out_dir, 'cr_otus')
    param_fp = join(out_dir, 'cr_params.txt')

    reference_fp = parameters.pop('reference-seq')
    taxonomy_fp = parameters.pop('reference-tax')

    write_parameters_file(param_fp, parameters)

    cmd = str("pick_closed_reference_otus.py -i %s -r %s -o %s -p %s -t %s"
              % (seqs_fp, reference_fp, output_dir, param_fp, taxonomy_fp))
    return cmd, output_dir


def generate_sortmerna_tgz(out_dir):
    """Generates the sortmerna failures tgz command

    Parameters
    ----------
    out_dir : str
        The job output directory

    Returns
    -------
    str
        The sortmerna failures tgz command
    """
    to_tgz = join(out_dir, 'sortmerna_picked_otus')
    tgz = to_tgz + '.tgz'
    with taropen(tgz, "w:gz") as tar:
        tar.add(to_tgz, arcname=basename(to_tgz))


def generate_artifact_info(pick_out):
    """Creates the artifact information to attach to the payload

    Parameters
    ----------
    pick_out : str
        Path to the pick otus directory

    Returns
    -------
    list
        The artifacts information
    """
    path_builder = partial(join, pick_out)
    filepaths = [(path_builder('otu_table.biom'), 'biom'),
                 (path_builder('sortmerna_picked_otus'), 'directory'),
                 (path_builder('sortmerna_picked_otus.tgz'), 'tgz'),
                 (glob(path_builder('log_*.txt'))[0], 'log')]
    return [ArtifactInfo('OTU table', 'BIOM', filepaths)]


def pick_closed_reference_otus(qclient, job_id, parameters, out_dir):
    """Run split libraries fastq with the given parameters

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

    Raises
    ------
    ValueError
        If there is any error gathering the information from the server
    """
    qclient.update_job_step(job_id, "Step 1 of 4: Collecting information")
    artifact_id = parameters['input_data']
    a_info = qclient.get("/qiita_db/artifacts/%s/" % artifact_id)
    fps = a_info['files']

    qclient.update_job_step(job_id, "Step 2 of 4: Generating command")
    command, pick_out = generate_pick_closed_reference_otus_cmd(
        fps, out_dir, parameters)

    qclient.update_job_step(job_id, "Step 3 of 4: Executing OTU picking")
    std_out, std_err, return_value = system_call(command)
    if return_value != 0:
        error_msg = ("Error running OTU picking:\nStd out: %s\nStd err: %s"
                     % (std_out, std_err))
        return False, None, error_msg

    qclient.update_job_step(job_id,
                            "Step 4 of 4: Generating tgz sortmerna folder")
    try:
        generate_sortmerna_tgz(pick_out)
    except Exception as e:
        error_msg = ("Error while tgz failures:\nError: %s" % str(e))
        return False, None, error_msg

    artifacts_info = generate_artifact_info(pick_out)

    return True, artifacts_info, ""
