# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from os.path import join, basename
import re

import pandas as pd

from qiita_client.util import system_call
from .util import (get_artifact_information, split_mapping_file,
                   generate_demux_file, generate_artifact_info)


def generate_parameters_string(parameters):
    """Generates the parameters string from the parameters dictionary

    Parameters
    ----------
    parameters : dict
        The parameter values, keyed by parameter name

    Returns
    -------
    str
        A string with the parameters to the CLI call
    """
    flag_params = ['rev_comp_barcode', 'rev_comp_mapping_barcodes', 'rev_comp']
    str_params = ['max_bad_run_length', 'min_per_read_length_fraction',
                  'sequence_max_n', 'phred_quality_threshold', 'barcode_type',
                  'max_barcode_errors']
    result = ["--%s %s" % (sp, parameters[sp]) for sp in str_params]
    if parameters['phred_offset'] != 'auto':
        result.append("--phred_offset %s" % parameters['phred_offset'])
    for fp in flag_params:
        if parameters[fp]:
            result.append("--%s" % fp)
    return ' '.join(result)


def get_sample_names_by_run_prefix(mapping_file):
    """Generates a dictionary of run_prefix and sample names

    Parameters
    ----------
    mapping_file : str
        The mapping file

    Returns
    -------
    dict
        Dict mapping run_prefix to sample id

    Raises
    ------
    ValueError
        If there is more than 1 sample per run_prefix
    """
    qiime_map = pd.read_csv(mapping_file, delimiter='\t', dtype=str,
                            encoding='utf-8')
    qiime_map.set_index('#SampleID', inplace=True)

    samples = {}
    errors = []
    for prefix, df in qiime_map.groupby('run_prefix'):
        len_df = len(df)
        if len_df != 1:
            errors.append('%s has %d samples (%s)' % (prefix, len_df,
                                                      ', '.join(df.index)))
        else:
            samples[prefix] = df.index.values[0]

    if errors:
        raise ValueError("You have run_prefix values with multiple "
                         "samples: %s" % ' -- '.join(errors))

    return samples


def generate_per_sample_fastq_command(forward_seqs, reverse_seqs, barcode_fps,
                                      mapping_file, output_dir, params_str):
    """Generates the per-sample FASTQ split_libraries_fastq.py command

    Parameters
    ----------
    forward_seqs : list of str
        The list of forward seqs filepaths
    reverse_seqs : list of str
        The list of reverse seqs filepaths
    barcode_fps : list of str
        The list of barcode filepaths
    mapping_file : str
        The path to the mapping file
    output_dir : str
        The path to the split libraries output directory
    params_str : str
        The string containing the parameters to pass to
        split_libraries_fastq.py

    Returns
    -------
    str
        The CLI to execute

    Raises
    ------
    ValueError
        - If barcode_fps is not an empty list
        - If there are run prefixes in the mapping file that do not match
        the sample names
    """
    if barcode_fps:
        raise ValueError('per_sample_FASTQ can not have barcodes: %s'
                         % (', '.join(basename(b) for b in barcode_fps)))
    sn_by_rp = get_sample_names_by_run_prefix(mapping_file)
    samples = []
    errors = []
    for fname in forward_seqs:
        fn = basename(fname)

        # removing extentions: fastq or fastq.gz
        if 'fastq' in fn.lower().rsplit('.', 2):
            f = fn[:fn.lower().rindex('.fastq')]
        else:
            f = fn
        m = [v for v in sn_by_rp if f.startswith(v)]

        # removing study_id, in case it's present
        if re.match(r"^[0-9]+\_.*", f):
            f = basename(fn).split('_', 1)[1]
        mi = [v for v in sn_by_rp if f.startswith(v)]

        # the matches is the largest between m/mi, if they are the same size
        # we are gonna use m
        matches = m if len(m) > len(mi) else mi

        if matches:
            len_matches = len(matches)
            if len_matches != 1:
                errors.append('%s has %s matches.' % (f, len_matches))
            for m in matches:
                samples.append(sn_by_rp[m])
                del sn_by_rp[m]
        else:
            errors.append('%s has NO matches' % f)

    if errors:
        raise ValueError('Errors found:\n%s' % '\n'.join(errors))

    cmd = str("split_libraries_fastq.py --store_demultiplexed_fastq "
              "-i %s --sample_ids %s -o %s %s"
              % (','.join(forward_seqs), ','.join(samples),
                 output_dir, params_str))
    return cmd


def generate_split_libraries_fastq_cmd(filepaths, mapping_file, atype,
                                       out_dir, parameters):
    """Generates the split_libraries_fastq.py command

    Parameters
    ----------
    filepaths : dict of {str: list of str}
        The artifact filepaths keyed by type
    mapping_file : str
        The artifact QIIME-compliant mapping file
    atype : str
        The artifact type
    out_dir : str
        The job output directory

    Returns
    -------
    str
        The CLI to execute

    Raises
    ------
    ValueError
        If the number of barcode files and the number of sequence files do not
        match
    """
    forward_seqs = filepaths.get('raw_forward_seqs', [])
    reverse_seqs = filepaths.get('raw_reverse_seqs', [])
    barcode_fps = filepaths.get('raw_barcodes', [])

    # We need to sort the filepaths to make sure that each lane's file is in
    # the same order, so they match when passed to split_libraries_fastq.py
    # All files should be prefixed with run_prefix, so the ordering is
    # ensured to be correct
    forward_seqs = sorted(forward_seqs)
    reverse_seqs = sorted(reverse_seqs)
    barcode_fps = sorted(barcode_fps)

    output_dir = join(out_dir, "sl_out")

    params_str = generate_parameters_string(parameters)

    if atype == "per_sample_FASTQ":
        cmd = generate_per_sample_fastq_command(
            forward_seqs, reverse_seqs, barcode_fps, mapping_file,
            output_dir, params_str)
    else:
        if len(barcode_fps) != len(forward_seqs):
            raise ValueError("The number of barcode files and the number of "
                             "sequence files should match: %d != %s"
                             % (len(barcode_fps), len(forward_seqs)))

        map_out_dir = join(out_dir, 'mappings')
        mapping_files = sorted(split_mapping_file(mapping_file, map_out_dir))

        cmd = str("split_libraries_fastq.py --store_demultiplexed_fastq -i %s "
                  "-b %s -m %s -o %s %s"
                  % (','.join(forward_seqs), ','.join(barcode_fps),
                     ','.join(mapping_files), output_dir, params_str))

    return cmd, output_dir


def split_libraries_fastq(qclient, job_id, parameters, out_dir):
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
    """
    # Step 1 get the rest of the information need to run split libraries
    qclient.update_job_step(job_id, "Step 1 of 4: Collecting information")
    artifact_id = parameters['input_data']
    filepaths, mapping_file, atype = get_artifact_information(
        qclient, artifact_id)

    # Step 2 generate the split libraries fastq command
    qclient.update_job_step(job_id, "Step 2 of 4: Generating command")
    command, sl_out = generate_split_libraries_fastq_cmd(
        filepaths, mapping_file, atype, out_dir, parameters)

    # Step 3 execute split libraries
    qclient.update_job_step(
        job_id, "Step 3 of 4: Executing demultiplexing and quality control")
    std_out, std_err, return_value = system_call(command)
    if return_value != 0:
        raise RuntimeError(
            "Error processing files:\nStd output: %s\n Std error:%s"
            % (std_out, std_err))

    # Step 4 generate the demux file
    qclient.update_job_step(job_id, "Step 4 of 4: Generating demux file")
    generate_demux_file(sl_out)

    artifacts_info = generate_artifact_info(sl_out)

    return True, artifacts_info, ""
