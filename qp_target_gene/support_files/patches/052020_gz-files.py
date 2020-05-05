# We need to gz all the preprocessed_fasta/preprocessed_fastq from the
# Demultiplexed artifacts. This needs to be ran on the main Qiita install
from os.path import getsize

from qiita_db.study import Study
from qiita_db.processing_job import _system_call as system_call
from qiita_db.sql_connection import TRN
from qiita_db.util import compute_checksum


artifacts = [a for s in Study.iter() for a in
             s.artifacts(artifact_type='Demultiplexed')]

len_artifacts = len(artifacts)
print("===> We are going to process %d artifacts" % len_artifacts)

sql = """UPDATE qiita.filepath SET (filepath, checksum, fp_size) = (%s, %s, %s)
         WHERE filepath_id = %s"""
command = "pigz -p20 %s"

for i, a in enumerate(artifacts):
    if (i+1) % 10:
        print("Processing artifact: %d/%d" % ((i+1), len_artifacts))

    fps_to_gz = [fp for fp in a.filepaths if fp['fp_type'] in (
        'preprocessed_fastq', 'preprocessed_fasta')
        and not fp['fp'].endswith('.gz')]

    for fp in fps_to_gz:
        with TRN:
            pout, perr, rv = system_call(command % fp['fp'])
            if rv != 0:
                raise ValueError("Erro: %s -- %s" % (pout, perr))
            new_fp = '%s.gz' % fp['fp']
            checksum = compute_checksum(new_fp)
            size = getsize(new_fp)
            TRN.add(sql, [new_fp, checksum, size, fp['fp_id']])
            TRN.execute()
