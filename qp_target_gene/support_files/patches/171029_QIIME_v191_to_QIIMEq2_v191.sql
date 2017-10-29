-- With the QIIMEq2 release, the parameters for the command "Pick closed-reference
-- OTUs" changed. The parameter "reference" has ben substituted for 2 parameters
-- "reference-tax" and "reference-seq" which contains the path to the actual
-- reference taxonomy file and the reference sequence file

DO $do$
DECLARE
    potu_cmd_id     BIGINT;
    new_potu_cmd_id BIGINT;
    a_info          RECORD;
    parameters      JSON;
    ref_seqs        VARCHAR;
    ref_taxa        VARCHAR;
    pj_id           UUID;
    input_data      VARCHAR;
BEGIN
    SELECT command_id INTO potu_cmd_id
        FROM qiita.software_command sc
            JOIN qiita.software s USING (software_id)
        WHERE s.name = 'QIIME' AND s.version = '1.9.1' AND sc.name = 'Pick closed-reference OTUs';

    SELECT command_id INTO new_potu_cmd_id
        FROM qiita.software_command sc
            JOIN qiita.software s USING (software_id)
        WHERE s.name = 'QIIMEq2' AND s.version = '1.9.1' AND sc.name = 'Pick closed-reference OTUs';

    -- Modify those artifacts and jobs that have been picked against Greengenes
    FOR a_info IN
        SELECT * FROM qiita.artifact WHERE command_id = potu_cmd_id
    LOOP
        IF a_info.command_parameters->>'reference' = '1' THEN
            -- Greengenes
            ref_seqs := '/databases/gg/13_8/rep_set/97_otus.fasta';
            ref_taxa := '/databases/gg/13_8/taxonomy/97_otu_taxonomy.txt';
        ELSIF a_info.command_parameters->>'reference' = '2' THEN
            -- Silva
            ref_seqs := '/projects/qiita_data/reference/silva_119_Silva_119_rep_set97.fna';
            ref_taxa := '/projects/qiita_test_data/reference/silva_119_taxonomy_97_7_levels.txt';
        ELSIF a_info.command_parameters->>'reference' = '3' THEN
            -- UNITE
            ref_seqs := '/projects/qiita_test_data/reference/unite_7_sh_refs_qiime_ver7_97_s_02.03.2015.fasta';
            ref_taxa := '/projects/qiita_test_data/reference/unite_7_sh_taxonomy_qiime_ver7_97_s_02.03.2015.txt';
        ELSE
            -- Unknown reference
            RAISE NOTICE 'Artifact with ID % contains an unknown reference id: %', a_info.artifact_id, a_info.command_parameters->>'reference';
        END IF;

        parameters := ('{"input_data": "' || (a_info.command_parameters->>'input_data')::varchar || '", '
                       '"reference-seq": "' || ref_seqs || '", '
                       '"reference-tax": "' || ref_taxa || '", '
                       '"similarity": "' || (a_info.command_parameters->>'similarity')::varchar || '", '
                       '"sortmerna_coverage": "' || (a_info.command_parameters->>'sortmerna_coverage')::varchar || '", '
                       '"sortmerna_e_value": "' || (a_info.command_parameters->>'sortmerna_e_value')::varchar || '", '
                       '"sortmerna_max_pos": "' || (a_info.command_parameters->>'sortmerna_max_pos')::varchar || '", '
                       '"threads": "' || (a_info.command_parameters->>'threads')::varchar || '"}')::json;

        SELECT processing_job_id INTO pj_id
            FROM qiita.processing_job
                JOIN qiita.artifact_output_processing_job USING (processing_job_id)
            WHERE artifact_id = a_info.artifact_id;

        UPDATE qiita.processing_job
            SET command_parameters = parameters, command_id = new_potu_cmd_id
            WHERE processing_job_id = pj_id;

        UPDATE qiita.artifact
            SET command_parameters = parameters, command_id = new_potu_cmd_id
            WHERE artifact_id = a_info.artifact_id;

    END LOOP;
END $do$
