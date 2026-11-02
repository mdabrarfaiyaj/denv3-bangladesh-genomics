#!/usr/bin/env nextflow
/*
 * denv3-bangladesh-genomics
 * ------------------------------------------------------------------
 * Genomic epidemiology of DENV-3 in Bangladesh, 2017-2025.
 *
 * Every process wraps a script from scripts/ rather than reimplementing it.
 * Those scripts are the analysis; this file only describes how they connect,
 * what each one needs, and what it produces. That means the pipeline and the
 * step-by-step runs cannot drift apart - there is one implementation.
 *
 * Two entry points:
 *   --input_mode raw        start from the NCBI download, run Nextclade
 *                           (needs network: the Nextclade dataset is fetched)
 *   --input_mode nextclade  start from an existing Nextclade result
 *                           (no network; this is what CI and -profile test use)
 *
 * The download itself is deliberately NOT a pipeline step. Fetching 9,815
 * sequences from NCBI on every run would be slow, impolite, and irreproducible
 * because the database changes daily. scripts/00_download.sh does it once and
 * records a checksum, exactly as nf-core pipelines handle reference data.
 */

nextflow.enable.dsl = 2

// ---------------------------------------------------------------- helpers ---

def banner() {
    log.info """
    ==========================================================
     DENV-3 Bangladesh  |  genomic epidemiology
    ==========================================================
     input mode   : ${params.input_mode}
     datasets     : ${params.sets}
     seed         : ${params.seed}   context cap: ${params.context_per_group}
     outgroup     : ${params.outgroup_lineages}
     outdir       : ${params.outdir}
    ==========================================================
    """.stripIndent()
}

/* Each script expects its inputs at fixed relative paths (results/...), which is
 * what makes them runnable on their own. Rather than rewrite six scripts, every
 * process rebuilds that small directory layout from its staged inputs. The
 * staging is explicit and visible in each script block below. */

// --------------------------------------------------------------- processes ---

process PARSE_METADATA {
    label 'light'
    publishDir "${params.outdir}/metadata", mode: 'copy'

    input:
    path report

    output:
    path 'metadata.tsv', emit: metadata
    path 'parse.log'

    script:
    """
    mkdir -p data/raw/denv3_all/ncbi_dataset/data results/metadata
    cp ${report} data/raw/denv3_all/ncbi_dataset/data/data_report.jsonl
    python ${projectDir}/scripts/01_parse_metadata.py | tee parse.log
    cp results/metadata/metadata.tsv .
    """
}

process NEXTCLADE_DATASET {
    label 'light'
    publishDir "${params.outdir}/nextclade", mode: 'copy'

    output:
    path 'dataset', emit: dataset

    script:
    """
    nextclade dataset get \\
        --name '${params.nextclade_name}' \\
        --tag '${params.nextclade_tag}' \\
        --output-dir dataset
    """
}

process NEXTCLADE_RUN {
    label 'medium'
    publishDir "${params.outdir}/nextclade", mode: 'copy'

    input:
    path fasta
    path dataset

    output:
    path 'nextclade.tsv', emit: tsv
    path 'nextclade.aligned.fasta', emit: aligned

    script:
    """
    nextclade run \\
        --input-dataset ${dataset} \\
        --output-tsv nextclade.tsv \\
        --output-fasta nextclade.aligned.fasta \\
        --jobs ${task.cpus} \\
        ${fasta}
    """
}

process SUMMARISE_CLADES {
    label 'light'
    publishDir "${params.outdir}/metadata", mode: 'copy', pattern: '*.tsv'

    input:
    path nextclade_tsv
    path metadata

    output:
    path 'metadata_clades.tsv', emit: clades
    path 'summarise.log'

    script:
    """
    mkdir -p results/nextclade results/metadata
    cp ${nextclade_tsv} results/nextclade/nextclade.tsv
    cp ${metadata} results/metadata/metadata.tsv
    python ${projectDir}/scripts/02b_summarise_nextclade.py | tee summarise.log
    cp results/metadata/metadata_clades.tsv .
    """
}

process LINEAGE_GEOGRAPHY {
    label 'light'
    publishDir "${params.outdir}/metadata", mode: 'copy'

    input:
    path clades

    output:
    path 'lineage_geography.tsv'
    path 'geography.log'

    script:
    """
    mkdir -p results/metadata
    cp ${clades} results/metadata/metadata_clades.tsv
    python ${projectDir}/scripts/03_lineage_geography.py | tee geography.log
    cp results/metadata/lineage_geography.tsv .
    """
}

process CURATION_AUDIT {
    label 'light'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path metadata
    path clades
    path nextclade_tsv

    output:
    path 'audit/*', optional: true
    path 'audit.log'

    script:
    """
    mkdir -p results/metadata results/nextclade
    cp ${metadata} results/metadata/metadata.tsv
    cp ${clades} results/metadata/metadata_clades.tsv
    cp ${nextclade_tsv} results/nextclade/nextclade.tsv
    python ${projectDir}/scripts/08_curation_audit.py | tee audit.log
    mkdir -p audit && cp -r results/audit/* audit/ 2>/dev/null || true
    """
}

process MUTATION_PROFILE {
    label 'light'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path nextclade_tsv
    path clades
    path gff

    output:
    path 'mutations/*', optional: true
    path 'mutations.log'

    script:
    """
    mkdir -p results/nextclade/dataset results/metadata
    cp ${nextclade_tsv} results/nextclade/nextclade.tsv
    cp ${clades} results/metadata/metadata_clades.tsv
    cp ${gff} results/nextclade/dataset/genome_annotation.gff3
    python ${projectDir}/scripts/05_mutation_profile.py | tee mutations.log
    mkdir -p mutations && cp -r results/mutations/* mutations/ 2>/dev/null || true
    """
}

process SELECT_SEQUENCES {
    label 'medium'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path clades
    path aligned
    path gff

    output:
    path 'datasets', emit: datasets
    path 'datasets/selection_log.tsv', emit: log

    script:
    """
    mkdir -p results/metadata results/nextclade/dataset
    cp ${clades} results/metadata/metadata_clades.tsv
    cp ${aligned} results/nextclade/nextclade.aligned.fasta
    cp ${gff} results/nextclade/dataset/genome_annotation.gff3
    # tee opens its output file immediately, before the script creates the
    # directory, so the directory has to exist first
    mkdir -p datasets
    python ${projectDir}/scripts/04_select_sequences.py \\
        --outdir datasets \\
        --seed ${params.seed} \\
        --context-per-group ${params.context_per_group} \\
        --min-cov-wg ${params.min_cov_wg} \\
        --min-cov-eg ${params.min_cov_eg} \\
        --outgroup-lineages '${params.outgroup_lineages}' \\
        --outgroup-max ${params.outgroup_max} \\
        | tee datasets/select.log
    """
}

process BUILD_TREE {
    label 'heavy'
    tag "${set}"
    // the output pattern already names trees/, so publish at the top level or
    // the files land in results/trees/trees/
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(set), path(datasets)

    output:
    tuple val(set), path("trees/${set}.treefile"), emit: tree
    path "trees/${set}.*"

    script:
    """
    DATA=${datasets} OUT=trees SETS=${set} \\
      THREADS=${task.cpus} MODEL='${params.iqtree_model}' BOOT=${params.bootstrap} \\
      bash ${projectDir}/scripts/06_build_trees.sh
    """
}

process INSPECT_TREES {
    label 'light'
    // publish the log too: the pattern was '*.tsv' and silently dropped
    // inspect.log, which is where the support values are reported
    publishDir "${params.outdir}/trees", mode: 'copy'

    input:
    path datasets
    path trees

    output:
    path '*_key_nodes.tsv', optional: true
    path 'inspect.log'

    script:
    """
    mkdir -p results/trees
    cp -r ${datasets} results/datasets
    cp ${trees} results/trees/
    python ${projectDir}/scripts/06b_inspect_trees.py | tee inspect.log
    cp results/trees/*_key_nodes.tsv . 2>/dev/null || true
    """
}

process CLOCK_CHECK {
    label 'medium'
    publishDir "${params.outdir}/clock", mode: 'copy'

    input:
    path datasets
    path trees

    output:
    path 'summary.tsv', optional: true
    path 'clock.log'

    script:
    """
    mkdir -p results/trees
    cp -r ${datasets} results/datasets
    cp ${trees} results/trees/
    python ${projectDir}/scripts/07_clock_check.py | tee clock.log
    cp results/clock/summary.tsv . 2>/dev/null || true
    """
}

process DATE_AND_INTRODUCTIONS {
    label 'medium'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path datasets
    path trees

    output:
    path 'dating/*', optional: true
    path 'dating.log'

    script:
    """
    mkdir -p results/trees
    cp -r ${datasets} results/datasets
    cp ${trees} results/trees/
    python ${projectDir}/scripts/11_date_and_introductions.py \\
        --sets ${params.sets.replaceAll(',', ' ')} | tee dating.log
    mkdir -p dating && cp -r results/dating/* dating/ 2>/dev/null || true
    """
}

// ---------------------------------------------------------------- workflow ---

workflow {

    banner()

    if (params.input_mode == 'raw') {
        report_ch  = Channel.fromPath(params.data_report, checkIfExists: true)
        fasta_ch   = Channel.fromPath(params.genomic_fasta, checkIfExists: true)

        PARSE_METADATA(report_ch)
        NEXTCLADE_DATASET()
        NEXTCLADE_RUN(fasta_ch, NEXTCLADE_DATASET.out.dataset)

        metadata_ch  = PARSE_METADATA.out.metadata
        nextclade_ch = NEXTCLADE_RUN.out.tsv
        aligned_ch   = NEXTCLADE_RUN.out.aligned
        gff_ch       = NEXTCLADE_DATASET.out.dataset.map { it.resolve('genome_annotation.gff3') }
    }
    else if (params.input_mode == 'nextclade') {
        metadata_ch  = Channel.fromPath(params.metadata, checkIfExists: true)
        nextclade_ch = Channel.fromPath(params.nextclade_tsv, checkIfExists: true)
        aligned_ch   = Channel.fromPath(params.aligned_fasta, checkIfExists: true)
        gff_ch       = Channel.fromPath(params.gff, checkIfExists: true)
    }
    else {
        error "--input_mode must be 'raw' or 'nextclade', got '${params.input_mode}'"
    }

    SUMMARISE_CLADES(nextclade_ch, metadata_ch)
    clades_ch = SUMMARISE_CLADES.out.clades

    LINEAGE_GEOGRAPHY(clades_ch)
    CURATION_AUDIT(metadata_ch, clades_ch, nextclade_ch)
    MUTATION_PROFILE(nextclade_ch, clades_ch, gff_ch)

    SELECT_SEQUENCES(clades_ch, aligned_ch, gff_ch)
    datasets_ch = SELECT_SEQUENCES.out.datasets

    // one tree per dataset, in parallel
    sets_ch = Channel.fromList(params.sets.split(/[ ,]+/) as List)
    BUILD_TREE(sets_ch.combine(datasets_ch))

    // collect every tree before the downstream steps, which need them together
    trees_ch = BUILD_TREE.out.tree.map { it[1] }.collect()

    INSPECT_TREES(datasets_ch, trees_ch)
    CLOCK_CHECK(datasets_ch, trees_ch)
    DATE_AND_INTRODUCTIONS(datasets_ch, trees_ch)
}
