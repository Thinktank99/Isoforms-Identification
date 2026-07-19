# Snakefile
#
# Drop into the same folder as your existing scripts, genome.fna, genomic.gtf,
# and *_1.fastq.gz / *_2.fastq.gz files. Assumes your usual environment
# (fastqc, cutadapt, picard, STAR, Drop-seq tools, R/Seurat) is already
# active on PATH -- activate it the same way env_conda.sh did, before
# running snakemake.
#
# Only the real production path is included (raw fastq -> tag -> align ->
# DGE -> Seurat -> QC/doublets -> integrate -> markers). 03/07/07A/08_relaxed/
# average_mapping.sh/09_summary.sh were one-off benchmarking scripts, not
# part of this DAG -- run them standalone if you need them.
#
# Usage:
#   snakemake -n                    # dry run
#   snakemake --cores 16            # local run
#   snakemake --cores 16 --sbatch ... / --workflow-profile ...  # SLURM

SAMPLES = ["ERR9871746", "ERR9871747", "ERR9871748", "ERR9871749"]
GENOME = "genome.fna"
GTF = "genomic.gtf"


rule all:
    input:
        expand("results/00_qc/{sample}_1_fastqc.html", sample=SAMPLES),
        expand("results/00_qc/{sample}_2_fastqc.html", sample=SAMPLES),
        "results/10_markers/all_markers.tsv",


rule fastqc:
    input:
        r1="{sample}_1.fastq.gz",
        r2="{sample}_2.fastq.gz",
    output:
        "results/00_qc/{sample}_1_fastqc.html",
        "results/00_qc/{sample}_2_fastqc.html",
    threads: 4
    shell:
        """
        mkdir -p results/00_qc
        fastqc -t {threads} -o results/00_qc {input.r1} {input.r2}
        """


rule build_reference:
    input:
        genome=GENOME,
        gtf=GTF,
    output:
        dict="genome.dict",
        refflat="dropseq_ref/genes.refFlat",
        reduced="dropseq_ref/genes.reduced.gtf",
        intervals="dropseq_ref/genes.intervals",
        star_index=directory("star_index"),
    threads: 16
    shell:
        """
        mkdir -p dropseq_ref {output.star_index}
        picard CreateSequenceDictionary R={input.genome} O={output.dict}
        ConvertToRefFlat ANNOTATIONS_FILE={input.gtf} SEQUENCE_DICTIONARY={output.dict} \
            OUTPUT={output.refflat}
        ReduceGtf SEQUENCE_DICTIONARY={output.dict} GTF={input.gtf} OUTPUT={output.reduced}
        CreateIntervalsFiles SEQUENCE_DICTIONARY={output.dict} REDUCED_GTF={output.reduced} \
            PREFIX=my OUTPUT=dropseq_ref
        STAR --runThreadN {threads} \
          --runMode genomeGenerate \
          --genomeDir {output.star_index} \
          --genomeFastaFiles {input.genome} \
          --sjdbGTFfile {input.gtf} \
          --sjdbOverhang 99 \
          --limitGenomeGenerateRAM 75000000000
        """


rule fastq_to_ubam:
    input:
        r1="{sample}_1.fastq.gz",
        r2="{sample}_2.fastq.gz",
    output:
        "results/05_tagged/{sample}.unmapped.bam",
    shell:
        """
        mkdir -p results/05_tagged
        picard FastqToSam \
          F1={input.r1} F2={input.r2} O={output} \
          SM={wildcards.sample} RG=RG_{wildcards.sample} PL=ILLUMINA \
          TMP_DIR=$(mktemp -d)
        """


rule tag_umi:
    input:
        "results/05_tagged/{sample}.unmapped.bam",
    output:
        bam="results/05_tagged/{sample}.xm.bam",
        summary="results/05_tagged/{sample}.xm.summary.txt",
    shell:
        """
        TagBamWithReadSequenceExtended \
          INPUT={input} OUTPUT={output.bam} SUMMARY={output.summary} \
          BASE_RANGE=1-10 BASE_QUALITY=10 BARCODED_READ=2 DISCARD_READ=false \
          TAG_NAME=XM NUM_BASES_BELOW_QUALITY=1
        """


rule tag_cell:
    input:
        "results/05_tagged/{sample}.xm.bam",
    output:
        bam="results/05_tagged/{sample}.xm_xc.bam",
        summary="results/05_tagged/{sample}.xc.summary.txt",
    shell:
        """
        TagBamWithReadSequenceExtended \
          INPUT={input} OUTPUT={output.bam} SUMMARY={output.summary} \
          BASE_RANGE=87-94:49-56:11-18 BASE_QUALITY=10 BARCODED_READ=2 DISCARD_READ=true \
          TAG_NAME=XC NUM_BASES_BELOW_QUALITY=1
        """


rule filter_tagged_bam:
    input:
        "results/05_tagged/{sample}.xm_xc.bam",
    output:
        "results/05_tagged/{sample}.xm_xc.filtered.bam",
    shell:
        "FilterBam INPUT={input} OUTPUT={output} TAG_REJECT=XQ"


rule trim_tso:
    input:
        "results/05_tagged/{sample}.xm_xc.filtered.bam",
    output:
        bam="results/06_align_count/{sample}/{sample}.trim_tso.bam",
        summary="results/06_align_count/{sample}/{sample}.tso.summary.txt",
    shell:
        """
        mkdir -p results/06_align_count/{wildcards.sample}
        TrimStartingSequence \
          INPUT={input} OUTPUT={output.bam} OUTPUT_SUMMARY={output.summary} \
          SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG MISMATCHES=0 NUM_BASES=5
        """


rule trim_polya:
    input:
        "results/06_align_count/{sample}/{sample}.trim_tso.bam",
    output:
        bam="results/06_align_count/{sample}/{sample}.trim_polya.bam",
        summary="results/06_align_count/{sample}/{sample}.polya.summary.txt",
    shell:
        """
        PolyATrimmer \
          INPUT={input} OUTPUT={output.bam} OUTPUT_SUMMARY={output.summary} \
          MISMATCHES=0 NUM_BASES=6
        """


rule bam_to_fastq:
    input:
        "results/06_align_count/{sample}/{sample}.trim_polya.bam",
    output:
        "results/06_align_count/{sample}/{sample}.cdna.fastq.gz",
    shell:
        """
        fq=results/06_align_count/{wildcards.sample}/{wildcards.sample}.cdna.fastq
        picard SamToFastq I={input} F=$fq TMP_DIR=$(mktemp -d) VALIDATION_STRINGENCY=SILENT
        gzip -f $fq
        """


rule star_align:
    input:
        fastq="results/06_align_count/{sample}/{sample}.cdna.fastq.gz",
        index="star_index",
    output:
        sam="results/06_align_count/{sample}/star_Aligned.out.sam",
        log_final="results/06_align_count/{sample}/star_Log.final.out",
    threads: 16
    shell:
        """
        STAR --runThreadN {threads} \
          --genomeDir {input.index} \
          --readFilesIn {input.fastq} \
          --readFilesCommand zcat \
          --outFileNamePrefix results/06_align_count/{wildcards.sample}/star_ \
          --outSAMtype SAM \
          --outSAMunmapped Within \
          --outSAMattributes NH HI AS nM
        """


rule sort_aligned_bam:
    input:
        "results/06_align_count/{sample}/star_Aligned.out.sam",
    output:
        "results/06_align_count/{sample}/{sample}.aligned.queryname.bam",
    shell:
        "picard SortSam I={input} O={output} SORT_ORDER=queryname TMP_DIR=$(mktemp -d)"


rule merge_bam_alignment:
    input:
        unmapped="results/06_align_count/{sample}/{sample}.trim_polya.bam",
        aligned="results/06_align_count/{sample}/{sample}.aligned.queryname.bam",
        genome=GENOME,
    output:
        "results/06_align_count/{sample}/{sample}.merged.bam",
    shell:
        """
        picard MergeBamAlignment \
          REFERENCE_SEQUENCE={input.genome} \
          UNMAPPED_BAM={input.unmapped} \
          ALIGNED_BAM={input.aligned} \
          OUTPUT={output} \
          INCLUDE_SECONDARY_ALIGNMENTS=false \
          PAIRED_RUN=false \
          TMP_DIR=$(mktemp -d)
        """


rule tag_gene_function:
    input:
        bam="results/06_align_count/{sample}/{sample}.merged.bam",
        refflat="dropseq_ref/genes.refFlat",
    output:
        "results/06_align_count/{sample}/{sample}.gene_function_tagged.bam",
    shell:
        "TagReadWithGeneFunction INPUT={input.bam} OUTPUT={output} ANNOTATIONS_FILE={input.refflat}"


rule digital_expression:
    input:
        "results/06_align_count/{sample}/{sample}.gene_function_tagged.bam",
    output:
        dge="results/06_align_count/{sample}/{sample}.DGE_with_introns.txt.gz",
        summary="results/06_align_count/{sample}/{sample}.DGE.summary.txt",
    shell:
        """
        DigitalExpression \
          INPUT={input} OUTPUT={output.dge} SUMMARY={output.summary} \
          READ_MQ=10 EDIT_DISTANCE=1 MIN_NUM_GENES_PER_CELL=100 \
          LOCUS_FUNCTION_LIST=CODING \
          LOCUS_FUNCTION_LIST=UTR \
          LOCUS_FUNCTION_LIST=INTRONIC
        """


rule build_seurat_objects:
    input:
        expand("results/06_align_count/{sample}/{sample}.DGE_with_introns.txt.gz", sample=SAMPLES),
    output:
        expand("results/07_seurat_raw/{sample}.seurat.raw.rds", sample=SAMPLES),
    shell:
        "Rscript 06_build_seurat_objects.R results/06_align_count results/07_seurat_raw"


rule qc_doubletfinder_sct:
    input:
        "results/07_seurat_raw/{sample}.seurat.raw.rds",
    output:
        "results/08_seurat_qc/{sample}.qc_sct.rds",
    shell:
        """
        mkdir -p results/08_seurat_qc
        Rscript 08_qc_doubletfinder_sct.R {input} {output}
        """


rule integrate_cluster:
    input:
        expand("results/08_seurat_qc/{sample}.qc_sct.rds", sample=SAMPLES),
    output:
        "results/09_integrated/integrated_seurat.rds",
    shell:
        """
        mkdir -p results/09_integrated
        Rscript 09_integrate_cluster.R results/08_seurat_qc {output}
        """


rule find_markers:
    input:
        "results/09_integrated/integrated_seurat.rds",
    output:
        "results/10_markers/all_markers.tsv",
    shell:
        """
        mkdir -p results/10_markers
        Rscript 10_find_marker.R {input} {output}
        """
