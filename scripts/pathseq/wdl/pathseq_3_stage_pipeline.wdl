########################################################################################################################
## PathSeq Pipeline WDL
########################################################################################################################
##
## Runs the PathSeq pipeline
##
## For further info see the GATK Documentation for the PathSeqPipelineSpark tool:
##   https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_spark_pathseq_PathSeqPipelineSpark.php
##
########################################################################################################################
##
## Input requirements :
## - Sequencing data in BAM format
## - Host and microbe references files available in the GATK Resource Bundle (available on FTP):
##     https://software.broadinstitute.org/gatk/download/bundle
##
## - Input BAM files must comply with the following requirements:
## - - file must pass validation by ValidateSamFile
## - - all reads must have an RG tag
## - - one or more read groups all belong to a single sample (SM)
##
## Output:
## - BAM file containing microbe-mapped reads and reads of unknown sequence
## - Tab-separated value (.tsv) file of taxonomic abundance scores
## - Picard-style metrics files for the filter and scoring phases of the pipeline
##
########################################################################################################################


workflow PathSeqPipeline {

  String sample_name
  File input_bam_or_cram
  File? input_bam_or_cram_index

  # If enabled, filter metrics will be estimated using a downsampled bam with this many reads (recommended)
  # If disabled, no filter metrics will be generated
  Boolean estimate_filter_metrics_with_downsampling = true
  Int estimate_filter_metrics_reads = 1000000

  # Required if the input is a cram
  File? cram_reference_fasta
  File? cram_reference_fasta_index
  File? cram_reference_dict

  File kmer_file
  File filter_bwa_image
  File microbe_bwa_image
  File microbe_dict
  File taxonomy_file

  Boolean is_host_aligned
  Boolean? filter_duplicates
  Boolean? skip_pre_bwa_repartition
  Int? filter_bam_partition_size
  Boolean? divide_by_genome_length
  Int? filter_bwa_seed_length
  Int? host_min_identity
  Int? min_clipped_read_length
  Float? min_score_identity
  Float? identity_margin

  File? gatk4_jar_override

  # Runtime parameters
  String gatk_docker

  Int? cram_to_bam_preemptible_attempts
  Int? downsample_preemptible_attempts
  Int? filter_preemptible_attempts
  Int? align_preemptible_attempts
  Int? score_preemptible_attempts

  Int? cram_to_bam_cpu
  Int? filter_cpu
  Int? align_cpu
  Int? score_cpu

  Float? cram_to_bam_mem_gb
  Int? filter_mem_gb
  Int? align_mem_gb
  Int? score_mem_gb

  Boolean? filter_ssd
  Boolean? align_ssd
  Boolean? score_ssd

  Int? cram_to_bam_max_retries

  # Optional input to increase all disk sizes in case of outlier sample with strange size behavior
  Int? cram_to_bam_additional_disk_gb
  Int? downsample_additional_disk_gb
  Int? filter_additional_disk_gb
  Int? align_additional_disk_gb
  Int? score_additional_disk_gb

  Boolean is_bam = basename(input_bam_or_cram, ".bam") + ".bam" == basename(input_bam_or_cram)
  String index_ext = if is_bam then ".bai" else ".crai"
  File bam_or_cram_index = if defined(input_bam_or_cram_index) then select_first([input_bam_or_cram_index]) else input_bam_or_cram + index_ext

  # Convert to BAM if we have a CRAM
  if (!is_bam) {
    call CramToBam {
      input:
        cram_file = input_bam_or_cram,
        cram_index = bam_or_cram_index,
        reference_fasta = select_first([cram_reference_fasta]),
        reference_index = select_first([cram_reference_fasta_index]),
        docker = gatk_docker,
        cpu = cram_to_bam_cpu,
        mem_gb = cram_to_bam_mem_gb,
        extra_disk_gb = cram_to_bam_additional_disk_gb,
        preemptible = cram_to_bam_preemptible_attempts,
        max_retries = cram_to_bam_max_retries,
    }
  }

  File bam_file = select_first([CramToBam.bam_file, input_bam_or_cram])
  File bam_index = select_first([CramToBam.bam_index, bam_or_cram_index])

  if (estimate_filter_metrics_with_downsampling) {
    # Downsample bam for filter metrics estimation
    call Downsample as DownsampleFilterMetricsBam {
      input:
        input_bam_file=bam_file,
        input_bam_index_file=bam_index,
        downsampled_bam_filename="${sample_name}.downsampled.bam",
        reads_after_downsampling=estimate_filter_metrics_reads,
        additional_disk_gb=downsample_additional_disk_gb,
        preemptible_tries=downsample_preemptible_attempts
    }

    call PathSeqFilter as GetDownsampledFilterMetrics {
      input:
        sample_name=sample_name,
        input_bam_or_cram=DownsampleFilterMetricsBam.output_bam_file,
        kmer_file=kmer_file,
        filter_bwa_image=filter_bwa_image,
        gather_metrics=true,
        is_host_aligned=is_host_aligned,
        filter_duplicates=filter_duplicates,
        min_clipped_read_length=min_clipped_read_length,
        bam_partition_size=filter_bam_partition_size,
        host_min_identity=host_min_identity,
        filter_bwa_seed_length=filter_bwa_seed_length,
        skip_pre_bwa_repartition=skip_pre_bwa_repartition,
        gatk4_jar_override=gatk4_jar_override,
        mem_gb=filter_mem_gb,
        gatk_docker=gatk_docker,
        preemptible_attempts=filter_preemptible_attempts,
        additional_disk_gb=filter_additional_disk_gb,
        cpu=filter_cpu,
        use_ssd=filter_ssd
    }
  }

  call PathSeqFilter as PathSeqFilter {
    input:
      sample_name=sample_name,
      input_bam_or_cram=bam_file,
      kmer_file=kmer_file,
      filter_bwa_image=filter_bwa_image,
      gather_metrics=false,
      is_host_aligned=is_host_aligned,
      filter_duplicates=filter_duplicates,
      min_clipped_read_length=min_clipped_read_length,
      bam_partition_size=filter_bam_partition_size,
      host_min_identity=host_min_identity,
      filter_bwa_seed_length=filter_bwa_seed_length,
      skip_pre_bwa_repartition=skip_pre_bwa_repartition,
      gatk4_jar_override=gatk4_jar_override,
      mem_gb=filter_mem_gb,
      gatk_docker=gatk_docker,
      preemptible_attempts=filter_preemptible_attempts,
      additional_disk_gb=filter_additional_disk_gb,
      cpu=filter_cpu,
      use_ssd=filter_ssd
  }

  call PathSeqAlign {
    input:
      sample_name=sample_name,
      input_paired_bam=PathSeqFilter.paired_bam_out,
      input_unpaired_bam=PathSeqFilter.unpaired_bam_out,
      microbe_bwa_image=microbe_bwa_image,
      microbe_dict=microbe_dict,
      gatk4_jar_override=gatk4_jar_override,
      mem_gb=align_mem_gb,
      gatk_docker=gatk_docker,
      preemptible_attempts=align_preemptible_attempts,
      additional_disk_gb=align_additional_disk_gb,
      cpu=align_cpu,
      use_ssd=align_ssd
  }

  call PathSeqScore {
    input:
      sample_name=sample_name,
      input_paired_bam=PathSeqAlign.paired_bam_out,
      input_unpaired_bam=PathSeqAlign.unpaired_bam_out,
      taxonomy_file=taxonomy_file,
      divide_by_genome_length=divide_by_genome_length,
      min_score_identity=min_score_identity,
      identity_margin=identity_margin,
      gatk4_jar_override=gatk4_jar_override,
      mem_gb=score_mem_gb,
      gatk_docker=gatk_docker,
      preemptible_attempts=score_preemptible_attempts,
      additional_disk_gb=score_additional_disk_gb,
      cpu=score_cpu,
      use_ssd=score_ssd
  }

  output {
    File final_bam = PathSeqScore.bam_out
    File scores = PathSeqScore.scores
    File filter_metrics = select_first([GetDownsampledFilterMetrics.filter_metrics, PathSeqFilter.filter_metrics])
    File score_metrics = PathSeqScore.score_metrics
  }
}

task CramToBam {
  File cram_file
  File? cram_index
  File reference_fasta
  File? reference_index

  String docker

  Int? cpu = 4
  Float? mem_gb = 15
  Int? extra_disk_gb = 10
  Int? preemptible = 3
  Int? max_retries = 1

  String bam_file_name = basename(cram_file, ".cram") + ".bam"

  File cram_index_file = select_first([cram_index, cram_file + ".crai"])
  File reference_index_file = select_first([reference_index, reference_fasta + ".fai"])

  Float cram_inflate_ratio = 3.0
  Float cram_size = size(cram_file, "GiB")
  Float cram_index_size = size(cram_index_file, "GiB")
  Float bam_size = cram_inflate_ratio * cram_size
  Float bam_index_size = cram_index_size
  Float ref_size = size(reference_fasta, "GiB")
  Float ref_index_size = size(reference_index_file, "GiB")
  Int vm_disk_size = ceil(cram_size + cram_index_size + bam_size + bam_index_size + ref_size + ref_index_size + extra_disk_gb)

  output {
    File bam_file = bam_file_name
    File bam_index = bam_file_name + ".bai"
  }
  command <<<

        set -Eeuo pipefail

        # covert cram to bam
        samtools view  -@ ${cpu} -h -b -T "${reference_fasta}" -o "${bam_file_name}" "${cram_file}"

        # index bam file
        samtools index -@ ${cpu} "${bam_file_name}"

  >>>
  runtime {
    cpu: 1
    memory: mem_gb + " GiB"
    disks: "local-disk " + vm_disk_size + " HDD"
    bootDiskSizeGb: 10
    docker: docker
    preemptible: preemptible
    maxRetries: max_retries
  }
}

# Downsamples BAM to a specified number of reads
task Downsample {
  File input_bam_file
  File input_bam_index_file
  String downsampled_bam_filename

  Int reads_after_downsampling

  Int additional_disk_gb = 20
  Int preemptible_tries = 3

  Int disk_size = ceil(size(input_bam_file, "GB")*2 + additional_disk_gb)

  command <<<
    set -e
    set -o pipefail
    NUM_READS=`samtools idxstats ${input_bam_file} | awk '{s+=$3+$4} END {print s}'`
    P_DOWNSAMPLE=`python -c "print ${reads_after_downsampling}/float($NUM_READS)"`
    RESULT=`python -c "print $P_DOWNSAMPLE > 1"`
    if [ "$RESULT" == "True" ]
    then
        P_DOWNSAMPLE="1"
    fi
    echo $NUM_READS > num_reads.txt
    java -Xmx2000m -jar /usr/gitc/picard.jar \
      DownsampleSam \
      INPUT=${input_bam_file} \
      OUTPUT=${downsampled_bam_filename} \
      P=$P_DOWNSAMPLE \
      VALIDATION_STRINGENCY=SILENT
  >>>
  output {
    File output_bam_file = "${downsampled_bam_filename}"
    Int total_reads = read_int("num_reads.txt")
  }
  runtime {
    preemptible: "${preemptible_tries}"
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.2.5-1486412288"
    memory: "3.75 GiB"
    cpu: "1"
    disks: "local-disk ${disk_size} HDD"
  }
}

task PathSeqFilter {

  # Inputs for this task
  String sample_name
  File input_bam_or_cram

  File kmer_file
  File filter_bwa_image

  # Required if cram is provided
  File? cram_reference_fasta
  File? cram_reference_fasta_index
  File? cram_reference_dict

  Boolean is_host_aligned
  Boolean skip_quality_filters = false
  Boolean skip_pre_bwa_repartition = false
  Boolean filter_duplicates = true
  Boolean gather_metrics = true

  Int filter_bwa_seed_length = 19
  Int host_min_identity = 30
  Int min_clipped_read_length = 60
  Int bam_partition_size = 4000000

  String paired_bam_output_path = "${sample_name}.non_host.paired.bam"
  String unpaired_bam_output_path = "${sample_name}.non_host.unpaired.bam"
  String filter_metrics_output_path = "${sample_name}.pathseq.filter_metrics"

  File? gatk4_jar_override

  # Default to WARNING which will avoid excessive Spark logging
  String verbosity = "WARNING"

  # Runtime parameters
  String gatk_docker
  Int mem_gb = 32
  Int preemptible_attempts = 3
  Float additional_disk_gb = 50
  Int cpu = 8
  Boolean use_ssd = false

  Int disk_size = ceil(size(input_bam_or_cram, "GB") + size(kmer_file, "GB") + size(filter_bwa_image, "GB") + additional_disk_gb)

  # Mem is in units of GB but our command and memory runtime values are in MB
  Int machine_mem = mem_gb * 1000
  Int command_mem = ceil((machine_mem - size(filter_bwa_image, "MB")) * 0.8)

  command <<<
    set -euo pipefail
    export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk4_jar_override}
    touch ${filter_metrics_output_path}
    gatk --java-options "-Xmx${command_mem}m" \
      PathSeqFilterSpark \
      --input ${input_bam_or_cram} \
      --paired-output ${paired_bam_output_path} \
      --unpaired-output ${unpaired_bam_output_path} \
      ${if defined(cram_reference_fasta) then "--reference ${cram_reference_fasta}" else ""} \
      ${if gather_metrics then "--filter-metrics ${filter_metrics_output_path}" else ""} \
      --kmer-file ${kmer_file} \
      --filter-bwa-image ${filter_bwa_image} \
      --bam-partition-size ${bam_partition_size} \
      --is-host-aligned ${is_host_aligned} \
      --skip-quality-filters ${skip_quality_filters} \
      --min-clipped-read-length ${min_clipped_read_length} \
      --filter-bwa-seed-length ${filter_bwa_seed_length} \
      --host-min-identity ${host_min_identity} \
      --filter-duplicates ${filter_duplicates} \
      --skip-pre-bwa-repartition ${skip_pre_bwa_repartition} \
      --verbosity ${verbosity}

    if [ ! -f "${paired_bam_output_path}" ]; then
    	echo "File ${paired_bam_output_path} not found, creating empty BAM"
    	printf "@HD     VN:1.5  SO:queryname\n@RG     ID:A    SM:${sample_name}\n" | samtools view -hb - > ${paired_bam_output_path}
    fi

    if [ ! -f "${unpaired_bam_output_path}" ]; then
    	echo "File ${unpaired_bam_output_path} not found, creating empty BAM"
    	printf "@HD     VN:1.5  SO:queryname\n@RG     ID:A    SM:${sample_name}\n" | samtools view -hb - > ${unpaired_bam_output_path}
    fi
  >>>
  runtime {
    docker: gatk_docker
    memory: machine_mem + " MB"
    # Note that the space before SSD and HDD should be included.
    disks: "local-disk " + disk_size + if use_ssd then " SSD" else " HDD"
    preemptible: preemptible_attempts
    cpu: cpu
  }
  output {
    File paired_bam_out = "${paired_bam_output_path}"
    File unpaired_bam_out = "${unpaired_bam_output_path}"
    File filter_metrics = "${sample_name}.pathseq.filter_metrics"
  }
}

task PathSeqAlign {

  # Inputs for this task
  String sample_name
  File input_paired_bam
  File input_unpaired_bam

  File microbe_bwa_image
  File microbe_dict

  String paired_bam_output_path = "${sample_name}.microbe_aligned.paired.bam"
  String unpaired_bam_output_path = "${sample_name}.microbe_aligned.unpaired.bam"

  File? gatk4_jar_override

  # Default to WARNING which will avoid excessive Spark logging
  String verbosity = "WARNING"

  # Runtime parameters
  String gatk_docker
  Int mem_gb = 140
  Int preemptible_attempts = 3
  Float additional_disk_gb = 20
  Int cpu = 32
  Boolean use_ssd = false

  Int disk_size = ceil(2*size(input_paired_bam, "GB") + 2*size(input_unpaired_bam, "GB") + size(microbe_bwa_image, "GB") + additional_disk_gb)

  # Mem is in units of GB but our command and memory runtime values are in MB
  Int machine_mem = mem_gb * 1000
  Int command_mem = ceil((machine_mem - size(microbe_bwa_image, "MB")) * 0.8)

  command <<<
    set -e
    export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk4_jar_override}
    gatk --java-options "-Xmx${command_mem}m" \
      PathSeqBwaSpark \
      --paired-input ${input_paired_bam} \
      --unpaired-input ${input_unpaired_bam} \
      --paired-output ${paired_bam_output_path} \
      --unpaired-output ${unpaired_bam_output_path} \
      --microbe-bwa-image ${microbe_bwa_image} \
      --microbe-dict ${microbe_dict} \
      --verbosity ${verbosity}
  >>>
  runtime {
    docker: gatk_docker
    memory: machine_mem + " MB"
    # Note that the space before SSD and HDD should be included.
    disks: "local-disk " + disk_size + if use_ssd then " SSD" else " HDD"
    preemptible: preemptible_attempts
    cpu: cpu
  }
  output {
    File paired_bam_out = "${paired_bam_output_path}"
    File unpaired_bam_out = "${unpaired_bam_output_path}"
  }
}

task PathSeqScore {

  # Inputs for this task
  String sample_name
  File input_paired_bam
  File input_unpaired_bam

  File taxonomy_file

  Boolean divide_by_genome_length = true
  Float min_score_identity = 0.9
  Float identity_margin = 0.02

  String bam_output_path = "${sample_name}.pathseq.bam"
  String scores_output_path = "${sample_name}.pathseq.tsv"
  String score_metrics_output_path = "${sample_name}.pathseq.score_metrics"

  File? gatk4_jar_override

  # Runtime parameters
  String gatk_docker
  Int mem_gb = 8
  Int preemptible_attempts = 3
  Float additional_disk_gb = 10
  Int cpu = 2
  Boolean use_ssd = false

  Int disk_size = ceil(size(input_paired_bam, "GB") + size(input_unpaired_bam, "GB") + additional_disk_gb)

  # Mem is in units of GB but our command and memory runtime values are in MB
  Int machine_mem = mem_gb * 1000
  Int command_mem = ceil(machine_mem * 0.8)

  command <<<
    set -e
    export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk4_jar_override}
    gatk --java-options "-Xmx${command_mem}m" \
      PathSeqScoreSpark \
      --paired-input ${input_paired_bam} \
      --unpaired-input ${input_unpaired_bam} \
      --output ${bam_output_path} \
      --scores-output ${scores_output_path} \
      --score-metrics ${score_metrics_output_path} \
      --taxonomy-file ${taxonomy_file} \
      --min-score-identity ${min_score_identity} \
      --identity-margin ${identity_margin} \
      --divide-by-genome-length ${divide_by_genome_length}
  >>>
  runtime {
    docker: gatk_docker
    memory: machine_mem + " MB"
    # Note that the space before SSD and HDD should be included.
    disks: "local-disk " + disk_size + if use_ssd then " SSD" else " HDD"
    preemptible: preemptible_attempts
    cpu: cpu
  }
  output {
    File bam_out = "${sample_name}.pathseq.bam"
    File scores = "${sample_name}.pathseq.tsv"
    File score_metrics = "${sample_name}.pathseq.score_metrics"
  }
}
