########################################################################################################################
## PathSeq Estimate Filter Metrics WDL
########################################################################################################################
##
## Estimates filter metrics for PathSeq
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
## - Picard-style metrics files for the filter phase of the pipeline at the downsampled coverage
##
########################################################################################################################

import "pathseq_3_stage_pipeline.wdl" as pathseq_pipeline

# Downsamples BAM to a specified number of reads
task Downsample {
  File input_bam_file
  File input_bam_index_file
  String downsampled_bam_filename

  Int reads_after_downsampling

  Int additional_disk_gb = 20
  Int preemptible_tries

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
      P=$P_DOWNSAMPLE
  >>>
  output {
    File output_bam_file = "${downsampled_bam_filename}"
    Int total_reads = read_tsv("num_reads.txt")[0][0]
  }
  runtime {
    preemptible: "${preemptible_tries}"
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.2.5-1486412288"
    memory: "3.75 GiB"
    cpu: "1"
    disks: "local-disk ${disk_size} HDD"
  }
}

workflow PathSeqFilterWithDownsampling {

  String sample_name
  File input_bam
  File input_bam_index
  Int reads_after_downsampling

  File kmer_file
  File filter_bwa_image

  Boolean? gather_filter_metrics
  Boolean is_host_aligned
  Boolean? filter_duplicates
  Boolean? skip_pre_bwa_repartition
  Int? filter_bam_partition_size
  Int? filter_bwa_seed_length
  Int? host_min_identity
  Int? min_clipped_read_length

  File? gatk4_jar_override

  # Runtime parameters
  String gatk_docker

  Int? downsample_preemptible_attempts
  Int? filter_preemptible_attempts
  Int? filter_cpu
  Int? filter_mem_gb
  Boolean? filter_ssd

  # Optional input to increase all disk sizes in case of outlier sample with strange size behavior
  Int? downsample_additional_disk_gb
  Int? filter_additional_disk_gb

  # Downsample input bam
  call Downsample {
    input:
      input_bam_file=input_bam,
      input_bam_index_file=input_bam_index,
      downsampled_bam_filename="${sample_name}.downsampled.bam",
      reads_after_downsampling=reads_after_downsampling,
      additional_disk_gb=downsample_additional_disk_gb,
      preemptible_tries=downsample_preemptible_attempts
  }

  call pathseq_pipeline.PathSeqFilter {
    input:
      sample_name=sample_name,
      input_bam_or_cram=Downsample.output_bam_file,
      kmer_file=kmer_file,
      filter_bwa_image=filter_bwa_image,
      gather_metrics=gather_filter_metrics,
      is_host_aligned=is_host_aligned,
      filter_duplicates=filter_duplicates,
      min_clipped_read_length=min_clipped_read_length,
      bam_partition_size=filter_bam_partition_size,
      host_min_identity=host_min_identity,
      filter_bwa_seed_length=filter_bwa_seed_length,
      gatk4_jar_override=gatk4_jar_override,
      gatk_docker=gatk_docker,
      preemptible_attempts=filter_preemptible_attempts,
      mem_gb=filter_mem_gb,
      cpu=filter_cpu,
      use_ssd=filter_ssd,
      additional_disk_gb=filter_additional_disk_gb
  }

  output {
    File filter_metrics = PathSeqFilter.filter_metrics
    Int original_total_reads = Downsample.total_reads
  }
}