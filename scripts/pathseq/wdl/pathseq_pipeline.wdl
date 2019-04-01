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

import "pathseq_3_stage_pipeline.wdl" as pathseq_pipeline
import "pathseq_with_downsampling.wdl" as pathseq_downsample

workflow PathSeqPipeline {

  String sample_name
  File input_bam
  File input_bam_index = input_bam + ".bai"

  # Number of reads to downsample to for estimating filter metrics
  Int filter_metrics_reads = 1000000

  File kmer_file
  File filter_bwa_image
  File microbe_bwa_image
  File microbe_dict
  File taxonomy_file

  Boolean is_host_aligned
  Boolean? filter_duplicates
  Boolean? skip_pre_bwa_repartition
  Int? filter_bam_partition_size
  Int? filter_bwa_seed_length
  Int? host_min_identity
  Int? min_clipped_read_length

  Float? min_score_identity
  Float? identity_margin
  Boolean? divide_by_genome_length

  File? gatk4_jar_override

  # Runtime parameters
  String gatk_docker

  Int? downsample_preemptible_attempts
  Int? downsample_filter_preemptible_attempts
  Int? filter_preemptible_attempts
  Int? align_preemptible_attempts
  Int? score_preemptible_attempts

  Int? downsample_filter_cpu
  Int? filter_cpu
  Int? align_cpu
  Int? score_cpu

  Int? downsample_filter_mem_gb
  Int? filter_mem_gb
  Int? align_mem_gb
  Int? score_mem_gb

  Boolean? downsample_filter_ssd
  Boolean? filter_ssd
  Boolean? align_ssd
  Boolean? score_ssd

  # Optional input to increase all disk sizes in case of outlier sample with strange size behavior
  Int? downsample_additional_disk_gb
  Int? downsample_filter_additional_disk_gb
  Int? filter_additional_disk_gb
  Int? align_additional_disk_gb
  Int? score_additional_disk_gb

  call pathseq_downsample.PathSeqFilterWithDownsampling {
    input:
      sample_name=sample_name,
      input_bam=input_bam,
      reads_after_downsampling=filter_metrics_reads,
      kmer_file=kmer_file,
      filter_bwa_image=filter_bwa_image,
      is_host_aligned=is_host_aligned,
      gather_filter_metrics=true,
      filter_duplicates=filter_duplicates,
      min_clipped_read_length=min_clipped_read_length,
      filter_bam_partition_size=filter_bam_partition_size,
      host_min_identity=host_min_identity,
      filter_bwa_seed_length=filter_bwa_seed_length,
      gatk4_jar_override=gatk4_jar_override,
      filter_mem_gb=filter_mem_gb,
      gatk_docker=gatk_docker,
      filter_preemptible_attempts=downsample_filter_preemptible_attempts,
      filter_cpu=downsample_filter_cpu,
      filter_mem_gb=downsample_filter_mem_gb,
      filter_ssd=downsample_filter_ssd,
      filter_additional_disk_gb=downsample_filter_additional_disk_gb
  }

  call pathseq_pipeline.PathSeqThreeStageWorkflow {
    input:
      sample_name=sample_name,
      input_bam_or_cram=input_bam,
      kmer_file=kmer_file,
      filter_bwa_image=filter_bwa_image,
      microbe_bwa_image=microbe_bwa_image,
      microbe_dict=microbe_dict,
      taxonomy_file=taxonomy_file,
      gather_filter_metrics=false,
      is_host_aligned=is_host_aligned,
      filter_duplicates=filter_duplicates,
      min_clipped_read_length=min_clipped_read_length,
      filter_bam_partition_size=filter_bam_partition_size,
      host_min_identity=host_min_identity,
      filter_bwa_seed_length=filter_bwa_seed_length,
      min_score_identity=min_score_identity,
      identity_margin=identity_margin,
      divide_by_genome_length=divide_by_genome_length,
      gatk4_jar_override=gatk4_jar_override,
      filter_mem_gb=filter_mem_gb,
      gatk_docker=gatk_docker,
      filter_preemptible_attempts=filter_preemptible_attempts,
      align_preemptible_attempts=align_preemptible_attempts,
      score_preemptible_attempts=score_preemptible_attempts,
      filter_cpu=filter_cpu,
      align_cpu=align_cpu,
      score_cpu=score_cpu,
      filter_mem_gb=filter_mem_gb,
      align_mem_gb=align_mem_gb,
      score_mem_gb=score_mem_gb,
      filter_ssd=filter_ssd,
      align_ssd=align_ssd,
      score_ssd=score_ssd,
      filter_additional_disk_gb=filter_additional_disk_gb,
      align_additional_disk_gb=align_additional_disk_gb,
      score_additional_disk_gb=score_additional_disk_gb
  }

  output {
    File final_bam = PathSeqThreeStageWorkflow.final_bam
    File scores = PathSeqThreeStageWorkflow.scores
    File score_metrics = PathSeqThreeStageWorkflow.score_metrics

    File downsampled_filter_metrics = PathSeqFilterWithDownsampling.filter_metrics
    Int original_total_reads = PathSeqFilterWithDownsampling.original_total_reads
  }
}