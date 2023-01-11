version 1.0

import "../common/structs.wdl"
import "../common/tasks/samtools_fasta.wdl" as SamtoolsFasta
import "../assemble_genome/assemble_genome.wdl" as AssembleGenome
import "../common/tasks/zip_index_vcf.wdl" as ZipIndexVcf
import "../common/tasks/bcftools_stats.wdl" as BcftoolsStats

workflow de_novo_assembly_sample {
	input {
		Sample sample

		ReferenceData reference

		RuntimeAttributes spot_runtime_attributes
	}

	scatter (movie_bam in sample.movie_bams) {
		call SamtoolsFasta.samtools_fasta {
			input:
				bam = movie_bam,
				runtime_attributes = spot_runtime_attributes
		}
	}

	call AssembleGenome.assemble_genome {
		input:
			sample_id = sample.sample_id,
			reads_fastas = samtools_fasta.reads_fasta,
			reference = reference,
			spot_runtime_attributes = spot_runtime_attributes
	}

	call htsbox {
		input:
			bam = assemble_genome.asm_bam.data,
			bam_index = assemble_genome.asm_bam.data_index,
			reference = reference.fasta.data,
			runtime_attributes = spot_runtime_attributes
	}

	call ZipIndexVcf.zip_index_vcf {
		input:
			vcf = htsbox.htsbox_vcf,
			runtime_attributes = spot_runtime_attributes
	}

	call BcftoolsStats.bcftools_stats {
		input:
			vcf = zip_index_vcf.zipped_vcf,
			bam = assemble_genome.asm_bam.data,
			reference = reference.fasta.data,
			runtime_attributes = spot_runtime_attributes
	}

	output {
		Array[File] assembly_noseq_gfas = assemble_genome.assembly_noseq_gfas
		Array[File] assembly_lowQ_beds = assemble_genome.assembly_lowQ_beds
		Array[File] zipped_assembly_fastas = assemble_genome.zipped_assembly_fastas
		Array[File] assembly_stats = assemble_genome.assembly_stats
		IndexData asm_bam = assemble_genome.asm_bam
		IndexData htsbox_vcf = {"data": zip_index_vcf.zipped_vcf, "data_index": zip_index_vcf.zipped_vcf_index}
		File htsbox_vcf_stats = bcftools_stats.stats
	}

	parameter_meta {
		sample: {help: "Sample information and associated data files"}
		reference: {help: "Reference genome data"}
		spot_runtime_attributes: {help: "RuntimeAttributes for spot (preemptible) tasks"}
	}
}

task htsbox {
	input {
		File bam
		File bam_index

		File reference

		RuntimeAttributes runtime_attributes
	}

	String bam_basename = basename(bam, ".bam")
	Int threads = 4
	Int disk_size = ceil((size(bam, "GB") + size(reference, "GB")) * 3 + 200)

	command <<<
		set -euo pipefail

		htsbox pileup \
			-q20 \
			-c \
			-f \
			~{reference} \
			~{bam} \
		> ~{bam_basename}.htsbox.vcf
	>>>

	output {
		File htsbox_vcf = "~{bam_basename}.htsbox.vcf"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/htsbox:b1a46c6"
		cpu: threads
		memory: "14 GB"
		disk: disk_size + " GB"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: runtime_attributes.preemptible_tries
		maxRetries: runtime_attributes.max_retries
		awsBatchRetryAttempts: runtime_attributes.max_retries
		queueArn: runtime_attributes.queue_arn
		zones: runtime_attributes.zones
	}
}
