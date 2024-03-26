version 1.0

# Run for each sample in the cohort. Aligns reads from each movie to the reference genome, then calls and phases small and structural variants.

import "../humanwgs_structs.wdl"
import "../wdl-common/wdl/tasks/pbsv_discover.wdl" as PbsvDiscover
import "../wdl-common/wdl/tasks/pbsv_call.wdl" as PbsvCall
import "../wdl-common/wdl/tasks/samtools_merge.wdl" as SamtoolsMerge
import "../wdl-common/wdl/workflows/deepvariant/deepvariant.wdl" as DeepVariant
import "../wdl-common/wdl/tasks/trgt.wdl" as Trgt
import "../wdl-common/wdl/tasks/mosdepth.wdl" as Mosdepth
import "../wdl-common/wdl/tasks/concat_vcf.wdl" as ConcatVcf
import "../wdl-common/wdl/tasks/paraphase.wdl" as Paraphase
import "../wdl-common/wdl/tasks/hificnv.wdl" as Hificnv
import "../wdl-common/wdl/workflows/hiphase/hiphase.wdl" as HiPhase
import "../wdl-common/wdl/tasks/cpg_pileup.wdl" as CpgPileup

workflow sample_analysis {
	input {
		Sample sample

		ReferenceData reference

		String deepvariant_version
		File? custom_deepvariant_model_tar

		RuntimeAttributes default_runtime_attributes
	}

	Array[Array[String]] pbsv_splits = read_json(reference.pbsv_splits)

	scatter (movie_bam in sample.movie_bams) {
		call pbmm2_align {
			input:
				sample_id = sample.sample_id,
				bam = movie_bam,
				reference = reference.fasta.data,
				reference_index = reference.fasta.data_index,
				reference_name = reference.name,
				runtime_attributes = default_runtime_attributes
		}

		call PbsvDiscover.pbsv_discover {
			input:
				aligned_bam = pbmm2_align.aligned_bam,
				aligned_bam_index = pbmm2_align.aligned_bam_index,
				reference_tandem_repeat_bed = reference.tandem_repeat_bed,
				runtime_attributes = default_runtime_attributes
		}

		IndexData aligned_bam = {
			"data": pbmm2_align.aligned_bam,
			"data_index": pbmm2_align.aligned_bam_index
		}
	}

	scatter (shard_index in range(length(pbsv_splits))) {
		Array[String] region_set = pbsv_splits[shard_index]

		call PbsvCall.pbsv_call {
			input:
				sample_id = sample.sample_id,
				svsigs = pbsv_discover.svsig,
				reference = reference.fasta.data,
				reference_index = reference.fasta.data_index,
				reference_name = reference.name,
				shard_index = shard_index,
				regions = region_set,
				runtime_attributes = default_runtime_attributes
		}
	}

	# concatenate pbsv vcfs
	call ConcatVcf.concat_vcf {
		input:
			vcfs = pbsv_call.pbsv_vcf,
			vcf_indices = pbsv_call.pbsv_vcf_index,
			output_vcf_name = "~{sample.sample_id}.~{reference.name}.pbsv.vcf.gz",
			runtime_attributes = default_runtime_attributes
	}

	# merge aligned bams if there are multiple
	if (length(aligned_bam) > 1) {
		scatter (bam_object in aligned_bam) {
			File bam_to_merge = bam_object.data
		}
		call SamtoolsMerge.merge_bams {
			input:
				bams = bam_to_merge,
				output_bam_name = "~{sample.sample_id}.~{reference.name}.bam",
				runtime_attributes = default_runtime_attributes
		}
	}

	# select the merged bam if it exists, otherwise select the first (only) aligned bam
	File aligned_bam_data = select_first([merge_bams.merged_bam, aligned_bam[0].data])
	File aligned_bam_index = select_first([merge_bams.merged_bam_index, aligned_bam[0].data_index])

	call Mosdepth.mosdepth {
		input:
			aligned_bam = aligned_bam_data,
			aligned_bam_index = aligned_bam_index,
			runtime_attributes = default_runtime_attributes
	}

	call DeepVariant.deepvariant {
		input:
			sample_id = sample.sample_id,
			aligned_bams = aligned_bam,
			reference_fasta = reference.fasta,
			reference_name = reference.name,
			deepvariant_version = deepvariant_version,
			custom_deepvariant_model_tar = custom_deepvariant_model_tar,
			default_runtime_attributes = default_runtime_attributes
	}

	call bcftools {
		input:
			vcf = deepvariant.vcf.data,
			stats_params = "--apply-filters PASS --samples ~{sample.sample_id}",
			reference = reference.fasta.data,
			runtime_attributes = default_runtime_attributes
	}

	call Trgt.trgt {
		input:
			sample_id = sample.sample_id,
			sex = sample.sex,
			bam = aligned_bam_data,
			bam_index = aligned_bam_index,
			reference = reference.fasta.data,
			reference_index = reference.fasta.data_index,
			tandem_repeat_bed = reference.trgt_tandem_repeat_bed,
			output_prefix = "~{sample.sample_id}.~{reference.name}",
			runtime_attributes = default_runtime_attributes
	}

	call HiPhase.hiphase {
		# vcfs order: small variants, SVs, TRGT
		input:
			id = sample.sample_id,
			refname = reference.name,
			sample_ids = [sample.sample_id],
			vcfs = [
				deepvariant.vcf,
				{"data": concat_vcf.concatenated_vcf, "data_index": concat_vcf.concatenated_vcf_index}, 
				{"data": trgt.repeat_vcf, "data_index": trgt.repeat_vcf_index}
				],
			bams = [{"data": aligned_bam_data, "data_index": aligned_bam_index}],
			haplotag = true,
			reference_fasta = reference.fasta,
			default_runtime_attributes = default_runtime_attributes
	}

	IndexData haplotagged_bam = {
		"data": hiphase.haplotagged_bams[0].data,
		"data_index": hiphase.haplotagged_bams[0].data_index
	}

	call CpgPileup.cpg_pileup {
		input:
			bam = haplotagged_bam.data,
			bam_index = haplotagged_bam.data_index,
			output_prefix = "~{sample.sample_id}.~{reference.name}",
			reference = reference.fasta.data,
			reference_index = reference.fasta.data_index,
			runtime_attributes = default_runtime_attributes
	}

	call Paraphase.paraphase {
		input:
			sample_id = sample.sample_id,
			bam = haplotagged_bam.data,
			bam_index = haplotagged_bam.data_index,
			reference = reference.fasta.data,
			reference_index = reference.fasta.data_index,
			out_directory = "~{sample.sample_id}.paraphase",
			runtime_attributes = default_runtime_attributes
	}

	call Hificnv.hificnv {
		input:
			sample_id = sample.sample_id,
			sex = sample.sex,
			bam = haplotagged_bam.data,
			bam_index = haplotagged_bam.data_index,
			phased_vcf = hiphase.phased_vcfs[0].data,
			phased_vcf_index = hiphase.phased_vcfs[0].data_index,
			reference = reference.fasta.data,
			reference_index = reference.fasta.data_index,
			exclude_bed = reference.hificnv_exclude_bed.data,
			exclude_bed_index = reference.hificnv_exclude_bed.data_index,
			expected_bed_male = reference.hificnv_expected_bed_male,
			expected_bed_female = reference.hificnv_expected_bed_female,
			output_prefix = "hificnv",
			runtime_attributes = default_runtime_attributes
	}

	output {
		# per movie stats, alignments
		Array[File] bam_stats = pbmm2_align.bam_stats
		Array[File] read_length_summary = pbmm2_align.read_length_summary
		Array[File] read_quality_summary = pbmm2_align.read_quality_summary
		Array[IndexData] aligned_bams = aligned_bam

		# phased_vcfs ouput order from HiPhase: small variants, SVs, TRGT

		# per sample structural variant signatures and calls
		IndexData phased_sv_vcf = hiphase.phased_vcfs[1]
		Array[File] svsigs = pbsv_discover.svsig

		# per sample small variant calls
		IndexData phased_small_variant_vcf = hiphase.phased_vcfs[0]
		IndexData small_variant_gvcf = deepvariant.gvcf
		File small_variant_vcf_stats = bcftools.stats
		File small_variant_roh_out = bcftools.roh_out
		File small_variant_roh_bed = bcftools.roh_bed

		# per sample phasing stats and haplotagged alignments
		File hiphase_stats = hiphase.hiphase_stats
		File hiphase_blocks = hiphase.hiphase_blocks
		File hiphase_haplotags = select_first([hiphase.hiphase_haplotags])
		IndexData merged_haplotagged_bam = haplotagged_bam
		File mosdepth_summary = mosdepth.summary
		File mosdepth_region_bed = mosdepth.region_bed

		# per sample trgt outputs
		IndexData trgt_repeat_vcf = hiphase.phased_vcfs[2]
		IndexData trgt_spanning_reads = {"data": trgt.spanning_reads, "data_index": trgt.spanning_reads_index}

		# per sample cpg outputs
		Array[File] cpg_pileup_beds = cpg_pileup.pileup_beds
		Array[File] cpg_pileup_bigwigs = cpg_pileup.pileup_bigwigs

		# per sample paraphase outputs
		File paraphase_output_json = paraphase.output_json
		IndexData paraphase_realigned_bam = {"data": paraphase.realigned_bam, "data_index": paraphase.realigned_bam_index}
		File? paraphase_vcfs = paraphase.paraphase_vcfs

		# per sample hificnv outputs
		IndexData hificnv_vcf = {"data": hificnv.cnv_vcf, "data_index": hificnv.cnv_vcf_index}
		File hificnv_copynum_bedgraph = hificnv.copynum_bedgraph
		File hificnv_depth_bw = hificnv.depth_bw
		File hificnv_maf_bw = hificnv.maf_bw
	}

	parameter_meta {
		sample: {help: "Sample information and associated data files"}
		reference: {help: "Reference genome data"}
		deepvariant_version: {help: "Version of deepvariant to use"}
		custom_deepvariant_model_tar: {help: "Optional deepvariant model to use"}
		default_runtime_attributes: {help: "Default RuntimeAttributes; spot if preemptible was set to true, otherwise on_demand"}
	}
}

task pbmm2_align {
	input {
		String sample_id
		File bam

		File reference
		File reference_index
		String reference_name

		RuntimeAttributes runtime_attributes
	}

	String movie = basename(bam, ".bam")

	Int threads = 24
	Int mem_gb = ceil(threads * 4)
	Int disk_size = ceil((size(bam, "GB") + size(reference, "GB")) * 4 + 20)

	command <<<
		set -euo pipefail

		pbmm2 --version

		pbmm2 align \
			--num-threads ~{threads} \
			--sort-memory 4G \
			--preset HIFI \
			--sample ~{sample_id} \
			--log-level INFO \
			--sort \
			--unmapped \
			~{reference} \
			~{bam} \
			~{sample_id}.~{movie}.~{reference_name}.aligned.bam

		# movie stats
		extract_read_length_and_qual.py \
			~{bam} \
		> ~{sample_id}.~{movie}.read_length_and_quality.tsv

		awk '{{ b=int($2/1000); b=(b>39?39:b); print 1000*b "\t" $2; }}' \
			~{sample_id}.~{movie}.read_length_and_quality.tsv \
			| sort -k1,1g \
			| datamash -g 1 count 1 sum 2 \
			| awk 'BEGIN {{ for(i=0;i<=39;i++) {{ print 1000*i"\t0\t0"; }} }} {{ print; }}' \
			| sort -k1,1g \
			| datamash -g 1 sum 2 sum 3 \
		> ~{sample_id}.~{movie}.read_length_summary.tsv

		awk '{{ print ($3>50?50:$3) "\t" $2; }}' \
				~{sample_id}.~{movie}.read_length_and_quality.tsv \
			| sort -k1,1g \
			| datamash -g 1 count 1 sum 2 \
			| awk 'BEGIN {{ for(i=0;i<=60;i++) {{ print i"\t0\t0"; }} }} {{ print; }}' \
			| sort -k1,1g \
			| datamash -g 1 sum 2 sum 3 \
		> ~{sample_id}.~{movie}.read_quality_summary.tsv
	>>>

	output {
		File aligned_bam = "~{sample_id}.~{movie}.~{reference_name}.aligned.bam"
		File aligned_bam_index = "~{sample_id}.~{movie}.~{reference_name}.aligned.bam.bai"
		File bam_stats = "~{sample_id}.~{movie}.read_length_and_quality.tsv"
		File read_length_summary = "~{sample_id}.~{movie}.read_length_summary.tsv"
		File read_quality_summary = "~{sample_id}.~{movie}.read_quality_summary.tsv"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/pbmm2@sha256:ed9dcb4db98c81967fff15f50fca89c8495b1f270eee00e9bec92f46d14d7e2f"
		cpu: threads
		memory: mem_gb + " GB"
		disk: disk_size + " GB"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: runtime_attributes.preemptible_tries
		maxRetries: runtime_attributes.max_retries
		awsBatchRetryAttempts: runtime_attributes.max_retries
		queueArn: runtime_attributes.queue_arn
		zones: runtime_attributes.zones
	}
}

task bcftools {
	input {
		File vcf

		String? stats_params
		File reference

		RuntimeAttributes runtime_attributes
	}

	String vcf_basename = basename(vcf, ".vcf.gz")

	Int threads = 2
	Int reference_size = if (defined(reference)) then ceil(size(reference, "GB")) else 0
	Int disk_size = ceil((size(vcf, "GB") + reference_size) * 2 + 20)

	command <<<
		set -euo pipefail

		bcftools --version

		bcftools stats \
			--threads ~{threads - 1} \
			~{stats_params} \
			~{"--fasta-ref " + reference} \
			~{vcf} \
		> ~{vcf_basename}.vcf.stats.txt

		bcftools roh \
			--threads ~{threads - 1} \
			--AF-dflt 0.4 \
			~{vcf} \
		> ~{vcf_basename}.bcftools_roh.out

		echo -e "#chr\\tstart\\tend\\tqual" > ~{vcf_basename}.roh.bed
		awk -v OFS='\t' '$1=="RG" {{ print $3, $4, $5, $8 }}' \
			~{vcf_basename}.bcftools_roh.out \
		>> ~{vcf_basename}.roh.bed
	>>>

	output {
		File stats = "~{vcf_basename}.vcf.stats.txt"
		File roh_out = "~{vcf_basename}.bcftools_roh.out"
		File roh_bed = "~{vcf_basename}.roh.bed"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/bcftools@sha256:46720a7ab5feba5be06d5269454a6282deec13060e296f0bc441749f6f26fdec"
		cpu: threads
		memory: "4 GB"
		disk: disk_size + " GB"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: runtime_attributes.preemptible_tries
		maxRetries: runtime_attributes.max_retries
		awsBatchRetryAttempts: runtime_attributes.max_retries
		queueArn: runtime_attributes.queue_arn
		zones: runtime_attributes.zones
	}
}
