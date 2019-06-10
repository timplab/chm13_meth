#!/bin/bash

# for version 4 mapping that came in cram format

ref=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v4_assembly/chm13.draft_v0.4.fasta
cram=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/rel2_to_v0.4.cram
outdir=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v4_assembly

if [ "$1" == "install" ]; then
	wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
	tar -xjf samtools-1.9.tar.bz2
	cd samtools-1.9
	./configure --enable-plugins --enable-libcurl
	make all all-htslib
	sudo make install install-htslib
	sudo cp samtools /usr/local/bin
fi

if [ "$1" == "subset" ]; then
samtools view -b -T chm13.draft_v0.4.fasta rel2_to_v0.4.cram "chrX_fixedBionanoSV_centromereV3:57012643-61034379" > cenx.bam
samtools sort cenx.bam -o cenx.bam
samtools index cenx.bam
samtools view -b -T chm13.draft_v0.4.fasta rel2_to_v0.4.cram "chrX_fixedBionanoSV_centromereV3:112000000-114000000" > DXZ4.bam
samtools sort DXZ4.bam -o DXZ4.bam
samtools index DXZ4.bam
fi

if [ "$1" == "meth" ]; then
	python3 /home/gmoney/repos/nanopore-methylation-utilities/mtsv2bedGraph.py -i methylation.tsv |\
		sort -k1,1 -k2,2n | bgzip > methylation.bed.gz
	tabix -p bed methylation.bed.gz

fi

if [ "$1" == "bam" ]; then
	python3 /home/gmoney/repos/nanopore-methylation-utilities/convert_bam_for_methylation.py --verbose -b cenx.bam \
		-c methylation.bed.gz -f chm13.draft_v0.4.fasta -w chrX_fixedBionanoSV_centromereV3:57012643-61034379 |\
		samtools sort -o rel2_to_v0.4_cenX_meth.bam
	samtools index rel2_to_v0.4_cenX_meth.bam
        python3 /home/gmoney/repos/nanopore-methylation-utilities/convert_bam_for_methylation.py --verbose -b DXZ4.bam \
                -c methylation.bed.gz -f chm13.draft_v0.4.fasta -w chrX_fixedBionanoSV_centromereV3:112000000-114000000 |\
                samtools sort -o rel2_to_v0.4_DXZ4_meth.bam
        samtools index rel2_to_v0.4_DXZ4_meth.bam
fi
if [ "$1" == "all" ]; then
#samtools view -b -T $ref $cram "chrX_fixedBionanoSV_centromereV3" > ${outdir}/chrx.bam
#samtools sort ${outdir}/chrx.bam -o ${outdir}/chrx.bam
#samtools index ${outdir}/chrx.bam
    python3 /home/gmoney/repos/nanopore-methylation-utilities/convert_bam_for_methylation.py --verbose -b ${outdir}/chrx.bam \
                -c  ${outdir}/methylation.bed.gz -f $ref |\
                samtools sort -o ${outdir}/rel2_to_v0.4_new_chrx_meth.bam
        samtools index ${outdir}/rel2_to_v0.4_chrx_new_meth.bam
fi

if [ "$1" == "bismark" ]; then
	python3 /home/gmoney/repos/nanopore-methylation-utilities/parseMethylbed.py frequency -i ${outdir}/methylation.bed.gz

fi
if [ "$1" == "bedgraph" ]; then
	bedtools genomecov -ibam ${outdir}/chrx.bam -bg > ${outdir}/v4_all.chrX.bg
	bedtools genomecov -ibam ${outdir}/rel2_to_v0.4_chrx_meth.bam -bg > ${outdir}/v4_meth.chrX.bg

fi

