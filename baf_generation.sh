#!/bin/bash

#SBATCH -N1 -n4 -t 2-00:00:00 --job-name RNA-Seq --mail-type=END --mail-user=hthirima@fredhutch.org -A XXXX

module load SAMtools/1.16.1-GCC-11.2.0

#generate a text file with a list of samplenames and save it as accessions.txt

if [ "${accessionId}" == "" ]
then
  echo "ERROR: No accession ID given- exiting" >&2
  exit 1
fi

sampleName=${accessionId}

echo $sampleName

bamfiles=~/MeningData2/GSE85133
outputdir=~/CASPER_BAF/GSE85133

BAFExtract=~/Analysis17_CNV/BAFExtract-master/bin/BAFExtract
hg38list=~/Analysis17_CNV/hg38.list
hg38_genome_fasta_pileup_dir=~/Analysis17_CNV/hg38_genome_fasta_pileup_dir/hg38

cd $outputdir
mkdir ${sampleName}
cd $outputdir/${sampleName}
finaloutput=$outputdir/${sampleName}
#output_baf_file=$outputdir/$dataset/${sampleName} 

echo "Extract BAF values from RNA-Seq bam files"

cd $bamfiles/${sampleName}_hg38

samtools view ${sampleName}.bam | $BAFExtract -generate_compressed_pileup_per_SAM stdin $hg38list $finaloutput 50 0; $BAFExtract -get_SNVs_per_pileup $hg38list $finaloutput $hg38_genome_fasta_pileup_dir 20 4 0.1 ${sampleName}.baf

echo "Job completed"
