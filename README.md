# Prerequisites

This program was developed for Python 3.X versions. Python prerequisites include Biopython and skbio. The skbio package is only available on Linux operating systems which limits the operation of this program to this platform.

# About

deGRIT attempts to remove (or de-grit) indel errors from a genome assembly utilising RNAseq transcript alignments. This program utilises information from a previous gene annotation run from PASA in order to reduce the false positive rate by focusing changes specifically on established exonic regions. An additional "gene rescuing" module is optionally available to allow correction of indels in regions that are not part of currently annotated genes; theoretically this may be more error prone, but in practice more strict measures are enacted to ensure that any changes should be valid (and, because these weren't predicted as exons, changes here won't adversely impact your annotation).

# Program inputs

This program is specifically designed to work with outputs from a few programs. 

Firstly, the gene annotation gff3 should be in a format produced by PASA. In my project, this file was the result of EVidenceModeler combination of BRAKER and PASA annotations which were subsequently improved with PASA twice, then augmented with tRNA and rRNA annotations using scripts as part of another of my repositories (https://github.com/zkstewart/Genome_analysis_scripts). 

Secondly, the nucleotide transcript alignment file is expected to be produced by GMAP using the settings that PASA automatically uses when running GMAP itself. In practice, this can be easily emulated with a command that resembles this.

```/home/../PASApipeline/scripts/..//scripts/run_spliced_aligners.pl --aligners gmap --genome $GENOME --transcripts $TRANSCRIPTOME -I 500000 -N 1 --CPU 2```

Unlike when running PASA, however, you should only align predicted CDS regions to the genome, not the full transcripts including UTR. You may extract CDS regions from a transcriptome using TransDecoder, EvidentialGene, or whatever pipeline you most prefer. BLAT alignment will work, but from my testing it does not work as well as GMAP. It is probable that some of the values BLAT produces (like identity or something else) differ from GMAP, and since I've built and tested my program around GMAP, the values I have are fine tuned to accept GMAP alignments instead of BLAT.

Finally, you should also provide the FASTA files for your genome assembly and the CDS nucleotide file you used for GMAP alignment. 
