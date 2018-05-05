# Prerequisites

This program was developed for Python 3.X versions. Python prerequisites include Biopython and skbio. The skbio package is currently only available on Unix/POSIX/MacOS/OS X operating systems which limits the operation of this program to these platforms.

# About

deGRIT attempts to remove (or de-grit) indel errors from a genome assembly utilising RNAseq transcript alignments. This program uses information from a previous gene annotation run from PASA in order to reduce the false positive rate by focusing changes specifically on established exonic regions. An additional "gene rescuing" module is optionally available to allow correction of indels in regions that are not part of currently annotated genes; theoretically this may be more error prone, but in practice more strict measures are enacted to ensure that any changes should be valid (and, because these weren't predicted as exons, changes here won't adversely impact your annotation).

# Program inputs

This program is specifically designed to work with outputs from a few programs. 

Firstly, the gene annotation GFF3 should be in a format produced by PASA. In my project, this file was the result of EVidenceModeler combination of BRAKER and PASA annotations which were subsequently improved with PASA twice, then augmented with tRNA and rRNA annotations using scripts as part of another of my repositories (https://github.com/zkstewart/Genome_analysis_scripts). However, any output that is formatted similar to how PASA formats a GFF3 file should work. An example files is provided in this repository. Pay attention to the comment line formatting where gene names are found, as well as the way in which multiple mRNA isoforms are presented and the ordering of positive and negative orientation lines.

Secondly, the nucleotide transcript alignment file is expected to be produced by GMAP using the settings that PASA automatically uses when running GMAP itself. In practice, this can be easily emulated with a command that resembles this.

```/home/dir_to/PASApipeline/scripts/..//scripts/run_spliced_aligners.pl --aligners gmap --genome $GENOME --transcripts $TRANSCRIPTOME -I 500000 -N 1 --CPU 2```

Unlike when running PASA, however, you should only align predicted CDS regions to the genome, not the full transcripts including UTR. You may extract CDS regions from a transcriptome using TransDecoder, EvidentialGene, or whatever pipeline you most prefer. BLAT alignment will work, but from my testing it does not work as well as GMAP. It is probable that some of the values BLAT produces (like identity or something else) differ from GMAP, and since I've built and tested my program around GMAP, the values I have are fine tuned to accept GMAP alignments instead of BLAT.

Finally, you should also provide the FASTA files for your genome assembly and the CDS nucleotide file you used for GMAP alignment. 

# Program use

This program accepts arguments on the command-line. When calling the program with the help argument

```python3 deGRIT.py -h```

The following information will be presented.

```
usage: deGRIT.py [-h] [-an GFF3FILE] [-gen GENOMEFILE] [-gm GMAPFILE]
                 [-tr TRANSCRIPTOMEFILE] [-o OUTPUTFILENAME] [-r] [-fo] [-v]
                 [-l]

deGRIT.py aims to improve the ability to reannotate gene models. In order to
work, this program requires a gff3 file of gene annotations alongside its
respective genome fasta file in addition to a gff3 file of transcript
alignments with its respective transcriptome fasta file. These files will be
used to compare the genome sequence to the aligned transcript sequence to
identify any occurrences of indel errors. By correcting these indels,
reannotation of gene models can take place which will provide more accurate
results. Note: This program is designed to work with CDS regions from
transcripts; this reduces the chance of falsely interrupting a reading frame
with an edit. You can predict the CDS region using TransDecoder or
EvidentialGene.

optional arguments:
  -h, --help            show this help message and exit
  -an GFF3FILE, --annotation GFF3FILE
                        Input gff3 gene annotation file name
  -gen GENOMEFILE, --genomefile GENOMEFILE
                        Input genome contig fasta file name
  -gm GMAPFILE, --gmap GMAPFILE
                        Input gff3 gmap transcript alignment file name
  -tr TRANSCRIPTOMEFILE, --transfile TRANSCRIPTOMEFILE
                        Input nucleotide transcriptome fasta file name (this
                        is the same transcript file used for GMAP alignment)
  -o OUTPUTFILENAME, --output OUTPUTFILENAME
                        Output results file name
  -r, --rescue_genes    Optionally perform extended gene model rescue module
                        (this is recommended)
  -fo, --force          Optionally allow the program overwrite existing files
                        at your own risk
  -v, --verbose         Print program details to terminal
  -l, --log             Additionally produce a detailed logging file as output
  ```
  
As mentioned above, you must provide your input file arguments (-an, -gen, -gm, -tr) in addition to an output file name (-o) wherein VCF-like formatted information of genome edit positions will be presented. This file can be used by the deGRIT_edit.py program to produce an output modified genome with the insertions and deletions suggested by deGRIT.

Optional tags can be used to change program behaviour.

## Gene rescuing module
The "gene rescuing" module can be activated by -r or --rescue_genes. The main part of this program uses predicted exon boundaries to narrow down spurious transcript matches (such as by closely related genes which don't originate from the exon currently being analysed) which enables greater confidence in the selection of transcripts for alignment. This module does not have this benefit which theoretically may make it more error prone. To counteract this, I use a slightly more stringent identity cutoff and modified behaviour when detecting indels which I believe remedies this problem. Thus, while this module is not activated by default it is recommended since it has been shown to fix previously unannotated exons which will allow these exons to be incorporated into PASA models in subsequent reannotation.

## Force
This argument allows the program to overwrite previous VCF-style results files. This should not be used unless when testing the program. This program handles this scenario by temporarily moving and renaming the file to the directory where this script is being run, then it will delete the file when the program successfully completes. Thus, if the user makes a mistake they can kill the process and retrieve the original file before program completion.

## Verbose
As the name suggests, this makes the program produce ongoing information in the terminal window. The results for each gene model are printed and overall program progress updates are provided.

## Log
This argument will result in a log file being produced for the current run. The log is a tabular file providing exon-by-exon information of the best transcript alignment, its alignment coordinates against the genome, the start and stop positions of the transcript segment which aligns, and any indels suggested by the transcript alignment. This is useful if you want to manually validate any changes suggested by the program. Logging does result in a slight performance hit as we need to calculate and reformat some of the information for human-readability.
