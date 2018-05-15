# Prerequisites

This program was developed for Python 3.X versions. Python prerequisites include Biopython, skbio, and ncls (https://github.com/hunt-genes/ncls). The skbio package is currently only available on Unix/POSIX/MacOS/OS X operating systems which limits the operation of this program to these platforms.

# About

deGRIT attempts to remove (or de-grit) indel errors from a genome assembly utilising RNAseq transcript alignments. Genome assemblies using PacBio reads even after polishing using Arrow and Pilon are likely to contain many indel errors which interrupt reading frames and result in gene model fragmentation during annotation. This program was built to address this issue within my own genome data, and is made available for other researchers who may find themselves in a similar situation when working with PacBio assemblies.

deGRIT uses information from a previous gene annotation run from PASA in order to reduce the false positive rate by focusing changes specifically on established exonic regions. An additional "gene rescuing" module is optionally available to allow correction of indels in regions that are not part of currently annotated genes; theoretically this may be more error prone, but in practice the strict measures enforced by deGRIT ensure that almost all changes should specifically fix indel errors only (and, because these weren't predicted as exons, changes here won't adversely impact your annotation).

While it may be possible to address indel errors in a genome by introducing N's to produce a consistent reading frame, I was not satisfied with this option and, to my knowledge, a program capable of doing this for genomic sequences is not currently available. By utilising RNAseq transcripts, we can identify the exact locations where indels occur and, in the case where an insertion must be made into the genome to fix an error, we can use a biologically "real" nucleotide. One drawback is that if your transcriptomic data does not originate from the exact organism you sequenced the genome from (it MUST be from the same species, however), the result may be somewhat of a mosaic whereby SNPs from multiple organisms may be present in the resultant genome. If the purpose of your genome assembly is to produce error-free gene models in a PacBio assembly, this should not be a serious problem, and the benefits of using this program far outweigh this drawback.

# Program inputs

This program is specifically designed to work with outputs from a few programs. 

Firstly, the gene annotation GFF3 should be in a format produced by EVidenceModeler with subsequent PASA updating. In my project, this file was the result of EVidenceModeler combination of BRAKER and PASA annotations which were subsequently improved with PASA twice, then augmented with tRNA and rRNA annotations using scripts as part of another of my repositories (https://github.com/zkstewart/Genome_analysis_scripts). However, any output that is formatted similar to how PASA formats a GFF3 file should work. An example files is provided in this repository. Pay attention to the comment line formatting where gene names are found, as well as the way in which multiple mRNA isoforms are presented and the ordering of positive and negative orientation lines.

Secondly, the nucleotide transcript alignment file is expected to be produced by GMAP using the settings that PASA automatically uses when running GMAP itself. In practice, this can be easily emulated with a command that resembles this.

```/home/dir_to/PASApipeline/scripts/..//scripts/run_spliced_aligners.pl --aligners gmap --genome $GENOME --transcripts $TRANSCRIPTOME -I 500000 -N 1 --CPU 2```

Unlike when running PASA, however, you should only align predicted CDS regions to the genome, not the full transcripts including UTR. You may extract CDS regions from a transcriptome using TransDecoder, EvidentialGene, or whatever pipeline you most prefer. BLAT alignment will work, but from my testing it does not work as well as GMAP. It is probable that some of the values BLAT produces (like identity or something else) differ from GMAP, and since I've built and tested my program around GMAP, the values I have are fine tuned to accept GMAP alignments instead of BLAT.

Finally, you should also provide the FASTA files for your genome assembly and the CDS nucleotide file you used for GMAP alignment.

# Program benchmark

Test 1 - scallop gene model prediction

Scallop predicts gene models based solely on RNAseq alignments against the genome sequence and thus any errors in this sequence will be present in the FASTA file produced when converting gtf to fasta. Within the worst affected genome sequence, deGRIT made the following modifications.

```
Genome editing complete. Basic information provided below.
Edits include gene module rescue positions: Yes
Number of deletion events: 8071
Number of substitution events: 148
Number of insertion events: 976
Total number of contigs: 829
Number of edited contigs: 409
```

BUSCO before: ```Before: C:94.1%[S:39.7%,D:54.4%],F:4.3%,M:1.6%,n:978```

BUSCO after: ```C:96.6%[S:41.3%,D:55.3%],F:1.7%,M:1.7%,n:978```

As demonstrated, in a genome that contains many indel errors, BUSCO scores can be improved noticeably by converting fragmented genes into complete genes.

Test 2 - EVM combination of BRAKER & PASA with PASA updates.

Forthcoming.

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
The "gene rescuing" module can be activated by -r or --rescue_genes. The main part of this program uses predicted exon boundaries to narrow down spurious transcript matches (such as by closely related genes which don't originate from the exon currently being analysed) which enables greater confidence in the selection of transcripts for alignment. 

The gene rescue module instead analyses the input GMAP alignment file and identifies probable gene models using a few heuristic checks. Specifically, it will identify gene models from transcript alignments where all exons have alignment identity >= 97%, the putative model has more than 1 exon, and where the transcript near-completely aligns to the genome (95% of the CDS must align with only small gaps allowed). Any models identified through this that contain exons not covered by annotated models will be extracted and subjected to the main module's processing to identify any potential indel errors that may have prevented this from being annotated originally by BRAKER/PASA/EVM.

While this module is not activated by default (simply because it may be slightly more error prone), it is recommended to use since it has been shown to fix previously unannotated exons which will allow these exons to be incorporated into current PASA models and it will likely result in other gene models that were originally discarded being identified in subsequent reannotation

## Force
This argument allows the program to overwrite previous VCF-style results files. This should not be used unless when testing the program. This program handles this scenario by temporarily moving and renaming the file to the directory where this script is being run, then it will delete the file when the program successfully completes. Thus, if the user makes a mistake they can kill the process and retrieve the original file before program completion.

## Verbose
As the name suggests, this makes the program produce ongoing information in the terminal window. The results for each gene model are printed and overall program progress updates are provided.

## Log
This argument will result in a log file being produced for the current run. The log is a tabular file providing exon-by-exon information of the best transcript alignment, its alignment coordinates against the genome, the start and stop positions of the transcript segment which aligns, and any indels suggested by the transcript alignment. This is useful if you want to manually validate any changes suggested by the program. Logging does result in a slight performance hit as we need to calculate and reformat some of the information for human-readability.

## Verbose and log
Both options will result in the calculation of probable gene overlaps. This process simply detects any gene models that, after editing, have overlapping terminal exons. During reannotation these genes will likely be merged. Other genes that will be merged with intervening introns are not detected by this program since determining this goes beyond the effort this program aims to employ when tracking how gene models have been modified. This simply serves as a useful validation to see that the program is fixing gene models that were previously fragmented.
