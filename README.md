# Prerequisites

This program was developed for Python 3.X versions. Python prerequisites include Biopython, skbio, and ncls (https://github.com/hunt-genes/ncls). The skbio package is currently only available on Unix/POSIX/MacOS/OS X operating systems which limits the operation of this program to these platforms.

# About

deGRIT attempts to remove (or de-grit) indel errors from a genome assembly utilising RNAseq transcript alignments. Genome assemblies using PacBio reads even after polishing using Arrow and Pilon are likely to contain many indel errors which interrupt reading frames and result in gene model fragmentation during annotation. This program was built to address this issue within my own genome data, and is made available for other researchers who may find themselves in a similar situation when working with PacBio assemblies.

deGRIT uses information from a previous gene annotation run from PASA in order to reduce the false positive rate by focusing changes specifically on established exonic regions. An additional "gene rescuing" module is optionally available to allow correction of indels in regions that are not part of currently annotated genes; theoretically this may be more error prone, but in practice the strict measures enforced by deGRIT ensure that almost all changes should specifically fix indel errors only (and, because these weren't predicted as exons, changes here won't adversely impact your annotation). In-frame stop codons can also be fixed by deGRIT, and this is detailed below.

While it may be possible to address indel errors in a genome by introducing N's to produce a consistent reading frame, I was not satisfied with this option and, to my knowledge, a program capable of doing this for genomic sequences is not currently available. By utilising RNAseq transcripts, we can identify the exact locations where indels occur and, in the case where an insertion must be made into the genome to fix an error, we can use a biologically "real" nucleotide. One drawback is that if your transcriptomic data does not originate from the exact organism you sequenced the genome from (it MUST be from the same species, however), the result may be somewhat of a mosaic whereby SNPs from multiple organisms may be present in the resultant genome. If the purpose of your genome assembly is to produce error-free gene models in a PacBio assembly, this should not be a serious problem, and the benefits of using this program far outweigh this drawback.

See https://flowersoftheocean.wordpress.com/2018/07/20/degrit-a-program-to-polish-genome-assemblies-using-transcriptome-data/ for a detailed run through of the way deGRIT works and many other things relating to this program.

# Program inputs

This program is specifically designed to work with outputs from a few programs. 

Firstly, the gene annotation GFF3 should be in a format produced by EVidenceModeler with subsequent PASA updating. In my project, this file was the result of EVidenceModeler combination of BRAKER and PASA annotations which were subsequently improved with PASA twice, then augmented with tRNA and rRNA annotations using scripts as part of another of my repositories (https://github.com/zkstewart/Genome_analysis_scripts). However, any output that is formatted similar to how PASA formats a GFF3 file should work. Example files are provided in this repository (FORTHCOMING). Pay attention to the comment line formatting where gene names are found, as well as the way in which multiple mRNA isoforms are presented and the ordering of positive and negative orientation lines.

Secondly, the nucleotide transcript alignment file is expected to be produced by GMAP using the settings that PASA automatically uses when running GMAP itself. In practice, this can be easily emulated with a command that resembles this.

```/home/dir_to/PASApipeline/scripts/run_spliced_aligners.pl --aligners gmap --genome $GENOME --transcripts $TRANSCRIPTOME -I 500000 -N 1 --CPU 2```

Note that the above PASA command assumes version 2.3.0 of the program is being used. How PASA runs GMAP may change with updates to the software.

Unlike when running PASA, however, you should only align predicted CDS regions to the genome, not the full transcripts including UTR. You may extract CDS regions from a transcriptome using TransDecoder, EvidentialGene, or whatever pipeline you most prefer. BLAT alignment will work, but from my testing it does not work as well as GMAP. It is probable that some of the values BLAT produces (like identity or something else) differ from GMAP, and since I've built and tested my program around GMAP, the values I have are fine tuned to accept GMAP alignments instead of BLAT.

Finally, you should also provide the FASTA files for your genome assembly and the CDS nucleotide file you used for GMAP alignment.

# Program benchmark

## Test 1 - scallop gene model prediction

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

BUSCO before: ```C:94.1%[S:39.7%,D:54.4%],F:4.3%,M:1.6%,n:978```

BUSCO after: ```C:96.6%[S:41.3%,D:55.3%],F:1.7%,M:1.7%,n:978```

As demonstrated, in a genome that contains many indel errors, BUSCO scores can be improved noticeably by converting fragmented genes into complete genes.

## Test 2 - EVM combination of BRAKER & PASA with PASA updates.

Using the above-mentioned edits, I re-ran BRAKER & PASA, used EVM to combine these results, then used PASA to update the gene models two times iteratively. BUSCO scores were computed using the representative isoform of each gene model. These scores from the first annotation to the second after deGRIT are below.

BUSCO before: ```C:94.4%[S:90.2%,D:4.2%],F:2.7%,M:2.9%,n:978```

BUSCO after: ```C:96.1%[S:91.7%,D:4.4%],F:1.5%,M:2.4%,n:978```

As with scallop, fragmented genes have been fixed and, in this case, we have also reduced the number of missing models compared to the first annotation.

Additional gene model statistics were calculated using the gff3_statistics.py script part of this repository. These are presented below.

## Before
```
# Statistics output for cal_smart.rnam-trna.predeGRIT.sorted.gff3 and cal_smart_cds.gmap.spliced_alignments.predeGRIT.gff3
# Exons covered by transcripts
Total num of exons part of longest isoforms	193920
Exons with transcriptional support	167219	86.231% of total
Exons with perfect support for boundaries	144402	86.355% of exons with support
# Unannotated exons and genes suggested by transcript alignment
Total num of exons indicated by transcript alignment absent in the annotation file	15235
Total num of good gene models indicated by transcript alignment absent in the annotation file	3009
# Exon statistics per gene
Most exons in a gene=162	Median exon count per gene=4.0	Mean exon count per gene=6.468
Longest exon length (bp)=14591	Median exon length (bp)=117.0	Mean exon length (bp)=202.943
# Intergenic distance between genes of the same orientation
Largest distance between genes (bp)=403376	Smallest distance between genes (bp)=8	Median distance between genes (bp)=5527	Mean distance between genes (bp)=11727.595
Extra info of smallest distance	Contig ID=utg380_pilon_pilon	Intergenic coordinates=443118-443127
```

## After
```
# Statistics output for cal_smart.rnam-trna.final.sorted.gff3 and gmap.spliced_alignments.gff3
# Exons covered by transcripts
Total num of exons part of longest isoforms	190620
Exons with transcriptional support	164470	86.282% of total
Exons with perfect support for boundaries	151184	91.922% of exons with support
# Unannotated exons and genes suggested by transcript alignment
Total num of exons indicated by transcript alignment absent in the annotation file	13331
Total num of good gene models indicated by transcript alignment absent in the annotation file	2818
# Exon statistics per gene
Most exons in a gene=221	Median exon count per gene=4.0	Mean exon count per gene=6.794
Longest exon length (bp)=14910	Median exon length (bp)=119.0	Mean exon length (bp)=209.929
# Intergenic distance between genes of the same orientation
Largest distance between genes (bp)=299025	Smallest distance between genes (bp)=7	Median distance between genes (bp)=6195.0	Mean distance between genes (bp)=12490.949
Extra info of smallest distance	Contig ID=utg34	Intergenic coordinates=301604-301612
```

Notably, we have reduced the total number of exons part of longest isoforms as expected when combining genes which were split by fragmentation within (what should be) a shared exon. 

Our number of perfectly supported exon boundaries has increased quite a bit and we have less absent exons and gene models.

Our largest distance between genes has decreased indicating the presence of a new gene model within this previous region, and the median and mean distances between genes has increased which shows that we have less gene models in the immediate vicinity of one another (which is often indicative of fragmentation).

The smallest distance between genes remained largely the same, which suggests that this gene or another gene was not able to be fixed by deGRIT. In the case of the "Before" instance, manual inspection revealed that there was no transcript CDS which aligned over the small gap separating these two "fragmented" gene models. Two separate CDS alignments perfectly supported the exon boundaries of these gene models separately, however. Based on a cursory overview of BLAST results, it is likely that these genes should be part of a much larger model which is not present in my transcriptome. In the case suggested by the "After" instance, this gene model entirely lacked transcriptional support, and it is probable that some pseudo-random system in the BRAKER pipeline resulted in another gene model being annotated just upstream of a gene that was part of the original annotation. Manual inspection indicates that a single bp deletion can merge these two models, but the fact that this is a single-exon gene with no transcriptional support suggests some need for caution when suggesting that this is "real".

Nonetheless, these statistics point towards a markedly better annotation with the only limitation being that deGRIT cannot fix errors present in genomic sequence which lack transcript alignments.

# Computional resource requirements

This program is single-threaded and uses minimal memory. From my own use, about 2GB of memory and 25 minutes was required when using the --rescue_genes module (explained below) on genomes approximately 300Mb in size. Not using this functionality speeds up program execution by approximately 20%, but the program is still quite quick (thanks to NCLS!) so this should not be a major concern. Larger input files may use slightly more resources, but it should never get to be more than a standard 8GB laptop could handle.

# Program use

This program accepts arguments on the command-line. When calling the program with the help argument

```python3 deGRIT.py -h```

The following information will be presented.

```
usage: deGRIT.py [-h] [-an GFF3FILE] [-gen GENOMEFILE] [-gm GMAPFILE]
                 [-tr TRANSCRIPTOMEFILE] [-o OUTPUTFILENAME] [-r] [-s] [-fo]
                 [-v] [-l]

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
  -s, --stop_codons     Optionally identify in-frame stop codons and correct
                        these if indicated by transcript evidence.
  -fo, --force          Optionally allow the program overwrite existing files
                        at your own risk
  -v, --verbose         Print program details to terminal
  -l, --log             Additionally produce a detailed logging file as output
  ```
  
As mentioned above, you must provide your input file arguments (-an, -gen, -gm, -tr) in addition to an output file name (-o) wherein VCF-like formatted information of genome edit positions will be presented. This file can be used by the deGRIT_edit.py program to produce an output modified genome with the insertions and deletions suggested by deGRIT.

Optional tags can be used to change program behaviour.

## Gene rescuing module
The "gene rescuing" module can be activated by -r or --rescue_genes. The main part of this program uses predicted exon boundaries to narrow down spurious transcript matches (such as by closely related genes which don't originate from the exon currently being analysed) which enables greater confidence in the selection of transcripts for alignment. 

The gene rescue module additionally analyses the input GMAP alignment file and identifies probable gene models (which weren't originally annotated) using a few heuristic checks. Specifically, it will identify gene models from transcript alignments where all exons have alignment identity >= 97%, the putative model has more than 1 exon, and where the transcript near-completely aligns to the genome (95% of the CDS must align with only small gaps allowed). Any models identified through this that contain exons not covered by annotated models will be extracted and subjected to the main module's processing to identify any potential indel errors that may have prevented this from being annotated originally by BRAKER/PASA/EVM.

While this module is not activated by default (simply because it may be slightly more error prone), it is recommended to use since it has been shown to fix previously unannotated exons which will allow these exons to be incorporated into current PASA models and it will likely result in other gene models that were originally discarded being identified in subsequent reannotation.

## In-frame stop codon fixing
This module can be activated by -s or --stop_codons. By default this program does not make substitutions to the genome that aren't indicated by a gap in alignment (in the case of high identity alignments, this typically means there is an indel, and that is what deGRIT checks for). In cases where your genome sequence incorrectly substitutes a C for a G, you might end up with an in-frame 'TAG' codon that results in a gene model's ORF being falsely terminated. Turning this function on means that transcript alignments will be checked against the genome, and where a substitution is evident in the genome causing an in-frame stop codon which is not present within your transcript, this will be treated as an error and fixed by deGRIT.

As with the gene rescuing module, it is expected that the error rate of this process will be slightly higher than the basic functionality of deGRIT; for example, RNA editing is a process which does occur and an in-frame codon in the genome may be corrected within mature mRNA transcript by the organism's biological machinery. In other words, sometimes in-frame stop codons are really part of the genome sequence. However, from utilising this process on my own genome, enough such errors existed that it was highly unlikely that these substitutions in the genome were real. Instead, it is likely that errors existed within the genome sequence which terminated ORFs earlier than they should, and deGRIT was able to remedy this. 

As such, I do recommend you use this function, but it is a choice to be made when considering your circumstances and whether you actually do need to use this.

## Force
This argument allows the program to overwrite previous VCF-style results files. This should not be used unless when testing the program. This program handles this scenario by temporarily moving and renaming the file to the directory where this script is being run, then it will delete the file when the program successfully completes. Thus, if the user makes a mistake they can kill the process and retrieve the original file before program completion.

## Verbose
As the name suggests, this makes the program produce ongoing information in the terminal window. The results for each gene model are printed and overall program progress updates are provided.

## Log
This argument will result in a log file being produced for the current run. The log is a tabular file providing exon-by-exon information of the best transcript alignment, its alignment coordinates against the genome, the start and stop positions of the transcript segment which aligns, and any indels suggested by the transcript alignment. This is useful if you want to manually validate any changes suggested by the program. Logging does result in a slight performance hit as we need to calculate and reformat some of the information for human-readability, but this makes little difference in reality.

I would recommend that you use this function so you can check where in the genome edits are being made and ensure that the program is working correctly. If it isn't, let me know - it works for me, but my genomes might be quite different to your own.

## Verbose and log
Both options will result in the calculation of probable gene overlaps. This process simply detects any gene models that, after editing, have overlapping terminal exons. During reannotation these genes will likely be merged. Other genes that will be merged with intervening introns are not detected by this program since determining this goes beyond the effort this program aims to employ when tracking how gene models have been modified. This simply serves as a useful validation to see that the program is fixing gene models that were previously fragmented.
