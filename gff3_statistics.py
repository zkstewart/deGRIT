#! python3

# gff3_statistics.py

# This program provides comprehensive statistical output including histograms
# for graphic purposes to analyse the current gene annotation. This is
# primarily intended for comparison of gene annotations before and after
# running the deGRIT program to see how indel correction has affected the
# gene annotation. Ideally, in a genome with many indel errors, if deGRIT
# has improved the genome we would expect to see (probably) less gene models
# in total with a longer average gene length, more exons and introns per gene,
# and longer gaps inbetween each gene model. This program will facilitate this
# comparison and act as a validation that deGRIT has improved your genome.

# Load packages
import os, re, argparse, copy
import pandas as pd
from Bio import SeqIO
from ncls import NCLS
from statistics import median, mean

### Various functions to perform operations throughout the program
def gmap_exon_finder(ncls, gmapLoc, model, coordIndex):
        start = int(model[0][coordIndex].split('-')[0])
        stop = int(model[0][coordIndex].split('-')[1])
        overlaps = ncls.find_overlap(start, stop)                                                       # Although our ncls is 1-based, find_overlap acts as a range and is thus 0-based. We need to +1 to the stop to offset this.
        # Figure out if we have any hits on the same contig / have perfect matching boundaries
        hit = 'n'
        same = 'n'
        for entry in overlaps:
                contigId = gmapLoc[entry[2]]
                if contigId == model[2]:
                        hit = 'y'
                        if entry[0] == start and entry[1] == stop + 1:                                  # Remember that ncls was built to be 1-based, so entry[1] is +1 to the original position.
                                same = 'y'
                                break
        # Return our same value
        return hit, same

## CORE FUNCTIONS ##
def gmap_parse_ncls(gmapFile, cutoff):
        gmapLoc = {}
        starts = []
        ends = []
        ids = []
        ongoingCount = 0
        with open(gmapFile, 'r') as fileIn:
                for line in fileIn:
                        # Skip unneccessary lines
                        if line.startswith('#'):
                                continue
                        sl = line.split('\t')
                        if sl[2] != 'cDNA_match':                                                        # I don't think any other type of line is present in a GMAP gff3 file produced with PASA's settings, but this could potentially future proof the script?
                                continue
                        # Get details from line including start, stop, and orientation
                        contigID = sl[0]
                        contigStart = int(sl[3])
                        contigStop = int(sl[4])
                        identity = float(sl[5])
                        if identity < cutoff:                                                           # Speed up program by only holding onto hits that will pass our cutoff check.
                                continue
                        # Add to our NCLS                                                               # We index using ranges since it provides an easy way to retrieve GMAP matches by coordinates. Since these coordinates aren't unique, we filter any results returned by their contig ID.
                        starts.append(contigStart)
                        ends.append(contigStop+1)                                                       # NCLS indexes 0-based, so +1 to make this more logically compliant with gff3 1-based system.
                        ids.append(ongoingCount)
                        gmapLoc[ongoingCount] = contigID
                        ongoingCount += 1
        # Build the NCLS object
        starts = pd.Series(starts)
        ends = pd.Series(ends)
        ids = pd.Series(ids)
        ncls = NCLS(starts.values, ends.values, ids.values)
        return ncls, gmapLoc

def exon_count(nuclDict):
        exonCount = []
        exonSize = []
        for key, model in nuclDict.items():
                exonCount.append(len(model[0]))
                for coord in model[0]:
                        coordSplit = coord.split('-')
                        exonSize.append(int(coordSplit[1]) - int(coordSplit[0]) + 1)                    # +1 since both numbers are 1-based; 1-1 == 0, but in reality it indicates an exon length of 1.
        # Basic statistics
        maxExonCount = max(exonCount)
        medExonCount = median(exonCount)
        meanExonCount = round(mean(exonCount), 3)
        maxExonSize = max(exonSize)
        medExonSize = median(exonSize)
        meanExonSize = round(mean(exonSize), 3)
        # Return values
        return maxExonCount, medExonCount, meanExonCount, maxExonSize, medExonSize, meanExonSize

def intergenic_spacing(nuclDict):
        # Parse the nuclDict and make ordered lists for contig gene annotations
        nuclLists = {}
        for key, model in nuclDict.items():
                # Derive coordinates
                if model[1] == '+':
                        start = int(model[0][0].split('-')[0])
                        stop = int(model[0][-1].split('-')[1])
                else:
                        start = int(model[0][-1].split('-')[0])
                        stop = int(model[0][0].split('-')[1])
                # Add to our list
                if model[2] not in nuclLists:
                        nuclLists[model[2]] = [[start, stop, model[1]]]
                else:
                        nuclLists[model[2]].append([start, stop, model[1]])
        # Prune out overlaps so we can look specifically at intergenic distances
        intergenDist = []
        for key, lst in nuclLists.items():
                lst.sort(key = lambda x: (x[0], -x[1]))
                # Separate lists into their respective orientations
                fwdList = []
                revList = []
                for entry in lst:
                        if entry[2] == '+':
                                fwdList.append(entry)
                        else:
                                revList.append(entry)
                # Remove overlaps
                breakCond = False
                while True:
                        if breakCond:
                                break
                        breakCond = True
                        for lst in [fwdList, revList]:
                                for i in range(len(lst)-1,-1,-1):
                                        if i != 0:
                                                # Remove redundant coords immediately
                                                if lst[i] == lst[i-1]:
                                                        del lst[i]
                                                        continue
                                                # Compare current to future coord and look for overlap in the same orientation
                                                if int(lst[i][1]) >= int(lst[i-1][0]) and int(lst[i][0]) <= int(lst[i-1][1]):       # If the end of seq1 > start of seq2 AND the start of seq1 < end of seq2, then we know there's overlap regardless of orientation.
                                                        breakCond = False
                                                        # Make a combined longer model
                                                        lst[i-1] = [min([lst[i][0], lst[i-1][0]]), max([lst[i][1], lst[i-1][1]]), lst[i][2]]
                                                        del lst[i]
                                                        continue
                # Calculate the intergenic spacing for this contig
                for lst in [fwdList, revList]:
                        for i in range(len(lst)):
                                if i != len(lst)-1:
                                        intergenDist.append(lst[i+1][0] - lst[i][1] - 1)
        # Calculate statistics
        maxDist = max(intergenDist)
        minDist = min(intergenDist)
        medDist = median(intergenDist)
        meanDist = round(mean(intergenDist), 3)
        # Return results
        return maxDist, minDist, medDist, meanDist

def cdna_parser(gffFile):                                                                               # I've essentially crammed the gff3_to_fasta.py script in here since we need to parse the gff3 file to get the CDS regions to perform the CDS merging and find out if we get a proper gene model.
        def group_process(currGroup):
                full_mrnaGroup = []                                                                     # This will hold processed mRNA positions.
                mrnaGroup = []                                                                          # This will be a temporary storage for mRNA lines.
                for entry in currGroup:
                        # Handle the first line in the group: we just want the gene ID
                        if entry[2] == 'gene':
                                geneID = idRegex.search(entry[8]).group(1)
                        # Handle mRNA lines: this will start a subgroup corresponding to the mRNA
                        elif entry[2] == 'mRNA':
                                if mrnaGroup == []:                                                     # i.e., if this is the first mRNA line in this gene group, we just need to start building it.
                                        mrnaGroup.append(entry)
                                else:                                                                   # i.e., there is more than one mRNA in this gene group, so we need to process the group we've built then initiate a new one.
                                        # Process current mrnaGroup
                                        for subentry in mrnaGroup:
                                                if subentry[2] == 'mRNA':
                                                        full_mrnaGroup.append([idRegex.search(subentry[8]).group(1), []])
                                                elif subentry[2] == 'CDS':
                                                        coords = subentry[3] + '-' + subentry[4]        # +1 here to make Python act 1-based like gff3 format.
                                                        full_mrnaGroup[-1][-1].append(coords)
                                        # Initiate new mrnaGroup
                                        full_mrnaGroup[-1] += [subentry[0],subentry[6]]                 # Append contig ID and orientation.
                                        mrnaGroup = [entry]
                        else:
                                mrnaGroup.append(entry)
                # Process the mrnaGroup that's currently sitting in the pipe (so to speak)
                for subentry in mrnaGroup:
                        if subentry[2] == 'mRNA':
                                full_mrnaGroup.append([idRegex.search(subentry[8]).group(1), []])
                        elif subentry[2] == 'CDS':
                                coords = subentry[3] + '-' + subentry[4]                                # +1 here to make Python act 1-based like gff3 format.
                                full_mrnaGroup[-1][-1].append(coords)
                full_mrnaGroup[-1] += [subentry[0],subentry[6]]                                         # Append contig ID and orientation.
                # Put info into the coordDict and move on
                gffCoordDict[geneID] = full_mrnaGroup
                
        idRegex = re.compile(r'ID=(.+?);')
        currGroup = []
        gffCoordDict = {}
        with open(gffFile, 'r') as fileIn:
                for line in fileIn:
                        # Skip filler lines
                        if line == '\n' or line.startswith('#'):
                                continue
                        # Get details
                        sl = line.rstrip('\n').split('\t')
                        lineType = sl[2]
                        # Building gene group/process it
                        if lineType == 'gene':
                                if currGroup == []:
                                        # First iteration: just play it cool, add the sl to the group
                                        currGroup.append(sl)
                                        continue
                                else:
                                        # Process group if we're encountering a new group
                                        group_process(currGroup)
                                        currGroup = [sl]
                        elif lineType == 'rRNA' or lineType == 'tRNA':                                  # Skip lines that aren't coding.
                                continue
                        else:
                                # Keep building group until we encounter another 'gene' lineType
                                currGroup.append(sl)
                # Process the last mrnaGroup
                group_process(currGroup)
        nuclDict = {}
        nuclDictNoIsos = {}
        for key, value in gffCoordDict.items():
                ## Full gene annotation dictionary
                for mrna in value:                                                                      # This just reorganises the results a bit. Only reason I'm doing this is because I'm reusing bits of old code, and it's easier to just take a small efficiency hit here to not bother rejigging the whole thing leading up to here.
                        nuclDict[mrna[0]] = [mrna[1], mrna[3], mrna[2], mrna[0]]
                ## Only representative isoform dictionary
                longestMrna = ['', 0]
                for mrna in value:
                        mrnaLen = 0
                        for pair in mrna[1]:
                                coords = pair.split('-')
                                mrnaLen += (int(coords[1]) - int(coords[0]) + 1)
                        if mrnaLen > longestMrna[1]:
                                longestMrna = [mrna, mrnaLen]
                value = [longestMrna[0]]
                for mrna in value:                                                                      # This just reorganises the results a bit. Only reason I'm doing this is because I'm reusing bits of old code, and it's easier to just take a small efficiency hit here to not bother rejigging the whole thing leading up to here.
                        nuclDictNoIsos[mrna[0]] = [mrna[1], mrna[3], mrna[2], mrna[0]]
        return nuclDict, nuclDictNoIsos

def gmap_parse_models(args, cutoff, transRecords, minSeqLen, nucleotideDict):
        # Re-index our nuclDict into a format capable of comparison to our gmapLoc dictionary
        nuclNcls, nuclModels = reindex_nucldict(nucleotideDict)
        # Parse the gmapFile and build models
        gmapModels = {}
        with open(args.gmapFile, 'r') as fileIn:
                for line in fileIn:
                        # Skip unneccessary lines
                        if line.startswith('#'):
                                continue
                        sl = line.split('\t')
                        if sl[2] != 'cDNA_match':                                                       # I don't think any other type of line is present in a GMAP gff3 file produced with PASA's settings, but this could potentially future proof the script?
                                continue
                        # Get details from line
                        contigID = sl[0]
                        pathID = sl[8].split(';')[0][3:]                                                # Start from 3 to get rid of 'ID='.
                        if not sl[8].startswith('ID=blat.proc'):                                        # May as well make the code compatible with BLAT.
                                geneID = sl[8].split(';')[1][5:]                                        # Start from 5 to get rid of 'Name='.
                        else:
                                geneID = sl[8].split(';')[1].split()[0][7:]                             # Start from 5 to get rid of 'Target='.
                        contigStart = sl[3]
                        contigStop = sl[4]
                        identity = float(sl[5])
                        orient = sl[6]
                        # Skip models that won't be accepted
                        if identity < cutoff:
                                gmapModels[pathID] = 'not_good'
                                continue
                        if pathID in gmapModels:
                                if gmapModels[pathID] == 'not_good':
                                        continue
                        # Get further details if this model is still OK
                        if not sl[8].startswith('ID=blat.proc'):
                                transStart = sl[8].split(';')[2].split()[1]
                                transStop = sl[8].split(';')[2].split()[2]
                        else:
                                transStart = sl[8].split(';')[1].split()[1]
                                transStop = sl[8].split(';')[1].split()[2]
                        # Add to our dictionary                                                         # We index using ranges since it provides an easy way to retrieve GMAP matches by coordinates. Since these coordinates aren't unique, we filter any results returned by their contig ID.
                        if pathID not in gmapModels:
                                gmapModels[pathID] = [[contigStart + '-' + contigStop], orient, contigID, geneID, [transStart + '-' + transStop]]
                        else:
                                gmapModels[pathID][0].append(contigStart + '-' + contigStop)
                                gmapModels[pathID][4].append(transStart + '-' + transStop)
        # Set up our values for maintaining memory of annotated exon regions
        cidDict = {}
        unrescuedExons = 0
        unrescuedGenes = 0
        for key, value in gmapModels.items():
                ## Check 1: Was this marked as 'not_good' earlier? (this means the gene model contained at least one exon which didn't meet out cutoff)
                if value == 'not_good':
                        continue
                ## Check 2: Does this model have more than 1 exon? (single exon genes may be pseudogenes)
                if len(value[0]) < 2:
                        continue
                ## Check 3: Is this model long enough to reasonably be annotated? (short genes are intentionally ignored by many gene prediction programs)
                # Get the transcript details
                record = transRecords[value[3]]
                seq = str(record.seq)
                if len(seq) < minSeqLen:
                        continue
                ## Check 4: Does the full CDS align to the genome? (fragmentary alignments mean this CDS doesn't belong to this genomic region or that the genome sequence itself is misassembled; in both cases, we can't annotate this gene model correctly)
                # Format the alignment details
                seqSet = set(range(1, len(seq)+1))
                alignSet = set()
                for coord in value[4]:
                        coordSplit = list(map(int, coord.split('-')))
                        alignSet = alignSet.union(set(range(coordSplit[0], coordSplit[1] + 1)))
                # Check for overall number of gaps, and check for long stretches of gaps
                gaps = list(seqSet - alignSet)
                gaps.sort()
                gapPerc = len(gaps) / len(seq)
                if gapPerc >= 0.05:                                                                     # Missing 5% of the sequence isn't much, but it's enough to raise concern.
                        continue
                gapStretch = 1
                check4pass = 'y'
                for x in range(len(gaps)-1):
                        if gaps[x] + 1 == gaps[x+1]:
                                gapStretch += 1
                        else:
                                gapStretch = 1
                        if gapStretch >= 5:                                                             # A 5bp gap also isn't much, but we don't want to make a "putative gene model" which doesn't match the transcript nearly perfectly.
                                check4pass = 'n'
                                break
                if check4pass == 'n':
                        continue
                ## Obtain information from this transcript, determing whether it is "unrescued" and to what extent
                entirelyUnrescued = 'y'
                for x in range(len(value[0])):
                        # Query the cdna ncls object for overlaps
                        coordSplit = value[0][x].split('-')
                        overlaps = nuclNcls.find_overlap(int(coordSplit[0]), int(coordSplit[1])+1)
                        nuclHits = []
                        for result in overlaps:
                                nuclHits.append(nuclModels[result[2]])
                        nuclHits = copy.deepcopy(nuclHits)
                        # Query the current ncls object for previously discovered exons
                        if value[2] in cidDict:
                                overlaps = cidDict[value[2]].intersection(set(range(int(coordSplit[0]), int(coordSplit[1])+1)))
                                if overlaps != set():
                                        nuclHits.append([value[2], ''])
                        # Narrow down our overlaps to hits on the same contig
                        for k in range(len(nuclHits)-1, -1, -1):
                                if nuclHits[k][0] != value[2]:
                                        del nuclHits[k]
                        # Save any exons that lack overlaps
                        if nuclHits != []:                                                              # This list won't be empty if it overlaps an established gene model.
                                entirelyUnrescued = 'n'
                                continue
                        else:
                                unrescuedExons += 1
                                # Build the ongoing set object
                                if value[2] not in cidDict:
                                        cidDict[value[2]] = set(range(int(coordSplit[0]), int(coordSplit[1])+1))
                                else:
                                        cidDict[value[2]] = cidDict[value[2]].union(set(range(int(coordSplit[0]), int(coordSplit[1])+1)))
                # If this transcript has entirely novel exons, +1 to this count
                if entirelyUnrescued == 'y':
                        unrescuedGenes += 1
        return unrescuedExons, unrescuedGenes

def reindex_nucldict(nuclDict):
        # Re-index our nuclDict into a NCLS format
        nuclModels = {}
        starts = []
        ends = []
        ids = []
        ongoingCount = 0
        for key, value in nuclDict.items():
                for coord in value[0]:
                        coordSplit = coord.split('-')
                        starts.append(int(coordSplit[0]))
                        ends.append(int(coordSplit[1])+1)
                        ids.append(ongoingCount)
                        nuclModels[ongoingCount] = [value[2], value[3]]
                        ongoingCount += 1
        # Build the NCLS object
        starts = pd.Series(starts)
        ends = pd.Series(ends)
        ids = pd.Series(ids)
        nuclNcls = NCLS(starts.values, ends.values, ids.values)
        # Return this object
        return nuclNcls, nuclModels

def validate_args(args):
        # Validate input file locations
        if not os.path.isfile(args.gff3File):
                print('I am unable to locate the input gff3 gene annotation file (' + args.gff3File + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        elif not os.path.isfile(args.gmapFile):
                print('I am unable to locate the input GMAP transcript alignment gff3 file (' + args.gmapFile + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        # Handle file overwrites
        if os.path.isfile(args.outputFileName):
                print(args.outputFileName + ' already exists. Delete/move/rename this file and run the program again.')
                quit()

### USER INPUT
usage = """%(prog)s analyses a gene annotation gff3 file alongside CDS alignments from GMAP
and the corresponding transcriptome file and calculates some basic statistics of the 
annotation. These include the number of exons and their amount of transcript support, the amount
of transcripts not represented in the annotation file, as well as statistics of exon count/length
in addition to intergenic distances. These statistics should help when evaluating the impact that
deGRIT has had on the annotation by comparing the original to the updated file.
"""
# Reqs
p = argparse.ArgumentParser(description=usage)

p.add_argument("-an", "--annotation", dest="gff3File",
               help="Input gff3 gene annotation file name")
p.add_argument("-gm", "--gmap", dest="gmapFile",
               help="Input gff3 gmap transcript alignment file name")
p.add_argument("-tr", "--transfile", dest="transcriptomeFile",
               help="Input nucleotide transcriptome fasta file name (this is the same transcript file used for GMAP alignment)")
p.add_argument("-o", "--output", dest="outputFileName",
               help="Output results file name")

args = p.parse_args()

# Validate arguments
tmpFileName = validate_args(args)

# Declare values needed for processing
gmapCutoff = 97

# Parse the gff3 file
nuclDict, nuclDictNoIso = cdna_parser(args.gff3File)

# Parse the gmap alignment file for transcript alignment locations
gmapNcls, gmapLoc = gmap_parse_ncls(args.gmapFile, gmapCutoff)

# Parse the transcriptome file
transRecords = SeqIO.to_dict(SeqIO.parse(open(args.transcriptomeFile, 'r'), 'fasta'))

## STATISTIC: Exons covered by transcripts
exonCover = [0, 0, 0]                                                           # Format is [num exons, num covered by transcripts, num perfect boundary matches]
for key, model in nuclDictNoIso.items():                                        # Look at representative isoforms only here since this means our total num exons value can provide a more objective measure of the amount of exons captured in this annotation.
        # Scan through each individual model's exons
        for i in range(len(model[0])):
                # Find GMAP matches and determine whether it is an exact match
                hit, same = gmap_exon_finder(gmapNcls, gmapLoc, model, i)
                # Update exonCover values
                exonCover[0] += 1
                if hit != 'n':
                        exonCover[1] += 1
                if same != 'n':
                        exonCover[2] += 1
percCover = round((exonCover[1] / exonCover[0]) * 100, 3)
percPerfect = round((exonCover[2] / exonCover[1]) * 100, 3)

## STATISTIC: Unannotated exons and genes suggested by transcript alignment
unrescuedExons, unrescuedGenes = gmap_parse_models(args, gmapCutoff, transRecords, 100, nuclDict)               # Most gene annotation programs won't annotate ORFs shorter than 100 bp. This is a reasonable default value to specify here.
                                                                                                                # We also use the full nuclDict with isoform variants since, otherwise, our unrescued gene values wouldn't be entirely correct (we'd classify isoform variants as "unrescued").
## STATISTIC: Exon statistics per gene
maxExonCount, medExonCount, meanExonCount, maxExonSize, medExonSize, meanExonSize = exon_count(nuclDictNoIso)   # Representative isoforms only is more valid here since we'll otherwise get a lot of repeated exons.

## STATISTIC: Intergenic distance between genes
maxDist, minDist, medDist, meanDist = intergenic_spacing(nuclDict)                                              # Use the full isoform dict since we handle overlaps by merging gene models - this provides a slightly more "realistic" view of the intergenic dist if a shorter isoform has an extra exon at its terminal, for example.

# Present results in output file
with open(args.outputFileName, 'w') as fileOut:
        fileOut.write('# Statistics output for ' + os.path.basename(args.gff3File) + ' and ' + os.path.basename(args.gmapFile))
        # Exons covered by transcripts
        fileOut.write('\n# Exons covered by transcripts')
        fileOut.write('\nTotal num of exons part of longest isoforms\t' + str(exonCover[0]))
        fileOut.write('\nExons with transcriptional support\t' + str(exonCover[1]) + '\t' + str(percCover) + '% of total')
        fileOut.write('\nExons with perfect support for boundaries\t' + str(exonCover[2]) + '\t' + str(percPerfect) + '% of exons with support')
        # Unannotated exons and genes suggested by transcript alignment
        fileOut.write('\n# Unannotated exons and genes suggested by transcript alignment')
        fileOut.write('\nTotal num of exons indicated by transcript alignment absent in the annotation file\t' + str(unrescuedExons))
        fileOut.write('\nTotal num of good gene models indicated by transcript alignment absent in the annotation file\t' + str(unrescuedGenes))
        # Exon statistics per gene
        fileOut.write('\n# Exon statistics per gene')
        fileOut.write('\nMost exons in a gene=' + str(maxExonCount) + '\tMedian exon count per gene=' + str(medExonCount) + '\tMean exon count per gene=' + str(meanExonCount))
        fileOut.write('\nLongest exon length (bp)=' + str(maxExonSize) + '\tMedian exon length (bp)=' + str(medExonSize) + '\tMean exon length (bp)=' + str(meanExonSize))
        # Intergenic distance between genes
        fileOut.write('\n# Intergenic distance between genes of the same orientation')
        fileOut.write('\nLargest distance between genes (bp)=' + str(maxDist) + '\tSmallest distance between genes (bp)=' + str(minDist) + '\tMedian distance between genes (bp)=' + str(medDist) + '\tMean distance between genes (bp)=' + str(meanDist))

#### SCRIPT ALL DONE
