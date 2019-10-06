#! python3

# deGRIT (DEtection and Genomic Rectification of Indels using Transcripts)

# This program attempts to improve (or de-grit) a genome assembly through rectification of indel errors
# utilising RNAseq transcripts. In order to reduce the false positive rate, changes are 
# limited specifically to regions that represent currently annotated genes. An additional
# "gene rescuing" module is optionally available to allow correction of indels in regions 
# that are not part of currently annotated genes; this behaviour is expected to be slightly
# more error prone, but in practice I have not found this module to incorrectly modify the genome
# and thus it may be considered part of the default behavior of this program. Through alignment
# of the transcript and genome sequences, we can find the precise locations where indel errors
# occur and correct these in the genomic sequence, enabling reannotation to occur.

# Load packages
import re, os, argparse, copy, warnings, shutil, parasail
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna, generic_protein
from ncls import NCLS

### Various functions to perform operations throughout the program
def reverse_comp(seq):
        reversedSeq = seq[::-1].lower()
        # Complement characters
        reversedSeq = reversedSeq.replace('a', 'T')
        reversedSeq = reversedSeq.replace('t', 'A')
        reversedSeq = reversedSeq.replace('c', 'G')
        reversedSeq = reversedSeq.replace('g', 'C')
        return reversedSeq

def gmap_exon_finder(ncls, gmapLoc, model, coordIndex, processType):
        start = int(model[0][coordIndex].split('-')[0])
        stop = int(model[0][coordIndex].split('-')[1])
        overlaps = ncls.find_overlap(start, stop)                                                       # Although our ncls is 1-based, find_overlap acts as a range and is thus 0-based. We need to +1 to the stop to offset this.
        dictEntries = []
        for result in overlaps:
                if processType == 'boundary':
                        hitRange = range(result[0], result[1])
                        if start in hitRange or stop in hitRange:                                       # Remember that we made the ncls act 1-based for the purpose of indexing. However, when you retrieve results, you get a tuple like (start, stop, index).
                                dictEntries.append(gmapLoc[result[2]])                                  # In this tuple, stop is equal to the number you used to build the ncls. Thus, it actually equals stop + 1 now. It's a bit annoying, but it's how we need to do it to get true ranges for intersection.
                else:
                        hitRange = range(result[0], result[1])                                          # We already made the ncls 1-based, so we don't need to do anything here since range will naturally "remove" the last number.
                        if start in hitRange and stop in hitRange:
                                dictEntries.append(gmapLoc[result[2]])
        dictEntries = copy.deepcopy(dictEntries)                                                        # Any time we're deleting things from a section of a dictionary we need to build a deepcopy to keep the original dictionary intact.
        # Narrow down our dictEntries to hits on the same contig
        for k in range(len(dictEntries)-1, -1, -1):                                                     # Remember: gmapLoc = [[contigStart, contigStop, transStart, transStop, orient, geneID, identity, contigID]]. Loop through in reverse so we can delete entries without messing up the list.
                if dictEntries[k][7] != model[2]:
                        del dictEntries[k]
        # Provide extra information for boundary matches for reprioritisation later in the program
        if processType == 'boundary':
                encompass = []
                exact = []
                # Narrow down this list further to make sure we're only holding onto perfect boundary matches OR fully encompassing matches
                for n in range(len(dictEntries)-1, -1, -1):
                        if dictEntries[n][0] == start or dictEntries[n][1] == stop:
                                exact.append(dictEntries[n] + ['exact'])
                        elif dictEntries[n][0] <= start and dictEntries[n][1] >= stop:                  # This wasn't previously in the program, but it was noticed that boundary and internal handling was too different and causing problems. This should help to equalise their behaviour for the most part.
                                encompass.append(dictEntries[n] + ['encompass'])
                outList = exact + encompass                                                             # We can use this extra info of 'exact' or 'encompass' for prioritising exact matches, but allowing encompassing matches to be accepted if we find no exact matches/for memory across genes. This was needed to handle a few scenarios.
        else:
                outList = []
                for entry in dictEntries:
                        outList.append(entry + ['.'])                                                   # This is to accommodate the 'exact'/'encompass' behaviour handling. Boundary matches are really annoying to handle but the behaviour is necessary, unfortunately.
        # Sort and return list
        outList.sort(key = lambda x: (int(x[6]), x[1] - x[0]), reverse = True)                          # Provides a sorted list where, at the top, we have the longest and best matching hits.
        return outList

def check_future(ncls, gmapMatches, model, i):
        # Process if there is actually another downstream exon
        if i != len(model[0]) - 1:
                # Grab the GMAP matches for the future exon
                if i+1 == 0 or i+1 == len(model[0]) - 1:
                        futureMatches = gmap_exon_finder(ncls, gmapLoc, model, i+1, "boundary")
                else:
                        futureMatches = gmap_exon_finder(ncls, gmapLoc, model, i+1, "internal")
                        if futureMatches == []:                                                         # If we can't find an encompassing match for an internal exon, there's a good chance that PASA has artificially extended this exon to correct for an indel.
                                futureMatches = gmap_exon_finder(ncls, gmapLoc, model, i+1, "boundary")
                # Get futureMatches IDs
                futureIDs = []
                for match in futureMatches:
                        futureIDs.append(match[5])
                # Cull any current matches that don't show up in the future
                outMatches = []
                for match in gmapMatches:
                        if match[5] in futureIDs:
                                outMatches.append(match)
                return outMatches
        return gmapMatches

def gmap_curate(gmapMatches, model, coordIndex):
        coords = model[0][coordIndex].split('-')
        bestMatches = []
        # Detect perfect exon boundary matches
        for x in range(len(gmapMatches)-1,-1,-1):
                """By putting this check in we can ensure that, in the case that we have perfect boundary matches but with
                lower identity we can prioritise these above 'better' matches which don't respect exon boundaries"""
                if gmapMatches[x][0] == int(coords[0]) and gmapMatches[x][1] == int(coords[1]):
                        bestMatches.append(gmapMatches[x])
        # Return the best match if we can or, if not, just return any matches that meet our curation conditions
        if bestMatches == []:
                return gmapMatches
        else:
                """I added in this condition due to situation where we had two gmapMatches, one was "perfect" (100% identity, exact exon boundary alignment)
                but the other, despite 98% identity, still had a better SSW score simply because it was longer. I could normalise SSW score to handle this, but this is probably
                just as good and it should reduce the computational time of the script.
                Also important consideration: if there is a GMAP alignment which matches this exon's boundaries perfectly, it provides solid evidence that this exon boundary
                should not be changed, so we can limit our consideration to these matches."""
                return bestMatches

def patch_seq_extract(match, model):
        # Transcriptomic records
        transcriptRecord = copy.deepcopy(transRecords[match[5]])
        if match[4] != '+':                                                                             # Put it in the same orientation as the genome sequence [the genome sequence is always + orientation].
                transcriptRecord = transcriptRecord.reverse_complement()
        # Genomic patch (correlating to transcript alignment positions)
        genomePatchRec = genomeRecords[model[2]][int(match[0])-1:int(match[1])]
        return transcriptRecord, genomePatchRec

def ssw(genomePatchRec, transcriptRecord):
        # Perform SSW with parasail implementation
        profile = parasail.profile_create_16(str(genomePatchRec.seq), parasail.blosum62)
        alignment = parasail.sw_trace_striped_profile_16(profile, str(transcriptRecord.seq), 10, 1)
        genomeAlign = alignment.traceback.query
        transcriptAlign = alignment.traceback.ref
        # Figure out where we're starting in the genome with this alignment
        startIndex = str(genomePatchRec.seq).find(genomeAlign.replace('-', ''))
        # Figure out if we need downstream processing to identify an indel
        hyphen = 'n'
        if '-' in genomeAlign:
                hyphen = 'y'
        elif '-' in transcriptAlign:
                hyphen = 'y'
        return [transcriptAlign, genomeAlign, hyphen, startIndex, alignment.optimal_alignment_score]

def indel_location(transcriptAlign, genomeAlign, matchStart, model, startIndex, exonIndex, stop_codons):             # This function will check hyphens in the transcript (== deletions in the genome) and hyphens in the genome (== insertion from the transcript).
        # Check if this is likely to be worth bothering
        badChars = ['---', 'n']
        for char in badChars:                                                                           # This is a rough metric, but gap opens larger than three make us wonder whether this transcript does actually originate from the alignment position (maybe it's a paralogue?); in almost all cases, a real indel has a length of one.
                if char in transcriptAlign.lower() or char in genomeAlign.lower():
                        return 0, '.'                                                                   # We return 0 since that tells the main part of the script that this hit isn't good enough and to stick to the current model coordinates, the tmpVcf value doesn't matter in this case so just return '.'.
        # Process the alignment to find differences
        identical = 0
        gapCorrection = 0                                                                               # If we have a gap in our genomeAlign, every position after that needs to be minused 1 to correct for this. This value will hold onto this value and apply it to our index calculation.
        tmpVcf = {}                                                                                     # We want to add results into a temporary dictionary because, for sequences which mysteriously do not have good identity, we don't want to save their edit positions.
        literalVcf = {}                                                                                 # This is for the incorrect stop codon identification, we just grab the coordinates for the genome segment itself directly.
        for x in range(len(transcriptAlign)):
                genomeIndex = matchStart + startIndex + x - gapCorrection                               # This will correspond to the genomic contig index [note that we add startIndex to match[0] because we may have trimmed some of the 5' sequence during SW alignment].
                pair = transcriptAlign[x] + genomeAlign[x]                                              # Note that these are 1-based, so we'll need to account for this behaviour later.
                if pair[0] == pair[1]:
                        identical += 1
                elif pair[0] == '-':
                        if model[2] not in tmpVcf:
                                tmpVcf[model[2]] = {genomeIndex: ['.']}
                                literalVcf[x-gapCorrection] = ['.']
                        else:
                                if genomeIndex not in tmpVcf[model[2]]:
                                        tmpVcf[model[2]][genomeIndex] = ['.']
                                        literalVcf[x-gapCorrection] = ['.']
                                else:                                                                   # By doing this, if we have an insertion followed by a deletion, we'll handle both cases by just substituting the deleted base with the inserted one.
                                        tmpVcf[model[2]][genomeIndex][0] += '*'
                                        literalVcf[x-gapCorrection][0] += '*'
                elif pair[1] == '-':
                        if model[2] not in tmpVcf:
                                tmpVcf[model[2]] = {genomeIndex: [pair[0]]}
                                literalVcf[x-gapCorrection] = [pair[0]]
                        else:
                                if genomeIndex not in tmpVcf[model[2]]:
                                        tmpVcf[model[2]][genomeIndex] = [pair[0]]
                                        literalVcf[x-gapCorrection] = [pair[0]]
                                else:
                                        tmpVcf[model[2]][genomeIndex][0] += pair[0]
                                        literalVcf[x-gapCorrection][0] += pair[0]
                        gapCorrection += 1
        # Stop codon identification
        if stop_codons: # stop_codons refers to the user-specified argument; if it is true, this module was turned on, and so we'll call this function here
                tmpVcf = stop_codon_identify(literalVcf, genomeAlign, transcriptAlign, startIndex, exonIndex, matchStart, model, tmpVcf)
        # Calculate the (rough) identity score between the alignments
        pctIdentity = (identical / len(transcriptAlign)) * 100
        return pctIdentity, tmpVcf

def stop_codon_identify(literalVcf, genomeAlign, transcriptAlign, startIndex, exonIndex, matchStart, model, tmpVcf):
        if args.stop_codons and exonIndex != len(model[0]) - 1:                                         # Don't look at our final exon since it should have a stop codon in it. 
                # Edit out genomeAlign sequence to look like it will in the final genome output
                genomeAlign = genomeAlign.replace('-', '')
                outAlign = ''
                genomeIndices = []
                for x in range(len(genomeAlign)):
                        if x in literalVcf:
                                if literalVcf[x] == ['.']:
                                        continue
                                elif '*' in literalVcf[x][0]:
                                        outAlign += literalVcf[x][0][:-1]
                                        genomeIndices.append(x)
                                else:
                                        outAlign += literalVcf[x][0]
                                        outAlign += genomeAlign[x]
                                        genomeIndices += ['-']*len(literalVcf[x][0])
                                        genomeIndices.append(x)
                        else:
                                genomeIndices.append(x)
                                outAlign += genomeAlign[x]
                # Rescue the last index if there's weird stuff going on in the end of the alignment (this should never happen naturally, but if genomeAlign ends with hyphens then we wouldn't end up adding the last bit from the literalVcf)
                if x+1 in literalVcf:
                        if literalVcf[x+1] == ['.']:
                                doNothing = 1
                        elif '*' in literalVcf[x+1][0]:
                                outAlign += literalVcf[x+1][0][:-1]
                                genomeIndices.append(x+1)
                        else:
                                outAlign += literalVcf[x+1][0]
                                genomeIndices += ['-']*len(literalVcf[x+1][0])
                # Check to see if there are any frames which lack a stop codon
                stopCodonsPos = ['tag', 'taa', 'tga']
                stopCodonsNeg = ['cta', 'tta', 'tca']
                if model[1] == '-':
                        stopCodons = stopCodonsNeg                                      # Just pick out our stop codons we're looking for now so we don't need to check what orientation the genomeAlign sequence is in future lines.
                else:
                        stopCodons = stopCodonsPos                                      # We want to treat everything in the + orientation since it's much easier to just stay in this mode versus switching between and recalculating the proper genomic index.
                firstStop = [-1,-1,-1]
                for frame in range(3):
                        frameNuc = outAlign[frame:].lower()
                        codons = [frameNuc[i:i+3] for i in range(0, len(frameNuc), 3)]
                        for c in range(len(codons)):
                                if codons[c] in stopCodons:
                                        firstStop[frame] = c
                                        break
                # If we lack any frames with stop codons, compare these stop codon locations to the transcript to see if the transcript lacks a stop codon in this position
                if -1 not in firstStop:
                        # Compare codons
                        transcriptAlign = transcriptAlign.replace('-', '')
                        for frame in range(3):
                                codonIndex = firstStop[frame]
                                genomeCodon = outAlign[frame:][codonIndex*3:(codonIndex*3) + 3]
                                codonIndices = genomeIndices[frame:][codonIndex*3:(codonIndex*3) + 3]
                                transcriptCodon = transcriptAlign[frame:][codonIndex*3:(codonIndex*3) + 3]
                                if genomeCodon != transcriptCodon:
                                        # Compare each base that is part of this codon for differences
                                        for x in range(len(genomeCodon)):
                                                if codonIndices[x] == '-':                                                      # If this base was inserted into the genome sequence then we don't consider it.
                                                        continue
                                                genomeIndex = matchStart + startIndex + codonIndices[x]                         # Because we tracked the indices of each base in the genomeAlign sequence, we can directly translate this back to the overall genomic coordinates.
                                                # Compare the codons and make any modifications suggested
                                                if genomeCodon[x] != transcriptCodon[x]:
                                                        if model[2] not in tmpVcf:
                                                                tmpVcf[model[2]] = {genomeIndex: [transcriptCodon[x] + '*']}    # We'll use asterisks to mark substitutions.
                                                        else:
                                                                tmpVcf[model[2]][genomeIndex] = [transcriptCodon[x] + '*']
        return tmpVcf

def prot_identity(prot1, prot2):
        identical = 0
        for x in range(len(prot1)):
                if prot1[x] == prot2[x]:
                        identical += 1
        # Calculate the (rough) identity score between the alignments
        pctIdentity = (identical / len(prot1)) * 100
        return pctIdentity

def vcf_edit(tmpVcf, contigID, coordRange):
        genomeSeq = str(genomeRecords[contigID].seq)[min(coordRange)-1:max(coordRange)]                         # -1 for 0-based.
        if contigID in tmpVcf:                                                                                  # We used to always have edits when using this function but because of the stop codon check sometimes that isn't true. We'll just directly return the genome sequence if there are no edits.
                # Extract edit positions
                subVcfDict = tmpVcf[contigID]
                tmpVcfList = []
                for key, value in subVcfDict.items():
                        if key in coordRange:
                                tmpVcfList.append([key, value[0]])
                tmpVcfList.sort(reverse=True)
                # Edit the genome sequence
                for pair in tmpVcfList:
                        indelIndex = pair[0] - min(coordRange)                                                  # pair[0] refers to the actual genomic index, but we want to find the location in this particular genome section, so we just minus the start coordinate.
                        if pair[1] == '.':
                                genomeSeq = genomeSeq[:indelIndex] + genomeSeq[indelIndex+1:]                   # Since pair[0] and coordRange are 1-based, minusing these results in an index that is, essentially, 0-based.
                        elif '*' in pair[1]:
                                genomeSeq = genomeSeq[:indelIndex] + pair[1][:-1] + genomeSeq[indelIndex+1:]      # This makes a substitution in our genome (as marked by the asterisk).
                        else:
                                genomeSeq = genomeSeq[:indelIndex] + pair[1] + genomeSeq[indelIndex:]           # Because of this, we +1 to the second bit to skip the indelIndex, and leave this neutral to simply insert a base at the indel index.
        return genomeSeq    

def cds_build(origCoords, newCoords, contigID, orientation, tmpVcf):
        # Build the original gene model
        origCDS = []
        for coord in origCoords:
                splitCoord = coord.split('-')
                cdsBit = str(genomeRecords[contigID].seq)[int(splitCoord[0])-1:int(splitCoord[1])]
                if orientation == '-':
                        cdsBit = reverse_comp(cdsBit)
                origCDS.append(cdsBit)
        # Resolve overlaps in the new gene model coords
        finalCDSCoords = []
        for i in range(len(newCoords)-1,-1,-1):
                if i != 0:
                        # Remove redundant coords immediately
                        if newCoords[i] == newCoords[i-1]:
                                del newCoords[i]                                                        # We get redundant coords when we naturally join two exons.
                                continue
                        # Compare current to future coord and look for overlap
                        split1 = newCoords[i].split('-')
                        split2 = newCoords[i-1].split('-')
                        if int(split1[1]) >= int(split2[0]) and int(split1[0]) <= int(split2[1]):       # If the end of seq1 > start of seq2 AND the start of seq1 < end of seq2, then we know there's overlap regardless of orientation.
                                # Pick out the longest coordinate
                                seq1Len = int(split1[1]) - int(split1[0])
                                seq2Len = int(split2[1]) - int(split2[0])
                                if seq2Len > seq1Len:
                                        del newCoords[i]
                                else:
                                        del newCoords[i-1]
                        else:
                                finalCDSCoords.append(newCoords[i])
                else:
                        finalCDSCoords.append(newCoords[i])
        finalCDSCoords.reverse()                                                                        # This is because we loop through newCoords in reverse - we need it to go back to its original order for the CDS orientation to be mainted.
        # Build the new gene model
        newCDS = []
        for coord in finalCDSCoords:
                splitCoord = coord.split('-')
                coordRange = range(int(splitCoord[0]), int(splitCoord[1])+1)                            # Our VCF dictionary is 1-based at this point, so we want our range to act like this, too.
                cdsBit = vcf_edit(tmpVcf, contigID, coordRange)
                if orientation == '-':
                        cdsBit = reverse_comp(cdsBit)
                newCDS.append(cdsBit)
        # Joing our CDS bits together
        origCDS = ''.join(origCDS)
        newCDS = ''.join(newCDS)
        return origCDS, newCDS   

def vcf_merge(vcf1, vcf2):                                                                              # This will merge vcf2 into vcf1 (currently this direction doesn't matter, but I might change this later).
        # Make deep copies
        vcf1 = copy.deepcopy(vcf1)
        vcf2 = copy.deepcopy(vcf2)
        # Merge the temporary vcf into the main one
        delKeys = []                                                                                    # In the rare scenario we have a serious disagreement (one transcript says insert, another says delete) we'll hold onto these keys and delete them from the dictionary [this will likely never happen, but just in case].
        delContigs = []                                                                                 # We just want to hold onto the corresponding contig ID for any index keys we delete.
        for key, value in vcf2.items():
                if key not in vcf1:
                        vcf1[key] = value
                else:
                        value2 = vcf1[key]
                        for k2, v2 in value.items():
                                if k2 in value2:
                                        if value2[k2] == v2:
                                                continue                                                # We don't care if it's identical, we just want to find situations that don't agree.
                                        elif value2[k2] != ['.'] and v2 != ['.']:                       # I haven't seen this scenario occur, but if both transcripts agree that an insertion should occur here then we'll just stick with what we found first.
                                                continue
                                        else:
                                                delKeys.append(k2)
                                                delContigs.append(key)
                                else:
                                        value2[k2] = v2
        # Delete any conflicting indel locations
        for i in range(len(delKeys)):
                del vcf1[delContigs[i]][delKeys[i]]
                # Clean up any empty dictionary keys
                if vcf1[delContigs[i]] == {}:
                        del vcf1[delContigs[i]]
        return vcf1

def vcf_output(outFileName, vcf, comment):
        if not os.path.isfile(outFileName):
                with open(outFileName, 'w') as fileOut:
                        fileOut.write('#contig_id\tposition\treplacement\n')
        with open(outFileName, 'a') as fileOut:
                if comment != '.':
                        fileOut.write(comment + '\n')
                for key, value in vcf.items():
                        value = list(value.items())
                        value.sort()
                        for pair in value:
                                fileOut.write('\t'.join([key, str(pair[0]), str(pair[1][0])]) + '\n')

def double_polish_check(modVcf, mainVcf, vcf, rangeDict, model, result):
        # Deep copy everything to stop annoying problems...
        modVcf = copy.deepcopy(modVcf)                                                                          # Need deepcopies so we can delete these freely in our tmpVcf without changing the original modelVcf/vcfDict dictionaries.
        mainVcf = copy.deepcopy(mainVcf)
        # Check for double polishing
        currentRange = set(range(int(result[5]), int(result[6]) + 1))
        ovlRanges = set()
        if model[2] in rangeDict:                                                                               # This function looks through out rangeDict (i.e., our previously polished positions) and looks at our current range that we intend to polish to see if there is overlap.
                for value in rangeDict[model[2]]:                                                               # If there is overlap, we want to record the exact positions where overlap occurs.
                        ovl = currentRange.intersection(value)                                                  # This lets us subsequently remove these positions from our input vcf (which will be the tmpVcf) to prevent double polishing.
                        if ovl != set():
                                ovlRanges = ovlRanges.union(ovl)
        # Remove any vcf positions covered by ovlRanges
        delPositions = []
        for vcfPos in vcf[model[2]]:
                if vcfPos in ovlRanges:
                        delPositions.append(vcfPos)
        for entry in delPositions:
                del vcf[model[2]][entry]
        # Add in any vcf positions covered by ovlRanges from our overall vcfDict
        if model[2] in modVcf:
                for key, value in modVcf[model[2]].items():
                        if key in currentRange:
                                vcf[model[2]][key] = value
        if mainVcf != '.':                                                                                      # This lets us re-use this function in the gene rescue module.
                if model[2] in mainVcf:
                        for key, value in mainVcf[model[2]].items():
                                if key in currentRange:
                                        vcf[model[2]][key] = value
        # Add to our rangeDict if we ended up with any remaining edits
        if vcf != {}:
                if model[2] not in rangeDict:
                        rangeDict[model[2]] = [range(int(result[5]), int(result[6]) + 1)]                       # Make the range 1-based.
                else:
                        if range(int(result[5]), int(result[6]) + 1) not in rangeDict[model[2]]:
                                rangeDict[model[2]].append(range(int(result[5]), int(result[6]) + 1))
        return vcf, rangeDict

def translate_cds(seq1, seq2):
        # Translate into ORFs and grab the longest bits inbetween stop codons
        records = [Seq(seq1, generic_dna), Seq(seq2, generic_dna)]
        longest = ['','']
        for i in range(len(records)):
                tmpLongest = ''
                for frame in range(3):
                        with warnings.catch_warnings():
                                warnings.simplefilter('ignore')                                                 # This is just to get rid of BioPython warnings about len(seq) not being a multiple of three. We know that in two of these frames that will be true so it's not a problem.
                                frameProt = str(records[i][frame:].translate(table=1))
                        frameProt = frameProt.split('*')
                        frameProt.sort(key = len, reverse = True)
                        frameOrf = frameProt[0]
                        if len(frameOrf) > len(tmpLongest):
                                tmpLongest = frameOrf
                longest[i] = tmpLongest
        return longest

def geneblocks_update(geneBlocksDict, model, modelCoords):
        orientation = model[1]
        contigID = model[2]
        # Derive coordinates
        if orientation == '+':
                start = int(modelCoords[0].split('-')[0])
                stop = int(modelCoords[-1].split('-')[1])
        else:
                start = int(modelCoords[-1].split('-')[0])
                stop = int(modelCoords[0].split('-')[1])
        # Update geneBlocksDict
        if contigID not in geneBlocksDict:
                geneBlocksDict[contigID] = [[start, stop, model[3], orientation]]
        else:
                geneBlocksDict[contigID].append([start, stop, model[3], orientation])
        return geneBlocksDict

def gene_overlap_validation(geneBlocks):
        outlist = []
        for key, value in geneBlocks.items():
                value.sort()
                for i in range(len(value)-1):
                        if value[i][1] >= value[i+1][0]:
                                basename1 = isoRegex.search(value[i][2]).group(1)
                                basename2 = isoRegex.search(value[i+1][2]).group(1)
                                if basename1 != basename2 and value[i][3] == value[i+1][3]:             # i.e., if these aren't isoforms (basenames will be identical if they are) and the end of gene1 overlaps the start of gene 2.
                                        outlist.append(value[i][2] + '\t' + value[i+1][2])
        return outlist

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
                        if not sl[8].startswith('ID=blat.proc'):                                        # May as well make the code compatible with BLAT.
                                geneID = sl[8].split(';')[1][5:]                                        # Start from 5 to get rid of 'Name='.
                        else:
                                geneID = sl[8].split(';')[1].split()[0][7:]                             # Start from 7 to get rid of 'Target='.
                        contigStart = int(sl[3])
                        contigStop = int(sl[4])
                        identity = float(sl[5])
                        if identity < cutoff:                                                           # Speed up program by only holding onto hits that will pass our cutoff check.
                                continue
                        orient = sl[6]
                        if not sl[8].startswith('ID=blat.proc'):
                                transStart = int(sl[8].split(';')[2].split()[1])
                                transStop = int(sl[8].split(';')[2].split()[2])
                        else:
                                transStart = int(sl[8].split(';')[1].split()[1])
                                transStop = int(sl[8].split(';')[1].split()[2])
                        # Add to our NCLS                                                               # We index using ranges since it provides an easy way to retrieve GMAP matches by coordinates. Since these coordinates aren't unique, we filter any results returned by their contig ID.
                        starts.append(contigStart)
                        ends.append(contigStop+1)                                                       # NCLS indexes 0-based, so +1 to make this more logically compliant with gff3 1-based system.
                        ids.append(ongoingCount)
                        gmapLoc[ongoingCount] = [contigStart, contigStop, transStart, transStop, orient, geneID, identity, contigID]
                        ongoingCount += 1
        # Build the NCLS object
        starts = pd.Series(starts)
        ends = pd.Series(ends)
        ids = pd.Series(ids)
        ncls = NCLS(starts.values, ends.values, ids.values)
        return ncls, gmapLoc

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
        for key, value in gffCoordDict.items():
                for mrna in value:                                                                      # This just reorganises the results a bit. Only reason I'm doing this is because I'm reusing bits of old code, and it's easier to just take a small efficiency hit here to not bother rejigging the whole thing leading up to here.
                        nuclDict[mrna[0]] = [mrna[1], mrna[3], mrna[2], mrna[0]]
        return nuclDict

def validate_args(args):
        # Validate input file locations
        if not os.path.isfile(args.gff3File):
                print('I am unable to locate the input gff3 gene annotation file (' + args.gff3File + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        elif not os.path.isfile(args.genomeFile):
                print('I am unable to locate the input genome fasta file (' + args.genomeFile + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        elif not os.path.isfile(args.gmapFile):
                print('I am unable to locate the input GMAP transcript alignment gff3 file (' + args.gmapFile + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        elif not os.path.isfile(args.transcriptomeFile):
                print('I am unable to locate the input transcriptome fasta file (' + args.transcriptomeFile + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        # Handle file overwrites
        tmpFileName = None
        if os.path.isfile(args.outputFileName):
                if args.force:
                        # Temporarily move the file to the current directory and delete the file at the end of program run - this acts as a safety mechanism if someone actually ends up not wanting to overwrite the output file.
                        tmpFileName = file_name_gen('deGRIT_backup', '_' + os.path.basename(args.outputFileName))
                        shutil.move(args.outputFileName, tmpFileName)                                   # I'm going to do this before I alert the user since they might immediately cause a KeyboardInterrupt and I don't know what happens if you do this during shutil.move().
                        print('You\'ve specified that you want to overwrite ' + args.outputFileName)
                        print('Is that right? I\'m going to rename this file to "' + tmpFileName + '" and hold onto the file in the current directory until this program exits.')
                        print('If you don\'t want to delete this file, kill this process and you can retrieve the file.')
                else:
                        print(args.outputFileName + ' already exists. Either provide the -fo argument to this program or delete/move/rename this file and run the program again.')
                        quit()
        return tmpFileName

## Functions for optional arguments
def file_name_gen(prefix, suffix):
        ongoingCount = 2
        while True:
                if not os.path.isfile(prefix + '1' + suffix):
                        return prefix + '1' + suffix
                elif os.path.isfile(prefix + str(ongoingCount) + suffix):
                        ongoingCount += 1
                else:
                        return prefix + str(ongoingCount) + suffix

def trans_pos(sswResult):                                                                               # Function designed for assisting the logging process by providing transcript start - stop positions.
        transcriptRecord = copy.deepcopy(transRecords[sswResult[7]])
        if sswResult[8] != '+':                                                                         # Put it in the same orientation as the genome sequence [the genome sequence is always + orientation].
                transcriptRecord = transcriptRecord.reverse_complement()
        transcriptBit = sswResult[0].replace('-', '')                                                   # remove hyphens since these aren't in the original record.
        startpos = str(transcriptRecord.seq).find(transcriptBit) + 1                                    # +1 to make it 1-based.
        endpos = str(startpos + len(transcriptBit) - 1)                                                 # -1 since our start pos is already 1-based.
        transPos = str(startpos) + '-' + str(endpos)
        return transPos

def log_update(args, logName, inputList):
        if args.log:
                if not os.path.isfile(logName):
                        with open(logName, 'w') as logFile:
                                logFile.write('#contig_id\tgene_name\toriginal_exon_coords\tbest_transcript_match\ttranscript_coords_against_genome\taligned_region_of_transcript\tmodified_locations\n')
                with open(logName, 'a') as logFile:
                        # Pull out GMAP match names and coordinates
                        if inputList[3] == []:
                                matchName = '_'
                                matchCoord = '_'
                        else:
                                matchName = inputList[3][0][7]
                                matchCoord = str(inputList[3][0][5]) + '-' + str(inputList[3][0][6])
                        # Pull out the transcript alignment positions if relevant
                        transPos = '_'
                        if inputList[3] != [] and inputList[5] != '.':                                  # i.e., if inputList[3] is a sswResult list and not '.', and if inputList[5] is 'hit' and not '.' [we just specify 'hit' for the purpose of this statement].
                                transPos = trans_pos(inputList[3][0])
                        # Format the edit positions for this exon if relevant
                        editPos = '_'
                        if inputList[4] != '.':
                                editPos = vcf_line_format(inputList[4])
                        # Write to log file
                        logFile.write('\t'.join([inputList[1][2], inputList[1][3], inputList[1][0][inputList[2]], matchName, matchCoord, transPos, editPos]) + '\n')

def vcf_line_format(inputVcf):
        editPos = ''
        for key, value in inputVcf.items():
                for k2, v2 in value.items():
                        editPos += str(k2) + ':' + v2[0] + ','
        editPos = editPos[:-1]                                                                          # Remove the last comma.
        if editPos == '':
                editPos = '_'                                                                           # This can happen if, after double polish checking, our tmpVcf is empty.
        return editPos

def log_comment(args, logName, text):
        if args.log:
                if not os.path.isfile(logName):
                        with open(logName, 'w') as logFile:
                                logFile.write('#contig_id\tgene_name\toriginal_exon_coords\tbest_transcript_match\ttranscript_coords_against_genome\taligned_region_of_transcript\tmodified_locations\n')
                with open(logName, 'a') as logFile:
                        logFile.write(text + '\n')

def verbose_print(args, text):
        if args.verbose:
                print(text)

## Gene model rescue functions
def gmap_parse_models(args, cutoff, transRecords):
        if args.rescue_genes == False:
                return {}
        # Re-index our nuclDict into a format capable of comparison to our gmapLoc dictionary
        nuclNcls, nuclModels = reindex_nucldict(nuclDict)
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
        # Curate our potential gene models through various checks
        validModels = {}
        for key, value in gmapModels.items():
                ## Check 1: Was this marked as 'not_good' earlier? (this means the gene model contained at least one exon which didn't meet out cutoff)
                if value == 'not_good':
                        continue
                ## Check 2: Does this model have more than 1 exon? (single exon genes may be pseudogenes)
                if len(value[0]) < 2:
                        continue
                ## Check 3: Does the full CDS align to the genome? (fragmentary alignments mean this CDS doesn't belong to this genomic region or that the genome sequence itself is misassembled; in both cases, we can't annotate this gene model correctly)
                # Get the transcript details
                record = transRecords[value[3]]
                seq = str(record.seq)
                seqSet = set(range(1, len(seq)+1))
                # Format the alignment details
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
                check3pass = 'y'
                for x in range(len(gaps)-1):
                        if gaps[x] + 1 == gaps[x+1]:
                                gapStretch += 1
                        else:
                                gapStretch = 1
                        if gapStretch >= 5:                                                             # A 5bp gap also isn't much, but we don't want to make a "putative gene model" which doesn't match the transcript nearly perfectly.
                                check3pass = 'n'
                                break
                if check3pass == 'n':
                        continue
                ## Check 4: Can we find any exons that do not match against the current gene annotation? (we don't want to re-process exons that are incorporated into the annotation, we just want to look at any exons that currently lack annotation)
                for x in range(len(value[0])):
                        # Query the ncls object for overlaps
                        coordSplit = value[0][x].split('-')
                        overlaps = nuclNcls.find_overlap(int(coordSplit[0]), int(coordSplit[1])+1)
                        nuclHits = []
                        for result in overlaps:
                                nuclHits.append(nuclModels[result[2]])
                        nuclHits = copy.deepcopy(nuclHits)
                        # Narrow down our overlaps to hits on the same contig
                        for k in range(len(nuclHits)-1, -1, -1):
                                if nuclHits[k][0] != value[2]:
                                        del nuclHits[k]
                        # Save any exons that lack overlaps
                        if nuclHits != []:                                                              # This list won't be empty if it overlaps an established gene model.
                                continue
                        else:
                                validModels[key] = value[0:3] + [key]                                   # This means we've found a transcript which contains an exon not covered by the annotation, so it's worth us inspecting it.
                                break
        return validModels

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

# Build regex for later use
isoRegex = re.compile(r'(evm\.model\.utg\d{1,10}(_pilon_pilon)?\.\d{1,10})')

### USER INPUT
usage = """%(prog)s aims to improve the ability to reannotate gene models. In order to work, this program
requires a gff3 file of gene annotations alongside its respective genome fasta file in addition to a gff3 file
of transcript alignments with its respective transcriptome fasta file. These files will be used to compare the
genome sequence to the aligned transcript sequence to identify any occurrences of indel errors. By correcting these
indels, reannotation of gene models can take place which will provide more accurate results.
Note: This program is designed to work with CDS regions from transcripts; this reduces the chance of falsely 
interrupting a reading frame with an edit. You can predict the CDS region using TransDecoder or EvidentialGene.
"""
# Reqs
p = argparse.ArgumentParser(description=usage)

p.add_argument("-an", "--annotation", dest="gff3File",
               help="Input gff3 gene annotation file name")
p.add_argument("-gen", "--genomefile", dest="genomeFile",
               help="Input genome contig fasta file name")
p.add_argument("-gm", "--gmap", dest="gmapFile",
               help="Input gff3 gmap transcript alignment file name")
p.add_argument("-tr", "--transfile", dest="transcriptomeFile",
               help="Input nucleotide transcriptome fasta file name (this is the same transcript file used for GMAP alignment)")
p.add_argument("-o", "--output", dest="outputFileName",
               help="Output results file name")
# Opts
p.add_argument('-r', '--rescue_genes', dest="rescue_genes", action='store_true',
               help="Optionally perform extended gene model rescue module (this is recommended)", default=False)
p.add_argument('-s', '--stop_codons', dest="stop_codons", action='store_true',
               help="Optionally identify in-frame stop codons and correct these if indicated by transcript evidence.", default=False)
p.add_argument("-fo", "--force", dest="force", action='store_true',
               help="Optionally allow the program overwrite existing files at your own risk", default=False)
p.add_argument('-v', '--verbose', dest="verbose", action='store_true',
               help="Print program details to terminal", default=False)
p.add_argument('-l', '--log', dest="log", action='store_true',
               help="Additionally produce a detailed logging file as output", default=False)

args = p.parse_args()

# Validate arguments and get log file name
tmpFileName = validate_args(args)
logName = file_name_gen('deGRIT_' + os.path.basename(args.genomeFile).rsplit('.', maxsplit=1)[0] + '_run', '.log')

# Declare values needed for processing
minCutoff = 97                                                                                          # I don't think this value should be modifiable - the program is built around this value, increasing it will result in finding very few results, and decreasing it will likely result in false changes.

# Load genome file as a dictionary
genomeRecords = SeqIO.to_dict(SeqIO.parse(open(args.genomeFile, 'r'), 'fasta'))
verbose_print(args, 'Loaded genome fasta file')

# Parse the gff3 file
nuclDict = cdna_parser(args.gff3File)
verbose_print(args, 'Parsed the annotations gff3 file')

# Parse the gmap alignment file for transcript alignment locations
ncls, gmapLoc = gmap_parse_ncls(args.gmapFile, minCutoff)
verbose_print(args, 'Parsed GMAP gff3 file')

# Parse the transcriptome file
transRecords = SeqIO.to_dict(SeqIO.parse(open(args.transcriptomeFile, 'r'), 'fasta'))
verbose_print(args, 'Loaded transcriptome fasta file')

# Declare values needed to retain results across loops and into the gene rescue module
vcfDict = {}                                                                                            # This dictionary will hold onto values in a style that is similar to VCF, making output and parsing easier.
geneBlocks = {}                                                                                         # This dictionary serves as a form of validation. By holding onto model starts/stops, we'll be able to check for overlap which will tell us if gene models will end up merged in the reannotation.
polishedRanges = {}
prevBestResult = ''                                                                                     # This will let us remember our previous best match. This will reduce the chance of small differences of opinion ruining gene models.
                                                                                                        # It's worth noting that we also want to hold this value across genes, hence why I want to declare a blank value before we start our core loop.
### CORE PROCESS
if args.rescue_genes:
        gmapRescueDict = gmap_parse_models(args, minCutoff, transRecords)
        verbose_print(args, 'Loaded GMAP gff3 file for gene rescue')
else:
        gmapRescueDict = {}     # Set up an empty dictionary so iteration stops during the below loop after processing nuclDict
for z in range(len([nuclDict, gmapRescueDict])):
        if z == 0:
                verbose_print(args, '### Main gene improvement module start ###')
        elif z == 1 and args.rescue_genes:
                verbose_print(args, '### Gene rescue module start ###')
        for key, model in [nuclDict, gmapRescueDict][z].items():
                # Hold onto both the original gene model, as well as the new gene model resulting from indel correction/exon boundary modification
                origModelCoords = []
                newModelCoords = []
                modelVcf = {}                                                                                   # This will hold onto the VCF-like dictionary for this model; we'll incorporate it into the main one if we accept these modifications.
                # Scan through each individual model's exons
                for i in range(len(model[0])):
                        # Find GMAP matches that align over the exon region of this coordinate
                        if i == 0 or i == len(model[0]) - 1:
                                gmapMatches = gmap_exon_finder(ncls, gmapLoc, model, i, "boundary")
                        else:
                                gmapMatches = gmap_exon_finder(ncls, gmapLoc, model, i, "internal")
                                if gmapMatches == []:                                                           # If we can't find an encompassing match for an internal exon, there's a good chance that PASA has artificially extended this exon to correct for an indel, in which case boundary matching should be more appropriate.
                                        gmapMatches = gmap_exon_finder(ncls, gmapLoc, model, i, "boundary")
                        if gmapMatches == []:
                                origModelCoords.append(model[0][i])                                             # If there is no transcript support for this exon, it might be a spurious attempt by PASA/EVM to keep the gene inframe. Thus, we'll only save these coords under the origModel.
                                # Log
                                log_update(args, logName, [key, model, i, gmapMatches, '.', '.'])
                                continue
                        gmapMatches = gmap_curate(gmapMatches, model, i)
                        # Look into the future: cull any GMAP matches if they don't show up as a match in the next exon too
                        gmapMatches = check_future(ncls, gmapMatches, model, i)
                        # Continue if no GMAP matches
                        if gmapMatches == []:
                                origModelCoords.append(model[0][i])
                                newModelCoords.append(model[0][i])                                              # In this case, there IS transcript support for this exon, but it's not good enough for us to make edits with. Thus, we'll just hold onto the original coordinates for our new model.
                                # Log
                                log_update(args, logName, [key, model, i, gmapMatches, '.', '.'])
                                continue
                        # Find the best GMAP match by SSW alignment score
                        sswResults = []
                        for match in gmapMatches:
                                # Grab the sequences for alignment                                              # Note that we're going to compare the portion of the genome which the transcript hits (from GMAP) to the full transcript since GMAP handles N's weirdly and thus its transcript coordinates cannot be used.
                                transcriptRecord, genomePatchRec = patch_seq_extract(match, model)
                                # Perform SSW alignment
                                sswResults.append(ssw(genomePatchRec, transcriptRecord) + [match[0], match[1], match[5], match[4], match[8]])  # SSW returns [transcriptAlign, genomeAlign, hyphen, startIndex, alignment.optimal_alignment_score), and we also + [matchStart, matchEnd, matchName, matchOrientation, 'exact/encompass'] to this.
                        sswResults.sort(key = lambda x: (-x[4], x[2], x[3]))                                    # i.e., sort so score is maximised, then sort by presence of hyphens then by the startIndex.
                        # Change our cutoff to reflect different levels of support (only 1 alignment is more suspect, multiple alignments means the best hit has a high chance of actually originating here)
                        if len(sswResults) == 1 and z == 0:                                                     # We can give the gene model rescue a bit more leeway since we can simply reject any changes that reduce ORF length - for the main module, ORF length reductions can be a good thing, so we need to maintain strictness.
                                tmpCutoff = 98
                        else:
                                tmpCutoff = minCutoff                                                           # Technically I could just declare the minCutoff at this part of the program and not before this loop begins, but I think it helps to understand the decision making process of this program.
                        # Loop through our sswResults and find any good alignments
                        acceptedResults = []
                        for j in range(len(sswResults)):
                                sswIdentity, tmpVcf = indel_location(sswResults[j][0], sswResults[j][1], sswResults[j][5], model, sswResults[j][3], i, args.stop_codons)     # This will update our vcfDict with indel locations.
                                if sswIdentity >= tmpCutoff:
                                        acceptedResults.append([sswIdentity, tmpVcf, sswResults[j]])
                        # End checking if no results are trustworthy
                        if acceptedResults == []:
                                origModelCoords.append(model[0][i])
                                newModelCoords.append(model[0][i])                                              # Like above after gmap_curate, there is transcript support for this exon. Here, we chose not to make any changes, so we'll stick to the original coordinates.
                                # Log
                                log_update(args, logName, [key, model, i, sswResults, '.', 'hit'])
                        else:
                                # Check to see if one of these results is identical to the accepted one from the last exon [This holds memory across genes]
                                sswIdentity, tmpVcf, result = acceptedResults[0]                                # Hold onto the best result and overwrite these values if we find the result used for the previous exon/an exact boundary match if relevant.
                                for l in range(len(acceptedResults)):
                                        if i == 0 or i == len(model[0]) - 1:
                                                if result[9] != 'exact' and acceptedResults[l][2][9] == 'exact':
                                                        sswIdentity, tmpVcf, result = acceptedResults[l]        # This allows us to reprioritise exact matches over encompassing matches if we're looking at a gene boundary.
                                        if acceptedResults[l][2][7] == prevBestResult:
                                                if acceptedResults[l][2][2] == 'n' or result[2] != 'n':
                                                        sswIdentity, tmpVcf, result = acceptedResults[l]        # If our prevBestMatch has no hyphen, choose this. If our highest scoring result has no hyphen and our prevBestMatch does, then stick with the best result to minimise potentially false changes.
                                                        break
                                prevBestResult = result[7]
                                # Does this result have any edits?
                                if tmpVcf == {}:
                                        origModelCoords.append(model[0][i])
                                        newModelCoords.append(str(result[5]) + '-' + str(result[6]))            # Like above after gmap_curate, there is transcript support for this exon. Here, there were simply no changes to make, so we'll actually use the newer coordinates.
                                        # Log
                                        log_update(args, logName, [key, model, i, [result], '.', 'hit'])
                                        continue
                                # Double polishing check [This also holds memory across genes]
                                tmpVcf, polishedRanges = double_polish_check(modelVcf, vcfDict, tmpVcf, polishedRanges, model, result)
                                # Modify our modelVcf
                                modelVcf = vcf_merge(modelVcf, tmpVcf)
                                origModelCoords.append(model[0][i])
                                newModelCoords.append(str(result[5]) + '-' + str(result[6]))
                                # Log
                                log_update(args, logName, [key, model, i, [result], tmpVcf, 'hit'])             # Put result into a list since our logging function expects a list of lists, for which it retrieves the first entry.
                                
                # Hold onto any indel positions and provide logging information about this
                if modelVcf == {}:
                        if z == 0:
                                geneBlocks = geneblocks_update(geneBlocks, model, origModelCoords)
                        # Verbose and log
                        verbose_print(args, 'Found no edits [' + model[3] + ']')
                        log_comment(args, logName, '#' + model[3] + '\tNo edits found')
                else:
                        if z == 1:                                                                              # i.e., if this is gene rescue module
                                # For the gene rescue module we don't have a true "original" gene model, and the transcript alignments don't tell us the true orientation w/r/t to the genome, so we need to figure this out for ourselves
                                origCDS, newCDS = cds_build(origModelCoords, origModelCoords, model[2], '+', modelVcf)
                                origRevCDS, newRevCDS = cds_build(origModelCoords, origModelCoords, model[2], '-', modelVcf)
                                origProt, newProt = translate_cds(origCDS, newCDS)
                                origRevProt, newRevProt = translate_cds(origRevCDS, newRevCDS)
                                if len(newRevProt) >= len(newProt) or len(origRevProt) >= len(origProt):        # This should work for almost all scenarios. In the rare scenario that it doesn't, all that happens is that we report logging output that looks a bit strange - the user can validate any indel positions themselves.
                                        model[1] = '-'
                                else:
                                        model[1] = '+'
                        origCDS, newCDS = cds_build(origModelCoords, newModelCoords, model[2], model[1], modelVcf)
                        origProt, newProt = translate_cds(origCDS, newCDS)
                        sizeDiff = abs(round((len(origProt) - len(newProt)) / len(newProt) * 100, 3))
                        # Is the newCDS at least as long as the original CDS without internal stop codons?
                        if len(newProt) > len(origProt):
                                vcfDict = vcf_merge(vcfDict, modelVcf)
                                if z == 0:
                                        geneBlocks = geneblocks_update(geneBlocks, model, newModelCoords)
                                # Verbose and log
                                verbose_print(args, 'Looks like we improved this model! [' + model[3] + ']')
                                log_comment(args, logName, '#' + model[3] + '\tModel length increased\tOriginal model is ' + str(sizeDiff) + '% shorter than new\tOld=' + origProt + '\tNew=' + newProt + '\tAccepted edits=' + vcf_line_format(modelVcf))
                        elif len(newProt) == len(origProt):
                                # How much did we change this model?
                                protAlign = ssw(SeqRecord(Seq(origProt, generic_protein)), SeqRecord(Seq(newProt, generic_protein)))
                                protIdentity = prot_identity(protAlign[0], protAlign[1])
                                lengthDiff = abs(round((len(origProt) - len(protAlign[0])) / len(protAlign[0]) * 100, 3))       # This lets us know how much of the sequence end is being chopped off. A high identity alignment might still occur in the 5' region but after any modifications it might stop aligning.
                                # Don't save changes if it was just a variant
                                if protIdentity >= 98.0 and lengthDiff <= 10.0:                                                 # This choice is arbitrary. Small changes in the sequence are likely to be variants from another member of the same species, large changes might indicate a correction of a frameshift resulting in the same eventual stop position.
                                        # Verbose and log
                                        verbose_print(args, 'Found no edits [' + model[3] + ']')
                                        log_comment(args, logName, '#' + model[3] + '\tNo edits found')
                                else:
                                        vcfDict = vcf_merge(vcfDict, modelVcf)
                                        if z == 0:
                                                geneBlocks = geneblocks_update(geneBlocks, model, newModelCoords)
                                        # Verbose and log
                                        verbose_print(args, 'Model length is the same but sequence differs a bit. [' + model[3] + ']')
                                        log_comment(args, logName, '#' + model[3] + '\tModel length is the same\tOld=' + origProt + '\tNew=' + newProt + '\tAccepted edits=' + vcf_line_format(modelVcf))
                        else:
                                # Should we save this modification?
                                ## Check 1: Single substitution that made things worse? (This likely occurs for genes that lack good transcriptional support. Exon boundaries aren't supported well and a transcript alignment proposes an exon that has no frame that lacks stop codons)
                                if len(modelVcf[model[2]]) == 1 and '*' in list(modelVcf[model[2]].items())[0][1][0]:
                                        verbose_print(args, 'Found no edits [' + model[3] + ']')
                                        log_comment(args, logName, '#' + model[3] + '\tNo edits found')
                                # Log how much shorter the new model is
                                else:
                                        vcfDict = vcf_merge(vcfDict, modelVcf)
                                        if z == 0:
                                                geneBlocks = geneblocks_update(geneBlocks, model, newModelCoords)
                                        # Verbose and log
                                        if sizeDiff <= 10.0:
                                                verbose_print(args, 'I shortened this model, but not by much. It\'s probably skipping an exon [' + model[3] + ']')
                                                log_comment(args, logName, '#' + model[3] + '\tModel length decreased _slightly_\tOriginal model is ' + str(sizeDiff) + '% longer than new\tOld=' + origProt + '\tNew=' + newProt + '\tAccepted edits=' + vcf_line_format(modelVcf))
                                        elif sizeDiff <= 60.0:
                                                verbose_print(args, 'I shortened this model a lot. Was this gene a chimer, or does it involve ab initio predicted exons? [' + model[3] + ']')
                                                log_comment(args, logName, '#' + model[3] + '\tModel length decreased _quite a bit_\tOriginal model is ' + str(sizeDiff) + '% longer than new\tOld=' + origProt + '\tNew=' + newProt + '\tAccepted edits=' + vcf_line_format(modelVcf))
                                        else:
                                                verbose_print(args, 'I shortened this model _dramatically_. Check this out manually [' + model[3] + ']')
                                                log_comment(args, logName, '#' + model[3] + '\tModel length decreased _dramatically_\tOriginal model is ' + str(sizeDiff) + '% longer than new\tOld=' + origProt + '\tNew=' + newProt + '\tAccepted edits=' + vcf_line_format(modelVcf))
        # Create output VCF-like file
        if z == 0:
                vcf_output(args.outputFileName, vcfDict, '.')
                stage1Vcf = copy.deepcopy(vcfDict)
        else:
                if gmapRescueDict != {}:
                        # Pull out the new indel positions so we can distinguish between indels found in the main module versus the gene rescue module
                        newVcf = {}
                        for key, value in vcfDict.items():
                                if key in stage1Vcf:
                                        st1ContigVcf = stage1Vcf[key]
                                        for pair in value.items():
                                                if pair[0] not in st1ContigVcf:
                                                        if key not in newVcf:
                                                                newVcf[key] = {pair[0]: pair[1]}
                                                        else:
                                                                newVcf[key][pair[0]] = pair[1]
                                else:
                                        newVcf[key] = value
                        vcf_output(args.outputFileName, newVcf, '# Gene model rescue module indel predictions')

# Check for probable gene joins
if args.verbose or args.log:
        joins = gene_overlap_validation(geneBlocks)
        verbose_print(args, '### Probable gene merges ###\n' + '\n'.join(joins))
        log_comment(args, logName, '### Probable gene merges ###\n' + '\n'.join(joins))
        """This function is only capable of finding gene models that join over shared exons. Gene models that join 
        through introns should be discovered by EVM/PASA but it's beyond the capacity of this program to try to find
        these scenarios. This at least provides a nice way to validate these changes manually to see if things
        are going according to plan."""

# Remove the output file that was being replaced if relevant
if tmpFileName != None and os.path.isfile(tmpFileName):
        os.remove(tmpFileName)

#### SCRIPT ALL DONE
verbose_print(args, 'Program completed successfully!')
