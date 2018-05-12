#! python3

# deGRIT_edit (DEtection and Genomic Rectification of Indels using Transcripts - edit)

# After running the main deGRIT program to produce a VCF-like file, this script will
# produce a modified genome FASTA file.

# Load packages
import os, argparse
from Bio import SeqIO

# Define VCF-like file handling
def vcf_parser(args):
        vcfDict = {}
        with open(args.editFile, 'r') as fileIn:
                for line in fileIn:
                        # Skip unnecessary lines or stop parsing if rescue_genes indels are not wanted
                        if line.startswith('#contig_id'):
                                continue
                        elif line.startswith('# Gene rescue module indel predictions'):
                                if args.rescue_genes == False:
                                        break
                                else:
                                        continue
                        # Pull out line details
                        sl = line.rstrip('\n').split('\t')
                        if sl[0] not in vcfDict:
                                vcfDict[sl[0]] = [[int(sl[1]), sl[2]]]
                        else:
                                vcfDict[sl[0]].append([int(sl[1]), sl[2]])
        return vcfDict

def vcf_edit(inputVcf, contigID, genomeSeq, editEvents):
        # Format edit positions list
        vcfList = inputVcf[contigID]
        vcfList.sort(reverse=True)
        # Edit the genome sequence
        for pair in vcfList:
                indelIndex = pair[0] - 1                                        # - 1 to make this act 0-based (in the main program we instead minused coordRange which accomplished the same goal of making the index 0-based).
                if pair[1] == '.':
                        genomeSeq = genomeSeq[:indelIndex] + genomeSeq[indelIndex+1:]
                        editEvents[0] += 1
                elif '*' in pair[1]:                                            # This makes a substitution in our genome by deleting the original base and replacing it with pair[1][:-1] (the last character is an asterisk so we remove it).
                        genomeSeq = genomeSeq[:indelIndex] + pair[1][:-1] + genomeSeq[indelIndex+1:]
                        editEvents[1] += 1
                else:
                        genomeSeq = genomeSeq[:indelIndex] + pair[1] + genomeSeq[indelIndex:]
                        editEvents[2] += 1
        return genomeSeq, editEvents

# Define file output function
def genome_output(args, seqid, sequence):
        if not os.path.isfile(args.outputFileName):
                open(args.outputFileName, 'w').close()
        with open(args.outputFileName, 'a') as fileOut:
                # Format multiline output if necessary
                if args.multilineLength != -1:
                        sequence = '\n'.join([sequence[i:i+args.multilineLength] for i in range(0, len(sequence), args.multilineLength)])
                # Write to file
                fileOut.write('>' + seqid + '\n' + sequence + '\n')

# Define argument validation
def validate_args(args):
        # Validate input file locations
        if not os.path.isfile(args.genomeFile):
                print('I am unable to locate the input genome fasta file (' + args.genomeFile + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        # Validate multiline length if provided
        if args.multilineLength != -1:
                if args.multilineLength <= 0:
                        print('You specified a number less than or equal to 0 for the multiline length. Make this an integer greater than 0 or leave this argument blank, and try again.')
                        quit()
        # Handle file overwrites
        if os.path.isfile(args.outputFileName):
                print(args.outputFileName + ' already exists. Either provide the -fo argument to this program or delete/move/rename this file and run the program again.')
                quit()

# Define final results info presentation
def results_info(args, editEvents, numContigs):
        # Present information of edit events
        print('Genome editing complete. Basic information provided below.')
        if args.rescue_genes == True:
                print('Edits include gene module rescue positions: Yes')
        else:
                print('Edits include gene module rescue positions: No')
        print('Number of deletion events: ' + str(editEvents[0]))
        print('Number of substitution events: ' + str(editEvents[1]))
        print('Number of insertion events: ' + str(editEvents[2]))
        # Present information of contigs edited
        print('Total number of contigs: ' + str(numContigs[0]))
        print('Number of edited contigs: ' + str(numContigs[1]))

### USER INPUT
usage = """%(prog)s is a utility program to perform the final genome editing process
using the output of the main deGRIT program.
"""
# Reqs
p = argparse.ArgumentParser(description=usage)

p.add_argument("-gen", "--genomefile", dest="genomeFile",
               help="Input genome contig fasta file name")
p.add_argument("-e", "--editFile", dest="editFile",
               help="Input VCF-like edit file name")
p.add_argument("-o", "--output", dest="outputFileName",
               help="Output edited genome file name")
# Opts
p.add_argument("-ml", "--multi-line", dest="multilineLength", type = int,
               help="Optionally produce a multi-line formatted FASTA file (specify the length of each line)", default = -1)
p.add_argument('-r', '--skip_rescue_genes', dest="rescue_genes", action='store_false',
               help="Optionally turn OFF the inclusion of edits derived from gene rescue module (do not provide this argument if you want all changes to be made)", default=True)

args = p.parse_args()

# Validate arguments
validate_args(args)

# Load genome file for iteration
genomeRecords = SeqIO.parse(open(args.genomeFile, 'r'), 'fasta')

# Parse the VCF-like file
vcfDict = vcf_parser(args)

### CORE PROCESS
# Perform the genome editing
editEvents = [0,0,0]                                                            # Provide some information to the user here for convenience [format is deletions, substitutions, insertions].
numContigs = [0,0]                                                              # Format is [total contig num, edited contig num]
for record in genomeRecords:
        seqid = record.id
        sequence = str(record.seq)
        if seqid in vcfDict:
                sequence, editEvents = vcf_edit(vcfDict, seqid, sequence, editEvents)
                genome_output(args, seqid, sequence)
                numContigs[0] += 1
                numContigs[1] += 1
        else:
                genome_output(args, seqid, sequence)
                numContigs[0] += 1

# Give some results info
results_info(args, editEvents, numContigs)

#### SCRIPT ALL DONE