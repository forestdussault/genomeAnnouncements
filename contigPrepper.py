#! /usr/env/python

import os, glob, errno
from Bio import SeqIO


def make_path(inPath):
    """from: http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary \
    does what is indicated by the URL"""
    try:
        os.makedirs(inPath)
        # os.chmod(inPath, 0777)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

# Either save the python file in the directory with the files, or provide the appropriate path below
# path = "/media/nas/akoziol/GeneSipping/NCBI_submission/Sequences/allSequences"
path = "/media/nas/akoziol/GeneSipping/NCBI_submission/Sequences/renamingContigs"
# path = os.getcwd()

# Go to the path
os.chdir(path)

# Grab all the contigs
fastaFile = glob.glob("*.fa*")


def prepper(fastaFile):
    """All the contig prepper stuff in one place"""
    for fasta in fastaFile:
        # Initialise variables for use in the function
        # MAYBE in the future, I will use argparse to get arguments for the organism name.
        # until then, edit the organism name below as appropriate using the following format:
        # organism = "[organism=Enterobacter cloacae]"
        organism = ""
        updated = []
        # Strip off the file extension from the file name
        name = fasta.split(".")[0]
        print name
        strain = '[strain=%s]' % name
        # This block is custom for my GeneSippr NCBI submission
        # Ignore it for now. Or not. Whatever
        if not "OLC-1682" in name and not "OLC-1683" in name:
            organism = "[organism=Escherichia coli]"
        else:
            organism = "[organism=Enterobacter cloacae]"

        # Using BioPython - open fasta file, and extract each record (header + sequence)
        for record in SeqIO.parse(open(fasta, "rU"), "fasta"):
            # Headers should have a format similar to the following:
            # >2014-SEQ-900_10_length_10013_cov_83.2596_ID_19
            # Splitting at "_"s allow for the extraction of the contig number (header[1])
            header = record.id.split("_")

            # Make a new ID with the contig number formatted with up to three floating zeroes, as well as the organism
            newID = "Cont%04d %s %s" % (float(header[1]), organism, strain)
            # Replace the old record.id
            record.id = newID
            # These need to be left blank (or properly filled), or else the script gets angry
            record.name = ''
            record.description = ''
            # Append all the updated records to updated
            updated.append(record)
        # Make the appropriate path
        make_path("%s/updatedHeaders" % path)
        # Open a new file for the updated headers
        fileName = "%s/updatedHeaders/%s.fa" % (path, name)
        formatted = open(fileName, "wb")
        # Write to the file using SeqIO - this formats the output into proper fasta format
        SeqIO.write(updated, formatted, "fasta")
        # Close the file
        formatted.close()

# Run the function
prepper(fastaFile)