#! /usr/env/python

import glob, errno, os, json, re, shutil
# Argument parser for user-inputted values, and a nifty help menu
from argparse import ArgumentParser

#Parser for arguments
parser = ArgumentParser(description='Automate tbl2asn to create beautiful .sqn files for your submission')
parser.add_argument('-v', '--version', action='version', version='%(prog)s v1.0')
parser.add_argument('-i', '--input', required=True, help='Specify input directory')

# Get the arguments into a list
args = vars(parser.parse_args())

# Define variables from the arguments - there may be a more streamlined way to do this
path = args['input']


def make_path(inPath):
    """from: http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary \
    does what is indicated by the URL"""
    try:
        os.makedirs(inPath)
        # os.chmod(inPath, 0777)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

# In order to make this script work properly, there are a few accessory files required. Examples of these files should
#  be included in the tbl2asnExamples directory
biosamples = open("biosamples.txt").readlines()
strainBiosample = {}

# Parse the biosamples into a dictionary
for sample in biosamples:
    data = sample.replace("\xc2\xa0", "").rstrip().split("\t")
    strainBiosample[data[0]] = data[1]


#Grab the sequence files
fasta = glob.glob("*.fsa")

genomeCoverage = open("genomeCoverages.txt").readlines()
coverages = {}

for coverage in genomeCoverage:
    cov = coverage.rstrip().split("\t")
    coverages[cov[0]] = cov[1]

comments = open("genome.asm").readlines()

# Read the template file into memory
template = open("template.sbt").readlines()


for fastaFile in fasta:
    name = fastaFile.split(".")[0]
    print strainBiosample[name]
    strainTemplate = template
    strainComments = comments
    customTemplate = open("%s.sbt" % name, "wb")
    customComment = open("%s.cmt" % name, "wb")
    for line in strainTemplate:
        # print line
        if not re.search("SUB800882", line):
            customTemplate.write(line)
            # print "!"
        else:
            spine = re.sub("SUB800882", str(strainBiosample[name]), line)
            customTemplate.write(spine)
    customTemplate.close()
    for com in strainComments:
        if not "Genome Coverage" in com:
            customComment.write(com)
        else:
            comm = "Genome Coverage\t%sx\n" % coverages[name]
            customComment.write(comm)
    make_path(name)
    shutil.move(fastaFile, name)
    shutil.move("%s.sbt" % name, name)
    shutil.move("%s.cmt" % name, name)

folders = glob.glob("*")
make_path("%s/sqnFiles/" % path)

print "Running tbl2asn"

for folder in folders:
    if os.path.isdir(folder) and not "zholding" in folder and not "tbl2asnOut" in folder and not "sqnFiles" in folder:
        print folder
        os.chdir("%s/%s" % (path, folder))
        make_path("%s/%s/tbl2asnOut" % (path, folder))
        command = "tbl2asn -p %s/%s -r %s/%s/tbl2asnOut/ -t %s.sbt -a s -X C " \
                  "-y 'Source DNA available from Burton Blais, 960 Carling Ave., Bldg. 22, Ottawa, Ontario," \
                  " Canada, K1C 0C6' -V vb -Z discrepancies.tx" % (path, folder, path, folder, folder)
        os.system(command)
        shutil.copy("%s/%s/tbl2asnOut/%s.sqn" % (path, folder, folder), "%s/sqnFiles" % path)
        os.chdir(path)



