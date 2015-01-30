#! /usr/env/python

import glob, sys, os, errno, shutil

path = "/media/nas/akoziol/GeneSipping/NCBI_submission/Sequences/2014-11-07"
# path = "/media/nas/akoziol/Collaborations/Ashley"

os.chdir(path)

def make_path(inPath):
    """from: http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary \
    does what is indicated by the URL"""
    try:
        os.makedirs(inPath)
        # os.chmod(inPath, 0777)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

files = glob.glob("*.fa*")

# print files

for fasta in files:
    fileName = fasta.split(".")[0]
    print fileName
    make_path("%s/%s" % (path, fileName))
    # os.chdir("%s/%s" % (path, fileName)) -e 700
    blastCall = "blastall -p blastn -i %s -d /media/nas/akoziol/UniVec/IlluminaPrimers.fa " \
                "-q -5 -G 3 -E 3 -F 'm D' -e 0.1 -Y 1.75e12 -m 8 -o %s/%s/%s_Illuminafmt8.tsv" % (fasta, path, fileName, fileName)
    os.system(blastCall)
    make_path("blastOutputs")
    shutil.copy("%s/%s/%s_Illuminafmt8.tsv" % (path, fileName, fileName), "%s/blastOutputs" % path)

    # os.chdir(path)