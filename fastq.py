#!/usr/bin/env python
from accessoryfiles import *
from glob import glob
__author__ = 'adamkoziol'


class Fastq(object):

    def objectprep(self):
        """
        Create objects for each file to be renamed
        """
        # Open the file containing the files to be renamed
        with open(self.renamefile, 'rb') as renamefile:
            for row in renamefile:
                # Extract the current, and desired names from the row
                sample, strain = row.rstrip().split(',')
                # Create and populate a metadata object
                metadata = MetadataObject()
                metadata.name = sample
                metadata.strain = strain
                # Append the metadata to the list of sample metadata
                self.samples.append(metadata)
        # Rename the files
        self.fastq()

    def fastq(self):
        """Renames fastq files to match assembly names"""
        import re
        for sample in self.samples:
            # Get the names of the fastq files
            sample.fastq = sorted(glob('{}{}*'.format(self.sequencepath, sample.name)))
            try:
                # Find the extension using a regex search; should be able to find both .fastq and .fastq.gz extensions
                sample.ext = re.search(r'.+?(\..+)', sample.fastq[0]).group(1)
                # Set the forward and reverse names - essentially ensures that both files are present. If not, then
                # the try/except loop will be triggered
                sample.forward, sample.reverse = sample.fastq[0], sample.fastq[1]
                # Map the file renaming function with the old name and the new name, which includes the 'R1/R2'
                # directionality as the second input of the lambda function
                map(lambda x, y: os.rename(x, '{}{}_{}{}'.format(self.sequencepath, sample.strain, y, sample.ext)),
                    sample.fastq, ['R1', 'R2'])
            # Except missing files
            except IndexError:
                # Try to find the files if they were already renamed
                sample.fastq = sorted(glob('{}{}*'.format(self.sequencepath, sample.name)))
                try:
                    # Test to see if the list contains two strings
                    sample.fastq[0], sample.fastq[1]
                except IndexError:
                    # If this strain is truly missing or misnamed, raise the error
                    print 'I can\'t find {}. Have the files already been renamed?'.format(sample.name)
                    raise

    def __init__(self, args):
        from time import time
        self.start = time()
        # Initialise the variables
        self.path = os.path.join(args.path, '')
        # If path information is not supplied in the -f argument, assume that the file is in the path
        self.renamefile = args.renamefile
        if '/' not in self.renamefile:
            self.renamefile = self.path + self.renamefile
        assert os.path.isfile(self.renamefile), u'Invalid file supplied for renaming file: {}' \
            .format(self.renamefile)
        self.sequencepath = os.path.join(args.sequencepath, '')
        assert os.path.isdir(self.sequencepath), u'Sequence path  is not a valid directory {0!r:s}' \
            .format(self.sequencepath)
        # A list to store metadata samples
        self.samples = list()
        # Run the analyses
        self.objectprep()

# If the script is called from the command line, then call the argument parser
if __name__ == '__main__':
    from argparse import ArgumentParser
    # Parser for arguments
    parser = ArgumentParser(description='Renames all .fastq files with desired names')
    parser.add_argument('path',
                        help='Folder in which analyses are to be performed')
    parser.add_argument('-s', '--sequencepath',
                        required=True,
                        help='Specify path of folder containing your fastq files')
    parser.add_argument('-f', '--renamefile',
                        required=True,
                        help='A comma-separated list with the sample and strain names '
                             '(one sample name/strain name/organism pair per row): '
                             '2015-SEQ-0947,OLF15251\n'
                             '2015-SEQ-0948,OLC2218\n'
                             'If the file is in the path, just provide the file name, but if it is in a different '
                             'folder, the path and file name must be supplied.')
    # Get the arguments into an object
    arguments = parser.parse_args()
    # Run the function
    Fastq(arguments)
