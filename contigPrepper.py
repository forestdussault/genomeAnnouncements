#!/usr/bin/env python
from glob import glob
from accessoryfiles import *
from tbl2asn import Tbl2asn

__author__ = 'adamkoziol'


def strainname(fullname):
    """
    Removes the path (os.path.split(fullname)[1] and the extension (.split(".")[0]
    :param fullname: full path + file name + file extension of a file
    """
    name = os.path.split(fullname)[1].split(".")[0]
    return name


class ContigPrepper(object):

    def populateobject(self):
        """Populate :self.strainorganism"""
        from csv import DictReader, Sniffer
        # Grab all the contigs
        self.fastafiles = sorted(glob('{}*.fa*'.format(self.path)))
        # Set the header information for reading the organism file into a dictionary
        fieldnames = ['sample', 'strain', 'organism', 'serotype', 'coverage']
        with open(self.organismfile, 'rb') as orgfile:
            # Detect the delimiter used in the file
            dialect = Sniffer().sniff(orgfile.read(1024))
            # Go back to the beginning of the file
            orgfile.seek(0)
            # Open the sequence profile file as a dictionary
            orgdict = DictReader(orgfile, fieldnames=fieldnames, dialect=dialect)
            # Iterate through the rows in the dictionary
            for row in orgdict:
                for fastafile in self.fastafiles:
                    # Create the metadata object
                    metadata = MetadataObject()
                    # Remove the path and the extension from the fasta file
                    filename = strainname(fastafile)
                    # If the name of the file matches the supplied file name
                    if row['sample'] == filename:
                        # Set the name, file name and path, as well as the organism
                        metadata.name = filename
                        metadata.fastafile = fastafile
                        metadata.strain = row['strain']
                        metadata.organism = row['organism']
                        metadata.serotype = row['serotype']
                        metadata.coverage = row['coverage'] if 'x' in row['coverage'].lower() else row['coverage'] + 'X'
                        # Append the metadata to the list of sample metadata
                        self.samples.append(metadata)
            # Assert that all the files are present in the organism file, and that the names match
            missingfiles = [fastafile for fastafile in self.fastafiles
                            if fastafile not in [sample.fastafile for sample in self.samples]]
            assert missingfiles == [], 'The following files differ from the organism file: \n{}'\
                .format('\n'.join(missingfiles))

    def prepcontigs(self):
        """Rename the contigs to NCBI format with on the organism"""
        from Bio import SeqIO
        #
        for sample in self.samples:
            # Reset the contig count to 1 for each sample
            contigcount = 1
            # Create the path for the reformatted file
            sample.reformattedpath = self.path + sample.strain
            make_path(sample.reformattedpath)
            sample.reformattedfile = '{}/{}.fsa'.format(sample.reformattedpath, sample.strain)
            if not os.path.isfile(sample.reformattedfile):
                print sample.name
                with open(sample.reformattedfile, 'wb') as formattedfile:
                    for record in SeqIO.parse(open(sample.fastafile, 'rU'), 'fasta'):
                        # Set the new record.id
                        record.id = '{}_Cont{:04d} [organism={}] [serotype={}] [strain={}]'.format(sample.strain,
                                        contigcount, sample.organism, sample.serotype, sample.strain)
                        # Clear the name and description attributes of the record
                        record.name = ''
                        record.description = ''
                        # Write each record to the combined file
                        SeqIO.write(record, formattedfile, 'fasta')
                        # Increment the contig count
                        contigcount += 1

    def fastq(self):
        """Renames fastq files to match assembly names"""
        import re
        print
        for sample in self.samples:
            # Get the names of the fastq files
            sample.fastq = sorted(glob('{}{}*'.format(self.sra, sample.name)))
            try:
                # Find the extension using a regex search; should be able to find both .fastq and .fastq.gz extensions
                sample.extension = re.search(r'.+?(\..+)', sample.fastq[0]).group(1)
                # Set the forward and reverse names - essentially ensures that both files are present. If not, then
                # the try/except loop will be triggered
                sample.forward, sample.reverse = sample.fastq[0], sample.fastq[1]
                # Map the file renaming function with the old name and the new name, which includes the 'R1/R2'
                # directionality as the second input of the lambda function
                map(lambda x, y: os.rename(x, '{}{}_{}{}'.format(self.sra, sample.strain, y, sample.extension)),
                    sample.fastq, ['R1', 'R2'])
            # Except missing files
            except IndexError:
                print sample.name
                # Try to find the files if they were already renamed
                sample.fastq = sorted(glob('{}{}*'.format(self.sra, sample.strain)))
                try:
                    # Test to see if the list contains two strings
                    sample.fastq[0], sample.fastq[1]
                except IndexError:
                    # If this strain is truly missing or misnamed, raise the error
                    print sample.strain
                    raise

    def __init__(self, args):
        from time import time
        self.start = time()
        # Initialise the variables
        self.path = os.path.join(args.path, '')
        self.organismfile = args.organismfile
        # If path information is not supplied in the -f argument, assume that the file is in the path
        if '/' not in self.organismfile:
            self.organismfile = self.path + self.organismfile
        assert os.path.isfile(self.organismfile), 'Invalid file supplied for organism file: {}'\
            .format(self.organismfile)
        self.commentfile = args.commentfile
        assert os.path.isfile(self.commentfile), 'Invalid file supplied for comment file: {}' \
            .format(self.commentfile)
        self.templatefile = args.templatefile
        assert os.path.isfile(self.commentfile), 'Invalid file supplied for template file: {}' \
            .format(self.template)
        # Initialise variables
        self.fastafiles = []
        self.samples = []
        # Populate :self.samples with strain name, file name, and organism
        self.populateobject()
        # Prep the contigs
        self.prepcontigs()
        # Run the tbl2asn script
        Tbl2asn(self)
        if args.sra:
            self.sra = os.path.join(args.sra, '')
            assert os.path.isdir(self.sra), 'Invalid path supplied for fastq files: {}' \
                .format(self.sra)
            self.fastq()
        # Print an exit statement
        print "\nElapsed Time: {:.2f} seconds".format(time() - self.start)


# If the script is called from the command line, then call the argument parser
if __name__ == '__main__':
    from argparse import ArgumentParser
    # Parser for arguments
    parser = ArgumentParser(description='Reformats all the headers in fasta files in a directory to be compatible for'
                                        'NCBI submissions. The organism must be supplied for each sample.'
                                        '>2015-SEQ-0947_19_length_6883_cov_11.9263_ID_37 becomes'
                                        '>Cont0001 [organism=Listeria monocytogenes] [strain=OLF15251]'
                                        'Reformatted files will be created in: path/reformatted')
    parser.add_argument('path',  help='Specify path of folder containing your fasta files')
    parser.add_argument('-f', '--organismfile',
                        required=True,
                        help='A comma-separated list with the sample, strain, and organism names, and the coverage  '
                             '(one sample name/strain name/organism pair per row):'
                             '2015-SEQ-0947,OLF15251,Listeria monocytogenes,27.25\n'
                             '2015-SEQ-0948,OLC2218,Escherichia coli,18.70\n'
                             'If the file is in the path, just provide the file name, but if it is in a different '
                             'folder, the path and file name must be supplied.')
    parser.add_argument('-c', '--commentfile',
                        required=True,
                        help='A file of comments used in preparing .sqn files. The file looks something like this:'
                             'Source DNA available from :name :address '
                             'StructuredCommentPrefix	##Genome-Assembly-Data-START## '
                             'Assembly Method	:method '
                             'Genome Coverage '
                             'Sequencing Technology	:platform')
    parser.add_argument('-t', '--templatefile',
                        required=True,
                        help='A template for tbl2asn to use. Generate the template at: '
                             'https://submit.ncbi.nlm.nih.gov/genbank/template/submission/')
    parser.add_argument('-s', '--sra',
                        help='Path to a folder containing the fastq files used to make the assemblies. These fastq'
                             'files will be renamed to match assembly names')
    # Get the arguments into an object
    arguments = parser.parse_args()
    # Run the function
    ContigPrepper(arguments)
