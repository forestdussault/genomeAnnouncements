#!/usr/bin/env python
from accessoryfiles import *
__author__ = 'adamkoziol'


class Remove(object):

    def runner(self):
        """
        Run the required methods in the appropriate order
        """
        self.report()
        self.remove()
        self.dotter.globalcounter()

    def report(self):
        """
        Parses the default format of NCBI contamination reports to extract sequences that should be excluded from the
        assemblies, and sequences that need to be trimmed
        """
        printtime('Parsing report', self.start)
        # Parse the required information from the report
        # from https://stackoverflow.com/questions/7559397/python-read-file-from-and-to-specific-lines-of-text
        with open(self.contaminationreport, 'rb') as report:
            # Since the reports have a consistent format, anything before a line of 'Exclude:' can be ignored
            for line in report:
                if line.rstrip() == 'Exclude:':
                    break
            # The block of text of with the sequences to include is separated from other sections by an empty line
            # when this line is encountered, stop populating the exclude dictionary
            for line in report:
                if line == '\n':
                    break
                # There is a header for the table of contigs to exclude - ignore this it consists of the following:
                # Sequence name, length, apparent source
                if 'Sequence' not in line:
                    # Split the data on tabs. It should look like this:
                    # CA_CFIA-1709_Cont0124	1439	pBR322
                    data = line.split('\t')
                    # Populate the dictionary with the strain name (data[0].split('_Cont')[0]) as the key, and the
                    # contig name ('Cont' + data[0].split('_Cont')[1]) as the value
                    try:
                        self.exclude[data[0].split('_Cont')[0]].update('Cont' + data[0].split('_Cont')[1])
                    # If the strain name hasn't been initialised in the dictionary, then it cannot be updated.
                    # Initialise it now
                    except KeyError:
                        self.exclude[data[0].split('_Cont')[0]] = 'Cont' + data[0].split('_Cont')[1]
            # Similar as excluding contigs above, except more detail from the tables are extracted. These data are
            # not yet used, but will hopefully be useful in later iterations of this script
            for line in report:
                if line.rstrip() == 'Trim:':
                    break
            # Stop populating dictionaries after an empty line
            for line in report:
                if line == '\n':
                    break
                # This table also has a header, as well as contigs that are safe to keep. Only the contigs marked
                # as '(trim)' should either be masked (future versions of the script), or trimmed.
                """After you remove the contamination, trim any Ns at the ends of the sequence and remove any sequences
                that are shorter than 200 nt and not part of a multi-component scaffold."""
                # Ignore any line that does not contain '(trim)' for the parsing
                # Sequence name, length, span(s), apparent source
                # CA_CFIA-728_Cont0027	69226	54803..55115	Rattus norvegicus (ignore)
                # D2435_Cont0079	1487	1..428,1410..1487	Gallus gallus (trim)
                if '(trim)' in line:
                    # Split on tabs
                    data = line.split('\t')
                    strainname = data[0].split('_Cont')[0]
                    try:
                        # Populate two separate dictionaries with the strain name. The trim dictionary contains the
                        # contig name, and a string of sequence to trim (e.g. 1..428,1410..1487), while the contigs
                        # dictionary contains the total length of the contig
                        self.trim[strainname].append(('Cont' + data[0].split('_Cont')[1], data[2]))
                        self.contigs[strainname] = data[1]
                    # Initialise and populate the list of tuples as required
                    except KeyError:
                        self.trim[strainname] = list()
                        self.trim[strainname].append(('Cont' + data[0].split('_Cont')[1], data[2]))

    def remove(self):
        """
        Remove the contaminated contigs from the assemblies
        """
        from Bio import SeqIO
        import re
        printtime('Removing contigs as required', self.start)
        for sample in self.samples:
            # Initialise lists to store contigs to exclude from the final assembly, and SeqIO records to be used
            # in preparing these final assemblies
            include = list()
            exclude = list()
            for name, data in self.exclude.items():
                if sample.name == name:
                    #
                    for record in SeqIO.parse(open(sample.reformattedfile, "rU"), "fasta"):
                        #
                        if data == 'Cont' + record.id.split('_Cont')[1]:
                            # print record.id, 'exc'
                            exclude.append(record.id)
            # Perform similar analyses on sequenced identified as requiring trimming. Right now, this only removes the
            # contigs as in the excluded section above, but it may become necessary to mask sequences if the contig
            # is a large one
            for name, data in self.trim.items():
                #
                if sample.name == name:
                    #
                    for contigs, locations in data:
                        #
                        for record in SeqIO.parse(open(sample.reformattedfile, "rU"), "fasta"):
                            #
                            if contigs == 'Cont' + record.id.split('_Cont')[1]:
                                #
                                exclude.append(record.id)
            # Variable to store the contig number - incremented so that each contig has a unique value. Since
            # contigs may be deleted from the file, this is designed to make the contig numbers sequential again
            contig = 1
            # Parse the reformatted file again
            for record in SeqIO.parse(open(sample.reformattedfile, "rU"), "fasta"):
                # Only process the contigs that should not be excluded
                if record.id not in exclude:
                    # Use regex to swap out the old contig number with the new, properly formatted number
                    description = re.sub("Cont\d+", 'Cont{}'.format('{:04d}'.format(contig)), record.description)
                    # Replace the id attribute with the modified description
                    record.id = description
                    # Clear the name and description attributes of the record
                    record.name = ''
                    record.description = ''
                    # Add the record to the list of all records to include for this strain
                    include.append(record)
                    # Increment the contig number
                    contig += 1
            # If the list of contigs passing is empty (extremely unlikely), then something has gone wrong
            if include:
                # Open the assembly file
                with open(sample.reformattedfile, 'wb') as removed:
                    # Write the records in the list to the file
                    SeqIO.write(include, removed, 'fasta')
            self.dotter.dotter()

    def __init__(self, inputobject):
        self.path = inputobject.path
        self.start = inputobject.start
        self.samples = inputobject.samples
        self.contaminationreport = inputobject.contaminationreport
        self.dotter = inputobject.dotter
        self.exclude = dict()
        self.trim = dict()
        self.contigs = dict()
        self.runner()
