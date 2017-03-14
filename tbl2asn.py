#! /usr/env/python
from accessoryfiles import *
import shutil


class Tbl2asn(object):

    def populatecomments(self):
        """
        Load comment and template files. Create a comment file for each strain, and modify this file with the
        calculated average coverage
        """
        printtime('Populating comments', self.start)
        # Load the accessory files into memory
        self.comments = open(self.commentfile).readlines()
        self.template = open(self.templatefile).readlines()
        for sample in self.samples:
            # Define the strain-specific comment and template files
            sample.commentfile = '{}/{}.cmt'.format(sample.reformattedpath, sample.name)
            # Write the strain-specific comment file
            with open(sample.commentfile, 'wb') as commentfile:
                for comment in self.comments:
                    # The genome coverage needs to be populated
                    if "Genome Coverage" in comment:
                        commentfile.write('Genome Coverage\t{}\n'.format(sample.coverage))
                    else:
                        commentfile.write(comment)

    def tbl2asnthreads(self):
        """
        Run tbl2asn on each of the strains in order to create .sqn files
        """
        from threading import Thread
        printtime('Creating .sqn files', self.start)
        for i in range(len(self.samples)):
            threads = Thread(target=self.tbl2asn, args=())
            threads.setDaemon(True)
            threads.start()
        for sample in self.samples:
            self.queue.put(sample)
        self.queue.join()

    def tbl2asn(self):
        from subprocess import call
        import os
        while True:
            sample = self.queue.get()
            # Set the output file name, and the system call
            sample.sqnfile = '{}sqnfiles/{}.sqn'.format(self.path, sample.name)
            # Include the generic comment regarding the availability of the source DNA
            sample.tbl2asncommand = 'tbl2asn -p {} -r {} -t {} -a s -X C -V v -Z {}/discrepancies.txt -y ' \
                                    '"Source DNA available from Burton Blais, 960 Carling Ave., Bldg. 22, Ottawa, ' \
                                    'Ontario, Canada, K1A 0C6"'.format(sample.reformattedpath, self.sqnpath,
                                                                       self.templatefile, sample.reformattedpath)
            # If the .sqn file has not already been created, run tbl2asn now
            if not os.path.isfile(sample.sqnfile):
                call(sample.tbl2asncommand, shell=True, stdout=self.devnull, stderr=self.devnull)
                # Move the .val file to the appropriate subfolder
                shutil.move('{}{}.val'.format(self.sqnpath, sample.name), sample.reformattedpath)
                try:
                    # Remove the errorsummary file from the sequence path
                    os.remove('{}errorsummary.val'.format(self.sqnpath))
                except OSError:
                    pass
            self.dotter.dotter()
            self.queue.task_done()

    def __init__(self, inputobject):
        from Queue import Queue
        # Initialise the variables
        self.start = inputobject.start
        self.path = inputobject.path
        self.commentfile = inputobject.commentfile
        self.templatefile = inputobject.templatefile
        self.samples = inputobject.samples
        self.dotter = inputobject.dotter
        self.comments = str()
        self.template = str()
        self.sqnpath = os.path.join(self.path, 'sqnfiles', '')
        make_path(self.sqnpath)
        self.queue = Queue()
        self.devnull = open(os.devnull, 'wb')
        # Run the analyses
        self.populatecomments()
        self.tbl2asnthreads()
