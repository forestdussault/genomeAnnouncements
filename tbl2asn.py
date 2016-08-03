#! /usr/env/python
from accessoryfiles import *
from glob import glob
import json
import re
import shutil


class Tbl2asn(object):

    def populatecomments(self):
        # Load the accessory files into memory
        self.comments = open(self.commentfile).readlines()
        self.template = open(self.templatefile).readlines()
        for sample in self.samples:
            # Define the strain-specific comment and template files
            sample.commentfile = '{}/{}.cmt'.format(sample.reformattedpath, sample.strain)
            # sample.templatefile = '{}reformatted/{}.sbt'.format(self.path, sample.strain)
            # Write the strain-specific comment file
            with open(sample.commentfile, 'wb') as commentfile:
                for comment in self.comments:
                    #
                    if "Genome Coverage" in comment:
                        commentfile.write('Genome Coverage\t{}\n'.format(sample.coverage))
                    else:
                        commentfile.write(comment)
            # Write the strain-specific template file
            # with open(sample.templatefile, 'wb') as templatefile:

    def tbl2asnthreads(self):
        from threading import Thread
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
        # for sample in self.samples:
        while True:
            sample = self.queue.get()
            sample.sqnfile = '{}sqnfiles/{}.sqn'.format(self.path, sample.strain)
            sample.tbl2asncommand = 'tbl2asn -p {} -r {} -t {} -a s -X C -V v -Z {}/discrepancies.txt -y ' \
                                    '"Source DNA available from Burton Blais, 960 Carling Ave., Bldg. 22, Ottawa, ' \
                                    'Ontario, Canada, K1A 0C6"'.format(sample.reformattedpath, self.sqnpath,
                                                                       self.templatefile, sample.reformattedpath)
            if not os.path.isfile(sample.sqnfile):
                call(sample.tbl2asncommand, shell=True, stdout=self.devnull, stderr=self.devnull)
                shutil.move('{}{}.val'.format(self.sqnpath, sample.strain), sample.reformattedpath)
                try:
                    os.remove('{}errorsummary.val'.format(self.sqnpath))
                except OSError:
                    pass
            dotter()
            self.queue.task_done()

    def __init__(self, inputobject):
        from time import time
        from Queue import Queue
        # Initialise the variables
        self.start = inputobject.start
        self.path = inputobject.path
        self.commentfile = inputobject.commentfile
        self.templatefile = inputobject.templatefile
        self.samples = inputobject.samples
        # Initialise variables
        self.comments = ''
        self.template = ''
        self.sqnpath = self.path + 'sqnfiles/'
        make_path(self.sqnpath)
        self.queue = Queue()
        self.devnull = open(os.devnull, 'wb')
        self.populatecomments()
        self.tbl2asnthreads()
