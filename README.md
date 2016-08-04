# genomeAnnouncements

## Requirements

* tbl2asn
* Biopython

 Installation
 
 `sudo apt-get update`
 `sudo apt-get install -y build-essential git python-dev python-pip wget` 

Clone the repository:

`https://github.com/adamkoziol/genomeAnnouncements.git`

Install python dependencies:

	
```
cd genomeAnnouncements/
python setup.py install
```


Install tbl2asn

`wget ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/linux64.tbl2asn.gz`
`gunzip linux64.tbl2asn.gz`
`sudo mv linux64.tbl2asn /usr/bin/tbl2asn`


## contigPrepper

Requires accessory files in order to run.

1) A template for tbl2asn to use. Generate the template at: https://submit.ncbi.nlm.nih.gov/genbank/template/submission/
2) A comma-separated list containing the sample name, strain name, file name (no extension), organism, and average coverage

Excerpt from example organism file:
organismdata.csv
```
D2653,OLC996,2013-SEQ-0123_2014-SEQ-0279,Escherichia coli,O157:H7,93.8052X
```
3) A comment file

Example comments file:
comments.asm
```
StructuredCommentPrefix	##Genome-Assembly-Data-START##
Assembly Method	SPAdes v. 3.7.1
Genome Coverage
Sequencing Technology	Illumina MiSeq
```

### usage

Example command (from within the genomeAnnouncements folder):

`python contigPrepper.py -f /path/to/organismfile.csv -c /path/to/commentfile.asm -t /path/to/templatefile.sbt -s /path/to/fastqfiles /path`

```
usage: contigPrepper.py [-h] -f ORGANISMFILE -c COMMENTFILE -t TEMPLATEFILE
                        [-s SRA]
                        path

Reformats all the headers in fasta files in a directory to be compatible
for NCBI submissions. The organism must be supplied for each sample.
>2015-SEQ-0947_19_length_6883_cov_11.9263_ID_37 becomes
>Cont0001 [organism=Listeria monocytogenes] [strain=OLF15251]
Reformatted files will be created in:
path/reformatted

positional arguments:
  path                  Specify path of folder containing your fasta files

optional arguments:
  -h, --help            show this help message and exit
  -f ORGANISMFILE, --organismfile ORGANISMFILE
                        A comma-separated list with the sample, strain, and
                        organism names, and the coverage (one sample
                        name/strain name/organism pair per
                        row):2015-SEQ-0947,OLF15251,Listeria
                        monocytogenes,27.25 2015-SEQ-0948,OLC2218,Escherichia
                        coli,18.70 If the file is in the path, just provide
                        the file name, but if it is in a different folder, the
                        path and file name must be supplied.
  -c COMMENTFILE, --commentfile COMMENTFILE
                        A file of comments used in preparing .sqn files. The
                        file looks something like this:Source DNA available
                        from :name :address StructuredCommentPrefix ##Genome-
                        Assembly-Data-START## Assembly Method :method Genome
                        Coverage Sequencing Technology :platform
  -t TEMPLATEFILE, --templatefile TEMPLATEFILE
                        A template for tbl2asn to use. Generate the template
                        at: https://submit.ncbi.nlm.nih.gov/genbank/template/s
                        ubmission/
  -s SRA, --sra SRA     Path to a folder containing the fastq files used to
                        make the assemblies. These fastq files will be renamed
                        to match assembly names
```
