BOTUX
=====

Prerequisites:
-----
 Python version 2.7 (version 3.0 or greater may not work)  
 Python module HTSeq .. installation instructions:http://www-huber.embl.de/users/anders/HTSeq/doc/install.html  

Usage:
-----

usage: runBOTUX.py [-h] [-v] [-fasta <in.fasta>] [-fastq <read1.fastq>]
                   [-2 <read2.fastq>] [-wl <8>] -project <P1>
                   [-loadModel <inModel.pkl>] [-ptl <0.75>] [-t <0.75>]

BOTUX: Bayesian-like OTU eXaminer

	Author: 	Vishal N Koparde, Ph. D.
	Created:	130315
	Modified:	130330

optional arguments:
  -h, --help            show this help message and exit  
  -v, --version         show program's version number and exit  
  -fasta <in.fasta>     input fasta file  
  -fastq <read1.fastq>  input fastq file(read1)  
  -2 <read2.fastq>      input fastq file(read2)  
  -wl <8>               word lenth (default 8)  
  -project <P1>         project name .. prefix used to name output files  
  -loadModel <inModel.pkl>
                        load OTU model from python pickle  
  -ptl <0.75>           percentile trimming length  
  -t <0.75>             threshold score .. advanced usage  
