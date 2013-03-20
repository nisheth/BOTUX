BOTUX
=====

usage: runBOTUX.py [-h] [-v] [-fasta <in.fasta>] [-fastq <read1.fastq>]
                   [-2 <read2.fastq>] [-wl <8>] [-otuFasta <OTUs.fasta>]
                   [-loadModel <inModel.pkl>] [-saveModel <outModel.pkl>]
                   [-ptl <0.9>] [-pa <OTU_assignments.txt>]
                   [-pp <OTU_profiles.txt>] [-t <0.75>]

BOTUX: Bayesian-like OTU eXaminer
	
	Author: 	Vishal N Koparde, Ph. D. 
	Created:	130315
	Modified:	130320

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -fasta <in.fasta>     input fasta file
  -fastq <read1.fastq>  input fastq file(read1)
  -2 <read2.fastq>      input fastq file(read2)
  -wl <8>               word lenth (default 8)
  -otuFasta <OTUs.fasta>
                        output fasta file containing OTU seeds (default
                        OTUs.fasta)
  -loadModel <inModel.pkl>
                        load OTU model from python pickle
  -saveModel <outModel.pkl>
                        save OTU model as python pickle
  -ptl <0.9>            percentile trimming length
  -pa <OTU_assignments.txt>
                        print detailed assignments
  -pp <OTU_profiles.txt>
                        print profile
  -t <0.75>             threshold score .. advanced usage

