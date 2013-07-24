#!/home/vnkoparde/opt/Python-2.7.2/bin/python
import HTSeq
import BOTUX
import os,sys,argparse,itertools,math,pickle
#import gc
#gc.disable()
import profile

def percentile(N, P):
	n = int(round(P * len(N) + 0.5))
	return N[n-1]

def geomMean(a,b):
	return math.sqrt(a*b)

def getBestOtu(otus,seqObj,thresholdScore,thresholdScore2):
	bestScore=0
	bestScore2=0
	bests=0
	bests2=0
	bestOTU=-1
	for otuIndex in range(len(otus)):
		s,s2=otus[otuIndex].getSeqScore(seqObj,thresholdScore)
		score=0
		product=0
		if s2==0:
			score=s
			if score > thresholdScore and score > bestScore:
				bestScore=score
				bests=s
				bests2=s2
				bestOTU=otuIndex
		else:
			product=s*s2
			if product > thresholdScore2 and product > bestScore2:
				bestScore2=product
				bests=s
				bests2=s2
				bestOTU=otuIndex
	return bests,bests2,bestOTU
		

#if __name__ == "__main__":
def main():
	parser=argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,description=
	"""BOTUX: Bayesian-like OTU eXaminer
	
	Author: 	Vishal N Koparde, Ph. D. 
	Created:	130315
	Modified:	130330""",
	version="3.0")
	parser.add_argument('-fasta',help='input fasta file',dest='fasta',required=False,metavar="<in.fasta>")
	parser.add_argument('-fastq',help='input fastq file(read1)',dest='fastq1',required=False,metavar="<read1.fastq>")
	parser.add_argument('-2',help='input fastq file(read2)',dest='fastq2',required=False,metavar="<read2.fastq>")
	parser.add_argument('-wl',help='word lenth (default 8)',dest='wl',required=False,metavar="<8>",type=int,default=8)
	parser.add_argument('-project',help='project name .. prefix used to name output files',dest='proj',required=True,metavar="<P1>")
	parser.add_argument('-loadModel',help='load OTU model from python pickle',dest='inModel',required=False,metavar="<inModel.pkl>")
	parser.add_argument('-ptl',help='percentile trimming length',dest='ptl',required=False,metavar="<0.75>",default=0.75,type=float)
	parser.add_argument('-t',help='threshold score .. advanced usage',dest='thresholdScore',required=False,default=0.65,metavar="<0.65>",type=float)
	args=vars(parser.parse_args())

        args['otuFasta']=args['proj']+"_OTUs.fasta"
        args['prnDetailedAssignments']=args['proj']+"_assignments.txt"
        args['prnProfile']=args['proj']+"_profiles.txt"
        args['outModel']=args['proj']+"_model.pkl"
	args['thresholdScore2']=args['thresholdScore']*args['thresholdScore']
	
	fasta=0
	fastq=0
	
	if not args['fasta'] and not args['fastq1']:
		print "Fasta or Fastq file must be supplied!"
		sys.exit()
		
	if args['fasta']:
		fasta=1
		if (not os.path.exists(args['fasta'])):
			print args['fasta']+" file not found!"
			sys.exit()
	
	if args['fastq1']:
		fastq+=1
		if (not os.path.exists(args['fastq1'])):
			print args['fastq1']+" file not found!"
			sys.exit()
	
	if args['fastq2']:
		fastq+=1
		if (not os.path.exists(args['fastq2'])):
			print args['fastq2']+" file not found!"
			sys.exit()
	
	
	seqlengths=[]
	seqstrings=dict()
	seqstrings2=dict()
	seqObjDict=dict()
	seqObjList=[]
	wordLength=args['wl']
	
	if fasta==1:
		seqstrings = dict( (s.name, s) for s in HTSeq.FastaReader(args['fasta']) )
	
	if fastq==1:
		seqstrings = dict( (s.name, s) for s in HTSeq.FastqReader(args['fastq1']) )
	elif fastq==2:
		for read1,read2 in itertools.izip(HTSeq.FastqReader(args['fastq1']),HTSeq.FastqReader(args['fastq2'])):
			seqstrings[read1.name]=read1
			seqstrings2[read1.name]=read2
		
	if fastq==2:
		seqlengths2=[]
	for seqid in seqstrings.iterkeys():
		seqlengths.append(len(seqstrings[seqid]))
		if fastq==2:
			seqlengths2.append(len(seqstrings[seqid]))
			
	maxLength=percentile(sorted(seqlengths),args['ptl'])
	if fastq==2:
		maxLength2=percentile(sorted(seqlengths2),args['ptl'])
	
	p=[int(round((i+1)*.1*len(seqlengths))) for i in range(10)]
	
	nseq=len(seqlengths)
	avgSeqLength=int(sum(seqlengths)*1.0/nseq)
	print "Total sequences: "+str(nseq)
	print "Average sequence length: "+str(avgSeqLength)
	print "Sequence trimming length: "+str(maxLength)
		
	del seqlengths
	if fastq==2:
		del seqlengths2
	
	
	count=0
	for seqid in seqstrings.keys():
		count+=1
		if count in p:
			q=(p.index(count)+1)*10
			print str(q)+"% sequences loaded"
		seqstr=str(seqstrings[seqid])
		if fastq==2:
			seqstr2=str(seqstrings2[seqid])
			try:
				seqObjDict[seqstr+seqstr2].append(seqid)
			except KeyError:
				seqObjDict[seqstr+seqstr2]=BOTUX.Seq(seqstr,seqid,wordLength,maxLength,seqstr2,maxLength2)
		else:
			try:
				seqObjDict[seqstr].append(seqid)
			except KeyError:
				seqObjDict[seqstr]=BOTUX.Seq(seqstr,seqid,wordLength,maxLength)
	
	for seqid in seqObjDict.keys():
		seqObjList.append(seqObjDict[seqid])
	
	del seqstrings
	del seqObjDict
	
	seqObjList.sort()
	seqObjList.reverse()
	
	print str(len(seqObjList))+" sequence objects created!"
	
	#for s in seqObjList:
	#	s.prn()

	#sys.exit()
	
	otus=[]
	if args['inModel']:
		IM=open(args['inModel'],'rb')
		otus=pickle.load(IM)
		IM.close()
	p=[int(round((i+1)*.1*len(seqObjList))) for i in range(10)]
	
	count=0
	for s in seqObjList:
		count+=1
		if count in p:
			q=(p.index(count)+1)*10
			print str(q)+"% sequences processed .. "+str(len(otus))+" distinct OTUs created."
		bestScore,bestScore2,bestOTU=getBestOtu(otus,s,args['thresholdScore'],args['thresholdScore2'])
		#print "count:"+str(count)+" otus:"+str(len(otus))+" bestScore:"+str(bestScore)
		if bestOTU==-1:
			otus.append(BOTUX.Otu(s))
		else:
			otus[bestOTU].addSeqObj(s,bestScore,bestScore2)
				
	otus.sort()
	otus.reverse()
	notus=len(otus)
	rangenotus=range(notus)
	print str(notus)+" distinct OTUs created!"
	nRealOtus=0
	for o in otus:
		if o.freq>5:
			nRealOtus+=1
	print str(nRealOtus)+" OTUs have frequency > 5."
	
	
	otuFasta="OTUs.fasta"
	if args['otuFasta']:
		otuFasta=args['otuFasta']
	OUTFILE=open(otuFasta,"w")
	for i in rangenotus:
		otus[i].prn2fasta(OUTFILE,i+1)
	OUTFILE.close()
	
	print "OTUs written to file "+otuFasta+"."
	
	if args['outModel']:
		OM=open(args['outModel'],'wb')
		pickle.dump(otus,OM)
		OM.close()
		
	if args['prnDetailedAssignments']:
		OUTFILE=open(args['prnDetailedAssignments'],'w')
		for i in rangenotus:
			otus[i].prnDetailedAssignments(OUTFILE,i+1)
		OUTFILE.close()

	if args['prnProfile']:		
		OUTFILE=open(args['prnProfile'],'w')
		for i in rangenotus:
			otus[i].prnProfile(OUTFILE,i+1,nseq)
		OUTFILE.close()
	sys.exit()

main()
#profile.run("main()", "profile_tmp")
	

