#!/home/vnkoparde/opt/Python-2.7.2/bin/python
import Seq,HTSeq
import sys,math

class Otu:
    def __init__(self,seqObj):
	self.pe=seqObj.pe
	self.seed=seqObj.seqstr
	self.seed2=seqObj.seqstr2
	self.seedLen=seqObj.len
	self.seedLen2=seqObj.len2
	self.freq=seqObj.freq
	self.seqids=seqObj.seqids[:]
	self.worddict=seqObj.worddict.copy()
	self.worddict2=seqObj.worddict2.copy()
	self.totalReadLen=self.seedLen*self.freq
	self.totalReadLen2=self.seedLen2*self.freq
	if self.totalReadLen!=0:
	    self.avgReadLen=self.totalReadLen/len(self.seqids)
	else:
	    self.avgReadLen=0
	if self.totalReadLen2!=0:
	    self.avgReadLen2=self.totalReadLen2/len(self.seqids)
	else:
	    self.avgReadLen2=0
	self.scores=[-1]*self.freq
	if self.pe==1:
	    self.scores2=[-1]*self.freq
	else:
	    self.scores2=[]
	self.sumOfAllFreq=sum(self.worddict.values())
	self.sumOfAllFreq2=sum(self.worddict2.values())
    def getSeqScore(self,seqObj,thresholdScore):
        score=0
        score2=0
        for w,f in seqObj.worddict.iteritems():
            if w in self.worddict:
                score += (self.worddict[w] * f * 1.0 / self.sumOfAllFreq)
        score = score * self.seedLen / seqObj.len
	if score < thresholdScore:
	    return score,score2
	if self.pe==1:
	    for w,f in seqObj.worddict2.iteritems():
		if w in self.worddict2:
		    score2 += (self.worddict2[w] * f * 1.0 / self.sumOfAllFreq2)
        if score2 != 0: # to avoid divide by zero
            score2 = score2 * self.seedLen2 / seqObj.len2        
        return score,score2        
    def addSeqObj(self,seqObj,score,score2=0):
        self.freq+=seqObj.freq
        for seqid in seqObj.seqids:
            self.seqids.append(seqid)
            self.totalReadLen+=seqObj.len
            self.scores.append(score)
            if self.pe==1:
                self.totalReadLen2+=seqObj.len2
                self.scores2.append(score2)
	self.avgReadLen=self.totalReadLen*1.0/self.freq
	self.avgReadLen2=self.totalReadLen2*1.0/self.freq
        for c,f in seqObj.worddict.iteritems():
	    x=f*seqObj.freq
	    self.sumOfAllFreq+=x
            try:
                self.worddict[c]+=x
            except KeyError:
                self.worddict[c]=x
	if self.pe==1:
	    for c,f in seqObj.worddict2.iteritems():
		x=f*seqObj.freq
		self.sumOfAllFreq2+=x
                try:
		    self.worddict2[c]+=x
                except KeyError:
		    self.worddict2[c]=x

    def __eq__(self,other):
        return self.seedLen == other.seedLen and self.seedLen2 == other.seedLen2 and self.freq == other.freq
    def __ne__(self,other):
        return not self == other
    def __lt__(self,other):
        if self.seedLen == other.seedLen and self.seedLen2 == other.seedLen2:
            return self.freq < other.freq
        if self.seedLen2 == 0:
            return self.seedLen < other.seedLen
        else:
            return math.sqrt(self.seedLen * self.seedLen2) < math.sqrt(other.seedLen * other.seedLen2)
    def __le__(self,other):
        if self.seedLen == other.seedLen and self.seedLen2 == other.seedLen2:
            return self.freq <= other.freq
        if self.seedLen2 == 0:
            return self.seedLen <= other.seedLen
        else:
            return math.sqrt(self.seedLen * self.seedLen2) <= math.sqrt(other.seedLen * other.seedLen2)
    def __gt__(self,other):
        if self.seedLen == other.seedLen and self.seedLen2 == other.seedLen2:
            return self.freq > other.freq
        if self.seedLen2 == 0:
            return self.seedLen > other.seedLen
        else:
            return math.sqrt(self.seedLen * self.seedLen2) > math.sqrt(other.seedLen * other.seedLen2)
    def __ge__(self,other):
        if self.seedLen == other.seedLen and self.seedLen2 == other.seedLen2:
            return self.freq >= other.freq
        if self.seedLen2 == 0:
            return self.seedLen >= other.seedLen
        else:
            return math.sqrt(self.seedLen * self.seedLen2) >= math.sqrt(other.seedLen * other.seedLen2)
    def prn(self,details=0):
        print "Seed:",self.seed," len.:",self.seedLen 
        print "Freq:",self.freq," avg. seq. len.:",self.avgReadLen
        print "Seqids:",self.seqids
        print "Scores:",self.scores
        if details==1:
            for w,f in self.worddict.iteritems():
                print w,f
        if self.pe==1:
            print "Scores2:",self.scores2
            if details==1:
                for w,f in self.worddict2.iteritems():
                    print w,f
    def prn2fasta(self,FH,otuid):
        seqid="OTU_%05d|AvgLen=%d|Freq=%d"%(otuid,self.avgReadLen,self.freq)
        seq=HTSeq.Sequence(self.seed,seqid)
        seq.write_to_fasta_file(FH)
        if self.pe==1:
            seqid="OTU_%05d_read2|AvgLen=%d|Freq=%d"%(otuid,self.avgReadLen2,self.freq)
            seq=HTSeq.Sequence(self.seed2,seqid)
            seq.write_to_fasta_file(FH)
    def prnDetailedAssignments(self,FH,otuid):
        #FH.write("OTU_%05d\nSeed:%s\n"%(otuid,self.seed))
        #if self.pe==1:
        #    FH.write("Seed2:%s\n"%(self.seed2))
        #FH.write("SeedLength:%d\n"%(self.seedLen))
        #if self.pe==1:
        #    FH.write("SeedLength2:%s\n"%(self.seedLen2))
        #FH.write("Frequency:%d\n"%(self.freq))
        #FH.write("AverageSequenceLength:%d\n"%(self.avgReadLen))
        #if self.pe==1:
        #    FH.write("AverageSequenceLength2:%s\n"%(self.avgReadLen2))
	if self.pe==1:
	    for i in range(self.freq):
		FH.write("%s\tOTU_%05d\t%7.2f\t%7.2f\n"%(self.seqids[i],otuid,self.scores[i],self.scores2[i]))
	else:
	    for i in range(self.freq):
		FH.write("%s\tOTU_%05d\t%7.2f\n"%(self.seqids[i],otuid,self.scores[i]))
    def prnProfile(self,FH,otuid,nseq):
	perc=self.freq*100.0/nseq
	FH.write("OTU_%05d\t%d\t%6.2f\n"%(otuid,self.freq,perc))

