#!/home/vnkoparde/opt/Python-2.7.2/bin/python

class Seq:
	def __init__(self,seqstr,seqid,wordLength,maxLength,seqstr2="",maxLength2=0):
		if seqstr2 != "":
			self.pe=1	#	is paired end
		else:
			self.pe=0	#	is not paired end
		length=len(seqstr)
		if length>maxLength:
			self.seqstr=seqstr[:maxLength].upper()
			self.len=maxLength
		else:
			self.seqstr=seqstr.upper()
			self.len=length
		self.freq=1
		self.seqids=[]
		self.seqids.append(seqid)
		self.worddict=dict()
		for word in [seqstr[i:i+wordLength] for i in range(self.len-wordLength+1)]:
			if word in self.worddict:
				self.worddict[word]+=1
			else:
				self.worddict[word]=1
		if self.pe==0:
			self.seqstr2=seqstr2
			self.len2=0
			self.worddict2=dict()
		else:
			length2=len(seqstr2)
			if length2>maxLength2:
				self.seqstr2=seqstr2[:maxLength2].upper()
				self.len2=maxLength2
			else:
				self.seqstr2=seqstr2.upper()
				self.len2=length2
			self.worddict2=dict()
			for word in [seqstr[i:i+wordLength] for i in range(self.len2-wordLength+1)]:
				if word in self.worddict2:
					self.worddict2[word]+=1
				else:
					self.worddict2[word]=1
	def __eq__(self,other):
		return (self.len+self.len2) == (other.len+other.len2)
	def __ne__(self,other):
		return (self.len+self.len2) != (other.len+other.len2)
	def __ge__(self,other):
		return (self.len+self.len2) >= (other.len+other.len2)
	def __le__(self,other):
		return (self.len+self.len2) <= (other.len+other.len2)
	def __gt__(self,other):
		return (self.len+self.len2) > (other.len+other.len2)
	def __lt__(self,other):
		return (self.len+self.len2) < (other.len+other.len2)
	def append(self,seqid):
		self.freq+=1
		self.seqids.append(seqid)
		for word in self.worddict.iterkeys():
			self.worddict[word]+=1
		for word in self.worddict2.iterkeys():
			self.worddict2[word]+=1
	def prn(self):
		if self.pe==0:
			print "Len:",self.len," Freq:",self.freq," Seq:",self.seqstr
			print "Seqids:",self.seqids
		else:
			print self.freq,self.len,self.seqstr,self.len2,self.seqstr2
			for seqid in self.seqids:
				print seqid
	
