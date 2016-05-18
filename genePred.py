##Bacterial Gene Predictor
##Author: Martin Basterrechea
##Python3
##This program takes a fasta nucleotide file and outputs the list of possible ORF's with info in gff3 format.
##It includes: GC content in all frames, Shine-Dalgarno sequence presence, and codon bias.
##NOTE: There is no support in calculations for non-canonical bases, but they are tolerated in input files.


import sys
import re
import getopt
import math




#Python throws an exception when sending the output on a pipe, the following two lines solve this.
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL)

##GLOBAL

#Start & Stop codons to consider
STARTCODONS = ['ATG']
STOPCODONS = ['TAA','TAG','TGA']



#Output file name
OUTFILE = ''

#Default ORF allowed sizes
LONGESTORF = 3000
SHORTESTORF = 200

#Verbose?
VERB = False

#Remove all sequences with Shine-Dalgarno probability of 0?
FILTERSD = True



#List with all ORF
ORFList = []

#Start & Stop regexp
START = ")|(".join(STARTCODONS)
START = re.compile("((" + START + "))")
STOP = ")|(".join(STOPCODONS)
STOP = re.compile("((" + STOP + "))")

#Shine-Dalgarno weight matrix, taken from http://www.ics.uci.edu/~kibler/pubs/Metmbs02.pdf
SD = ({'A':0.29,'T':0.47,'G':0.11,'C':0.13},
      {'A':0.42,'T':0.00,'G':0.15,'C':0.43},
      {'A':0.81,'T':0.11,'G':0.00,'C':0.08},
      {'A':0.00,'T':0.00,'G':1.0,'C':0.00},
      {'A':0.00,'T':0.00,'G':1.0,'C':0.00},
      {'A':0.97,'T':0.00,'G':0.02,'C':0.01},
      {'A':0.23,'T':0.07,'G':0.66,'C':0.04},)

GCcontent = 0

aaTable = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
       "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
       "TAT":"Y", "TAC":"Y", "TAA":"*", "TAG":"*",
       "TGT":"C", "TGC":"C", "TGA":"*", "TGG":"W",
       "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
       "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
       "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
       "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
       "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
       "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
       "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
       "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
       "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
       "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
       "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
       "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",}

#Store codon count to calculate frequency
globalTrTable = {"TTT":0, "TTC":0, "TTA":0, "TTG":0,
       "TCT":0, "TCC":0, "TCA":0, "TCG":0,
       "TAT":0, "TAC":0, "TAA":0, "TAG":0,
       "TGT":0, "TGC":0, "TGA":0, "TGG":0,
       "CTT":0, "CTC":0, "CTA":0, "CTG":0,
       "CCT":0, "CCC":0, "CCA":0, "CCG":0,
       "CAT":0, "CAC":0, "CAA":0, "CAG":0,
       "CGT":0, "CGC":0, "CGA":0, "CGG":0,
       "ATT":0, "ATC":0, "ATA":0, "ATG":0,
       "ACT":0, "ACC":0, "ACA":0, "ACG":0,
       "AAT":0, "AAC":0, "AAA":0, "AAG":0,
       "AGT":0, "AGC":0, "AGA":0, "AGG":0,
       "GTT":0, "GTC":0, "GTA":0, "GTG":0,
       "GCT":0, "GCC":0, "GCA":0, "GCG":0,
       "GAT":0, "GAC":0, "GAA":0, "GAG":0,
       "GGT":0, "GGC":0, "GGA":0, "GGG":0,}

#Store global aminoacid frequency
globalAaCount = {}


class ORF:
  def __init__(self, seq, strand, b20, seqId,start,stop, seqIdLength):
    self.seq = seq #The sequence of the ORF itself
    self.seqId = seqId #Name of the sequence containing this ORF (from input file)
    self.seqIdLength = seqIdLength #Length of the sequence containing this ORF (from input file)
    self.start = start #Starting position (always in the + strand). This would be the stop codon if strand is -
    self.stop = stop
    self.length = len(seq)
    self.strand = strand # '+' or '-'
    self.b20 = b20 #20 bases upstream of the gene, to calculate SD presence
    self.SDscore = self.checkSD()
    ##SD filter is here:
    if FILTERSD and self.SDscore > 0.0:
      ORFList.append(self) #Won't append to ORFList if filtered, and stop calculating
    else:
      return
    self.GC1 = self.GC(1)
    self.GC2 = self.GC(2)
    self.GC3 = self.GC(3)
    self.trTable, self.aaCount = self.createTransTable()

  def GC(self,frame): #Return GC content in frame
    ##NOTE, the non-canonical bases here are counted as CG
    seq = self.seq[frame-1::3]
    seqTot = seq.replace('N','')
    seqNoAT = re.sub('[ATN]', '', seq)
    return len(seqNoAT)/len(seqTot)

  def checkSD(self): #Return max probability of SD in previous 20 bases, as suggested in http://www.ics.uci.edu/~kibler/pubs/Metmbs02.pdf
    score = 0
    for n in range(len(self.b20)-7):
      sub = self.b20[n:n+7]
      subscore = 1
      for l in range(len(sub)):
        if sub[l] in SD[l]:
          subscore *= SD[l][sub[l]]
        else:
          subscore *= 0
      score = max(subscore, score)
    return score


  def createTransTable(self):
    #Creates self codon table count and adds codon count to global as well.
    #Creates aminoacid count and updates global aminoacid count as well (to calculate codon frequency per amino acid).
    global globalTrTable
    trTable = globalTrTable.copy()
    aaCount = {}
    for n in range(0,len(self.seq),3):
      sub = self.seq[n:n+3]
      if sub in trTable: #Ignores non-canonical codons
        trTable[sub] += 1
        globalTrTable[sub] += 1
        if aaTable[sub] in globalAaCount:
          globalAaCount[aaTable[sub]] += 1
        else:
          globalAaCount[aaTable[sub]] = 1
        if aaTable[sub] in aaCount:
          aaCount[aaTable[sub]] += 1
        else:
          aaCount[aaTable[sub]] = 1
    return (trTable, aaCount)



  def scoreTr(self): #Score self Trtable against global tr table (compare frequencies of each codon)
    score = 0
    for key in self.trTable.keys():
      if aaTable[key] in self.aaCount:
        self.trTable[key] = self.trTable[key]/self.aaCount[aaTable[key]]
      else:
        self.trTable[key] = 0

    for key in self.trTable:
      score += (self.trTable[key]-globalTrTable[key])**2 #Score is sum of the squares of the differences
    self.TRscore = score


  def printToFile(self, filename):
    with open(filename,'a') as o:
      o.write("\t".join([self.seqId,'MartinsGenePred','ORF', str(self.start), str(self.stop),".",self.strand,".","; ".join(["GC1=" + str(self.GC1),"GC2=" + str(self.GC2),"GC3=" + str(self.GC3),"SDscore="+str(self.SDscore),"TRscore="+str(self.TRscore),"SeqLength="+str(self.seqIdLength)])])+'\n')

  def remove(self):
    ORFList.remove(self)


class Seq:
  def __init__(self,seqId,seq1):
    self.seqId = seqId
    self.seq1 = seq1
    self.seq2 = revComp(seq1) #Reverse complement

  def findORFs(self): #Find ORFs in self.
    ORFs = []

    #Strand +
    starts = []
    stops = []
    #Get location of all start & stop codons
    for m in re.finditer(START, self.seq1):
      starts.append(m.start())
    for m in re.finditer(STOP,self.seq1):
      stops.append(m.start())

    starts.sort()
    stops.sort()

    for sr in starts:
      done = False
      i = 0
      while not done and i< len(stops): #Did we find a suitable stop codon?
        sp = stops[i]
        if (sp - sr) >= LONGESTORF: #Longer than max, stop with this start codon.
          break
        if (sp - sr) >= SHORTESTORF-3 and (sp - sr) % 3 == 0: #-3 is because sp is the first base of the codon, not the third
          ORFs.append(ORF(self.seq1[sr:sp+3], '+', self.seq1[sr-20:sr],self.seqId,sr,sp,len(self.seq1)))
          done = True
        i += 1


    #Strand -
    starts = []
    stops = []
    for m in re.finditer(START, self.seq2):
      starts.append(m.start())
    for m in re.finditer(STOP,self.seq2):
      stops.append(m.start())
    starts.sort()
    stops.sort()
    for sr in starts:
      done = False
      i = 0
      while not done and i< len(stops):
        sp = stops[i]
        if (sp - sr) >= LONGESTORF:
          break
        if (sp - sr) >= SHORTESTORF-3 and (sp - sr) % 3 == 0:
          ORFs.append(ORF(self.seq2[sr:sp+3], '-',self.seq2[sr-20:sr],self.seqId,len(self.seq2)-sp,len(self.seq2)-sr,len(self.seq1)))
          done = True
        i += 1
    self.ORFs = ORFs


def revComp(seq): #Reverse complement
  seq2 = ''
  d = {'A':'T','T':'A','C':'G','G':'C','M':'K','K':'M','R':'Y','Y':'R','S':'S','W':'W','N':'N'}
  for i in range(len(seq)):
    seq2 = d[seq[i]] + seq2
  return seq2



def parsefile(f): #Return all sequences on a fasta file.
  seqName = ''
  sequence = ''
  seqs = []
  first = True
  for line in f:
    line = line.strip()
    #Look for seq name
    if line[0] == ">":
      m = re.search("(?<=\>).*$",line)

      if not first:
        seqs.append( Seq(seqName,sequence))
      first = False
      seqName = m.group(0)
      sequence = ''
    else:
      sequence += line.upper()
  seqs.append(Seq(seqName,sequence))
  return seqs

def getGC(seqs): #Calculate global GC content
  total = 0
  GC = 0
  for seq in seqs:
    seq1 = seq.seq1
    seqTot = seq1.replace('N','')
    seqNoAT = re.sub('[ATN]', '', seq1)
    total += len(seqTot)
    GC += len(seqNoAT)

  return GC/total

def updateTrTable(): #Turn global codon count into frequencies.
  for key in globalTrTable.keys():
    if aaTable[key] in globalAaCount:
      globalTrTable[key] = globalTrTable[key]/globalAaCount[aaTable[key]]
    else:
      globalTrTable[key] = 0


def analyse_seqs(seqs):
  oprint("Searching for ORFs in sequences...")
  ORFcount = 0
  for i in range(len(seqs)):
    seqs[i].findORFs()
    oprint("Sequence " + str(i) + " analysed, " + str(len(seqs[i].ORFs)) + ' ORFs found')
    ORFcount += len(seqs[i].ORFs)
  if FILTERSD:
    oprint("{0} ORFs found. {1} were filtered out because of 0 Shine-Dalgarno probability.".format(ORFcount, ORFcount-len(ORFList)))
  else:
    oprint("{0} ORFs found.".format(ORFcount))

  updateTrTable() #Turn global codon count into frequencies.
  oprint("Scoring ORFs and writing to file...")
  scoreAndWriteORF()
  i = 0
  for key in globalTrTable.keys():
    if i < len(globalTrTable.keys()):
      oprint("{0}:{1} ".format(key,math.floor(globalTrTable[key]*100)/100), '')
    if i % 4 == 3:
      oprint('')
    i += 1
  oprint('')


def scoreAndWriteORF():
  c = 0
  for orf in ORFList:
    if (c % (len(ORFList)/20) <= 1):
      oprint("{0}% completed.".format(math.floor(100*c/len(ORFList))),'\r')
    orf.scoreTr()
    orf.printToFile(OUTFILE)
    c += 1
  oprint("100% completed")


def oprint(pStr='', pEnd='\n'):
  if VERB:
    print(pStr, end=pEnd)



def main(argv):
  global VERB, OUTFILE, GCcontent, LONGESTORF, SHORTESTORF, FILTERSD
  inputfile = ''
  helpText = """###Bacterial Gene Predictor
###Author: Martin Basterrechea
###Usage: python3 blastParser.py [options] -i FILE
Options:
    -h : print this help
    -v : verbose mode
    -i [FILE] : path to input file (in fasta)
    -o [FILE] : name of the output file (gff format)
    -m [INT] : Shortest ORF threshold (default = 200)
    -n [INT] : Longest ORF threshold (default = 3000)
    -s : Don't filter ORFs with Shine-Dalgarno probability of 0"""
  try:
    opts, args = getopt.getopt(argv,"hvi:o:m:n:s")
  except getopt.GetoptError:
    print(helpText)
    quit()
  for opt, arg in opts:
    if opt == '-h':
      print(helpText)
      quit()
    elif opt == '-v':
      VERB = True
    elif opt == '-i':
      inputfile = arg
    elif opt == '-o':
      OUTFILE = arg
    elif opt == '-n':
      try:
        LONGESTORF = int(arg)
      except ValueError:
        print("Error: Longest ORF size must be integer")
        quit()
    elif opt == '-m':
      try:
        SHORTESTORF = int(arg)
      except ValueError:
        print("Error: Shortest ORF size must be integer")
        quit()
    elif opt == '-s':
      FILTERSD = False
  if OUTFILE == '':
    OUTFILE = inputfile.split(".")[0] + ".gff3"
  if inputfile == '':
    print(helpText)
  if SHORTESTORF > LONGESTORF:
    print("Error: Longest ORF size must be higher than shortest ORF size")
    quit()
  try:
    fasta = open(inputfile,'r')
  except FileNotFoundError:
    print("Error: File not found")
    quit()
  oprint("Reading File")
  seqs = parsefile(fasta)
  fasta.close()
  oprint(str(len(seqs)) + " sequences read.")
  GCcontent = getGC(seqs)
  oprint('Global GC content is: ' + str(GCcontent))
  o = open(OUTFILE,'w')
  o.close()
  analyse_seqs(seqs)


if __name__ == "__main__":
  main(sys.argv[1:])
