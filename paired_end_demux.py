#!/usr/bin/python

#################################
# E. Roberson                   #
# Deindexes Illumina sequencing #
#################################

###########
# Imports #
###########
import argparse
from collections import namedtuple
from datetime import datetime
import sys

####################
# Version and name #
####################
SCRIPT_PATH = sys.argv[0]
SCRIPT_NAME = SCRIPT_PATH.split('/')[-1].split('\\')[-1]
VERSION = '1.0.0'

#######################
# Constants & globals #
#######################
SequenceInterpretation = namedtuple( 'SequenceInterpretation', 'totalBp adjustBp extensionName' )

dnaMismatchDict = {'A':('T','C','G','N'), 'T':('A','C','G','N'), 'C':('A','G','T','N'), 'G':('A','C','T','N')} # Dictionary of mutations for index mismatch

EXT = 'fqs' # fqs used rather than FASTQ to signify +33 scale output

############
# Counters #
############
readCount = 0
ambigCount = 0
matchCount = 0
highNsCount = 0

########
# Fxns #
########
def errMsg( s ):
	print "Error: " + s
	sys.exit(1)

def nonZeroDivisor( count ):
	return count if count > 0 else 1

def isIntString( s ):
	try:
		int(s)
	except ValueError:
		return False
	return True

def boolString( boolVal ):
	if boolVal == True:
		return "True"
	else:
		return "False"

def getSequenceAmount( totalSeqBp ):
	if totalSeqBp > 1E15:
		ext = "Pbp"
		adjBp = totalSeqBp / 1E15
	elif totalSeqBp > 1E12:
		ext = "Tbp"
		adjBp = totalSeqBp / 1E12
	elif totalSeqBp > 1E9:
		ext= "Gbp"
		adjBp = totalSeqBp / 1E9
	elif totalSeqBp > 1E6:
		ext = "Mbp"
		adjBp = totalSeqBp / 1E6
	elif totalSeqBp > 1E3:
		ext = "Kbp"
		adjBp = totalSeqBp / 1E3
	else:
		ext = "bp"
		adjBp = totalSeqBp
	return SequenceInterpretation( totalSeqBp, adjBp, ext )

def getIndexSize( filename, readLength=None ):
	try:
		INFH = open(filename, 'r')
	except:
		errMsg( "Error when trying to open file [%s] to determine index length" % (filename) )
	
	# FASTQ formatted input
	# Assumes all indexes are equal length!
	line = INFH.readline() # identifier
	line = INFH.readline().rstrip("\r\n") # sequence
	INFH.close()
	
	if readLength == None:
		return len(line.rstrip("\r\n"))
	else:
		return len(line.rstrip("\r\n")) - readLength

def getIdSep( filename ):
	try:
		INFH = open(filename, 'r')
	except:
		errMsg( "Error when trying to open file [%s] to determine index length" % (filename) )
	
	splitString = ""
	
	line = INFH.readline() # identifier
	INFH.close()
	
	line = line.rstrip('\r\n')
	
	if len(line.split('#',1)) == 2:
		return "#"
	elif len(line.split(' ', 1)) == 2:
		return " "
	else:
		errMsg("No recognized split character in ID field.")

def getQualScale( filename ):
	try:
		INFH = open( filename, 'r' )
	except:
		errMsg( "Error when trying to open file [%s] to determine quality scale" % (filename) )
	
	# phred		 0 to 93, ASCII 33 to 126; offset 33; phred scale
	# sol1.0	-5 to 62, ASCII 59 to 126; offset 64; log-odds
	# sol1.3+	 0 to 62, ASCII 64 to 126; offset 64; phred scale
	# Qsanger = -10 * log10(p)
	# Qsol1.0 = -10 * log10( p / 1-p )
	# 33 - 58 == phred
	# 59 - 63 == sol1.0
	# sol1.3 if no evidence for either
	
	maxReads = 5000
	readCount = 0
	minQual = 10000
	localIlmn1p0Min = 59
	localIlmn1p3Min = 64
	localPhredMin = 33
	
	buffer = ''
	
	while INFH and readCount < maxReads:
		buffer = INFH.readline() # @ identifier
		buffer = INFH.readline() # seq
		buffer = INFH.readline() # + identifier
		buffer = INFH.readline().rstrip("\r\n") # quality
		
		if len( buffer ) == 0:
			break
		
		readCount += 1
		
		for qScore in buffer:
			if ord(qScore) < minQual:
				minQual = ord(qScore)
				
			if minQual < localIlmn1p0Min:
				readCounts = maxReads
				break
	
	INFH.close()
	
	if minQual < localIlmn1p0Min:
		return "phred"
	elif minQual < localIlmn1p3Min:
		return "illumina1.0"
	else:
		return "illumina1.3"
	
def buildMismatches( seq, mismatch, mutationDictionary ):
	#########################
	# Clean up the sequence #
	#########################
	seq = seq.upper().rstrip("\r\n")
	
	#########################
	# Build with mismatches #
	#########################
	mismatchIndexDict= {seq:0}

	if mismatch == True:
		for base in mutationDictionary[seq[-1]]:
			newIndex = seq[:-1] + base
			
			# check not in dictionary
			if newIndex not in mismatchIndexDict:
				mismatchIndexDict[newIndex] = 0
	
	return mismatchIndexDict

#####################
# class Definitions #
#####################
class Index:	
	def __init__( self, index, basename, fileExtension ):
		self.index = index.upper().rstrip("\r\n")
		self.matches = 0
		self.filename = []
		self.filehandle = []
		self.indexMatches = None
		
		self.filename.append( "%s_%s_1.%s" % ( basename, index, fileExtension ) )
		self.filename.append( "%s_%s_2.%s" % ( basename, index, fileExtension ) )
		
		for fname in self.filename:
			try:
				self.filehandle.append( open( fname, 'w' ) )
			except:
				errMsg( "Couldn't open output file [%s]" % (fname) )
	
	def write( self, fhIndex, s ):
		try:
			self.filehandle[fhIndex].write( s )
		except:
			errMsg( "Error writing to filehandle %s of %s open filehandles" % (fhIndex, len(self.filehandle) ) )
	
	def __del__( self ):
		for fhandle in self.filehandle:
			fhandle.close()
			
		if self.indexMatches != None:
			print ("\nIndex %s: %s matches" % (self.index, self.matches))
			for tags in self.indexMatches.keys():
				print ("\t%s %s" % (tags, self.indexMatches[tags]))

#############
# arg parse #
#############
parser = argparse.ArgumentParser(prog=SCRIPT_NAME, epilog="%s v%s" % (SCRIPT_NAME, VERSION))

parser.add_argument('fastq1', help="Name of read 1 fastq")
parser.add_argument('fastq2', help="Name of read 2 fastq")
parser.add_argument('indexes', help="Comma-separated list of indexes")
parser.add_argument('readLength', help="Number of cycles per end, e.g. 100", type=int)
parser.add_argument('outBase', help="Basename of output files")
parser.add_argument('--indexRead', help="If index read is separate file specify here", required=False, default=None)
parser.add_argument('--nPercentCutoff', help="Adjust the N-base percent to trigger a 'bad' read", required=False, default=0.85, type=float)
parser.add_argument('--noIndexMismatch', help="Disallows terminal mismatch of index if active", required=False, dest='indexMismatch', default=True, action="store_false")

args = parser.parse_args()

if args.nPercentCutoff < 0.0 or args.nPercentCutoff > 1.0:
	errMsg( "N content cutoff must be between 0.0 and 1.0!" )

INDEXES = [val.upper().rstrip("\n\r") for val in args.indexes.split(',')] 

####################
# Set index length #
####################
if args.indexRead != None:
	INDEXLENGTH = getIndexSize( args.indexRead )
else:
	INDEXLENGTH = getIndexSize( args.fastq1, args.readLength )

#####################
# Set quality scale #
#####################
QUALITYSCALE = getQualScale( args.fastq1 )
if QUALITYSCALE != 'phred':
	errMsg( "Input quality scale appears to be %s. Please convert to Phred (+33) scale before proceeding." % (QUALITYSCALE) )

#########################
# print settings to log #
#########################
print "%s v%s" % (SCRIPT_PATH, VERSION)
print "Options set"
print "==========="
print "Input-scale: %s" % (QUALITYSCALE)
print "Read length: %s" % (args.readLength)
print "Index length: %s" % (INDEXLENGTH)
print "Indexes: %s" % (', '.join(INDEXES))
print "Mismatches allowed: %s" % (boolString(args.indexMismatch))
print "Percent N content cutoff: %s" % (args.nPercentCutoff)
print "Input files: %s & %s" % (args.fastq1, args.fastq2)
if args.indexRead != None:
	print "Index read: %s" % (args.indexRead)
print "Output basename: %s" % (args.outBase)
print
sys.stdout.flush() # must flush to get early feedback

#######################################################
# Set custom identifier string return based on format #
#######################################################
IDSEPSTRING = getIdSep( args.fastq1 )

if IDSEPSTRING == "#":
	def identifierString( readId, indexSequence ):
		idVals = readId.split( "#", 1 )
		return '%s#%s/%s' % (idVals[0], indexSequence, idVals[1][-1])
elif IDSEPSTRING == " ":
	def identifierString( readId, indexSequence ):
		idVals = readId.split(' ', 1)
		pairYqc = idVals[1].split(":")
		return '%s %s:%s:%s:%s' % (idVals[0], pairYqc[0], pairYqc[1], pairYqc[2], indexSequence)
else:
	errMsg( "Identifier string [%s] not recognized" % (IDSEPSTRING) )

#############################
# Build index dictionary    #
# Populate information      #
# Build subindexes for each #
#############################
indexDict = {}
indexList = []
for i in range( len(INDEXES) ):
	INDEXES[i] = INDEXES[i][:INDEXLENGTH]
	
	##################################################
	# Warn and ignore if index listed more than once #
	##################################################
	if INDEXES[i] in indexDict:
		print "%s was already listed - ignoring" % (INDEXES[i])
		continue
	
	######################
	# Populate the index #
	######################
	indexDict[INDEXES[i]] = Index( INDEXES[i], args.outBase, EXT )
	
	#################
	# Build matches #
	#################
	indexDict[INDEXES[i]].indexMatches = buildMismatches( INDEXES[i], args.indexMismatch, dnaMismatchDict )

###########################
# Good sequence, no index #
###########################
indexDict['noindex'] = Index('noindex', args.outBase, EXT )

indexDict['noindex'].filename.append( "%s_%s_index.%s" % ( args.outBase, 'noindex', EXT ) )
try:
	indexDict['noindex'].filehandle.append( open( indexDict['noindex'].filename[-1], 'w' ) )
except:
	errMsg( "Can't make index file for unindexed reads" )

#########################
# Poor quality sequence #
#########################
indexDict['junk'] = Index('junk', args.outBase, EXT )

indexDict['junk'].filename.append( "%s_%s_index.%s" % ( args.outBase, 'junk', EXT ) )
try:
	indexDict['junk'].filehandle.append( open( indexDict['junk'].filename[-1], 'w' ) )
except:
	errMsg( "Can't make index file for junk reads" )

###################
# start the clock #
###################
analysisStartTime = datetime.now()

###############
# Open input1 #
###############
try:
	inputFile1 = open(args.fastq1, 'r')
except:
	errMsg( "Could not open input file %s" % (args.fastq1) )

###############
# Open input2 #
###############
try:
	inputFile2 = open(args.fastq2, 'r')
except:
	errMsg( "Could not open input file %s" % (args.fastq2) )

if args.indexRead != None:
	###################
	# Open index file #
	###################
	try:
		indexFile = open(args.indexRead, 'r')
	except:
		errMsg( "Could not open index file %s" % (args.indexRead) )

######################
# Parse line by line #
######################
loopIteration = 1
while inputFile1 and inputFile2:
	#############
	# first end #
	#############
	currentId1 = inputFile1.readline().rstrip('\r\n')
	currentSeq1 = inputFile1.readline().rstrip('\r\n')
	repeatId1 = inputFile1.readline().rstrip('\r\n')
	currentQuality1 = inputFile1.readline().rstrip('\r\n')
	
	##############
	# second end #
	##############
	currentId2 = inputFile2.readline().rstrip('\r\n')
	currentSeq2 = inputFile2.readline().rstrip('\r\n')
	repeatId2 = inputFile2.readline().rstrip('\r\n')
	currentQuality2 = inputFile2.readline().rstrip('\r\n')
	
	if args.indexRead != None:
		#########
		# index #
		#########
		indexId = indexFile.readline().rstrip('\r\n')
		indexRead = indexFile.readline().rstrip('\r\n')
		indexRepeatId = indexFile.readline().rstrip('\r\n')
		indexQuality = indexFile.readline().rstrip('\r\n')
	else:
		indexId = currentId1
		indexRead = currentSeq1[-INDEXLENGTH:]
		currentSeq1 = currentSeq1[:-INDEXLENGTH]
		indexQuality = currentQuality1[-INDEXLENGTH:]
		currentQuality1 = currentQuality1[:-INDEXLENGTH]
		
	##########################
	# Check empty & id syncs #
	##########################
	if len(currentId1) == 0 and len(currentId2) == 0 and len(indexId) == 0:
		break
	elif len(currentId1) == 0 or len(currentId2) == 0 or len(indexId) == 0:
		errMsg( "Whoa! Something is empty but not all reads are in sync.\nRead 1: %s\nRead 2: %s\nIndex: %s\n" % (currentId1, currentId2, indexId) )
	elif currentId1.split(IDSEPSTRING, 1)[0] != currentId2.split(IDSEPSTRING, 1)[0] or currentId1.split(IDSEPSTRING, 1)[0] != indexId.split(IDSEPSTRING, 1)[0]:
		errMsg( "Reads out of sync!!!Read1: [%s]\nRead2: [%s]\nIndex Read: [%s]\n" % (currentId1, currentId2, indexId ))
	
	#############################
	# Now process the read data #
	#############################
	match = False
	
	#################################
	# Get N percent from seq length #
	#################################
	n_perc1 = float(currentSeq1.count('N')) / float(len(currentSeq1))
	n_perc2 = float(currentSeq2.count('N')) / float(len(currentSeq2))
	
	if n_perc1 >= args.nPercentCutoff and n_perc2 >= args.nPercentCutoff:
		indexDict['junk'].write(0, "%s\n%s\n+\n%s\n" % (currentId1, currentSeq1, currentQuality1))
		indexDict['junk'].write(1, "%s\n%s\n+\n%s\n" % (currentId2, currentSeq2, currentQuality2))
		indexDict['junk'].write(2, "%s\n%s\n+\n%s\n" % (indexId, indexRead, indexQuality))
		
		highNsCount += 1
		match = True
	
	else:
		##########################
		# Search for index match #
		##########################
		for indexSeq in INDEXES:
			if match == True:
				break
			
			##############################
			# For each possible match of #
			# the primary index          #
			##############################
			for indexTag in indexDict[indexSeq].indexMatches.keys():
				if indexTag == indexRead:
					###################################
					# Write line to proper index file #
					###################################
					indexDict[indexSeq].write(0,"%s\n%s\n+\n%s\n" % ( identifierString(currentId1, indexTag), currentSeq1, currentQuality1))
					
					indexDict[indexSeq].write(1,"%s\n%s\n+\n%s\n" % ( identifierString(currentId2, indexTag), currentSeq2, currentQuality2))
					
					####################
					# Increment counts #
					####################
					indexDict[indexSeq].indexMatches[indexTag] += 1
					indexDict[indexSeq].matches += 1
					matchCount += 1
					
					##########################
					# Confirm index match    #
					# Short-circuit the loop #
					##########################
					match = True
					break
	
	if match != True:
		##########################################
		# If no match it is by default noindex #
		##########################################
		# Write to catch-all file #
		###########################
		indexDict['noindex'].write(0, "%s\n%s\n+\n%s\n" % (currentId1, currentSeq1, currentQuality1))
		indexDict['noindex'].write(1, "%s\n%s\n+\n%s\n" % (currentId2, currentSeq2, currentQuality2))
		indexDict['noindex'].write(2, "%s\n%s\n+\n%s\n" % (indexId, indexRead, indexQuality))
		
		ambigCount += 1
	
	readCount += 1
	
inputFile1.close()	
inputFile2.close()

if args.indexRead != None:
	indexFile.close()
	
#############################
# Print summary information #
#############################
nominalSequence = getSequenceAmount( readCount * args.readLength * 2 )
indexedSequence = getSequenceAmount( matchCount * args.readLength * 2 )

print """
Deindexing completed in %s seconds

###########
# Summary #
###########
""" % (datetime.now() - analysisStartTime)
# Total X end reads processed
print "%s paired-end reads processed" % (readCount)

# Total matching specified indexes
print "%s/%s (%.1f%%) reads matched specified indexes" % (matchCount, readCount, float(matchCount)/float( nonZeroDivisor(readCount) )*100.0)

print "%s reads (%.1f%%) ignored due to high 'N' content" % (highNsCount, float(highNsCount)/float( nonZeroDivisor(readCount) )*100.0)

print "%.1f %s total sequence" % (nominalSequence.adjustBp, nominalSequence.extensionName)
print "%.1f %s indexed sequence" % (indexedSequence.adjustBp, indexedSequence.extensionName)
