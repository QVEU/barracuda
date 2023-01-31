#
# Barracuda.py
#
# Usage: python Barracuda <coding sequence> [<start position> <length of window>]
#
#

#Imports
import sys
import numpy as np
import pandas as pd

#Functions
'''
Function: parseArgs()

Description:

Arguments:
    - sArgV: list, input vector from command line

'''

def parseArgs(sArgv):
    inputsequence=sArgv[1]
    if len(sArgv)>2:
        start=int(sArgv[2])
        if len(sArgv)>3:
            length=int(sArgv[3])
        else: length=len(inputsequence)
    else: 
    	start=1
    	length=0
    return(inputsequence,start,length)
'''
Function: readCodonsFromFiles()

Description:

Arguments:
	none

'''

def readCodonsFromFile():
    codons=pd.read_csv("/hpcdata/lvd_qve/QVEU_Code/bioinfo/CodonTable.csv",delimiter="\t")
    myDict=dict()
    myDict={X:[] for (X) in np.unique(codons.Letter)}
    for i in range(len(codons)):
        myDict[codons.Letter[i]].append(codons.Codon[i])
    return(codons,myDict)

'''
Function: readInput(inputsequence,myDict,start,lenMax)

Description:

Arguments: 
    - inputsequence : string, input string
    - myDict : dict, codons from input sequence
    - start: int, starting position (1 indexed)
    - lenMax: int, length of window from which to generate barcodes (potential DANGER ZONE)
'''
def readInput(inputsequence0,myDict,start=0,lenMax=30):
    left=""
    if start!=0:
        left=inputsequence0[0:(start-1)]
        inputsequence=inputsequence0[(start-1):len(inputsequence0)]
        print("Starting from position "+str(start))
    lenMod=len(inputsequence)%3
    if len(inputsequence)>lenMax:
        lenMod=len(inputsequence)-lenMax
        print("WARNING: Trimmed to "+str(lenMax)+" bases: -"+str(lenMod))
        inputsequence=inputsequence[:-lenMod]
        lenMod=len(inputsequence)%3
    if lenMod!=0:
        print("WARNING: Trimmed 3' End: -"+str(lenMod))
        inputsequence=inputsequence[:-lenMod]
    lenT=len(inputsequence)
    right=inputsequence0[(len(left)+lenT):len(inputsequence0)]
    outfile="barracuda_"+inputsequence+".fasta"
    print("Input Sequence: "+inputsequence)
    posDict=dict()
    for i in range(len(inputsequence.replace("U","T"))):    #For position in query, convert RNA to DNA in query and phase into codons
        if i%3==0:
            posDict[i]=myDict[codons.Letter[codons.Codon==inputsequence[i:i+3]].values[0]]  # Parse into codons.
    return(posDict,inputsequence,outfile,right,left) #return the dict of codons and the input sequence (trimmed).

'''
Function: generateVariants()

Description:

Arguments:
    - posDict : dict, dict of all positions
    - inputsequence : string, input nucleotide sequence to mutate
    - baseSeq : a prefix to put on all barcodes

'''
def generateVariants(posDict,inputsequence,baseSeq=""):
    baseSequence=baseSeq
    baseName=""
    sequenceList=[""]
    idList=[""]
    countDict=dict()
    prevCount=1
    for k in posDict.keys(): #for each codon position
        newSequences=sequenceList
        newNames=idList
        sequenceList=[]
        idList=[]
        count=0
        c=0
        for i in posDict[k]: #for each alternative codon
            for x in range(len(newSequences)):
                if i==inputsequence[k:(k+3)]:
                    jointSeq=baseSequence+newSequences[x]+i.lower()
                else:
                    jointSeq=baseSequence+newSequences[x]+i
                jointName=baseName+newNames[x]+str(c)
                sequenceList.append(jointSeq)
                idList.append(jointName)
                count+=1
            c+=1
        countDict.update({k:(count/prevCount)})
        prevCount=count
    print("Total Number of Barcoded Sequences: " +str(len(sequenceList)))
    return(sequenceList,idList,countDict)

'''
Function: plotVariantDist()

Description:

Arguments:
    - countDict : dict, counts of variants at each codon position

'''
def plotVariantDist(countDict):
    for K in countDict.keys():
        print(str(K)+":\t"+"*"*int(countDict[K]))
    print ("- Variants/Site ->")
'''
Function: plotVariantDist()

Description:

Arguments:
	- seqList,idList,right,left

'''
def writeOutBarcodedSeqs(seqList,idList,right,left,outfile="barcode.fasta"):
    with open(outfile,'w') as OUT:
        for INDEX in range(len(seqList)-1):
            OUT.write(">VarSeq_"+str(idList[INDEX])+"\n"+left.lower()+seqList[INDEX]+right.lower()+"\n")
    print("Wrote "+str(outfile))

if __name__=="__main__":
    if len(sys.argv)<2:
        print("\nUsage: python Barracuda <coding sequence> [<start position> <length of window>]")
        exit()
    inputsequence,start,length = parseArgs(sys.argv)
    codons,myDict = readCodonsFromFile()
    if length==0:
        print(length)
        posDict,trimInput,outfileF,right,left = readInput(inputsequence,myDict,start)
    else:
    	posDict,trimInput,outfileF,right,left = readInput(inputsequence,myDict,start,length)
    seqList,idList,countDict = generateVariants(posDict=posDict,inputsequence=trimInput)
    plotVariantDist(countDict)
    writeOutBarcodedSeqs(seqList,idList,right,left,outfile=outfileF)
