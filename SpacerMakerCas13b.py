#purpose of script is to parse phage genbank and rank the top protospacer choices given to rules published by Smargon et al.
from subprocess import Popen, PIPE
import numpy as np
from Bio import SeqIO

def PFSRules(seq):
    rule1, rule2, = False, False
    if seq[0] in ['A','U','G']:
        rule1 = True
    if 'A' in seq[-2:]:
        rule2 = True
    if rule1 and rule2:
        return True
    else:
        return False

def runRNAplfold(seq, idx):
    seq = seq[1:len(seq)-2]
    echo = Popen(['echo', seq], stdout=PIPE)
    RNAplfold = Popen(['RNAplfold','-W', '240', '-L', '180', '-u', '15'], stdin=echo.stdout, stdout=PIPE)
    echo.stdout.close()
    out, err = RNAplfold.communicate()
    data = np.genfromtxt(fname="plfold_lunp", delimiter="\t", skip_header=2,) #filling_values=1)  # change filling_values as req'd to fill in missing values
    results = []
    for i in idx:
        try:
            results.append([data[i][-1], i]) #last 15 nt window of transcript , need to change to last 15nt of potential protopacer
        except IndexError:
            pass
    return results

def parseGenbank(filename):
    transcripts = []
    utr_lengh = 10
    for gb_record in SeqIO.parse(open(filename,"r"), "genbank"):
        for feature in gb_record.features:
            if feature.type == 'CDS':
                start, end = feature.location.start, feature.location.end
                if feature.location.strand == 1:
                    mrna = gb_record.seq[start-utr_lengh:end+utr_lengh].transcribe()
                else:
                    mrna = gb_record.seq[start-utr_lengh:end+utr_lengh].reverse_complement().transcribe()
                transcripts.append(Transcript(mrna, feature.location))
    return transcripts

def findProtospacers(transcripts):
    for mrna in transcripts:
        windows = [(mrna.seq[i:i+34],i) for i in range(len(mrna.seq)-33)] #filter to abide by both rules
        protospacers = list(filter(lambda x: PFSRules(x[0]), windows))
        idx = [x[1]+34 for x in protospacers]
        results = runRNAplfold(str(mrna.seq),idx)
        final = list(zip(results,protospacers))
        a = sorted(final, key=lambda x: float(x[0][0]), reverse=True)
        for i in a:
            print('{}, {}, {}, {}'.format(i[0][0], i[0][1], i[1][0].back_transcribe(), i[1][1]))

class Transcript:
    def __init__(self, seq, location):
        self.seq = seq
        self.location = location

filename = 'JBD30.gb'
transcripts = parseGenbank(filename)
findProtospacers(transcripts)
#for mrna in transcripts:
#    runRNAplfold(str(mrna.seq))
#seq = 'AUGAUGAUGAGUAGGAUAUGAAGCUAGACAUAA'
#print(PFSRules(seq))
#runRNAplfold(seq)
