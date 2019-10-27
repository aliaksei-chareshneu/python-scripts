import os
from Bio import SeqIO
from Bio.Align.Applications import ClustalwCommandline

clustalw_exe = r"C:\ClustalW2\clustalw2.exe"
directory = ".\\"

pathToAlignmentInput = "alignment_input.fasta"

def runAlignment(alignment_input):
    #record_dict = SeqIO.index(alignment_input, "fasta")
    sequences_A = []
    sequences_B = []
    sequences_C = []
    sequences_E = []
    #D is absent because it is effectively A in PDB encoding

    for record in SeqIO.parse(alignment_input, "fasta"):
        if record.id.endswith('_A') or record.id.endswith('_D'):
            sequences_A.append(record)
        elif record.id.endswith('_B'):
            sequences_B.append(record)
        elif record.id.endswith('_C'):
            sequences_C.append(record)
        elif record.id.endswith('_E'):
            sequences_E.append(record)
    
    SeqIO.write(sequences_A, "sequences_A.fasta", "fasta")
    SeqIO.write(sequences_B, "sequences_B.fasta", "fasta")
    SeqIO.write(sequences_C, "sequences_C.fasta", "fasta")
    SeqIO.write(sequences_E, "sequences_E.fasta", "fasta")
    

    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.startswith("sequences_") and filename.endswith(".fasta"):        
            in_file = filename
            clustalw_cline = ClustalwCommandline(clustalw_exe, infile=in_file)
            assert os.path.isfile(clustalw_exe), "Clustal W executable missing"
            stdout, stderr = clustalw_cline()
            print(clustalw_cline)
        else:
            continue        

runAlignment(pathToAlignmentInput)