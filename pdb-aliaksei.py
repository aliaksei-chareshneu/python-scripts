#!/usr/bin/env python
# coding: utf-8

### Parses PDB using Bio.PDB module
import os
import re
from Bio.PDB import *

def numbersInCOMPNDMolecule(s): #checks if string contains numbers, and returns the number (or 1, if none found)
    #right now is not necessary, but can help to handle other receptor types such as a4b2 etc.
    if any(i.isdigit() for i in s) == True:
        return (str(int(re.search(r'\d+', s).group())))
    else:
        return ('1')

def summaryLengthOfPDBMolecule(structure):
    summary_length = 0
    for chain in structure.get_chains():
        summary_length = summary_length + len([_ for _ in chain.get_residues() if is_aa(_)])
    return (summary_length)

def chainNameDefinition(chain_id): #returns string with Chain_Type depending on the Chain_ID
    #cannot unambiguosly classify epsilon chain, and cannot deal with non-muscle type structures
    #however, there is no epsilon chains in the tested PDB files
    if chain_id == 'A':
        return ('Alpha')
    elif chain_id == 'B':
        return ('Beta')
    elif chain_id == 'C':
        return ('Delta')
    elif chain_id == 'D':
        return ('Alpha')
    elif chain_id == 'E':
        return ('Gamma')
    else:
        print("Error!")
        raise

def chainNameFromCompound(structure, chain): #Alternative of the above function
    #returns string with Chain_Type depending on the content of COMPND fields of a PDB file
    #seems more reliable and can retrieve chain type from more different types of PDB files tested
    for f in structure.header['compound']:
        if re.search(chain.id, structure.header['compound'][f]['chain'], re.IGNORECASE):
            if re.search("alpha", structure.header['compound'][f]['molecule'], re.IGNORECASE):
                return ('Alpha' + numbersInCOMPNDMolecule(structure.header['compound'][f]['molecule']))
            elif re.search("beta", structure.header['compound'][f]['molecule'], re.IGNORECASE):
                return ('Beta' + numbersInCOMPNDMolecule(structure.header['compound'][f]['molecule']))
            elif re.search("delta", structure.header['compound'][f]['molecule'], re.IGNORECASE):
                return ('Delta')
            elif re.search("gamma", structure.header['compound'][f]['molecule'], re.IGNORECASE):
                return ('Gamma')
            elif re.search("epsilon", structure.header['compound'][f]['molecule'], re.IGNORECASE):
                return ('Epsilon')
            else:
                print("Error!")
                raise
        
def iterateThroughPDBFiles(): #iterates through all PDB files in the current directory, and prints the retrieved data
    directory = os.fsencode(os.getcwd())
    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith(".pdb"):
            parser = PDBParser()
            structure = parser.get_structure(filename, filename) # the 1st argument is some custom "user-defined name", the 2nd - actual name of a PDB file 
            #Below is a code which retrieves the data for different fields in Molecule entry (based on sample XML database file)
            Database_Name = 'PDB'
            Database_ID = (filename.split('.pdb', 1)[0].lower().rstrip())
            Name = 'nicotinic acetylcholine receptor (neuromuscular type)' #for now; there is no uniform protein name in PDB files
            Organism = (structure.header['source']['1']['organism_scientific'])
            #Tissue: not always present and not in a consistent way
            #therefore we may include 'organ' to make it more clear (but still not uniform)
            #In tested PDB, tissue and organ data are listed five times, probably for each chain
            if 'source' in structure.header:
                if ('organ' or 'tissue') in structure.header['source']['1']:
                    Tissue = (structure.header['source']['1']['organ'] + ", " + structure.header['source']['1']['tissue'])
            Length = summaryLengthOfPDBMolecule(structure)
            Resolution = structure.header['resolution']
            Source = structure.header['structure_method']
            
            print(Database_Name, Database_ID, Name, Organism, Tissue, Length, Resolution, Source, '\n')
            
            for model in structure:
                #none of the tested PDB files has multiple models
                for chain in model:
                    #print(chain.id, chainNameDefinition(chain.id))# alternative version
                    Chain_ID = chain.id
                    Chain_Type = chainNameFromCompound(structure, chain)
                    print(Chain_ID, Chain_Type, '\n')
                    #The following retrieves AA_Code (residue.get_resname) and Residue_Position (residue.id)
                    #Residue_Position_Protein would be the same
                    for residue in chain: #some pdb produce warning such as 'discontinuous chain' at certain lines of PDB file
                    #but still handles the sequence
                    #for those mentioned in the database, it works fine without any warnings
                        if is_aa(residue):
                            print(residue.id[1], residue.get_resname()) #residue ID is a tuple with three elements
                    print('\n')
        else:
            continue
            
iterateThroughPDBFiles()