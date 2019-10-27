#!/usr/bin/env python
# coding: utf-8

### Parses mmCIF
### Issue 1: molecule name (alternative: can be set as "nAChR" by default)
import os
import re
from Bio.PDB import *
from Bio.PDB import MMCIF2Dict

def numbersInPdbxDescription(s): #checks if string contains numbers, and returns the number (or 1, if none found)
    #right now is not necessary, but can help to handle other receptor types such as a4b2 etc.
    if any(i.isdigit() for i in s) == True:
        return (str(int(re.search(r'\d+', s).group())))
    else:
        return ('1')

def summaryLengthOfMmcifMolecule(structure):
    summary_length = 0
    for chain in structure.get_chains():
        summary_length = summary_length + len([_ for _ in chain.get_residues() if is_aa(_)])
    return (summary_length)

def chainNameFromPdbxDescription(mmcif_dict, chain):
    #returns string with Chain_Type depending on the content of the values of _entity.pdbx_description key
    DictChainEntityId = dict(zip(mmcif_dict['_struct_asym.id'], mmcif_dict['_struct_asym.entity_id']))#dictionary pairing id and entity_id
    DictChainIdPdbxDescription = dict(zip(mmcif_dict['_entity.id'], mmcif_dict['_entity.pdbx_description']))#dictionary pairing entity.id and pdbx_description
    
    if chain.id in DictChainEntityId.keys():
        if DictChainEntityId[chain.id] in DictChainIdPdbxDescription.keys(): #DictChainIdPdbxDescription[DictChainEntityId[chain.id]] 
            if re.search("alpha", DictChainIdPdbxDescription[DictChainEntityId[chain.id]], re.IGNORECASE):
                return ('Alpha' + numbersInPdbxDescription(DictChainIdPdbxDescription[DictChainEntityId[chain.id]]))
            elif re.search("beta", DictChainIdPdbxDescription[DictChainEntityId[chain.id]], re.IGNORECASE):
                return ('Beta' + numbersInPdbxDescription(DictChainIdPdbxDescription[DictChainEntityId[chain.id]]))
            elif re.search("delta", DictChainIdPdbxDescription[DictChainEntityId[chain.id]], re.IGNORECASE):
                return ('Delta')
            elif re.search("gamma", DictChainIdPdbxDescription[DictChainEntityId[chain.id]], re.IGNORECASE):
                return ('Gamma')
            elif re.search("epsilon", DictChainIdPdbxDescription[DictChainEntityId[chain.id]], re.IGNORECASE):
                return ('Epsilon')
            else:
                print("Error!")
                raise

def iterateThroughMmcifFiles(): #iterates through all PDB files in the current directory, and prints the retrieved data
    directory = os.fsencode(os.getcwd())
    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith(".cif"):
            parser = MMCIFParser()
            #some pdb produce warning such as 'discontinuous chain' at certain lines of mmCIF file
            #but still handles the sequence
            #for those mentioned in the database, it works fine without any warnings
            structure = parser.get_structure(filename, filename)
            mmcif_dict = MMCIF2Dict.MMCIF2Dict(filename)
            
            #Below is a code which retrieves the data for different fields in Molecule entry (based on sample XML database file)
            Database_Name = 'PDB'
            Database_ID = mmcif_dict['_entry.id']
            Name = 'nicotinic acetylcholine receptor (neuromuscular type)' #for now; there is no uniform protein name in mmCIF files (same as with PDB files)
            if '_entity_src_nat.pdbx_organism_scientific' in mmcif_dict:
                Organism = mmcif_dict['_entity_src_nat.pdbx_organism_scientific'][0]
            elif '_entity_src_gen.pdbx_host_org_scientific_name' in mmcif_dict:
                Organism = mmcif_dict['_entity_src_gen.pdbx_host_org_scientific_name'][0]
            elif '_entity_src_gen.pdbx_gene_src_scientific_name' in mmcif_dict:
                Organism = mmcif_dict['_entity_src_gen.pdbx_gene_src_scientific_name'][0]
            
            if '_entity_src_nat.tissue' in mmcif_dict:
                if '_entity_src_nat.pdbx_organ' in mmcif_dict:
                    Tissue = (mmcif_dict['_entity_src_nat.pdbx_organ'][0] + ', ' + mmcif_dict['_entity_src_nat.tissue'][0])
                else:
                    Tissue = (mmcif_dict['_entity_src_nat.tissue'][0])
            else:
                Tissue = 'not available'
            
            Length = summaryLengthOfMmcifMolecule(structure)
            Source = mmcif_dict['_exptl.method'].lower()
            if (Source).lower() == ('Electron Microscopy').lower():  #add conditions for X-ray crystallography and NMR
                #if it returns '?' instead of resolution, it means that in mmCIF file the resolution was not specified correctly ('?' instead of number in the corresponding field)
                Resolution = mmcif_dict['_em_3d_reconstruction.resolution'] # key depends on the method; it makes sense to detect the method first
                
            print(Database_Name, Database_ID, Name, Organism, Tissue, Length, Resolution, Source, '\n')
            
            for model in structure:
                #none of the tested mmCIF files has multiple models
                for chain in model:
                    Chain_ID = chain.id
                    Chain_Type = chainNameFromPdbxDescription(mmcif_dict, chain)
                    print(Chain_ID, Chain_Type, '\n')
                    #The following retrieves AA_Code (residue.get_resname) and Residue_Position (residue.id)
                    #Residue_Position_Protein would be the same
                    for residue in chain: #some pdb produce warning such as 'discontinuous chain' at certain lines of mmCIF file
                    #but still handles the sequence
                    #for those mentioned in the database, it works fine without any warnings
                        if is_aa(residue):
                            print(residue.id[1], residue.get_resname()) #residue ID is a tuple with three elements
                    print('\n')
        else:
            continue

iterateThroughMmcifFiles()


# In[ ]:





# In[ ]:




