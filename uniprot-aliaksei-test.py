#!/usr/bin/env python
# coding: utf-8

### Iterate through one big UniProt XML file with many UniProt entries, parsing them

import os
import re
from Bio import SeqIO
import xml.etree.ElementTree as ET

def SigPepLenFromFeatures(record): #determines the length of signal peptide based on the content under 'features = signal peptide' tag
    for f in range(len(record.features)):
        if record.features[f].type == 'signal peptide':
            return(record.features[f].location.end)
        else:
            continue
    return('signal peptide is not available')

def numbersInCOMPNDMolecule(s): #checks if string contains numbers, and returns the number (or 1, if none found)
    #right now is not necessary, but can help to handle other receptor types such as a4b2 etc.
    if any(i.isdigit() for i in s) == True:
        return (str(int(re.search(r'\d+', s).group())))
    else:
        return ('')

def sigPepLen(uniprotFilename): #not suitable for multiple-entry XML #returns the length of the signal peptide in the corresponding UniProt entry
    uniprot_xml_namespace = "{http://uniprot.org/uniprot}"
    tree = ET.parse(uniprotFilename)
    get_feature = tree.findall(".//"+uniprot_xml_namespace+"feature[@type='signal peptide']")
    if get_feature:
        feature_location_position = []
        for eachitem in get_feature:
            for item in list(eachitem)[0]:
                feature_location_position.append(item.attrib.get('position'))
            return (int(feature_location_position[1]))    
    else:
        print("no signal peptide found")
        return None

def chainSearch(record): #returns tuple with Chain_ID and Chain_Type; works properly only with vertebrate typology (17 subunits)
    if "vertebrata" in map(str.lower, record.annotations['taxonomy']):
        if re.search("alpha", record.description, re.IGNORECASE):
            return ('A', 'Alpha' + numbersInCOMPNDMolecule(record.description))
        elif re.search("beta", record.description, re.IGNORECASE):
            return ('B', 'Beta' + numbersInCOMPNDMolecule(record.description))
        elif re.search("delta", record.description, re.IGNORECASE):
            return ('C', 'Delta' + numbersInCOMPNDMolecule(record.description))
        elif re.search("gamma", record.description, re.IGNORECASE):
            return ('E', 'Gamma')
        elif re.search("epsilon", record.description, re.IGNORECASE):
            return ('E', 'Epsilon')
        elif re.search("uncharacterized", record.description, re.IGNORECASE):
            return ('uncharacterized protein', 'not available')
        else:
            return ('not available', 'not available')
    else:
        return ('not available', 'not available')

def iterateThroughUniprotXMLFiles(): #iterates through all XML Uniprot files in the current directory, and prints the retrieved data
    directory = os.fsencode(os.getcwd())
    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith("format.xml"):
            continue
        elif filename.endswith(".xml"):
            #j = 0
            for record in SeqIO.parse(filename, 'uniprot-xml'):
                if re.search("toxin", record.description, re.IGNORECASE):
                    continue
                #Below is a code which retrieves the data for different fields in Molecule entry (based on sample XML database file)
                #j = j + 1
                Database_Name = 'UNIPROT'
                Database_ID = record.id
                Name = record.name
                Organism = (record.annotations['organism'].split('(', 1)[0].rstrip())
                #Tissue: for some entries it is present, but not in a consistent way
                Length = (record.annotations['sequence_length'])
                Chain_ID = chainSearch(record)[0]
                Chain_Type = chainSearch(record)[1]
                Uniprot_Protein_Name = record.description
                Signal_peptide_length = SigPepLenFromFeatures(record)

                i = 10 #variable for iteration through Uniprot sequence
                AA_Code = record.seq[i]
                Residue_Position = i + 1
                if str(Signal_peptide_length).isdigit(): #checks if the signal peptide length is available
                    Residue_Position_Protein = i + 1 + Signal_peptide_length
                else: #if not, it is a string 'signal peptide is not available', then Residue_Position_Protein equals Residue_Position
                    Residue_Position_Protein = Residue_Position

                #print(j, Database_Name, Database_ID, Name, Organism, Length, Chain_ID, Chain_Type, AA_Code, Residue_Position, '\n', record.description, '\n')
                print(Database_Name, Database_ID, Name, Organism, Length, Chain_ID, Chain_Type, Uniprot_Protein_Name, AA_Code, Residue_Position, Residue_Position_Protein)
        else:
            continue

iterateThroughUniprotXMLFiles()