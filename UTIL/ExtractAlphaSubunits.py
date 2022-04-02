import yaml

# Open pdb id config file
with open('../CONFIG/config_data.yaml', 'r') as f:
    pdb_ids = yaml.load(f, Loader=yaml.FullLoader)

# Open the data files and loop over the channels, gene names, and structure ids
for channel in pdb_ids.keys():
    for gene in pdb_ids[channel]:
        for structure in list(gene.values())[0]:

            # Open the selected pdb file and return an array of the lines in the pdb file
            with open('../DATA/RAW/' + channel + '/' + structure + '.pdb') as pdb_structure:
                lines = pdb_structure.readlines()
            
            # Loop over the lines in the pdb structure and find the chain ids within the DBREF lines of the pdb files
            chain_ids = []
            for line in lines:
                gene_name = list(gene.keys())[0]
                if 'DBREF' in line and gene_name in line:
                    chain_ids.append(line.split()[2])
            
            # If unable to extract any chain ids for the structure, print a warning
            if len(chain_ids) < 1:
                print("Warning: Unable to extract any chain ids for " + structure + " of " + list(gene.keys())[0])
            
            # Loop over the lines in the pdb structure and find the ATOM lines that are associated with the structure's 
            # chain ids. Then, write that line to the file within the SUBUNITS data folder.
            with open('../DATA/SUBUNITS/' + channel + '/' + structure + '.pdb', 'w') as extracted_subunit_pdb:
                for line in lines:
                    if not 'REMARK' in line and 'ATOM' in line:
                        if line.split()[4][0] in chain_ids:
                            extracted_subunit_pdb.write(line)