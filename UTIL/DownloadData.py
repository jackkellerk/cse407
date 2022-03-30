import yaml
import urllib.request

# Open pdb id config file
with open('../CONFIG/config_data.yaml', 'r') as f:
    pdb_ids = yaml.load(f, Loader=yaml.FullLoader)

# Download the data files from pdb
for channel in pdb_ids.keys():
    for gene in pdb_ids[channel]:
        for structure in list(gene.values())[0]:
            try:
                urllib.request.urlretrieve('https://files.rcsb.org/download/' + structure + '.pdb', '../DATA/RAW/' + channel + '/' + structure + '.pdb')
            except:
                print('Unable to download structure: ', structure)