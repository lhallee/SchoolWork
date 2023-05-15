import requests
import pandas as pd
from io import StringIO
from Bio import SeqIO
from tqdm import tqdm

# Get the list of available datasets

df = pd.read_csv('martservice.txt', delimiter='\t', header=None)
datasets = df[1].tolist()


url = 'http://www.ensembl.org/biomart/martservice?type=datasets&mart=ENSEMBL_MART_ENSEMBL'
response = requests.get(url)
print(response)
'''
datasets = response.text.split('\n')
datasets = [dataset.split('\t')[0] for dataset in datasets if dataset]
'''

'''
# Initialize an empty DataFrame
coding_sequences_df = pd.DataFrame(columns=['ID', 'Sequence'])

# Loop through the datasets and download the Coding Sequences
for dataset in tqdm(datasets, desc='Downloading Coding Sequences'):
	fasta_url = f'http://www.ensembl.org/biomart/martservice?type=fasta&mart=ensembl&dataset={dataset}&attributes=coding'
	fasta_response = requests.get(fasta_url)
	print(fasta_response.text)
	# Parse the downloaded FASTA data
	fasta_data = StringIO(fasta_response.text)
	sequences = SeqIO.parse(fasta_data, 'fasta')

	# Add the sequences to the DataFrame
	for seq in sequences:
		coding_sequences_df = coding_sequences_df.append(
			{'ID': seq.id, 'Sequence': str(seq.seq)}, ignore_index=True
		)

# Save the DataFrame to a CSV file
coding_sequences_df.to_csv('coding_sequences.csv', index=False)

'''

