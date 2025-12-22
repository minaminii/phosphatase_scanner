from Bio import SeqIO
import pandas as pd

sequences = []
for record in SeqIO.parse("/content/TmF4.fasta.general.hit_sequence.fasta", "fasta"):
    sequences.append({'id': record.id, 'sequence': str(record.seq)})

fasta_file = "/content/TmF4.fasta.general.hit_sequence.fasta"
count = 0

for record in SeqIO.parse(fasta_file, "fasta"):
    count += 1

print("Number of sequences:", count)

df = pd.DataFrame(sequences)
df.head()
df.to_csv("TmF4Phosphatome.csv", index=False)

df = pd.read_csv('/content/TmF4.fasta.summary.tab', sep='\t')

# Save as a .csv file
df.to_csv('TmF4Phosphatome2.csv', index=False)

# Join/merge to .csv files # Change name
# Load your two CSV files
df1 = pd.read_csv("/content/TmF4Phosphatome2.csv")  # contains 'Query'
df2 = pd.read_csv("/content/TmF4Phosphatome.csv")  # contains 'id'

# Merge based on Query and id
merged_df = pd.merge(df1, df2, left_on='Query', right_on='id', how='inner')

# Drop original columns and create one unified ID column
merged_df['Prot_ID'] = merged_df['Query']  # or merged_df['id']
merged_df = merged_df.drop(columns=['Query', 'id'])

# Save or preview the result
merged_df.to_csv("merged_TmF4Phosphatome.csv", index=False)
print(merged_df.head())