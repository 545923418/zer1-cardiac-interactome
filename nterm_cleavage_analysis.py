import pandas as pd
import requests
import time
import os

FILE_PATH = "trans.xlsx"

def get_protein_sequence(accession):
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
    try:
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            fasta_data = response.text.split('\n', 1)
            if len(fasta_data) > 1:
                return fasta_data[1].replace('\n', '')
    except Exception as e:
        print(f"Failed to fetch {accession}: {e}")
    return None

def process_n_term(seq):
    if not seq or len(seq) < 2:
        return seq, "Sequence Too Short"

    first_aa = seq[0]
    second_aa = seq[1]
    
    loose_set = {'A', 'C', 'G', 'P', 'S', 'T', 'V'}
    strict_set = {'G', 'A', 'S', 'T', 'C'}
    
    if first_aa == 'M' and second_aa in loose_set:
        processed_seq = seq[1:]
        nt1 = processed_seq[0]
        
        if nt1 in strict_set:
            return processed_seq, "Strict"
        else:
            return processed_seq, "Loose"
    
    return seq, "Not Processed"

if os.path.exists(FILE_PATH):
    df = pd.read_excel(FILE_PATH)
    print(f"Loaded {len(df)} entries")
else:
    print(f"File not found: {FILE_PATH}")
    exit()

print("Fetching protein sequences...")
sequences = []
for index, row in df.iterrows():
    acc = row['Accession']
    if index % 10 == 0:
        print(f"Progress: {index}/{len(df)}")
        
    seq = get_protein_sequence(acc)
    sequences.append(seq)
    time.sleep(0.2)

df['Original_Sequence'] = sequences

print("Processing N-terminal cleavage...")
processed_seqs = []
groups = []

for seq in df['Original_Sequence']:
    p_seq, group = process_n_term(seq)
    processed_seqs.append(p_seq)
    groups.append(group)

df['Processed_Sequence'] = processed_seqs
df['Group'] = groups
df['Length'] = df['Processed_Sequence'].apply(lambda x: len(x) if x else 0)

df_final = df[df['Group'].isin(['Strict', 'Loose'])].copy()

print("\n" + "="*30)
print("Processing Summary")
print("="*30)
print(df_final[['Accession', 'Gene Symbol', 'Group', 'Processed_Sequence']].head())
print(f"\nTotal input: {len(df)}")
print(f"Met cleavage criteria: {len(df_final)}")
print(f"  - Strict (G/A/S/T/C): {len(df_final[df_final['Group'] == 'Strict'])}")
print(f"  - Loose (V/P): {len(df_final[df_final['Group'] == 'Loose'])}")

output_path = FILE_PATH.replace(".xlsx", "_output.xlsx")

try:
    with pd.ExcelWriter(output_path) as writer:
        df_final.to_excel(writer, sheet_name='Filtered_All', index=False)
        df_final[df_final['Group'] == 'Strict'].to_excel(writer, sheet_name='Strict_Group', index=False)
        df_final[df_final['Group'] == 'Loose'].to_excel(writer, sheet_name='Loose_Group', index=False)
        df.to_excel(writer, sheet_name='Raw_Data_With_Seq', index=False)
        
    print(f"\nResults saved to: {output_path}")
except Exception as e:
    print(f"Failed to save file: {e}")