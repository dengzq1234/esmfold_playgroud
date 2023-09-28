import argparse
import requests
import biotite.structure.io as bsio
from Bio import SeqIO

def fold_sequence(sequence):
    headers = {
        'Content-Type': 'application/x-www-form-urlencoded',
    }
    response = requests.post('https://api.esmatlas.com/foldSequence/v1/pdb/', headers=headers, data=sequence)
    pdb_string = response.content.decode('utf-8')
    return pdb_string

def parse_fasta(fastafile):
    fasta_dict = {}
    with open(fastafile,'r') as f:
        head = ''
        seq = ''
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if seq != '':
                    fasta_dict[head] = seq
                    seq = ''
                    head = line[1:]
                else:
                    head = line[1:]
            else:
                seq += line
    fasta_dict[head] = seq
    return fasta_dict

def main():
    parser = argparse.ArgumentParser(description='Predict protein structure from a FASTA file using ESMFold and save as PDB.')
    parser.add_argument('input_fasta', type=str, help='Input FASTA file.')
    #parser.add_argument('output_pdb', type=str, help='Output PDB file.')
    
    args = parser.parse_args()

    # Parse the first sequence from the input FASTA file
    fasta_dict = parse_fasta(args.input_fasta)
    
    for name,seq in fasta_dict.items():

        # Predict the structure
        pdb_string = fold_sequence(seq)

        # Save the predicted structure
        with open(name+'.pdb', 'w') as pdb_file:
            pdb_file.write(pdb_string)
    
        # Load structure and calculate pLDDT
        struct = bsio.load_structure(name+'.pdb', extra_fields=["b_factor"])
        b_value = round(struct.b_factor.mean(), 4)
    
        print('pLDDT for sequence:', name, 'is', b_value)

if __name__ == "__main__":
    main()
