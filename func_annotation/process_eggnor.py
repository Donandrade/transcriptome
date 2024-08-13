import csv


import csv

def process_eggnog_tsv(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        
        # Ignora as linhas de comentário iniciais
        for line in reader:
            if line[0].startswith('#query'):
                header = line
                break
        
        reader = csv.DictReader(infile, delimiter='\t', fieldnames=header)
        writer = csv.writer(outfile, delimiter='\t')
        
        for row in reader:
            protein_id = row['#query']
            go_terms = row['GOs'].split(',') if row['GOs'] else ['-']
            
            for go_term in go_terms:
                writer.writerow([protein_id, go_term])

# Especifique o caminho para o arquivo TSV de entrada e o arquivo de saída desejado
input_file = 'MM_3jfpbrke.emapper.annotations.tsv'
output_file = 'output_go_terms.tsv'

# Processa o arquivo
process_eggnog_tsv(input_file, output_file)

