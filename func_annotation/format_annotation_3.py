import sys

# Processa o arquivo eggnog diretamente, linha por linha
dict_annot = {}

with open(sys.argv[1], 'r') as eggnog_file:
    for line in eggnog_file:
        if line.startswith("#"):
            continue

        cols = line.strip().split("\t")
        key = cols[0][:-7]

        if cols[8] != "-":
            go_terms = cols[9].split(",") if isinstance(cols[9], str) else cols[9]
            dict_annot[key] = [cols[8], cols[7], go_terms]
        else:
            go_terms = cols[9].split(",") if isinstance(cols[9], str) else cols[9]
            dict_annot[key] = [key, cols[7], go_terms]

# Cria um dicionário de busca para o arquivo pannzer2
pannzer2_dict = {}

with open(sys.argv[2], 'r') as pannzer2_file:
    for line in pannzer2_file:
        cols = line.split("\t")
        key = cols[0][:-7]
        go_term = "GO:" + str(cols[2])
        if key in pannzer2_dict:
            pannzer2_dict[key].append(go_term)
        else:
            pannzer2_dict[key] = [go_term]

# Preenche o dict_final_annot
dict_final_annot = {}

for k, v in dict_annot.items():
    if k in pannzer2_dict:
        for go_term in pannzer2_dict[k]:
            if go_term not in v[2]:  # Verifica se o termo GO já está na lista dentro de v[2]
                v[2].append(go_term)  # Adiciona o novo termo GO na lista que está na posição 2
        dict_final_annot.setdefault(k, []).extend(v)  # Armazena ou atualiza a lista completa em dict_final_annot

# Exemplo de uso do dict_final_annot
# for k, v in dict_final_annot.items():
#     print(k, v)

count = 1

with open('go_annotation_2.tsv', 'w') as go_file:
    go_file.write("GID\tGO\tEVIDENCE\n")
    for k, v in dict_final_annot.items():
        for go in v[2]:
            if go != "-":
                go_file.write(f"{count}\t{go}\tIEA\n")
        count += 1

count = 1

with open('gene_info_2.txt', 'w') as name_file:
    name_file.write("GID\tSYMBOL\tGENENAME\n")
    for k, v in dict_final_annot.items():
        if v[0] != "-":
            name_file.write(f"{count}\t{k}\t{v[1]}\n")
        count += 1
import sys

# Processa o arquivo eggnog diretamente, linha por linha
dict_annot = {}

with open(sys.argv[1], 'r') as eggnog_file:
    for line in eggnog_file:
        if line.startswith("#"):
            continue

        cols = line.strip().split("\t")
        key = cols[0][:-7]

        if cols[8] != "-":
            go_terms = cols[9].split(",") if isinstance(cols[9], str) else cols[9]
            dict_annot[key] = [cols[8], cols[7], go_terms]
        else:
            go_terms = cols[9].split(",") if isinstance(cols[9], str) else cols[9]
            dict_annot[key] = [key, cols[7], go_terms]

# Cria um dicionário de busca para o arquivo pannzer2
pannzer2_dict = {}

with open(sys.argv[2], 'r') as pannzer2_file:
    for line in pannzer2_file:
        cols = line.split("\t")
        key = cols[0][:-7]
        go_term = "GO:" + str(cols[2])
        if key in pannzer2_dict:
            pannzer2_dict[key].append(go_term)
        else:
            pannzer2_dict[key] = [go_term]

# Preenche o dict_final_annot
dict_final_annot = {}

for k, v in dict_annot.items():
    if k in pannzer2_dict:
        for go_term in pannzer2_dict[k]:
            if go_term not in v[2]:  # Verifica se o termo GO já está na lista dentro de v[2]
                v[2].append(go_term)  # Adiciona o novo termo GO na lista que está na posição 2
        dict_final_annot.setdefault(k, []).extend(v)  # Armazena ou atualiza a lista completa em dict_final_annot

# Exemplo de uso do dict_final_annot
# for k, v in dict_final_annot.items():
#     print(k, v)

count = 1

with open('go_annotation_2.tsv', 'w') as go_file:
    go_file.write("GID\tGO\tEVIDENCE\n")
    for k, v in dict_final_annot.items():
        for go in v[2]:
            if go != "-":
                go_file.write(f"{count}\t{go}\tIEA\n")
        count += 1

count = 1

with open('gene_info_2.txt', 'w') as name_file:
    name_file.write("GID\tSYMBOL\tGENENAME\n")
    for k, v in dict_final_annot.items():
        if v[0] != "-":
            name_file.write(f"{count}\t{k}\t{v[1]}\n")
        count += 1
import sys

# Processa o arquivo eggnog diretamente, linha por linha
dict_annot = {}

with open(sys.argv[1], 'r') as eggnog_file:
    for line in eggnog_file:
        if line.startswith("#"):
            continue

        cols = line.strip().split("\t")
        key = cols[0][:-7]

        if cols[8] != "-":
            go_terms = cols[9].split(",") if isinstance(cols[9], str) else cols[9]
            dict_annot[key] = [cols[8], cols[7], go_terms]
        else:
            go_terms = cols[9].split(",") if isinstance(cols[9], str) else cols[9]
            dict_annot[key] = [key, cols[7], go_terms]

# Cria um dicionário de busca para o arquivo pannzer2
pannzer2_dict = {}

with open(sys.argv[2], 'r') as pannzer2_file:
    for line in pannzer2_file:
        cols = line.split("\t")
        key = cols[0][:-7]
        go_term = "GO:" + str(cols[2])
        if key in pannzer2_dict:
            pannzer2_dict[key].append(go_term)
        else:
            pannzer2_dict[key] = [go_term]

# Preenche o dict_final_annot
dict_final_annot = {}

for k, v in dict_annot.items():
    if k in pannzer2_dict:
        for go_term in pannzer2_dict[k]:
            if go_term not in v[2]:  # Verifica se o termo GO já está na lista dentro de v[2]
                v[2].append(go_term)  # Adiciona o novo termo GO na lista que está na posição 2
        dict_final_annot.setdefault(k, []).extend(v)  # Armazena ou atualiza a lista completa em dict_final_annot

# Exemplo de uso do dict_final_annot
# for k, v in dict_final_annot.items():
#     print(k, v)

count = 1

with open('go_annotation_2.tsv', 'w') as go_file:
    go_file.write("GID\tGO\tEVIDENCE\n")
    for k, v in dict_final_annot.items():
        for go in v[2]:
            if go != "-":
                go_file.write(f"{count}\t{go}\tIEA\n")
        count += 1

count = 1

with open('gene_info_2.txt', 'w') as name_file:
    name_file.write("GID\tSYMBOL\tGENENAME\n")
    for k, v in dict_final_annot.items():
        if v[0] != "-":
            name_file.write(f"{count}\t{k}\t{v[1]}\n")
        count += 1
import sys

# Processa o arquivo eggnog diretamente, linha por linha
dict_annot = {}

with open(sys.argv[1], 'r') as eggnog_file:
    for line in eggnog_file:
        if line.startswith("#"):
            continue

        cols = line.strip().split("\t")
        key = cols[0][:-7]

        if cols[8] != "-":
            go_terms = cols[9].split(",") if isinstance(cols[9], str) else cols[9]
            dict_annot[key] = [cols[8], cols[7], go_terms]
        else:
            go_terms = cols[9].split(",") if isinstance(cols[9], str) else cols[9]
            dict_annot[key] = [key, cols[7], go_terms]

# Cria um dicionário de busca para o arquivo pannzer2
pannzer2_dict = {}

with open(sys.argv[2], 'r') as pannzer2_file:
    for line in pannzer2_file:
        cols = line.split("\t")
        key = cols[0][:-7]
        go_term = "GO:" + str(cols[2])
        if key in pannzer2_dict:
            pannzer2_dict[key].append(go_term)
        else:
            pannzer2_dict[key] = [go_term]

# Preenche o dict_final_annot
dict_final_annot = {}

for k, v in dict_annot.items():
    if k in pannzer2_dict:
        for go_term in pannzer2_dict[k]:
            if go_term not in v[2]:  # Verifica se o termo GO já está na lista dentro de v[2]
                v[2].append(go_term)  # Adiciona o novo termo GO na lista que está na posição 2
        dict_final_annot.setdefault(k, []).extend(v)  # Armazena ou atualiza a lista completa em dict_final_annot

# Exemplo de uso do dict_final_annot
# for k, v in dict_final_annot.items():
#     print(k, v)

count = 1

with open('go_annotation_2.tsv', 'w') as go_file:
    go_file.write("GID\tGO\tEVIDENCE\n")
    for k, v in dict_final_annot.items():
        for go in v[2]:
            if go != "-":
                go_file.write(f"{count}\t{go}\tIEA\n")
        count += 1

count = 1

with open('gene_info_2.txt', 'w') as name_file:
    name_file.write("GID\tSYMBOL\tGENENAME\n")
    for k, v in dict_final_annot.items():
        if v[0] != "-":
            name_file.write(f"{count}\t{k}\t{v[1]}\n")
        count += 1
import sys

# Processa o arquivo eggnog diretamente, linha por linha
dict_annot = {}

with open(sys.argv[1], 'r') as eggnog_file:
    for line in eggnog_file:
        if line.startswith("#"):
            continue

        cols = line.strip().split("\t")
        key = cols[0][:-7]

        if cols[8] != "-":
            go_terms = cols[9].split(",") if isinstance(cols[9], str) else cols[9]
            dict_annot[key] = [cols[8], cols[7], go_terms]
        else:
            go_terms = cols[9].split(",") if isinstance(cols[9], str) else cols[9]
            dict_annot[key] = [key, cols[7], go_terms]

# Cria um dicionário de busca para o arquivo pannzer2
pannzer2_dict = {}

with open(sys.argv[2], 'r') as pannzer2_file:
    for line in pannzer2_file:
        cols = line.split("\t")
        key = cols[0][:-7]
        go_term = "GO:" + str(cols[2])
        if key in pannzer2_dict:
            pannzer2_dict[key].append(go_term)
        else:
            pannzer2_dict[key] = [go_term]

# Preenche o dict_final_annot
dict_final_annot = {}

for k, v in dict_annot.items():
    if k in pannzer2_dict:
        for go_term in pannzer2_dict[k]:
            if go_term not in v[2]:  # Verifica se o termo GO já está na lista dentro de v[2]
                v[2].append(go_term)  # Adiciona o novo termo GO na lista que está na posição 2
        dict_final_annot.setdefault(k, []).extend(v)  # Armazena ou atualiza a lista completa em dict_final_annot

# Exemplo de uso do dict_final_annot
# for k, v in dict_final_annot.items():
#     print(k, v)

count = 1

with open('go_annotation_2.tsv', 'w') as go_file:
    go_file.write("GID\tGO\tEVIDENCE\n")
    for k, v in dict_final_annot.items():
        for go in v[2]:
            if go != "-":
                go_file.write(f"{count}\t{go}\tIEA\n")
        count += 1

count = 1

with open('gene_info_2.txt', 'w') as name_file:
    name_file.write("GID\tSYMBOL\tGENENAME\n")
    for k, v in dict_final_annot.items():
        if v[0] != "-":
            name_file.write(f"{count}\t{k}\t{v[1]}\n")
        count += 1
import sys

# Processa o arquivo eggnog diretamente, linha por linha
dict_annot = {}

with open(sys.argv[1], 'r') as eggnog_file:
    for line in eggnog_file:
        if line.startswith("#"):
            continue

        cols = line.strip().split("\t")
        key = cols[0][:-7]

        if cols[8] != "-":
            go_terms = cols[9].split(",") if isinstance(cols[9], str) else cols[9]
            dict_annot[key] = [cols[8], cols[7], go_terms]
        else:
            go_terms = cols[9].split(",") if isinstance(cols[9], str) else cols[9]
            dict_annot[key] = [key, cols[7], go_terms]

# Cria um dicionário de busca para o arquivo pannzer2
pannzer2_dict = {}

with open(sys.argv[2], 'r') as pannzer2_file:
    for line in pannzer2_file:
        cols = line.split("\t")
        key = cols[0][:-7]
        go_term = "GO:" + str(cols[2])
        if key in pannzer2_dict:
            pannzer2_dict[key].append(go_term)
        else:
            pannzer2_dict[key] = [go_term]

# Preenche o dict_final_annot
dict_final_annot = {}

for k, v in dict_annot.items():
    if k in pannzer2_dict:
        for go_term in pannzer2_dict[k]:
            if go_term not in v[2]:  # Verifica se o termo GO já está na lista dentro de v[2]
                v[2].append(go_term)  # Adiciona o novo termo GO na lista que está na posição 2
        dict_final_annot.setdefault(k, []).extend(v)  # Armazena ou atualiza a lista completa em dict_final_annot

# Exemplo de uso do dict_final_annot
# for k, v in dict_final_annot.items():
#     print(k, v)

count = 1

with open('go_annotation_2.tsv', 'w') as go_file:
    go_file.write("GID\tGO\tEVIDENCE\n")
    for k, v in dict_final_annot.items():
        for go in v[2]:
            if go != "-":
                go_file.write(f"{count}\t{go}\tIEA\n")
        count += 1

count = 1

with open('gene_info_2.txt', 'w') as name_file:
    name_file.write("GID\tSYMBOL\tGENENAME\n")
    for k, v in dict_final_annot.items():
        if v[0] != "-":
            name_file.write(f"{count}\t{k}\t{v[1]}\n")
        count += 1
import sys

# Processa o arquivo eggnog diretamente, linha por linha
dict_annot = {}

with open(sys.argv[1], 'r') as eggnog_file:
    for line in eggnog_file:
        if line.startswith("#"):
            continue

        cols = line.strip().split("\t")
        key = cols[0][:-7]

        if cols[8] != "-":
            go_terms = cols[9].split(",") if isinstance(cols[9], str) else cols[9]
            dict_annot[key] = [cols[8], cols[7], go_terms]
        else:
            go_terms = cols[9].split(",") if isinstance(cols[9], str) else cols[9]
            dict_annot[key] = [key, cols[7], go_terms]

# Cria um dicionário de busca para o arquivo pannzer2
pannzer2_dict = {}

with open(sys.argv[2], 'r') as pannzer2_file:
    for line in pannzer2_file:
        cols = line.split("\t")
        key = cols[0][:-7]
        go_term = "GO:" + str(cols[2])
        if key in pannzer2_dict:
            pannzer2_dict[key].append(go_term)
        else:
            pannzer2_dict[key] = [go_term]

# Preenche o dict_final_annot
dict_final_annot = {}

for k, v in dict_annot.items():
    if k in pannzer2_dict:
        for go_term in pannzer2_dict[k]:
            if go_term not in v[2]:  # Verifica se o termo GO já está na lista dentro de v[2]
                v[2].append(go_term)  # Adiciona o novo termo GO na lista que está na posição 2
        dict_final_annot.setdefault(k, []).extend(v)  # Armazena ou atualiza a lista completa em dict_final_annot

# Exemplo de uso do dict_final_annot
# for k, v in dict_final_annot.items():
#     print(k, v)

count = 1

with open('go_annotation_2.tsv', 'w') as go_file:
    go_file.write("GID\tGO\tEVIDENCE\n")
    for k, v in dict_final_annot.items():
        for go in v[2]:
            if go != "-":
                go_file.write(f"{count}\t{go}\tIEA\n")
        count += 1

count = 1

with open('gene_info_2.txt', 'w') as name_file:
    name_file.write("GID\tSYMBOL\tGENENAME\n")
    for k, v in dict_final_annot.items():
        if v[0] != "-":
            name_file.write(f"{count}\t{k}\t{v[1]}\n")
        count += 1
import sys

# Processa o arquivo eggnog diretamente, linha por linha
dict_annot = {}

with open(sys.argv[1], 'r') as eggnog_file:
    for line in eggnog_file:
        if line.startswith("#"):
            continue

        cols = line.strip().split("\t")
        key = cols[0][:-7]

        if cols[8] != "-":
            go_terms = cols[9].split(",") if isinstance(cols[9], str) else cols[9]
            dict_annot[key] = [cols[8], cols[7], go_terms]
        else:
            go_terms = cols[9].split(",") if isinstance(cols[9], str) else cols[9]
            dict_annot[key] = [key, cols[7], go_terms]

# Cria um dicionário de busca para o arquivo pannzer2
pannzer2_dict = {}

with open(sys.argv[2], 'r') as pannzer2_file:
    for line in pannzer2_file:
        cols = line.split("\t")
        key = cols[0][:-7]
        go_term = "GO:" + str(cols[2])
        if key in pannzer2_dict:
            pannzer2_dict[key].append(go_term)
        else:
            pannzer2_dict[key] = [go_term]

# Preenche o dict_final_annot
dict_final_annot = {}

for k, v in dict_annot.items():
    if k in pannzer2_dict:
        for go_term in pannzer2_dict[k]:
            if go_term not in v[2]:  # Verifica se o termo GO já está na lista dentro de v[2]
                v[2].append(go_term)  # Adiciona o novo termo GO na lista que está na posição 2
        dict_final_annot.setdefault(k, []).extend(v)  # Armazena ou atualiza a lista completa em dict_final_annot

# Exemplo de uso do dict_final_annot
# for k, v in dict_final_annot.items():
#     print(k, v)

count = 1

with open('go_annotation_2.tsv', 'w') as go_file:
    go_file.write("GID\tGO\tEVIDENCE\n")
    for k, v in dict_final_annot.items():
        for go in v[2]:
            if go != "-":
                go_file.write(f"{count}\t{go}\tIEA\n")
        count += 1

count = 1

with open('gene_info_2.txt', 'w') as name_file:
    name_file.write("GID\tSYMBOL\tGENENAME\n")
    for k, v in dict_final_annot.items():
        if v[0] != "-":
            name_file.write(f"{count}\t{k}\t{v[1]}\n")
        count += 1
import sys

# Processa o arquivo eggnog diretamente, linha por linha
dict_annot = {}

with open(sys.argv[1], 'r') as eggnog_file:
    for line in eggnog_file:
        if line.startswith("#"):
            continue

        cols = line.strip().split("\t")
        key = cols[0][:-7]

        if cols[8] != "-":
            go_terms = cols[9].split(",") if isinstance(cols[9], str) else cols[9]
            dict_annot[key] = [cols[8], cols[7], go_terms]
        else:
            go_terms = cols[9].split(",") if isinstance(cols[9], str) else cols[9]
            dict_annot[key] = [key, cols[7], go_terms]

# Cria um dicionário de busca para o arquivo pannzer2
pannzer2_dict = {}

with open(sys.argv[2], 'r') as pannzer2_file:
    for line in pannzer2_file:
        cols = line.split("\t")
        key = cols[0][:-7]
        go_term = "GO:" + str(cols[2])
        if key in pannzer2_dict:
            pannzer2_dict[key].append(go_term)
        else:
            pannzer2_dict[key] = [go_term]

# Preenche o dict_final_annot
dict_final_annot = {}

for k, v in dict_annot.items():
    if k in pannzer2_dict:
        for go_term in pannzer2_dict[k]:
            if go_term not in v[2]:  # Verifica se o termo GO já está na lista dentro de v[2]
                v[2].append(go_term)  # Adiciona o novo termo GO na lista que está na posição 2
        dict_final_annot.setdefault(k, []).extend(v)  # Armazena ou atualiza a lista completa em dict_final_annot

# Exemplo de uso do dict_final_annot
# for k, v in dict_final_annot.items():
#     print(k, v)

count = 1

with open('go_annotation_2.tsv', 'w') as go_file:
    go_file.write("GID\tGO\tEVIDENCE\n")
    for k, v in dict_final_annot.items():
        for go in v[2]:
            if go != "-":
                go_file.write(f"{count}\t{go}\tIEA\n")
        count += 1

count = 1

with open('gene_info_2.txt', 'w') as name_file:
    name_file.write("GID\tSYMBOL\tGENENAME\n")
    for k, v in dict_final_annot.items():
        if v[0] != "-":
            name_file.write(f"{count}\t{k}\t{v[1]}\n")
        count += 1

