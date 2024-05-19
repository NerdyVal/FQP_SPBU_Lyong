import os
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib_venn import venn3, venn3_circles
from sklearn.metrics import cohen_kappa_score

# принимаем на вход аннотации GTDB,Greengenes2 и SILVA
def get_csv(file_name: str | bytes):
    return pd.read_csv(file_name) if os.path.isfile(file_name) else None

gtdb = get_csv("taxa_sample_gtdb.csv")
gg2 = get_csv("taxa_sample_gg2.csv")
silva = get_csv("taxa_sample_silva.csv")

gg2_cleaned = pd.DataFrame(gg2)
gtdb_cleaned = pd.DataFrame(gtdb)
silva = pd.DataFrame(silva)


# стандартизация формата аннотаций для каждой базы данных
def remove_first_three_chars(s):
    if pd.notna(s) and isinstance(s, str):
        return s[3:]
    else:
        return s

def remove_split(s):
    if pd.notna(s) and isinstance(s, str):
        return s.split("_")[0]
    else:
        return s


columns_to_process = [col for col in gg2_cleaned.columns if col != 'Sequence']

gg2_cleaned[columns_to_process] = gg2_cleaned[columns_to_process].applymap(remove_first_three_chars)
gg2_cleaned[columns_to_process] = gg2_cleaned[columns_to_process].applymap(remove_split)
gtdb_cleaned[columns_to_process] = gtdb_cleaned[columns_to_process].applymap(remove_split)

gg2 = gg2_cleaned
gtdb = gtdb_cleaned


replace_dictionary = {
    "Bacillota": "Firmicutes",
    "Pseudomonadota": "Proteobacteria",
    "Epsilonproteobacteria": "Gammaproteobacteria",
    "Escherichia-Shigella": "Escherihia"
}


def replace_values(df):
    return df.replace(replace_dictionary)


gtdb = replace_values(gtdb)
gg2 = replace_values(gg2)
silva = replace_values(silva)


taxonomy_levels = ['Phylum', 'Class', 'Order', 'Family', 'Genus']
# вычисляем коэффициент Жаккара для каждого ранга и каждой пары баз данных
jaccard_scores = {}
for rank in taxonomy_levels:
    for db1, df1 in {'gtdb': gtdb, 'silva': silva, 'gg2': gg2}.items():
        for db2, df2 in {'gtdb': gtdb, 'silva': silva, 'gg2': gg2}.items():
            if db1 != db2:
                set1 = set(df1[rank])
                set2 = set(df2[rank])
                jaccard = len(set1.intersection(set2)) / len(set1.union(set2))
                jaccard_scores[(rank, db1, db2)] = jaccard

for (rank, db1, db2), score in sorted(jaccard_scores.items()):
    print(f'{rank} - {db1} vs. {db2}: {score:.4f}')


for rank in taxonomy_levels:
    gtdb[rank] = gtdb[rank].astype('category').cat.codes
    silva[rank] = silva[rank].astype('category').cat.codes
    gg2[rank] = gg2[rank].astype('category').cat.codes

# вычисляем коэффициент Коэна для каждого ранга и каждой пары баз данных
kappa_scores = {}
for rank in taxonomy_levels:
    for db1, df1 in {'gtdb': gtdb, 'silva': silva, 'gg2': gg2}.items():
        for db2, df2 in {'gtdb': gtdb, 'silva': silva, 'gg2': gg2}.items():
            if db1 != db2:
                kappa = cohen_kappa_score(df1[rank], df2[rank])
                kappa_scores[(rank, db1, db2)] = kappa

for (rank, db1, db2), score in sorted(kappa_scores.items()):
    print(f'Каппа {rank} - {db1} vs. {db2}: {score:.4f}')



for taxonomy_level in taxonomy_levels:
    # извлечение уникальных значений из колонки 'Sequence' и колонки текущего уровня таксономии каждой базы данных
    gtdb_values = set(gtdb.loc[:, [taxonomy_level, 'Sequence']].dropna().drop_duplicates().iloc[:, 0])
    gg2_values = set(gg2.loc[:, [taxonomy_level, 'Sequence']].dropna().drop_duplicates().iloc[:, 0])
    silva_values = set(silva.loc[:, [taxonomy_level, 'Sequence']].dropna().drop_duplicates().iloc[:, 0])

    # вычисление пересечений между базами данных
    gtdb_gg2 = gtdb_values & gg2_values
    gtdb_silva = gtdb_values & silva_values
    gg2_silva = gg2_values & silva_values
    gtdb_gg2_silva = gtdb_values & gg2_values & silva_values

    # создание словаря для хранения данных для диаграммы Венна
    venn_data = {
        'GTDB': len(gtdb_values),
        'GG2': len(gg2_values),
        'Silva': len(silva_values),
        'GTDB & GG2': len(gtdb_gg2),
        'GTDB & Silva': len(gtdb_silva),
        'GG2 & Silva': len(gg2_silva),
        'GTDB & GG2 & Silva': len(gtdb_gg2_silva),
    }

    # вывод числовых значений пересечений
    print(f"Пересечения на уровне {taxonomy_level}:")
    for key, value in venn_data.items():
        print(f"{key}: {value}")

    # создание диаграммы Венна
    plt.figure(figsize=(10, 7))
    v = venn3([gtdb_values, gg2_values, silva_values], set_labels=('GTDB', 'GG2', 'SILVA'))
    for text in v.subset_labels:
        text.set_fontsize(14)  
    plt.title(f'Перекрываемость таксономии на уровне {taxonomy_level} между базами данных GTDB, GG2 и SILVA ', fontdict={'fontsize': 25})
    plt.show()

