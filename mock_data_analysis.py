import os
import pandas as pd
from matplotlib import pyplot as plt
from sklearn.metrics import f1_score
from matplotlib_venn import venn3, venn3_circles

# принимаем на вход аннотации GTDB, Greengenes 2 и SILVA, 
# а также истинную аннотацию известных микроорганизмов (ожидаемые микроорганизмы)
def get_csv(file_name: str | bytes):
    return pd.read_csv(file_name) if os.path.isfile(file_name) else None


gtdb = get_csv("taxa_mock_gtdb.csv")
gg2 = get_csv("taxa_mock_gg2.csv")
silva = get_csv("taxa_mock_silva.csv")
true_data = get_csv('true_data_mock.csv')

gg2_cleaned = pd.DataFrame(gg2)
gtdb_cleaned = pd.DataFrame(gtdb)
silva = pd.DataFrame(silva)
true_data = pd.DataFrame(true_data)

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

common_true_gtdb = 0
common_true_gg2 = 0
common_true_silva = 0
print("Length", len(true_data["Sequence"].unique()))
true_data = true_data.sort_values(by='Sequence')
true_data.index = true_data["Sequence"]
gtdb.index = gtdb["Sequence"]
gg2.index = gg2["Sequence"]
silva.index = silva["Sequence"]

for true_data_index, _ in true_data.iterrows():
    for gtdb_index, _ in gtdb.iterrows():
        if true_data_index == gtdb_index:
            common_true_gtdb += 1
    for gg2_index, _ in gg2.iterrows():
        if true_data_index == gg2_index:
            common_true_gg2 += 1
    for silva_index, _ in silva.iterrows():
        if true_data_index == silva_index:
            common_true_silva += 1
print(common_true_gtdb)
print(common_true_gg2)
print(common_true_silva)

gtdb = gtdb.reindex(true_data.index)
gg2 = gg2.reindex(true_data.index)
silva = silva.reindex(true_data.index)

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


score_gtdb = 0
score_gg2 = 0
score_silva = 0
score_gtdb_silva = 0
score_gtdb_gg2 = 0
score_gg2_silva = 0
score_gtdb_gg2_silva = 0

taxonomy_levels = ['Phylum', 'Class', 'Order', 'Family', 'Genus']
#  вычисляем количество пересечений в аннотациях трех баз данных,
# сравниваем аннотации исследуемых баз данных с истинной аннотацией
score_dict = {}
for level in taxonomy_levels:
    score_dict[level] = {
        'gtdb': 0,
        'gg2': 0,
        'silva': 0,
        'gtdb_silva': 0,
        'gtdb_gg2': 0,
        'gg2_silva': 0,
        'gtdb_gg2_silva': 0
    }
true_data = true_data.drop(columns=["Sequence", "Species"])
gtdb = gtdb.drop(columns=["Sequence", "Species"])
gg2 = gg2.drop(columns=["Sequence", "Species"])
silva = silva.drop(columns=["Sequence"])

for level in taxonomy_levels:
    for idx, row in true_data.iterrows():
        try:
            true_value = true_data.loc[idx][level]
            gtdb_value = gtdb.loc[idx][level]
            gg2_value = gg2.loc[idx][level]
            silva_value = silva.loc[idx][level]

            if true_value == gtdb_value == gg2_value == silva_value:
                score_dict[level]['gtdb_gg2_silva'] += 1
                continue
            if true_value == gtdb_value == silva_value:
                score_dict[level]['gtdb_silva'] += 1
                continue
            if true_value == gg2_value == silva_value:
                score_dict[level]['gg2_silva'] += 1
                continue
            if true_value == gtdb_value == gg2_value:
                score_dict[level]['gtdb_gg2'] += 1
                continue
            if true_value == gtdb_value:
                score_dict[level]['gtdb'] += 1
                continue
            if true_value == gg2_value:
                score_dict[level]['gg2'] += 1
                continue
            if true_value == silva_value:
                score_dict[level]['silva'] += 1
                continue
        except ValueError:
            print("Same values:", true_value, gtdb_value, gg2_value, silva_value)
            continue

for level in taxonomy_levels:
    print(level)
    print(score_dict[level])

# создание диаграммы Венна
for level in taxonomy_levels:
    plt.figure(figsize=(10, 7))
    venn_labels = ('GTDB', 'GG2', 'Silva')
    venn_data = (
        score_dict[level]['gtdb'],
        score_dict[level]['gg2'],
        score_dict[level]['gtdb_gg2'],
        score_dict[level]['silva'],
        score_dict[level]['gtdb_silva'],
        score_dict[level]['gg2_silva'],
        score_dict[level]['gtdb_gg2_silva']
    )

    venn = venn3(subsets=venn_data, set_labels=venn_labels)
    plt.title(f'Перекрываемость таксономии на уровне {level} между базами данных GTDB, GG2 и SILVA ', fontdict={'fontsize': 25})
    plt.show()

# вычисление метрики F1
taxonomy_levels = ['Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']

def compute_f1_scores(reference, test, levels):
    f1_scores = {}
    for level in levels:
        if level in reference.columns and level in test.columns:
            true_labels = reference[level].tolist()
            predicted_labels = test[level].tolist()
            f1 = f1_score(true_labels, predicted_labels, average='weighted')
            f1_scores[level] = f1
    return f1_scores

f1_scores_gtdb = compute_f1_scores(true_data, gtdb, taxonomy_levels)
f1_scores_gg2 = compute_f1_scores(true_data, gg2, taxonomy_levels)
f1_scores_silva = compute_f1_scores(true_data, silva, taxonomy_levels)

print("F1 Scores - true vs. gtdb (Common Sequence):")
print(f1_scores_gtdb)

print("\nF1 Scores - true vs. gg2 (Common Sequence):")
print(f1_scores_gg2)

print("\nF1 Scores - true vs. silva (Common Sequence):")
print(f1_scores_silva)
