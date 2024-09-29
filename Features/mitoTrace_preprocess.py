import pandas as pd

def preprocess_columns(columns):
    return [col.lower().replace(' ', '').replace('-', '').replace('_', '').replace('(', '').replace(')', '') for col in columns]
def preprocess_tissues_file(tissues_file_path):
    tissues_file = pd.read_excel(tissues_file_path,'mitoTrace')
    split_features = tissues_file.columns[1].split('|')
    for feature in split_features:
            tissues_file[feature] = tissues_file.iloc[:,1]
    tissues_file.drop(columns=[tissues_file.columns[1],tissues_file.columns[12]], inplace=True)
    tissues_file.columns = preprocess_columns(tissues_file.columns)
    tissues_file = tissues_file.applymap(lambda x: x.lower().replace(' ', '').replace('-', '').replace('_', '').replace('(', '').replace(')','') if isinstance(x, str) else x)
    tissues_file.to_excel('tissues_file_features_split_tissues_preprocessed.xlsx')
    return tissues_file



def create_mitoTrace(Trace_file_path,mito_genes):
    Trace_file = pd.read_csv(Trace_file_path)
    mitoTrace = Trace_file[Trace_file['Gene_ID'].isin(mito_genes)]
    mitoTrace.columns = preprocess_columns(mitoTrace.columns)
    mitoTrace.to_excel('mitoTrace.xlsx')
    return mitoTrace
all_mito_genes = pd.read_excel(r"C:\Users\liorr\Dropbox\Lior\mitochondrial genes lists\mitoCarta plus review\merged_mitochondrial_genes_MitoCarta and review.xlsx")
tissues_file = preprocess_tissues_file(r"C:\Users\liorr\Dropbox\Lior\mitoModel tissues associations.xlsx")

create_mitoTrace(r"C:\Users\liorr\Documents\Ben Gurion University\Yeger-Lotem lab\projects\Mitochondrial diseases\Files\Aneuploidy\5_TRACE_Dataset_with_RNAi.csv", all_mito_genes['First ENSG ID'])

