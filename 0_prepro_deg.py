# ------------------
# Script to normalize DEG database CSV file
# ------------------
import pandas as pd

raw_path = "genes/deg_annotation_p.csv"
norm_path = "genes/deg.csv"

column_names = [
    "general", "deg_id", "gene_name", "gene_ref", "cog",
    "function_class", "description", "organism", "refseq", "condition",
    "unnamed_placeholder",
    "gene_ontology", "uniprot"
]

try:
    df = pd.read_csv(
        raw_path,
        sep=';',
        header=None,
        quotechar='"'
    )

    if df.shape[1] > len(column_names):
        df = df.iloc[:, :-1]

    df.columns = column_names
    print(df.head())
    print(df.info())

    df.to_csv(norm_path, header=True, index=False)

except Exception as e:
    print(f"An error occurred: {e}")