import pandas as pd
from pathlib import Path
import argparse


# merges three dataframes of contact_frequencies,
# and compute the mean and stdev of the contact_frequencies

def main(freq_paths:list[Path]) -> pd.DataFrame :
    """
    each freq_paths contains a dataframe of contact_frequencies
    whose cols are
    [residue_1_bw, residue_2_bw, residue_1, residue_2, contact_frequency]

    compute the mean and stdev of the contact_frequencies
    """

    dfs=[]
    for freq_path in freq_paths:
        df=pd.read_csv(freq_path, sep='\t')
        df['identifier'] = df.apply(lambda x: str(sorted([x['residue_1'], x['residue_2']])), axis=1)
        dfs.append(df)

    # Merge all dataframes on the identifier column
    merged = dfs[0]
    for df in dfs[1:]:
        merged = merged.merge(df, on="identifier", how="outer", suffixes=(None, None))

    # Preserve only one copy of residue_1 and residue_2 columns
    merged[['residue_1', 'residue_2']] = merged[['residue_1_x', 'residue_2_x']].fillna(
        merged[['residue_1_y', 'residue_2_y']]
    )
    merged = merged.drop(columns=[col for col in merged.columns if '_x' in col or '_y' in col])

    # Fill missing contact frequencies with 0
    contact_cols = [col for col in merged.columns if 'contact_frequency' in col]
    merged[contact_cols] = merged[contact_cols].fillna(0)

    # Calculate average and standard deviation for contact frequencies
    merged['freq_avg'] = merged[contact_cols].mean(axis=1)
    merged['freq_stdev'] = merged[contact_cols].std(axis=1)

    # Keep only required columns in the output
    result = merged[['residue_1', 'residue_2', 'freq_avg', 'freq_stdev']]

    return result

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="merge frequency tables")
    parser.add_argument(
        "-i", type=str, help="tsv files to merge", required=True, nargs='+')
    parser.add_argument(
        "-o", type=argparse.FileType('w'), help="output tsv file", required=True) 
    
    args = parser.parse_args()
    files=args.i
    output_path = args.o

    freq_paths=[]
    for file in files:
        if not Path(file).exists():
            raise FileNotFoundError(f"{file} does not exist")
        else:
            freq_paths.append(Path(file))

    merged_df=main(freq_paths)
