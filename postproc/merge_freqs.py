import pandas as pd
from pathlib import Path
import argparse
from functools import reduce
import re
# merges three dataframes of contact_frequencies,
# and compute the mean and stdev of the contact_frequencies

def main(freq_paths:list[Path], freq_cutoff:float = 0.1) -> pd.DataFrame :
    """
    each freq_paths contains a dataframe of contact_frequencies
    whose cols are
    [residue_1_bw, residue_2_bw, residue_1, residue_2, contact_frequency]

    compute the mean and stdev of the contact_frequencies
    """

    dataframes=[]
    for freq_path in freq_paths:
        df=pd.read_csv(freq_path, sep='\t', index_col=None)
        df['identifier'] = df.apply(lambda x: str(sorted([x['residue_1'], x['residue_2']])), axis=1)
        dataframes.append(df)

    # Use functools.reduce to merge all DataFrames on the identifier column
    merged = reduce(lambda left, right: pd.merge(left, right, on="identifier", how="outer"), dataframes)

    # Reconcile residue_1 and residue_2 columns using regex
    residue_1_cols = [col for col in merged.columns if re.match(r'residue_1(_\w)?$', col)]
    residue_2_cols = [col for col in merged.columns if re.match(r'residue_2(_\w)?$', col)]

    merged['residue_1'] = merged[residue_1_cols].bfill(axis=1).iloc[:, 0]
    merged['residue_2'] = merged[residue_2_cols].bfill(axis=1).iloc[:, 0]

    # Reconcile residue_1 and residue_2 columns using regex
    residue_1_bw_cols = [col for col in merged.columns if re.match(r'residue_1_bw(_\w)?$', col)]
    residue_2_bw_cols = [col for col in merged.columns if re.match(r'residue_2_bw(_\w)?$', col)]

    merged['residue_1_bw'] = merged[residue_1_bw_cols].bfill(axis=1).iloc[:, 0]
    merged['residue_2_bw'] = merged[residue_2_bw_cols].bfill(axis=1).iloc[:, 0]

    # Fill missing contact frequencies with 0
    contact_cols = [col for col in merged.columns if 'contact_frequency' in col]
    merged[contact_cols] = merged[contact_cols].fillna(0)

    # Calculate average and standard deviation for contact frequencies
    merged['freq_avg'] = merged[contact_cols].mean(axis=1)
    merged['freq_stdev'] = merged[contact_cols].std(axis=1)

    # Keep all necessary columns in the output, including _bw columns
    result = merged[['residue_1_bw', 'residue_2_bw', 'residue_1', 'residue_2', 'freq_avg', 'freq_stdev']]
    result = result[result['freq_avg'] > freq_cutoff]

    return result

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="merge frequency tables")
    parser.add_argument(
        "-i", type=str, help="tsv files to merge", required=True, nargs='+')
    parser.add_argument(
        "-o", type=argparse.FileType('w'), help="output tsv file", required=True) 
    parser.add_argument(
        "-c", type=float, help="frequency cutoff", default=0.1) 
    
    args = parser.parse_args()
    files=args.i
    output_path = args.o
    cutoff = args.c

    freq_paths=[]
    for file in files:
        if not Path(file).exists():
            raise FileNotFoundError(f"{file} does not exist")
        else:
            freq_paths.append(Path(file))

    merged_df=main(freq_paths, cutoff)
    merged_df.to_csv(output_path, sep='\t', index=False)

