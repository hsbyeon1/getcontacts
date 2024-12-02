import pandas as pd
from pathlib import Path
import argparse

# converts native numbering to bw numbering
# drops the rows with diff_resnum < cutoff 

AMBER_NONCAN = {
    "ASH": "ASP",
    "GLH": "GLU",
    "CYX": "CYS",
    "HID": "HIS",
    "HIE": "HIS",
    "HIP": "HIS",
    "LYN": "LYS",
    "CYM": "CYS",
    "CYF": "CYS",
    "CYR": "CYS",
    "CYT": "CYS",
}

RESIDUES_3_TO_1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "GLN": "Q",
    "SER": "S",
    "SEC": "U",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}

def converter(query:str, label_dict:dict[str,str]) -> str:
    query_sp = query.split(':')
    if query_sp[1] in AMBER_NONCAN:
        query_sp = [query_sp[0], AMBER_NONCAN[query_sp[1]], query_sp[2]]
    elif query_sp[1] in RESIDUES_3_TO_1:
        pass
    else:
        raise ValueError(f"Unknown residue: {query_sp[1]}")

    if query in label_dict:
        return label_dict[query]
    else:
        return f"{RESIDUES_3_TO_1[query_sp[1]]}{query_sp[2]}"

def resno_diff(res1:str, res2:str) -> int:
    """
    inputs are in form of X:ARG:1
    """

    resno1 = int(res1.split(':')[-1])
    resno2 = int(res2.split(':')[-1])

    return abs(resno1-resno2)

def main(label_path:Path, freq_path:Path, output_path:Path, cutoff:int=4) -> pd.DataFrame:
    label_df = pd.read_csv(label_path, sep='\t', header=None, names=['native','bw','color'],index_col=None)
    label_dict = label_df.set_index('native')['bw'].to_dict()

    freq_df = pd.read_csv(freq_path, sep='\t', comment='#', names=['residue_1','residue_2','contact_frequency'])

    freq_df['close_res'] = freq_df.apply(lambda x: resno_diff(x['residue_1'], x['residue_2']) < cutoff, axis=1)
    freq_df = freq_df[~freq_df['close_res']]
    freq_df.drop(columns=['close_res'], inplace=True)

    freq_df['residue_1_bw'] = freq_df['residue_1'].apply(lambda x: converter(x, label_dict))
    freq_df['residue_2_bw'] = freq_df['residue_2'].apply(lambda x: converter(x, label_dict))

    freq_df = freq_df.sort_values(by=['contact_frequency'], ascending=False)

    freq_df = freq_df[['residue_1_bw','residue_2_bw', 'residue_1', 'residue_2', 'contact_frequency']]

    return freq_df

if __name__ == '__main__':
    # label_path = Path('/home/hsbyeon1/proj/ago_class/data/water_md/md_nolig_interaction/label/b1ar_melga.tsv')
    # freq_path = Path('/home/hsbyeon1/proj/ago_class/data/water_md/md_nolig_interaction/frequency/2ycw_A_0_frequencies_wb.tsv')
    # output_path = Path('/home/hsbyeon1/proj/ago_class/data/water_md/md_nolig_interaction/frequency_bw/2ycw_A_0_frequencies_wb_bw.tsv')
    parser = argparse.ArgumentParser(description="converts frequency table to bw numbering")
    parser.add_argument(
        "-i", type=Path, help="Path to the frequency table", required=True
    )
    parser.add_argument(
        "-l", type=Path, help="Path to the bw label mapping table", required=True
    )
    parser.add_argument(
        "-o", type=Path, help="Path to the output tsv file", required=True
    )
    parser.add_argument(
        "-c", type=int, help="Minimum resno distance", default=4
    )
    args = parser.parse_args()

    label_path = args.l
    freq_path = args.i
    output_path = args.o
    cutoff = args.c

    freq_df=main(label_path, freq_path, output_path, cutoff)
    freq_df.to_csv(output_path, sep='\t', index=False)