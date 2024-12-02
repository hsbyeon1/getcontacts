# reads residue table and generates label tsv
# ex.
# "A:ALA:4    A4  white\n",
# "A:CYS:5    C5  yellow\n",
# "A:TRP:6    W6\n"
# Reasidue table can be downloaded from GPCRdb
import pandas as pd
from pathlib import Path
import argparse

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
    "-": "-",
}

RESIDUES_1_TO_3 = {v: k for k, v in RESIDUES_3_TO_1.items()}

COLOR_TABLE = {
    1: "red",
    2: "orange",
    3: "yellow",
    4: "lightgreen",
    5: "green",
    6: "blue",
    7: "navy",
    8: "purple",
}


def convert_res(resinfo: str, offset: int) -> tuple[str, int]:
    """
    converts 'R20' into ('ARG', 20 + offset), and '-' into ('-', -1)
    """
    # offset b1ar_melga : -32

    if resinfo == "-":
        return ("-", -1)
    else:
        return (RESIDUES_1_TO_3[resinfo[0]], int(resinfo[1:]) + offset)


def fetch_receptor_data(
    df_path: Path, receptor_name: str, chainid: str = "A", offset: int = 0
) -> list[tuple[str, str, str]]:
    """
    Fetches the receptor data from the residue table, and returns a list of a form
    [("A:ALA:4", "1x28"),
    ("A:CYS:5", "1x29"), ...]

    Args:
    df_path: Path to the residue table
    receptor_name: Name of the receptor, e.g. "5-HT1A receptor Human"
    chainid: Chain ID of the receptor, default is "A"
    offset: Offset to be added to the residue number, default is 0

    Returns:
    list[tuple[str,str,str]]: List of tuples of the form [("A:ALA:4", "1x28", 'red'), ...]
    """

    # columns : ["GPCRdb(A)", "5-HT1A receptor Human", "", ...]
    df = pd.read_csv(df_path, index_col=None)

    # get the column index of the receptor
    receptor_idx = df.columns.get_loc(receptor_name)

    # get the receptor data
    # 0 for the first col (GPCRdb numberings)
    receptor_df = df.iloc[:, [0, receptor_idx]]
    receptor_df.reset_index(drop=True, inplace=True)

    receptor_df.iloc[:, 0] = receptor_df.iloc[:, 0].apply(lambda x: x.split()[0])

    # convert into a zip object
    receptor_data_1:list[tuple[str,str]] = [
        tup for tup in receptor_df.itertuples(index=False, name=None) if tup[1] != "-"
    ]  # ("1x28", "T32")
    receptor_data_3 = list(
        map(lambda x: (x[0], convert_res(x[1], offset)), receptor_data_1)
    )  # ("1x28", ("THR", 32))

    final_receptor_data = list(
        map(
            lambda x: (
                f"{chainid}:{x[1][0]}:{x[1][1]}",
                x[0],
                COLOR_TABLE.get(int(x[0].split("x")[0]), "white"),
            ),
            receptor_data_3,
        )
    )  # ("A:THR:32", "1x28", "red")

    return final_receptor_data


def main(
    input_df_path,
    output_csv_path,
    receptor_name: str = "5-HT1A receptor Human",
    chainid: str = "A",
    offset: int = 0,
):
    # input_df_path=Path("/home/hsbyeon1/proj/ago_class/data/gpcrdb_data/residue_table.csv" )
    # output_csv_path=Path("/home/hsbyeon1/proj/ago_class/data/gpcrdb_data/5-HT1A_human_label.tsv")
    receptor_data = fetch_receptor_data(input_df_path, receptor_name, chainid, offset)

    with open(output_csv_path, "w") as f:
        for tup in receptor_data:
            t = "\t".join(tup) + "\n"
            f.write(t)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate label tsv file")
    parser.add_argument(
        "-i", type=Path, help="Path to the residue table", required=True
    )
    parser.add_argument(
        "-o", type=Path, help="Path to the output label tsv file", required=True
    )
    parser.add_argument("-n", type=str, help="Name of the receptor", required=True)
    parser.add_argument("-c", type=str, help="Chain ID of the receptor", default="A")
    parser.add_argument(
        "--offset", type=int, help="Offset to be added to the residue number", default=0
    )

    args = parser.parse_args()

    main(args.i, args.o, args.n, args.c, args.offset)
