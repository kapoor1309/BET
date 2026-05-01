import sys
import pandas as pd
import re
import os

ss_file = sys.argv[1]
output_file = sys.argv[2]
target_aa = sys.argv[3]

pdb_id = os.path.basename(ss_file).split(".")[0]

aa3_to_aa1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

records = []

with open(ss_file) as f:
    for line in f:
        if line.startswith("ASG"):
            parts = line.split()
            if len(parts) < 10:
                continue

            m = re.match(r"\d+", parts[3])
            if not m:
                continue

            try:
                records.append({
                    "res3": parts[1],
                    "chain": parts[2],
                    "resnum": int(m.group()),
                    "idx": int(parts[4]),
                    "ss_code": parts[5],
                    "ss_name": parts[6],
                    "phi": float(parts[7]),
                    "psi": float(parts[8]),
                    "area": float(parts[9]),
                })
            except Exception:
                continue

df = pd.DataFrame(records)

rows = []

if not df.empty:
    for i in range(1, len(df) - 1):
        center = df.iloc[i]
        if center["res3"] != target_aa:
            continue

        prev = df.iloc[i - 1]
        next_ = df.iloc[i + 1]

        tri1 = (
            aa3_to_aa1.get(prev["res3"], "X")
            + aa3_to_aa1.get(center["res3"], "X")
            + aa3_to_aa1.get(next_["res3"], "X")
        )
        ss_tri = prev["ss_code"] + center["ss_code"] + next_["ss_code"]

        pos_str = (
            f"{center['chain']}:{prev['resnum']},"
            f"{center['chain']}:{center['resnum']},"
            f"{center['chain']}:{next_['resnum']}"
        )

        for r in [prev, center, next_]:
            rows.append([
                r["res3"],
                r["chain"],
                r["resnum"],
                r["idx"],
                r["ss_code"],
                r["ss_name"],
                r["phi"],
                r["psi"],
                r["area"],
                aa3_to_aa1.get(r["res3"], "X"),
                tri1,
                ss_tri,
                pdb_id,
                pos_str,
            ])

os.makedirs(os.path.dirname(output_file), exist_ok=True)
out_df = pd.DataFrame(rows)
out_df.to_csv(output_file, sep="\t", index=False, header=False)
