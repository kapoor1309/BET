import gzip
import multiprocessing as mp
import os
import sys
from functools import partial

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from tqdm import tqdm


SIZE_MAP = {
    "K": "Large", "R": "Bulky", "D": "Intermediate", "E": "Large",
    "G": "Tiny", "A": "Tiny", "V": "Small", "L": "Intermediate", "I": "Intermediate",
    "M": "Large", "P": "Small", "S": "Small", "T": "Small", "N": "Intermediate",
    "Q": "Large", "C": "Small", "F": "Bulky", "Y": "Bulky", "W": "Bulky", "H": "Large",
}


def signed_angle_3d(v1, v2, axis, degrees=True):
    #print ([v1, v2])
    v1 = np.array(v1, dtype=float)
    v2 = np.array(v2, dtype=float)
    axis = np.array(axis, dtype=float)

    v1 /= np.linalg.norm(v1)
    v2 /= np.linalg.norm(v2)
    axis /= np.linalg.norm(axis)

    cross = np.cross(v1, v2)
    dot = np.dot(v1, v2)

    angle = np.arctan2(
        np.dot(cross, axis),
        dot
    )

    return np.degrees(angle) if degrees else angle


def load_structure(pdb_id, pdb_dir):
    file_path = os.path.join(pdb_dir, f"{pdb_id}.pdb.gz")
    if not os.path.exists(file_path):
        return None
    parser = PDBParser(QUIET=True)
    try:
        with gzip.open(file_path, "rt") as f:
            return parser.get_structure(pdb_id, f)
    except Exception:
        return None


def get_residue(structure, chain_id, resnum):
    try:
        for model in structure:
            chain = model[chain_id]
            for res in chain:
                if res.id[1] == resnum:
                    return res
    except Exception:
        return None
    return None


def get_ca(res):
    if res and "CA" in res:
        return res["CA"].get_coord()
    return None


def get_centroid(res):
    coords = []
    for atom in res:
        if atom.get_name() not in ["N", "CA", "C", "O"]:
            coords.append(atom.get_coord())
    if not coords:
        return None
    return np.mean(coords, axis=0)


def process_context(filename, contexts_dir, pdb_dir, target_aa):
    """Process one context TSV. Returns (rows, pdb_id_or_None)."""
    path = os.path.join(contexts_dir, filename)

    if os.path.getsize(path) == 0:
        return ([], None)

    try:
        df = pd.read_csv(path, sep="\t", header=None)
    except Exception:
        return ([], None)

    if df.empty:
        return ([], None)

    pdb_id = df.iloc[0, 12]
    structure = load_structure(pdb_id, pdb_dir)
    if structure is None:
        return ([], None)

    rows = []
    for i in range(0, len(df), 3):
        try:
            prev = df.iloc[i]
            center = df.iloc[i + 1]
            next_ = df.iloc[i + 2]
        except Exception:
            continue

        if center[0] != target_aa:
            continue
        if center[11] != "HHH":
            continue

        chain = center[1]
        r_prev = int(prev[2])
        r_cent = int(center[2])
        r_next = int(next_[2])

        res_prev = get_residue(structure, chain, r_prev)
        res_cent = get_residue(structure, chain, r_cent)
        res_next = get_residue(structure, chain, r_next)
        if not res_prev or not res_cent or not res_next:
            continue

        ca_prev = get_ca(res_prev)
        ca_cent = get_ca(res_cent)
        ca_next = get_ca(res_next)
        cen_prev = get_centroid(res_prev)
        cen_cent = get_centroid(res_cent)

        if any(x is None for x in [ca_prev, ca_cent, ca_next, cen_prev, cen_cent]):
            continue

        v1 = cen_prev - ca_prev
        v2 = cen_cent - ca_cent
        axis = ca_cent - ca_prev
        if np.linalg.norm(axis) == 0:
            continue

        try:
            angle = signed_angle_3d(v1, v2, axis)
        except Exception:
            continue

        prev_aa1 = prev[9]
        size_class = SIZE_MAP.get(prev_aa1, "Unknown")
        if size_class == "Unknown":
            continue

        rows.append([pdb_id, prev_aa1, size_class, angle])

    return (rows, pdb_id if rows else None)


def main():
    contexts_dir = sys.argv[1]
    output_file = sys.argv[2]
    target_aa = sys.argv[3]
    pdb_dir = sys.argv[4] if len(sys.argv) > 4 else "pdbs"
    workers = int(sys.argv[5]) if len(sys.argv) > 5 else 12

    files = sorted(f for f in os.listdir(contexts_dir) if f.endswith(".tsv"))

    worker = partial(
        process_context,
        contexts_dir=contexts_dir,
        pdb_dir=pdb_dir,
        target_aa=target_aa,
    )

    all_rows = []
    valid_pdbs = set()

    with mp.Pool(workers) as pool:
        for rows, pdb_id in tqdm(
            pool.imap_unordered(worker, files, chunksize=16),
            total=len(files),
            desc="Processing contexts",
        ):
            if rows:
                all_rows.extend(rows)
            if pdb_id is not None:
                valid_pdbs.add(str(pdb_id))

    out_df = pd.DataFrame(all_rows, columns=["pdb", "left_aa", "size_class", "angle"])
    out_dir = os.path.dirname(output_file) or "."
    os.makedirs(out_dir, exist_ok=True)
    out_df.to_csv(output_file, sep="\t", index=False)

    valid_pdbs_path = os.path.join(out_dir, "valid_pdbs.txt")
    with open(valid_pdbs_path, "w") as f:
        for pdb_id in sorted(valid_pdbs):
            f.write(f"{pdb_id}\n")

    print(f"Saved: {output_file}")
    print(f"Saved: {valid_pdbs_path}")
    print(f"Total angles: {len(out_df)}")
    print(f"Valid PDB files: {len(valid_pdbs)}")


if __name__ == "__main__":
    main()
