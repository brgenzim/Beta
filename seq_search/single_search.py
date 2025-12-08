#!/usr/bin/env python3


from __future__ import annotations

import argparse
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple


def read_fasta_sequence(path: Path) -> Tuple[str, str]:
    qid = ""
    parts: List[str] = []
    with path.open() as inf:
        for line in inf:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                qid = line[1:].strip()
            else:
                parts.append(line.strip())
    return qid, "".join(parts)


def run_psiblast(psi_blast_binary: Path, input_fas: Path, pdb_blast_db: Path) -> None:

    cmd = [
        str(psi_blast_binary),
        "-query",
        str(input_fas),
        "-db",
        str(pdb_blast_db),
        "-evalue",
        "0.0001",
        "-max_target_seqs=50000",
        "-num_iterations",
        "3",
        "-out_ascii_pssm",
        "query.pssm",
        "-out",
        "query.blast",
    ]
    subprocess.run(cmd, check=True)


def load_pdb_summary(path: Path) -> Dict[str, Tuple[str, str, str]]:
    meta: Dict[str, Tuple[str, str, str]] = {}
    with path.open() as inf:
        for line in inf:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            pdb_id = cols[0]
            date_, method, resolution = cols[1], cols[2], cols[3].strip()
            meta[pdb_id] = (date_, method, resolution)
    return meta


def _extras_for_entry(date_str: str, method: str) -> Tuple[str, str]:
    extraspace = "\t" if ("NMR" in method) else ""
    extramsg = ""

    try:
        year = int(date_str[:4])
        month = int(date_str[5:7])
    except Exception:
        return extraspace, extramsg  # fail-soft

    if year < 2021:
        extramsg = "AF3 Server and AF2 Database may seen this structure"
    elif year == 2021:
        if month < 3:
            extramsg = "AF3 Server and AF2 Database may seen this structure"
        elif month < 10:
            extramsg = "AF3 Server may seen this structure"

    return extraspace, extramsg


def parse_blast_collect(
    blast_path: Path, query_seq: str, pdb_meta: Dict[str, Tuple[str, str, str]]
) -> Dict[str, str]:
    used: Dict[str, str] = {}

    with blast_path.open() as inf:
        while True:
            line = inf.readline()
            if line == "":
                break
            if not line.startswith(">"):
                continue

            hit_id = line[1:7]  # preserve exact behavior

            # Skip 4 lines, then locate the 'Identities' line
            for _ in range(4):
                line = inf.readline()

            if not line.split():
                while True:
                    line = inf.readline()
                    try:
                        if line.split()[0] == "Identities":
                            break
                    except IndexError:
                        pass

            if line.split()[0] != "Identities":
                while True:
                    line = inf.readline()
                    try:
                        if line.split()[0] == "Identities":
                            break
                    except IndexError:
                        pass

            parts = line.split()
            try:
                identity_pct = int(parts[3].replace("(", "").replace("%)", "").replace(",", ""))
                aln_len = int(parts[2].split("/")[1].replace(",", ""))
            except (IndexError, ValueError):
                continue  # fail-soft

            if identity_pct <= 20 or aln_len <= 10:
                continue

            # Parse alignment block
            q_aln = ""
            h_aln = ""
            qs = 1000000
            qe = -100000
            hs = 1000000
            he = -100000
            blanks_in_row = 0

            while True:
                line = inf.readline()
                if line == "":
                    break
                toks = line.split()
                if len(toks) == 0:
                    blanks_in_row += 1
                    if blanks_in_row == 2:
                        break
                    continue
                blanks_in_row = 0

                tag = toks[0]
                if tag == "Query":
                    q_aln += toks[2]
                    try:
                        qs = min(qs, int(toks[1]))
                        qe = max(qe, int(toks[3].strip()))
                    except (ValueError, IndexError):
                        pass
                elif tag == "Sbjct":
                    h_aln += toks[2]
                    try:
                        hs = min(hs, int(toks[1]))
                        he = max(he, int(toks[3].strip()))
                    except (ValueError, IndexError):
                        pass

            # Build spacing exactly like the original
            pre = query_seq[:qs]
            post = query_seq[qe:]
            pre_spaces = " " * max(len(pre) - 1, 0)
            post_spaces = " " * len(post)

            # Remove gaps in hit using query alignment as mask
            h_compact_chars: List[str] = []
            hi = 0
            for qc in q_aln:
                if qc != "-":
                    h_compact_chars.append(h_aln[hi])
                hi += 1
            h_compact = "".join(h_compact_chars)

            pdb_code = hit_id[:4]
            try:
                date_, method, resolution = pdb_meta[pdb_code]  # KeyError if missing (preserved)

                extraspace, extramsg = _extras_for_entry(date_, method)
                res_clean = resolution.replace("unknown", "NA")[:3]

            # Match your exact column ordering/spacing
                line_out = (
                    f"{hit_id}\t{pre_spaces}{h_compact}{post_spaces}"
                    f"\t{qs}\t{qe}\t{hs}\t{he}\t{identity_pct}\t\t{date_}\t{method}"
                    f"{extraspace}\tresolution:{res_clean} A    \t{extramsg}\n"
                )

                used[hit_id] = line_out  # keep last occurrence
            except KeyError:
                line_out = (
                    f"{hit_id}\t{pre_spaces}{h_compact}{post_spaces}"
                    f"\t{qs}\t{qe}\t{hs}\t{he}\t{identity_pct}\t\tNA\tNA"
                    f"\tNA\tNA\n"
                )

    return used


def write_output(out_path: Path, query_seq: str, used: Dict[str, str]) -> None:
    empty = " " * len(query_seq)
    with out_path.open("w") as outf:
        outf.write(
            "ID\t"
            f"{empty}\tq_from\tq_to\th_from\th_to\tidentity\tdate\t\tmethod\t\t\tresolution\t\tnote\n"
        )
        outf.write(f"query\t{query_seq}\n")
        for key in used:
            outf.write(used[key])


def main() -> None:
    parser = argparse.ArgumentParser(
        description=""
    )
    parser.add_argument("input_fas", type=Path, help="Query FASTA file")
    parser.add_argument("pdb_blast_db", type=Path, help="PDB BLAST database")
    parser.add_argument("psi_blast_binary", type=Path, help="Path to psiblast binary")
    parser.add_argument("pdb_summary", type=Path, help="PDB summary TSV")
    parser.add_argument("out", type=Path, help="Output file")
    args = parser.parse_args()

    # Read query sequence
    _, qseq = read_fasta_sequence(args.input_fas)

    # Run PSI-BLAST (produces query.blast and query.pssm)
    run_psiblast(args.psi_blast_binary, args.input_fas, args.pdb_blast_db)

    # Load PDB metadata
    pdb_meta = load_pdb_summary(args.pdb_summary)

    # Parse BLAST and collect output lines by ID
    used = parse_blast_collect(Path("query.blast"), qseq, pdb_meta)

    # Write final output
    write_output(args.out, qseq, used)


if __name__ == "__main__":
    print("python3 search.py <INPUT_FAS(SINGLE)> <BLAST_DB> <PSIB_BLAST COMMAND> <PDB_SUMMARY> <OUTFILE>")
    main()
