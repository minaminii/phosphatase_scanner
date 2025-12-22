from argparse import ArgumentParser
from os import PathLike
from typing import Tuple
import pandas as pd
from collections import Counter
from jinja2 import Environment, FileSystemLoader
from Bio import SeqIO

# ======= Mapping "General HMM Name" to the correct "Fold". =========
FOLD_DICT = {
        'Alk_phosphatase': 'ALK',
        'CPDc_SM00577': 'HAD',
        'Chronophin_1': 'HAD',
        'Chronophin_2': 'HAD',
        'DSPc': 'CC1',
        'DSPc_SM00195': 'CC1',
        'DUF442': 'Unknown',
        'His_Phos_1': 'HP',
        'His_Phos_2': 'HP',
        'Init_tRNA_PT': 'CC1',
        'LMWPc': 'CC2',
        'LMWPc_SSU72': 'CC2',
        'LMWPc_SSU72_ArsC': 'CC2',
        'MDP_1': 'HAD',
        'MDP_2': 'HAD',
        'Myotub-related': 'CC1',
        'NIF': 'HAD',
        'PP2Ac_SM00156': 'PPPL',
        'PP2C': 'PPM',
        'PP2C_2': 'PPM',
        'PP2C_SIG_SM00331': 'PPM',
        'PP2Cc_SM00332': 'PPM',
        'PPM_PTC7': 'PPM',
        'PPM_PTC7_SpoIIE': 'PPM',
        'PPM_PTC7_motif': 'PPM',
        'PPP_1': 'PPPL',
        'PPP_2': 'PPPL',
        'PPP_3': 'PPPL',
        'PPP_4': 'PPPL',
        'PPP_5_Dcr2': 'PPPL',
        'PTEN_1': 'CC1',
        'PTPc_DSPc_SM00012': 'CC1',
        'PTPc_motif_SM00404': 'CC1',
        'PTPlike_phytase': 'CC1',
        'RHOD_SM00450': 'CC3',
        'RPAP2_Rtr1': 'RTR1',
        'Ssu72': 'CC2',
        'Y_phosphatase': 'CC1',
        'Y_phosphatase2': 'CC1',
        'Y_phosphatase3': 'CC1',
        'PPM_PTC7_SpoIIE_motif': 'PPM',
        'SpoIIE': 'PPM',
    }

E_VALUE_CUTOFF = 1e-6

def get_args() -> Tuple[str | PathLike, str | PathLike, str | PathLike]:
    parser = ArgumentParser(description="Parse HMMER scan inputs")

    parser.add_argument(
        "-t", "--hmmscan-tblout",
        dest="args_hmmscan_tblout",
        required=True,
        help="Path to hmmscan --tblout file"
    )
    parser.add_argument(
        "-f", "--fasta-file",
        dest="args_fasta_file",
        required=True,
        help="Path to input FASTA file"
    )
    parser.add_argument(
        "-o", "--match-output-file",
        dest="args_match_output_file",
        required=True,
        help="Path to output match file"
    )

    args = parser.parse_args()

    return (
        args.args_hmmscan_tblout,
        args.args_fasta_file,
        args.args_match_output_file,
    )

def read_hmmscan_tblout(filepath: str| PathLike):
    """
    Read an hmmscan --tblout file into a pandas DataFrame.
    
    Skips comment lines (#) and parses the standard 18-column HMMER output.
    Returns a DataFrame with column names matching the hmmscan manual.
    """
    colnames = [
        "target_name", "target_accession",
        "query_name", "query_accession",
        "full_seq_Evalue", "full_seq_score", "full_seq_bias",
        "best_dom_Evalue", "best_dom_score", "best_dom_bias",
        "exp", "reg", "clu", "ov", "env", "dom", "rep", "inc",
        "description_of_target"
    ]
    numeric_cols = [
        "full_seq_Evalue", "full_seq_score", "full_seq_bias",
        "best_dom_Evalue", "best_dom_score", "best_dom_bias",
        "exp", "reg", "clu", "ov", "env", "dom", "rep", "inc"
    ]

    _df = pd.read_csv(
        filepath,
        comment="#",
        sep='\\s+',
        names=colnames,
        usecols=range(len(colnames))  # ignore any trailing extra whitespace
    )
    # convert numeric columns to floats/ints where appropriate
    _df[numeric_cols] = _df[numeric_cols].apply(pd.to_numeric, errors="coerce")

    return _df

if __name__ == '__main__':

    hmmscan_tblout_file, fasta_file, match_output_file = get_args()

    df = read_hmmscan_tblout(hmmscan_tblout_file)
    df = df.loc[df['full_seq_Evalue'] < E_VALUE_CUTOFF].reset_index(drop=True)

    df.insert(1, 'fold_family', df['target_name'].map(FOLD_DICT).fillna('unknown'))

    # Get best hit for each query.
    best_hits = df.groupby(
        'query_name', group_keys=True
    ).apply(
        lambda g:g.sort_values(['full_seq_Evalue', 'full_seq_score'],ascending=True)
        .head(1), include_groups=False
    )
    best_hits = best_hits.reset_index().drop('level_1', axis=1)

    sequences = SeqIO.parse(fasta_file, 'fasta')
    matched_sequences = []
    for seq in sequences:
        # For each sequence in fasta file, find a match hmmscan result.
        data = best_hits.loc[best_hits['query_name'] == seq.id]
        if len(data) > 0:
            data = data.to_dict('records')[0]
            # Modify the sequnece id to include predition result
            seq.description += f' fold_family={data.get("fold_family")}'
            matched_sequences.append(seq)

    print(matched_sequences[0])
    SeqIO.write(matched_sequences, match_output_file, 'fasta')

    ###################

    env = Environment(loader=FileSystemLoader("/mnt/sdb/SOFTWARE/phosphatase-scanner/scanner-v1.1/phosphatase_scanner/report/"))
    template = env.get_template("template.j2")

    fold_families = []

    for r in matched_sequences:
        # Split description into tokens
        parts = r.description.split()

        # Find the token starting with "fold_family="
        family_token = next((p for p in parts if p.startswith("fold_family=")), None)

        if family_token:
            family = family_token.split("=")[1]
            fold_families.append(family)
        else:
            fold_families.append("Unknown")

    fold_family_counts = Counter(fold_families)

    html = template.render(
        project_name="PhosphoScan_2025_Run01",
        timestamp="2025-12-02 14:33:55",
        database="Custom Phosphatase DB v3.4",
        sequences=matched_sequences,
        fold_family_counts=fold_family_counts
    )

    with open("report_.html", "w", encoding='utf-8') as f:
        f.write(html)
