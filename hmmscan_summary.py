import os
import re
from tqdm import tqdm
from argparse import ArgumentParser
from typing import Tuple
import pandas as pd
from collections import defaultdict
from jinja2 import Environment, FileSystemLoader
from datetime import datetime
from Bio import SeqIO

__author__ = "Sara Wattanasombat"
__date__ = "2025-01-21"
__version__ = "1.0.0"
__project__ = "Phosphatase scanner result summary and motif finder"

# ======= Mapping "General HMM Name" to the correct "Fold". =========
# E-value cutoff for HMM scan confidence.
E_VALUE_CUTOFF = 1e-6
# Map "HMM Name" to "Fold family".
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

MOTIF_PATTERNS = {
    "CC1":   [
        {
            "pattern": r"HC[A-Z]{5}R",
            "match_type": "exact",
            "motif": "CC1",
            "description": "Canonical CC1 motif",
        },
        {
            "pattern": r"HC[A-Z]{5}[A-Z]",
            "match_type": "relaxed",
            "motif": "CC1",
            "description": "Conserved canonical HC motif with mutated C-terminal arginine",
        },
        {
            "pattern": r"[A-Z]C[A-Z]{5}R",
            "match_type": "relaxed",
            "motif": "CC1",
            "description": "His mutated at position 1",
        },
        {
            "pattern": r"H[A-Z]{6}R", # only R conserved
            "match_type": "degenerate",
            "motif": "CC1",
            "description": "Mutated C at position 2",
        },
    ],
    "CC2": [
        {
            "pattern": r"(V|I)C[A-Z]{5}R",
            "match_type": "exact",
            "motif": "CC2",
            "description": "Canonical CC2 motif",
        },
        {
            "pattern": r"(V|I)C[A-Z]{5}K",
            "match_type": "degenerate",
            "motif": "CC2",
            "description": "Arg→Lys substitution at position 8",
        },
    ],
    "CC3": [
        {
            "pattern": r"HC[A-Z]{5}R",
            "match_type": "exact",
            "motif": "CC3",
            "description": "Canonical CC3 motif",
        },
        {
            "pattern": r"HC[A-Z]{5}[A-Z]",
            "match_type": "relaxed",
            "motif": "CC3",
            "description": "Conserved canonical HC motif with mutated C-terminal arginine",
        },
        {
            "pattern": r"[A-Z]C[A-Z]{5}R",
            "match_type": "relaxed",
            "motif": "CC3",
            "description": "N-terminal residue mutated",
        },
        {
            "pattern": r"H[A-Z]{6}R", # only R conserved
            "match_type": "degenerate",
            "motif": "CC3",
            "description": "Mutated C at position 2",
        },
        {
            "pattern": r"[A-Z]C[A-Z]{5}K", # only C conserved
            "match_type": "degenerate",
            "motif": "CC3",
            "description": "Conserved C with Arg→Lys substitution at position 8",
        },
    ],
    "PPM":   [
        {
            "pattern": r"(([NR][A-Z]{4}D).*?(D[ND]))",
            "match_type": "exact",
            "motif": "PPM",
            "description": "Canonical PPM connected motif",
            "connected": True,
            "part": "full",
        },
        {
            "pattern": r"((DG[A-Z]{2}G).*?(D[ND]))",
            "match_type": "relaxed",
            "motif": "PPM",
            "description": "Alternative N-terminal PPM anchor",
            "connected": True,
            "part": "full",
        },
        {
            "pattern": r"[NR][A-Z]{4}D",
            "match_type": "anchor",
            "motif": "PPM",
            "description": "PPM N-terminal anchor only",
            "connected": False,
            "part": "A",
        },
        {
            "pattern": r"DG[A-Z]{2}G",
            "match_type": "anchor",
            "motif": "PPM",
            "description": "PPM alternative N-terminal anchor",
            "connected": False,
            "part": "A",
        },
        {
            "pattern": r"(([ND][A-Z]{4}D).*?(D[ND]))",
            "match_type": "remnant",
            "motif": "PPM",
            "description": "PPM fold remnant: acidic anchors retained with degraded chemistry",
            "connected": True,
            "part": "full",
        },
        {
            "pattern": r"DGL[WFY][DE]",
            "match_type": "remnant",
            "motif": "PPM",
            "description": "Degenerate PPM signature",
            "connected": False,
            "part": "anchor",
        },
    ],
    "PPPL": [
        {
            "pattern": r"GD[STAV]HG",
            "match_type": "exact",
            "motif": "PPPL",
            "description": "Core PPPL catalytic GDxHG motif",
            "connected": False,
            "part": "full",
        },
        {
            "pattern": r"GD[A-Z]HG[A-Z]{2}D",
            "match_type": "exact",
            "motif": "PPPL",
            "description": "Extended PPPL metal-binding motif",
            "connected": False,
            "part": "full",
        },
        {
            "pattern": r"GD[A-Z]VDRG",
            "match_type": "relaxed",
            "motif": "PPPL",
            "description": "Extended PPPL variant",
            "connected": False,
            "part": "full",
        },
        {
            "pattern": r"G[DN][A-Z]HG",
            "match_type": "relaxed",
            "motif": "PPPL",
            "description": "Divergent PPPL catalytic loop",
            "connected": False,
            "part": "core",
        },
        {
            "pattern": r"GN[HQ]E",
            "match_type": "remnant",
            "motif": "PPPL",
            "description": "Degenerate PPPL remnant",
            "connected": False,
            "part": "anchor",
        },
        {
            "pattern": r"G[A-Z][A-Z]HG",
            "match_type": "remnant",
            "motif": "PPPL",
            "description": "Severely degenerate GDxHG motif",
            "connected": False,
            "part": "anchor",
        },
    ],
    "RTR1": [
        {
            "pattern": r"H[A-Z]{2}H[A-Z]{10,30}[DE]",
            "match_type": "canonical",
            "motif": "RTR1",
            "description": "Canonical RTR1 His–His–acid catalytic motif",
            "connected": False,
            "part": "core",
        },
        {
            "pattern": r"H[A-Z]H[A-Z]{5,40}[DE]",
            "match_type": "relaxed",
            "motif": "RTR1",
            "description": "Divergent but likely active RTR1 motif",
            "connected": False,
            "part": "core",
        },
        {
            "pattern": r"[HQ][A-Z]H[A-Z]{10,40}[DE]",
            "match_type": "relaxed",
            "motif": "RTR1",
            "description": "Conservatively substituted RTR1 motif",
            "connected": False,
            "part": "core",
        },
        {
            "pattern": r"H[A-Z]{10,50}[DE]",
            "match_type": "remnant",
            "motif": "RTR1",
            "description": "Degenerate RTR1 remnant",
            "connected": False,
            "part": "anchor",
        },
    ],
    "HAD":   [
        {
            "pattern": r"D[A-Z]D[A-Z]T",
            "match_type": "exact",
            "motif": "HAD",
            "description": "Canonical HAD catalytic motif",
            "connected": False,
            "part": "full",
        },
    ],
    "HP":    [
        {
            "pattern": r"((RHG[A-Z]R[A-Z]P).*?(HD))",
            "match_type": "exact",
            "motif": "HP",
            "description": "Canonical HP connected motif",
            "connected": True,
            "part": "full",
        },
        {
            "pattern": r"((RHG[A-Z]R[A-Z]P).*?(H[TQ]))",
            "match_type": "degenerate",
            "motif": "HP",
            "description": "Degenerate HP remnant H-G core preserved with degraded downstream motif (His followed by Thr / Gln)",
            "connected": True,
            "part": "full",
        },
        {
            "pattern": r"RHG[A-Z]R[A-Z]P",
            "match_type": "remnent",
            "motif": "HP",
            "description": "HP N-terminal anchor only",
            "connected": False,
            "part": "A",
        },
        {
            "pattern": r"R?HG[A-Z]{3}P?",
            "match_type": "remnant",
            "motif": "HP",
            "description": "Degenerate HP remnant without downstream motif",
        },
    ],
    "Alk":   [
        {
            "pattern": r"[IV][A-Z]DS[GAS][GASC][GAST][GA]T",
            "match_type": "exact",
            "motif": "Alk",
            "description": "Canonical Alk catalytic motif",
            "connected": False,
            "part": "full",
        },
    ],
}

def get_args() -> Tuple[str | os.PathLike, str | os.PathLike, str | os.PathLike]:
    parser = ArgumentParser(description="Parse HMMER scan inputs")

    parser.add_argument(
        "-t", "--title",
        dest="title",
        default="Untitled",
        help="Title of the run i.e. 'PhosphoScan_2025_Run01'"
    )
    parser.add_argument(
        "-tbl", "--hmmscan-tblout",
        dest="args_hmmscan_tblout",
        required=True,
        help="Path to hmmscan --tblout file"
    )
    parser.add_argument(
        "-fa", "--fasta-file",
        dest="args_fasta_file",
        required=True,
        help="Path to input FASTA file"
    )
    parser.add_argument(
        "-o", "--outdir",
        dest="outdir",
        required=True,
        help="Output directory"
    )

    args = parser.parse_args()

    return (
        args.args_hmmscan_tblout,
        args.args_fasta_file,
        args.outdir,
        args.title
    )

def read_hmmscan_tblout(filepath: str| os.PathLike):
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
        usecols=range(len(colnames))
    )
    _df[numeric_cols] = _df[numeric_cols].apply(pd.to_numeric, errors="coerce")
    return _df

def classify_hmmscan(hmmscan_tbldf: pd.DataFrame, min_evalue=1e-6):
    """
    Classify HMMER hmmscan hits into confidence categories.

    This function takes the tabular output from `read_hmmscan_tblout`, filters out
    hits with very weak E-values, and assigns each hit to one of the following
    categories based on heuristic thresholds:

        - "High-confidence"
        - "Probable"
        - "Weak"
        - "Fragmentary"
        - "Ambiguous"
        - "No match"

    Parameters
    ----------
    hmmscan_tbldf : pandas.DataFrame
        DataFrame containing hmmscan results, as formatted by
        `read_hmmscan_tblout`. Must include the following columns:
        ``["full_seq_Evalue", "best_dom_Evalue", "best_dom_score", "best_dom_bias",
        "exp", "clu", "ov", "env", "rep", "inc"]``.
    min_evalue : float, optional
        Minimum full-sequence E-value threshold for filtering hits (default = 1e-6).

    Returns
    -------
    pandas.DataFrame
        A copy of the input DataFrame with an additional column ``"confidence"``
        containing the assigned category for each hit.

    Raises
    ------
    BaseException
        If one or more of the required columns are missing from the input DataFrame.

    Notes
    -----
    The classification scheme is heuristic and may need to be adjusted for
    specific datasets or applications.
    """
    _df = hmmscan_tbldf.copy()

    _df = _df.loc[_df['full_seq_Evalue'] < min_evalue]
    expected_cols = [
        "full_seq_Evalue",
        "best_dom_Evalue", "best_dom_score", "best_dom_bias",
        "exp", "clu", "ov", "env", "rep", "inc",
    ]

    if set(expected_cols).intersection(set(_df.columns)) != set(expected_cols):
        raise BaseException('classify_hmmscan: Expected column(s) not found in dataframe.')

    def classify(row):
        full_e = row["full_seq_Evalue"]
        dom_e = row["best_dom_Evalue"]
        score = row["best_dom_score"]
        bias = row["best_dom_bias"]
        inc = row["inc"]
        rep = row["rep"]
        exp = row["exp"]
        env = row["env"]
        clu = row["clu"]
        ov = row["ov"]

        # Decision tree
        if full_e <= 1e-5:
            if inc >= 1:
                if dom_e <= 1e-3 and bias < 0.1 * score:
                    return "High-confidence"
                else:
                    return "Probable"
            else:
                if dom_e <= 1e-3 and rep >= 1:
                    return "Probable"
                elif exp >= 0.5 and env > 0:
                    return "Weak"
                else:
                    return "Fragmentary"
        else:
            if dom_e <= 1e-3 and inc >= 1:
                return "Probable"
            elif clu > 0 or ov > 0:
                return "Ambiguous"
            elif rep >= 1 and dom_e <= 1e-2:
                return "Weak"
            else:
                return "No match"

    _df["confidence"] = _df.apply(classify, axis=1)
    return _df

def find_motif(sequence, pattern_list):
    for i, motif in enumerate(pattern_list):
        matches = []
        pattern = motif.get('pattern')
        for m in re.finditer(pattern, str(sequence)):
            matches.append({
                "pattern": str(pattern),
                "match_type": motif.get('match_type'),
                "description": motif.get('description'),
                "match": m.group(),       # the actual substring
                "start": m.start(),       # 0-based start index
                "end": m.end()            # 0-based end index (exclusive)
            })
        if len(matches) >= 1:
            return matches
            break
        if i == len(pattern_list)-1:
            return None

def safe_int(x):
    if pd.isna(x):
        return None
    return int(x)

def highlight_sequence(sequence, spans):
    seq = list(sequence)
    
    for start, end in spans:
        start -= 1  # convert to 0-based
        end -= 1
        seq[start] = f"<span class='seq-highlight'>{seq[start]}"
        seq[end] = f"{seq[end]}</span>"

    return "".join(seq)

def render_report_html(template_dir, title, data):
    env = Environment(
        loader=FileSystemLoader(template_dir)
    )
    template = env.get_template("template.j2")

    grouped = (
        data.groupby(['motif', 'match_type'])
        .size()
        .reset_index(name='value')
    )

    fold_family_counts = {
        "name": "root",
        "children": []
    }

    motif_map = defaultdict(list)
    for row in grouped.itertuples(index=False):
        motif_map[row.motif].append({
            "name": row.match_type.capitalize(),
            "value": int(row.value)
        })

    for motif, children in motif_map.items():
        fold_family_counts["children"].append({
            "name": motif,
            "children": children
        })

    grouped_sequences = defaultdict(lambda: {
        "sequence_id": None,
        "sequence": None,
        "motif": None,
        "HMM confidence": None,
        "e-value": None,
        "sequence_html": None,
        "matches": [],
        "n_matches": 0
    })

    for _, row in data.iterrows():
        key = (row["sequence_id"], row["motif"])

        entry = grouped_sequences[key]
        entry["sequence_id"] = row["sequence_id"]
        entry["sequence"] = row["sequence"]
        entry["motif"] = row["motif"]
        entry["HMM confidence"] = row["HMM confidence"]
        entry["e-value"] = row["e-value"]

        entry["matches"].append({
            "pattern": row["pattern"],
            "match": row["match"],
            "start": safe_int(row["start"]),
            "end": safe_int(row["end"]),
            "match_type": row["match_type"].capitalize(),
            "description": row["description"]
        })
        if safe_int(row["start"]) is not None and safe_int(row["end"]) is not None:
            entry['n_matches'] += 1

    for entry in grouped_sequences.values():
        spans = [
            (m["start"], m["end"])
            for m in entry["matches"]
            if m["match_type"] is not None and m["end"] is not None
        ]
        entry["sequence_html"] = highlight_sequence(
            entry["sequence"],
            spans
        )

    sequences = list(grouped_sequences.values())

    html = template.render(
        project_name=title,
        timestamp=datetime.now(),
        database="Phosphatase Scanner DB",
        sequences=sequences,
        fold_family_counts=fold_family_counts
    )

    with open("report_.html", "w", encoding="utf-8") as f:
        f.write(html)

if __name__ == '__main__':

    hmmscan_tblout_file, fasta_file, outdir, project_title = get_args()

    df = read_hmmscan_tblout(hmmscan_tblout_file)

    # 2. HMM hit filtering and classification
    df = df.loc[df['full_seq_Evalue'] < E_VALUE_CUTOFF].reset_index(drop=True)
    df = classify_hmmscan(df)
    best_hits = (
        df.sort_values(
            ['query_name', 'full_seq_Evalue', 'full_seq_score'],
            ascending=[True, True, False]
        )
        .groupby('query_name', as_index=False)
        .head(1)
    )

    best_hits_map = best_hits.set_index('query_name').to_dict('index')

    sequences = SeqIO.parse(fasta_file, 'fasta')
    matched_sequences = []

    # 3. Motif pattern matching
    all_pattern_matches = []
    for seq in tqdm(sequences, desc="Processing sequences", mininterval=0.5):
        data = best_hits_map.get(seq.id)
        if data:
            fold_family = FOLD_DICT.get(data["target_name"], 'unknown')
            seq.description += f' fold_family={fold_family}'
            matched_sequences.append(seq)
            if fold_family != 'unknown':
                motif_pattern_match = find_motif(seq.seq, MOTIF_PATTERNS.get(fold_family, []))
                if motif_pattern_match:
                    for match in motif_pattern_match:
                        all_pattern_matches.append({
                            "sequence_id": seq.id,
                            "sequence": str(seq.seq),
                            "motif": fold_family,
                            "HMM confidence": data.get('confidence'),
                            "e-value": f"{data.get('full_seq_Evalue'):.2e}",
                            "Motif confidence e-value": f"{data.get('full_seq_Evalue'):.2e}",
                            **match
                        })
                else:
                    all_pattern_matches.append({
                        "sequence_id": seq.id,
                        "sequence": str(seq.seq),
                        "motif": fold_family,
                        "HMM confidence": data.get('confidence'),
                        "e-value": f"{data.get('full_seq_Evalue'):.2e}",
                        "pattern": None,
                        "match_type": "No Match",
                        "match": None,
                        "description": f"No evidence of canonical or remnent of {fold_family}",
                        "start": None,
                        "end": None
                    })

    # 4. Summerize result
    n_matched_sequnece = SeqIO.write(matched_sequences, f'{outdir}/match_sequences.fa', 'fasta')
    print(f'Saved {n_matched_sequnece} sequences to {outdir}/match_sequences.fa')
    pattern_match_summary = pd.DataFrame(all_pattern_matches)
    pattern_match_summary.to_csv(f'{outdir}/match_summary.tsv', sep='\t', index=None)

    render_report_html(f'{os.path.dirname(os.path.realpath(__file__))}/report/', project_title, pattern_match_summary)
