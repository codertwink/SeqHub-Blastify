"""
seqtools.py - core sequence analysis functions (no UI code here).

Covers:
  * input parsing (raw sequence or FASTA) and alphabet detection
  * DNA/RNA: GC / AT content, reverse, complement, reverse complement, DNA -> RNA
  * BLAST via NCBI (qblast), XML parsing, hit filtering, ortholog candidates
  * Protein: six-frame translation, ORF detection, amino acid composition
  * Export helpers (FASTA, CSV, text report)
"""

from __future__ import annotations

import re
from collections import Counter
from io import StringIO

import pandas as pd
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

NUCLEOTIDE_ALPHABET = set("ACGTUNRYKMSWBDHV-")  # includes IUPAC ambiguity codes
PROTEIN_ALPHABET = set("ACDEFGHIKLMNPQRSTVWYXBZJUO*")
STOP_CODONS = {"TAA", "TAG", "TGA"}
AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")

BLAST_PROGRAMS = {
    "nucleotide": {
        "program": "blastn",
        "databases": {
            "nt (all GenBank nucleotide)": "nt",
            "refseq_rna (RefSeq transcripts)": "refseq_rna",
            "core_nt (curated subset)": "core_nt",
        },
    },
    "protein": {
        "program": "blastp",
        "databases": {
            "nr (all non-redundant protein)": "nr",
            "refseq_protein (RefSeq)": "refseq_protein",
            "swissprot (curated)": "swissprot",
        },
    },
}


# ---------------------------------------------------------------------------
# Input parsing
# ---------------------------------------------------------------------------

class SequenceInputError(ValueError):
    """Raised when the user input cannot be parsed as a biological sequence."""


def parse_fasta_or_raw(text: str) -> tuple[str, str]:
    """
    Accept either raw sequence text or a single FASTA record.
    Returns (header, sequence). Header is "" for raw input.
    Whitespace, digits and line breaks are stripped; the result is uppercased.
    """
    if not text or not text.strip():
        raise SequenceInputError("The input is empty.")

    header = ""
    lines = text.strip().splitlines()
    body_lines = []
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header:
                # Only the first record is analysed; the rest are ignored.
                break
            header = line[1:].strip()
        elif line.startswith(";"):  # old-style FASTA comment
            continue
        else:
            body_lines.append(line)

    raw = "".join(body_lines)
    seq = re.sub(r"[\s\d]", "", raw).upper()
    if not seq:
        raise SequenceInputError("No sequence characters found after the header.")
    return header, seq


def detect_sequence_type(seq: str) -> str:
    """
    Returns "DNA", "RNA", or "protein".
    Raises SequenceInputError if the characters match neither alphabet.
    """
    chars = set(seq)
    if chars <= NUCLEOTIDE_ALPHABET:
        if "U" in chars and "T" not in chars:
            return "RNA"
        return "DNA"
    if chars <= PROTEIN_ALPHABET:
        return "protein"
    bad = sorted(chars - PROTEIN_ALPHABET)
    raise SequenceInputError(
        f"Unrecognised characters in sequence: {''.join(bad)}. "
        "Expected DNA (ACGT), RNA (ACGU) or a protein sequence."
    )


# ---------------------------------------------------------------------------
# Nucleotide analysis
# ---------------------------------------------------------------------------

def nucleotide_counts(seq: str) -> dict:
    """Counts of A, C, G, T (or U for RNA) and ambiguous/other characters."""
    counts = Counter(seq)
    base_keys = ["A", "C", "G", "T"] if "U" not in seq else ["A", "C", "G", "U"]
    result = {b: counts.get(b, 0) for b in base_keys}
    result["Other"] = len(seq) - sum(result.values())
    return result


def gc_content(seq: str) -> float:
    """GC% computed over unambiguous bases (A, C, G, T/U). Returns 0.0 if none."""
    counts = nucleotide_counts(seq)
    unamb = sum(counts[b] for b in ("A", "C", "G", "T", "U") if b in counts)
    if unamb == 0:
        return 0.0
    return 100.0 * (counts.get("G", 0) + counts.get("C", 0)) / unamb


def at_content(seq: str) -> float:
    """AT% (A + T, or A + U for RNA) over unambiguous bases. Equals 100 - GC%."""
    counts = nucleotide_counts(seq)
    unamb = sum(counts[b] for b in ("A", "C", "G", "T", "U") if b in counts)
    if unamb == 0:
        return 0.0
    return 100.0 * (counts.get("A", 0) + counts.get("T", 0) + counts.get("U", 0)) / unamb


def sliding_gc(seq: str, window: int = 100, step: int = 10) -> pd.DataFrame:
    """GC% along the sequence. Returns a DataFrame with position (window centre) and gc."""
    if len(seq) < window:
        window = max(1, len(seq))
    rows = []
    for start in range(0, len(seq) - window + 1, step):
        chunk = seq[start:start + window]
        rows.append({"position": start + window // 2, "gc": gc_content(chunk)})
    return pd.DataFrame(rows)


def reverse(seq: str) -> str:
    return seq[::-1]


def complement(seq: str) -> str:
    """Complement without reversing. Uses Biopython so IUPAC codes are handled."""
    return str(Seq(seq).complement())


def reverse_complement(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


def dna_to_rna(seq: str) -> str:
    """Transcribe DNA to RNA (T -> U). Leaves RNA unchanged."""
    return seq.replace("T", "U")


def transformations(seq: str) -> dict:
    """All nucleotide transformations shown in the dashboard, keyed by label."""
    dna = seq.replace("U", "T") if "U" in seq and "T" not in seq else seq
    return {
        "Original": seq,
        "Reverse": reverse(dna),
        "Complement": complement(dna),
        "Reverse complement": reverse_complement(dna),
        "DNA to RNA (transcript)": dna_to_rna(dna),
    }


# ---------------------------------------------------------------------------
# Protein analysis: six-frame translation, ORFs, composition
# ---------------------------------------------------------------------------

def _translate(nt: str) -> str:
    """Translate a nucleotide string (any length) to protein, trimming partial codons."""
    usable = len(nt) - (len(nt) % 3)
    if usable == 0:
        return ""
    return str(Seq(nt[:usable]).translate(table=1))


def six_frame_translation(seq: str) -> dict:
    """
    Translate all six reading frames of a DNA sequence.
    Keys: "+1", "+2", "+3" (forward) and "-1", "-2", "-3" (reverse complement).
    Stop codons are shown as '*'.
    """
    dna = seq.replace("U", "T")
    rc = reverse_complement(dna)
    frames = {}
    for offset in range(3):
        frames[f"+{offset + 1}"] = _translate(dna[offset:])
    for offset in range(3):
        frames[f"-{offset + 1}"] = _translate(rc[offset:])
    return frames


def find_orfs(seq: str, min_aa: int = 30, start_codons=("ATG",)) -> pd.DataFrame:
    """
    Find open reading frames on both strands.

    An ORF starts at the first start codon after a stop (so the longest ORF
    in each stretch is reported) and ends at the next stop codon. ORFs with
    no stop codon before the sequence ends are skipped.

    Coordinates are 1-based and always refer to the forward (input) strand.
    """
    dna = seq.replace("U", "T")
    L = len(dna)
    rc = reverse_complement(dna)
    records = []

    for strand, s in (("+", dna), ("-", rc)):
        for frame in range(3):
            codons = [s[i:i + 3] for i in range(frame, L - 2, 3)]
            start_idx = None
            for idx, codon in enumerate(codons):
                if start_idx is None:
                    if codon in start_codons:
                        start_idx = idx
                    continue
                if codon in STOP_CODONS:
                    nt_start = frame + 3 * start_idx          # 0-based on strand s
                    nt_end = frame + 3 * (idx + 1)            # exclusive, includes stop
                    n_aa = (nt_end - nt_start) // 3 - 1       # amino acids, stop excluded
                    if n_aa >= min_aa:
                        protein = _translate(s[nt_start:nt_end - 3])
                        if strand == "+":
                            f_start, f_end = nt_start + 1, nt_end
                        else:
                            # Position p on the reverse complement is L - 1 - p on the
                            # forward strand, so the span flips. Stop codon included.
                            f_start = L - nt_end + 1
                            f_end = L - nt_start
                        records.append({
                            "strand": strand,
                            "frame": f"{strand}{frame + 1}",
                            "start": f_start,
                            "end": f_end,
                            "length_nt": nt_end - nt_start,
                            "length_aa": n_aa,
                            "protein": protein,
                        })
                    start_idx = None

    columns = ["strand", "frame", "start", "end", "length_nt", "length_aa", "protein"]
    df = pd.DataFrame(records, columns=columns)
    if not df.empty:
        df = df.sort_values("length_aa", ascending=False).reset_index(drop=True)
        df.insert(0, "ORF", [f"ORF_{i + 1}" for i in range(len(df))])
    return df


def aa_composition(protein: str) -> pd.DataFrame:
    """Amino acid counts and percentages (stop codons and unknowns excluded)."""
    clean = protein.replace("*", "")
    total = len(clean)
    counts = Counter(clean)
    rows = []
    for aa in AMINO_ACIDS:
        n = counts.get(aa, 0)
        rows.append({
            "amino_acid": aa,
            "count": n,
            "percent": round(100.0 * n / total, 2) if total else 0.0,
        })
    return pd.DataFrame(rows)


def protein_properties(protein: str) -> dict:
    """Molecular weight, isoelectric point, GRAVY, and instability index."""
    clean = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", protein.upper())
    if not clean:
        return {}
    pa = ProteinAnalysis(clean)
    return {
        "Length (aa)": len(clean),
        "Molecular weight (Da)": round(pa.molecular_weight(), 2),
        "Isoelectric point (pI)": round(pa.isoelectric_point(), 2),
        "GRAVY (hydropathy)": round(pa.gravy(), 3),
        "Instability index": round(pa.instability_index(), 2),
    }


# ---------------------------------------------------------------------------
# BLAST
# ---------------------------------------------------------------------------

def run_blast(seq: str, program: str, database: str, evalue: float = 10.0,
              hitlist_size: int = 50, entrez_query: str | None = None) -> str:
    """
    Submit a sequence to NCBI BLAST (web service) and return the XML result.
    This is a slow network call (often 30 s to several minutes), so callers
    should cache the result.
    """
    from Bio.Blast import NCBIWWW

    kwargs = dict(expect=evalue, hitlist_size=hitlist_size, format_type="XML")
    if entrez_query:
        kwargs["entrez_query"] = entrez_query
    handle = NCBIWWW.qblast(program, database, seq, **kwargs)
    try:
        return handle.read()
    finally:
        handle.close()


_ORGANISM_RE = re.compile(r"\[([^\[\]]+)\]\s*$")


def parse_blast_xml(xml_text: str, query_length: int) -> pd.DataFrame:
    """
    Turn BLAST XML into one row per hit (best HSP per hit).
    Identity and coverage are percentages; coverage is relative to the query.
    """
    from Bio.Blast import NCBIXML

    rows = []
    for record in NCBIXML.parse(StringIO(xml_text)):
        for aln in record.alignments:
            if not aln.hsps:
                continue
            hsp = min(aln.hsps, key=lambda h: h.expect)
            title = getattr(aln, "hit_def", None) or aln.title or ""
            acc = getattr(aln, "accession", None) or aln.hit_id or ""
            org_match = _ORGANISM_RE.search(title)
            organism = org_match.group(1) if org_match else "Unknown"
            description = _ORGANISM_RE.sub("", title).strip()
            covered = hsp.query_end - hsp.query_start + 1
            rows.append({
                "accession": acc,
                "description": description,
                "organism": organism,
                "hit_length": aln.length,
                "identity_pct": round(100.0 * hsp.identities / hsp.align_length, 2),
                "coverage_pct": round(100.0 * covered / query_length, 2) if query_length else 0.0,
                "query_start": hsp.query_start,
                "query_end": hsp.query_end,
                "subject_start": hsp.sbjct_start,
                "subject_end": hsp.sbjct_end,
                "align_length": hsp.align_length,
                "gaps": hsp.gaps,
                "evalue": hsp.expect,
                "bit_score": hsp.score,
            })

    columns = ["accession", "description", "organism", "hit_length", "identity_pct",
               "coverage_pct", "query_start", "query_end", "subject_start", "subject_end",
               "align_length", "gaps", "evalue", "bit_score"]
    df = pd.DataFrame(rows, columns=columns)
    if not df.empty:
        df = df.sort_values(["evalue", "bit_score"], ascending=[True, False]).reset_index(drop=True)
        df.insert(0, "rank", range(1, len(df) + 1))
    return df


def filter_hits(df: pd.DataFrame, min_identity: float = 0.0, min_coverage: float = 0.0,
                max_evalue: float = 10.0) -> pd.DataFrame:
    if df.empty:
        return df
    mask = (
        (df["identity_pct"] >= min_identity)
        & (df["coverage_pct"] >= min_coverage)
        & (df["evalue"] <= max_evalue)
    )
    return df[mask].reset_index(drop=True)


def ortholog_candidates(df: pd.DataFrame) -> pd.DataFrame:
    """
    Best hit per organism. This gives candidate orthologs, not confirmed ones:
    confirming orthology needs a reciprocal best hit check or a phylogenetic
    analysis, which this tool does not do.
    """
    if df.empty:
        return df
    best = (df.sort_values(["evalue", "bit_score"], ascending=[True, False])
              .groupby("organism", sort=False)
              .head(1)
              .reset_index(drop=True))
    return best


# ---------------------------------------------------------------------------
# Export helpers
# ---------------------------------------------------------------------------

def to_fasta(records: list[tuple[str, str]], width: int = 60) -> str:
    """records: list of (header, sequence). Returns FASTA text with wrapped lines."""
    out = []
    for header, seq in records:
        out.append(f">{header}")
        for i in range(0, len(seq), width):
            out.append(seq[i:i + width])
    return "\n".join(out) + "\n"


def build_summary_report(name: str, seq: str, seq_type: str, orfs: pd.DataFrame | None,
                         blast_hits: pd.DataFrame | None = None) -> str:
    lines = [
        "SEQUENCE ANALYSIS REPORT",
        "=" * 40,
        f"Name:            {name or 'unnamed'}",
        f"Type:            {seq_type}",
        f"Length:          {len(seq)} bp" if seq_type != "protein" else f"Length:          {len(seq)} aa",
    ]
    if seq_type != "protein":
        lines += [
            f"GC content:      {gc_content(seq):.2f} %",
            f"AT content:      {at_content(seq):.2f} %",
        ]
        if orfs is not None:
            lines.append(f"ORFs found:      {len(orfs)} (min 30 aa)")
    if blast_hits is not None and not blast_hits.empty:
        lines.append(f"BLAST candidate hits: {len(blast_hits)}")
        top = blast_hits.iloc[0]
        lines.append(f"Top hit:         {top['accession']} | {top['organism']} | "
                     f"{top['identity_pct']}% identity | E={top['evalue']}")
    return "\n".join(lines) + "\n"
