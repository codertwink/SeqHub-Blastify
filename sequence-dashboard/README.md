# Sequence Analysis Dashboard

A Streamlit web app for nucleotide and protein sequence analysis, built on your
junior's original scripts (GC content, reverse, complement, reverse complement)
and extended with BLAST, six-frame translation, ORF detection, and exports.

## Features

**Input:** paste a raw sequence or FASTA, or upload a `.fasta` / `.fa` / `.txt` file.
DNA, RNA and protein are detected automatically.

**Nucleotide analysis**
- Length, GC %, AT %
- Base composition chart
- Sliding-window GC % plot (adjustable window and step)
- Reverse, complement, reverse complement, DNA-to-RNA

**BLAST search** (NCBI web service)
- blastn for nucleotide input, blastp for protein input
- Choice of database, E-value and hit limit
- Filters on identity %, query coverage % and E-value
- Candidate orthologs: best hit per organism
- Identity vs coverage plot

**Protein analysis**
- Six-frame translation (+1 to +3 and -1 to -3)
- ORF detection on both strands, with an ORF map
- Amino acid composition chart
- Molecular weight, pI, GRAVY and instability index

**Exports:** FASTA (transformations, ORF proteins, six-frame translations),
CSV (BLAST hits, ortholog candidates, ORF table, composition), TXT summary report.

## Project layout

```
app.py            Streamlit UI (tabs, charts, downloads)
seqtools.py       Core logic, no UI code (testable on its own)
requirements.txt  Python dependencies
```

## Run locally

```bash
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -r requirements.txt
streamlit run app.py
```

Then open http://localhost:8501.

## Hosting

**Streamlit Community Cloud (free):** push this folder to a GitHub repo, then at
share.streamlit.io choose the repo and set the main file to `app.py`.

**Render:** create a new Web Service from the repo with:
- Build command: `pip install -r requirements.txt`
- Start command: `streamlit run app.py --server.port $PORT --server.address 0.0.0.0`

## Important notes

- **BLAST needs internet access** to NCBI. Runs usually take 30 seconds to a few
  minutes, and NCBI may be slow or busy at peak times. Results are cached for 24 hours
  for identical queries.
- **Candidate orthologs are not confirmed orthologs.** The "best hit per organism"
  table is a starting point. Confirming orthology needs a reciprocal best-hit check
  or phylogenetic analysis.
- GC % and AT % are calculated over unambiguous bases only (A, C, G, T/U); N and
  other IUPAC codes are excluded from the denominator.
- ORF coordinates are 1-based and always refer to the input strand, including the
  stop codon. Reverse-strand ORFs are mapped back to these coordinates.
- Only the first FASTA record is analysed. Multi-FASTA batch processing is not
  included yet.

## Tests

The core logic was checked against the original `reverse_complement.txt` output and
against known synthetic ORFs (both strands). Run `python -m pytest` if you add tests.
