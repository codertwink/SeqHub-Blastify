"""
Sequence Analysis Dashboard
Run locally:  streamlit run app.py
"""

import math

import plotly.graph_objects as go
import streamlit as st

import seqtools as sq

st.set_page_config(page_title="SeqHub Blastify", page_icon="🧬", layout="wide")

EXAMPLE_SEQUENCE = """Paste your sequence here"""

FRAME_ORDER = ["+1", "+2", "+3", "-1", "-2", "-3"]


# ---------------------------------------------------------------------------
# Cached helpers
# ---------------------------------------------------------------------------

@st.cache_data(ttl=60 * 60 * 24, show_spinner=False)
def cached_blast(seq: str, program: str, database: str, evalue: float, hitlist: int) -> str:
    return sq.run_blast(seq, program=program, database=database,
                        evalue=evalue, hitlist_size=hitlist)


def fig_layout(fig: go.Figure, height: int = 380) -> go.Figure:
    fig.update_layout(height=height, margin=dict(l=10, r=10, t=40, b=10),
                      template="plotly_white")
    return fig


# ---------------------------------------------------------------------------
# Sidebar: input
# ---------------------------------------------------------------------------

st.sidebar.title("🧬 Sequence Input")
input_mode = st.sidebar.radio("Input method", ["Paste sequence", "Upload FASTA file"], horizontal=False)

raw_text = ""
if input_mode == "Paste sequence":
    raw_text = st.sidebar.text_area(
        "Sequence (raw or FASTA)",
        value=EXAMPLE_SEQUENCE,
        height=260,
        help="Paste DNA, RNA, or protein. Spaces, numbers and line breaks are ignored.",
    )
else:
    uploaded = st.sidebar.file_uploader("FASTA / TXT file", type=["fasta", "fa", "fna", "txt"])
    if uploaded is not None:
        raw_text = uploaded.getvalue().decode("utf-8", errors="ignore")

seq_name_input = st.sidebar.text_input("Sequence name (optional)", value="")

try:
    header, seq = sq.parse_fasta_or_raw(raw_text)
    seq_type = sq.detect_sequence_type(seq)
except sq.SequenceInputError as err:
    st.sidebar.error(str(err))
    st.title("Sequence Analysis Dashboard")
    st.info("Paste a sequence or upload a FASTA file in the sidebar to begin.")
    st.stop()

seq_name = seq_name_input or header or "unnamed_sequence"
st.sidebar.success(f"Detected: **{seq_type}**, length **{len(seq):,}**")
if len(seq) > 30000:
    st.sidebar.warning("Long sequence: BLAST will be slow and NCBI may reject very large queries.")

st.title(f"Sequence Analysis: {seq_name}")

m1, m2, m3, m4 = st.columns(4)
m1.metric("Type", seq_type)
m2.metric("Length", f"{len(seq):,} {'aa' if seq_type == 'protein' else 'bp'}")
if seq_type != "protein":
    m3.metric("GC content", f"{sq.gc_content(seq):.2f} %")
    m4.metric("AT content", f"{sq.at_content(seq):.2f} %")
else:
    m3.metric("Stop codons (*)", seq.count("*"))
    m4.metric("Unique residues", len(set(seq) - {"*"}))

tab_nt, tab_blast, tab_prot, tab_export = st.tabs(
    ["Nucleotide analysis", "BLAST search", "Protein analysis", "Exports"]
)

# ---------------------------------------------------------------------------
# Tab 1: nucleotide analysis
# ---------------------------------------------------------------------------

with tab_nt:
    if seq_type == "protein":
        st.info("This is a protein sequence, so nucleotide analysis does not apply. "
                "See the Protein analysis tab.")
    else:
        counts = sq.nucleotide_counts(seq)
        c1, c2 = st.columns([1, 1])

        with c1:
            st.subheader("Base composition")
            fig = go.Figure(go.Bar(
                x=list(counts.keys()), y=list(counts.values()),
                marker_color=["#2E86AB", "#A23B72", "#F18F01", "#C73E1D", "#8D99AE"][:len(counts)],
                text=list(counts.values()), textposition="outside",
            ))
            fig.update_layout(yaxis_title="Count")
            st.plotly_chart(fig_layout(fig, 340))

        with c2:
            st.subheader("GC content along the sequence")
            w1, w2 = st.columns(2)
            window = w1.slider("Window (bp)", 20, 500, 100, step=10)
            step = w2.slider("Step (bp)", 1, 100, 10)
            gc_df = sq.sliding_gc(seq, window=window, step=step)
            if gc_df.empty:
                st.info("Sequence is shorter than the window.")
            else:
                fig = go.Figure(go.Scatter(x=gc_df["position"], y=gc_df["gc"],
                                           mode="lines", line=dict(color="#2E86AB")))
                fig.add_hline(y=sq.gc_content(seq), line_dash="dash", line_color="grey",
                              annotation_text="overall GC")
                fig.update_layout(xaxis_title="Position (bp)", yaxis_title="GC %")
                st.plotly_chart(fig_layout(fig, 340))

        st.subheader("Sequence transformations")
        st.caption("GC% and AT% are calculated over unambiguous bases (A, C, G, T/U). "
                   "Ambiguous codes such as N are excluded.")
        for label, value in sq.transformations(seq).items():
            with st.expander(f"{label}  ({len(value):,} nt)", expanded=(label == "Reverse complement")):
                st.code(value, language=None, wrap_lines=True)


# ---------------------------------------------------------------------------
# Tab 2: BLAST
# ---------------------------------------------------------------------------

with tab_blast:
    blast_cfg = sq.BLAST_PROGRAMS["protein" if seq_type == "protein" else "nucleotide"]
    st.write(f"Program: **{blast_cfg['program']}**. Runs against NCBI's public BLAST service, "
             "which can take from under a minute to several minutes.")

    b1, b2, b3, b4 = st.columns(4)
    db_label = b1.selectbox("Database", list(blast_cfg["databases"].keys()))
    database = blast_cfg["databases"][db_label]
    evalue = b2.select_slider("Max E-value (search)", options=[1e-50, 1e-20, 1e-10, 1e-5, 1e-2, 1.0, 10.0], value=10.0)
    hitlist = b3.slider("Max hits from NCBI", 10, 500, 100, step=10)
    run_key = (seq, blast_cfg["program"], database, evalue, hitlist)

    if st.button("Run BLAST", type="primary"):
        with st.spinner("Submitting to NCBI BLAST... this can take a few minutes."):
            try:
                xml = cached_blast(seq, blast_cfg["program"], database, evalue, hitlist)
                st.session_state["blast_xml"] = xml
                st.session_state["blast_key"] = run_key
                st.session_state["blast_error"] = None
            except Exception as exc:  # network, NCBI busy, malformed response
                st.session_state["blast_xml"] = None
                st.session_state["blast_error"] = str(exc)

    if st.session_state.get("blast_error"):
        st.error(f"BLAST failed: {st.session_state['blast_error']}. "
                 "NCBI may be busy or unreachable. Try again in a minute.")

    xml = st.session_state.get("blast_xml")
    if xml:
        if st.session_state.get("blast_key") != run_key:
            st.warning("Settings or sequence have changed since the last BLAST run. "
                       "Results below are from the previous run. Click Run BLAST to update.")

        all_hits = sq.parse_blast_xml(xml, query_length=len(seq))
        st.session_state["blast_all"] = all_hits

        st.subheader("Filter hits")
        f1, f2, f3 = st.columns(3)
        min_id = f1.slider("Min identity %", 0.0, 100.0, 70.0, step=0.5)
        min_cov = f2.slider("Min query coverage %", 0.0, 100.0, 50.0, step=1.0)
        max_ev = f3.number_input("Max E-value", value=1e-5, format="%.1e")

        hits = sq.filter_hits(all_hits, min_identity=min_id, min_coverage=min_cov, max_evalue=max_ev)
        st.session_state["blast_filtered"] = hits

        if all_hits.empty:
            st.info("BLAST returned no hits.")
        else:
            k1, k2, k3 = st.columns(3)
            k1.metric("Total hits", len(all_hits))
            k2.metric("Passing filters", len(hits))
            k3.metric("Distinct organisms (passing)", hits["organism"].nunique() if not hits.empty else 0)

            if not hits.empty:
                fig = go.Figure(go.Scatter(
                    x=hits["coverage_pct"], y=hits["identity_pct"], mode="markers",
                    marker=dict(size=(hits["bit_score"] / hits["bit_score"].max() * 30).clip(6, 30),
                                color=hits["evalue"].apply(lambda e: -math.log10(max(e, 1e-300))),
                                colorscale="Viridis", colorbar=dict(title="-log10(E)"), showscale=True),
                    text=hits["accession"] + " | " + hits["organism"],
                    hovertemplate="%{text}<br>coverage %{x:.1f}%<br>identity %{y:.1f}%<extra></extra>",
                ))
                fig.update_layout(xaxis_title="Query coverage %", yaxis_title="Identity %")
                st.plotly_chart(fig_layout(fig, 380))

            st.subheader("Candidate orthologs (best hit per organism)")
            st.caption("These are candidates. Confirming orthology needs a reciprocal best-hit check "
                       "or phylogenetic analysis.")
            cands = sq.ortholog_candidates(hits)
            st.dataframe(cands, hide_index=True)

            st.subheader("All filtered hits")
            st.dataframe(hits, hide_index=True)
    else:
        st.info("No BLAST results yet. Click Run BLAST.")


# ---------------------------------------------------------------------------
# Tab 3: protein analysis
# ---------------------------------------------------------------------------

with tab_prot:
    if seq_type == "protein":
        protein_seq = seq
        st.info("Input is already a protein sequence, so translation is skipped.")
        orfs = None
        frames = None
    else:
        if seq_type == "RNA":
            st.caption("RNA input is converted to DNA for translation.")
        min_aa = st.slider("Minimum ORF length (amino acids)", 10, 300, 30, step=5)
        frames = sq.six_frame_translation(seq)
        orfs = sq.find_orfs(seq, min_aa=min_aa)
        protein_seq = None

        st.subheader("Six-frame translation")
        st.caption("Frames +1 to +3 read the input strand. Frames -1 to -3 read the reverse complement. "
                   "* marks a stop codon.")
        frame_choice = st.selectbox("Frame", FRAME_ORDER, index=0)
        st.code(frames[frame_choice], language=None, wrap_lines=True)

        st.subheader("ORF map")
        if orfs.empty:
            st.info(f"No ORFs of at least {min_aa} amino acids found. Try lowering the minimum length.")
        else:
            fig = go.Figure()
            for _, row in orfs.iterrows():
                y = FRAME_ORDER.index(row["frame"])
                fig.add_trace(go.Scatter(
                    x=[row["start"], row["end"]], y=[y, y], mode="lines",
                    line=dict(width=10, color="#2E86AB" if row["strand"] == "+" else "#F18F01"),
                    name=row["ORF"], showlegend=False,
                    hovertemplate=f"{row['ORF']} ({row['frame']})<br>{row['start']}-{row['end']} nt"
                                  f"<br>{row['length_aa']} aa<extra></extra>",
                ))
            fig.update_yaxes(tickmode="array", tickvals=list(range(6)), ticktext=FRAME_ORDER,
                             autorange="reversed")
            fig.update_layout(xaxis_title="Position on input strand (nt)")
            st.plotly_chart(fig_layout(fig, 340))

            st.subheader(f"ORFs ({len(orfs)})")
            st.dataframe(orfs.drop(columns=["protein"]), hide_index=True)

            orf_choice = st.selectbox("View protein for ORF", orfs["ORF"].tolist())
            protein_seq = orfs.loc[orfs["ORF"] == orf_choice, "protein"].iloc[0]

    if protein_seq:
        st.subheader("Protein properties")
        props = sq.protein_properties(protein_seq)
        if props:
            cols = st.columns(len(props))
            for col, (k, v) in zip(cols, props.items()):
                col.metric(k, v)

        comp = sq.aa_composition(protein_seq)
        c1, c2 = st.columns([2, 1])
        with c1:
            st.subheader("Amino acid composition")
            fig = go.Figure(go.Bar(x=comp["amino_acid"], y=comp["percent"],
                                   marker_color="#2E86AB",
                                   customdata=comp["count"],
                                   hovertemplate="%{x}: %{y:.2f}% (n=%{customdata})<extra></extra>"))
            fig.update_layout(yaxis_title="% of residues")
            st.plotly_chart(fig_layout(fig, 340))
        with c2:
            st.subheader("Protein sequence")
            st.code(protein_seq, language=None, wrap_lines=True)
            st.caption("Composition table")
            st.dataframe(comp, hide_index=True)


# ---------------------------------------------------------------------------
# Tab 4: exports
# ---------------------------------------------------------------------------

with tab_export:
    st.write("Download results. FASTA files are wrapped at 60 characters per line.")

    st.subheader("Sequences (FASTA)")
    if seq_type != "protein":
        tf = sq.transformations(seq)
        records = [(f"{seq_name}_{label.replace(' ', '_')}", value) for label, value in tf.items()]
        st.download_button("Download nucleotide transformations (FASTA)",
                           sq.to_fasta(records), file_name=f"{seq_name}_transformations.fasta",
                           mime="text/plain")

        orfs_exp = sq.find_orfs(seq, min_aa=30)
        if not orfs_exp.empty:
            prot_records = [
                (f"{seq_name}_{r.ORF}_{r.frame}_{r.start}-{r.end}", r.protein)
                for r in orfs_exp.itertuples()
            ]
            st.download_button("Download ORF proteins (FASTA, min 30 aa)",
                               sq.to_fasta(prot_records), file_name=f"{seq_name}_orfs.fasta",
                               mime="text/plain")

        frames_exp = sq.six_frame_translation(seq)
        st.download_button("Download six-frame translation (FASTA)",
                           sq.to_fasta([(f"{seq_name}_frame{f}", p) for f, p in frames_exp.items()]),
                           file_name=f"{seq_name}_six_frames.fasta", mime="text/plain")
    else:
        st.download_button("Download protein sequence (FASTA)",
                           sq.to_fasta([(seq_name, seq)]),
                           file_name=f"{seq_name}.fasta", mime="text/plain")

    st.subheader("Tables (CSV)")
    hits_df = st.session_state.get("blast_filtered")
    if hits_df is not None and not hits_df.empty:
        st.download_button("Download filtered BLAST hits (CSV)",
                           hits_df.to_csv(index=False), file_name=f"{seq_name}_blast_hits.csv",
                           mime="text/csv")
        st.download_button("Download ortholog candidates (CSV)",
                           sq.ortholog_candidates(hits_df).to_csv(index=False),
                           file_name=f"{seq_name}_ortholog_candidates.csv", mime="text/csv")
    else:
        st.caption("BLAST hit exports appear here after a BLAST run with hits.")

    if seq_type != "protein":
        st.download_button("Download ORF table (CSV)",
                           sq.find_orfs(seq, min_aa=30).drop(columns=["protein"]).to_csv(index=False),
                           file_name=f"{seq_name}_orfs.csv", mime="text/csv")
    prot_for_comp = seq if seq_type == "protein" else None
    if prot_for_comp:
        st.download_button("Download amino acid composition (CSV)",
                           sq.aa_composition(prot_for_comp).to_csv(index=False),
                           file_name=f"{seq_name}_aa_composition.csv", mime="text/csv")

    st.subheader("Summary report (TXT)")
    report_orfs = sq.find_orfs(seq, min_aa=30) if seq_type != "protein" else None
    report = sq.build_summary_report(seq_name, seq, seq_type, report_orfs, hits_df)
    st.code(report, language=None)
    st.download_button("Download summary report (TXT)", report,
                       file_name=f"{seq_name}_report.txt", mime="text/plain")
