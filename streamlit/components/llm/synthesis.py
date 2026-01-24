"""LLM Research Synthesis component for OncoMind."""

from typing import Any
import streamlit as st

from oncomind.config.constants import UI_MAX_RESEARCH_HYPOTHESES, UI_MAX_REFERENCES


def render_llm_synthesis(result: dict[str, Any]) -> None:
    """Render the LLM Research Synthesis section.

    Args:
        result: The full insight result dictionary containing 'insight' key
    """
    # Check for LLM error first
    llm_rationale = result['insight'].get('rationale', '')
    if llm_rationale and llm_rationale.startswith("LLM narrative generation failed:"):
        error_msg = llm_rationale.replace("LLM narrative generation failed: ", "")
        st.warning(f"⚠️ **LLM synthesis failed:** {error_msg}\n\nEvidence-only results are shown above.")
        return

    llm_narrative = result['insight'].get('llm_narrative')
    if not llm_narrative:
        return

    if llm_rationale and llm_rationale.startswith("LLM narrative generation failed:"):
        return

    with st.container(border=True):
        st.markdown("""### Experimental - LLM Research Synthesis <span style='font-size: 0.7rem; font-weight: normal; color: #666;'>AI-generated — may contain inaccuracies. Verify against primary sources.</span>""", unsafe_allow_html=True)
        st.markdown("""<div style='font-size: 0.95rem; color: #1976d2;'>
<b>Locus Match:</b> variant-level = evidence match at exact locus &nbsp;·&nbsp; codon-level = evidence match at other variants in this codon &nbsp;·&nbsp; gene-level = evidence match at other variants in this gene
</div>""", unsafe_allow_html=True)
        st.markdown("<hr style='margin: 0.5rem 0; border: none; border-top: 1px solid #ddd;'>", unsafe_allow_html=True)

        functional_summary = result['insight'].get('functional_summary')
        biological_context = result['insight'].get('biological_context')
        therapeutic_summary = result['insight'].get('therapeutic_summary')

        if functional_summary:
            st.markdown(f"**Functional Impact:** {functional_summary}")
        if biological_context:
            st.markdown(f"**Biological Context:** {biological_context}")
        if therapeutic_summary:
            st.markdown(f"**Therapeutic Landscape:** {therapeutic_summary}")

        # NOTE: Therapeutic Landscape (structured dict) removed from LLM output
        # This data is already shown in the Therapies tab with accurate source attribution
        # LLM was adding context from training knowledge (hallucination risk)
        # But we now show therapeutic_summary (prose) which is constrained to provided evidence

        if not any([functional_summary, biological_context, therapeutic_summary]):
            st.markdown(llm_narrative)

        conflicting_evidence = result['insight'].get('conflicting_evidence', [])
        if conflicting_evidence:
            st.markdown(f"**⚠️ Conflicting Evidence:** {' · '.join(conflicting_evidence)}")

        evidence_tags = [tag.title() for tag in result['insight'].get('evidence_tags', [])]
        if evidence_tags:
            st.markdown(f"**🏷️ Evidence Types:** {' · '.join(evidence_tags)}")

        research_implications = result['insight'].get('research_implications')
        if research_implications and research_implications != result['insight'].get('rationale'):
            st.markdown(f"**🔬 Research Implications:** {research_implications}")

        research_hypotheses = result['insight'].get('research_hypotheses', [])
        if research_hypotheses:
            st.markdown("**💡 Emerging Research Hypotheses:**")
            for i, hypothesis in enumerate(research_hypotheses[:UI_MAX_RESEARCH_HYPOTHESES], 1):
                st.markdown(f"  {i}. {hypothesis}")

        references = result['insight'].get('references', [])
        if references:
            clickable_refs = _format_references(references[:UI_MAX_REFERENCES])
            st.markdown(f"**📚 Key References:** {', '.join(clickable_refs)}")

        st.caption("_Synthesis incorporates established domain knowledge beyond queried databases._")


def _format_references(references: list) -> list[str]:
    """Format references as clickable links.

    Args:
        references: List of reference strings (PMIDs, NCT IDs, etc.)

    Returns:
        List of formatted markdown links
    """
    clickable_refs = []
    for ref in references:
        ref_str = str(ref).strip()
        if ref_str.startswith("PMID"):
            pmid = ref_str.replace("PMID", "").replace(":", "").strip()
            clickable_refs.append(f"[PMID {pmid}](https://pubmed.ncbi.nlm.nih.gov/{pmid}/)")
        elif ref_str.isdigit():
            clickable_refs.append(f"[PMID {ref_str}](https://pubmed.ncbi.nlm.nih.gov/{ref_str}/)")
        elif "cBioPortal" in ref_str or "cbioportal" in ref_str.lower():
            if ":" in ref_str:
                study_id = ref_str.split(":")[-1].strip()
                clickable_refs.append(f"[{ref_str}](https://www.cbioportal.org/study/summary?id={study_id})")
            else:
                clickable_refs.append(ref_str)
        elif ref_str.startswith("NCT"):
            clickable_refs.append(f"[{ref_str}](https://clinicaltrials.gov/study/{ref_str})")
        else:
            clickable_refs.append(ref_str)
    return clickable_refs
