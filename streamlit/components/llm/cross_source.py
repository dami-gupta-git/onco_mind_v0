"""LLM Cross-Source Drug Analysis component for OncoMind."""

from typing import Any
import streamlit as st


def render_cross_source_analysis(result: dict[str, Any]) -> None:
    """Render the LLM Cross-Source Drug Analysis section.

    Args:
        result: The full insight result dictionary containing 'cross_source_analysis' key
    """
    cross_source = result.get('cross_source_analysis')
    if not cross_source:
        return

    with st.container(border=True):
        st.markdown("### 🔬 LLM Cross-Source Drug Analysis")
        st.caption("AI synthesis of drug evidence across CGI, CIViC, VICC, and Literature sources")

        # Summary at top
        summary = cross_source.get('summary', '')
        if summary:
            st.markdown(f"**Summary:** {summary}")
            st.markdown("<hr style='margin: 0.5rem 0;'>", unsafe_allow_html=True)

        # Create columns for structured display
        col1, col2 = st.columns(2)

        # Strongest Evidence (corroborated across sources)
        with col1:
            _render_strongest_evidence(cross_source.get('strongest_evidence', []))

        # Conflicting Signals
        with col2:
            _render_conflicting_signals(cross_source.get('conflicting_signals', []))

        # Emerging Targets (single-source, needs validation)
        _render_emerging_targets(cross_source.get('emerging_targets', []))

        # Key Gaps
        _render_key_gaps(cross_source.get('key_gaps', []))


def _render_strongest_evidence(strongest: list) -> None:
    """Render the strongest evidence section."""
    if not strongest:
        return

    lines = ["**💪 Strongest Evidence**"]
    for item in strongest[:5]:  # Limit to top 5
        drug = item.get('drug', 'Unknown')
        signal = item.get('signal', '')
        sources = ", ".join(item.get('sources', []))
        level = item.get('evidence_level', '')
        rationale = item.get('rationale', '')
        signal_emoji = "🟢" if signal == "sensitivity" else "🔴" if signal == "resistance" else "⚪"
        lines.append(f"{signal_emoji} **{drug}** ({level})  ")
        lines.append(f"&nbsp;&nbsp;&nbsp;Sources: {sources}  ")
        if rationale:
            lines.append(f"&nbsp;&nbsp;&nbsp;_{rationale}_  ")
    st.markdown("<br>".join(lines), unsafe_allow_html=True)


def _render_conflicting_signals(conflicts: list) -> None:
    """Render the conflicting signals section."""
    if not conflicts:
        return

    lines = ["**⚠️ Conflicting Signals**"]
    for item in conflicts[:5]:
        drug = item.get('drug', 'Unknown')
        conflict = item.get('conflict', '')
        reason = item.get('likely_reason', '')
        question = item.get('research_question', '')
        lines.append(f"**{drug}**: {conflict}  ")
        if reason:
            lines.append(f"&nbsp;&nbsp;&nbsp;_Likely reason: {reason}_  ")
        if question:
            lines.append(f"&nbsp;&nbsp;&nbsp;❓ {question}  ")
    st.markdown("<br>".join(lines), unsafe_allow_html=True)


def _render_emerging_targets(emerging: list) -> None:
    """Render the emerging targets section."""
    if not emerging:
        return

    st.markdown("<hr style='margin: 0.5rem 0;'>", unsafe_allow_html=True)
    lines = ["**🌱 Emerging Targets** (single-source, needs validation)"]
    for item in emerging[:3]:  # Limit to top 3
        drug = item.get('drug', 'Unknown')
        source = item.get('source', '')
        level = item.get('evidence_level', '')
        rationale = item.get('biological_rationale', '')
        validation = item.get('validation_needed', '')
        lines.append(f"• **{drug}** ({source}, {level})  ")
        if rationale:
            lines.append(f"&nbsp;&nbsp;_{rationale}_  ")
        if validation:
            lines.append(f"&nbsp;&nbsp;📋 Needs: {validation}  ")
    st.markdown("<br>".join(lines), unsafe_allow_html=True)


def _render_key_gaps(gaps: list) -> None:
    """Render the key gaps section."""
    if not gaps:
        return

    st.markdown("<hr style='margin: 0.5rem 0;'>", unsafe_allow_html=True)
    lines = ["**🔍 Key Gaps**"]
    for gap in gaps[:3]:
        lines.append(f"• {gap}")
    st.markdown("<br>".join(lines), unsafe_allow_html=True)
