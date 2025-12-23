# Agentic Architecture for OncoMind

This document outlines how agents could be integrated into OncoMind for complex, multi-step reasoning tasks.

---

## When to Use Agents vs. Deterministic Pipelines

| Task Type | Best Approach | Why |
|-----------|---------------|-----|
| Variant annotation | Deterministic pipeline | Fixed steps, no branching logic |
| Literature search | Deterministic + LLM scoring | Search is API call, scoring is single LLM call |
| Evidence synthesis | Single LLM call | One-shot generation with context |
| **Trial eligibility** | **Agent** | Multi-step reasoning, conditional logic |
| **Clinical Q&A** | **Agent** | Retrieve → reason → maybe retrieve more |
| **Multi-variant reasoning** | **Agent** | Combinatorial, needs iterative search |

**Rule of thumb:** Use agents when the task requires conditional branching or iterative refinement based on intermediate results.

---

## Agent 1: Trial Eligibility Agent

### The Problem

Clinical trial eligibility criteria are free text:

```
Inclusion Criteria:
- Histologically confirmed NSCLC
- Documented EGFR exon 19 deletion or L858R mutation
- Disease progression on or after osimertinib
- ECOG performance status 0-1

Exclusion Criteria:
- Known EGFR T790M at time of enrollment
- Prior treatment with MET inhibitor
- Symptomatic brain metastases
```

Determining eligibility requires:
1. Fetching the criteria
2. Parsing into structured format
3. Matching against patient data
4. Handling ambiguity (ask for clarification or flag uncertainty)

### Agent Design

```python
from typing import Literal
from pydantic import BaseModel

class PatientProfile(BaseModel):
    """Input: what we know about the patient."""
    variants: list[str]                    # ["EGFR L858R", "MET amp"]
    tumor_type: str                        # "NSCLC"
    prior_therapies: list[str] | None      # ["osimertinib", "carboplatin"]
    ecog_status: int | None                # 1
    brain_mets: bool | None                # False
    other_factors: dict[str, str] | None   # {"age": "65", "histology": "adenocarcinoma"}

class ParsedEligibility(BaseModel):
    """Structured representation of trial criteria."""
    nct_id: str

    # Variant requirements
    required_variants: list[str]           # ["EGFR exon 19 del", "EGFR L858R"]
    variant_logic: Literal["any_of", "all_of"]
    excluded_variants: list[str]           # ["EGFR T790M"]

    # Therapy requirements
    required_prior_therapies: list[str]    # ["EGFR TKI"]
    excluded_prior_therapies: list[str]    # ["MET inhibitor"]

    # Other requirements
    other_inclusion: list[str]             # ["ECOG 0-1", "Measurable disease"]
    other_exclusion: list[str]             # ["Symptomatic brain mets"]

class EligibilityVerdict(BaseModel):
    """Output: is the patient eligible?"""
    nct_id: str
    verdict: Literal["eligible", "ineligible", "uncertain", "needs_info"]
    confidence: float

    # Detailed reasoning
    criteria_met: list[str]                # ["EGFR L858R matches inclusion"]
    criteria_not_met: list[str]            # ["Prior MET inhibitor excluded"]
    criteria_uncertain: list[str]          # ["Brain met status unknown"]

    # What's missing
    missing_info: list[str]                # ["Need brain met status"]

    # Full reasoning chain
    reasoning_steps: list[str]

class TrialEligibilityAgent:
    """Agent for determining clinical trial eligibility."""

    def __init__(self, llm_model: str = "gpt-4o"):
        self.llm_model = llm_model
        self.tools = [
            self.fetch_trial_criteria,
            self.parse_eligibility,
            self.check_variant_match,
            self.check_therapy_history,
            self.request_clarification,
        ]

    async def run(
        self,
        patient: PatientProfile,
        nct_id: str
    ) -> EligibilityVerdict:
        """
        Main agent loop.

        1. Fetch trial criteria
        2. Parse into structured format
        3. Check each criterion against patient
        4. Handle uncertainty (ask or flag)
        5. Return verdict with reasoning
        """
        ...

    # --- Tools ---

    async def fetch_trial_criteria(self, nct_id: str) -> str:
        """Fetch eligibility criteria text from ClinicalTrials.gov."""
        # Uses existing ClinicalTrialsClient
        ...

    async def parse_eligibility(self, criteria_text: str) -> ParsedEligibility:
        """LLM parses free text into structured criteria."""
        prompt = """
        Parse this clinical trial eligibility into structured format.

        {criteria_text}

        Return JSON with:
        - required_variants: variants that qualify patient
        - variant_logic: "any_of" or "all_of"
        - excluded_variants: variants that disqualify patient
        - required_prior_therapies: treatments patient must have had
        - excluded_prior_therapies: treatments that disqualify
        - other_inclusion: other inclusion criteria
        - other_exclusion: other exclusion criteria
        """
        ...

    async def check_variant_match(
        self,
        patient_variants: list[str],
        parsed: ParsedEligibility
    ) -> dict:
        """Check if patient variants match trial requirements."""
        # This might need LLM for fuzzy matching:
        # "EGFR L858R" should match "EGFR exon 21 L858R mutation"
        ...

    async def check_therapy_history(
        self,
        prior_therapies: list[str],
        parsed: ParsedEligibility
    ) -> dict:
        """Check if therapy history meets requirements."""
        ...

    async def request_clarification(self, missing: list[str]) -> None:
        """Flag what additional info is needed."""
        ...
```

### Agent Loop (Pseudocode)

```python
async def run(self, patient: PatientProfile, nct_id: str) -> EligibilityVerdict:
    reasoning_steps = []

    # Step 1: Fetch criteria
    criteria_text = await self.fetch_trial_criteria(nct_id)
    reasoning_steps.append(f"Fetched eligibility criteria for {nct_id}")

    # Step 2: Parse into structure
    parsed = await self.parse_eligibility(criteria_text)
    reasoning_steps.append(f"Parsed criteria: {len(parsed.required_variants)} required variants, {len(parsed.excluded_variants)} excluded")

    # Step 3: Check variants
    variant_result = await self.check_variant_match(patient.variants, parsed)
    reasoning_steps.append(f"Variant check: {variant_result['status']}")

    if variant_result['status'] == 'no_match':
        return EligibilityVerdict(
            nct_id=nct_id,
            verdict="ineligible",
            confidence=0.9,
            criteria_not_met=[variant_result['reason']],
            reasoning_steps=reasoning_steps,
            ...
        )

    # Step 4: Check therapy history
    therapy_result = await self.check_therapy_history(patient.prior_therapies, parsed)
    reasoning_steps.append(f"Therapy check: {therapy_result['status']}")

    if therapy_result['status'] == 'excluded':
        return EligibilityVerdict(
            nct_id=nct_id,
            verdict="ineligible",
            confidence=0.85,
            criteria_not_met=[therapy_result['reason']],
            reasoning_steps=reasoning_steps,
            ...
        )

    # Step 5: Check other criteria
    # ... similar pattern

    # Step 6: Handle missing info
    missing = []
    if patient.brain_mets is None and "brain" in str(parsed.other_exclusion).lower():
        missing.append("Brain metastasis status")

    if missing:
        return EligibilityVerdict(
            nct_id=nct_id,
            verdict="needs_info",
            confidence=0.5,
            missing_info=missing,
            reasoning_steps=reasoning_steps,
            ...
        )

    # Step 7: Return verdict
    return EligibilityVerdict(
        nct_id=nct_id,
        verdict="eligible",
        confidence=0.8,
        criteria_met=[...],
        reasoning_steps=reasoning_steps,
        ...
    )
```

### Example Interaction

```python
patient = PatientProfile(
    variants=["EGFR L858R", "MET amp"],
    tumor_type="NSCLC",
    prior_therapies=["erlotinib", "osimertinib"],
    ecog_status=1,
    brain_mets=False,
)

agent = TrialEligibilityAgent()
verdict = await agent.run(patient, "NCT04816214")

print(verdict)
# EligibilityVerdict(
#     nct_id="NCT04816214",
#     verdict="eligible",
#     confidence=0.85,
#     criteria_met=[
#         "EGFR L858R matches 'EGFR mutation' requirement",
#         "Prior osimertinib satisfies 'progression on EGFR TKI'",
#         "ECOG 1 within 0-1 range",
#         "No brain mets (exclusion avoided)",
#     ],
#     criteria_not_met=[],
#     criteria_uncertain=[
#         "MET amp not explicitly mentioned - may be relevant for stratification"
#     ],
#     reasoning_steps=[...]
# )
```

---

## Agent 2: Clinical Q&A Agent

### The Problem

User asks: "Should I try pembrolizumab for BRAF V600E melanoma?"

This requires:
1. Understanding the question (drug, variant, tumor)
2. Fetching relevant evidence
3. Reasoning about standard of care
4. Citing sources
5. Surfacing caveats

### Agent Design

```python
class ClinicalQuestion(BaseModel):
    """Parsed clinical question."""
    question_type: Literal["drug_choice", "prognosis", "resistance", "trial", "other"]
    drug: str | None
    variant: str | None
    tumor_type: str | None
    context: str | None  # "after BRAF inhibitor failure"

class ClinicalAnswer(BaseModel):
    """Grounded answer to clinical question."""
    question: str
    answer: str
    confidence: Literal["high", "moderate", "low"]

    supporting_evidence: list[AttributedClaim]
    contradicting_evidence: list[AttributedClaim]

    caveats: list[str]
    evidence_gaps: list[str]

    reasoning_steps: list[str]

class ClinicalQAAgent:
    """Agent for answering clinical questions about variants."""

    def __init__(self, llm_model: str = "gpt-4o"):
        self.llm_model = llm_model
        self.tools = [
            self.parse_question,
            self.fetch_variant_evidence,
            self.fetch_drug_evidence,
            self.search_literature,
            self.synthesize_answer,
        ]

    async def run(self, question: str) -> ClinicalAnswer:
        """
        1. Parse the question
        2. Fetch relevant evidence (may iterate)
        3. Synthesize answer with citations
        4. Surface caveats and gaps
        """
        ...

    async def parse_question(self, question: str) -> ClinicalQuestion:
        """Extract drug, variant, tumor, question type."""
        ...

    async def fetch_variant_evidence(self, variant: str, tumor: str) -> EvidencePanel:
        """Get OncoMind evidence panel for the variant."""
        # Uses existing get_insight()
        ...

    async def fetch_drug_evidence(self, drug: str, variant: str) -> list[dict]:
        """Get drug-specific evidence from FDA, trials."""
        ...

    async def search_literature(self, query: str) -> list[PubMedEvidence]:
        """Search for additional literature if needed."""
        ...

    async def synthesize_answer(
        self,
        question: ClinicalQuestion,
        evidence: dict
    ) -> ClinicalAnswer:
        """Generate grounded answer from evidence."""
        ...
```

### Agent Loop

```python
async def run(self, question: str) -> ClinicalAnswer:
    steps = []

    # Step 1: Parse question
    parsed = await self.parse_question(question)
    steps.append(f"Parsed: drug={parsed.drug}, variant={parsed.variant}, tumor={parsed.tumor_type}")

    # Step 2: Fetch variant evidence
    if parsed.variant:
        panel = await self.fetch_variant_evidence(parsed.variant, parsed.tumor_type)
        steps.append(f"Retrieved evidence from {len(panel.meta.sources_with_data)} sources")

    # Step 3: Check if we have enough evidence
    if not panel.clinical.fda_approvals and not panel.kb.civic_assertions:
        # Need to search literature
        lit = await self.search_literature(f"{parsed.variant} {parsed.drug} {parsed.tumor_type}")
        steps.append(f"Searched literature, found {len(lit)} papers")

    # Step 4: Synthesize answer
    answer = await self.synthesize_answer(parsed, {
        "panel": panel,
        "literature": lit,
    })
    answer.reasoning_steps = steps

    return answer
```

### Example Interaction

```python
agent = ClinicalQAAgent()
answer = await agent.run("Should I try pembrolizumab for BRAF V600E melanoma?")

print(answer)
# ClinicalAnswer(
#     question="Should I try pembrolizumab for BRAF V600E melanoma?",
#     answer="Pembrolizumab is not first-line for BRAF V600E melanoma. BRAF/MEK inhibitor
#            combination (dabrafenib + trametinib or vemurafenib + cobimetinib) is preferred
#            per NCCN guidelines. However, pembrolizumab may be considered after BRAF inhibitor
#            progression or if patient has contraindication to targeted therapy.",
#     confidence="high",
#     supporting_evidence=[
#         AttributedClaim(claim="BRAF/MEK combination is first-line", source="NCCN Guidelines"),
#         AttributedClaim(claim="Dabrafenib+trametinib FDA approved", source="FDA label"),
#     ],
#     caveats=[
#         "Immunotherapy may be preferred if high tumor mutational burden",
#         "Consider clinical trial if available",
#     ],
#     evidence_gaps=[
#         "Optimal sequencing of targeted vs immunotherapy not fully established",
#     ],
# )
```

---

## Agent 3: Multi-Variant Reasoning Agent

### The Problem

Patient has: EGFR L858R + TP53 R248W + MET amplification

Each variant alone has evidence. Together, they interact:
- TP53 mutation may reduce EGFR TKI efficacy
- MET amp may cause primary resistance to EGFR TKI
- Combined prognosis is worse than EGFR alone

### Agent Design

```python
class VariantInteraction(BaseModel):
    """Interaction between two variants."""
    variant_a: str
    variant_b: str
    interaction_type: Literal["synergistic", "antagonistic", "prognostic", "none", "unknown"]
    clinical_impact: str
    evidence_level: str
    source_pmids: list[str]

class MultiVariantAnalysis(BaseModel):
    """Combined analysis of multiple variants."""
    variants: list[str]
    tumor_type: str

    # Individual annotations
    individual_panels: dict[str, EvidencePanel]

    # Interactions
    interactions: list[VariantInteraction]

    # Combined implications
    combined_prognosis: str
    therapeutic_implications: list[str]
    resistance_risk: str

    # Literature specific to combination
    combination_papers: list[PubMedEvidence]

    reasoning_steps: list[str]

class MultiVariantAgent:
    """Agent for reasoning about multiple variants together."""

    async def run(
        self,
        variants: list[str],
        tumor_type: str
    ) -> MultiVariantAnalysis:
        """
        1. Annotate each variant individually
        2. Search for co-occurrence literature
        3. Identify interactions
        4. Synthesize combined implications
        """
        ...
```

---

## Implementation Approach

### Phase 1: Learn the Primitives
Start with raw Claude tool use to understand the mechanics:

```python
import anthropic

client = anthropic.Anthropic()

tools = [
    {
        "name": "fetch_trial_criteria",
        "description": "Fetch eligibility criteria for a clinical trial",
        "input_schema": {
            "type": "object",
            "properties": {
                "nct_id": {"type": "string", "description": "ClinicalTrials.gov ID"}
            },
            "required": ["nct_id"]
        }
    },
    # ... more tools
]

response = client.messages.create(
    model="claude-sonnet-4-20250514",
    max_tokens=4096,
    tools=tools,
    messages=[{"role": "user", "content": "Is a patient with EGFR L858R eligible for NCT04816214?"}]
)

# Handle tool calls in a loop
while response.stop_reason == "tool_use":
    tool_call = response.content[0]  # Get tool call
    result = execute_tool(tool_call)  # Run the tool
    response = client.messages.create(
        model="claude-sonnet-4-20250514",
        max_tokens=4096,
        tools=tools,
        messages=[
            {"role": "user", "content": "..."},
            {"role": "assistant", "content": response.content},
            {"role": "user", "content": [{"type": "tool_result", "tool_use_id": tool_call.id, "content": result}]}
        ]
    )
```

### Phase 2: Add Structure
Wrap in a class with proper state management:

```python
class AgentRunner:
    def __init__(self, tools: list[Tool]):
        self.tools = tools
        self.client = anthropic.Anthropic()

    async def run(self, task: str, max_steps: int = 10) -> AgentResult:
        messages = [{"role": "user", "content": task}]
        steps = []

        for _ in range(max_steps):
            response = await self.client.messages.create(...)

            if response.stop_reason == "end_turn":
                return AgentResult(answer=response.content, steps=steps)

            if response.stop_reason == "tool_use":
                tool_call = self.extract_tool_call(response)
                result = await self.execute_tool(tool_call)
                steps.append(Step(tool=tool_call.name, result=result))
                messages.append(...)

        return AgentResult(answer=None, steps=steps, error="Max steps reached")
```

### Phase 3: Consider a Framework
If the raw approach gets unwieldy, consider:

- **LangGraph** — Good for complex state machines, explicit control flow
- **Claude Agent SDK** — If Anthropic releases one (check current status)
- **PydanticAI** — Type-safe agent framework (https://github.com/pydantic/pydantic-ai)

---

## Where Agents Fit in OncoMind Roadmap

| Version | Feature | Agent? |
|---------|---------|--------|
| v0.2 | Trust Foundation | No |
| v0.3 | Structural Variants | No |
| v0.4 | Validation | No |
| v0.5 | Visualizations | No |
| v0.6 | Ensemble LLM | No (parallel calls, not agent) |
| v0.7 | Clinical Intelligence | **Maybe** (trial eligibility) |
| v0.8 | LLM Infrastructure | No |
| Future | Clinical Q&A | **Yes** |
| Future | Multi-Variant | **Yes** |
| Future | Contextual Q&A | **Yes** |

**Recommendation:** Build v0.2-v0.6 without agents. Introduce agents in v0.7 with Trial Eligibility as the learning project.

---

## Resources

- [Anthropic Tool Use Guide](https://docs.anthropic.com/claude/docs/tool-use)
- [Building Effective Agents](https://www.anthropic.com/research/building-effective-agents) (Anthropic blog)
- [LangGraph Documentation](https://langchain-ai.github.io/langgraph/)
- [PydanticAI](https://github.com/pydantic/pydantic-ai)

---

## Open Questions

1. **State persistence** — Do agents need memory across sessions? (Probably not for OncoMind use cases)
2. **Error recovery** — What happens when a tool fails mid-agent-loop?
3. **Cost control** — Agents can make many LLM calls. How to cap?
4. **Evaluation** — How to benchmark agent performance vs. single-call approaches?
