"""
FDA Label Parser for OncoMind
============================

Parses FDA drug labels to extract structured biomarker-indication associations,
properly handling:
- Negation patterns ("no EGFR mutation", "without ALK aberrations")
- Variant-level specificity (gene, codon, specific variant)
- Tumor type extraction
- Combination therapy detection
- Exclusion criteria vs inclusion criteria

Author: Dami (OncoMind)
"""

import logging
import re
from pathlib import Path
from typing import Optional

from oncomind.models.evidence.fda_biomarker import (
    BiomarkerRequirement,
    SpecificityLevel,
    FDABiomarkerEvidence,
    match_variant_to_indication,
)


class FDALabelParser:
    """
    Parser for FDA drug labels from OpenFDA API.

    Handles the full complexity of FDA label text including:
    - Multiple indications per drug
    - Negation patterns
    - Combination therapies
    - Variant specificity
    """

    # Negation patterns that indicate biomarker should be ABSENT
    NEGATION_PATTERNS = [
        r"with\s+no\s+(?:sensitizing\s+)?",
        r"without\s+(?:sensitizing\s+)?",
        r"(?:who\s+)?do\s+not\s+have\s+",
        r"negative\s+for\s+",
        r"wild[- ]?type\s+",
        r"lacking\s+",
        r"absence\s+of\s+",
        r"no\s+(?:known\s+)?(?:sensitizing\s+)?",
        r"excluding\s+(?:patients\s+with\s+)?",
        r"tumors?\s+(?:that\s+)?(?:do\s+)?not\s+(?:have|harbor|express)\s+",
    ]

    # Patterns that indicate the drug is NOT effective (different from negation)
    NEGATIVE_EFFICACY_PATTERNS = [
        r"has\s+not\s+been\s+(?:established|demonstrated|evaluated)",
        r"efficacy\s+(?:was|has)\s+not\s+(?:been\s+)?(?:established|demonstrated)",
        r"not\s+(?:indicated|approved|recommended)\s+for",
        r"should\s+not\s+be\s+used\s+(?:in|for)",
        r"is\s+not\s+effective\s+(?:in|for)",
        r"safety\s+and\s+efficacy.*not\s+(?:been\s+)?established",
    ]

    # Patterns indicating prior therapy requirements (not an indication FOR this biomarker)
    # e.g., "Patients with EGFR mutations should have disease progression on FDA-approved therapy"
    # This means: if you have EGFR, you must fail targeted therapy first - NOT that this drug treats EGFR
    PRIOR_THERAPY_PATTERNS = [
        r"should\s+have\s+disease\s+progression\s+on",
        r"must\s+have\s+(?:disease\s+)?progression\s+(?:on|following)",
        r"(?:have|has)\s+progressed\s+(?:on|after|following)\s+(?:FDA[- ]approved|targeted)\s+therapy",
        r"after\s+(?:failure|progression)\s+(?:of|on)\s+(?:FDA[- ]approved|targeted)\s+therapy",
        r"who\s+have\s+(?:disease\s+)?progression\s+on\s+FDA[- ]approved\s+therapy\s+for\s+these\s+aberrations",
        # Conditional prior therapy: "and if BRAF V600 positive, a BRAF inhibitor"
        # This means: IF you have this mutation, you must have tried this therapy first
        r"(?:and\s+)?if\s+(?:\w+\s+)?(?:V600|mutation)[- ]?positive,?\s+(?:a\s+)?(?:BRAF|EGFR|ALK)\s+inhibitor",
        r"if\s+(?:BRAF|EGFR|ALK)\s+(?:V600\s+)?(?:mutation[- ]?)?positive",
    ]

    # Patterns indicating mechanism-of-action (drug targets the protein, NOT a mutation requirement)
    # e.g., "EGFR antagonist" means the drug targets EGFR protein, not that patients need EGFR mutations
    # These patterns match text FOLLOWING a gene name
    MECHANISM_OF_ACTION_PATTERNS = [
        r"antagonist",
        r"inhibitor",
        r"antibody",
        r"blocker",
        r"targeting\s+(?:agent|drug|therapy)",
        r"directed\s+(?:therapy|treatment|agent)",
    ]

    # Biomarker/gene patterns
    GENE_PATTERNS = {
        "EGFR": r"(?:epidermal\s+growth\s+factor\s+receptor|EGFR)",
        "ALK": r"(?:anaplastic\s+lymphoma\s+kinase|ALK)",
        "BRAF": r"BRAF",
        "KRAS": r"KRAS",
        "NRAS": r"NRAS",
        "HER2": r"(?:HER2|ERBB2|human\s+epidermal\s+growth\s+factor\s+receptor\s+2)",
        "ROS1": r"ROS1",
        "RET": r"RET",
        "MET": r"\bMET\b(?!\s*(?:astatic|hod|ric|al))",  # Exclude "metastatic", "method", etc.
        "NTRK": r"(?:NTRK[123]?|neurotrophic\s+tyrosine\s+receptor\s+kinase)",
        "PIK3CA": r"PIK3CA",
        "AKT1": r"AKT1",
        "IDH1": r"(?:isocitrate\s+dehydrogenase[- ]?1|IDH1)",
        "IDH2": r"(?:isocitrate\s+dehydrogenase[- ]?2|IDH2)",
        "FGFR": r"(?:fibroblast\s+growth\s+factor\s+receptor|FGFR[1234]?)",
        "BRCA": r"BRCA[12]?",
        "KIT": r"(?:c[- ]?Kit|Kit\s*\(CD117\)|KIT|CD117)",
        "PDGFRA": r"(?:PDGFRA|PDGFR\s*[αa]?|platelet[- ]derived\s+growth\s+factor\s+receptor(?:\s+alpha)?)",
        "MSI": r"(?:microsatellite\s+instability|MSI[- ]?(?:H|high)?)",
        "TMB": r"(?:tumor\s+mutational\s+burden|TMB[- ]?(?:H|high)?)",
        "PD-L1": r"(?:PD[- ]?L1|CD274)",
    }

    # Variant patterns (gene-specific)
    VARIANT_PATTERNS = {
        "EGFR": [
            (r"exon\s+19\s+deletion", "exon 19 del", SpecificityLevel.CODON),
            (r"exon\s+21\s+(?:\()?L858R(?:\))?", "L858R", SpecificityLevel.VARIANT),
            (r"L858R", "L858R", SpecificityLevel.VARIANT),
            (r"T790M", "T790M", SpecificityLevel.VARIANT),
            (r"exon\s+20\s+insertion", "exon 20 ins", SpecificityLevel.CODON),
            (r"C797S", "C797S", SpecificityLevel.VARIANT),
        ],
        "BRAF": [
            (r"V600E", "V600E", SpecificityLevel.VARIANT),
            (r"V600K", "V600K", SpecificityLevel.VARIANT),
            (r"V600", "V600", SpecificityLevel.CODON),
        ],
        "KRAS": [
            (r"G12C", "G12C", SpecificityLevel.VARIANT),
            (r"G12D", "G12D", SpecificityLevel.VARIANT),
            (r"G12", "G12", SpecificityLevel.CODON),
        ],
        "IDH1": [
            (r"R132H", "R132H", SpecificityLevel.VARIANT),
            (r"R132C", "R132C", SpecificityLevel.VARIANT),
            (r"R132", "R132", SpecificityLevel.CODON),
        ],
        "IDH2": [
            (r"R172K", "R172K", SpecificityLevel.VARIANT),
            (r"R140Q", "R140Q", SpecificityLevel.VARIANT),
        ],
    }

    # Tumor type patterns (order matters - more specific first)
    # Note: Specific glioma subtypes are listed separately to preserve FDA indication specificity
    TUMOR_PATTERNS = {
        "NSCLC": r"(?:non[- ]?small\s+cell\s+lung\s+cancer|NSCLC(?![A-Z]))",
        "SCLC": r"(?:(?<!non[- ])(?<!non)small\s+cell\s+lung\s+cancer|(?<![N])SCLC(?![A-Z]))",
        "breast cancer": r"(?:breast\s+(?:cancer|carcinoma))",
        "colorectal cancer": r"(?:colorectal\s+(?:cancer|carcinoma)|CRC|colon\s+cancer)",
        "melanoma": r"melanoma",
        # Glioma subtypes - extract specific types when mentioned
        # Use word boundaries to avoid matching "glioma" within "oligodendroglioma"
        "astrocytoma": r"\bastrocytoma\b",
        "oligodendroglioma": r"\boligodendroglioma\b",
        "glioblastoma": r"\bglioblastoma\b",
        "glioma": r"\bglioma\b",  # Generic fallback (only matches standalone "glioma")
        "thyroid cancer": r"(?:thyroid\s+(?:cancer|carcinoma))",
        "gastric cancer": r"(?:gastric\s+(?:cancer|carcinoma)|stomach\s+cancer)",
        "hepatocellular carcinoma": r"(?:hepatocellular\s+carcinoma|HCC|liver\s+cancer)",
        "cholangiocarcinoma": r"(?:cholangiocarcinoma|bile\s+duct\s+cancer)",
        "ovarian cancer": r"(?:ovarian\s+(?:cancer|carcinoma))",
        "pancreatic cancer": r"(?:pancreatic\s+(?:cancer|carcinoma|adenocarcinoma))",
        "AML": r"(?:acute\s+myeloid\s+leukemia|AML)",
        "MDS": r"(?:myelodysplastic\s+syndrome|MDS)",
        "solid tumor": r"(?:solid\s+tumor|solid\s+malignancies)",
    }

    # Combination therapy patterns
    COMBINATION_PATTERNS = [
        r"in\s+combination\s+with\s+([^,\.]+?)(?:\s+for|\s+in|\s+is|\.|,|$)",
        r"combined\s+with\s+([^,\.]+?)(?:\s+for|\s+in|\.|,|$)",
        r"coadministered\s+with\s+([^,\.]+?)(?:\s+for|\s+in|\.|,|$)",
        r"(?:plus|and)\s+([a-z]+(?:mab|nib|lib|tib|tinib|ciclib))",
        r"with\s+([a-z]+(?:mab|nib|lib|tib|tinib|ciclib))\s+(?:for|in)",
    ]

    # Stage patterns
    STAGE_PATTERNS = [
        (r"metastatic", "metastatic"),
        (r"locally\s+advanced", "locally advanced"),
        (r"unresectable", "unresectable"),
        (r"advanced", "advanced"),
        (r"stage\s+(?:III|IV|3|4)", "advanced"),
        (r"grade\s+2", "grade 2"),
        (r"recurrent", "recurrent"),
    ]

    # Line of therapy patterns
    LINE_OF_THERAPY_PATTERNS = [
        (r"first[- ]?line", "first-line"),
        (r"1L", "first-line"),
        (r"second[- ]?line", "second-line"),
        (r"2L", "second-line"),
        (
            r"after\s+prior\s+therapy(?!\s+with)",
            "subsequent",
        ),  # "after prior therapy" (not followed by "with")
        (r"after\s+(?:prior\s+)?(?:treatment|therapy)\s+with", "subsequent"),
        (
            r"(?:whose\s+disease\s+)?(?:has\s+)?progress(?:ed|ion)\s+(?:on|after|following)",
            "post-progression",
        ),
        (r"previously\s+treated", "previously treated"),
        (r"treatment[- ]?na[ïi]ve", "first-line"),
    ]

    GENE_LIST_PATH = (
        Path(__file__).resolve().parents[3] / "data" / "genes" / "gene_list.txt"
    )

    def __init__(self):
        self.gene_patterns = dict(self.GENE_PATTERNS)
        self._load_gene_list()
        self._compile_patterns()

    def _load_gene_list(self):
        """Load genes from gene_list.txt and add any not already in GENE_PATTERNS."""
        logger = logging.getLogger(__name__)
        try:
            with open(self.GENE_LIST_PATH) as f:
                for line in f:
                    gene = line.strip()
                    if gene and gene not in self.gene_patterns:
                        self.gene_patterns[gene] = rf"\b{gene}\b"
        except FileNotFoundError:
            logger.warning(f"Gene list not found at {self.GENE_LIST_PATH}")

    def _compile_patterns(self):
        """Pre-compile regex patterns for performance."""
        self._negation_re = re.compile(
            "|".join(f"({p})" for p in self.NEGATION_PATTERNS), re.IGNORECASE
        )
        self._negative_efficacy_re = re.compile(
            "|".join(f"({p})" for p in self.NEGATIVE_EFFICACY_PATTERNS), re.IGNORECASE
        )
        self._prior_therapy_re = re.compile(
            "|".join(f"({p})" for p in self.PRIOR_THERAPY_PATTERNS), re.IGNORECASE
        )
        # Pre-compile per-gene patterns and build a single combined regex for fast screening
        self._compiled_gene_patterns = {
            gene: re.compile(pattern, re.IGNORECASE)
            for gene, pattern in self.gene_patterns.items()
        }
        # Combined pattern: match any gene symbol as a quick screen
        gene_symbols = sorted(self.gene_patterns.keys(), key=len, reverse=True)
        self._gene_screen_re = re.compile(
            r"\b(?:" + "|".join(re.escape(g) for g in gene_symbols) + r")\b",
            re.IGNORECASE,
        )

    def _normalize_text(self, text: str) -> str:
        """Normalize text by replacing ligatures and special characters."""
        # Common ligatures found in FDA labels (from PDF extraction)
        ligature_map = {
            "\ufb01": "fi",  # ﬁ -> fi
            "\ufb02": "fl",  # ﬂ -> fl
            "\ufb00": "ff",  # ﬀ -> ff
            "\ufb03": "ffi",  # ﬃ -> ffi
            "\ufb04": "ffl",  # ﬄ -> ffl
            "\u2019": "'",  # ' -> '
            "\u2018": "'",  # ' -> '
            "\u201c": '"',  # " -> "
            "\u201d": '"',  # " -> "
            "\u2013": "-",  # – -> -
            "\u2014": "-",  # — -> -
        }
        for ligature, replacement in ligature_map.items():
            text = text.replace(ligature, replacement)
        return text

    def parse_label(self, label_data: dict) -> list[FDABiomarkerEvidence]:
        """
        Parse a complete FDA label from OpenFDA API response.

        Args:
            label_data: Single result from OpenFDA /drug/label.json

        Returns:
            List of FDABiomarkerEvidence objects
        """
        indications = []

        # Extract drug names (normalized to uppercase for consistency)
        openfda = label_data.get("openfda", {})
        drug_name = openfda.get("generic_name", [None])[0]
        brand_name = openfda.get("brand_name", [None])[0]

        if not drug_name:
            drug_name = brand_name or "Unknown"
        # Normalize drug name to uppercase
        drug_name = drug_name.upper() if drug_name else drug_name

        # Parse only Indications and Usage section
        # DO NOT parse clinical_studies - it mentions variants in context of
        # resistance mechanisms, trial eligibility, etc. which are NOT approved indications
        sections_to_parse = [
            ("indications_and_usage", label_data.get("indications_and_usage", [])),
        ]

        for section_name, section_content in sections_to_parse:
            if not section_content:
                continue

            # OpenFDA returns lists of strings
            text = (
                " ".join(section_content)
                if isinstance(section_content, list)
                else section_content
            )

            # Normalize text to handle ligatures from PDF extraction
            text = self._normalize_text(text)

            # Split into indication blocks (numbered sections like "1.1", "1.2")
            indication_blocks = self._split_indication_blocks(text)

            for block in indication_blocks:
                block_indications = self._parse_indication_block(
                    block, drug_name, brand_name, section_name
                )
                indications.extend(block_indications)

        # Deduplicate indications (same gene + requirement + tumor + combination partners + line of therapy)
        # Different combination partners or lines of therapy represent distinct FDA approvals
        seen = set()
        unique_indications = []
        for ind in indications:
            key = (
                ind.gene,
                ind.requirement,
                ind.specificity,
                tuple(sorted(ind.specified_variants)),
                tuple(sorted(ind.tumor_types)),
                tuple(sorted(ind.combination_partners)),
                ind.line_of_therapy,
            )
            if key not in seen:
                seen.add(key)
                unique_indications.append(ind)

        return unique_indications

    def _split_indication_blocks(self, text: str) -> list[str]:
        """Split indication text into separate indication blocks."""

        # Split off "Limitation of Use" sections first - discard them
        limitation_split = re.split(
            r"Limitation(?:s)?\s+of\s+Use\s*:", text, flags=re.IGNORECASE
        )
        if len(limitation_split) > 1:
            text = limitation_split[0]

        # Pattern for numbered sections like "1.1", "( 1.2 )", etc.
        split_pattern = r"(?:^|\s)(?:\(\s*)?(\d+\.\d+)(?:\s*\))?(?:\s|$)"

        # First, try to split by section numbers
        parts = re.split(split_pattern, text)

        if len(parts) > 1:
            # Recombine section numbers with their content
            blocks = []
            i = 0
            while i < len(parts):
                if re.match(r"\d+\.\d+", parts[i].strip()):
                    # This is a section number, combine with next part
                    if i + 1 < len(parts):
                        blocks.append(f"{parts[i]} {parts[i+1]}")
                        i += 2
                    else:
                        i += 1
                else:
                    if parts[i].strip():
                        blocks.append(parts[i])
                    i += 1
            if blocks:
                return blocks

        # Try splitting by bullet points (•, ·, -, *)
        # Split on bullet that's followed by text (not just whitespace)
        bullet_split = re.split(r"(?:^|[\n\r])\s*[•·\-\*]\s+(?=[A-Z])", text)
        if len(bullet_split) > 1:
            # Filter out empty blocks and return
            return [b.strip() for b in bullet_split if b.strip()]

        # Try splitting by "for the treatment of" patterns
        treatment_blocks = re.split(
            r"(?:•|·|\n)\s*(?=for the treatment)", text, flags=re.IGNORECASE
        )
        if len(treatment_blocks) > 1:
            return [b.strip() for b in treatment_blocks if b.strip()]

        # Try splitting by "indicated for:" followed by multiple items
        if re.search(r"indicated\s+for\s*:", text, re.IGNORECASE):
            # Split after the colon by sentence-like boundaries
            post_colon = re.split(r"indicated\s+for\s*:", text, flags=re.IGNORECASE)
            if len(post_colon) > 1:
                items = post_colon[1]
                # Split by periods followed by capital letters (new sentences)
                sentences = re.split(r"\.\s+(?=[A-Z])", items)
                if len(sentences) > 1:
                    return [s.strip() for s in sentences if s.strip()]

        return [text]

    def _parse_indication_block(
        self, text: str, drug_name: str, brand_name: Optional[str], section: str
    ) -> list[FDABiomarkerEvidence]:
        """Parse a single indication block for biomarker associations."""
        indications = []

        # Check for negative efficacy statements first
        if self._is_negative_efficacy(text):
            return []  # Skip blocks that say the drug wasn't effective

        # Quick screen: find which genes are mentioned in this block
        screen_hits = {m.group().upper() for m in self._gene_screen_re.finditer(text)}
        if not screen_hits:
            return indications

        # Only run detailed patterns for genes found in screen
        candidate_genes = {
            gene
            for gene in self.gene_patterns
            if gene.upper() in screen_hits
            or any(s in gene.upper() for s in screen_hits)
        }
        # Also include genes whose custom patterns might match differently (e.g. EGFR matching "epidermal growth factor receptor")
        # Always check genes with custom patterns (non-word-boundary patterns from GENE_PATTERNS class var)
        candidate_genes.update(self.GENE_PATTERNS.keys())

        # Find all genes mentioned
        for gene in candidate_genes:
            compiled = self._compiled_gene_patterns[gene]
            gene_matches = list(compiled.finditer(text))

            for gene_match in gene_matches:
                # Check if this gene mention is in a prior therapy requirement context
                # e.g., "Patients with EGFR should have disease progression on..."
                if self._is_in_prior_therapy_context(text, gene_match.start()):
                    continue  # Skip this gene mention - it's not an indication

                # Check if this gene mention describes mechanism of action
                # e.g., "EGFR antagonist" means drug targets EGFR protein, not mutation
                if self._is_mechanism_of_action_context(text, gene_match.end()):
                    continue  # Skip - this is mechanism of action, not biomarker

                # Check if this gene mention is negated
                requirement = self._check_negation(text, gene_match.start())

                # Find variant specificity
                specificity, variants, codon = self._extract_variant_specificity(
                    text, gene
                )

                # Find tumor types
                tumor_types = self._extract_tumor_types(text)

                # Find tumor stage
                stage = self._extract_stage(text)

                # Find combination partners
                partners = self._extract_combination_partners(text)

                # Find line of therapy
                line = self._extract_line_of_therapy(text)

                indication = FDABiomarkerEvidence(
                    drug_name=drug_name,
                    brand_name=brand_name,
                    gene=gene,
                    requirement=requirement,
                    specificity=specificity,
                    specified_variants=variants,
                    codon=codon,
                    tumor_types=tumor_types,
                    tumor_stage=stage,
                    combination_partners=partners,
                    is_monotherapy=len(partners) == 0,
                    line_of_therapy=line,
                    indication_text=text[:500],  # Store first 500 chars for reference
                    section=section,
                )

                indications.append(indication)

        return indications

    def _check_negation(self, text: str, gene_position: int) -> BiomarkerRequirement:
        """
        Check if the gene mention is negated (e.g., "no EGFR mutation").

        Args:
            text: Full indication text
            gene_position: Character position of the gene mention

        Returns:
            BiomarkerRequirement indicating if biomarker should be present or absent
        """
        # Look at the 100 characters before the gene mention
        context_start = max(0, gene_position - 100)
        context = text[context_start:gene_position].lower()

        # Check for negation patterns
        if self._negation_re.search(context):
            return BiomarkerRequirement.REQUIRED_NEGATIVE

        # Also check for "without X or Y" patterns where gene might be Y
        # e.g., "without EGFR mutation or ALK aberrations"
        extended_context = text[context_start : gene_position + 50].lower()
        if re.search(
            r"without\s+\w+\s+(?:mutation|aberration|alteration)s?\s+or",
            extended_context,
        ):
            return BiomarkerRequirement.REQUIRED_NEGATIVE

        return BiomarkerRequirement.REQUIRED_POSITIVE

    def _is_negative_efficacy(self, text: str) -> bool:
        """Check if this block states the drug was NOT effective."""
        return bool(self._negative_efficacy_re.search(text))

    def _is_prior_therapy_requirement(self, text: str) -> bool:
        """
        Check if this block describes a prior therapy requirement rather than an indication.

        e.g., "Patients with EGFR or ALK genomic tumor aberrations should have disease
        progression on FDA-approved therapy for these aberrations prior to receiving KEYTRUDA"

        This is NOT saying the drug treats EGFR-positive patients - it's saying IF a patient
        has EGFR, they must have failed targeted therapy first.
        """
        return bool(self._prior_therapy_re.search(text))

    def _is_in_prior_therapy_context(self, text: str, gene_position: int) -> bool:
        """
        Check if a gene mention at a specific position is in a prior therapy requirement context.

        Looks at the sentence containing the gene mention to see if it's describing
        prior therapy requirements rather than an indication.
        """
        # Find the sentence containing this gene mention
        # Look backwards for sentence start
        sentence_start = gene_position
        for i in range(gene_position - 1, max(0, gene_position - 200), -1):
            if text[i] in ".!?" and i < gene_position - 1:
                sentence_start = i + 1
                break
        else:
            sentence_start = max(0, gene_position - 200)

        # Look forwards for sentence end
        sentence_end = gene_position
        for i in range(gene_position, min(len(text), gene_position + 200)):
            if text[i] in ".!?":
                sentence_end = i + 1
                break
        else:
            sentence_end = min(len(text), gene_position + 200)

        sentence = text[sentence_start:sentence_end]

        # Check if this sentence contains prior therapy patterns
        return bool(self._prior_therapy_re.search(sentence))

    def _is_mechanism_of_action_context(self, text: str, gene_match_end: int) -> bool:
        """
        Check if a gene mention is describing mechanism of action, not a biomarker requirement.

        e.g., "EGFR antagonist" or "epidermal growth factor receptor (EGFR) antibody"
        means the drug targets EGFR protein - NOT that patients need EGFR mutations.

        Args:
            text: Full indication text
            gene_match_end: Character position where the gene pattern match ends

        Returns:
            True if this is a mechanism-of-action description (should skip this gene)
        """
        # Look at the text following the gene mention (within 60 chars)
        following_text = text[gene_match_end : gene_match_end + 60].lower()

        # Handle patterns like "receptor (EGFR) antagonist" by removing parenthetical abbreviations
        # This regex removes patterns like "(EGFR)" or "( EGFR )"
        following_text = re.sub(
            r"\s*\([A-Z0-9]+\)\s*", " ", following_text, flags=re.IGNORECASE
        )

        # Strip leading whitespace and closing parentheses (for cases where gene is inside parens)
        # e.g., "(EGFR) antagonist" -> after matching EGFR, remaining is ") antagonist"
        following_text = following_text.lstrip(" \t)")

        for pattern in self.MECHANISM_OF_ACTION_PATTERNS:
            if re.match(pattern, following_text, re.IGNORECASE):
                return True
        return False

    def _extract_variant_specificity(
        self, text: str, gene: str
    ) -> tuple[SpecificityLevel, list[str], Optional[str]]:
        """
        Extract variant-level specificity from text.

        Returns:
            (specificity_level, list_of_variants, codon_number)
        """
        variants = []
        codon = None

        # Check gene-specific variant patterns
        if gene in self.VARIANT_PATTERNS:
            for pattern, variant_name, level in self.VARIANT_PATTERNS[gene]:
                if re.search(pattern, text, re.IGNORECASE):
                    variants.append(variant_name)

                    # Extract codon if it's a variant
                    codon_match = re.search(r"[A-Z]?(\d+)[A-Z]?", variant_name)
                    if codon_match:
                        codon = codon_match.group(1)

        if variants:
            # Check if we have full variants or just codons
            has_full_variant = any(
                re.match(r"^[A-Z]\d+[A-Z]$", v)
                or "del" in v.lower()
                or "ins" in v.lower()
                for v in variants
            )
            if has_full_variant:
                # Remove codon-only entries if we have full variants
                # e.g., remove "G12" if we have "G12C"
                full_variants = [
                    v
                    for v in variants
                    if re.match(r"^[A-Z]\d+[A-Z]$", v)
                    or "del" in v.lower()
                    or "ins" in v.lower()
                ]
                return SpecificityLevel.VARIANT, full_variants, codon
            else:
                return SpecificityLevel.CODON, variants, codon

        # Check for generic mutation mentions
        gene_pattern = self.gene_patterns.get(gene, gene)
        if re.search(
            rf"{gene_pattern}[- ](?:positive|mutant|mutated|altered|mutation)",
            text,
            re.IGNORECASE,
        ):
            return SpecificityLevel.GENE, [], None

        return SpecificityLevel.GENE, [], None

    def _extract_tumor_types(self, text: str) -> list[str]:
        """Extract tumor types from indication text."""
        tumor_types = []
        for tumor_name, tumor_pattern in self.TUMOR_PATTERNS.items():
            if re.search(tumor_pattern, text, re.IGNORECASE):
                tumor_types.append(tumor_name)
        return tumor_types

    def _extract_stage(self, text: str) -> Optional[str]:
        """Extract tumor stage/setting from indication text."""
        for pattern, stage_name in self.STAGE_PATTERNS:
            if re.search(pattern, text, re.IGNORECASE):
                return stage_name
        return None

    def _extract_combination_partners(self, text: str) -> list[str]:
        """Extract combination therapy partners from indication text."""
        partners = []
        for pattern in self.COMBINATION_PATTERNS:
            matches = re.findall(pattern, text, re.IGNORECASE)
            for match in matches:
                # Clean up the partner name
                partner = match.strip().lower()
                # Remove common suffixes/prefixes
                partner = re.sub(r"^(?:and|or|with)\s+", "", partner)
                partner = re.sub(r"\s+(?:for|in|after|until).*$", "", partner)
                if partner and len(partner) > 2:
                    # Split compound partners like "palbociclib and fulvestrant" into separate drugs
                    # Handle both "and" and "," separators
                    sub_partners = re.split(r"\s+and\s+|,\s*", partner)
                    for sub in sub_partners:
                        sub = sub.strip()
                        if sub and len(sub) > 2:
                            partners.append(sub)
        return list(set(partners))  # Deduplicate

    def _extract_line_of_therapy(self, text: str) -> Optional[str]:
        """Extract line of therapy from indication text."""
        for pattern, line_name in self.LINE_OF_THERAPY_PATTERNS:
            if re.search(pattern, text, re.IGNORECASE):
                return line_name
        return None


def match_variant_to_indications(
    indications: list[FDABiomarkerEvidence],
    query_gene: str,
    query_variant: str,
    query_tumor: Optional[str] = None,
) -> list[dict]:
    """
    Match a user's variant query against parsed FDA indications.

    Args:
        indications: List of parsed FDABiomarkerEvidence objects
        query_gene: Gene symbol (e.g., "EGFR")
        query_variant: Variant notation (e.g., "L858R", "T790M")
        query_tumor: Optional tumor type filter

    Returns:
        List of match results with full context
    """
    results = []

    for ind in indications:
        match_result = ind.matches_variant(query_gene, query_variant)

        # Skip indications for different genes - they're not relevant
        if (
            match_result["match_type"] is None
            and "Different gene" in match_result["reason"]
        ):
            continue

        # Filter by tumor if specified
        if query_tumor and ind.tumor_types:
            tumor_match = any(
                query_tumor.lower() in t.lower() or t.lower() in query_tumor.lower()
                for t in ind.tumor_types
            )
            if not tumor_match:
                continue

        results.append(
            {
                "drug": ind.drug_name,
                "brand": ind.brand_name,
                "matches": match_result["matches"],
                "match_type": match_result["match_type"],
                "reason": match_result["reason"],
                "requirement": ind.requirement.value,
                "specificity": ind.specificity.value,
                "specified_variants": ind.specified_variants,
                "tumor_types": ind.tumor_types,
                "tumor_stage": ind.tumor_stage,
                "combination_partners": ind.combination_partners,
                "is_monotherapy": ind.is_monotherapy,
                "line_of_therapy": ind.line_of_therapy,
            }
        )

    # Sort: exact matches first, then excluded (so user sees what NOT to use)
    results.sort(
        key=lambda x: (
            not x["matches"],  # Matches first
            x["match_type"] != "exact",  # Exact matches before partial
            x["match_type"] == "excluded",  # Exclusions last among non-matches
        )
    )

    return results


# ============================================================================
# OpenFDA Integration
# ============================================================================

import requests
from typing import Optional


class OpenFDAClient:
    """Client for fetching FDA labels from OpenFDA API."""

    BASE_URL = "https://api.fda.gov/drug/label.json"

    def fetch_label_by_drug(self, drug_name: str) -> Optional[dict]:
        """
        Fetch FDA label by drug name (generic or brand).

        Args:
            drug_name: Generic or brand name of the drug

        Returns:
            Label data dict or None if not found
        """
        search = (
            f'openfda.generic_name:"{drug_name}" OR openfda.brand_name:"{drug_name}"'
        )
        params = {"search": search, "limit": 1}

        try:
            resp = requests.get(self.BASE_URL, params=params, timeout=10)
            resp.raise_for_status()
            data = resp.json()

            if data.get("results"):
                return data["results"][0]
        except Exception as e:
            print(f"Error fetching label for {drug_name}: {e}")

        return None

    def fetch_labels_by_gene(self, gene: str, limit: int = 100) -> list[dict]:
        """
        Fetch all FDA labels that mention a specific gene/biomarker.

        Searches using both the gene symbol (e.g., "EGFR") and full names
        (e.g., "epidermal growth factor receptor") to ensure comprehensive coverage.

        Args:
            gene: Gene symbol (e.g., "EGFR", "BRAF")
            limit: Maximum number of labels to return

        Returns:
            List of label data dicts
        """
        from oncomind.config.constants import GENE_FULL_NAMES

        # Build search terms: gene symbol + full names
        search_terms = [gene]
        gene_upper = gene.upper()
        if gene_upper in GENE_FULL_NAMES:
            search_terms.extend(GENE_FULL_NAMES[gene_upper])

        # Build OR query for all search terms
        # Quote each term and join with OR
        search_parts = [f'indications_and_usage:"{term}"' for term in search_terms]
        search = " OR ".join(search_parts)
        params = {"search": search, "limit": limit}

        try:
            resp = requests.get(self.BASE_URL, params=params, timeout=30)
            resp.raise_for_status()
            data = resp.json()
            return data.get("results", [])
        except Exception as e:
            print(f"Error fetching labels for {gene}: {e}")

        return []


def get_fda_approved_drugs_for_variant(
    gene: str, variant: str, tumor_type: Optional[str] = None
) -> list[dict]:
    """
    High-level function to find FDA-approved drugs for a specific variant.

    Args:
        gene: Gene symbol (e.g., "EGFR")
        variant: Variant notation (e.g., "T790M", "L858R")
        tumor_type: Optional tumor type filter (e.g., "NSCLC")

    Returns:
        List of matching drug information with approval details
    """
    client = OpenFDAClient()
    parser = FDALabelParser()

    # Fetch all labels mentioning this gene
    labels = client.fetch_labels_by_gene(gene)

    all_matches = []

    for label_data in labels:
        # Parse the label
        indications = parser.parse_label(label_data)

        # Match against the query
        matches = match_variant_to_indications(indications, gene, variant, tumor_type)

        for match in matches:
            if match["matches"]:
                all_matches.append(match)
            elif match["match_type"] == "excluded":
                # Include exclusions so user knows what NOT to use
                match["_is_exclusion"] = True
                all_matches.append(match)

    # Deduplicate by drug name
    seen_drugs = set()
    unique_matches = []
    for match in all_matches:
        if match["drug"] not in seen_drugs:
            seen_drugs.add(match["drug"])
            unique_matches.append(match)

    return unique_matches
