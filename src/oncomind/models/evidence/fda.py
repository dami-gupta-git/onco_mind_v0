import re

from pydantic import BaseModel, Field

from oncomind.models.evidence.base import EvidenceItemBase


def extract_combination_partners(indication_text: str | None) -> list[str]:
    """Extract combination drug partners from indication text."""
    if not indication_text:
        return []

    patterns = [
        r"in combination with\s+([^,\.]+?)(?:\s+for|\s+in|\.|,|$)",
        r"combined with\s+([^,\.]+?)(?:\s+for|\s+in|\.|,|$)",
        r"coadministered with\s+([^,\.]+?)(?:\s+for|\s+in|\.|,|$)",
    ]

    partners = []
    for pattern in patterns:
        matches = re.findall(pattern, indication_text, re.IGNORECASE)
        partners.extend([m.strip() for m in matches])

    # Handle "X and Y" or "X or Y" in partner string
    expanded = []
    for p in partners:
        parts = re.split(r"\s+(?:and|or|\+|with)\s+", p, flags=re.IGNORECASE)
        expanded.extend([x.strip() for x in parts])

    # Clean up artifacts
    cleaned = []
    for p in expanded:
        # Remove leading "either", "without", "or without"
        p = re.sub(r"^(either|without|or without)\s+", "", p, flags=re.IGNORECASE)
        # Remove trailing "or without X"
        p = re.sub(r"\s+(or without|without).*$", "", p, flags=re.IGNORECASE)
        p = p.strip()
        if p:
            cleaned.append(p)

    # Filter to only drug-like names
    drug_suffixes = (
        'mab', 'nib', 'lib', 'tib', 'tin', 'zumab', 'ximab', 'mumab',
        'platin', 'taxel', 'mustine', 'rubicin', 'mycin', 'tinib',
        'ciclib', 'lisib', 'rafenib', 'metinib', 'sertib', 'parin',
        'mide', 'uracil', 'citabine', 'trexed', 'fosfamide', 'imod',
        'olimus', 'fulvestrant', 'pemetrexed', 'asertib',
    )

    filtered = []
    for p in cleaned:
        p_lower = p.lower()
        # Skip generic terms
        skip_terms = (
            'chemotherapy', 'radiotherapy', 'radiation', 'therapy',
            'platinum-based', 'platinum', 'chemoradiotherapy',
            'neoadjuvant treatment', 'fluoropyrimidine-based',
            'platinum-containing', 'single agent',
        )
        if any(term in p_lower for term in skip_terms):
            continue
        # Check suffix
        if any(p_lower.endswith(suffix) for suffix in drug_suffixes):
            filtered.append(p)

    # Deduplicate while preserving order
    seen = set()
    unique = []
    for p in filtered:
        p_lower = p.lower()
        if p_lower not in seen:
            seen.add(p_lower)
            unique.append(p)

    return unique


class BiomarkerMatch(BaseModel):
    """Result of matching a patient's biomarker against FDA indication.

    Populated by match_fda_approval() from fda_processor.py.
    """
    matched: bool = False  # Whether the variant is covered by FDA approval
    match_level: str | None = None  # "variant", "codon", "gene", "contraindicated", or None
    tumor_matched: bool | None = None  # Whether the tumor type matches FDA indication
    tumor_match_type: str | None = None  # "exact", "pan_cancer", or None
    combination_partners: list[str] = Field(default_factory=list)  # Partner drugs from indication


class ClinicalStudyEvidence(BaseModel):
    """Structured clinical study data from FDA label."""
    trial_name: str | None = None
    nct_id: str | None = None
    patients_n: int | None = None
    pfs_months_treatment: float | None = None
    pfs_months_control: float | None = None
    hazard_ratio: float | None = None
    hazard_ratio_ci: tuple[float, float] | None = None
    orr_treatment: float | None = None
    orr_control: float | None = None
    biomarker_breakdown: dict[str, float] | None = None


class MechanismEvidence(BaseModel):
    """Mechanism of action data from FDA label."""
    targets: list[str] = Field(default_factory=list)
    mechanism: str | None = None
    preclinical: str | None = None


class AdverseReactionsEvidence(BaseModel):
    """Adverse reactions data from FDA label."""
    common_toxicities: list[tuple[str, float]] = Field(default_factory=list)
    serious_rate: float | None = None
    discontinuation_rate: float | None = None


class FDALabelEvidence(EvidenceItemBase):
    """Complete FDA drug label data for display in FDA tab.

    This contains the structured data extracted from OpenFDA drug labels,
    including clinical trial results, mechanism of action, and adverse reactions.
    """
    drug: str
    gene: str
    set_id: str | None = None  # FDA label unique identifier (UUID)
    version: str | None = None  # FDA label version number
    brand_name: str | None = None
    generic_name: str | None = None
    manufacturer: str | None = None
    indications_and_usage: str | None = None
    initial_approval_date: str | None = None
    effective_time: str | None = None  # Label revision date (YYYY-MM-DD)
    approved_indications: list[str] = Field(default_factory=list)  # List of approved diseases

    # Structured data
    clinical_studies: ClinicalStudyEvidence | None = None
    mechanism_of_action: MechanismEvidence | None = None
    adverse_reactions: AdverseReactionsEvidence | None = None
    last_label_update: str | None = None
    update_reason: str | None = None

    # Raw text for advanced use
    clinical_studies_text: str | None = None
    mechanism_of_action_text: str | None = None
    adverse_reactions_text: str | None = None

    # Match result from match_fda_approval - populated when matching patient variant
    biomarker_match: BiomarkerMatch | None = None

    def to_approval(self) -> "FDAApproval":
        """Convert this label to an FDAApproval for LLM consumption.

        Returns a simplified FDAApproval object suitable for LLM prompts,
        preserving the locus_variant_match and biomarker_match from this label.
        """
        # Extract flattened biomarker_match fields
        biomarker_matched = False
        biomarker_match_level = None
        biomarker_tumor_matched = None
        biomarker_tumor_match_type = None
        combination_partners: list[str] = []

        if self.biomarker_match:
            biomarker_matched = self.biomarker_match.matched
            biomarker_match_level = self.biomarker_match.match_level
            biomarker_tumor_matched = self.biomarker_match.tumor_matched
            biomarker_tumor_match_type = self.biomarker_match.tumor_match_type
            combination_partners = self.biomarker_match.combination_partners

        return FDAApproval(
            drug_name=self.drug,  # Already has combination name (e.g., "capivasertib + fulvestrant")
            brand_name=self.brand_name,
            generic_name=self.generic_name,
            indication=self.indications_and_usage,
            approval_date=self.initial_approval_date,
            gene=self.gene,
            locus_variant_match=self.locus_variant_match,
            cancer_type_match=self.cancer_type_match,
            # Flattened biomarker match fields
            biomarker_matched=biomarker_matched,
            biomarker_match_level=biomarker_match_level,
            biomarker_tumor_matched=biomarker_tumor_matched,
            biomarker_tumor_match_type=biomarker_tumor_match_type,
            combination_partners=combination_partners,
        )


class FDAApproval(EvidenceItemBase):
    """FDA drug approval information.

    Inherits from EvidenceItemBase which provides:
    - locus_variant_match: EvidenceLevel tracking variant/gene level match
    - cancer_type_match: EvidenceLevel tracking cancer specificity
    - locus_match: Computed field returning 'variant', 'codon', or 'gene'
    - tumor_match: Computed field returning True/False/None for cancer match
    """

    drug_name: str | None = None
    brand_name: str | None = None
    generic_name: str | None = None
    indication: str | None = None
    approval_date: str | None = None
    gene: str | None = None
    variant_in_indications: bool = False
    variant_in_clinical_studies: bool = False
    # Additional fields for clinical context
    companion_diagnostic: str | None = None
    black_box_warning: str | None = None
    dosing_for_variant: str | None = None
    # Drug response association (Responsive, Resistant, etc.)
    # Note: FDA labels from OpenFDA don't contain this - it comes from CGI/CIViC
    # This field is None for FDA label-derived approvals
    association: str | None = None

    # Flattened BiomarkerMatch fields - populated from FDALabelEvidence.biomarker_match
    biomarker_matched: bool = False  # Whether the queried variant is covered by this approval
    biomarker_match_level: str | None = None  # "variant", "codon", "gene", or None
    biomarker_tumor_matched: bool | None = None  # Whether tumor type matches FDA indication
    biomarker_tumor_match_type: str | None = None  # "exact", "pan_cancer", or None
    combination_partners: list[str] = Field(default_factory=list)  # Partner drugs from indication

    def extract_approved_variant(self) -> str | None:
        """Extract the specific variant(s) the drug is approved for from indication text.

        Looks for patterns like 'KRAS G12C-mutated', 'EGFR exon 19 deletion', etc.
        Falls back to known drug-variant associations when text parsing fails.

        Returns:
            The variant string (e.g., 'G12C', 'exon 19 del') or None if not found
        """
        if not self.gene:
            return None

        gene_upper = self.gene.upper()

        # Try to extract from indication text first
        if self.indication:
            indication_upper = self.indication.upper()

            # Pattern 1: Gene + specific variant (e.g., "KRAS G12C", "BRAF V600E")
            # Matches: KRAS G12C, EGFR L858R, BRAF V600E, etc.
            variant_pattern = rf'{gene_upper}\s+([A-Z]\d+[A-Z])'
            match = re.search(variant_pattern, indication_upper)
            if match:
                return match.group(1)

            # Pattern 2: Gene + exon notation (e.g., "EGFR exon 19 deletion")
            exon_pattern = rf'{gene_upper}\s+(EXON\s*\d+\s*(?:DELETION|DEL|INSERTION|INS))'
            match = re.search(exon_pattern, indication_upper)
            if match:
                return match.group(1).lower().replace('deletion', 'del').replace('insertion', 'ins')

            # Pattern 3: Gene-mutated without specific variant (e.g., "BRCA-mutated")
            # This means any mutation in the gene, not a specific variant
            mutated_pattern = rf'{gene_upper}[\s-]*(MUTATED|MUTATION)'
            if re.search(mutated_pattern, indication_upper):
                return "any mutation"

        # Fallback: Known variant-specific drugs
        # These drugs are FDA-approved ONLY for specific variants, not the whole gene
        drug_name = (self.generic_name or self.drug_name or "").lower()
        known_variant_drugs = {
            # KRAS G12C inhibitors - only approved for G12C
            ("kras", "sotorasib"): "G12C",
            ("kras", "adagrasib"): "G12C",
        }

        for (gene, drug), variant in known_variant_drugs.items():
            if gene_upper == gene.upper() and drug in drug_name:
                return variant

        return None

    def parse_indication_for_tumor(self, tumor_type: str) -> dict:
        """Parse FDA indication text to extract line-of-therapy and approval type for a specific tumor."""
        if not self.indication or not tumor_type:
            return {
                'tumor_match': False,
                'line_of_therapy': 'unspecified',
                'approval_type': 'unspecified',
                'indication_excerpt': ''
            }

        indication_lower = self.indication.lower()
        tumor_lower = tumor_type.lower()

        # Check for tumor type match (flexible matching)
        tumor_keywords = {
            'colorectal': ['colorectal', 'colon', 'rectal', 'crc', 'mcrc'],
            'melanoma': ['melanoma'],
            'lung': ['lung', 'nsclc', 'non-small cell'],
            'breast': ['breast'],
            'thyroid': ['thyroid', 'atc', 'anaplastic thyroid'],
            'gist': ['gist', 'gastrointestinal stromal tumor', 'gastrointestinal stromal'],
            'gastrointestinal stromal tumor': ['gist', 'gastrointestinal stromal tumor', 'gastrointestinal stromal'],
            'bladder': ['bladder', 'urothelial', 'transitional cell'],
            'bladder cancer': ['bladder', 'urothelial', 'transitional cell', 'urothelial carcinoma'],
            'urothelial': ['urothelial', 'bladder', 'transitional cell'],
            'cholangiocarcinoma': ['cholangiocarcinoma', 'bile duct', 'biliary'],
            # Myeloproliferative neoplasms - these are DEFINED by MPL/JAK2/CALR mutations
            # The FDA labels say "myelofibrosis" or "polycythemia vera" but patients present
            # with a diagnosis of "myeloproliferative neoplasm" containing these mutations
            'myeloproliferative neoplasm': ['myelofibrosis', 'polycythemia vera', 'myeloproliferative', 'mpn'],
            'myeloproliferative': ['myelofibrosis', 'polycythemia vera', 'myeloproliferative', 'mpn'],
            'mpn': ['myelofibrosis', 'polycythemia vera', 'myeloproliferative', 'mpn'],
            'myelofibrosis': ['myelofibrosis', 'myeloproliferative'],
            'polycythemia vera': ['polycythemia vera', 'myeloproliferative'],
        }

        tumor_match = False
        matched_section = ""

        # Priority 0: Detect TUMOR-AGNOSTIC MSI-H/dMMR approvals
        # These apply to ANY solid tumor (endometrial, pancreatic, ovarian, etc.)
        # FDA label says "MSI-H or dMMR solid tumors" or "MSI-H or dMMR Cancer"
        # The [FDA APPROVED FOR MSI-H...] prefix indicates this is a tumor-agnostic approval
        msi_h_tumor_agnostic_patterns = [
            'fda approved for msi-h',
            'fda approved for dmmr',
            'microsatellite instability-high',
            'mismatch repair deficient',
        ]
        is_msi_h_approval = any(p in indication_lower for p in msi_h_tumor_agnostic_patterns)

        # For MSI-H/dMMR approvals, check if approval is tumor-agnostic (applies to all solid tumors)
        # vs tumor-specific (e.g., "MSI-H colorectal cancer" only applies to CRC)
        if is_msi_h_approval:
            # Check if this is a tumor-agnostic approval (no specific tumor mentioned with MSI-H)
            # Look for phrases like "MSI-H solid tumors" or "MSI-H Cancer" without a specific site
            tumor_agnostic_phrases = [
                'msi-h or dmmr cancer',
                'msi-h cancer',
                'dmmr cancer',
                'msi-h solid tumor',
                'dmmr solid tumor',
                'msi-h or mismatch repair deficient cancer',
                'microsatellite instability-high or mismatch repair deficient cancer',
            ]
            is_tumor_agnostic = any(p in indication_lower for p in tumor_agnostic_phrases)

            if is_tumor_agnostic:
                # This is a tumor-agnostic approval - applies to ANY solid tumor
                # Including endometrial, pancreatic, ovarian, gastric, etc.
                # Extract the MSI-H section as the matched section
                for pattern in ['[fda approved for msi-h', '[fda approved for dmmr']:
                    if pattern in indication_lower:
                        idx = indication_lower.find(pattern)
                        bracket_end = self.indication.find(']', idx)
                        if bracket_end > 0:
                            matched_section = self.indication[idx:bracket_end + 1]
                            tumor_match = True
                            break

                if not tumor_match:
                    # Fallback: find MSI-H mention in indication
                    for pattern in msi_h_tumor_agnostic_patterns:
                        if pattern in indication_lower:
                            idx = indication_lower.find(pattern)
                            start = max(0, idx - 50)
                            end = min(len(self.indication), idx + 300)
                            matched_section = self.indication[start:end]
                            tumor_match = True
                            break

        # Priority 1: If indication has a variant-specific section at the start (from fda.py),
        # use that section for line-of-therapy detection. This handles cases like TAGRISSO
        # where T790M has its own later-line indication separate from L858R/exon19del first-line.
        if not tumor_match and self.indication.startswith('[FDA APPROVED FOR'):
            # Extract the variant-specific section
            bracket_end = self.indication.find(']')
            if bracket_end > 0:
                variant_section = self.indication[:bracket_end + 1]
                # Check if this variant section mentions the tumor type
                variant_section_lower = variant_section.lower()
                for key, keywords in tumor_keywords.items():
                    if any(kw in tumor_lower for kw in keywords):
                        if any(kw in variant_section_lower for kw in keywords):
                            tumor_match = True
                            matched_section = variant_section
                            break
                if not tumor_match and tumor_lower in variant_section_lower:
                    tumor_match = True
                    matched_section = variant_section

        # Priority 2: Standard tumor type matching in full indication
        if not tumor_match:
            tumor_keys = []
            for key, keywords in tumor_keywords.items():
                if any(kw in tumor_lower for kw in keywords):
                    tumor_keys = keywords
                    break
            if not tumor_keys:
                tumor_keys = [tumor_lower]

            for kw in tumor_keys:
                if kw in indication_lower:
                    tumor_match = True
                    idx = indication_lower.find(kw)
                    start = max(0, idx - 50)
                    next_section_markers = [
                        'non-small cell lung cancer',
                        'nsclc)',
                        'melanoma •',
                        'breast cancer',
                        'thyroid cancer',
                        'limitations of use',
                        '1.1 braf',
                        '1.2 braf',
                        '1.3 braf',
                        '1.4 ',
                    ]
                    end = len(self.indication)
                    for next_sec in next_section_markers:
                        next_idx = indication_lower.find(next_sec, idx + len(kw) + 100)
                        if next_idx > idx and next_idx < end:
                            end = next_idx
                    matched_section = self.indication[start:end]
                    break

        if not tumor_match:
            return {
                'tumor_match': False,
                'line_of_therapy': 'unspecified',
                'approval_type': 'unspecified',
                'indication_excerpt': ''
            }

        later_line_phrases = [
            'after prior therapy',
            'after progression',
            'following progression',
            'following recurrence',
            'has progressed',  # "whose disease has progressed on or after"
            'progressed on or after',
            'second-line',
            'second line',
            'third-line',
            'third line',
            'previously treated',
            'refractory',
            'who have failed',
            'after failure',
            'following prior',
            'disease progression',
        ]

        first_line_phrases = [
            'first-line',
            'first line',
            'frontline',
            'initial treatment',
            'treatment-naive',
            'previously untreated',
        ]

        matched_lower = matched_section.lower()
        line_of_therapy = 'unspecified'

        for phrase in later_line_phrases:
            if phrase in matched_lower:
                line_of_therapy = 'later-line'
                break

        if line_of_therapy == 'unspecified':
            for phrase in first_line_phrases:
                if phrase in matched_lower:
                    line_of_therapy = 'first-line'
                    break

        approval_type = 'full'
        accelerated_phrases = [
            'accelerated approval',
            'approved under accelerated',
            'contingent upon verification',
            'confirmatory trial',
        ]

        for phrase in accelerated_phrases:
            if phrase in matched_lower:
                approval_type = 'accelerated'
                break

        return {
            'tumor_match': True,
            'line_of_therapy': line_of_therapy,
            'approval_type': approval_type,
            'indication_excerpt': matched_section[:300]
        }

    def extract_indication_cancer_type(self) -> str | None:
        """Extract the primary cancer type from the FDA indication text.

        Returns:
            The cancer type string (e.g., "ovarian cancer", "melanoma") or None if not found
        """
        if not self.indication:
            return None

        indication_lower = self.indication.lower()

        # Common cancer type patterns to look for
        # Order matters - more specific patterns first
        cancer_patterns = [
            # Specific cancer types with modifiers
            ("low-grade serous ovarian cancer", "ovarian cancer"),
            ("high-grade serous ovarian cancer", "ovarian cancer"),
            ("serous ovarian cancer", "ovarian cancer"),
            ("epithelial ovarian cancer", "ovarian cancer"),
            ("ovarian cancer", "ovarian cancer"),
            ("non-small cell lung cancer", "NSCLC"),
            ("nsclc", "NSCLC"),
            ("small cell lung cancer", "small cell lung cancer"),
            ("lung cancer", "lung cancer"),
            ("metastatic melanoma", "melanoma"),
            ("melanoma", "melanoma"),
            ("colorectal cancer", "colorectal cancer"),
            ("colon cancer", "colorectal cancer"),
            ("rectal cancer", "colorectal cancer"),
            ("breast cancer", "breast cancer"),
            ("pancreatic cancer", "pancreatic cancer"),
            ("pancreatic adenocarcinoma", "pancreatic cancer"),
            ("thyroid cancer", "thyroid cancer"),
            ("anaplastic thyroid cancer", "thyroid cancer"),
            ("hepatocellular carcinoma", "hepatocellular carcinoma"),
            ("liver cancer", "hepatocellular carcinoma"),
            ("renal cell carcinoma", "renal cell carcinoma"),
            ("kidney cancer", "renal cell carcinoma"),
            ("bladder cancer", "bladder cancer"),
            ("urothelial carcinoma", "bladder cancer"),
            ("gastric cancer", "gastric cancer"),
            ("stomach cancer", "gastric cancer"),
            ("esophageal cancer", "esophageal cancer"),
            ("glioblastoma", "glioblastoma"),
            ("brain cancer", "brain cancer"),
            ("prostate cancer", "prostate cancer"),
            ("endometrial cancer", "endometrial cancer"),
            ("uterine cancer", "endometrial cancer"),
            ("cervical cancer", "cervical cancer"),
            ("head and neck cancer", "head and neck cancer"),
            ("squamous cell carcinoma of the head and neck", "head and neck cancer"),
            ("cholangiocarcinoma", "cholangiocarcinoma"),
            ("bile duct cancer", "cholangiocarcinoma"),
            ("gastrointestinal stromal tumor", "GIST"),
            ("gist", "GIST"),
            ("acute myeloid leukemia", "AML"),
            ("aml", "AML"),
            ("chronic myeloid leukemia", "CML"),
            ("cml", "CML"),
            ("multiple myeloma", "multiple myeloma"),
            ("lymphoma", "lymphoma"),
            ("myelofibrosis", "myelofibrosis"),
            ("polycythemia vera", "polycythemia vera"),
            ("myeloproliferative", "myeloproliferative neoplasm"),
            # Tumor-agnostic
            ("solid tumor", "solid tumors (pan-cancer)"),
            ("msi-h", "MSI-H tumors (pan-cancer)"),
            ("dmmr", "dMMR tumors (pan-cancer)"),
        ]

        for pattern, cancer_type in cancer_patterns:
            if pattern in indication_lower:
                return cancer_type

        return None

