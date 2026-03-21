"""Integration tests for clinical trial NCT ID / title / phase alignment.

These tests verify that NCT IDs returned by search are not misaligned with their
titles, phases, or statuses. They were written in response to a report where
NCT03126916 and NCT05770037 were flagged as possibly having stale/outdated
ID→title/status mappings.

Ground truth for each NCT ID is confirmed via a direct live API call and locked
in here. If a trial is amended, these tests will catch the drift.

Tests here hit the live ClinicalTrials.gov v2 API.
"""

import pytest
import httpx

from oncomind.api.clinicaltrials import ClinicalTrialsClient


# ---------------------------------------------------------------------------
# Ground truth confirmed by direct API call on 2026-03-21
# ---------------------------------------------------------------------------

# NCT03126916 — amended trial (lorlatinib replaced crizotinib)
# Title: "Testing the Addition of 131I-MIBG or Lorlatinib to Intensive Therapy
#          in People With High-Risk Neuroblastoma (NBL)"
# Phase: PHASE3  Status: ACTIVE_NOT_RECRUITING
NCT03126916 = "NCT03126916"
NCT03126916_TITLE_FRAGMENT = "neuroblastoma"
NCT03126916_KNOWN_DRUG = "lorlatinib"
NCT03126916_EXPECTED_PHASE = "PHASE3"

# NCT05770037 — DETERMINE trial, Arm 01: Alectinib in ALK+ cancers
# Title: "DETERMINE Trial Treatment Arm 01: Alectinib in Adult, Paediatric and
#          Teenage/Young Adult Patients With ALK Positive Cancers"
# Phase: PHASE2 (phases[0]; trial is listed as PHASE2/PHASE3)
# Status: RECRUITING
NCT05770037 = "NCT05770037"
NCT05770037_TITLE_FRAGMENT = "alectinib"
NCT05770037_KNOWN_DRUG = "alectinib"
NCT05770037_EXPECTED_PHASE = "PHASE2"
NCT05770037_EXPECTED_STATUS = "RECRUITING"


def _fetch_nct_direct(nct_id: str) -> dict:
    """Fetch a single study directly from the ClinicalTrials.gov v2 API."""
    resp = httpx.get(
        f"https://clinicaltrials.gov/api/v2/studies/{nct_id}",
        params={"format": "json"},
        timeout=15,
    )
    resp.raise_for_status()
    return resp.json()


class TestNCT03126916DirectLookup:
    """Verify NCT03126916 title/phase against live API ground truth.

    Ground truth confirmed 2026-03-21:
      Title: "Testing the Addition of 131I-MIBG or Lorlatinib to Intensive
               Therapy in People With High-Risk Neuroblastoma (NBL)"
      Phase: PHASE3
    """

    @pytest.mark.integration
    def test_nct03126916_title_and_phase_direct(self):
        """Direct API lookup: title contains neuroblastoma + lorlatinib, phase is PHASE3."""
        data = _fetch_nct_direct(NCT03126916)
        proto = data["protocolSection"]
        title = (
            proto["identificationModule"].get("briefTitle")
            or proto["identificationModule"].get("officialTitle")
            or ""
        ).lower()
        phases = proto.get("designModule", {}).get("phases", [])
        phase = phases[0] if phases else None

        assert NCT03126916_TITLE_FRAGMENT in title, (
            f"{NCT03126916} title should contain {NCT03126916_TITLE_FRAGMENT!r}, got: {title!r}"
        )
        assert NCT03126916_KNOWN_DRUG in title, (
            f"{NCT03126916} title should mention {NCT03126916_KNOWN_DRUG!r}, got: {title!r}"
        )
        assert phase == NCT03126916_EXPECTED_PHASE, (
            f"{NCT03126916} expected phase {NCT03126916_EXPECTED_PHASE!r}, got {phase!r}"
        )

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_nct03126916_search_title_matches_direct_lookup(self):
        """Title returned by search_trials() must match title from direct lookup.

        This is the core alignment check: if search returns a different title for
        the same NCT ID, it indicates an ID/data mismatch in parsing.
        """
        direct = _fetch_nct_direct(NCT03126916)
        expected_title = (
            direct["protocolSection"]["identificationModule"].get("briefTitle")
            or direct["protocolSection"]["identificationModule"].get("officialTitle")
        )

        async with ClinicalTrialsClient() as client:
            trials = await client.search_trials(
                gene="ALK",
                tumor_type="Neuroblastoma",
                recruiting_only=False,
                max_results=50,
            )

        nct_ids = [t.nct_id for t in trials]
        if NCT03126916 not in nct_ids:
            pytest.skip(
                f"{NCT03126916} not returned in ALK Neuroblastoma search "
                f"(may be status-excluded). IDs returned: {nct_ids}"
            )

        trial = next(t for t in trials if t.nct_id == NCT03126916)
        assert trial.title == expected_title, (
            f"{NCT03126916} title from search ({trial.title!r}) does not match "
            f"direct lookup ({expected_title!r}) — ID/title misalignment in parsing"
        )


class TestNCT05770037DirectLookup:
    """Verify NCT05770037 title/phase/status against live API ground truth.

    Ground truth confirmed 2026-03-21:
      Title: "DETERMINE Trial Treatment Arm 01: Alectinib in Adult, Paediatric
               and Teenage/Young Adult Patients With ALK Positive Cancers"
      Phase: PHASE2/PHASE3 (phases[0] == PHASE2)
      Status: RECRUITING
    """

    @pytest.mark.integration
    def test_nct05770037_title_phase_status_direct(self):
        """Direct API lookup: title contains alectinib, phase includes PHASE2, status RECRUITING."""
        data = _fetch_nct_direct(NCT05770037)
        proto = data["protocolSection"]
        title = (
            proto["identificationModule"].get("briefTitle")
            or proto["identificationModule"].get("officialTitle")
            or ""
        ).lower()
        phases = proto.get("designModule", {}).get("phases", [])
        phase = phases[0] if phases else None
        status = proto["statusModule"].get("overallStatus")

        assert NCT05770037_TITLE_FRAGMENT in title, (
            f"{NCT05770037} title should contain {NCT05770037_TITLE_FRAGMENT!r}, got: {title!r}"
        )
        assert NCT05770037_KNOWN_DRUG in title, (
            f"{NCT05770037} title should mention {NCT05770037_KNOWN_DRUG!r}, got: {title!r}"
        )
        assert phase == NCT05770037_EXPECTED_PHASE, (
            f"{NCT05770037} expected phase {NCT05770037_EXPECTED_PHASE!r}, got {phase!r}"
        )
        assert status == NCT05770037_EXPECTED_STATUS, (
            f"{NCT05770037} expected status {NCT05770037_EXPECTED_STATUS!r}, got {status!r}"
        )

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_nct05770037_search_title_matches_direct_lookup(self):
        """Title returned by search_trials() must match title from direct lookup."""
        direct = _fetch_nct_direct(NCT05770037)
        expected_title = (
            direct["protocolSection"]["identificationModule"].get("briefTitle")
            or direct["protocolSection"]["identificationModule"].get("officialTitle")
        )

        async with ClinicalTrialsClient() as client:
            trials = await client.search_trials(
                gene="ALK",
                tumor_type="Neuroblastoma",
                recruiting_only=False,
                max_results=50,
            )

        nct_ids = [t.nct_id for t in trials]
        if NCT05770037 not in nct_ids:
            pytest.skip(
                f"{NCT05770037} not returned in ALK Neuroblastoma search. "
                f"IDs returned: {nct_ids}"
            )

        trial = next(t for t in trials if t.nct_id == NCT05770037)
        assert trial.title == expected_title, (
            f"{NCT05770037} title from search ({trial.title!r}) does not match "
            f"direct lookup ({expected_title!r}) — ID/title misalignment in parsing"
        )

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_nct05770037_status_matches_direct_lookup(self):
        """Status returned by search_trials() must match status from direct lookup."""
        direct = _fetch_nct_direct(NCT05770037)
        expected_status = direct["protocolSection"]["statusModule"].get("overallStatus")

        async with ClinicalTrialsClient() as client:
            trials = await client.search_trials(
                gene="ALK",
                tumor_type="Neuroblastoma",
                recruiting_only=False,
                max_results=50,
            )

        nct_ids = [t.nct_id for t in trials]
        if NCT05770037 not in nct_ids:
            pytest.skip(
                f"{NCT05770037} not returned in ALK Neuroblastoma search. "
                f"IDs returned: {nct_ids}"
            )

        trial = next(t for t in trials if t.nct_id == NCT05770037)
        assert trial.status == expected_status, (
            f"{NCT05770037} status from search ({trial.status!r}) does not match "
            f"direct lookup ({expected_status!r})"
        )


class TestNCTAlignmentALKNeuroblastoma:
    """Verify that all trials returned for ALK Neuroblastoma have consistent NCT/title/phase."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_each_nct_id_unique(self):
        """Search results must not have duplicate NCT IDs."""
        async with ClinicalTrialsClient() as client:
            trials = await client.search_all_trial_evidence(
                gene="ALK",
                variant="F1174L",
                tumor_type="Neuroblastoma",
                recruiting_only=False,
                max_results=20,
            )

        nct_ids = [t.nct_id for t in trials]
        assert len(nct_ids) == len(set(nct_ids)), (
            f"Duplicate NCT IDs found: "
            f"{[x for x in nct_ids if nct_ids.count(x) > 1]}"
        )

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_all_nct_ids_have_titles(self):
        """Every returned trial must have a non-empty title."""
        async with ClinicalTrialsClient() as client:
            trials = await client.search_all_trial_evidence(
                gene="ALK",
                variant="F1174L",
                tumor_type="Neuroblastoma",
                recruiting_only=False,
                max_results=20,
            )

        for trial in trials:
            assert trial.title and trial.title.strip(), (
                f"Trial {trial.nct_id} has an empty title — "
                f"ID and title may have been misaligned during parsing"
            )

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_all_nct_ids_start_with_nct_prefix(self):
        """All NCT IDs must follow the standard NCT + 8-digit format."""
        import re

        async with ClinicalTrialsClient() as client:
            trials = await client.search_all_trial_evidence(
                gene="ALK",
                variant="F1174L",
                tumor_type="Neuroblastoma",
                recruiting_only=False,
                max_results=20,
            )

        nct_pattern = re.compile(r"^NCT\d{8}$")
        for trial in trials:
            assert nct_pattern.match(trial.nct_id), (
                f"Malformed NCT ID: {trial.nct_id!r}"
            )

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_phase_values_are_known(self):
        """Phase values must be from the known ClinicalTrials.gov enum set.

        The API returns values like 'PHASE1', 'PHASE2', 'PHASE3', 'PHASE4',
        'EARLY_PHASE1', or None. Any other value suggests a parsing error.
        """
        known_phases = {
            "PHASE1", "PHASE2", "PHASE3", "PHASE4", "EARLY_PHASE1", None
        }

        async with ClinicalTrialsClient() as client:
            trials = await client.search_all_trial_evidence(
                gene="ALK",
                variant="F1174L",
                tumor_type="Neuroblastoma",
                recruiting_only=False,
                max_results=20,
            )

        for trial in trials:
            assert trial.phase in known_phases, (
                f"Trial {trial.nct_id} ({trial.title!r}) has unexpected phase "
                f"value: {trial.phase!r}. Known values: {known_phases}"
            )

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_status_values_are_known(self):
        """Status values must be from the known ClinicalTrials.gov enum set."""
        known_statuses = {
            "RECRUITING",
            "ENROLLING_BY_INVITATION",
            "NOT_YET_RECRUITING",
            "ACTIVE_NOT_RECRUITING",
            "COMPLETED",
            "SUSPENDED",
            "TERMINATED",
            "WITHDRAWN",
            "UNKNOWN",
        }

        async with ClinicalTrialsClient() as client:
            trials = await client.search_all_trial_evidence(
                gene="ALK",
                variant="F1174L",
                tumor_type="Neuroblastoma",
                recruiting_only=False,
                max_results=20,
            )

        for trial in trials:
            assert trial.status in known_statuses, (
                f"Trial {trial.nct_id} ({trial.title!r}) has unexpected status "
                f"value: {trial.status!r}"
            )

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_alk_gene_appears_in_title_or_eligibility(self):
        """Every trial in an ALK search must mention ALK in its title or eligibility.

        The search_trials() method already enforces this with a regex filter.
        This test confirms the filter is working correctly for this specific query.
        """
        import re

        async with ClinicalTrialsClient() as client:
            trials = await client.search_trials(
                gene="ALK",
                tumor_type="Neuroblastoma",
                recruiting_only=False,
                max_results=20,
            )

        alk_pattern = re.compile(r"\bALK\b")
        for trial in trials:
            searchable_text = f"{trial.title} {trial.eligibility_criteria or ''}".upper()
            assert alk_pattern.search(searchable_text), (
                f"Trial {trial.nct_id} ({trial.title!r}) was returned for an ALK "
                f"search but does not mention ALK in title or eligibility criteria"
            )
