# TODO

## Validation

- [ ] Build systematic validation pipeline with domain expert review
- [ ] Automate testing across representative variant sets (beyond the current 15 manual checks)
- [ ] Improve negation detection in FDA label parsing ("not demonstrated", "not approved")
- [ ] Address edge cases for rare variant-disease pairings with inconsistent evidence grading
- [ ] Validate LLM synthesis quality for understudied variants (where hallucination risk is highest)

## Structural Variant Support

- [ ] Add fusion variant parsing and normalization
- [ ] Add amplification and copy-number variant support
- [ ] Extend evidence aggregation to handle structural variants across all sources
- [ ] Update LLM prompts for structural variant context

## Code Quality

- [ ] Remove commented-out Evidence Specificity panel in `cli.py` (lines 278-309) -- pending validation
- [ ] Remove commented-out CBIOPORTAL_STUDY_MAPPINGS duplicate in `constants.py` (lines 974-1105)
- [ ] Replace remaining `print()` calls in `llm/service.py` (lines 257, 453, 579) with `logger.debug()`
- [ ] Replace `print()` calls in `cli.py` batch command (lines 542, 571, 583, etc.) with proper logging
- [ ] Audit `service.py` for the `test_llm_baseline_knowledge` method -- determine if it should remain or be test-only

## Output Formats

- [ ] Add `--format csv` output for CLI
- [ ] Add `--format markdown` output for CLI
- [ ] Add `--split` option for per-variant output in batch mode
