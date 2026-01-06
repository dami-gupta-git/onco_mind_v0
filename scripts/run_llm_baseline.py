#!/usr/bin/env python3
"""Test LLM baseline knowledge for a variant (no evidence provided).

Usage:
    python scripts/test_llm_baseline.py PIK3CA H1047R "Thyroid Cancer"
    python scripts/test_llm_baseline.py CHEK2 I157T "Breast Cancer"
    python scripts/test_llm_baseline.py BRAF V600E
"""

import asyncio
import sys

from oncomind.llm.service import LLMService


async def main():
    if len(sys.argv) < 3:
        print("Usage: python scripts/test_llm_baseline.py GENE VARIANT [TUMOR_TYPE]")
        print("Example: python scripts/test_llm_baseline.py PIK3CA H1047R 'Thyroid Cancer'")
        sys.exit(1)

    gene = sys.argv[1]
    variant = sys.argv[2]
    tumor_type = sys.argv[3] if len(sys.argv) > 3 else None

    # Use the default model (Claude)
    service = LLMService()

    result = await service.test_llm_baseline_knowledge(gene, variant, tumor_type)

    print("\n" + "=" * 80)
    print("PARSED RESULT:")
    print("=" * 80)
    import json
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    asyncio.run(main())
