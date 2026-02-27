

from civicpy import civic

# Search for the variant (returns a list; usually 0 or 1 match for precise names)
variant = civic.get_variant_by_id(931)

if variant:
    # variant = variants[0]  # Take the first (best) match
    print(f"Found variant: gene={variant.gene.name} variant={variant.name} ID={variant.id}")

    simple_mp = variant.single_variant_molecular_profile

    if simple_mp:
        print(f"Simple Molecular Profile:")
        print(f"  - ID: {simple_mp.id}")
        print(f"  - Name: {simple_mp.name}")  # e.g., "BRAF V600E"
        print(f"  - CIViC Score: {simple_mp.molecular_profile_score}")
        print(f"  - Description/Structure: {simple_mp.description or 'N/A'}")

        # Now access evidence attached to this simple MP
        evidence_items = simple_mp.evidence_items  # list of Evidence objects

        print(f"\nAssociated Evidence Items: {len(evidence_items)}")
        for ev in evidence_items[:5]:  # show first few
            print(f"  - EID {ev.id=}: {ev.evidence_type=} | Level {ev.evidence_level=}")
            print(f"  - {ev.evidence_direction=}, {ev.evidence_type=}, {ev.significance=}")
            print(f"    {ev.description[:100]}...")
            # print(f"    Significance: {ev.clinical_significance}")
            print(f"    Therapies: {', '.join(t.name for t in ev.therapies) if ev.therapies else 'None'}")
            print("-" * 60)
    else:
        print("No simple molecular profile found for this variant (unusual).")