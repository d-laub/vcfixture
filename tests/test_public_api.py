import vcfixture


def test_public_exports():
    for name in [
        "VcfBuilder",
        "Number",
        "Type",
        "Genotype",
        "GroundTruth",
        "Reference",
        "ReferenceBuilder",
        "ReferenceSpec",
        "RepeatFeature",
        "strategies",
    ]:
        assert hasattr(vcfixture, name), f"missing public export: {name}"
