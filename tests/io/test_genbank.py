import pytest

from plasmidlab.core import Feature, FeatureSegment, MoleculeTopology, SequenceRecord
from plasmidlab.io.genbank import (
    LOCATION_OPERATOR_QUALIFIER,
    LOCATION_WARNINGS_QUALIFIER,
    GenBankLossyWarning,
    dumps_genbank,
    loads_genbank,
)


def test_genbank_round_trip_synthetic_circular_plasmid() -> None:
    record = _synthetic_circular_plasmid()

    imported = loads_genbank(dumps_genbank(record))[0]
    reimported = loads_genbank(dumps_genbank(imported))[0]

    assert reimported.sequence == record.sequence
    assert reimported.topology is MoleculeTopology.CIRCULAR
    assert len(reimported.features) == len(record.features)
    assert _feature_signature(reimported) == _feature_signature(record)
    assert _qualifier_values(reimported.features[1], "gene") == ("synCds",)
    assert _qualifier_values(reimported.features[4], "note") == ("synthetic antibiotic marker",)


def test_genbank_round_trip_linear_pcr_product() -> None:
    record = SequenceRecord(
        id="pcr_linear",
        name="pcr_linear",
        description="Linear PCR product",
        sequence="ATGCGTACGTAGCTAGCTAGCGTACGAT",
        topology=MoleculeTopology.LINEAR,
        features=(
            Feature(
                type="misc_feature",
                start=4,
                end=20,
                strand=1,
                name="amplicon_body",
                qualifiers={"label": "amplicon_body"},
            ),
        ),
    )

    imported = loads_genbank(dumps_genbank(record))[0]

    assert imported.sequence == record.sequence
    assert imported.topology is MoleculeTopology.LINEAR
    assert _feature_signature(imported) == _feature_signature(record)


def test_genbank_round_trip_reverse_strand_feature() -> None:
    record = SequenceRecord(
        id="reverse_feature",
        sequence="ATGCGTACGTAGCTAGCTAGCGTACGAT",
        features=(
            Feature(
                type="CDS",
                start=8,
                end=24,
                strand=-1,
                name="rev_gene",
                qualifiers={"gene": "rev_gene", "product": "reverse strand protein"},
            ),
        ),
    )

    imported = loads_genbank(dumps_genbank(record))[0]

    assert imported.features[0].segments == (FeatureSegment(8, 24, strand=-1),)
    assert _qualifier_values(imported.features[0], "gene") == ("rev_gene",)
    assert _qualifier_values(imported.features[0], "product") == ("reverse strand protein",)


def test_genbank_round_trip_compound_cds_feature() -> None:
    record = SequenceRecord(
        id="split_cds",
        sequence="ATGCGTACGTAGCTAGCTAGCGTACGATATGCGTACGTAGCTAGCTAGCGTACGAT",
        features=(
            Feature(
                type="CDS",
                segments=(
                    FeatureSegment(5, 15, strand=1),
                    FeatureSegment(30, 45, strand=1),
                ),
                name="split_gene",
                qualifiers={"gene": "split_gene", "product": "split protein"},
            ),
        ),
    )

    imported = loads_genbank(dumps_genbank(record))[0]

    assert imported.features[0].type == "CDS"
    assert imported.features[0].segments == (
        FeatureSegment(5, 15, strand=1),
        FeatureSegment(30, 45, strand=1),
    )
    assert _qualifier_values(imported.features[0], "gene") == ("split_gene",)
    assert _qualifier_values(imported.features[0], LOCATION_OPERATOR_QUALIFIER) == ("join",)


def test_genbank_round_trip_reverse_strand_compound_cds_feature() -> None:
    record = SequenceRecord(
        id="rev_split_cds",
        sequence="ACGTACGTACGT",
        features=(
            Feature(
                type="CDS",
                segments=(
                    FeatureSegment(7, 10, strand=-1),
                    FeatureSegment(0, 3, strand=-1),
                ),
                name="rev_split",
                qualifiers={"gene": "rev_split", "translation": "MY"},
            ),
        ),
    )

    text = dumps_genbank(record)
    imported = loads_genbank(text)[0]
    reimported = loads_genbank(dumps_genbank(imported))[0]

    assert "complement(join(1..3,8..10))" in text.replace("\n", "")
    assert reimported.features[0].segments == (
        FeatureSegment(7, 10, strand=-1),
        FeatureSegment(0, 3, strand=-1),
    )
    assert _qualifier_values(reimported.features[0], LOCATION_OPERATOR_QUALIFIER) == ("join",)
    assert _qualifier_values(reimported.features[0], "translation") == ("MY",)


def test_genbank_multivalued_qualifiers_are_preserved() -> None:
    record = SequenceRecord(
        id="multi_qualifier",
        sequence="ATGCGTACGTAG",
        features=(
            Feature(
                type="misc_feature",
                start=0,
                end=6,
                qualifiers={
                    "note": ("first note", "second note"),
                    "db_xref": ("taxon:32630", "SO:0000001"),
                },
            ),
        ),
    )

    imported = loads_genbank(dumps_genbank(record))[0]
    reimported = loads_genbank(dumps_genbank(imported))[0]

    assert _qualifier_values(reimported.features[0], "note") == ("first note", "second note")
    assert _qualifier_values(reimported.features[0], "db_xref") == ("taxon:32630", "SO:0000001")


def test_genbank_translation_qualifier_is_preserved() -> None:
    translation = "MSTNPKPQRKTKRNTNRRPQDVKFPGG"
    record = SequenceRecord(
        id="translation_qualifier",
        sequence="ATG" * 12,
        features=(
            Feature(
                type="CDS",
                start=0,
                end=36,
                strand=1,
                qualifiers={"gene": "synthetic_cds", "translation": translation},
            ),
        ),
    )

    imported = loads_genbank(dumps_genbank(record))[0]

    assert _qualifier_values(imported.features[0], "translation") == (translation,)


def test_genbank_fuzzy_location_import_warns_and_keeps_numeric_boundaries() -> None:
    record = SequenceRecord(
        id="fuzzy_feature",
        sequence="ACGTACGTACGT",
        features=(Feature(type="misc_feature", start=2, end=8, qualifiers={"note": "fuzzy"}),),
    )
    fuzzy_text = dumps_genbank(record).replace("3..8", "<3..>8")

    with pytest.warns(GenBankLossyWarning, match="fuzzy location"):
        imported = loads_genbank(fuzzy_text)[0]

    feature = imported.features[0]
    assert feature.segments == (FeatureSegment(2, 8, strand=1),)
    assert "numeric boundaries" in _qualifier_values(feature, LOCATION_WARNINGS_QUALIFIER)[0]

    with pytest.warns(GenBankLossyWarning, match="fuzzy GenBank location"):
        exported = dumps_genbank(imported)
    assert "<3..>8" not in exported
    assert "3..8" in exported


def test_genbank_order_compound_operator_round_trip() -> None:
    record = SequenceRecord(
        id="ordered_feature",
        sequence="ACGTACGTACGT",
        features=(
            Feature(
                type="misc_feature",
                segments=(FeatureSegment(0, 3, strand=1), FeatureSegment(7, 10, strand=1)),
                qualifiers={LOCATION_OPERATOR_QUALIFIER: "order", "note": "ordered pieces"},
            ),
        ),
    )

    text = dumps_genbank(record)
    imported = loads_genbank(text)[0]

    assert "order(1..3,8..10)" in text.replace("\n", "")
    assert LOCATION_OPERATOR_QUALIFIER not in text
    assert imported.features[0].segments == (
        FeatureSegment(0, 3, strand=1),
        FeatureSegment(7, 10, strand=1),
    )
    assert _qualifier_values(imported.features[0], LOCATION_OPERATOR_QUALIFIER) == ("order",)


def test_genbank_source_feature_and_circular_topology_round_trip() -> None:
    record = SequenceRecord(
        id="source_topology",
        sequence="ACGTACGTACGT",
        topology=MoleculeTopology.CIRCULAR,
        features=(
            Feature(
                type="source",
                start=0,
                end=12,
                qualifiers={"organism": "synthetic construct", "mol_type": "other DNA"},
            ),
        ),
    )

    imported = loads_genbank(dumps_genbank(record))[0]

    assert imported.topology is MoleculeTopology.CIRCULAR
    assert imported.features[0].type == "source"
    assert imported.features[0].segments == (FeatureSegment(0, 12, strand=1),)
    assert _qualifier_values(imported.features[0], "organism") == ("synthetic construct",)
    assert _qualifier_values(imported.features[0], "mol_type") == ("other DNA",)


def test_genbank_remote_location_fails_clearly() -> None:
    record = SequenceRecord(
        id="remote_feature",
        sequence="ACGTACGTACGT",
        features=(Feature(type="misc_feature", start=0, end=3, qualifiers={"note": "remote"}),),
    )
    remote_text = dumps_genbank(record).replace("1..3", "J00194.1:1..3")

    with pytest.raises(ValueError, match="remote GenBank locations"):
        loads_genbank(remote_text)


def _synthetic_circular_plasmid() -> SequenceRecord:
    return SequenceRecord(
        id="pLAB1",
        name="pLAB1",
        description="Synthetic circular PlasmidLab test plasmid",
        sequence=("ATGCGTAC" * 15),
        topology=MoleculeTopology.CIRCULAR,
        features=(
            Feature(
                type="promoter",
                start=4,
                end=18,
                strand=1,
                name="P_syn",
                qualifiers={"label": "P_syn", "note": "synthetic promoter"},
            ),
            Feature(
                type="CDS",
                start=25,
                end=61,
                strand=1,
                name="synCds",
                qualifiers={"gene": "synCds", "product": "synthetic protein"},
            ),
            Feature(
                type="terminator",
                start=62,
                end=78,
                strand=1,
                name="T_syn",
                qualifiers={"label": "T_syn"},
            ),
            Feature(
                type="rep_origin",
                segments=(
                    FeatureSegment(112, 120, strand=1),
                    FeatureSegment(0, 6, strand=1),
                ),
                name="ori_syn",
                qualifiers={"label": "ori_syn"},
            ),
            Feature(
                type="misc_feature",
                start=82,
                end=108,
                strand=1,
                name="KanR_marker",
                qualifiers={"label": "KanR_marker", "note": "synthetic antibiotic marker"},
            ),
        ),
    )


def _feature_signature(record: SequenceRecord) -> tuple[tuple[str, tuple[tuple[int, int, int], ...]], ...]:
    return tuple(
        (
            feature.type,
            tuple((segment.start, segment.end, segment.strand) for segment in feature.segments),
        )
        for feature in record.features
    )


def _qualifier_values(feature: Feature, key: str) -> tuple[str, ...]:
    value = feature.qualifiers[key]
    if isinstance(value, str):
        return (value,)
    if isinstance(value, tuple | list):
        return tuple(str(item) for item in value)
    return (str(value),)
