import pytest
from hypothesis import given
from hypothesis import strategies as st

from plasmidlab.core import (
    EnzymeSite,
    Feature,
    FeatureSegment,
    MoleculeTopology,
    Primer,
    SequenceRecord,
    Topology,
)


def test_sequence_record_uses_zero_based_half_open_features() -> None:
    feature = Feature(start=0, end=3, type="CDS", label="tiny protein")

    record = SequenceRecord(
        id="synthetic-1",
        sequence="atgc",
        topology=Topology.CIRCULAR,
        features=(feature,),
    )

    assert record.sequence == "ATGC"
    assert record.length == 4
    assert record.is_circular
    assert record.features[0].start == 0
    assert record.features[0].end == 3


def test_feature_rejects_empty_interval() -> None:
    with pytest.raises(ValueError, match="greater than start"):
        Feature(start=3, end=3, type="misc_feature")


def test_record_rejects_feature_outside_sequence() -> None:
    feature = Feature(start=1, end=5, type="promoter")

    with pytest.raises(ValueError, match="exceeds sequence length"):
        SequenceRecord(id="bad", sequence="ATGC", features=(feature,))


def test_compound_feature_can_wrap_on_circular_record() -> None:
    feature = Feature(
        type="CDS",
        name="origin spanning gene",
        segments=(
            FeatureSegment(8, 10, strand=1),
            FeatureSegment(0, 2, strand=1),
        ),
    )

    record = SequenceRecord(
        id="circular",
        sequence="AACCGGTTAA",
        topology=MoleculeTopology.CIRCULAR,
        features=(feature,),
    )

    assert record.features[0].segments[0].start == 8
    assert record.features[0].segments[1].end == 2


def test_compound_feature_wrap_is_rejected_on_linear_record() -> None:
    feature = Feature(
        type="CDS",
        segments=(
            FeatureSegment(8, 10, strand=1),
            FeatureSegment(0, 2, strand=1),
        ),
    )

    with pytest.raises(ValueError, match="linear compound features"):
        SequenceRecord(id="linear", sequence="AACCGGTTAA", features=(feature,))


def test_sequence_and_feature_name_validation() -> None:
    with pytest.raises(ValueError, match="invalid DNA sequence"):
        SequenceRecord(id="bad-sequence", sequence="ATX")

    with pytest.raises(ValueError, match="feature name"):
        Feature(start=0, end=1, type="misc_feature", name=" bad")


def test_primer_and_enzyme_site_models_validate_against_record_length() -> None:
    primer = Primer(name="forward", sequence="gaat", start=0, end=4, strand=1)
    enzyme_site = EnzymeSite(
        enzyme_name="EcoRI",
        recognition_sequence="GAATTC",
        start=0,
        end=6,
        cut_index=1,
    )

    record = SequenceRecord(
        id="with-sites",
        sequence="GAATTC",
        primers=(primer,),
        enzyme_sites=(enzyme_site,),
    )

    assert record.primers[0].sequence == "GAAT"
    assert record.enzyme_sites[0].recognition_sequence == "GAATTC"


def test_insert_invalidates_primer_bindings_and_drops_enzyme_sites() -> None:
    record = _record_with_derived_annotations(id="insert-derived")

    inserted = record.insert(4, "GG")

    assert inserted.primers == (
        Primer(name="bound", sequence="CCGG", strand=0, target_id="insert-derived"),
        Primer(name="inventory", sequence="TTAA"),
    )
    assert inserted.enzyme_sites == ()
    summary = inserted.history[-1].parameters["annotation_semantics"]
    assert summary["primer_bindings_invalidated"] == 1
    assert summary["enzyme_sites_dropped"] == 1
    assert "primer binding(s) invalidated" in inserted.history[-1].description
    assert "enzyme site annotation(s) dropped" in inserted.history[-1].description


def test_delete_invalidates_primer_bindings_and_drops_enzyme_sites() -> None:
    record = _record_with_derived_annotations(id="delete-derived")

    deleted = record.delete(2, 5)

    assert deleted.primers[0].start is None
    assert deleted.primers[0].end is None
    assert deleted.primers[0].strand == 0
    assert deleted.enzyme_sites == ()
    summary = deleted.history[-1].parameters["annotation_semantics"]
    assert summary["primer_bindings_invalidated"] == 1
    assert summary["enzyme_sites_dropped"] == 1


def test_replace_invalidates_primer_bindings_and_drops_enzyme_sites() -> None:
    record = _record_with_derived_annotations(id="replace-derived")

    replaced = record.replace(3, 7, "A")

    assert replaced.primers[0].start is None
    assert replaced.primers[0].end is None
    assert replaced.enzyme_sites == ()
    summary = replaced.history[-1].parameters["annotation_semantics"]
    assert summary["primer_bindings_invalidated"] == 1
    assert summary["enzyme_sites_dropped"] == 1


def test_reverse_complement_remaps_primer_bindings_and_drops_enzyme_sites() -> None:
    record = _record_with_derived_annotations(id="rc-derived")

    reversed_record = record.reverse_complement()

    assert reversed_record.primers[0] == Primer(
        name="bound",
        sequence="CCGG",
        start=8,
        end=12,
        strand=-1,
        target_id="rc-derived",
    )
    assert reversed_record.primers[1] == Primer(name="inventory", sequence="TTAA")
    assert reversed_record.enzyme_sites == ()
    summary = reversed_record.history[-1].parameters["annotation_semantics"]
    assert summary["primer_bindings_remapped"] == 1
    assert summary["primer_bindings_invalidated"] == 0
    assert summary["enzyme_sites_dropped"] == 1


def test_circular_edit_crossing_origin_does_not_leave_stale_derived_annotations() -> None:
    record = _record_with_derived_annotations(
        id="circular-derived",
        topology=MoleculeTopology.CIRCULAR,
    )

    edited = record.delete(14, 2)

    assert edited.sequence == "CCGGTTAATTCC"
    assert edited.primers[0].start is None
    assert edited.primers[0].end is None
    assert edited.enzyme_sites == ()
    summary = edited.history[-1].parameters["annotation_semantics"]
    assert summary["primer_policy"] == "binding_coordinates_invalidated"
    assert summary["enzyme_site_policy"] == "dropped_recompute_from_selected_enzyme_set"


@given(st.text(alphabet="ACGTRYSWKMBDHVN", min_size=0, max_size=80))
def test_reverse_complement_twice_returns_original_sequence(sequence: str) -> None:
    record = SequenceRecord(id="property", sequence=sequence)

    restored = record.reverse_complement().reverse_complement()

    assert restored.sequence == record.sequence
    assert restored.topology == record.topology


def test_reverse_complement_maps_feature_coordinates_and_strand() -> None:
    feature = Feature(start=2, end=5, type="CDS", strand=1, name="gene")
    record = SequenceRecord(id="rc", sequence="AACCGGTT", features=(feature,))

    reversed_record = record.reverse_complement()

    assert reversed_record.sequence == "AACCGGTT"
    assert reversed_record.features[0].segments == (FeatureSegment(3, 6, strand=-1),)
    assert reversed_record.history[-1].operation == "reverse_complement"


def test_core_edits_create_new_output_version_ids() -> None:
    record = SequenceRecord(id="versioned", sequence="AACCGGTT")

    inserted = record.insert(2, "AA")
    deleted = inserted.delete(0, 1)
    replaced = deleted.replace(1, 3, "TTT")
    reversed_record = replaced.reverse_complement()

    records = (record, inserted, deleted, replaced, reversed_record)
    assert len({item.version_id for item in records}) == len(records)
    for previous, current in zip(records[:-1], records[1:], strict=True):
        event = current.history[-1]
        assert event.input_version_ids == (previous.version_id,)
        assert event.output_version_id == current.version_id
        assert event.output_record_id == current.id


def test_insert_shifts_and_expands_features() -> None:
    feature = Feature(start=2, end=4, type="CDS", strand=1)
    record = SequenceRecord(id="insert", sequence="AAAAAA", features=(feature,))

    inserted_before = record.insert(1, "GGG")
    inserted_inside = record.insert(3, "GG")

    assert inserted_before.sequence == "AGGGAAAAA"
    assert inserted_before.features[0].segments == (FeatureSegment(5, 7, strand=1),)
    assert inserted_before.history[-1].operation == "insert"
    assert inserted_inside.features[0].segments == (FeatureSegment(2, 6, strand=1),)


def test_delete_shifts_truncates_and_merges_features() -> None:
    spanning = Feature(start=2, end=8, type="CDS", strand=1, name="spanning")
    after = Feature(start=8, end=10, type="terminator", strand=1, name="after")
    record = SequenceRecord(id="delete", sequence="AAAAAAAAAA", features=(spanning, after))

    deleted = record.delete(4, 6)

    assert deleted.sequence == "AAAAAAAA"
    assert deleted.features[0].segments == (FeatureSegment(2, 6, strand=1),)
    assert deleted.features[1].segments == (FeatureSegment(6, 8, strand=1),)
    assert deleted.history[-1].operation == "delete"


def test_replace_preserves_overlapping_feature_and_adds_single_event() -> None:
    feature = Feature(start=1, end=5, type="CDS", strand=1)
    record = SequenceRecord(id="replace", sequence="AAAAAA", features=(feature,))

    replaced = record.replace(2, 4, "GGG")

    assert replaced.sequence == "AAGGGAA"
    assert replaced.features[0].segments == (FeatureSegment(1, 6, strand=1),)
    assert [event.operation for event in replaced.history] == ["replace"]


def test_circular_replace_not_crossing_origin_preserves_origin() -> None:
    record = SequenceRecord(
        id="circular-replace",
        sequence="AAAACCCC",
        topology=MoleculeTopology.CIRCULAR,
        features=(Feature(start=4, end=8, type="misc_feature", strand=1, name="right-half"),),
    )

    replaced = record.replace(2, 6, "GG")

    assert replaced.sequence == "AAGGCC"
    assert replaced.topology is MoleculeTopology.CIRCULAR
    assert replaced.sequence.startswith("AA")
    assert replaced.features[0].segments == (FeatureSegment(2, 6, strand=1),)
    assert replaced.history[-1].affected_ranges == ((2, 6),)


def test_circular_replace_crossing_origin_is_rejected_without_reordering_sequence() -> None:
    record = _record_with_derived_annotations(
        id="circular-replace-derived",
        topology=MoleculeTopology.CIRCULAR,
    )

    with pytest.raises(ValueError, match="origin-crossing circular replace is ambiguous"):
        record.replace(14, 2, "GG")

    assert record.sequence == "AACCGGTTAATTCCGG"
    assert record.primers[0].start == 4
    assert record.enzyme_sites[0].start == 8


def test_linear_replace_does_not_wrap() -> None:
    record = SequenceRecord(id="linear-replace", sequence="AAAACCCC")

    with pytest.raises(ValueError, match="linear range start"):
        record.replace(6, 2, "GG")


def test_circular_slice_crossing_origin_preserves_wrapping_feature() -> None:
    feature = Feature(
        type="misc_feature",
        name="origin feature",
        segments=(
            FeatureSegment(6, 8, strand=1),
            FeatureSegment(0, 2, strand=1),
        ),
    )
    record = SequenceRecord(
        id="origin-slice",
        sequence="AACCGGTT",
        topology=MoleculeTopology.CIRCULAR,
        features=(feature,),
    )

    sliced = record.slice(6, 2)

    assert sliced.sequence == "TTAA"
    assert sliced.topology is MoleculeTopology.LINEAR
    assert sliced.features[0].segments == (FeatureSegment(0, 4, strand=1),)
    assert sliced.history[-1].operation == "slice"


def test_circular_slice_same_start_end_is_zero_length() -> None:
    record = SequenceRecord(
        id="zero-slice",
        sequence="AACCGGTT",
        topology=MoleculeTopology.CIRCULAR,
        features=(Feature(start=2, end=6, type="misc_feature"),),
    )

    sliced = record.slice(3, 3)

    assert sliced.sequence == ""
    assert sliced.features == ()
    assert sliced.history[-1].affected_ranges == ()


def test_circularize_and_linearize_add_provenance() -> None:
    record = SequenceRecord(id="topology", sequence="AACCGGTT")

    circular = record.circularize()
    linear = circular.linearize(cut=4)

    assert circular.topology is MoleculeTopology.CIRCULAR
    assert circular.history[-1].operation == "circularize"
    assert linear.sequence == "GGTTAACC"
    assert linear.topology is MoleculeTopology.LINEAR
    assert linear.history[-1].operation == "linearize"


def _record_with_derived_annotations(
    *,
    id: str,
    topology: MoleculeTopology = MoleculeTopology.LINEAR,
) -> SequenceRecord:
    return SequenceRecord(
        id=id,
        sequence="AACCGGTTAATTCCGG",
        topology=topology,
        primers=(
            Primer(name="bound", sequence="CCGG", start=4, end=8, strand=1, target_id=id),
            Primer(name="inventory", sequence="TTAA"),
        ),
        enzyme_sites=(
            EnzymeSite(
                enzyme_name="EcoRI",
                recognition_sequence="GAATTC",
                start=8,
                end=14,
                strand=1,
                cut_index=1,
            ),
        ),
    )
