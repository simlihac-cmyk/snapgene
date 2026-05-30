from importlib import resources

from plasmidlab.core import (
    FeatureLibraryEntry,
    FeatureStrandBehavior,
    MoleculeTopology,
    SequenceRecord,
    apply_feature_annotations,
    detect_features,
    load_feature_library,
)
from plasmidlab.core.feature_annotation import DEFAULT_FEATURE_LIBRARY_PACKAGE


def test_exact_feature_detection_from_open_default_library() -> None:
    library = load_feature_library()
    record = SequenceRecord(id="t7", sequence="CCCC" + "TAATACGACTCACTATAGGG" + "AAAA")

    detections = detect_features(record, library)

    t7 = next(detection for detection in detections if detection.name == "T7 promoter fragment")
    assert t7.start == 4
    assert t7.end == 24
    assert t7.strand == 1
    assert t7.identity == 1.0
    assert t7.exact


def test_default_feature_library_loads_from_package_resources() -> None:
    library = load_feature_library()
    names = {entry.name for entry in library}
    resource = resources.files(DEFAULT_FEATURE_LIBRARY_PACKAGE).joinpath(
        "common_synthetic_features.json"
    )

    assert resource.is_file()
    assert "PlasmidLab common synthetic feature fragments" in resource.read_text(encoding="utf-8")
    assert {
        "lac promoter fragment",
        "T7 promoter fragment",
        "generic ori fragment",
        "generic antibiotic marker fragment",
    } <= names
    assert all(entry.source_path is not None for entry in library)
    assert all(str(entry.source_path).startswith("package:plasmidlab.data.features/") for entry in library)


def test_reverse_complement_detection() -> None:
    entry = FeatureLibraryEntry(
        name="synthetic feature",
        type="misc_feature",
        sequence="ATGAAACCC",
        strand_behavior=FeatureStrandBehavior.BOTH,
    )
    record = SequenceRecord(id="reverse", sequence="GGG" + "GGGTTTCAT" + "TTT")

    detections = detect_features(record, (entry,))

    assert len(detections) == 1
    assert detections[0].start == 3
    assert detections[0].end == 12
    assert detections[0].strand == -1
    assert detections[0].matched_sequence == "GGGTTTCAT"


def test_near_exact_feature_detection() -> None:
    entry = FeatureLibraryEntry(
        name="near feature",
        type="misc_feature",
        sequence="ATGAAACCC",
        minimum_identity=8 / 9,
    )
    record = SequenceRecord(id="near", sequence="GGGATGAAAGCCTTT")

    detections = detect_features(record, (entry,))

    assert len(detections) == 1
    assert detections[0].mismatches == 1
    assert detections[0].identity == 8 / 9


def test_circular_origin_crossing_feature_detection() -> None:
    entry = FeatureLibraryEntry(
        name="origin crossing",
        type="misc_feature",
        sequence="AAACCCGGG",
        strand_behavior=FeatureStrandBehavior.FORWARD_ONLY,
    )
    record = SequenceRecord(
        id="circular",
        sequence="CCCGGGTTTTAAA",
        topology=MoleculeTopology.CIRCULAR,
    )

    detections = detect_features(record, (entry,))

    assert len(detections) == 1
    detection = detections[0]
    assert detection.wraps_origin
    assert detection.start == 10
    assert detection.end == 6
    assert [(segment.start, segment.end) for segment in detection.segments] == [(10, 13), (0, 6)]


def test_duplicate_feature_resolution_keeps_best_site_proposal() -> None:
    entries = (
        FeatureLibraryEntry(name="duplicate", type="misc_feature", sequence="ATGAAACCC"),
        FeatureLibraryEntry(
            name="duplicate loose",
            type="misc_feature",
            sequence="ATGAAACCC",
            minimum_identity=8 / 9,
        ),
    )
    record = SequenceRecord(id="dup", sequence="GGGATGAAACCCTTT")

    detections = detect_features(record, entries)

    assert len(detections) == 1
    assert detections[0].name == "duplicate"
    assert detections[0].identity == 1.0


def test_apply_feature_annotations_adds_features_and_provenance() -> None:
    entry = FeatureLibraryEntry(name="apply me", type="promoter", sequence="ATGAAACCC")
    record = SequenceRecord(id="apply", sequence="GGGATGAAACCCTTT")
    detections = detect_features(record, (entry,))

    updated = apply_feature_annotations(record, detections)
    repeated = apply_feature_annotations(updated, detections)

    assert len(updated.features) == 1
    assert updated.features[0].name == "apply me"
    assert updated.features[0].qualifiers["detected_identity"] == "1.000"
    assert updated.history[-1].operation == "detect_features"
    assert updated.history[-1].parameters["applied_count"] == 1
    assert repeated is updated


def test_user_feature_library_json_can_be_loaded(tmp_path) -> None:
    library_path = tmp_path / "lab_features.json"
    library_path.write_text(
        """
        {
          "features": [
            {
              "name": "lab tag",
              "type": "misc_feature",
              "sequence": "ACACACAC",
              "aliases": ["LT"],
              "strand_behavior": "forward_only",
              "minimum_identity": 1.0,
              "notes": "Local lab feature"
            }
          ]
        }
        """,
        encoding="utf-8",
    )

    library = load_feature_library(library_path)

    assert len(library) == 1
    assert library[0].name == "lab tag"
    assert library[0].aliases == ("LT",)
    assert library[0].source_path == str(library_path)
