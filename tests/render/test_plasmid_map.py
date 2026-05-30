from pathlib import Path

from plasmidlab.core import (
    EnzymeSite,
    Feature,
    FeatureSegment,
    MoleculeTopology,
    Primer,
    PrimerBinding,
    SequenceRecord,
    analyze_restriction_sites,
)
from plasmidlab.render import MapStyle, prepare_map_overlays, render_plasmid_map, to_svg, write_svg


def test_circular_map_geometry_contains_expected_primitives() -> None:
    record = _circular_record()

    model = render_plasmid_map(record, style=MapStyle(width=600, height=500, circular_radius=120))

    assert model.map_kind == "circular"
    assert model.center is not None
    assert model.center.x == 300
    assert model.center.y == 250
    assert model.origin_marker.angle == -90
    assert len(model.feature_arcs) == 3
    assert len(model.feature_arrows) == 3
    assert len(model.enzyme_ticks) == 2
    assert len(model.primer_arrows) == 1
    assert model.feature_arcs[0].start_angle == -90
    assert model.feature_arcs[0].end_angle == -18
    assert model.feature_arcs[0].color == "#2f9e44"
    assert model.feature_arcs[1].lane != model.feature_arcs[2].lane


def test_linear_map_geometry_uses_scaled_coordinates() -> None:
    record = SequenceRecord(
        id="linear_map",
        sequence="A" * 100,
        topology=MoleculeTopology.LINEAR,
        features=(Feature(type="CDS", start=25, end=75, strand=-1, name="gene"),),
        enzyme_sites=(EnzymeSite("EcoRI", "GAATTC", 50, 56, cut_index=1),),
    )
    style = MapStyle(width=500, height=300, margin=50)

    model = render_plasmid_map(record, style=style)

    assert model.map_kind == "linear"
    assert model.backbone_start.x == 50
    assert model.backbone_end.x == 450
    assert model.feature_arcs[0].start_point.x == 150
    assert model.feature_arcs[0].end_point.x == 350
    assert model.enzyme_ticks[0].start_point.x == 250
    assert model.feature_arrows[0].strand == -1


def test_svg_export_smoke_for_synthetic_circular_plasmid(tmp_path: Path) -> None:
    model = render_plasmid_map(_circular_record())
    svg = to_svg(model)
    path = tmp_path / "map.svg"

    write_svg(model, path)

    assert svg.startswith('<?xml version="1.0"')
    assert "<svg" in svg
    assert "<path" in svg
    assert "P_test" in svg
    assert "EcoRI" in svg
    assert "Fwd" in svg
    assert path.read_text(encoding="utf-8") == svg


def test_feature_colors_are_configurable() -> None:
    record = _circular_record()

    model = render_plasmid_map(record, feature_colors={"promoter": "#123456"})

    assert model.feature_arcs[0].color == "#123456"


def test_svg_export_handles_full_circle_feature() -> None:
    record = SequenceRecord(
        id="full_feature",
        sequence="A" * 40,
        topology=MoleculeTopology.CIRCULAR,
        features=(Feature(type="misc_feature", start=0, end=40, name="whole"),),
    )

    svg = to_svg(render_plasmid_map(record))

    assert "whole" in svg
    assert svg.count("<circle") >= 2


def test_map_overlay_keeps_circular_primer_binding_across_origin() -> None:
    record = SequenceRecord(
        id="primer_wrap",
        sequence="A" * 100,
        topology=MoleculeTopology.CIRCULAR,
        primers=(Primer(name="wrap_primer", sequence="A" * 10),),
    )
    overlays = prepare_map_overlays(
        record,
        primer_bindings_by_primer_index={
            0: (
                PrimerBinding(
                    primer_sequence="A" * 10,
                    binding_sequence="A" * 10,
                    tail_sequence="",
                    start=95,
                    end=5,
                    strand=1,
                    wraps_origin=True,
                ),
            ),
        },
    )

    model = render_plasmid_map(record, overlays=overlays)

    assert len(model.primer_arrows) == 1
    assert model.primer_arrows[0].start == 95
    assert model.primer_arrows[0].end == 5
    assert model.primer_arrows[0].wraps_origin


def test_map_overlay_keeps_circular_enzyme_site_across_origin() -> None:
    record = SequenceRecord(
        id="enzyme_wrap",
        sequence="ATTC" + "C" * 14 + "GA",
        topology=MoleculeTopology.CIRCULAR,
    )
    analysis = analyze_restriction_sites(record, "EcoRI")
    overlays = prepare_map_overlays(record, restriction_analysis=analysis)

    model = render_plasmid_map(record, overlays=overlays)

    assert len(model.enzyme_ticks) == 1
    assert model.enzyme_ticks[0].enzyme_name == "EcoRI"
    assert model.enzyme_ticks[0].site_start == 18
    assert model.enzyme_ticks[0].site_end == 4
    assert model.enzyme_ticks[0].wraps_origin


def test_compound_feature_segments_are_rendered_without_flattening() -> None:
    record = SequenceRecord(
        id="compound",
        sequence="A" * 100,
        topology=MoleculeTopology.CIRCULAR,
        features=(
            Feature(
                type="CDS",
                name="split",
                segments=(FeatureSegment(80, 100, strand=1), FeatureSegment(0, 10, strand=1)),
            ),
        ),
    )

    model = render_plasmid_map(record)

    assert [(arc.start, arc.end) for arc in model.feature_arcs] == [(80, 100), (0, 10)]
    assert len(model.feature_arrows) == 2


def _circular_record() -> SequenceRecord:
    return SequenceRecord(
        id="pMap",
        sequence="A" * 100,
        topology=MoleculeTopology.CIRCULAR,
        features=(
            Feature(type="promoter", start=0, end=20, strand=1, name="P_test"),
            Feature(type="CDS", start=25, end=70, strand=1, name="geneA"),
            Feature(type="terminator", start=40, end=90, strand=-1, name="T_test"),
        ),
        enzyme_sites=(
            EnzymeSite("EcoRI", "GAATTC", 10, 16, cut_index=1),
            EnzymeSite("BamHI", "GGATCC", 75, 81, cut_index=1),
        ),
        primers=(Primer(name="Fwd", sequence="AAAAAA", start=5, end=11, strand=1),),
    )
