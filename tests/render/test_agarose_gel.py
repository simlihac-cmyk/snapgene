from plasmidlab.core import GelLaneInput, simulate_gel
from plasmidlab.render import gel_to_svg, write_gel_svg


def test_gel_svg_contains_lanes_bands_and_run_parameters(tmp_path) -> None:
    model = simulate_gel(
        (GelLaneInput(name="digest", fragments=(3000, 1000, 500)),),
        ladder="100_bp_dna_ladder",
        agarose_percentage=1.2,
        run_time_minutes=50,
        voltage=90,
    )
    path = tmp_path / "gel.svg"

    svg = gel_to_svg(model)
    write_gel_svg(model, path)

    assert svg.startswith('<?xml version="1.0"')
    assert "Agarose gel simulation" in svg
    assert "digest" in svg
    assert "3 kb" in svg
    assert "1.2% agarose" in svg
    assert path.read_text(encoding="utf-8") == svg
