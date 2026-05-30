from plasmidlab.core import (
    BandAmountMode,
    GelInputFragment,
    GelLaneInput,
    MoleculeTopology,
    SequenceRecord,
    builtin_ladders,
    digest,
    load_ladder,
    simulate_gel,
    simulate_pcr,
)


def test_fragment_sizes_order_by_migration() -> None:
    model = simulate_gel((GelLaneInput(name="digest", fragments=(10000, 3000, 500)),))
    bands = model.lanes[0].bands

    assert [band.size_bp for band in bands] == [10000, 3000, 500]
    assert bands[0].migration < bands[1].migration < bands[2].migration


def test_larger_fragments_migrate_less_than_smaller_fragments() -> None:
    model = simulate_gel((GelLaneInput(name="sizes", fragments=(5000, 1000)),))
    large, small = model.lanes[0].bands

    assert large.size_bp == 5000
    assert small.size_bp == 1000
    assert large.migration < small.migration


def test_equal_mass_and_equal_moles_intensity_are_explicit() -> None:
    equal_mass = simulate_gel(
        (GelLaneInput(name="mass", fragments=(1000, 500)),),
        amount_mode=BandAmountMode.EQUAL_MASS,
    )
    equal_moles = simulate_gel(
        (GelLaneInput(name="moles", fragments=(1000, 500)),),
        amount_mode=BandAmountMode.EQUAL_MOLES,
    )

    mass_bands = equal_mass.lanes[0].bands
    mole_bands = equal_moles.lanes[0].bands
    assert mass_bands[0].relative_intensity == mass_bands[1].relative_intensity
    assert mole_bands[0].mass_ng > mole_bands[1].mass_ng
    assert mole_bands[0].relative_intensity > mole_bands[1].relative_intensity


def test_explicit_fragment_mass_overrides_loading_mode() -> None:
    model = simulate_gel(
        (
            GelLaneInput(
                name="weighted",
                fragments=(
                    GelInputFragment(size_bp=1000, mass_ng=10),
                    GelInputFragment(size_bp=500, mass_ng=50),
                ),
            ),
        ),
        amount_mode=BandAmountMode.EQUAL_MOLES,
    )

    large, small = model.lanes[0].bands
    assert large.mass_ng == 10
    assert small.mass_ng == 50
    assert small.relative_intensity == 1.0


def test_builtin_ladders_load_from_open_json_definitions() -> None:
    ladders = builtin_ladders()
    ladder = load_ladder("1 kb DNA ladder")
    model = simulate_gel((GelLaneInput(name="sample", fragments=(750,)),), ladder=ladder)

    assert "1kb_dna_ladder" in {definition.id for definition in ladders}
    assert model.lanes[0].is_ladder
    assert model.lanes[0].name == "1 kb DNA ladder"
    assert model.lanes[0].bands[0].size_bp == 10000
    assert model.lanes[1].bands[0].size_bp == 750


def test_digest_and_pcr_products_can_be_loaded_as_fragments() -> None:
    digest_record = SequenceRecord(
        id="digest_template",
        sequence="A" * 30 + "GAATTC" + "C" * 40 + "GGATCC" + "T" * 50,
        topology=MoleculeTopology.LINEAR,
    )
    pcr_record = SequenceRecord(id="pcr_template", sequence="AAAAACCCCCGGGGGTTTTT")

    digest_fragments = digest(digest_record, "EcoRI,BamHI")
    pcr_product = simulate_pcr(pcr_record, "AAAAA", "AAAAA")
    model = simulate_gel(
        (
            ("Digest", digest_fragments),
            ("PCR", (pcr_product,)),
        )
    )

    assert [band.size_bp for band in model.lanes[0].bands] == sorted(
        (fragment.length for fragment in digest_fragments),
        reverse=True,
    )
    assert model.lanes[1].bands[0].size_bp == pcr_product.length
