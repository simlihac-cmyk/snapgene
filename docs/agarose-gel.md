# Agarose Gel Simulation

PlasmidLab includes an approximate agarose gel visualization model for digest and PCR
fragments. The model is intended for planning, reports, and quick sanity checks. It is
not an empirical, instrument-grade migration predictor, and it should not be expected
to exactly match any commercial software, gel box, ladder lot, stain, buffer, or imaging
workflow.

## Inputs

`plasmidlab.core.gel.simulate_gel()` accepts one or more lanes. Lane fragments may be:

- integer sizes in base pairs
- `DigestFragment` values from restriction digest simulation
- `PCRProduct` values from PCR simulation
- explicit `GelInputFragment` values with optional `mass_ng` or `moles_fmol`

Display parameters include agarose percentage, run time, voltage, and amount mode.
These parameters alter the approximate visual spacing but are not calibrated to a
specific instrument.

## Ladders

Built-in ladder definitions are open JSON files under
`src/plasmidlab/resources/ladders/`. They are synthetic reference definitions for
visualization, not proprietary product data.

Available bundled definitions include:

- `1kb_dna_ladder`
- `100_bp_dna_ladder`

Custom ladder JSON can be loaded with `load_ladder(path)` or passed directly as a
mapping.

## Migration Model

Band migration is modeled as a log-size relationship: larger fragments remain closer
to the wells, and smaller fragments migrate farther down the gel. Agarose percentage,
run time, and voltage scale the display in a deterministic, approximate way.

Internal migration values are relative distances from the well:

- `0.0` means near the well
- `1.0` means near the bottom of the gel

## Intensity Model

Band intensity is based on estimated DNA mass.

In `equal_mass` mode, fragments without explicit amounts are assigned the same mass,
so equal-size-independent loading produces equal relative intensities.

In `equal_moles` mode, fragments without explicit amounts are assigned the same molar
amount, so larger fragments have more DNA mass and appear brighter.

Explicit `mass_ng` takes precedence over the selected mode. Explicit `moles_fmol` is
converted using an approximate dsDNA mass of 660 g/mol per base pair.

## Exports

The GUI-independent renderer exports SVG with `gel_to_svg()` or `write_gel_svg()`.
The PySide6 `AgaroseGelWidget` can render the same model and export PNG files.
