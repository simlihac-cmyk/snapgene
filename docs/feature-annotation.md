# Automatic Feature Annotation

PlasmidLab can propose sequence annotations from open, user-editable feature libraries.
The detector is deterministic and lives in `plasmidlab.core.feature_annotation`; the GUI
only reviews and applies proposals.

## Library Format

The built-in open demo library is packaged under
`plasmidlab.data.features` from `src/plasmidlab/data/features/*.json`, so it is available
from source checkouts, editable installs, wheels, and PyInstaller bundles. User-provided
feature libraries are ordinary JSON files that can live anywhere on disk and be passed
to `load_feature_library(path_or_directory)`.

A file may contain a single feature object, a list of feature objects, or an object with
a `features` list.

Each feature entry supports:

```json
{
  "name": "T7 promoter fragment",
  "type": "promoter",
  "sequence": "TAATACGACTCACTATAGGG",
  "aliases": ["T7 promoter"],
  "strand_behavior": "both",
  "minimum_identity": 0.95,
  "notes": "Canonical short T7 promoter consensus fragment."
}
```

`strand_behavior` may be:

- `both`: search the sequence and reverse complement; annotate `+` or `-`
- `forward_only`: search only the sequence as written
- `reverse_only`: search only the reverse complement
- `non_directional`: search both orientations but annotate strand as unknown

`minimum_identity` is a fraction from `0` to `1`. Use `1.0` for exact matching, or a
lower threshold for near-exact matching.

## Detection Behavior

The detector supports linear and circular DNA records. Circular matches that cross the
origin are reported as compound features with two zero-based, half-open segments, such
as `[(10, 13), (0, 6)]`.

Duplicate proposal resolution keeps the best match for repeated detections at the same
site, preferring higher identity and longer matches.

Applying proposals returns a new `SequenceRecord` with a `detect_features` provenance
event. Existing identical annotations are not duplicated.

## Legal Boundary

The bundled `plasmidlab.data.features/common_synthetic_features.json` file is a small
open synthetic demo library. PlasmidLab does not import, copy, reverse engineer, or
depend on proprietary feature databases. Labs can maintain their own JSON files outside
the package and pass custom file or directory paths to `load_feature_library()`.
