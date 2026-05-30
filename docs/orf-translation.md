# ORF And Translation

PlasmidLab implements ORF and translation tools in `plasmidlab.core.orfs`. These
functions are pure core utilities and do not depend on GUI code.

## Coordinates

All returned ORF coordinates use PlasmidLab's internal zero-based, half-open coordinate
system. ORFs include the terminal stop codon in their nucleotide interval. Circular ORFs
that cross the origin are represented with multiple `FeatureSegment` values, for
example `[14, 18)` followed by `[0, 5)`.

## ORF Detection

`find_orfs()` searches DNA records in three plus-strand frames and, by default, three
minus-strand frames. It detects closed ORFs from an allowed start codon to the first
in-frame stop codon.

Configurable inputs include:

- Start codons, defaulting to `ATG`.
- Stop codons, defaulting to the selected genetic code table.
- Minimum amino acid length, excluding the terminal stop symbol.
- Genetic code table.
- Whether to include reverse-strand frames.

Circular records are searched across the origin without duplicating ORFs that start at
the same template coordinate.

## Translation

Translation helpers include:

- `translate_sequence()` for a DNA string.
- `translate_record()` for a whole record or selected region.
- `translate_feature()` for simple or compound features.
- `extract_region()` for selected-region extraction with circular support.

Translation accepts a genetic code table id or name supported by Biopython. Incomplete
trailing codons are ignored; ambiguous codons translate to `X`.

## CDS Validation

`validate_cds()` checks a candidate CDS for:

- Allowed start codon.
- Allowed terminal stop codon.
- Internal stop codons.
- Nucleotide length divisible by three.

The result is a `CDSValidationResult` containing a boolean `is_valid`, issues, start and
stop codons, internal stop positions, and nucleotide length.

## Reverse Translation

`reverse_translate()` back-translates a protein sequence deterministically. Without codon
usage data it chooses a stable codon from the selected genetic code table. With usage
weights, it chooses the highest-weight codon for each amino acid.

Codon usage can be passed directly as a mapping or loaded from:

- JSON: `{ "codons": { "GCC": 5.0, "GCT": 1.0 } }` or a direct codon-to-weight object.
- CSV: columns named `codon` and one of `weight`, `frequency`, or `usage`.
