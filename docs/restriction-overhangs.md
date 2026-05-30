# Restriction Cut Geometry and Overhang Semantics

PlasmidLab models restriction digest products as double-stranded DNA fragments with
explicit physical ends. Compatibility decisions must come from the actual template
sequence and the enzyme cut geometry, not from static enzyme overhang labels.

## Coordinates

All coordinates are zero-based, half-open positions on the stored top-strand sequence.
For each restriction site, PlasmidLab records:

- `recognition_start` and `recognition_end`: the matched recognition sequence range.
- `top_cut`: the nick position on the stored top strand.
- `bottom_cut`: the nick position on the opposite strand, expressed in the same
  top-strand coordinate system.
- `strand`: the recognition-site orientation, `+1` for the stored recognition sequence
  and `-1` for a reverse-complement match.

For circular molecules, public cut coordinates are normalized into the sequence length.
The cut metadata also retains the absolute cut calculation used to derive the nick-to-nick
sequence when a cut crosses the origin.

## Physical Fragment Ends

A digest fragment has a `left_end` and a `right_end`. These ends are not interchangeable.
An end can be:

- `none`: an uncut linear molecule edge.
- `blunt`: top and bottom strands are nicked at the same coordinate.
- `5_prime`: the protruding strand has a free 5-prime end.
- `3_prime`: the protruding strand has a free 3-prime end.

`FragmentEnd.overhang_sequence` is the protruding single-stranded sequence read in its
own physical 5-prime-to-3-prime direction. `top_strand_overhang_sequence` is the same
nick-to-nick interval read from the stored template top strand. For palindromic enzymes
such as EcoRI these are often the same string; for non-palindromic Type IIS overhangs
they can differ by reverse complement and by fragment side.

## Type IIS Enzymes

Type IIS enzymes, such as BsaI, cut outside their recognition sequence. Static metadata
may describe their overhang as an ambiguous placeholder such as `NNNN`. PlasmidLab does
not use that placeholder for ligation or Golden Gate decisions. Instead, it derives the
actual overhang from the template bases between `top_cut` and `bottom_cut`.

For example, a BsaI site with top/bottom nicks flanking template bases `ATGC` produces
physical ends derived from `ATGC`; it is never treated as a generic `NNNN` end.

## Compatibility

Two fragment ends are compatible only when:

- one is a left end and the other is a right end;
- neither end is an uncut `none` edge;
- both ends are blunt, or both have the same overhang type;
- sticky overhangs are physically complementary after accounting for strand direction;
- ambiguous bases such as `N` are rejected unless the caller explicitly enables
  ambiguous matching.

This means PlasmidLab compares physical complementarity, end side, and overhang type
together. It does not accept ligation based on overhang string equality alone.
