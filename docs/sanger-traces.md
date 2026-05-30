# Alignment and Sanger Trace Support

PlasmidLab provides pure core reference-alignment utilities and AB1 trace import for
Sanger verification workflows. These tools are independent of the GUI and can be
tested with synthetic sequences or mocked ABIF dictionaries.

## Pairwise Alignment

`plasmidlab.core.alignment` exposes:

- `pairwise_align(reference, query, mode=...)`
- `align_to_reference(reference, observed, qualities=...)`
- `discrepancy_report_tsv(result)`

Supported modes are `global`, `local`, and `semi_global`. Semi-global alignment is
implemented as global alignment with free terminal gaps, which is useful when aligning
a read or imported fragment against a longer reference.

Reference coordinates in discrepancy reports are zero-based, matching PlasmidLab's
internal coordinate convention. Insertions use the reference coordinate before which
the inserted observed base appears; an insertion after the final base is reported at
`len(reference)`.

## AB1 Trace Import

`plasmidlab.io.sanger.read_ab1(path)` reads ABI/AB1 files with Biopython. It extracts:

- called bases from the ABI record or `PBAS2`/`PBAS1`
- quality values from Biopython `phred_quality` annotations or `PCON2`/`PCON1`
- peak positions from `PLOC2`/`PLOC1`
- chromatogram signal channels from `DATA9` through `DATA12`, using `FWO_1` channel
  order when available

PlasmidLab does not perform de novo base calling. The called sequence comes from the
trace file's existing base calls.

## Trace Verification

`align_trace_to_reference(reference, trace)` aligns the trace's called bases to a
reference `SequenceRecord` or string and carries quality values into mismatch and
insertion discrepancies when available.

The PySide6 `SangerTraceWidget` renders the parsed chromatogram channels, called bases,
and quality values. It is a display widget only; biological logic and trace parsing
remain in `core` and `io`.

## Test Data Policy

The committed tests use mocked ABIF raw dictionaries and synthetic traces. To add real
`.ab1` fixtures, only use files with clear redistribution permission, place them under
`tests/data/`, document the license/source beside the fixture, and add a parser test
that calls `read_ab1()`.
