"""Restriction enzyme search and digest simulation."""

from __future__ import annotations

import re
from collections import defaultdict
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from enum import StrEnum
from typing import Any

from Bio.Restriction import AllEnzymes

from plasmidlab.core.models import MoleculeTopology, MoleculeType, SequenceRecord


IUPAC_DNA: Mapping[str, str] = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "AG",
    "Y": "CT",
    "S": "GC",
    "W": "AT",
    "K": "GT",
    "M": "AC",
    "B": "CGT",
    "D": "AGT",
    "H": "ACT",
    "V": "ACG",
    "N": "ACGT",
}

DNA_COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


class OverhangKind(StrEnum):
    """Restriction fragment end type."""

    NONE = "none"
    BLUNT = "blunt"
    FIVE_PRIME = "5_prime"
    THREE_PRIME = "3_prime"


class FragmentEndSide(StrEnum):
    """Physical side of a digested fragment."""

    LEFT = "left"
    RIGHT = "right"


def _coerce_overhang_kind(kind: OverhangKind | str) -> OverhangKind:
    try:
        return kind if isinstance(kind, OverhangKind) else OverhangKind(kind)
    except ValueError as error:
        msg = f"unsupported overhang kind: {kind!r}"
        raise ValueError(msg) from error


def _coerce_fragment_end_side(side: FragmentEndSide | str) -> FragmentEndSide:
    try:
        return side if isinstance(side, FragmentEndSide) else FragmentEndSide(side)
    except ValueError as error:
        msg = f"unsupported fragment end side: {side!r}"
        raise ValueError(msg) from error


@dataclass(frozen=True, slots=True)
class Overhang:
    """A fragment-end overhang.

    ``NONE`` is used for uncut linear molecule edges. Restriction-cut ends use
    ``BLUNT``, ``FIVE_PRIME``, or ``THREE_PRIME``.
    """

    kind: OverhangKind = OverhangKind.NONE
    sequence: str = ""

    def __post_init__(self) -> None:
        object.__setattr__(self, "kind", _coerce_overhang_kind(self.kind))
        object.__setattr__(self, "sequence", self.sequence.upper())
        if self.kind in (OverhangKind.NONE, OverhangKind.BLUNT) and self.sequence:
            msg = "none and blunt overhangs must not have a sequence"
            raise ValueError(msg)

    @property
    def length(self) -> int:
        """Return the overhang length."""

        return len(self.sequence)


NO_OVERHANG = Overhang()


@dataclass(frozen=True, slots=True)
class CutGeometry:
    """Actual top-strand and bottom-strand nick positions for a restriction site."""

    enzyme_name: str
    recognition_start: int
    recognition_end: int
    strand: int
    top_cut: int
    bottom_cut: int
    recognition_sequence: str
    cut_source_metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "enzyme_name", _validate_name(self.enzyme_name, "enzyme name"))
        if self.strand not in (-1, 1):
            msg = "cut geometry strand must be +1 or -1"
            raise ValueError(msg)
        for field_name in ("recognition_start", "recognition_end", "top_cut", "bottom_cut"):
            value = getattr(self, field_name)
            if not isinstance(value, int):
                msg = f"{field_name} must be an integer"
                raise TypeError(msg)
            if value < 0:
                msg = f"{field_name} must be non-negative"
                raise ValueError(msg)
        object.__setattr__(self, "recognition_sequence", self.recognition_sequence.upper())
        object.__setattr__(self, "cut_source_metadata", dict(self.cut_source_metadata))


@dataclass(frozen=True, slots=True)
class FragmentEnd:
    """A physical fragment end produced by a cut or by an uncut molecule edge.

    ``overhang_sequence`` is the protruding single-stranded sequence in its own
    physical 5-prime-to-3-prime direction. ``top_strand_overhang_sequence`` stores the
    same nick-to-nick interval in template top-strand coordinates for reporting and
    deterministic assembly calculations.
    """

    side: FragmentEndSide
    kind: OverhangKind
    overhang_sequence: str = ""
    top_cut: int | None = None
    bottom_cut: int | None = None
    source_enzyme: str | None = None
    source_site_coordinate: int | None = None
    orientation: int = 0
    top_strand_overhang_sequence: str = ""
    cut_geometry: CutGeometry | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "side", _coerce_fragment_end_side(self.side))
        object.__setattr__(self, "kind", _coerce_overhang_kind(self.kind))
        overhang_sequence = self.overhang_sequence.upper()
        top_strand_sequence = (self.top_strand_overhang_sequence or overhang_sequence).upper()
        object.__setattr__(self, "overhang_sequence", overhang_sequence)
        object.__setattr__(self, "top_strand_overhang_sequence", top_strand_sequence)
        if self.kind in (OverhangKind.NONE, OverhangKind.BLUNT):
            if overhang_sequence or top_strand_sequence:
                msg = "none and blunt fragment ends must not have overhang sequences"
                raise ValueError(msg)
        elif not overhang_sequence:
            msg = "sticky fragment ends must include an overhang sequence"
            raise ValueError(msg)
        for field_name in ("top_cut", "bottom_cut", "source_site_coordinate"):
            value = getattr(self, field_name)
            if value is not None and (not isinstance(value, int) or value < 0):
                msg = f"{field_name} must be a non-negative integer or None"
                raise ValueError(msg)
        if self.orientation not in (-1, 0, 1):
            msg = "fragment end orientation must be +1, -1, or 0"
            raise ValueError(msg)
        if self.source_enzyme is not None:
            object.__setattr__(self, "source_enzyme", _validate_name(self.source_enzyme, "source enzyme"))

    @property
    def sequence(self) -> str:
        """Return the physical overhang sequence for backward-compatible reporting."""

        return self.overhang_sequence

    @property
    def overhang(self) -> Overhang:
        """Return the end as the legacy lightweight overhang value."""

        return Overhang(self.kind, self.overhang_sequence)


UNCUT_LEFT_END = FragmentEnd(FragmentEndSide.LEFT, OverhangKind.NONE)
UNCUT_RIGHT_END = FragmentEnd(FragmentEndSide.RIGHT, OverhangKind.NONE)


@dataclass(frozen=True, slots=True)
class RestrictionEnzyme:
    """Open restriction enzyme metadata used by PlasmidLab."""

    name: str
    recognition_sequence: str
    top_cut_offset: int
    bottom_cut_offset: int
    overhang: Overhang
    methylation_sensitive: bool | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "name", _validate_name(self.name, "enzyme name"))
        recognition_sequence = self.recognition_sequence.upper()
        invalid = sorted(set(recognition_sequence) - set(IUPAC_DNA))
        if invalid:
            msg = f"invalid recognition sequence characters: {''.join(invalid)}"
            raise ValueError(msg)
        object.__setattr__(self, "recognition_sequence", recognition_sequence)
        if not isinstance(self.top_cut_offset, int) or not isinstance(self.bottom_cut_offset, int):
            msg = "cut offsets must be integers"
            raise TypeError(msg)
        if not isinstance(self.overhang, Overhang):
            msg = "overhang must be an Overhang"
            raise TypeError(msg)

    @property
    def recognition_length(self) -> int:
        """Return the recognition sequence length."""

        return len(self.recognition_sequence)

    @property
    def is_blunt(self) -> bool:
        """Return whether this enzyme produces blunt ends."""

        return self.overhang.kind is OverhangKind.BLUNT

    @property
    def is_sticky(self) -> bool:
        """Return whether this enzyme produces sticky ends."""

        return self.overhang.kind in (OverhangKind.FIVE_PRIME, OverhangKind.THREE_PRIME)


@dataclass(frozen=True, slots=True)
class RestrictionSite:
    """A restriction enzyme recognition site and its cut positions."""

    enzyme: RestrictionEnzyme
    start: int
    end: int
    strand: int
    top_cut: int
    bottom_cut: int
    wraps_origin: bool = False
    top_cut_absolute: int | None = None
    bottom_cut_absolute: int | None = None

    def __post_init__(self) -> None:
        if self.strand not in (-1, 1):
            msg = "restriction site strand must be +1 or -1"
            raise ValueError(msg)
        for field_name in ("start", "end", "top_cut", "bottom_cut"):
            value = getattr(self, field_name)
            if not isinstance(value, int):
                msg = f"{field_name} must be an integer"
                raise TypeError(msg)
            if value < 0:
                msg = f"{field_name} must be non-negative"
                raise ValueError(msg)
        for field_name in ("top_cut_absolute", "bottom_cut_absolute"):
            value = getattr(self, field_name)
            if value is not None and not isinstance(value, int):
                msg = f"{field_name} must be an integer or None"
                raise TypeError(msg)

    @property
    def enzyme_name(self) -> str:
        """Return the enzyme name."""

        return self.enzyme.name

    @property
    def recognition_sequence(self) -> str:
        """Return the enzyme recognition sequence."""

        return self.enzyme.recognition_sequence

    @property
    def cut_positions(self) -> tuple[int, int]:
        """Return top-strand and bottom-strand cut positions."""

        return (self.top_cut, self.bottom_cut)

    @property
    def overhang(self) -> Overhang:
        """Return enzyme-level overhang metadata.

        This is descriptive only. Digest and cloning compatibility use sequence-derived
        :class:`FragmentEnd` values instead.
        """

        return self.enzyme.overhang

    @property
    def cut_geometry(self) -> CutGeometry:
        """Return actual nick geometry for the site."""

        return CutGeometry(
            enzyme_name=self.enzyme_name,
            recognition_start=self.start,
            recognition_end=self.end,
            strand=self.strand,
            top_cut=self.top_cut,
            bottom_cut=self.bottom_cut,
            recognition_sequence=self.recognition_sequence,
            cut_source_metadata={
                "top_cut_absolute": self.top_cut_absolute,
                "bottom_cut_absolute": self.bottom_cut_absolute,
                "top_cut_offset": self.enzyme.top_cut_offset,
                "bottom_cut_offset": self.enzyme.bottom_cut_offset,
            },
        )

    @property
    def methylation_sensitive(self) -> bool | None:
        """Return methylation-sensitivity metadata when known."""

        return self.enzyme.methylation_sensitive


@dataclass(frozen=True, slots=True)
class DigestFragment:
    """A fragment produced by an in-silico digest."""

    start: int
    end: int
    length: int
    sequence: str
    left_end: FragmentEnd = UNCUT_LEFT_END
    right_end: FragmentEnd = UNCUT_RIGHT_END
    source_enzymes: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if self.start < 0 or self.end < 0 or self.length < 0:
            msg = "fragment coordinates and length must be non-negative"
            raise ValueError(msg)
        object.__setattr__(self, "sequence", self.sequence.upper())
        object.__setattr__(self, "source_enzymes", tuple(dict.fromkeys(self.source_enzymes)))
        if self.length != len(self.sequence):
            msg = "fragment length must match sequence length"
            raise ValueError(msg)
        if self.left_end.side is not FragmentEndSide.LEFT:
            msg = "left_end must have side='left'"
            raise ValueError(msg)
        if self.right_end.side is not FragmentEndSide.RIGHT:
            msg = "right_end must have side='right'"
            raise ValueError(msg)

    @property
    def source_start(self) -> int:
        """Return the source-record start coordinate."""

        return self.start

    @property
    def source_end(self) -> int:
        """Return the source-record end coordinate."""

        return self.end

    @property
    def left_overhang(self) -> Overhang:
        """Return the left physical end as the legacy overhang value."""

        return self.left_end.overhang

    @property
    def right_overhang(self) -> Overhang:
        """Return the right physical end as the legacy overhang value."""

        return self.right_end.overhang


@dataclass(frozen=True, slots=True)
class RestrictionAnalysis:
    """Restriction-site search result with cutter classifications."""

    record_id: str
    sequence_length: int
    topology: MoleculeTopology
    enzymes: tuple[RestrictionEnzyme, ...]
    sites: tuple[RestrictionSite, ...]

    @property
    def sites_by_enzyme(self) -> dict[str, tuple[RestrictionSite, ...]]:
        """Return sites grouped by enzyme name."""

        grouped: dict[str, list[RestrictionSite]] = {enzyme.name: [] for enzyme in self.enzymes}
        for site in self.sites:
            grouped.setdefault(site.enzyme_name, []).append(site)
        return {name: tuple(sites) for name, sites in grouped.items()}

    @property
    def non_cutters(self) -> tuple[str, ...]:
        """Return enzymes with no sites."""

        return tuple(name for name, sites in self.sites_by_enzyme.items() if len(sites) == 0)

    @property
    def single_cutters(self) -> tuple[str, ...]:
        """Return enzymes with exactly one site."""

        return tuple(name for name, sites in self.sites_by_enzyme.items() if len(sites) == 1)

    @property
    def multi_cutters(self) -> tuple[str, ...]:
        """Return enzymes with more than one site."""

        return tuple(name for name, sites in self.sites_by_enzyme.items() if len(sites) > 1)


EnzymeInput = str | RestrictionEnzyme | Any


def restriction_enzyme_from_biopython(enzyme: Any) -> RestrictionEnzyme:
    """Create PlasmidLab enzyme metadata from a Biopython enzyme class."""

    recognition_sequence = str(enzyme.site).upper()
    top_cut_offset = int(enzyme.fst5)
    bottom_cut_offset = len(recognition_sequence) + int(enzyme.fst3)
    return RestrictionEnzyme(
        name=enzyme.__name__,
        recognition_sequence=recognition_sequence,
        top_cut_offset=top_cut_offset,
        bottom_cut_offset=bottom_cut_offset,
        overhang=_overhang_from_offsets(recognition_sequence, top_cut_offset, bottom_cut_offset),
        methylation_sensitive=None,
    )


def resolve_enzymes(enzymes: EnzymeInput | Iterable[EnzymeInput]) -> tuple[RestrictionEnzyme, ...]:
    """Resolve names, Biopython enzyme classes, or custom definitions."""

    if isinstance(enzymes, str):
        raw_items: Iterable[EnzymeInput] = [item.strip() for item in enzymes.split(",") if item.strip()]
    elif isinstance(enzymes, RestrictionEnzyme):
        raw_items = (enzymes,)
    else:
        raw_items = enzymes

    resolved: list[RestrictionEnzyme] = []
    seen: set[str] = set()
    for item in raw_items:
        enzyme = _resolve_single_enzyme(item)
        if enzyme.name not in seen:
            resolved.append(enzyme)
            seen.add(enzyme.name)
    return tuple(resolved)


def analyze_restriction_sites(
    record: SequenceRecord,
    enzymes: EnzymeInput | Iterable[EnzymeInput],
) -> RestrictionAnalysis:
    """Find restriction enzyme sites on a linear or circular DNA record."""

    _require_dna(record)
    enzyme_definitions = resolve_enzymes(enzymes)
    sites = tuple(
        sorted(
            (
                site
                for enzyme in enzyme_definitions
                for site in _find_sites_for_enzyme(record, enzyme)
            ),
            key=lambda site: (site.start, site.enzyme_name, site.strand),
        )
    )
    return RestrictionAnalysis(
        record_id=record.id,
        sequence_length=record.length,
        topology=record.topology,
        enzymes=enzyme_definitions,
        sites=sites,
    )


def find_restriction_sites(
    record: SequenceRecord,
    enzymes: EnzymeInput | Iterable[EnzymeInput],
) -> tuple[RestrictionSite, ...]:
    """Return restriction sites for the requested enzymes."""

    return analyze_restriction_sites(record, enzymes).sites


def digest(
    record: SequenceRecord,
    enzymes: EnzymeInput | Iterable[EnzymeInput],
) -> tuple[DigestFragment, ...]:
    """Simulate a single, double, or multi-enzyme digest."""

    analysis = analyze_restriction_sites(record, enzymes)
    cuts_by_position = _cuts_by_top_position(analysis.sites, record.length, record.topology)
    if not cuts_by_position:
        return (
            DigestFragment(
                start=0,
                end=record.length,
                length=record.length,
                sequence=record.sequence,
                left_end=_uncut_fragment_end(FragmentEndSide.LEFT, 0),
                right_end=_uncut_fragment_end(FragmentEndSide.RIGHT, record.length),
            ),
        )

    if record.topology is MoleculeTopology.LINEAR:
        return _digest_linear(record, cuts_by_position)
    return _digest_circular(record, cuts_by_position)


def _find_sites_for_enzyme(
    record: SequenceRecord,
    enzyme: RestrictionEnzyme,
) -> tuple[RestrictionSite, ...]:
    sequence = record.sequence
    sequence_length = record.length
    recognition = enzyme.recognition_sequence
    recognition_length = enzyme.recognition_length
    if sequence_length == 0 or recognition_length == 0:
        return ()

    search_sequence = (
        sequence + sequence[: recognition_length - 1]
        if record.topology is MoleculeTopology.CIRCULAR and recognition_length > 1
        else sequence
    )
    patterns = ((recognition, 1),)
    reverse = _reverse_complement(recognition)
    if reverse != recognition:
        patterns = patterns + ((reverse, -1),)

    sites: list[RestrictionSite] = []
    seen: set[tuple[str, int, int]] = set()
    for pattern, strand in patterns:
        regex = _iupac_regex(pattern)
        for match in regex.finditer(search_sequence):
            start = match.start()
            if start >= sequence_length:
                continue
            end_abs = start + recognition_length
            wraps_origin = record.topology is MoleculeTopology.CIRCULAR and end_abs > sequence_length
            if record.topology is MoleculeTopology.LINEAR and end_abs > sequence_length:
                continue
            top_cut_abs, bottom_cut_abs = _cut_positions_for_match(enzyme, start, strand)
            top_cut = _normalize_cut(top_cut_abs, sequence_length, record.topology)
            bottom_cut = _normalize_cut(bottom_cut_abs, sequence_length, record.topology)
            if top_cut is None or bottom_cut is None:
                continue
            site_end = end_abs % sequence_length if wraps_origin else end_abs
            key = (enzyme.name, start, strand)
            if key in seen:
                continue
            sites.append(
                RestrictionSite(
                    enzyme=enzyme,
                    start=start,
                    end=site_end,
                    strand=strand,
                    top_cut=top_cut,
                    bottom_cut=bottom_cut,
                    wraps_origin=wraps_origin,
                    top_cut_absolute=top_cut_abs,
                    bottom_cut_absolute=bottom_cut_abs,
                )
            )
            seen.add(key)
    return tuple(sorted(sites, key=lambda site: (site.start, site.strand)))


def _digest_linear(
    record: SequenceRecord,
    cuts_by_position: Mapping[int, tuple[RestrictionSite, ...]],
) -> tuple[DigestFragment, ...]:
    cut_positions = sorted(position for position in cuts_by_position if 0 < position < record.length)
    boundaries = (0, *cut_positions, record.length)
    fragments: list[DigestFragment] = []
    for start, end in zip(boundaries, boundaries[1:], strict=False):
        if start == end:
            continue
        left_sites = cuts_by_position.get(start, ())
        right_sites = cuts_by_position.get(end, ())
        sequence = record.sequence[start:end]
        fragments.append(
            DigestFragment(
                start=start,
                end=end,
                length=len(sequence),
                sequence=sequence,
                left_end=_representative_fragment_end(
                    record,
                    left_sites,
                    FragmentEndSide.LEFT,
                    start,
                ),
                right_end=_representative_fragment_end(
                    record,
                    right_sites,
                    FragmentEndSide.RIGHT,
                    end,
                ),
                source_enzymes=_source_enzyme_names(left_sites, right_sites),
            )
        )
    return tuple(fragments)


def _digest_circular(
    record: SequenceRecord,
    cuts_by_position: Mapping[int, tuple[RestrictionSite, ...]],
) -> tuple[DigestFragment, ...]:
    cut_positions = sorted(cuts_by_position)
    if len(cut_positions) == 1:
        cut = cut_positions[0]
        sequence = record.sequence[cut:] + record.sequence[:cut]
        sites = cuts_by_position[cut]
        return (
            DigestFragment(
                start=cut,
                end=cut,
                length=len(sequence),
                sequence=sequence,
                left_end=_representative_fragment_end(record, sites, FragmentEndSide.LEFT, cut),
                right_end=_representative_fragment_end(record, sites, FragmentEndSide.RIGHT, cut),
                source_enzymes=_source_enzyme_names(sites),
            ),
        )

    fragments: list[DigestFragment] = []
    for index, start in enumerate(cut_positions):
        end = cut_positions[(index + 1) % len(cut_positions)]
        sequence = _circular_subsequence(record.sequence, start, end)
        left_sites = cuts_by_position[start]
        right_sites = cuts_by_position[end]
        fragments.append(
            DigestFragment(
                start=start,
                end=end,
                length=len(sequence),
                sequence=sequence,
                left_end=_representative_fragment_end(record, left_sites, FragmentEndSide.LEFT, start),
                right_end=_representative_fragment_end(record, right_sites, FragmentEndSide.RIGHT, end),
                source_enzymes=_source_enzyme_names(left_sites, right_sites),
            )
        )
    return tuple(fragments)


def _cuts_by_top_position(
    sites: Iterable[RestrictionSite],
    sequence_length: int,
    topology: MoleculeTopology,
) -> dict[int, tuple[RestrictionSite, ...]]:
    grouped: dict[int, list[RestrictionSite]] = defaultdict(list)
    for site in sites:
        cut = site.top_cut
        if topology is MoleculeTopology.LINEAR and not 0 <= cut <= sequence_length:
            continue
        if topology is MoleculeTopology.CIRCULAR:
            cut %= sequence_length
        grouped[cut].append(site)
    return {
        position: tuple(sorted(sites_at_position, key=lambda site: site.enzyme_name))
        for position, sites_at_position in grouped.items()
    }


def _source_enzyme_names(*site_groups: Iterable[RestrictionSite]) -> tuple[str, ...]:
    names: list[str] = []
    for site_group in site_groups:
        for site in site_group:
            if site.enzyme_name not in names:
                names.append(site.enzyme_name)
    return tuple(names)


def _representative_fragment_end(
    record: SequenceRecord,
    sites: Iterable[RestrictionSite],
    side: FragmentEndSide,
    coordinate: int,
) -> FragmentEnd:
    sorted_sites = sorted(sites, key=lambda site: site.enzyme_name)
    if not sorted_sites:
        return _uncut_fragment_end(side, coordinate)
    return fragment_end_for_site(record, sorted_sites[0], side)


def _uncut_fragment_end(side: FragmentEndSide, coordinate: int) -> FragmentEnd:
    return FragmentEnd(
        side=side,
        kind=OverhangKind.NONE,
        top_cut=coordinate,
        bottom_cut=coordinate,
    )


def fragment_end_for_site(
    record: SequenceRecord,
    site: RestrictionSite,
    side: FragmentEndSide | str,
) -> FragmentEnd:
    """Derive a physical fragment end from a site and template sequence."""

    resolved_side = _coerce_fragment_end_side(side)
    top_abs = site.top_cut_absolute if site.top_cut_absolute is not None else site.top_cut
    bottom_abs = site.bottom_cut_absolute if site.bottom_cut_absolute is not None else site.bottom_cut
    if top_abs == bottom_abs:
        return FragmentEnd(
            side=resolved_side,
            kind=OverhangKind.BLUNT,
            top_cut=site.top_cut,
            bottom_cut=site.bottom_cut,
            source_enzyme=site.enzyme_name,
            source_site_coordinate=site.start,
            orientation=site.strand,
            cut_geometry=site.cut_geometry,
        )

    region_start = min(top_abs, bottom_abs)
    region_end = max(top_abs, bottom_abs)
    top_strand_sequence = _subsequence_between_absolute_cuts(
        record.sequence,
        region_start,
        region_end,
        record.topology,
    )
    kind = OverhangKind.FIVE_PRIME if top_abs < bottom_abs else OverhangKind.THREE_PRIME
    if (kind is OverhangKind.FIVE_PRIME and resolved_side is FragmentEndSide.LEFT) or (
        kind is OverhangKind.THREE_PRIME and resolved_side is FragmentEndSide.RIGHT
    ):
        physical_sequence = top_strand_sequence
    else:
        physical_sequence = _reverse_complement(top_strand_sequence)

    return FragmentEnd(
        side=resolved_side,
        kind=kind,
        overhang_sequence=physical_sequence,
        top_cut=site.top_cut,
        bottom_cut=site.bottom_cut,
        source_enzyme=site.enzyme_name,
        source_site_coordinate=site.start,
        orientation=site.strand,
        top_strand_overhang_sequence=top_strand_sequence,
        cut_geometry=site.cut_geometry,
    )


def compatible_fragment_ends(
    first: FragmentEnd,
    second: FragmentEnd,
    *,
    allow_ambiguous: bool = False,
) -> bool:
    """Return whether two physical fragment ends can be ligated."""

    if first.side is second.side:
        return False
    if first.kind is OverhangKind.NONE or second.kind is OverhangKind.NONE:
        return False
    if first.kind is OverhangKind.BLUNT or second.kind is OverhangKind.BLUNT:
        return first.kind is OverhangKind.BLUNT and second.kind is OverhangKind.BLUNT
    if first.kind is not second.kind:
        return False
    if len(first.overhang_sequence) != len(second.overhang_sequence):
        return False
    if not allow_ambiguous and (
        _contains_ambiguous_bases(first.overhang_sequence)
        or _contains_ambiguous_bases(second.overhang_sequence)
    ):
        return False
    return _iupac_sequences_match(
        first.overhang_sequence,
        _reverse_complement(second.overhang_sequence),
        allow_ambiguous=allow_ambiguous,
    )


def _circular_subsequence(sequence: str, start: int, end: int) -> str:
    if start == end:
        return sequence[start:] + sequence[:start]
    if start < end:
        return sequence[start:end]
    return sequence[start:] + sequence[:end]


def _subsequence_between_absolute_cuts(
    sequence: str,
    start: int,
    end: int,
    topology: MoleculeTopology,
) -> str:
    if start >= end:
        return ""
    if topology is MoleculeTopology.LINEAR:
        return sequence[start:end]
    if not sequence:
        return ""
    return "".join(sequence[position % len(sequence)] for position in range(start, end))


def _cut_positions_for_match(
    enzyme: RestrictionEnzyme,
    start: int,
    strand: int,
) -> tuple[int, int]:
    if strand == 1:
        return (
            start + enzyme.top_cut_offset,
            start + enzyme.bottom_cut_offset,
        )
    recognition_length = enzyme.recognition_length
    return (
        start + recognition_length - enzyme.bottom_cut_offset,
        start + recognition_length - enzyme.top_cut_offset,
    )


def _normalize_cut(
    cut: int,
    sequence_length: int,
    topology: MoleculeTopology,
) -> int | None:
    if topology is MoleculeTopology.CIRCULAR:
        if sequence_length == 0:
            return None
        return cut % sequence_length
    if 0 <= cut <= sequence_length:
        return cut
    return None


def _resolve_single_enzyme(item: EnzymeInput) -> RestrictionEnzyme:
    if isinstance(item, RestrictionEnzyme):
        return item
    if isinstance(item, str):
        enzyme_class = _enzyme_lookup().get(item.lower())
        if enzyme_class is None:
            msg = f"unknown restriction enzyme: {item}"
            raise ValueError(msg)
        return restriction_enzyme_from_biopython(enzyme_class)
    if hasattr(item, "site") and hasattr(item, "fst5") and hasattr(item, "fst3"):
        return restriction_enzyme_from_biopython(item)
    msg = f"unsupported restriction enzyme input: {item!r}"
    raise TypeError(msg)


def _enzyme_lookup() -> dict[str, Any]:
    return {enzyme.__name__.lower(): enzyme for enzyme in AllEnzymes}


def _overhang_from_offsets(
    recognition_sequence: str,
    top_cut_offset: int,
    bottom_cut_offset: int,
) -> Overhang:
    if top_cut_offset == bottom_cut_offset:
        return Overhang(OverhangKind.BLUNT)
    kind = (
        OverhangKind.FIVE_PRIME
        if top_cut_offset < bottom_cut_offset
        else OverhangKind.THREE_PRIME
    )
    region_start = min(top_cut_offset, bottom_cut_offset)
    region_end = max(top_cut_offset, bottom_cut_offset)
    sequence = (
        recognition_sequence[region_start:region_end]
        if 0 <= region_start < region_end <= len(recognition_sequence)
        else ""
    )
    return Overhang(kind, sequence)


def _contains_ambiguous_bases(sequence: str) -> bool:
    return any(base not in "ACGT" for base in sequence.upper())


def _iupac_sequences_match(
    observed: str,
    expected: str,
    *,
    allow_ambiguous: bool,
) -> bool:
    if len(observed) != len(expected):
        return False
    if not allow_ambiguous:
        return observed == expected
    return all(
        bool(set(IUPAC_DNA[observed_base]) & set(IUPAC_DNA[expected_base]))
        for observed_base, expected_base in zip(observed.upper(), expected.upper(), strict=True)
    )


def _iupac_regex(pattern: str) -> re.Pattern[str]:
    regex = "".join(f"[{IUPAC_DNA[base]}]" for base in pattern.upper())
    return re.compile(f"(?=({regex}))")


def _reverse_complement(sequence: str) -> str:
    return sequence.upper().translate(DNA_COMPLEMENT)[::-1]


def _validate_name(value: str, field_name: str) -> str:
    if not isinstance(value, str):
        msg = f"{field_name} must be a string"
        raise TypeError(msg)
    if not value or value != value.strip():
        msg = f"{field_name} must be non-empty without leading or trailing whitespace"
        raise ValueError(msg)
    return value


def _require_dna(record: SequenceRecord) -> None:
    if record.molecule_type is not MoleculeType.DNA:
        msg = "restriction enzyme analysis requires a DNA sequence"
        raise ValueError(msg)
