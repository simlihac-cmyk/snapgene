"""FASTA import and export for PlasmidLab sequence records."""

from __future__ import annotations

from collections.abc import Iterable
from io import StringIO
from pathlib import Path
from typing import TextIO

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord as BioSeqRecord

from plasmidlab.core import MoleculeTopology, MoleculeType, SequenceRecord


def loads_fasta(
    text: str,
    *,
    molecule_type: MoleculeType | str = MoleculeType.DNA,
    topology: MoleculeTopology | str = MoleculeTopology.LINEAR,
) -> tuple[SequenceRecord, ...]:
    """Parse FASTA text into PlasmidLab records."""

    return parse_fasta(StringIO(text), molecule_type=molecule_type, topology=topology)


def dumps_fasta(records: SequenceRecord | Iterable[SequenceRecord]) -> str:
    """Serialize one or more PlasmidLab records to FASTA text."""

    output = StringIO()
    write_fasta(records, output)
    return output.getvalue()


def read_fasta(
    path: str | Path,
    *,
    molecule_type: MoleculeType | str = MoleculeType.DNA,
    topology: MoleculeTopology | str = MoleculeTopology.LINEAR,
) -> tuple[SequenceRecord, ...]:
    """Read FASTA records from a path."""

    with Path(path).open(encoding="utf-8") as handle:
        return parse_fasta(handle, molecule_type=molecule_type, topology=topology)


def write_fasta(records: SequenceRecord | Iterable[SequenceRecord], target: str | Path | TextIO) -> None:
    """Write one or more PlasmidLab records to FASTA."""

    bio_records = [_to_biopython_record(record) for record in _as_record_tuple(records)]
    if hasattr(target, "write"):
        SeqIO.write(bio_records, target, "fasta")
        return

    with Path(target).open("w", encoding="utf-8") as handle:
        SeqIO.write(bio_records, handle, "fasta")


def parse_fasta(
    handle: TextIO,
    *,
    molecule_type: MoleculeType | str = MoleculeType.DNA,
    topology: MoleculeTopology | str = MoleculeTopology.LINEAR,
) -> tuple[SequenceRecord, ...]:
    """Parse FASTA records from a text handle."""

    return tuple(
        SequenceRecord(
            id=bio_record.id,
            name=bio_record.name if bio_record.name and bio_record.name != "<unknown name>" else None,
            description=_description_from_biopython(bio_record),
            sequence=str(bio_record.seq),
            molecule_type=molecule_type,
            topology=topology,
        )
        for bio_record in SeqIO.parse(handle, "fasta")
    )


def _to_biopython_record(record: SequenceRecord) -> BioSeqRecord:
    description = record.description or record.name or record.id
    return BioSeqRecord(
        Seq(record.sequence),
        id=record.id,
        name=record.name or record.id,
        description=description,
    )


def _as_record_tuple(records: SequenceRecord | Iterable[SequenceRecord]) -> tuple[SequenceRecord, ...]:
    if isinstance(records, SequenceRecord):
        return (records,)
    return tuple(records)


def _description_from_biopython(bio_record: BioSeqRecord) -> str | None:
    description = bio_record.description
    if not description or description == "<unknown description>":
        return None
    prefix = f"{bio_record.id} "
    if description.startswith(prefix):
        description = description[len(prefix) :]
    return description or None
