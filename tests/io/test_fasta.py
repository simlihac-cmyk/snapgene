from pathlib import Path

from plasmidlab.core import MoleculeTopology, MoleculeType, SequenceRecord
from plasmidlab.io.fasta import dumps_fasta, loads_fasta, read_fasta, write_fasta


def test_fasta_round_trip_linear_pcr_product() -> None:
    record = SequenceRecord(
        id="pcr_product_1",
        name="PCR_product_1",
        description="Linear PCR product",
        sequence="ATGCGTACGTAGCTAGCTAG",
        molecule_type=MoleculeType.DNA,
        topology=MoleculeTopology.LINEAR,
    )

    imported = loads_fasta(dumps_fasta(record))[0]

    assert imported.id == record.id
    assert imported.sequence == record.sequence
    assert imported.description == record.description
    assert imported.molecule_type is MoleculeType.DNA
    assert imported.topology is MoleculeTopology.LINEAR
    assert imported.features == ()


def test_fasta_file_helpers(tmp_path: Path) -> None:
    record = SequenceRecord(id="rna_fragment", sequence="AUGCAU", molecule_type=MoleculeType.RNA)
    path = tmp_path / "fragment.fasta"

    write_fasta(record, path)
    imported = read_fasta(path, molecule_type=MoleculeType.RNA)

    assert imported[0].sequence == "AUGCAU"
    assert imported[0].molecule_type is MoleculeType.RNA
