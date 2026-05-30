from plasmidlab.cli import main
from plasmidlab.core import MoleculeTopology, SequenceRecord
from plasmidlab.io.genbank import write_genbank
from plasmidlab.io.fasta import read_fasta


def test_digest_cli_outputs_fragment_table(tmp_path, capsys) -> None:
    sequence = (
        "A" * 10
        + "GAATTC"
        + "C" * 24
        + "GGATCC"
        + "T" * 24
        + "AAGCTT"
        + "G" * 24
    )
    path = tmp_path / "input.gb"
    write_genbank(
        SequenceRecord(id="cli_digest", sequence=sequence, topology=MoleculeTopology.CIRCULAR),
        path,
    )

    exit_code = main(["digest", str(path), "--enzymes", "EcoRI,BamHI"])

    output = capsys.readouterr().out
    assert exit_code == 0
    assert "start\tend\tlength" in output
    assert "11\t41\t30\tEcoRI,BamHI" in output
    assert "41\t11\t70\tBamHI,EcoRI" in output


def test_pcr_cli_outputs_product_table(tmp_path, capsys) -> None:
    sequence = (
        "TTGACCGTAA"
        "GCTAGTCCATGA"
        "CCGTAAGGTTCCGATACCGGTTACGATCGGATATCGAA"
        "ATGACCTGATCG"
        "GGTTAACC"
    )
    path = tmp_path / "template.gb"
    write_genbank(
        SequenceRecord(id="cli_pcr", sequence=sequence, topology=MoleculeTopology.CIRCULAR),
        path,
    )

    exit_code = main(
        [
            "pcr",
            str(path),
            "--fwd",
            sequence[60:72],
            "--rev",
            _reverse_complement(sequence[10:22]),
            "--min-anneal-length",
            "12",
        ]
    )

    output = capsys.readouterr().out
    assert exit_code == 0
    assert "start\tend\tlength\twraps_origin\tsequence" in output
    assert f"60\t22\t42\ttrue\t{sequence[60:] + sequence[:22]}" in output


def test_convert_cli_writes_fasta_from_genbank(tmp_path, capsys) -> None:
    input_path = tmp_path / "input.gb"
    output_path = tmp_path / "output.fasta"
    write_genbank(SequenceRecord(id="convert_me", sequence="ATGCATGC"), input_path)

    exit_code = main(["convert", str(input_path), str(output_path)])

    output = capsys.readouterr().out
    converted = read_fasta(output_path)
    assert exit_code == 0
    assert "converted\t1\tfasta" in output
    assert converted[0].id == "convert_me"
    assert converted[0].sequence == "ATGCATGC"


def test_map_export_cli_writes_svg(tmp_path, capsys) -> None:
    input_path = tmp_path / "map.gb"
    output_path = tmp_path / "map.svg"
    write_genbank(
        SequenceRecord(id="map_me", sequence="ATGC" * 20, topology=MoleculeTopology.CIRCULAR),
        input_path,
    )

    exit_code = main(["map-export", str(input_path), str(output_path)])

    output = capsys.readouterr().out
    assert exit_code == 0
    assert "map_svg\tmap_me" in output
    assert output_path.read_text(encoding="utf-8").startswith('<?xml version="1.0"')


def test_plasmidlab_without_subcommand_prints_help(capsys) -> None:
    exit_code = main([])

    output = capsys.readouterr().out
    assert exit_code == 0
    assert "digest" in output
    assert "convert" in output


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]
