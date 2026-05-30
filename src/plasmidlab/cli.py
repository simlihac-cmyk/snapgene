"""Command line interface for PlasmidLab."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

from plasmidlab.core.primers import PCRSimulationError, simulate_pcr
from plasmidlab.core.restriction import digest
from plasmidlab.io.fasta import read_fasta, write_fasta
from plasmidlab.io.genbank import read_genbank
from plasmidlab.io.genbank import write_genbank
from plasmidlab.render import render_plasmid_map, write_svg


def main(argv: Sequence[str] | None = None) -> int:
    """Run the PlasmidLab CLI."""

    parser = _build_parser()
    args = parser.parse_args(argv)
    if not hasattr(args, "func"):
        parser.print_help()
        return 0
    return args.func(args)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="plasmidlab")
    subparsers = parser.add_subparsers()

    digest_parser = subparsers.add_parser("digest", help="simulate restriction digest")
    digest_parser.add_argument("input", type=Path, help="input FASTA or GenBank file")
    digest_parser.add_argument(
        "--enzymes",
        required=True,
        help="comma-separated restriction enzymes, e.g. EcoRI,BamHI",
    )
    digest_parser.set_defaults(func=_digest_command)

    pcr_parser = subparsers.add_parser("pcr", help="simulate PCR")
    pcr_parser.add_argument("template", type=Path, help="input FASTA or GenBank template")
    pcr_parser.add_argument("--fwd", required=True, help="forward primer sequence")
    pcr_parser.add_argument("--rev", required=True, help="reverse primer sequence")
    pcr_parser.add_argument("--mismatches", type=int, default=0, help="allowed mismatches per primer")
    pcr_parser.add_argument(
        "--min-anneal-length",
        type=int,
        default=12,
        help="minimum 3-prime annealing length used for primers with 5-prime tails",
    )
    pcr_parser.add_argument(
        "--max-anneal-length",
        type=int,
        default=None,
        help="maximum 3-prime annealing length to consider; defaults to full primer length",
    )
    pcr_parser.add_argument(
        "--target-region",
        type=_parse_region,
        default=None,
        metavar="START:END",
        help="prefer a product spanning this zero-based half-open template interval",
    )
    pcr_parser.add_argument(
        "--pair-index",
        type=int,
        default=None,
        help="select one product when multiple primer binding pairs are possible",
    )
    pcr_parser.set_defaults(func=_pcr_command)

    convert_parser = subparsers.add_parser("convert", help="convert FASTA and GenBank files")
    convert_parser.add_argument("input", type=Path, help="input FASTA or GenBank file")
    convert_parser.add_argument("output", type=Path, help="output FASTA or GenBank file")
    convert_parser.add_argument(
        "--to",
        choices=("fasta", "genbank"),
        default=None,
        help="output format; inferred from extension when omitted",
    )
    convert_parser.set_defaults(func=_convert_command)

    map_parser = subparsers.add_parser("map-export", help="export a sequence map as SVG")
    map_parser.add_argument("input", type=Path, help="input FASTA or GenBank file")
    map_parser.add_argument("output", type=Path, help="output SVG file")
    map_parser.add_argument(
        "--record-index",
        type=int,
        default=0,
        help="zero-based record index to export from multi-record files",
    )
    map_parser.set_defaults(func=_map_export_command)
    return parser


def _digest_command(args: argparse.Namespace) -> int:
    records = _read_sequence_file(args.input)
    if not records:
        raise SystemExit(f"No sequence records found in {args.input}")
    record = records[0]
    fragments = digest(record, args.enzymes)
    print("start\tend\tlength\tsource_enzymes\tleft_overhang\tright_overhang")
    for fragment in fragments:
        print(
            "\t".join(
                (
                    str(fragment.start),
                    str(fragment.end),
                    str(fragment.length),
                    ",".join(fragment.source_enzymes) or "-",
                    _format_overhang(fragment.left_overhang),
                    _format_overhang(fragment.right_overhang),
                )
            )
        )
    return 0


def _format_overhang(overhang: object) -> str:
    kind = getattr(overhang, "kind")
    sequence = getattr(overhang, "sequence")
    return f"{kind.value}:{sequence}" if sequence else kind.value


def _pcr_command(args: argparse.Namespace) -> int:
    records = _read_sequence_file(args.template)
    if not records:
        raise SystemExit(f"No sequence records found in {args.template}")
    try:
        product = simulate_pcr(
            records[0],
            args.fwd,
            args.rev,
            mismatches=args.mismatches,
            min_anneal_length=args.min_anneal_length,
            max_anneal_length=args.max_anneal_length,
            target_region=args.target_region,
            pair_index=args.pair_index,
        )
    except PCRSimulationError as error:
        raise SystemExit(str(error)) from error

    print("start\tend\tlength\twraps_origin\tsequence")
    print(
        "\t".join(
            (
                str(product.start),
                str(product.end),
                str(product.length),
                str(product.wraps_origin).lower(),
                product.sequence,
            )
        )
    )
    return 0


def _convert_command(args: argparse.Namespace) -> int:
    records = _read_sequence_file(args.input)
    if not records:
        raise SystemExit(f"No sequence records found in {args.input}")
    output_format = args.to or _format_from_path(args.output)
    if output_format == "fasta":
        write_fasta(records, args.output)
    elif output_format == "genbank":
        write_genbank(records, args.output)
    else:
        raise SystemExit(f"Cannot infer output format from {args.output}; use --to")
    print(f"converted\t{len(records)}\t{output_format}\t{args.output}")
    return 0


def _map_export_command(args: argparse.Namespace) -> int:
    records = _read_sequence_file(args.input)
    if not records:
        raise SystemExit(f"No sequence records found in {args.input}")
    if args.record_index < 0 or args.record_index >= len(records):
        raise SystemExit(f"record-index {args.record_index} is outside {len(records)} records")
    record = records[args.record_index]
    write_svg(render_plasmid_map(record), args.output)
    print(f"map_svg\t{record.id}\t{args.output}")
    return 0


def _read_sequence_file(path: Path):
    suffix = path.suffix.lower()
    if suffix in {".fa", ".fasta", ".fna"}:
        return read_fasta(path)
    return read_genbank(path)


def _format_from_path(path: Path) -> str | None:
    suffix = path.suffix.lower()
    if suffix in {".fa", ".fasta", ".fna"}:
        return "fasta"
    if suffix in {".gb", ".gbk", ".genbank"}:
        return "genbank"
    return None


def _parse_region(value: str) -> tuple[int, int]:
    try:
        start_text, end_text = value.split(":", maxsplit=1)
        return (int(start_text), int(end_text))
    except ValueError as error:
        msg = "target regions must use START:END zero-based half-open syntax"
        raise argparse.ArgumentTypeError(msg) from error


if __name__ == "__main__":
    raise SystemExit(main())
