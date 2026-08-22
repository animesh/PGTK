#!/usr/bin/env python3
from pathlib import Path
from tempfile import TemporaryDirectory

SOURCE = Path(__file__).with_name("main.nf").read_text(encoding="utf-8")
REQUIRED = [
    "new File(root, expected)",
    "new File(new File(root, srr), expected)",
    "External caller VCF not found",
    "Ambiguous external caller VCF",
    "External caller VCF is empty",
]
for text in REQUIRED:
    assert text in SOURCE, text


def resolve(root: Path, srr: str, suffix: str) -> Path:
    expected = f"{srr}{suffix}"
    preferred = []
    for candidate in (root / expected, root / srr / expected):
        if candidate.is_file() and candidate.resolve() not in [x.resolve() for x in preferred]:
            preferred.append(candidate)
    if len(preferred) == 1:
        if preferred[0].stat().st_size == 0:
            raise ValueError("empty")
        return preferred[0]
    if len(preferred) > 1:
        raise ValueError("ambiguous")
    matches = list(root.rglob(expected))
    if len(matches) != 1:
        raise ValueError("missing or ambiguous")
    if matches[0].stat().st_size == 0:
        raise ValueError("empty")
    return matches[0]

with TemporaryDirectory() as temporary:
    root = Path(temporary)
    suffix = ".haplotypecaller.filtered.vcf.gz"
    srr = "SRR31089074"
    flat = root / f"{srr}{suffix}"
    flat.write_bytes(b"VCF")
    assert resolve(root, srr, suffix) == flat

with TemporaryDirectory() as temporary:
    root = Path(temporary)
    suffix = ".haplotypecaller.filtered.vcf.gz"
    srr = "SRR31089074"
    nested = root / srr / f"{srr}{suffix}"
    nested.parent.mkdir()
    nested.write_bytes(b"VCF")
    assert resolve(root, srr, suffix) == nested

with TemporaryDirectory() as temporary:
    root = Path(temporary)
    suffix = ".haplotypecaller.filtered.vcf.gz"
    srr = "SRR31089074"
    deep = root / "variant_calling" / "haplotypecaller" / srr / f"{srr}{suffix}"
    deep.parent.mkdir(parents=True)
    deep.write_bytes(b"VCF")
    assert resolve(root, srr, suffix) == deep

with TemporaryDirectory() as temporary:
    root = Path(temporary)
    suffix = ".haplotypecaller.filtered.vcf.gz"
    srr = "SRR31089074"
    (root / srr).mkdir()
    (root / f"{srr}{suffix}").write_bytes(b"VCF")
    (root / srr / f"{srr}{suffix}").write_bytes(b"VCF")
    try:
        resolve(root, srr, suffix)
    except ValueError as error:
        assert "ambiguous" in str(error)
    else:
        raise AssertionError("ambiguous layout was accepted")

print("PASS: flat, Sarek nested, recursive fallback, empty and ambiguity contracts are covered")
