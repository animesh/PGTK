import csv
from pathlib import Path

fusion_dir = Path("resultsbkup/rna_validation/fusions")

if not fusion_dir.is_dir():
    raise SystemExit(f"ERROR: directory not found: {fusion_dir.resolve()}")

files = sorted(fusion_dir.glob("*.fusion.validated.tsv"))

if not files:
    raise SystemExit(f"ERROR: no validated fusion files found in: {fusion_dir.resolve()}")

total_interchromosomal = 0

for path in files:
    with path.open(encoding="utf-8", errors="replace", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    interchromosomal = []

    for row in rows:
        breakpoint1 = row.get("breakpoint1", row.get("Breakpoint1", "")).strip()
        breakpoint2 = row.get("breakpoint2", row.get("Breakpoint2", "")).strip()

        chromosome1 = breakpoint1.split(":", 1)[0]
        chromosome2 = breakpoint2.split(":", 1)[0]

        if chromosome1 and chromosome2 and chromosome1 != chromosome2:
            interchromosomal.append(row)

    total_interchromosomal += len(interchromosomal)

    print(
        path.name,
        f"validated={len(rows)}",
        f"interchromosomal={len(interchromosomal)}",
    )

    for row in interchromosomal:
        split_reads = (
            int(row.get("split_reads1", "0") or 0)
            + int(row.get("split_reads2", "0") or 0)
        )

        print(
            " ",
            row.get("#gene1", row.get("gene1", "")),
            row.get("gene2", ""),
            row.get("breakpoint1", ""),
            row.get("breakpoint2", ""),
            f"split_reads={split_reads}",
            f"discordant_mates={row.get('discordant_mates', '')}",
            f"confidence={row.get('confidence', '')}",
            f"reading_frame={row.get('reading_frame', '')}",
        )

print(f"total_interchromosomal_sample_events={total_interchromosomal}")
