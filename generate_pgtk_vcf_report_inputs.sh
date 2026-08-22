#!/usr/bin/env bash
set -Eeuo pipefail

PROJECT="/cluster/projects/nn9036k/scrbkup/pgtk/checkport"
RESULTS="$PROJECT/results"
IMAGE="/cluster/projects/nn9036k/scrbkup/pgtk/singularity_cache/quay.io-biocontainers-bcftools-1.21--h8b25389_0.img"
OUTPUT="$RESULTS/report_inputs"
ARCHIVE="$HOME/pgtk_report_inputs_19246217.tgz"
LOG="$HOME/pgtk_report_inputs_19246217.log"

exec > >(tee "$LOG") 2>&1

cd "$PROJECT"

test -f "$IMAGE"
test -d "$RESULTS"

rm -rf "$OUTPUT"
mkdir -p "$OUTPUT/vcf" "$OUTPUT/go" "$OUTPUT/navigation" "$OUTPUT/provenance"

apptainer exec \
    --no-home \
    --pid \
    -B /cluster \
    -B "$PROJECT:$PROJECT" \
    --pwd "$PROJECT" \
    "$IMAGE" \
    python3 - <<'PY'
import csv
import gzip
import json
import re
import subprocess
from collections import Counter
from pathlib import Path

project = Path("/cluster/projects/nn9036k/scrbkup/pgtk/checkport")
results = project / "results"
output = results / "report_inputs"
vcf_output = output / "vcf"

def run_bcftools(arguments):
    completed = subprocess.run(
        ["bcftools", *arguments],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    )
    return completed.stdout

def infer_sample(path):
    match = re.search(r"(TK[0-9]+)", path.name)
    return match.group(1) if match else "unknown"

def infer_stage(path):
    relative = path.relative_to(results).as_posix()

    mappings = [
        ("gvcf/", "gvcf"),
        ("vcf_raw/", "raw"),
        ("vcf_normalized/", "normalized"),
        ("vcf_filtered/", "filtered"),
        ("vcf_pass/", "pass"),
        ("vep/", "vep"),
        ("rna_validation/variants/", "rna_validated"),
        ("progression_vcf/", "progression"),
    ]

    for prefix, stage in mappings:
        if relative.startswith(prefix):
            return stage

    return relative.split("/", 1)[0]

def variant_type(ref, alt):
    if alt in {".", "*", "<NON_REF>"} or alt.startswith("<"):
        return "SYMBOLIC"
    if len(ref) == 1 and len(alt) == 1:
        return "SNP"
    if len(ref) < len(alt):
        return "INSERTION"
    if len(ref) > len(alt):
        return "DELETION"
    if len(ref) == len(alt) and len(ref) > 1:
        return "MNP"
    return "COMPLEX"

def parse_csq_format(header):
    match = re.search(
        r'##INFO=<ID=CSQ[^>]*Format: ([^">]+)',
        header,
    )
    if not match:
        return []
    return [field.strip() for field in match.group(1).split("|")]

all_rows = []
all_chromosome_rows = []
all_type_rows = []
all_filter_rows = []
all_consequence_rows = []
all_impact_rows = []
all_gene_rows = []
inventory_rows = []

vcfs = sorted(results.rglob("*.vcf.gz"))

for index, vcf in enumerate(vcfs, start=1):
    relative = vcf.relative_to(results).as_posix()
    sample = infer_sample(vcf)
    stage = infer_stage(vcf)

    print(f"[{index}/{len(vcfs)}] {relative}", flush=True)

    header = run_bcftools(["view", "-h", str(vcf)])
    csq_fields = parse_csq_format(header)
    csq_index = {name: position for position, name in enumerate(csq_fields)}

    chromosome_counts = Counter()
    type_counts = Counter()
    chromosome_type_counts = Counter()
    filter_counts = Counter()
    consequence_counts = Counter()
    impact_counts = Counter()
    gene_counts = Counter()

    records = 0
    alleles = 0

    query = run_bcftools(
        [
            "query",
            "-f",
            r"%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/CSQ\n",
            str(vcf),
        ]
    )

    for line in query.splitlines():
        fields = line.split("\t", 5)
        if len(fields) < 6:
            continue

        chrom, position, ref, alt_field, filter_value, csq_value = fields
        records += 1
        chromosome_counts[chrom] += 1
        filter_counts[filter_value or "."] += 1

        for alt in alt_field.split(","):
            current_type = variant_type(ref, alt)
            alleles += 1
            type_counts[current_type] += 1
            chromosome_type_counts[(chrom, current_type)] += 1

        if csq_fields and csq_value not in {"", "."}:
            seen_annotations = set()

            for annotation in csq_value.split(","):
                values = annotation.split("|")

                def value(name):
                    position = csq_index.get(name)
                    if position is None or position >= len(values):
                        return ""
                    return values[position].strip()

                consequence = value("Consequence")
                impact = value("IMPACT")
                gene = (
                    value("SYMBOL")
                    or value("Gene")
                    or value("HGNC_ID")
                )

                annotation_key = (consequence, impact, gene)
                if annotation_key in seen_annotations:
                    continue
                seen_annotations.add(annotation_key)

                for term in consequence.split("&"):
                    if term:
                        consequence_counts[term] += 1

                if impact:
                    impact_counts[impact] += 1

                if gene:
                    gene_counts[gene] += 1

    file_id = re.sub(r"[^A-Za-z0-9_.-]+", "_", relative)

    inventory_rows.append(
        {
            "Sample": sample,
            "Stage": stage,
            "RelativePath": relative,
            "Records": records,
            "Alleles": alleles,
            "Bytes": vcf.stat().st_size,
            "CSQAvailable": bool(csq_fields),
            "CSQFields": "|".join(csq_fields),
        }
    )

    for chrom, count in chromosome_counts.items():
        all_chromosome_rows.append(
            {
                "Sample": sample,
                "Stage": stage,
                "VCF": relative,
                "Chromosome": chrom,
                "Records": count,
            }
        )

    for current_type, count in type_counts.items():
        all_type_rows.append(
            {
                "Sample": sample,
                "Stage": stage,
                "VCF": relative,
                "VariantType": current_type,
                "Alleles": count,
            }
        )

    for (chrom, current_type), count in chromosome_type_counts.items():
        all_rows.append(
            {
                "Sample": sample,
                "Stage": stage,
                "VCF": relative,
                "Chromosome": chrom,
                "VariantType": current_type,
                "Alleles": count,
            }
        )

    for filter_value, count in filter_counts.items():
        all_filter_rows.append(
            {
                "Sample": sample,
                "Stage": stage,
                "VCF": relative,
                "Filter": filter_value,
                "Records": count,
            }
        )

    for consequence, count in consequence_counts.items():
        all_consequence_rows.append(
            {
                "Sample": sample,
                "Stage": stage,
                "VCF": relative,
                "Consequence": consequence,
                "Annotations": count,
            }
        )

    for impact, count in impact_counts.items():
        all_impact_rows.append(
            {
                "Sample": sample,
                "Stage": stage,
                "VCF": relative,
                "Impact": impact,
                "Annotations": count,
            }
        )

    for gene, count in gene_counts.most_common():
        all_gene_rows.append(
            {
                "Sample": sample,
                "Stage": stage,
                "VCF": relative,
                "Gene": gene,
                "Annotations": count,
            }
        )

    summary = {
        "sample": sample,
        "stage": stage,
        "vcf": relative,
        "records": records,
        "alleles": alleles,
        "variant_types": dict(type_counts),
        "filters": dict(filter_counts),
        "csq_fields": csq_fields,
        "top_consequences": consequence_counts.most_common(50),
        "impacts": dict(impact_counts),
        "top_genes": gene_counts.most_common(100),
    }

    (vcf_output / f"{file_id}.summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )

def write_tsv(path, rows, columns):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)

write_tsv(
    vcf_output / "vcf_inventory.tsv",
    inventory_rows,
    [
        "Sample",
        "Stage",
        "RelativePath",
        "Records",
        "Alleles",
        "Bytes",
        "CSQAvailable",
        "CSQFields",
    ],
)

write_tsv(
    vcf_output / "records_by_chromosome.tsv",
    all_chromosome_rows,
    ["Sample", "Stage", "VCF", "Chromosome", "Records"],
)

write_tsv(
    vcf_output / "variant_types.tsv",
    all_type_rows,
    ["Sample", "Stage", "VCF", "VariantType", "Alleles"],
)

write_tsv(
    vcf_output / "variant_types_by_chromosome.tsv",
    all_rows,
    [
        "Sample",
        "Stage",
        "VCF",
        "Chromosome",
        "VariantType",
        "Alleles",
    ],
)

write_tsv(
    vcf_output / "filters.tsv",
    all_filter_rows,
    ["Sample", "Stage", "VCF", "Filter", "Records"],
)

write_tsv(
    vcf_output / "vep_consequences.tsv",
    all_consequence_rows,
    ["Sample", "Stage", "VCF", "Consequence", "Annotations"],
)

write_tsv(
    vcf_output / "vep_impacts.tsv",
    all_impact_rows,
    ["Sample", "Stage", "VCF", "Impact", "Annotations"],
)

write_tsv(
    vcf_output / "vep_genes.tsv",
    all_gene_rows,
    ["Sample", "Stage", "VCF", "Gene", "Annotations"],
)

print(f"Summarized {len(vcfs)} VCF files", flush=True)
PY

# Copy complete GO and progression tables needed for report redesign.
find results/expression results/progression_biology results/progression_vcf \
    -type f \
    \( -name '*.tsv' -o -name '*.txt' -o -name '*.md' \) \
    -exec cp -aL --parents {} "$OUTPUT/go/" \;

# Create stable navigation to every published result.
{
    printf 'Section\tSample\tStage\tRelativePath\tBytes\tModified\n'

    find results -type f \
        ! -path 'results/report_inputs/*' \
        -exec stat --printf='%n\t%s\t%y\n' {} \; |
    sort |
    awk -F '\t' '
        BEGIN { OFS="\t" }
        {
            path=$1
            sample=""

            if (match(path, /TK[0-9]+/)) {
                sample=substr(path, RSTART, RLENGTH)
            }

            relative=path
            sub(/^results\//, "", relative)

            split(relative, parts, "/")
            section=parts[1]
            stage=(length(parts) > 1 ? parts[2] : parts[1])

            print section, sample, stage, relative, $2, $3
        }
    '
} > "$OUTPUT/navigation/report_index.tsv"

# Full current source manifest with hashes and timestamps.
{
    printf 'Path\tBytes\tModified\tSHA256\n'

    find . -maxdepth 1 -type f \
        \( -name '*.nf' -o -name '*.config' -o -name '*.py' \
        -o -name '*.sh' -o -name '*.slurm' -o -name '*.csv' \
        -o -name '*.yaml' -o -name '*.yml' -o -name '*.json' \
        -o -name '*.md' -o -name '*.xml' \) \
        -print0 |
    sort -z |
    while IFS= read -r -d '' file
    do
        printf '%s\t%s\t%s\t%s\n' \
            "$file" \
            "$(stat -c '%s' "$file")" \
            "$(stat -c '%y' "$file")" \
            "$(sha256sum "$file" | awk '{print $1}')"
    done
} > "$OUTPUT/provenance/script_manifest.tsv"

# Record the exact analysis container.
{
    printf 'Tool\tVersion\tContainerPath\tContainerSHA256\tModified\n'

    version=$(
        apptainer exec \
            --no-home \
            --pid \
            -B /cluster \
            "$IMAGE" \
            bcftools --version |
        head -n 1
    )

    printf 'bcftools\t%s\t%s\t%s\t%s\n' \
        "$version" \
        "$IMAGE" \
        "$(sha256sum "$IMAGE" | awk '{print $1}')" \
        "$(stat -c '%y' "$IMAGE")"
} > "$OUTPUT/provenance/bcftools_runtime.tsv"

{
    printf 'GeneratedUTC\t%s\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
    printf 'GeneratedLocal\t%s\n' "$(date '+%Y-%m-%dT%H:%M:%S%z')"
    printf 'PipelineJobID\t19246217\n'
    printf 'PipelineRevision\t835e76644e\n'
    printf 'NextflowVersion\t26.04.6\n'
} > "$OUTPUT/provenance/report_generation.tsv"

# Package only compact report inputs.
rm -f "$ARCHIVE" "$ARCHIVE.sha256" "$ARCHIVE.contents.txt"

tar -C "$RESULTS" -czf "$ARCHIVE" report_inputs
gzip -t "$ARCHIVE"
sha256sum "$ARCHIVE" > "$ARCHIVE.sha256"
tar -tzf "$ARCHIVE" > "$ARCHIVE.contents.txt"

printf '\nGenerated report inputs:\n'
find "$OUTPUT" -type f | sort

printf '\nArchive:\n'
ls -lh "$ARCHIVE" "$ARCHIVE.sha256" "$ARCHIVE.contents.txt"

printf '\nChecksum:\n'
cat "$ARCHIVE.sha256"
