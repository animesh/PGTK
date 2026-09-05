#!/usr/bin/env python3
import gzip
import json
import sys
from pathlib import Path

findings_path = Path(sys.argv[1])
geometry_path = findings_path.parents[1] / "event_geometry.json"
geometry = json.loads(geometry_path.read_text()) if geometry_path.is_file() else {}
errors = []
total = 0
structural = 0
with gzip.open(findings_path, "rt") as handle:
    for total, line in enumerate(handle, 1):
        row = json.loads(line)
        event_id = row["EventID"]
        alt = int(row["ExactAltReads"]); ref = int(row["CleanReferenceReads"]); excluded = int(row["ExcludedReads"])
        callable_count = int(row.get("CallableReads", row["CallableAlignments"]))
        examined = int(row.get("TotalReadsExamined", row["UniqueAlignments"]))
        expected_status = "MIXED_ALT_AND_REFERENCE" if alt and ref else "ALT_SUPPORTED" if alt else "NO_EXACT_ALT_SUPPORT" if ref else "NO_CALLABLE_READS"
        if callable_count != alt + ref or examined != callable_count + excluded or row["ReadValidationStatus"] != expected_status:
            errors.append(event_id + ":count_or_status")
        fraction = row.get("AltFractionAmongClean")
        if callable_count == 0 and fraction not in (None, "NA", "N/A", ""):
            errors.append(event_id + ":undefined_fraction")
        if geometry:
            item = geometry.get(event_id)
            if not item or not item.get("regions"):
                errors.append(event_id + ":missing_geometry"); continue
            event_type = item.get("event_type")
            if row.get("EventType") != event_type: errors.append(event_id + ":event_type")
            for region in item["regions"]:
                if not region.get("chrom") or int(region["start0"]) < 0 or int(region["end0"]) <= int(region["start0"]): errors.append(event_id + ":invalid_region")
            if event_type == "FUSION":
                structural += 1
                if len(item["regions"]) < 2 or {region["role"] for region in item["regions"]} != {"BREAKPOINT_1", "BREAKPOINT_2"}: errors.append(event_id + ":fusion_breakpoints")
            elif event_type == "SPLICE_JUNCTION":
                structural += 1
                if any(region["role"] != "JUNCTION" for region in item["regions"]): errors.append(event_id + ":splice_geometry")
            if event_type in {"FUSION", "SPLICE_JUNCTION", "CONTEXT_EVENT"}:
                expected_visual = "CONTEXT_ALIGNMENTS_AVAILABLE" if int(row.get("ContextAlignments", 0)) else "NO_CONTEXT_ALIGNMENTS"
                if row.get("VisualEvidenceStatus") != expected_visual: errors.append(event_id + ":visual_status")
if geometry and len(geometry) != total: errors.append(f"geometry_count:{len(geometry)}!={total}")
if errors: raise SystemExit(f"FAIL: {len(errors)} invariant errors; first={errors[:10]}")
print(f"published finding validation: PASS ({total} records; structural={structural}; geometry={len(geometry)})")
