import csv
import gzip
from pathlib import Path
from contextlib import ExitStack

def get_lane(header: str) -> str:
    """Extracts lane number from SRA-prepended Illumina header."""
    # Example: @SRR31089073.1 ST-E00192:445:H2L22CCXY:7:1101:2503:1749/1
    tokens = header.strip().split()
    if len(tokens) > 1:
        parts = tokens[1].split(":")
        if len(parts) >= 4 and parts[3].isdigit():
            return parts[3]
    return "UNKNOWN"

def split_fastqs_by_lane(input_csv: str, output_csv: str):
    new_rows = []
    
    with open(input_csv, mode="r", encoding="utf-8") as f_in:
        reader = csv.DictReader(f_in)
        fieldnames = reader.fieldnames
        
        for row in reader:
            patient = row["patient"]
            sample = row["sample"]
            fq1_path = Path(row["fastq_1"])
            fq2_path = Path(row["fastq_2"])
            
            print(f"Processing {sample}...")
            
            if not (fq1_path.exists() and fq2_path.exists()):
                print(f"  Missing files for {sample}, skipping.")
                continue

            # Dictionary to store open file handles for each detected lane
            # Structure: { '7': (fq1_out_handle, fq2_out_handle) }
            out_handles = {}
            
            with ExitStack() as stack:
                in_fq1 = stack.enter_context(gzip.open(fq1_path, "rt"))
                in_fq2 = stack.enter_context(gzip.open(fq2_path, "rt"))
                
                while True:
                    # Read 4 lines per record for R1
                    r1_h = in_fq1.readline()
                    if not r1_h: break  # EOF
                    r1_s = in_fq1.readline()
                    r1_p = in_fq1.readline()
                    r1_q = in_fq1.readline()
                    
                    # Read 4 lines per record for R2
                    r2_h = in_fq2.readline()
                    r2_s = in_fq2.readline()
                    r2_p = in_fq2.readline()
                    r2_q = in_fq2.readline()
                    
                    lane = get_lane(r1_h)
                    
                    # If this lane hasn't been seen yet, open new FASTQ output files
                    if lane not in out_handles:
                        out_fq1_path = fq1_path.parent / f"{sample}_L00{lane}_1.fastq.gz"
                        out_fq2_path = fq2_path.parent / f"{sample}_L00{lane}_2.fastq.gz"
                        
                        # Open in binary/compress mode
                        handle1 = stack.enter_context(gzip.open(out_fq1_path, "wt", compresslevel=4))
                        handle2 = stack.enter_context(gzip.open(out_fq2_path, "wt", compresslevel=4))
                        out_handles[lane] = (handle1, handle2)
                        
                        # Register the new lane in the sample sheet
                        new_rows.append({
                            "patient": patient,
                            "sample": sample,
                            "fastq_1": str(out_fq1_path),
                            "fastq_2": str(out_fq2_path),
                            "lane": lane
                        })
                        print(f"  -> Discovered Lane {lane}")

                    # Write records to the respective lane files
                    out_handles[lane][0].writelines([r1_h, r1_s, r1_p, r1_q])
                    out_handles[lane][1].writelines([r2_h, r2_s, r2_p, r2_q])

    # Write the new sample sheet
    with open(output_csv, mode="w", newline="", encoding="utf-8") as f_out:
        writer = csv.DictWriter(f_out, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(new_rows)
        
    print(f"\nDone. Updated sample sheet saved to {output_csv}")

if __name__ == "__main__":
    split_fastqs_by_lane("PGTK/sample_laneNA.csv", "PGTK/samples_lane_split.csv")