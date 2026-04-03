#!/usr/bin/env bash
# Quick status check for the bedcov background job
BEDCOV_DIR="/mnt/z/Download/TK/results/coverage/bedcov"
DONE=$(ls "$BEDCOV_DIR"/*.bedcov.tsv 2>/dev/null | wc -l)
NONZERO=$(ls -s "$BEDCOV_DIR"/*.bedcov.tsv 2>/dev/null | awk '$1>0' | wc -l)
RUNNING=$(ps aux | grep "samtools bedcov" | grep -v grep | wc -l)
echo "[$( date '+%H:%M:%S')] BAMs: $DONE/28 started  |  $NONZERO/28 complete  |  $RUNNING samtools processes running"
if [ -f "$BEDCOV_DIR/DONE" ]; then
    echo "  ✅ ALL DONE — run: python3 /mnt/z/Download/TK/plot_windowed_cn.py"
else
    tail -3 /mnt/z/Download/TK/results/coverage/bedcov_run.log
fi
