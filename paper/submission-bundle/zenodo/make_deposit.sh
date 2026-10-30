#!/bin/bash
# Build the Zenodo tarball from the live project and validate the payload offline.
# Deposit (sandbox/production) and publication are separate, explicit steps — never run here.
set -e
HERE=$(cd "$(dirname "$0")" && pwd); ROOT=$(cd "$HERE/../../.." && pwd)
cd "$ROOT"
tar czf "$HERE/skinheat-v0.1.0.tar.gz" \
  --exclude='*.pkl' --exclude='__pycache__' --exclude='*.egg-info' --exclude='*.bak' --exclude='results/NOTES.md' \
  --exclude='paper/submission-bundle' --exclude='paper/*.aux' --exclude='paper/*.log' --exclude='paper/*.out' --exclude='paper/*.blg' --exclude='paper/*Notes.bib' --exclude='paper/*.compile.log' \
  LICENSE README.md Makefile pyproject.toml skinheat scripts tests results paper
ls -la "$HERE/skinheat-v0.1.0.tar.gz"
cd "$HERE" && /opt/anaconda3/bin/zenodo-deposit --config submission.json --dry-run
