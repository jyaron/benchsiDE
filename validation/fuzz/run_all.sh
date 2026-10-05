#!/bin/sh
# Randomized comparison with limma/edgeR for one seed: python3/Rscript on PATH or RSCRIPT set.
# Usage: sh validation/fuzz/run_all.sh <seed>   (from the repository root, after tests/engine/fetch_data.py)
set -e
RS=${RSCRIPT:-Rscript}; SEED=$1
python3 validation/fuzz/gen_configs.py $SEED
python3 -c "import json;print('\n'.join(l['name'] for l in json.load(open('validation/fuzz/work/configs.json'))['loads']))" > validation/fuzz/work/loads.txt
xargs -P 8 -I{} sh -c 'python3 validation/fuzz/run_app.py {} > validation/fuzz/work/{}.app.log 2>&1' < validation/fuzz/work/loads.txt
xargs -P 8 -I{} sh -c 'python3 validation/fuzz/run_ds.py {} > validation/fuzz/work/{}.ds.log 2>&1' < validation/fuzz/work/loads.txt
xargs -P 8 -I{} sh -c 'python3 validation/fuzz/run_sets.py {} > validation/fuzz/work/{}.sets.log 2>&1' < validation/fuzz/work/loads.txt
xargs -P 8 -I{} sh -c "$RS validation/fuzz/reference.R {} > validation/fuzz/work/{}.ref.log 2>&1" < validation/fuzz/work/loads.txt
xargs -P 8 -I{} sh -c "$RS validation/fuzz/reference.R {} sets > validation/fuzz/work/{}.setsref.log 2>&1" < validation/fuzz/work/loads.txt
echo "== seed $SEED"; python3 validation/fuzz/compare.py 2>&1 | grep -v Warning; python3 validation/fuzz/compare_ds.py; python3 validation/fuzz/compare_sets.py | head -4
