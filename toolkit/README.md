# Toolkit CLI

## Required user inputs

- `--targets` CSV/FASTA (`target_id` + `sequence`)
- `--assignments` CSV (`file,target_id,clone_id,primer`)
- `--reads-dir` with `.seq` / `.ab1`

See `../examples/`.

## Commands

```bash
pip install -r requirements.txt

python stage1_seq_matrix.py --targets T.csv --assignments A.csv --reads-dir READS --out-dir stage1_out
python stage2_ab1_matrix.py --targets T.csv --assignments A.csv --reads-dir READS --worth-csv stage1_out/stage1_files_worth_ab1.csv --out-dir stage2_out --peak-min 0.75
python export_final_sample_matrix.py --stage1-clones stage1_out/stage1_clone_union.csv --stage2-clones stage2_out/stage2_clone_union.csv --targets T.csv --out-dir final_out
```

Optional HER2-style filename helper only:

```bash
python draft_her2_inputs.py --insert-csv pcr_insert_flanks.csv --reads-dir READS --aliases aliases.json --out-dir drafted_inputs
```

## Insertion thresholds

See `insertion.py` (`LOCAL_IDENTITY_MIN`, `FLANK_MATCH_MIN_BP`, `EDGE_MARGIN_BP`, …).
Base union and insertion evidence are separate tracks.

## Tests

```bash
python -m unittest tests.test_insertion -v
```
