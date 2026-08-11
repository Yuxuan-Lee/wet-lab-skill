# sanger-sequencing-qc

**Repository:** [github.com/Yuxuan-Lee/sanger-sequencing-qc-skill](https://github.com/Yuxuan-Lee/sanger-sequencing-qc-skill)

[English](#english) | [中文](#中文)

---

## English

A general **two-stage Sanger sequencing QC** skill (Cursor/Codex) plus a local Python CLI.

You tell the tool which DNA intervals to verify and which reads belong to which sample. It screens identity / indels / SNVs, checks chromatogram peak purity, and produces a **return-plasmid** list.

### Two core inputs (you do not have to hand-write the CSVs)

The pipeline needs:

1. **`targets`** — sequences to verify  
2. **`assignments`** — which read files belong to which target / clone / primer  

**Easiest path with an agent:** describe the situation in plain language (and point at folders/files). Ask the agent to draft `targets.csv` / `assignments.csv` (or FASTA) for you, then **confirm** before trusting results.

Tell the agent things like:

- where the sequencing reads live (`.seq` / `.ab1`)
- which sequences you want checked (paste them, or point to a FASTA/CSV/SnapGene export)
- which files belong together as one clone/sample (e.g. “these three files are T7 / MID / T7T for clone-1 of geneA”)
- which reads are off-target / wrong locus and should be **excluded**
- that a long insert may need **2–4 primers**; same `clone_id`, multiple rows

**Hard requirement:** every `targets` sequence must be an interval your assigned primers can actually cover. Do not include unsequenced flanks and then demand full coverage.

#### `targets` format (CSV or FASTA)

```csv
target_id,sequence
geneA,ATGC...
geneB,ATGC...
```

#### `assignments` format

```csv
file,target_id,clone_id,primer
xxx_T7_A01.seq,geneA,clone-1,T7
xxx_mid_A02.seq,geneA,clone-1,MID
xxx_T7T_A03.seq,geneA,clone-1,T7T
```

- Off-target reads → omit from that target’s rows  
- Multi-primer coverage → same `clone_id`, different `primer` (union)  
- Same-primer replicates → separate files are kept as separate evidence units  

Examples: [`examples/`](examples/).

### Pipeline

| Stage | Input | Role |
|-------|--------|------|
| Stage-1 | `.seq` | Identity + indel/SNV (primary) |
| Stage-2 | `.ab1` (worth list from stage-1) | Correct-base peak purity ≥ 75% + insertion chromatogram check |
| Final | clone tables | target × **sample** (all primers merged) → return list |

Within one `clone_id`:

- **Base correctness** = reference-**position** union (any reliable read calling the correct base is enough)
- **Insertion** = reference-**boundary** event (local flank gate; **not** whole-target 90% coverage)

Details: [SKILL.md](SKILL.md) · [reference.md](reference.md) · [toolkit/README.md](toolkit/README.md)

### Install

```bash
pip install -r toolkit/requirements.txt
```

Cursor skill: copy this repo to `~/.cursor/skills/sanger-sequencing-qc/`  
(`SKILL.md` `name` must match the folder name.)

### Run

```bash
# 1) seq
python toolkit/stage1_seq_matrix.py \
  --targets targets.csv \
  --assignments assignments.csv \
  --reads-dir /path/to/reads \
  --out-dir stage1_out

# 2) ab1
python toolkit/stage2_ab1_matrix.py \
  --targets targets.csv \
  --assignments assignments.csv \
  --reads-dir /path/to/reads \
  --worth-csv stage1_out/stage1_files_worth_ab1.csv \
  --out-dir stage2_out \
  --peak-min 0.75

# 3) return list
python toolkit/export_final_sample_matrix.py \
  --stage1-clones stage1_out/stage1_clone_union.csv \
  --stage2-clones stage2_out/stage2_clone_union.csv \
  --targets targets.csv \
  --out-dir final_out
```

### License

MIT — see [LICENSE](LICENSE).

---

## 中文

通用的 **Sanger 测序两阶段 QC**（Cursor/Codex skill + 本地 Python CLI）。

你说明「要验证哪些序列」和「哪些测序结果属于哪个样品」，工具完成比对筛选，并给出可返样质粒清单。

### 两份核心输入（不必自己手搓 CSV）

流水线需要：

1. **`targets`**：待验证序列  
2. **`assignments`**：测序文件归属（target / clone / primer）

**最省事的用法：** 用自然语言把情况告诉 agent（并指出目录/文件），让 agent 帮你整理成 `targets.csv` / `assignments.csv`（或 FASTA），**你确认后再跑**。

你可以这样说：

- 测序结果目录在哪（`.seq` / `.ab1`）
- 想检查哪些序列（直接粘贴，或指向 FASTA/CSV 等）
- 哪些文件属于同一个克隆/样品（例如 geneA 的 clone-1 用了 T7 / MID / T7T 三条）
- 哪些结果测飞了、应**排除**
- 长片段可能需要 **2–4 条引物**：同一 `clone_id` 多行即可

**硬性要求：** `targets` 里的序列必须是测序引物实际能覆盖到的区间；不要塞进没人测到的区段却要求全覆盖。

#### `targets` 格式（CSV 或 FASTA）

```csv
target_id,sequence
geneA,ATGC...
geneB,ATGC...
```

#### `assignments` 格式

```csv
file,target_id,clone_id,primer
xxx_T7_A01.seq,geneA,clone-1,T7
xxx_mid_A02.seq,geneA,clone-1,MID
xxx_T7T_A03.seq,geneA,clone-1,T7T
```

- 测飞的 read → 不要写进该 target  
- 多引物覆盖 → 同一 `clone_id`、不同 `primer`（并集）  
- 同引物复测 → 不同文件会作为独立 evidence 保留  

示例见 [`examples/`](examples/)。

### 流水线

| 阶段 | 输入 | 作用 |
|------|------|------|
| Stage-1 | `.seq` | 一致性 + indel/SNV（主判） |
| Stage-2 | `.ab1`（stage-1 值得继续的文件） | 正确碱基峰纯度 ≥ 75% + insertion 峰验证 |
| Final | 克隆表 | target × **样品**（多引物合并）→ 返样清单 |

同一 `clone_id` 下：

- **碱基正确性**：按 reference **位点**并集（任一可靠 read 读对即可）
- **Insertion**：按 reference **边界**单独合并（局部侧翼门控；**不再**用整段 target 90% 覆盖率）

详见 [SKILL.md](SKILL.md) · [reference.md](reference.md) · [toolkit/README.md](toolkit/README.md)

### 安装

```bash
pip install -r toolkit/requirements.txt
```

Cursor skill：拷到 `~/.cursor/skills/sanger-sequencing-qc/`（目录名须与 `SKILL.md` 的 `name` 一致）。

### 运行

```bash
# 1) seq
python toolkit/stage1_seq_matrix.py \
  --targets targets.csv \
  --assignments assignments.csv \
  --reads-dir /path/to/reads \
  --out-dir stage1_out

# 2) ab1
python toolkit/stage2_ab1_matrix.py \
  --targets targets.csv \
  --assignments assignments.csv \
  --reads-dir /path/to/reads \
  --worth-csv stage1_out/stage1_files_worth_ab1.csv \
  --out-dir stage2_out \
  --peak-min 0.75

# 3) 返样表
python toolkit/export_final_sample_matrix.py \
  --stage1-clones stage1_out/stage1_clone_union.csv \
  --stage2-clones stage2_out/stage2_clone_union.csv \
  --targets targets.csv \
  --out-dir final_out
```

### 许可

MIT（见 [LICENSE](LICENSE)）。
