# sanger-stage1-seq-matrix

[中文](#中文) | [English](#english)

---

## 中文

Cursor / Codex 可用的 **Sanger 测序两阶段 QC skill**，也可作为本地 Python CLI。面向 binder / 引物合成批次：用 `flank + CDS + flank` 参考序列筛克隆，并输出**可返样质粒**清单。

### 能做什么

1. **Stage-1（`.seq`）**：全长一致性 + **indel / SNV**（主判，更简单）
2. **Stage-2（`.ab1`）**：正确碱基**峰纯度 ≥ 75%**
3. **Final**：基因 × **样品**矩阵（T7/T7T 合并为同一样品）→ 收集返样质粒

每阶段维护同一套矩阵视图：哪些测序文件值得继续；哪些目标基因尚未整条失败。

### 仓库结构

```text
sanger-stage1-seq-matrix/
├── README.md
├── LICENSE
├── SKILL.md              # Agent skill 说明（Cursor / Codex）
├── reference.md          # 判定规则细节
└── toolkit/              # 可独立运行的脚本
    ├── requirements.txt
    ├── aliases.example.json
    ├── stage1_seq_matrix.py
    ├── stage2_ab1_matrix.py
    ├── export_final_sample_matrix.py
    └── README.md
```

### 安装依赖

```bash
pip install -r toolkit/requirements.txt
```

需要：Python 3.9+、`biopython`。

### 安装为 Cursor skill

把本仓库拷到 Cursor skills 目录（名称保持文件夹名即可）：

```powershell
# Windows PowerShell
Copy-Item -Recurse .\sanger-stage1-seq-matrix $env:USERPROFILE\.cursor\skills\
```

```bash
# macOS / Linux
mkdir -p ~/.cursor/skills
cp -r sanger-stage1-seq-matrix ~/.cursor/skills/
```

新开对话后，用自然语言描述「Sanger / seq / ab1 / 返样质粒筛选」即可触发；也可直接运行 `toolkit/` 下脚本。

### 本地 CLI 用法

准备：

- 测序目录：含 `*.seq` 与对应 `*.ab1`（文件名含样品、T7/T7T）
- `pcr_insert_flanks.csv`：至少含 `vendor_id`、`full_insert_dna`
- 可选 `aliases.json`：样品前缀笔误 → `vendor_id`（见 `toolkit/aliases.example.json`）

```bash
# 1) seq 初筛（identity + indel/SNV）
python toolkit/stage1_seq_matrix.py \
  --seq-dir /path/to/reads \
  --insert-csv /path/to/pcr_insert_flanks.csv \
  --aliases /path/to/aliases.json \
  --out-dir stage1_out

# 2) ab1 峰纯度（仅 stage1 合格文件）
python toolkit/stage2_ab1_matrix.py \
  --ab1-dir /path/to/reads \
  --insert-csv /path/to/pcr_insert_flanks.csv \
  --worth-csv stage1_out/stage1_files_worth_ab1.csv \
  --aliases /path/to/aliases.json \
  --out-dir stage2_out \
  --peak-min 0.75

# 3) 最终返样表（基因 × 样品）
python toolkit/export_final_sample_matrix.py \
  --stage1-clones stage1_out/stage1_clone_union.csv \
  --stage2-clones stage2_out/stage2_clone_union.csv \
  --genes-csv /path/to/pcr_insert_flanks.csv \
  --out-dir final_out
```

打开 `final_out/final_gene_sample_matrix.html` 或 `final_plasmid_return_list.csv` 即可按样品收质粒。

### 判定要点（摘要）

- 参考序列 = **侧翼 + CDS + 侧翼**（`full_insert_dna`）
- 正反向并集：任一端读对即该位点成功；未覆盖或读错 = 失败
- Stage-1 额外抓**可信内部插入**（高置信比对：match≥90% 参考长度且 identity≥0.95）
- Stage-2：正确碱基通道高度 / 四通道之和 ≥ **0.75**
- 详情见 [reference.md](reference.md)、[SKILL.md](SKILL.md)

### 别人如何使用

1. `git clone` 本仓库（或 Download ZIP）
2. `pip install -r toolkit/requirements.txt`
3. 按上面 CLI 跑，或拷到 `~/.cursor/skills/` 当 skill 用

### 许可

MIT（见 [LICENSE](LICENSE)）。

---

## English

A **two-stage Sanger QC skill** (Cursor / Codex) and local Python CLI for binder / primer synthesis batches. Screens clones against `flank + CDS + flank` references and exports a **plasmid-return** list.

### Pipeline

1. **Stage-1 (`.seq`)**: full-insert identity + **indel/SNV** (primary, simpler)
2. **Stage-2 (`.ab1`)**: correct-base **peak fraction ≥ 75%**
3. **Final**: gene × **sample** matrix (T7/T7T merged) for return plasmids

### Install

```bash
pip install -r toolkit/requirements.txt
```

Copy this folder to `~/.cursor/skills/sanger-stage1-seq-matrix/` to use as a Cursor skill.

### CLI

See the Chinese section above for the three commands (`stage1_seq_matrix.py` → `stage2_ab1_matrix.py` → `export_final_sample_matrix.py`). Details: [toolkit/README.md](toolkit/README.md), [SKILL.md](SKILL.md), [reference.md](reference.md).

### License

MIT — see [LICENSE](LICENSE).
