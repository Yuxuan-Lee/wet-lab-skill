# sanger-stage1-seq-matrix / wet-lab-skill

[中文](#中文) | [English](#english)

---

## 中文

通用的 **Sanger 测序两阶段 QC**（Cursor/Codex skill + 本地 Python CLI）。

你提供「想验证的序列」和「哪些测序文件属于哪个样品」，工具做比对筛选，并给出可返样质粒清单。

### 你需要准备的两份核心输入

1. **`targets`（待测序列）**  
   CSV 或 FASTA。内容必须是**测序引物实际能覆盖到的区间**（不要丢一段没人测到的序列却要求全覆盖）。

   CSV 示例：

   ```csv
   target_id,sequence
   geneA,ATGC...
   geneB,ATGC...
   ```

2. **`assignments`（测序归属表）**  
   告诉 agent/脚本：哪几个结果文件属于哪个 target、哪个克隆/样品、哪条引物。  
   - 有的 read 可能测到别处 → **不要写进该 target 的 assignment**  
   - 整段可能需要 2 / 3 / 4 条引物 → **同一 `clone_id` 下挂多行即可**（并集覆盖）

   ```csv
   file,target_id,clone_id,primer
   xxx_T7_A01.seq,geneA,clone-1,T7
   xxx_mid_A02.seq,geneA,clone-1,MID
   xxx_T7T_A03.seq,geneA,clone-1,T7T
   ```

### 流水线分工

| 阶段 | 输入 | 作用 |
|------|------|------|
| Stage-1 | `.seq` | 全长一致性 + indel/SNV（主判） |
| Stage-2 | `.ab1`（仅 stage-1 合格） | 正确碱基峰纯度 ≥ 75% |
| Final | 克隆表 | target × **样品**（多引物合并）→ 返样清单 |

同一 `clone_id` 下：

- **碱基正确性**：按 reference **位点**做并集（任一可靠 read 读对即可）
- **Insertion**：按 reference **边界**单独合并（局部侧翼质量门控；**不再**用整段 target 90% 覆盖率）

详见 [reference.md](reference.md)。

### 安装

```bash
pip install -r toolkit/requirements.txt
```

Cursor skill：把本仓库拷到 `~/.cursor/skills/sanger-stage1-seq-matrix/`（目录名须与 `SKILL.md` 的 `name` 一致）。

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

示例文件见 `examples/`。

### Agent 使用方式

用自然语言说明：

1. targets 文件路径（或粘贴序列让 agent 写成 CSV/FASTA）  
2. 测序结果目录  
3. 哪些文件对应哪个 target/clone（可先让 agent 根据文件名起草 `assignments.csv`，再由你确认）

详细规则： [SKILL.md](SKILL.md) · [reference.md](reference.md) · [toolkit/README.md](toolkit/README.md)

### 许可

MIT（见 [LICENSE](LICENSE)）。

---

## English

General **two-stage Sanger QC** skill/CLI.

Provide:

1. **targets** — sequences to verify (must be coverable by the assigned primers/reads)
2. **assignments** — which read files belong to which `target_id` / `clone_id` / `primer` (2–N primers per clone supported; omit off-target reads)

Then run stage-1 (`.seq` identity + indel/SNV) → stage-2 (`.ab1` peak ≥75%) → final target×sample return list.

See commands above and `examples/`.
