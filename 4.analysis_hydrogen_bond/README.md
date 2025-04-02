# 氢键后处理分析脚本（适用于 gmx hbond 输出）

本脚本用于对 GROMACS 中 `gmx hbond` 命令输出的氢键数据进行可视化与统计分析。

---

## 💡 使用前提

请确保已正确使用如下命令生成氢键信息：

```bash
gmx hbond -s Structure.tpr -f trajectory.xtc -n index.ndx \
  -num prefix.xvg -hbm prefix.xpm -hbn prefix.ndx
```

并将 `prefix.xpm` 与 `prefix.ndx` 放在当前工作目录下。

---

## 📦 依赖库

- pandas：用于读取数据表格  
- numpy：处理数值型数组  
- matplotlib：绘图（热图、条形图等）  
- seaborn（可选）：美化图像  
- argparse：解析命令行参数  
- os / sys：文件与路径处理  

安装命令：

```bash
pip install pandas numpy matplotlib seaborn
```

---

## 🚀 脚本使用方式

```bash
python analyze_hbond_pairs.py \
  --prefixes a,b \
  --Struc_file structure.pdb \
  --summary_occupation \
  --hydrogen_distribution \
  --hydrogen_occupation \
  --hydrogen_occupation_single \
  --height_range 45,67,89
```

---

## ⚙️ 参数说明

| 参数 | 类型 | 说明 |
|------|------|------|
| `--prefixes` | str | 输入待分析的多个前缀，逗号分隔，例如 `"a,b,c"` |
| `--Struc_file` | str | 参考结构文件（.pdb），需与原始 tpr 对应 |
| `--summary_occupation` | flag | 绘制所有轨迹中单残基氢键热图 |
| `--hydrogen_distribution` | flag | 每对氢键随时间分布图 |
| `--hydrogen_occupation` | flag | 每对氢键在整个轨迹中的占有率 |
| `--hydrogen_occupation_single` | flag | 统计单残基参与氢键的频率 |
| `--height_range` | str | 高亮的残基编号，逗号分隔，例如 `"123,456"` |

---

## 📊 输出结果

根据你选择的功能，将生成以下图像（可能包含）：

- `heatmap.png`：多轨迹下单残基氢键热图
- `a_distribution.png`：轨迹 `a` 中每对氢键随时间出现频率
- `a_hydrogen_occupation.png`：氢键占有率条形图
- `a_single_residue_hydrogen_occupation.png`：每个残基的氢键参与情况

---

## 📬 联系方式
作者：zhouyq
邮箱：zhouyq@shanghaitech.edu.cn  
更新日期：2025年4月

---
   
   
   
---

# Post-processing Script for `gmx hbond` Output

This script provides post-analysis and visualization tools for hydrogen bond data produced by the `gmx hbond` command in GROMACS.

---

## 💡 Prerequisites

Make sure you have run the following `gmx hbond` command:

```bash
gmx hbond -s Structure.tpr -f trajectory.xtc -n index.ndx \
  -num prefix.xvg -hbm prefix.xpm -hbn prefix.ndx
```

Ensure `prefix.xpm` and `prefix.ndx` are in the working directory.

---

## 📦 Dependencies

- `pandas`: for handling tabular data  
- `numpy`: numerical computations  
- `matplotlib`: plotting (heatmaps, bar plots)  
- `seaborn` (optional): enhanced plot styling  
- `argparse`: parsing command-line arguments  
- `os` / `sys`: file and system operations  

Install with:

```bash
pip install pandas numpy matplotlib seaborn
```

---

## 🚀 Usage

```bash
python analyze_hbond_pairs.py \
  --prefixes a,b \
  --Struc_file structure.pdb \
  --summary_occupation \
  --hydrogen_distribution \
  --hydrogen_occupation \
  --hydrogen_occupation_single \
  --height_range 45,67,89
```

---

## ⚙️ Argument Descriptions

| Argument | Type | Description |
|----------|------|-------------|
| `--prefixes` | str | Comma-separated file prefixes, e.g., `"a,b,c"` |
| `--Struc_file` | str | Reference `.pdb` structure file (should match `tpr`) |
| `--summary_occupation` | flag | Plot heatmap of single-residue hydrogen bond occupation |
| `--hydrogen_distribution` | flag | Time-resolved hydrogen bond distribution per pair |
| `--hydrogen_occupation` | flag | Overall occupancy of each hydrogen bond |
| `--hydrogen_occupation_single` | flag | Frequency of each residue’s involvement in hydrogen bonding |
| `--height_range` | str | Highlight specific residues on Y-axis, e.g., `"123,456"` |

---

## 📊 Output

Depending on the selected options, the following files may be generated:

- `heatmap.png`: heatmap of single-residue H-bond occupation across trajectories
- `a_distribution.png`: time-resolved H-bond presence in trajectory `a`
- `a_hydrogen_occupation.png`: bar chart of H-bond occupancy
- `a_single_residue_hydrogen_occupation.png`: summary by individual residue

---

## 📬 Contact

Author: zhouyq  
Email: zhouyq@shanghaitech.edu.cn  
Last Updated: April 2025
