import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os,sys
from argparse import RawTextHelpFormatter
def clean_string(string):
    string = string.strip()
    while "  " in string:
        string = string.replace("  "," ")
    return string.split(" ")
def analysis_pdb(pdb_file):
    aa_dict = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
        'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
        'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
        'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
        'HIE': 'H', 'HIP': 'H', 'HID': 'H', 'CYM': 'C', 'CYX': 'C',
        'ASH': 'D', 'GLH': 'E'
    }
    atom_dict = {}
    with open(pdb_file, "r") as f:
        lines = f.readlines()
    for line in lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            atom_serial = int(line[6:11].strip())
            res_name = line[17:20].strip()
            res_id = line[22:26].strip()
            one_letter = aa_dict.get(res_name, res_name)
            atom_dict[atom_serial] = (f"{res_id}{one_letter}")
    return atom_dict

def analysis_hbond_index(index_file,pdb_file):
    hbond = False
    hb_id = 0
    hb_dic = {}
    atom_dict = analysis_pdb(pdb_file)
    #print(atom_dict)
    with open(index_file,"r") as f:
        lines = f.readlines()
    for line in lines:
        if "hbonds_" in line:
            hbond = True
            continue
        if hbond:
            line = clean_string(line)
            Residue1 = atom_dict[int(line[0])]
            Residue2 = atom_dict[int(line[-1])]
            #print(f"{Residue1} vs. {Residue2}")
            hb_dic[hb_id] = f"{Residue1} vs. {Residue2}"
            hb_id += 1
    return hb_dic

def plot_hydrogen_occupation(labels,values,prefix,Title,height_range):
    #print(height_range)
    plt.figure(figsize=(10, max(6, len(labels) * 0.3)))
    bars = plt.barh(labels, values, color='skyblue')
    plt.xlabel("Occupancy (%)", fontsize=12)
    plt.ylabel("Hydrogen Bond", fontsize=12)
    plt.title(Title, fontsize=14, weight='bold')
    plt.gca().invert_yaxis()
    plt.grid(axis='x', linestyle='--', alpha=0.5)
    # 设置红色label
    ax = plt.gca()
    yticks = ax.get_yticklabels()
    for i, label_text in enumerate(labels):
        # 提取两个残基编号（假设格式是 "123Q vs. 456E"）
        try:
            res1, res2 = label_text.strip().split(" vs. ")
            id1 = int(''.join(filter(str.isdigit, res1)))
            id2 = int(''.join(filter(str.isdigit, res2)))
            if id1 in height_range or id2 in height_range:
                yticks[i].set_color('red')
        except:
            try:
                res = int(label_text[:3])
                #print(res)
                if res in height_range:
                    yticks[i].set_color('red')
            except:
                pass  # 如果格式不对就跳过
    # 添加右侧占有率文字
    for bar, value in zip(bars, values):
        plt.text(bar.get_width() + 1, bar.get_y() + bar.get_height()/2, f"{value:.1f}%", 
                 va='center', fontsize=10)
    plt.tight_layout()
    plt.savefig(f'{prefix}.png', dpi=300, bbox_inches='tight')
    #plt.show()
    plt.close()

def plot_hydrogen_distribution(matrix,hb_dic,x_label,y_label,prefix,Title,height_range):
    total_hbonds_per_frame = np.sum(matrix, axis=0)  # shape = (num_frames,)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 6 + 0.2 * len(hb_dic)), gridspec_kw={'height_ratios': [40, 1]}, sharex=True)
    im = ax1.imshow(matrix, aspect='auto', cmap='Reds', interpolation='nearest', origin='lower')
    ax1.set_ylabel(y_label, fontsize=12, weight='bold')
    ax1.set_title(Title, fontsize=14, weight='bold')
    ax1.set_yticks(np.arange(len(hb_dic)))
    labels = [hb_dic[i] for i in range(len(hb_dic))]
    ax1.set_yticklabels(labels, fontsize=8)
    # 获取当前 ytick 对象，并设置红色 label
    yticks = ax1.get_yticklabels()
    for label, tick in zip(labels, yticks):
        try:
            # 尝试解析为残基对
            res1, res2 = label.strip().split(" vs. ")
            id1 = int(''.join(filter(str.isdigit, res1)))
            id2 = int(''.join(filter(str.isdigit, res2)))
            if id1 in height_range or id2 in height_range:
                tick.set_color('red')
        except:
            try:
                res = int(label[:3])
                if res in height_range:
                    tick.set_color('red')
            except:
                pass  # 忽略解析失败的 label
    # 子图：1行热图表示每帧氢键总数
    hb_total_matrix = np.expand_dims(total_hbonds_per_frame, axis=0)
    im2 = ax2.imshow(hb_total_matrix, aspect='auto', cmap='Greys', interpolation='nearest', origin='lower')
    ax2.set_yticks([0])
    ax2.set_yticklabels(["Hydrogen Bond\nCount"], fontsize=10, weight='bold')
    ax2.set_xlabel(x_label, fontsize=12, weight='bold')
    ax2.tick_params(axis='x', labelsize=10)
    plt.tight_layout()
    plt.savefig(f'{prefix}.png', dpi=300, bbox_inches='tight')
    #plt.show()
    plt.close()

def xpm2png(xpm_file,index_file,pdb_file):
    with open(xpm_file,"r") as f:
        lines = f.readlines()
    Matrix=False
    num = 0
    data = []
    for line in lines:
        if "x-label" in line:
            x_label = line.split("\"")[1]
            #print(x_label)
        if "y-label" in line:
            y_label = line.split("\"")[1]
            #print(y_label)
        if "\"" in line and not Matrix and "*" not in line:
            Matrix = True
            temp = line.split("\"")[1]
            temp = clean_string(temp)
            x_len = int(temp[0])
            y_len = int(temp[1])
            continue
        if Matrix and "*" not in line:
            line = line.split("\"")[1]
            lenth = len(line)
            if lenth == x_len:
               num += 1
               line = line.replace(" ","0")
               line = line.replace("o","1")
               data.append([int(c) for c in line])
    hb_dic = analysis_hbond_index(index_file,pdb_file)
    labels=[hb_dic[i] for i in range(len(hb_dic))]
    if num != y_len:
        print("Error during read xpm data. Exit")
        exit(1)
    matrix = np.flipud(np.array(data))
    # 打印每一对氢键的占有率
    occurrences = np.sum(matrix, axis=1)  # 每一行对应一个氢键的总出现次数
    total_frames = matrix.shape[1]
    #print("Hydrogen Bond Occupancy:")
    Occupation = {}
    Occupation_single_residue = {}
    for i, count in enumerate(occurrences):
        occupancy = 100 * count / total_frames
        index = hb_dic[i]
        index1 = index.split(" ")[0]
        index2 = index.split(" ")[-1]
        temp1 = f'{index1} vs. {index2}'
        temp2 = f'{index2} vs. {index1}'
        if temp1 not in Occupation.keys() and temp2 not in Occupation.keys():
            Occupation[temp1] = occupancy
        else:
            try:
                Occupation[temp1] = occupancy + Occupation[temp1]
            except:
                try:
                    Occupation[temp2] = occupancy + Occupation[temp2]
                except:
                    print("Warning: hydrogen bond index may error.")
        #print(f"{hb_dic[i]}: {occupancy:.1f}%")
        if index1 not in Occupation_single_residue.keys():
            Occupation_single_residue[index1] = occupancy
        else:
            Occupation_single_residue[index1] = occupancy + Occupation_single_residue[index1]
        if index2 not in Occupation_single_residue.keys():
            Occupation_single_residue[index2] = occupancy
        else:
            Occupation_single_residue[index2] = occupancy + Occupation_single_residue[index2]
    prefix = xpm_file.split(".xpm")[0]
    single_residue_sorted_items = sorted(Occupation_single_residue.items(), key=lambda x: x[1], reverse=True)
    pair_residues_sorted_items = sorted(Occupation.items(), key=lambda x: x[1], reverse=True)
    return single_residue_sorted_items, pair_residues_sorted_items, matrix, hb_dic, x_label, y_label

def main(prefix_list, Struc_file, summary_occupation, hydrogen_distribution, hydrogen_occupation, hydrogen_occupation_single, height_range):
    # 这里是你主要的逻辑部分
    if summary_occupation:
        data_dict = {}
    # 在这里写你实际要执行的操作
    for prefix in prefix_list:
        single_residue_sorted_items, pair_residues_sorted_items, matrix, hb_dic, x_label, y_label = xpm2png(f"{prefix}.xpm",f"{prefix}.ndx",Struc_file)
        if summary_occupation:
            data_dict[prefix] = single_residue_sorted_items
        if hydrogen_distribution:
            Title = f"{prefix} hydrogen bond\nexistence map"
            prefix1 = f"{prefix}_distribution"
            plot_hydrogen_distribution(matrix,hb_dic,x_label,y_label,prefix1,Title,height_range)
            print(f"Hydrogen bond distribution over time ({prefix1}.png) generated successfully.")
        if hydrogen_occupation:
            labels = [item[0] for item in pair_residues_sorted_items]
            values = [item[1] for item in pair_residues_sorted_items]
            Title = f"{prefix} hydrogen bond occupation"
            prefix2 = f"{prefix}_hydrogen_occupation"
            plot_hydrogen_occupation(labels, values, prefix2, Title,height_range)
            print(f"Hydrogen bond occupation ratio for each residue pair ({prefix2}.png) generated successfully.")
        if hydrogen_occupation_single:
            labels = [item[0] for item in single_residue_sorted_items]
            values = [item[1] for item in single_residue_sorted_items]
            Title = f"{prefix} single residue\nhydrogen bond occupation"
            prefix3 = f"{prefix}_single_residue_hydrogen_occupation"
            plot_hydrogen_occupation(labels, values, prefix3, Title,height_range)
            print(f"Hydrogen bond occupation ratio for individual residue ({prefix3}.png) generated successfully.")
    if summary_occupation:
        # 构建 DataFrame
        df = pd.DataFrame()
        for label, pairs in data_dict.items():
            series = pd.Series(dict(pairs), name=label)
            df = pd.concat([df, series], axis=1)
        # 替换 NaN 为 0
        df = df.fillna(0)
        # 按残基编号排序（提取数字部分）
        def sort_key(x):
            return int(''.join(filter(str.isdigit, x)))
        df = df.sort_index(key=lambda x: [sort_key(i) for i in x])
        # 绘图
        plt.figure(figsize=(10, max(6, 0.3 * len(df))))
        sns.heatmap(df, cmap="Reds", linewidths=0.5, linecolor='gray', annot=True, fmt=".1f")
        plt.title("Residue Occupancy Across Multiple replications", fontsize=14, weight='bold')
        plt.ylabel("Residue", fontsize=12)
        plt.tight_layout()
        plt.savefig(f'heatmap.png', dpi=300, bbox_inches='tight')
        plt.close()
        print("Hydrogen bond distribution heatmap (heatmap.png) generated successfully.")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "该脚本用于 gmx hbond 的后续分析。\n"
            "通常情况下使用如下 gmx hbond 命令：\n"
            "    gmx hbond -s Structure.tpr -f trajectory.xtc -n index.ndx "
            "-num prefix.xvg -hbm prefix.xpm -hbn prefix.ndx\n"
            "将产生的 prefix.xpm 与 prefix.ndx 放置在同一目录下，通过下列选项进行分析。\n"
            "This script is used for post-analysis of gmx hbond output.\n"
            "Usually run gmx hbond as:\n"
            "    gmx hbond -s Structure.tpr -f trajectory.xtc -n index.ndx "
            "-num prefix.xvg -hbm prefix.xpm -hbn prefix.ndx\n"
            "Make sure prefix.xpm and prefix.ndx are in the same directory for analysis."
        ),
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        "--prefixes",
        type=str,
        help="如果有多个待分析的轨迹，可以使用逗号分割的字符串输入前缀，例如 \"a,b,c\"。\n"
             "If you have multiple trajectories to analyze, use comma-separated prefixes, e.g., \"a,b,c\"."
    )

    parser.add_argument(
        "--Struc_file",
        type=str,
        help="用于判断氢键残基类型的参考结构文件 <.pdb>，应与分析氢键时的 tpr 文件一致。\n"
             "Reference structure file <.pdb> used to determine residue types, must match the tpr used in gmx hbond."
    )

    parser.add_argument(
        "--summary_occupation",
        action="store_true",
        help="用热图展示不同轨迹中单个氨基酸的氢键占有率，大于100%%表示每帧可能出现多个氢键。\n"
             "Generate a heatmap of single residue H-bond occupation across trajectories (>100%% means >1 bond per frame). "
             "<heatmap.png>"
    )

    parser.add_argument(
        "--hydrogen_distribution",
        action="store_true",
        help="统计每对残基之间氢键在时间上的分布情况。\n"
             "Plot time-resolved H-bond presence for each residue pair. <prefix_distribution.png>"
    )

    parser.add_argument(
        "--hydrogen_occupation",
        action="store_true",
        help="统计每对氢键在轨迹中的占有率。\n"
             "Calculate total occupation ratio of each H-bond over the trajectory. <prefix_hydrogen_occupation.png>"
    )

    parser.add_argument(
        "--hydrogen_occupation_single",
        action="store_true",
        help="根据单个氨基酸统计其氢键占有率。\n"
             "Compute occupation ratios of H-bonds grouped by individual residues. <prefix_single_residue_hydrogen_occupation.png>"
    )

    parser.add_argument(
        "--height_range",
        type=str,
        help="如果想高亮标红某些氨基酸编号，可以逗号分隔输入，例如 \"123,456\"。\n"
             "To highlight certain residues on Y-axis (in red), provide comma-separated residue numbers, e.g., \"123,456\"."
    )
    args = parser.parse_args()
    prefix_list = args.prefixes.split(",")
    temp = ""
    for prefix in prefix_list:
        temp = f"{prefix} {temp}".strip()
        if not os.path.exists(f"{prefix}.xpm") or not os.path.exists(f"{prefix}.ndx"):
            print(f"Error: {prefix}.xpm or {prefix}.ndx does not exist.")
            sys.exit(1)
    height_range = [int(i) for i in args.height_range.split(",")] if args.height_range else []
    print(f"height_range: {height_range}.")
    print(f"Your systems are:{temp}.")
    print(f"reference structure: {args.Struc_file}.")
    main(prefix_list, args.Struc_file, args.summary_occupation, args.hydrogen_distribution, args.hydrogen_occupation, args.hydrogen_occupation_single, height_range)
