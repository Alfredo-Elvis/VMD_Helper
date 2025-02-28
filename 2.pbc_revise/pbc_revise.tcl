# ===================================================================================
# 用途： 修正由于GROMACS周期性边界条件处理失败导致分子在轨迹中跳跃的情况（通过移除跳跃帧）。
# 安装方法： 将此文件放置到VMD安装目录，并在vmd.rc文件中设置如下命令：
# source <Your path to VMD>/pbc_revise.tcl
# 使用方法：
# ① 打开vmd → 导入轨迹 
# ② 菜单栏 → Extensions → Tk Console，并在面板中输入：
# %  pbc_revise -h (%不用输入)
# 根据帮助信息使用
# 最后更新时间：2025年2月28日
# ===================================================================================

# ===================================================================================
# Purpose: Revise instances where molecules jump in the trajectory due to PBC handling failures by GROMACS (by removing jump frames).
# Installation Method: Place this file in the VMD installation directory and set the following command in the vmd.rc file:
# source <Your path to VMD>/pbc_revise.tcl
# Usage:
# ① Open VMD → Import the trajectory
# ② VMD main interface → Extensions → Tk Console, and enter the following command in the panel:
# % pbc_revise -h  (No need to input the % symbol)
# You can use the function as described in the help.
# Last update: 2025-02-28
# ===================================================================================

proc pbc_revise {args} {
    if {[lsearch -exact $args "-h"] != -1} {
        puts "Info: "
        puts "Select frames from the last input model where the RMSD is below a specified threshold."
        puts "Usage: pbc_revise -rmsd_threshold value -reference_frame value -rmsd_cal_region value"
        puts "Arguments:"
        puts "  -rmsd_threshold (10): The threshold value for RMSD filtering."
        puts "  -reference_frame(0) : The frame index used as the reference frame for RMSD calculation."
        puts "  -rmsd_cal_region    : The selection string to specify which atoms to calculate RMSD."
        puts "Example:"
        puts "  pbc_revise -rmsd_threshold 20 -reference_frame 10000 -rmsd_cal_region 'resname MOL and noh'"
        puts "The output file is: selected_frames.trr"
        return
    }
    # 初始化默认参数
    set rmsd_threshold 10
    set reference_frame 0
    if {[llength $args] == 0} {
        puts "Please set rmsd_cal_region, use '-h' for help."
        return
    }
    # 解析传入的参数
    set i 0
    while {$i < [llength $args]} {
        set arg [lindex $args $i]
        if {$arg == "-rmsd_threshold"} {
            incr i
            set rmsd_threshold [lindex $args $i]
        } elseif {$arg == "-reference_frame"} {
            incr i
            set reference_frame [lindex $args $i]
        } elseif {$arg == "-rmsd_cal_region"} {
            incr i
            set rmsd_cal_region [lindex $args $i]
        } else {
            puts "Unknown option: $arg. Use '-h' for help."
            return
        }
        incr i
    }
    # 输出参数确认
    puts "RMSD threshold: $rmsd_threshold"
    puts "Reference frame: $reference_frame"
    puts "RMSD calculation region: $rmsd_cal_region"
    set out_file "selected_frames.trr"
    set num_models [molinfo list]
    set mol_id [expr {[llength $num_models] - 1}]
    set out_id [expr [llength $num_models]]
    set selection [atomselect $mol_id $rmsd_cal_region]
    set selected_frames {}
    set total_frames [molinfo $mol_id get numframes]
    animate write pdb $out_file beg 0 end 0 $mol_id
    mol new $out_file type pdb
    # 遍历每一帧
    for {set i 0} {$i < $total_frames} {incr i} {
        # 获取当前帧的 RMSD，相对于参考帧
        # 第一个 atomselect 是参考帧上的原子，第二个 atomselect 是当前帧上的原子
        set rmsd [measure rmsd [atomselect $mol_id $rmsd_cal_region frame $reference_frame] \
                          [atomselect $mol_id $rmsd_cal_region frame $i]]
        # 如果 RMSD 小于设定阈值，加入到帧列表中
        if {$rmsd < $rmsd_threshold} {
            lappend selected_frames $i
            #puts "RMSD:$rmsd, frame:$i"
            break
        }
    }
    # 遍历选中的帧
    foreach frame $selected_frames {
        if {$first_frame} {
            animate write trr $out_file beg $frame end $frame $mol_id
            mol addfile $out_file type trr first 0 last -1 step 1 waitfor -1 $out_id
        }
    }
    #写出结果
    animate write trr $out_file beg 0 end -1 $out_id
    puts "Total select frames: [llength $selected_frames]"
}