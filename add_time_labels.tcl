# =======================================================================
# 用途： 自动根据帧设置时间戳，适用于同时载入多条model的情况。
# 安装方法： 将此文件放置到VMD安装目录，并在vmd.rc文件中设置如下命令：
# set add_time_labels "Your path to VMD/add_time_labels.tcl" (引号需要输入)
# 使用方法：
# ① 打开vmd → 导入轨迹 
# ② 菜单栏 → Extensions → Tk Console，并在面板中输入：
# % source $add_time_labels (%不用输入)
# 在VMD main界面中会出现"Time Lable"，将时间调至左下角，
# 双击F固定，并按T or R调整其余视角到合适位置
# =======================================================================

# 最后更新时间：2024年8月29日

# 每帧的时间跨度，单位在下面修改，默认为 1 ns
set md_time_step 1

color Display Background white ;# 设置背景为白色
axes location Off ;# 取消坐标轴

# 寻找时间最长的model
set longest_frames 0
set num_models [molinfo num]
for {set i 0} {$i < $num_models} {incr i} {
    set frames [molinfo $i get numframes]
    if {$frames > $longest_frames} {
        set longest_frames $frames
        set longest_model $i
    }
}
# 绘制时间标签的函数
proc draw_time_label {args} {
    global vmd_frame
    global md_time_step
    global longest_model
    global num_models
    graphics $num_models delete all
    graphics $num_models color black ;# 调节字体颜色
    graphics $num_models text {0 0 0} "[expr $vmd_frame($longest_model) * $md_time_step]ns" size 6 ;# 标签内容位置及字体大小，可以在这里修改单位ns为ps等
}

trace add variable vmd_frame($longest_model) write draw_time_label
mol new
mol rename $num_models "Time Label"