
#===============================================================================  
用途： 自动导入蛋白-配体系统的GROMACS轨迹，并优化显示方式。  
安装方法： 将此python文件放置到轨迹所在目录。  
使用方法：  
① 打开python文件。   
② 修改xtc_file, pdb_file, model_ID, Show_Residues_Near_MOL为对应数值，并运行代码。  
③ 打开vmd → 菜单栏 → Extensions → Tk Console，并在面板中输入：  
% source VMD.tcl (%不用输入)  
最后更新时间：2025年2月28日  
#===============================================================================  
Purpose: Automatically imports GROMACS trajectories of protein-ligand systems and and refines the display style.  
Installation Method: Place this Python file in the directory containing the trajectory files.  
Usage:  
① Open the Python file.  
② Modify xtc_file, pdb_file, model_ID, and Show_Residues_Near_MOL with the corresponding values, and run it.  
③ VMD main interface → Extensions → Tk Console, and input the following command in the panel:  
% source VMD.tcl (No need to input the % symbol)  
Last update: 2025-02-28  
