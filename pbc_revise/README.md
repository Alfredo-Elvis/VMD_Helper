
#=======================================================================  
用途： 修正由于GROMACS周期性边界条件处理失败导致分子在轨迹中跳跃的情况。（通过移除跳跃帧）  
安装方法： 将此文件放置到VMD安装目录，并在vmd.rc文件中设置如下命令：  
source <Your path to VMD>/pbc_revise.tcl  
使用方法：  
① 打开vmd → 导入轨迹   
② 菜单栏 → Extensions → Tk Console，并在面板中输入：  
%  pbc_revise -h (%不用输入)  
根据帮助信息使用  
最后更新时间：2025年2月28日  
#=======================================================================  
Purpose: Revise instances where molecules jump in the trajectory due to PBC handling failures by GROMACS (by removing jump frames).  
Installation Method: Place this file in the VMD installation directory and set the following command in the vmd.rc file:  
source <Your path to VMD>/pbc_revise.tcl  
Usage:  
① Open VMD → Import the trajectory  
② VMD main interface → Extensions → Tk Console, and enter the following command in the panel:  
% pbc_revise -h  (No need to input the % symbol)  
You can use the function as described in the help.  
Last update: 2025-02-28  
