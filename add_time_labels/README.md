
#=======================================================================  
用途： 自动根据帧设置时间戳，适用于同时载入多条model的情况。  
安装方法： 将此文件放置到VMD安装目录，并在vmd.rc文件中设置如下命令：  
set add_time_labels "Your path to VMD/add_time_labels.tcl" (引号需要输入)  
使用方法：  
① 打开vmd → 导入轨迹   
② 菜单栏 → Extensions → Tk Console，并在面板中输入：  
% source $add_time_labels (%不用输入)  
在VMD main界面中会出现"Time Lable"，将时间调至左下角，  
双击F固定，并按T or R调整其余蛋白到合适位置  
注意调整每帧代表的时间长度  
最后更新时间：2024年08月29日  
#=======================================================================  
Purpose: Automatically set time label according to frames, suitable for loading multiple models simultaneously.  
Installation Method: Place this file in the VMD installation directory and set the following command in the vmd.rc file:  
set add_time_labels "Your path to VMD/add_time_labels.tcl" (quotation marks are required)  
Usage:  
① Open VMD → Import trajectory  
② VMD main interface → Extensions → Tk Console, and input the following command in the panel:  
% source $add_time_labels (No need to input the % symbol)  
"Time Label" will appear in the VMD main interface, adjust the time to the lower-left corner,  
double-click F to fix it, and press T or R to adjust other proteins to appropriate positions.  
Note to adjust the time length represented by each frame.  
Last update: 2024-08-29  
