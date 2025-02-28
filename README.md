# 1. add_time_labels:  
  Automatically set time label according to frames, suitable for loading multiple models simultaneously.  
  自动根据帧设置时间戳，适用于同时载入多条model的情况。  
# 2. pbc_revise:   
  Revise instances where molecules jump in the trajectory due to PBC handling failures by GROMACS (by removing jump frames).  
  修正由于GROMACS周期性边界条件处理失败导致分子在轨迹中跳跃的情况（通过移除跳跃帧）。  