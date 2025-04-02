# ===================================================================================
# Purpose: Fix molecule jumping in trajectories caused by GROMACS PBC handling failure
#          by removing frames with large deviations.
# Logic  : Use frame 0 as the reference, and compare each subsequent frame with
#          the last retained frame based on RMSD.
#          If RMSD exceeds the threshold, the frame is skipped; otherwise, it is kept.
# Installation:
# Place this file in the VMD installation directory and add the following line
# to your vmd.rc file:
# source pbc_revise.tcl
# Usage:
# ① Open VMD → Load your trajectory
# ② From the menu bar: Extensions → Tk Console, then enter:
# % pbc_revise -h (omit the '%' symbol)
# Follow the help instructions for usage.
# ===================================================================================

# ===================================================================================
# Purpose: Revise instances where molecules jump in the trajectory due to PBC handling failures by GROMACS (by removing jump frames).
# Installation Method: Place this file in the VMD installation directory and set the following command in the vmd.rc file:
# source pbc_revise.tcl
# Usage:
# ① Open VMD → Import the trajectory
# ② VMD main interface → Extensions → Tk Console, and enter the following command in the panel:
# % pbc_revise -h  (No need to input the % symbol)
# You can use the function as described in the help.
# Last update: 2025-03-26
# ===================================================================================

proc pbc_revise {args} {
    # Default parameters
    set rmsd_threshold 10
    set rmsd_cal_region ""
    set out_file ""
    set window_size 5

    # Argument parsing
    set i 0
    while {$i < [llength $args]} {
        set arg [lindex $args $i]
        switch -- $arg {
            "-rmsd_threshold" {
                incr i
                set rmsd_threshold [lindex $args $i]
            }
            "-rmsd_cal_region" {
                incr i
                set rmsd_cal_region [lindex $args $i]
            }
            "-out_file" {
                incr i
                set out_file [lindex $args $i]
            }
            "-window_size" {
                incr i
                set window_size [lindex $args $i]
            }
            "-h" {
                puts "Usage: pbc_revise -rmsd_threshold <value> -rmsd_cal_region <atom selection> ?-window_size <int>? ?-out_file <filename>?"
                puts "Example: pbc_revise -rmsd_threshold 10 -rmsd_cal_region \"resname MOL and noh\" -window_size 5"
                return
            }
            default {
                puts "Unknown option: $arg. Use '-h' for help."
                return
            }
        }
        incr i
    }

    if {$rmsd_cal_region eq ""} {
        puts "Error: -rmsd_cal_region must be specified."
        return
    }

    if {$out_file eq ""} {
        set timestamp [clock format [clock seconds] -format "%Y%m%d_%H%M%S"]
        set out_file "smart_window_cleaned_$timestamp.trr"
    }

    if {[catch {molinfo top} mol_id]} {
        puts "Error: No molecule loaded. Please load a trajectory first."
        return
    }

    set total_frames [molinfo $mol_id get numframes]
    set natoms [molinfo $mol_id get numatoms]

    puts "Smart Window PBC Frame Revision"
    puts "-------------------------------"
    puts "Model ID          : $mol_id"
    puts "Total frames      : $total_frames"
    puts "RMSD threshold    : $rmsd_threshold"
    puts "Window size       : $window_size"
    puts "Atom region       : $rmsd_cal_region"
    puts "Output file       : $out_file"

    # Initialize progress tracker
    set progress_step [expr {int($total_frames * 0.05)}]
    if {$progress_step < 1} { set progress_step 1 }
    set next_progress $progress_step

    set selected_frames {0}
    set ref_index 0
    set i 1

    while {$i < $total_frames} {
        if {$i >= $next_progress} {
            set percent [expr {int(100.0 * $i / $total_frames)}]
            puts "Progress: $percent% completed..."
            incr next_progress $progress_step
        }

        set ref_sel [atomselect $mol_id "$rmsd_cal_region" frame $ref_index]
        set cur_sel [atomselect $mol_id "$rmsd_cal_region" frame $i]
        set rmsd [measure rmsd $ref_sel $cur_sel]
        $ref_sel delete
        $cur_sel delete

        if {$rmsd < $rmsd_threshold} {
            lappend selected_frames $i
            set ref_index $i
            incr i
        } else {
            # Jump detected at frame i → window search
            set jump_start $i
            set jump_end -1
            set j [expr {$i + 1}]
            while {$j < $total_frames && $j <= $i + $window_size} {
                set prev_sel [atomselect $mol_id "$rmsd_cal_region" frame [expr {$j - 1}]]
                set test_sel [atomselect $mol_id "$rmsd_cal_region" frame $j]
                set rmsd_j [measure rmsd $prev_sel $test_sel]
                $prev_sel delete
                $test_sel delete

                if {$rmsd_j > $rmsd_threshold} {
                    set jump_end $j
                    break
                }
                incr j
            }

            if {$jump_end != -1} {
                #puts "Jump confirmed: removing frames $jump_start to [expr {$jump_end - 1}]"
                set i [expr {$jump_end + 1}]
                # Reference frame remains unchanged
            } else {
                # No second jump → assume A is okay, keep A ~ j-1
                #puts "Jump tolerated: keeping frames $jump_start to [expr {$j - 1}]"
                for {set k $jump_start} {$k < $j} {incr k} {
                    lappend selected_frames $k
                }
                set ref_index [expr {$j - 1}]
                set i $j
            }
        }
    }

    # Build new trajectory
    puts "Creating new molecule and copying selected frames..."
    set out_mol [mol new atoms $natoms]
    set all_sel_orig [atomselect $mol_id all]
    set all_sel_new [atomselect $out_mol all]

    foreach frame $selected_frames {
        $all_sel_orig frame $frame
        $all_sel_new frame [molinfo $out_mol get numframes]
        $all_sel_new set {x y z} [$all_sel_orig get {x y z}]
        animate dup $out_mol
    }

    $all_sel_orig delete
    $all_sel_new delete

    puts "Writing [llength $selected_frames] frames to $out_file ..."
    animate write trr $out_file beg 0 end -1 $out_mol

    puts "Done. Total selected frames: [llength $selected_frames]"
}