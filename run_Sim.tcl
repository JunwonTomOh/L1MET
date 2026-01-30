array set opt {
  nevents    10
  write      1
  read       0
  debug      1
  verbose    1
}

foreach arg $::argv {
  foreach o [lsort [array names opt]] {
    regexp "$o=+(\\w+)" $arg unused opt($o)
  }
}

# Configuration
set cflags "-std=c++14"
if {$opt(nevents)} {
    append cflags " -DNEVENTS=$opt(nevents)"
}
if {$opt(write)} {
    append cflags " -DMET_WRITE_TB_FILE"
}
if {$opt(read)} {
    append cflags " -DMET_READ_TB_FILE"
}
if {$opt(debug)} {
    append cflags " -DMET_DEBUG"
}

if {$opt(verbose)} {
    append cflags " -DMET_VERBOSE"
}

# Project
open_project -reset "MET_Sim"

set_top puppimet_xy

add_files -tb tb_puppimet.cpp -cflags "${cflags}"
add_files -tb ./firmware/puppimet.cpp -cflags "${cflags}"

add_files -tb ./data/ParticlesIn.txt
add_files -tb ./data/DeregionizerIn.txt
add_files -tb ./data/METsOut.txt

# Solution
open_solution -reset "solution"
set_part {xcvu13p-flga2577-2-e}
create_clock -period 3.0 -name default

# csim_design
csim_design -ldflags "-B /usr/lib/x86_64-linux-gnu"

exit
