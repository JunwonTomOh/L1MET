# Project
open_project -reset "MET_Sim"

set_top Get_LUT
set cflags "-std=c++14"

add_files -tb tb_LUT.cpp -cflags "${cflags}"
add_files -tb ./firmware/puppimet.cpp -cflags "${cflags}"

# Solution
open_solution -reset "solution"
set_part {xcvu13p-flga2577-2-e}
create_clock -period 3.0 -name default

# csim_design
csim_design -ldflags "-B /usr/lib/x86_64-linux-gnu"

exit
