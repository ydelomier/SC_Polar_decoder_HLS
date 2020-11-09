############################################################
## This file is generated automatically by Vivado HLS.
## Please DO NOT edit it.
## Copyright (C) 1986-2018 Xilinx, Inc. All Rights Reserved.
############################################################
open_project SC_RTL_Simulation
set_top my_module
add_files src/module/my_module.cpp
add_files src/module/my_module.h
add_files -tb src/module/config.h
add_files -tb ../shared/src/functions.h
add_files -tb ../shared/src/library.h
add_files -tb src/rtl_simu_testbench/main.cpp
add_files -tb src/module/polar_parameters.h
add_files -tb src/rtl_simu_testbench/sc_generator/sc_generator.h
add_files -tb src/rtl_simu_testbench/sc_monitor/sc_monitor.h
add_files -tb src/rtl_simu_testbench/sc_top_module.h
add_files -tb ../shared/src/scalar.h
add_files -tb ../shared/src/vector.h

open_solution "solution_v7_100MHz"
set_part {xc7vx690tffg1761-2} -tool vivado
create_clock -period 10 -name default
config_compile -no_signed_zeros=0 -unsafe_math_optimizations=0
set_clock_uncertainty 1
#source "./SC_RTL_Simulation/solution1/directives.tcl"

#csim_design
#csynth_design
#cosim_design -rtl vhdl -tool xsim
export_design -flow impl -rtl vhdl -format ip_catalog

exit