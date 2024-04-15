set ddr3dir ../../../dram/vc707/
set libdir ../../../../bluelib/src/coregen/vc707/

############# DDR3 Stuff
read_ip $ddr3dir/core/ddr3_0/ddr3_0.xci
read_verilog [ glob $ddr3dir/*.v ]
read_xdc $ddr3dir/dram.xdc
############# end Flash Stuff

############# Floating Point Stuff
source $libdir/../fp_import.tcl
############# Floating Point Stuff
