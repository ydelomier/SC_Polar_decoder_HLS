#!/bin/bash
set -x

declare -a NBIT
declare -a KBIT
declare -a QUANT

NBIT=(65536 32768 32768 32768)
KBIT=(58982 29492 29492 29492)
QUANT=(6 9 7 8)

for ELAG in 1
do

	for REF in 0 1 2 3
	do

		for PAR in 8 16
		do

			# Generate Proper parameters file
			cd /c/Polar_decoder_HLS/Polar_Code_Generator/Frozen_Bit_Generator/cmake-build-debug/   
			./FB_Generator.exe ${NBIT[$REF]} ${KBIT[$REF]} $PAR 0 ../../Generated_Frozen_Bit/frozen_n_${NBIT[$REF]}_k_${KBIT[$REF]}.txt 1 /c/Polar_decoder_HLS/SC_Polar_decoder_v5_Full_Elag/src/module/
					
			cd /c/Polar_decoder_HLS/SC_Polar_decoder_v5_Full_Elag/ 
					
			# Generate Proper config file
			echo " 
#define LLR_BITS		${QUANT[$REF]} 
#define LLR_MAXV    	(+127)
#define LLR_MINV    	(-127)

#define LLR 			sc_bigint<LLR_BITS>
#define BIT 			sc_biguint<1>
#define TYPE_LLRS   	sc_bigint<PAR * LLR_BITS>
#define TYPE_BITS   	sc_biguint<PAR>

#define SIGMAG 				// CA2 or SIGMAG , define the representation format
						// CA2 : 2's complement , SIGMAG : sign and magnitude

#define EXTENDED		1	// 1: on, 0: off. Enable the extension of the quantization format for the internal signals of the operators

#define PRUNING_LEVEL 	$ELAG 	// 0: No Pruning, 1: Leaf Pruning, 2: Full Pruning enable (R0 node pruned by default when enable)

// Choose the type of node to pruned. 1 : enable, 0 : off
#define ELAG_R1			1	// Node of Rate 1 						   , by default only for G_function
#define ELAG_REP		0	// Repetition Node 			(ex: 00000001) , by default only for F_function
#define ELAG_SPC		0	// Single Parity Check Node (ex: 01111111) , by default only for G_function
#define ELAG_REP2		0	// REP2 Node 				(ex: 00000011) , by default only for F_function
#define ELAG_SPC2		0	// SPC2 Node 				(ex: 00111111) , by default only for G_function

// These options enable pruning for the rariest cases  ( F_R1, G_R0, F_SPC, ...)
#define ELAG_RARE		0	// 1 enable, 0 off

#define ELAG_H0			0	// When pruning a node R0 from a F function, the F_R0 state can be bypass at the cost of a dedicated H function : H0

#define _MONITORING_ " > src/module/config.h
					
			# Launch Vivado_HLS Script
			vivado_hls -f script_Syn_Lat.tcl 
			vivado_hls -f script_Syn_Lat2.tcl 
			vivado_hls -f script_Syn_Lat3.tcl 

			# Copy the csim log
			mkdir -p Results_Implem_Latency/Comparaison_with_literature/Comp_BLG_Ref_$REF-N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-Elag_$ELAG/PAR_$PAR/csim 
			cp SC_RTL_Simulation/solution_v7_100MHz/csim/report/my_module_csim.log Results_Implem_Latency/Comparaison_with_literature/Comp_BLG_Ref_$REF-N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-Elag_$ELAG/PAR_$PAR/csim 

			# Copy the synthesis report files for Solution 100MHz
			mkdir -p Results_Implem_Latency/Comparaison_with_literature/Comp_BLG_Ref_$REF-N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-Elag_$ELAG/PAR_$PAR/Solution_100MHz/syn 
			cp SC_RTL_Simulation/solution_v7_100MHz/syn/report/*.rpt Results_Implem_Latency/Comparaison_with_literature/Comp_BLG_Ref_$REF-N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-Elag_$ELAG/PAR_$PAR/Solution_100MHz/syn 

			# Copy the vhdl files for Solution 100MHz
			cp -R SC_RTL_Simulation/solution_v7_100MHz/syn/vhdl Results_Implem_Latency/Comparaison_with_literature/Comp_BLG_Ref_$REF-N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-Elag_$ELAG/PAR_$PAR/Solution_100MHz/syn 

			# Copy the cosim log for Solution 100MHz
			mkdir -p Results_Implem_Latency/Comparaison_with_literature/Comp_BLG_Ref_$REF-N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-Elag_$ELAG/PAR_$PAR/Solution_100MHz/cosim 
			cp SC_RTL_Simulation/solution_v7_100MHz/sim/report/vhdl/my_module.log Results_Implem_Latency/Comparaison_with_literature/Comp_BLG_Ref_$REF-N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-Elag_$ELAG/PAR_$PAR/Solution_100MHz/cosim 

			# Copy the implementation report for Solution 100MHz
			mkdir -p Results_Implem_Latency/Comparaison_with_literature/Comp_BLG_Ref_$REF-N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-Elag_$ELAG/PAR_$PAR/Solution_100MHz/impl 
			cp SC_RTL_Simulation/solution_v7_100MHz/impl/report/vhdl/my_module_export.rpt Results_Implem_Latency/Comparaison_with_literature/Comp_BLG_Ref_$REF-N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-Elag_$ELAG/PAR_$PAR/Solution_100MHz/impl 


			# Copy the synthesis report files for Solution 300MHz
			mkdir -p Results_Implem_Latency/Comparaison_with_literature/Comp_BLG_Ref_$REF-N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-Elag_$ELAG/PAR_$PAR/Solution_300MHz/syn 
			cp SC_RTL_Simulation/solution_v7_300MHz/syn/report/*.rpt Results_Implem_Latency/Comparaison_with_literature/Comp_BLG_Ref_$REF-N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-Elag_$ELAG/PAR_$PAR/Solution_300MHz/syn 

			# Copy the vhdl files for Solution 300MHz
			cp -R SC_RTL_Simulation/solution_v7_300MHz/syn/vhdl Results_Implem_Latency/Comparaison_with_literature/Comp_BLG_Ref_$REF-N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-Elag_$ELAG/PAR_$PAR/Solution_300MHz/syn 

			# Copy the cosim log for Solution 300MHz
			mkdir -p Results_Implem_Latency/Comparaison_with_literature/Comp_BLG_Ref_$REF-N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-Elag_$ELAG/PAR_$PAR/Solution_300MHz/cosim 
			cp SC_RTL_Simulation/solution_v7_300MHz/sim/report/vhdl/my_module.log Results_Implem_Latency/Comparaison_with_literature/Comp_BLG_Ref_$REF-N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-Elag_$ELAG/PAR_$PAR/Solution_300MHz/cosim 

			# Copy the implementation report for Solution 300MHz
			mkdir -p Results_Implem_Latency/Comparaison_with_literature/Comp_BLG_Ref_$REF-N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-Elag_$ELAG/PAR_$PAR/Solution_300MHz/impl 
			cp SC_RTL_Simulation/solution_v7_300MHz/impl/report/vhdl/my_module_export.rpt Results_Implem_Latency/Comparaison_with_literature/Comp_BLG_Ref_$REF-N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-Elag_$ELAG/PAR_$PAR/Solution_300MHz/impl 

			# Erase the Solution for new test
			cd SC_RTL_Simulation/ 
			rm -r solution_v7_100MHz/ 
			rm -r solution_v7_300MHz/ 
				
		done
	done
done


NBIT=(32768)
KBIT=(29492)
QUANT=(6)

############### FULL Elag ###################

for ELAG in 2
do

	for REF in 0
	do

		for PAR in 4 8 16 32 64
		do

			# Generate Proper parameters file
			cd /c/Polar_decoder_HLS/Polar_Code_Generator/Frozen_Bit_Generator/cmake-build-debug/   
			./FB_Generator.exe ${NBIT[$REF]} ${KBIT[$REF]} $PAR 0 ../../Generated_Frozen_Bit/frozen_n_${NBIT[$REF]}_k_${KBIT[$REF]}.txt 1 /c/Polar_decoder_HLS/SC_Polar_decoder_v5_Full_Elag/src/module/
					
			# Generate Proper config file
			cd /c/Polar_decoder_HLS/SC_Polar_decoder_v5_Full_Elag/ 
			echo " 
#define LLR_BITS		${QUANT[$REF]} 
#define LLR_MAXV    	(+127)
#define LLR_MINV    	(-127)

#define LLR 			sc_bigint<LLR_BITS>
#define BIT 			sc_biguint<1>
#define TYPE_LLRS   	sc_bigint<PAR * LLR_BITS>
#define TYPE_BITS   	sc_biguint<PAR>

#define CA2 				// CA2 or SIGMAG , define the representation format
							// CA2 : 2's complement , SIGMAG : sign and magnitude

#define EXTENDED		1	// 1: on, 0: off. Enable the extension of the quantization format for the internal signals of the operators

#define PRUNING_LEVEL 	$ELAG 	// 0: No Pruning, 1: Leaf Pruning, 2: Full Pruning enable (R0 node pruned by default when enable)

// Choose the type of node to pruned. 1 : enable, 0 : off
#define ELAG_R1			1	// Node of Rate 1 						   , by default only for G_function
#define ELAG_REP		1	// Repetition Node 			(ex: 00000001) , by default only for F_function
#define ELAG_SPC		1	// Single Parity Check Node (ex: 01111111) , by default only for G_function
#define ELAG_REP2		0	// REP2 Node 				(ex: 00000011) , by default only for F_function
#define ELAG_SPC2		0	// SPC2 Node 				(ex: 00111111) , by default only for G_function

// These options enable pruning for the rariest cases  ( F_R1, G_R0, F_SPC, ...)
#define ELAG_RARE		0	// 1 enable, 0 off

#define ELAG_H0			1	// When pruning a node R0 from a F function, the F_R0 state can be bypass at the cost of a dedicated H function : H0

#define _MONITORING_ " > src/module/config.h
					
			# Launch Vivado_HLS Script
			vivado_hls -f script_Syn_Lat.tcl 
			vivado_hls -f script_Syn_Lat2.tcl 
			vivado_hls -f script_Syn_Lat3.tcl 

			# Copy the csim log
			mkdir -p Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-FULL_Elag/PAR_$PAR/csim 
			cp SC_RTL_Simulation/solution_v7_100MHz/csim/report/my_module_csim.log Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-FULL_Elag/PAR_$PAR/csim 

			# Copy the synthesis report files for Solution 100MHz
			mkdir -p Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-FULL_Elag/PAR_$PAR/Solution_100MHz/syn 
			cp SC_RTL_Simulation/solution_v7_100MHz/syn/report/*.rpt Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-FULL_Elag/PAR_$PAR/Solution_100MHz/syn 

			# Copy the vhdl files for Solution 100MHz
			cp -R SC_RTL_Simulation/solution_v7_100MHz/syn/vhdl Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-FULL_Elag/PAR_$PAR/Solution_100MHz/syn 

			# Copy the cosim log for Solution 100MHz
			mkdir -p Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-FULL_Elag/PAR_$PAR/Solution_100MHz/cosim 
			cp SC_RTL_Simulation/solution_v7_100MHz/sim/report/vhdl/my_module.log Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-FULL_Elag/PAR_$PAR/Solution_100MHz/cosim 

			# Copy the implementation report for Solution 100MHz
			mkdir -p Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-FULL_Elag/PAR_$PAR/Solution_100MHz/impl 
			cp SC_RTL_Simulation/solution_v7_100MHz/impl/report/vhdl/my_module_export.rpt Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-FULL_Elag/PAR_$PAR/Solution_100MHz/impl 


			# Copy the synthesis report files for Solution 300MHz
			mkdir -p Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-FULL_Elag/PAR_$PAR/Solution_300MHz/syn 
			cp SC_RTL_Simulation/solution_v7_300MHz/syn/report/*.rpt Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-FULL_Elag/PAR_$PAR/Solution_300MHz/syn 

			# Copy the vhdl files for Solution 300MHz
			cp -R SC_RTL_Simulation/solution_v7_300MHz/syn/vhdl Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-FULL_Elag/PAR_$PAR/Solution_300MHz/syn 

			# Copy the cosim log for Solution 300MHz
			mkdir -p Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-FULL_Elag/PAR_$PAR/Solution_300MHz/cosim 
			cp SC_RTL_Simulation/solution_v7_300MHz/sim/report/vhdl/my_module.log Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-FULL_Elag/PAR_$PAR/Solution_300MHz/cosim 

			# Copy the implementation report for Solution 300MHz
			mkdir -p Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-FULL_Elag/PAR_$PAR/Solution_300MHz/impl 
			cp SC_RTL_Simulation/solution_v7_300MHz/impl/report/vhdl/my_module_export.rpt Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-FULL_Elag/PAR_$PAR/Solution_300MHz/impl 

			# Erase the Solution for new test
			cd SC_RTL_Simulation/ 
			rm -r solution_v7_100MHz/ 
			rm -r solution_v7_300MHz/ 
				
		done
	done
done

############### ELAG R0R1 ###################

for ELAG in 2
do

	for REF in 0
	do

		for PAR in 4 8 16 32 64
		do

			# Generate Proper parameters file
			cd /c/Polar_decoder_HLS/Polar_Code_Generator/Frozen_Bit_Generator/cmake-build-debug/   
			./FB_Generator.exe ${NBIT[$REF]} ${KBIT[$REF]} $PAR 0 ../../Generated_Frozen_Bit/frozen_n_${NBIT[$REF]}_k_${KBIT[$REF]}.txt 1 /c/Polar_decoder_HLS/SC_Polar_decoder_v5_Full_Elag/src/module/
					
			# Generate Proper config file
			cd /c/Polar_decoder_HLS/SC_Polar_decoder_v5_Full_Elag/ 
			echo " 
#define LLR_BITS		${QUANT[$REF]} 
#define LLR_MAXV    	(+127)
#define LLR_MINV    	(-127)

#define LLR 			sc_bigint<LLR_BITS>
#define BIT 			sc_biguint<1>
#define TYPE_LLRS   	sc_bigint<PAR * LLR_BITS>
#define TYPE_BITS   	sc_biguint<PAR>

#define CA2 				// CA2 or SIGMAG , define the representation format
							// CA2 : 2's complement , SIGMAG : sign and magnitude

#define EXTENDED		1	// 1: on, 0: off. Enable the extension of the quantization format for the internal signals of the operators

#define PRUNING_LEVEL 	$ELAG 	// 0: No Pruning, 1: Leaf Pruning, 2: Full Pruning enable (R0 node pruned by default when enable)

// Choose the type of node to pruned. 1 : enable, 0 : off
#define ELAG_R1			1	// Node of Rate 1 						   , by default only for G_function
#define ELAG_REP		0	// Repetition Node 			(ex: 00000001) , by default only for F_function
#define ELAG_SPC		0	// Single Parity Check Node (ex: 01111111) , by default only for G_function
#define ELAG_REP2		0	// REP2 Node 				(ex: 00000011) , by default only for F_function
#define ELAG_SPC2		0	// SPC2 Node 				(ex: 00111111) , by default only for G_function

// These options enable pruning for the rariest cases  ( F_R1, G_R0, F_SPC, ...)
#define ELAG_RARE		0	// 1 enable, 0 off

#define ELAG_H0			1	// When pruning a node R0 from a F function, the F_R0 state can be bypass at the cost of a dedicated H function : H0

#define _MONITORING_ " > src/module/config.h
					
			# Launch Vivado_HLS Script
			vivado_hls -f script_Syn_Lat.tcl 
			vivado_hls -f script_Syn_Lat2.tcl 
			vivado_hls -f script_Syn_Lat3.tcl 

			# Copy the csim log
			mkdir -p Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-Elag_R0R1/PAR_$PAR/csim 
			cp SC_RTL_Simulation/solution_v7_100MHz/csim/report/my_module_csim.log Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-Elag_R0R1/PAR_$PAR/csim 

			# Copy the synthesis report files for Solution 100MHz
			mkdir -p Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-Elag_R0R1/PAR_$PAR/Solution_100MHz/syn 
			cp SC_RTL_Simulation/solution_v7_100MHz/syn/report/*.rpt Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-Elag_R0R1/PAR_$PAR/Solution_100MHz/syn 

			# Copy the vhdl files for Solution 100MHz
			cp -R SC_RTL_Simulation/solution_v7_100MHz/syn/vhdl Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-Elag_R0R1/PAR_$PAR/Solution_100MHz/syn 

			# Copy the cosim log for Solution 100MHz
			mkdir -p Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-Elag_R0R1/PAR_$PAR/Solution_100MHz/cosim 
			cp SC_RTL_Simulation/solution_v7_100MHz/sim/report/vhdl/my_module.log Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-Elag_R0R1/PAR_$PAR/Solution_100MHz/cosim 

			# Copy the implementation report for Solution 100MHz
			mkdir -p Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-Elag_R0R1/PAR_$PAR/Solution_100MHz/impl 
			cp SC_RTL_Simulation/solution_v7_100MHz/impl/report/vhdl/my_module_export.rpt Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-Elag_R0R1/PAR_$PAR/Solution_100MHz/impl 


			# Copy the synthesis report files for Solution 300MHz
			mkdir -p Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-Elag_R0R1/PAR_$PAR/Solution_300MHz/syn 
			cp SC_RTL_Simulation/solution_v7_300MHz/syn/report/*.rpt Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-Elag_R0R1/PAR_$PAR/Solution_300MHz/syn 

			# Copy the vhdl files for Solution 300MHz
			cp -R SC_RTL_Simulation/solution_v7_300MHz/syn/vhdl Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-Elag_R0R1/PAR_$PAR/Solution_300MHz/syn 

			# Copy the cosim log for Solution 300MHz
			mkdir -p Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-Elag_R0R1/PAR_$PAR/Solution_300MHz/cosim 
			cp SC_RTL_Simulation/solution_v7_300MHz/sim/report/vhdl/my_module.log Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-Elag_R0R1/PAR_$PAR/Solution_300MHz/cosim 

			# Copy the implementation report for Solution 300MHz
			mkdir -p Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-Elag_R0R1/PAR_$PAR/Solution_300MHz/impl 
			cp SC_RTL_Simulation/solution_v7_300MHz/impl/report/vhdl/my_module_export.rpt Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-Elag_R0R1/PAR_$PAR/Solution_300MHz/impl 

			# Erase the Solution for new test
			cd SC_RTL_Simulation/ 
			rm -r solution_v7_100MHz/ 
			rm -r solution_v7_300MHz/ 
				
		done
	done
done

############### NO Elag ###################

for ELAG in 0
do

	for REF in 0
	do

		for PAR in 4 8 16 32 64
		do

			# Generate Proper parameters file
			cd /c/Polar_decoder_HLS/Polar_Code_Generator/Frozen_Bit_Generator/cmake-build-debug/   
			./FB_Generator.exe ${NBIT[$REF]} ${KBIT[$REF]} $PAR 0 ../../Generated_Frozen_Bit/frozen_n_${NBIT[$REF]}_k_${KBIT[$REF]}.txt 1 /c/Polar_decoder_HLS/SC_Polar_decoder_v5_Full_Elag/src/module/
					
			# Generate Proper config file
			cd /c/Polar_decoder_HLS/SC_Polar_decoder_v5_Full_Elag/ 
			echo " 
#define LLR_BITS		${QUANT[$REF]} 
#define LLR_MAXV    	(+127)
#define LLR_MINV    	(-127)

#define LLR 			sc_bigint<LLR_BITS>
#define BIT 			sc_biguint<1>
#define TYPE_LLRS   	sc_bigint<PAR * LLR_BITS>
#define TYPE_BITS   	sc_biguint<PAR>

#define CA2 				// CA2 or SIGMAG , define the representation format
							// CA2 : 2's complement , SIGMAG : sign and magnitude

#define EXTENDED		1	// 1: on, 0: off. Enable the extension of the quantization format for the internal signals of the operators

#define PRUNING_LEVEL 	$ELAG 	// 0: No Pruning, 1: Leaf Pruning, 2: Full Pruning enable (R0 node pruned by default when enable)

// Choose the type of node to pruned. 1 : enable, 0 : off
#define ELAG_R1			0	// Node of Rate 1 						   , by default only for G_function
#define ELAG_REP		0	// Repetition Node 			(ex: 00000001) , by default only for F_function
#define ELAG_SPC		0	// Single Parity Check Node (ex: 01111111) , by default only for G_function
#define ELAG_REP2		0	// REP2 Node 				(ex: 00000011) , by default only for F_function
#define ELAG_SPC2		0	// SPC2 Node 				(ex: 00111111) , by default only for G_function

// These options enable pruning for the rariest cases  ( F_R1, G_R0, F_SPC, ...)
#define ELAG_RARE		0	// 1 enable, 0 off

#define ELAG_H0			0	// When pruning a node R0 from a F function, the F_R0 state can be bypass at the cost of a dedicated H function : H0

#define _MONITORING_ " > src/module/config.h
					
			# Launch Vivado_HLS Script
			vivado_hls -f script_Syn_Lat.tcl 
			vivado_hls -f script_Syn_Lat2.tcl 
			vivado_hls -f script_Syn_Lat3.tcl 

			# Copy the csim log
			mkdir -p Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-NO_Elag/PAR_$PAR/csim 
			cp SC_RTL_Simulation/solution_v7_100MHz/csim/report/my_module_csim.log Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-NO_Elag/PAR_$PAR/csim 

			# Copy the synthesis report files for Solution 100MHz
			mkdir -p Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-NO_Elag/PAR_$PAR/Solution_100MHz/syn 
			cp SC_RTL_Simulation/solution_v7_100MHz/syn/report/*.rpt Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-NO_Elag/PAR_$PAR/Solution_100MHz/syn 

			# Copy the vhdl files for Solution 100MHz
			cp -R SC_RTL_Simulation/solution_v7_100MHz/syn/vhdl Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-NO_Elag/PAR_$PAR/Solution_100MHz/syn 

			# Copy the cosim log for Solution 100MHz
			mkdir -p Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-NO_Elag/PAR_$PAR/Solution_100MHz/cosim 
			cp SC_RTL_Simulation/solution_v7_100MHz/sim/report/vhdl/my_module.log Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-NO_Elag/PAR_$PAR/Solution_100MHz/cosim 

			# Copy the implementation report for Solution 100MHz
			mkdir -p Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-NO_Elag/PAR_$PAR/Solution_100MHz/impl 
			cp SC_RTL_Simulation/solution_v7_100MHz/impl/report/vhdl/my_module_export.rpt Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-NO_Elag/PAR_$PAR/Solution_100MHz/impl 


			# Copy the synthesis report files for Solution 300MHz
			mkdir -p Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-NO_Elag/PAR_$PAR/Solution_300MHz/syn 
			cp SC_RTL_Simulation/solution_v7_300MHz/syn/report/*.rpt Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-NO_Elag/PAR_$PAR/Solution_300MHz/syn 

			# Copy the vhdl files for Solution 300MHz
			cp -R SC_RTL_Simulation/solution_v7_300MHz/syn/vhdl Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-NO_Elag/PAR_$PAR/Solution_300MHz/syn 

			# Copy the cosim log for Solution 300MHz
			mkdir -p Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-NO_Elag/PAR_$PAR/Solution_300MHz/cosim 
			cp SC_RTL_Simulation/solution_v7_300MHz/sim/report/vhdl/my_module.log Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-NO_Elag/PAR_$PAR/Solution_300MHz/cosim 

			# Copy the implementation report for Solution 300MHz
			mkdir -p Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-NO_Elag/PAR_$PAR/Solution_300MHz/impl 
			cp SC_RTL_Simulation/solution_v7_300MHz/impl/report/vhdl/my_module_export.rpt Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-NO_Elag/PAR_$PAR/Solution_300MHz/impl 

			# Erase the Solution for new test
			cd SC_RTL_Simulation/ 
			rm -r solution_v7_100MHz/ 
			rm -r solution_v7_300MHz/ 
				
		done
	done
done

exit 0