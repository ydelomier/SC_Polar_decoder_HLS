#!/bin/bash
set -x

declare -a NBIT
declare -a KBIT
declare -a QUANT


NBIT=(65536 262144 32768 32768 32768 131072)
KBIT=(58982 235930 29492 29492 29492 117964)
QUANT=(6 6 9 7 8 5)

for ELAG in 0 2
do

	for REF in 1 5
	do

		for PAR in 8 16 64
		do

			# Generate Proper parameters file
			cd /c/Polar_decoder_HLS/Polar_Code_Generator/Frozen_Bit_Generator/cmake-build-debug/   
			./FB_Generator.exe ${NBIT[$REF]} ${KBIT[$REF]} $PAR 0 ../../Generated_Frozen_Bit/frozen_n_32768_k_16384.txt 1 /c/Polar_decoder_HLS/SC_Polar_decoder_v5_Full_Elag/src/module/
					
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
	#define ELAG_REP		1	// Repetition Node 			(ex: 00000001) , by default only for F_function
	#define ELAG_SPC		1	// Single Parity Check Node (ex: 01111111) , by default only for G_function
	#define ELAG_REP2		0	// REP2 Node 				(ex: 00000011) , by default only for F_function
	#define ELAG_SPC2		0	// SPC2 Node 				(ex: 00111111) , by default only for G_function

	// These options enable pruning for the rariest cases  ( F_R1, G_R0, F_SPC, ...)
	#define ELAG_RARE		0	// 1 enable, 0 off

	#define ELAG_H0			1	// When pruning a node R0 from a F function, the F_R0 state can be bypass at the cost of a dedicated H function : H0

	#define _MONITORING_ " > src/module/config.h

			# Copy the csim log
			mkdir -p Results_Implem_Latency/Ref_$REF-N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-Elag_$ELAG/PAR_$PAR/csim 
	
		done
	done
done

exit 0