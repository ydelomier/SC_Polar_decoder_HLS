#!/bin/bash
set -x

rm parser.txt

declare -a NBIT
declare -a KBIT
declare -a QUANT
declare -a ELAG

NBIT=(1024)
KBIT=(512 922)
QUANT=(8 6)
ELAG=("NULL" "R0" "R0R1" "R0R1REP" "ROR1REPSPC")

for KREF in 0 1
do
echo "/*********************************/" >> parser.txt
echo "/******** K = ${KBIT[$KREF]} ********/" >> parser.txt
echo "/*********************************/" >> parser.txt
	
	for QREF in 0 1
	do
echo "/******** Q = ${QUANT[$QREF]} ********/" >> parser.txt
echo "/*********************************/" >> parser.txt

		for PARA in 2 4 8 16 32 64
		do
		echo "$$$$$ P = $PARA $$$$$" >> parser.txt
		
			for ELREF in 0 1 2 3 4		
			do
			echo "%%% ELAG_${ELAG[$ELREF] %%%" >> parser.txt
			echo "   # 100MHz" >> parser.txt
	# Impl Result
	cat /c/THESIS/Polar_decoder_HLS/SC_Polar_decoder/RESULTAT_TECS_POLAR/N_${NBIT[0]}-K_${KBIT[$KREF]}-Q_${QUANT[$QREF]}/PAR_$PARA/ELAG_${ELAG[$ELREF]}/Solution_100MHz/impl/my_module_export.rpt | grep SLICE: | awk '{ print $2 }' >> parser.txt
	cat /c/THESIS/Polar_decoder_HLS/SC_Polar_decoder/RESULTAT_TECS_POLAR/N_${NBIT[0]}-K_${KBIT[$KREF]}-Q_${QUANT[$QREF]}/PAR_$PARA/ELAG_${ELAG[$ELREF]}/Solution_100MHz/impl/my_module_export.rpt | grep LUT: | awk '{ print $2 }' >> parser.txt
	cat /c/THESIS/Polar_decoder_HLS/SC_Polar_decoder/RESULTAT_TECS_POLAR/N_${NBIT[0]}-K_${KBIT[$KREF]}-Q_${QUANT[$QREF]}/PAR_$PARA/ELAG_${ELAG[$ELREF]}/Solution_100MHz/impl/my_module_export.rpt | grep FF: | awk '{ print $2 }' >> parser.txt
	cat /c/THESIS/Polar_decoder_HLS/SC_Polar_decoder/RESULTAT_TECS_POLAR/N_${NBIT[0]}-K_${KBIT[$KREF]}-Q_${QUANT[$QREF]}/PAR_$PARA/ELAG_${ELAG[$ELREF]}/Solution_100MHz/impl/my_module_export.rpt | grep DSP: | awk '{ print $2 }' >> parser.txt
	cat /c/THESIS/Polar_decoder_HLS/SC_Polar_decoder/RESULTAT_TECS_POLAR/N_${NBIT[0]}-K_${KBIT[$KREF]}-Q_${QUANT[$QREF]}/PAR_$PARA/ELAG_${ELAG[$ELREF]}/Solution_100MHz/impl/my_module_export.rpt | grep BRAM: | awk '{ print $2 }' >> parser.txt
	cat /c/THESIS/Polar_decoder_HLS/SC_Polar_decoder/RESULTAT_TECS_POLAR/N_${NBIT[0]}-K_${KBIT[$KREF]}-Q_${QUANT[$QREF]}/PAR_$PARA/ELAG_${ELAG[$ELREF]}/Solution_100MHz/impl/my_module_export.rpt | grep post-implementation: | awk '{ print $4 }' >> parser.txt
	# Latence Mesurer
	cat /c/THESIS/Polar_decoder_HLS/SC_Polar_decoder/RESULTAT_TECS_POLAR/N_${NBIT[0]}-K_${KBIT[$KREF]}-Q_${QUANT[$QREF]}/PAR_$PARA/ELAG_${ELAG[$ELREF]}/Solution_100MHz/cosim/my_module.log | grep PROCESS | awk '{ print $5 }' >> parser.txt
	cat /c/THESIS/Polar_decoder_HLS/SC_Polar_decoder/RESULTAT_TECS_POLAR/N_${NBIT[0]}-K_${KBIT[$KREF]}-Q_${QUANT[$QREF]}/PAR_$PARA/ELAG_${ELAG[$ELREF]}/Solution_100MHz/cosim/my_module.log | grep STORE | awk '{ print $5 }' >> parser.txt
			echo "   # 200MHz" >> parser.txt
	# Impl Result
	cat /c/THESIS/Polar_decoder_HLS/SC_Polar_decoder/RESULTAT_TECS_POLAR/N_${NBIT[0]}-K_${KBIT[$KREF]}-Q_${QUANT[$QREF]}/PAR_$PARA/ELAG_${ELAG[$ELREF]}/Solution_200MHz/impl/my_module_export.rpt | grep SLICE: | awk '{ print $2 }' >> parser.txt
	cat /c/THESIS/Polar_decoder_HLS/SC_Polar_decoder/RESULTAT_TECS_POLAR/N_${NBIT[0]}-K_${KBIT[$KREF]}-Q_${QUANT[$QREF]}/PAR_$PARA/ELAG_${ELAG[$ELREF]}/Solution_200MHz/impl/my_module_export.rpt | grep LUT: | awk '{ print $2 }' >> parser.txt
	cat /c/THESIS/Polar_decoder_HLS/SC_Polar_decoder/RESULTAT_TECS_POLAR/N_${NBIT[0]}-K_${KBIT[$KREF]}-Q_${QUANT[$QREF]}/PAR_$PARA/ELAG_${ELAG[$ELREF]}/Solution_200MHz/impl/my_module_export.rpt | grep FF: | awk '{ print $2 }' >> parser.txt
	cat /c/THESIS/Polar_decoder_HLS/SC_Polar_decoder/RESULTAT_TECS_POLAR/N_${NBIT[0]}-K_${KBIT[$KREF]}-Q_${QUANT[$QREF]}/PAR_$PARA/ELAG_${ELAG[$ELREF]}/Solution_200MHz/impl/my_module_export.rpt | grep DSP: | awk '{ print $2 }' >> parser.txt
	cat /c/THESIS/Polar_decoder_HLS/SC_Polar_decoder/RESULTAT_TECS_POLAR/N_${NBIT[0]}-K_${KBIT[$KREF]}-Q_${QUANT[$QREF]}/PAR_$PARA/ELAG_${ELAG[$ELREF]}/Solution_200MHz/impl/my_module_export.rpt | grep BRAM: | awk '{ print $2 }' >> parser.txt
	cat /c/THESIS/Polar_decoder_HLS/SC_Polar_decoder/RESULTAT_TECS_POLAR/N_${NBIT[0]}-K_${KBIT[$KREF]}-Q_${QUANT[$QREF]}/PAR_$PARA/ELAG_${ELAG[$ELREF]}/Solution_200MHz/impl/my_module_export.rpt | grep post-implementation: | awk '{ print $4 }' >> parser.txt
	# Latence Mesurer
	cat /c/THESIS/Polar_decoder_HLS/SC_Polar_decoder/RESULTAT_TECS_POLAR/N_${NBIT[0]}-K_${KBIT[$KREF]}-Q_${QUANT[$QREF]}/PAR_$PARA/ELAG_${ELAG[$ELREF]}/Solution_200MHz/cosim/my_module.log | grep PROCESS | awk '{ print $5 }' >> parser.txt
	cat /c/THESIS/Polar_decoder_HLS/SC_Polar_decoder/RESULTAT_TECS_POLAR/N_${NBIT[0]}-K_${KBIT[$KREF]}-Q_${QUANT[$QREF]}/PAR_$PARA/ELAG_${ELAG[$ELREF]}/Solution_200MHz/cosim/my_module.log | grep STORE | awk '{ print $5 }' >> parser.txt
			echo "   # 300MHz" >> parser.txt
	# Impl Result
	cat /c/THESIS/Polar_decoder_HLS/SC_Polar_decoder/RESULTAT_TECS_POLAR/N_${NBIT[0]}-K_${KBIT[$KREF]}-Q_${QUANT[$QREF]}/PAR_$PARA/ELAG_${ELAG[$ELREF]}/Solution_300MHz/impl/my_module_export.rpt | grep SLICE: | awk '{ print $2 }' >> parser.txt
	cat /c/THESIS/Polar_decoder_HLS/SC_Polar_decoder/RESULTAT_TECS_POLAR/N_${NBIT[0]}-K_${KBIT[$KREF]}-Q_${QUANT[$QREF]}/PAR_$PARA/ELAG_${ELAG[$ELREF]}/Solution_300MHz/impl/my_module_export.rpt | grep LUT: | awk '{ print $2 }' >> parser.txt
	cat /c/THESIS/Polar_decoder_HLS/SC_Polar_decoder/RESULTAT_TECS_POLAR/N_${NBIT[0]}-K_${KBIT[$KREF]}-Q_${QUANT[$QREF]}/PAR_$PARA/ELAG_${ELAG[$ELREF]}/Solution_300MHz/impl/my_module_export.rpt | grep FF: | awk '{ print $2 }' >> parser.txt
	cat /c/THESIS/Polar_decoder_HLS/SC_Polar_decoder/RESULTAT_TECS_POLAR/N_${NBIT[0]}-K_${KBIT[$KREF]}-Q_${QUANT[$QREF]}/PAR_$PARA/ELAG_${ELAG[$ELREF]}/Solution_300MHz/impl/my_module_export.rpt | grep DSP: | awk '{ print $2 }' >> parser.txt
	cat /c/THESIS/Polar_decoder_HLS/SC_Polar_decoder/RESULTAT_TECS_POLAR/N_${NBIT[0]}-K_${KBIT[$KREF]}-Q_${QUANT[$QREF]}/PAR_$PARA/ELAG_${ELAG[$ELREF]}/Solution_300MHz/impl/my_module_export.rpt | grep BRAM: | awk '{ print $2 }' >> parser.txt
	cat /c/THESIS/Polar_decoder_HLS/SC_Polar_decoder/RESULTAT_TECS_POLAR/N_${NBIT[0]}-K_${KBIT[$KREF]}-Q_${QUANT[$QREF]}/PAR_$PARA/ELAG_${ELAG[$ELREF]}/Solution_300MHz/impl/my_module_export.rpt | grep post-implementation: | awk '{ print $4 }' >> parser.txt
	# Latence Mesurer
	cat /c/THESIS/Polar_decoder_HLS/SC_Polar_decoder/RESULTAT_TECS_POLAR/N_${NBIT[0]}-K_${KBIT[$KREF]}-Q_${QUANT[$QREF]}/PAR_$PARA/ELAG_${ELAG[$ELREF]}/Solution_300MHz/cosim/my_module.log | grep PROCESS | awk '{ print $5 }' >> parser.txt
	cat /c/THESIS/Polar_decoder_HLS/SC_Polar_decoder/RESULTAT_TECS_POLAR/N_${NBIT[0]}-K_${KBIT[$KREF]}-Q_${QUANT[$QREF]}/PAR_$PARA/ELAG_${ELAG[$ELREF]}/Solution_300MHz/cosim/my_module.log | grep STORE | awk '{ print $5 }' >> parser.txt

			done
		done	
	done
done

exit 0