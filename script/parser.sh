#!/bin/bash
set -x

rm parser.txt

declare -a NBIT
declare -a KBIT
declare -a QUANT

NBIT=(32768 32768)
KBIT=(29492 29492)
QUANT=(8 6)

echo "/*********************************/" >> parser.txt
echo "/************ CA2 ************/" >> parser.txt
echo "/*********************************/" >> parser.txt

for REF in 0 1
do
	echo "/*****  QUANT ${QUANT[$REF]}  *****/" >> parser.txt
	for PAR in 4 8 16 32 64 
	do
		echo "%%% [PAR $PAR] %%%" >> parser.txt

		# # Latence Théorique
		# cat /c/Polar_decoder_HLS/SC_Polar_decoder_v5_Full_Elag/Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-NO_Elag/PAR_$PAR/csim/my_module_csim.log | grep PROCESS | awk '{ print $5 }' >> parser.txt
		
		# Impl Result
		cat /c/Polar_decoder_HLS/SC_Polar_decoder_v5_Full_Elag/Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-FULL_Elag/PAR_$PAR/Solution_100MHz/impl/my_module_export.rpt | grep SLICE: | awk '{ print $2 }' >> parser.txt
		cat /c/Polar_decoder_HLS/SC_Polar_decoder_v5_Full_Elag/Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-FULL_Elag/PAR_$PAR/Solution_100MHz/impl/my_module_export.rpt | grep LUT: | awk '{ print $2 }' >> parser.txt
		cat /c/Polar_decoder_HLS/SC_Polar_decoder_v5_Full_Elag/Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-FULL_Elag/PAR_$PAR/Solution_100MHz/impl/my_module_export.rpt | grep FF: | awk '{ print $2 }' >> parser.txt
		cat /c/Polar_decoder_HLS/SC_Polar_decoder_v5_Full_Elag/Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-FULL_Elag/PAR_$PAR/Solution_100MHz/impl/my_module_export.rpt | grep DSP: | awk '{ print $2 }' >> parser.txt
		cat /c/Polar_decoder_HLS/SC_Polar_decoder_v5_Full_Elag/Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-FULL_Elag/PAR_$PAR/Solution_100MHz/impl/my_module_export.rpt | grep BRAM: | awk '{ print $2 }' >> parser.txt
		cat /c/Polar_decoder_HLS/SC_Polar_decoder_v5_Full_Elag/Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-FULL_Elag/PAR_$PAR/Solution_100MHz/impl/my_module_export.rpt | grep post-implementation: | awk '{ print $4 }' >> parser.txt
		# Latence Mesurer
		cat /c/Polar_decoder_HLS/SC_Polar_decoder_v5_Full_Elag/Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-FULL_Elag/PAR_$PAR/Solution_100MHz/cosim/my_module.log | grep PROCESS | awk '{ print $5 }' >> parser.txt
		cat /c/Polar_decoder_HLS/SC_Polar_decoder_v5_Full_Elag/Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-FULL_Elag/PAR_$PAR/Solution_100MHz/cosim/my_module.log | grep STORE | awk '{ print $5 }' >> parser.txt
	
	done
done

echo "/*********************************/" >> parser.txt
echo "/************ SIGMAG ************/" >> parser.txt
echo "/*********************************/" >> parser.txt

for REF in 0 1
do
	echo "/*****  QUANT ${QUANT[$REF]}  *****/" >> parser.txt
	for PAR in 4 8 16 32 64 
	do
		echo "%%% [PAR $PAR] %%%" >> parser.txt

		# # Latence Théorique
		# cat /c/Polar_decoder_HLS/SC_Polar_decoder_v5_Full_Elag/Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-NO_Elag/PAR_$PAR/csim/my_module_csim.log | grep PROCESS | awk '{ print $5 }' >> parser.txt
		
		# Impl Result
		cat /c/Polar_decoder_HLS/SC_Polar_decoder_v5_Full_Elag/Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-SIGMAG-Full_Elag/PAR_$PAR/Solution_100MHz/impl/my_module_export.rpt | grep SLICE: | awk '{ print $2 }' >> parser.txt
		cat /c/Polar_decoder_HLS/SC_Polar_decoder_v5_Full_Elag/Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-SIGMAG-Full_Elag/PAR_$PAR/Solution_100MHz/impl/my_module_export.rpt | grep LUT: | awk '{ print $2 }' >> parser.txt
		cat /c/Polar_decoder_HLS/SC_Polar_decoder_v5_Full_Elag/Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-SIGMAG-Full_Elag/PAR_$PAR/Solution_100MHz/impl/my_module_export.rpt | grep FF: | awk '{ print $2 }' >> parser.txt
		cat /c/Polar_decoder_HLS/SC_Polar_decoder_v5_Full_Elag/Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-SIGMAG-Full_Elag/PAR_$PAR/Solution_100MHz/impl/my_module_export.rpt | grep DSP: | awk '{ print $2 }' >> parser.txt
		cat /c/Polar_decoder_HLS/SC_Polar_decoder_v5_Full_Elag/Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-SIGMAG-Full_Elag/PAR_$PAR/Solution_100MHz/impl/my_module_export.rpt | grep BRAM: | awk '{ print $2 }' >> parser.txt
		cat /c/Polar_decoder_HLS/SC_Polar_decoder_v5_Full_Elag/Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-SIGMAG-Full_Elag/PAR_$PAR/Solution_100MHz/impl/my_module_export.rpt | grep post-implementation: | awk '{ print $4 }' >> parser.txt
		# Latence Mesurer
		cat /c/Polar_decoder_HLS/SC_Polar_decoder_v5_Full_Elag/Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-SIGMAG-Full_Elag/PAR_$PAR/Solution_100MHz/cosim/my_module.log | grep PROCESS | awk '{ print $5 }' >> parser.txt
		cat /c/Polar_decoder_HLS/SC_Polar_decoder_v5_Full_Elag/Results_Implem_Latency/N_${NBIT[$REF]}-K_${KBIT[$REF]}-Q_${QUANT[$REF]}-SIGMAG-Full_Elag/PAR_$PAR/Solution_100MHz/cosim/my_module.log | grep STORE | awk '{ print $5 }' >> parser.txt
	
	done
done

# echo "/*********************************/" >> parser.txt
# echo "/************ NO ELAG ************/" >> parser.txt
# echo "/*********************************/" >> parser.txt

# echo "/*****  Solution 100MHZ  *****/" >> parser.txt
# for PAR in 4 8 1632 64 
# do
	# echo "%%% [PAR $PAR] %%%" >> parser.txt

	# # Impl Result
	# cat /c/Polar_decoder_HLS/SC_Polar_decoder_v5_Full_Elag/Results_Implem_Latency/N_$NBIT-K_$KBIT-Q_$QUANT-NO_Elag/PAR_$PAR/Solution_100MHz/impl/my_module_export.rpt | grep SLICE: | awk '{ print $2 }' >> parser.txt
	# cat /c/Polar_decoder_HLS/SC_Polar_decoder_v5_Full_Elag/Results_Implem_Latency/N_$NBIT-K_$KBIT-Q_$QUANT-NO_Elag/PAR_$PAR/Solution_100MHz/impl/my_module_export.rpt | grep LUT: | awk '{ print $2 }' >> parser.txt
	# cat /c/Polar_decoder_HLS/SC_Polar_decoder_v5_Full_Elag/Results_Implem_Latency/N_$NBIT-K_$KBIT-Q_$QUANT-NO_Elag/PAR_$PAR/Solution_100MHz/impl/my_module_export.rpt | grep FF: | awk '{ print $2 }' >> parser.txt
	# cat /c/Polar_decoder_HLS/SC_Polar_decoder_v5_Full_Elag/Results_Implem_Latency/N_$NBIT-K_$KBIT-Q_$QUANT-NO_Elag/PAR_$PAR/Solution_100MHz/impl/my_module_export.rpt | grep DSP: | awk '{ print $2 }' >> parser.txt
	# cat /c/Polar_decoder_HLS/SC_Polar_decoder_v5_Full_Elag/Results_Implem_Latency/N_$NBIT-K_$KBIT-Q_$QUANT-NO_Elag/PAR_$PAR/Solution_100MHz/impl/my_module_export.rpt | grep BRAM: | awk '{ print $2 }' >> parser.txt
	# cat /c/Polar_decoder_HLS/SC_Polar_decoder_v5_Full_Elag/Results_Implem_Latency/N_$NBIT-K_$KBIT-Q_$QUANT-NO_Elag/PAR_$PAR/Solution_100MHz/impl/my_module_export.rpt | grep post-implementation: | awk '{ print $4 }' >> parser.txt
	# # Latence Mesurer
	# cat /c/Polar_decoder_HLS/SC_Polar_decoder_v5_Full_Elag/Results_Implem_Latency/N_$NBIT-K_$KBIT-Q_$QUANT-NO_Elag/PAR_$PAR/Solution_100MHz/cosim/my_module.log | grep PROCESS | awk '{ print $5 }' >> parser.txt
	# cat /c/Polar_decoder_HLS/SC_Polar_decoder_v5_Full_Elag/Results_Implem_Latency/N_$NBIT-K_$KBIT-Q_$QUANT-NO_Elag/PAR_$PAR/Solution_100MHz/cosim/my_module.log | grep STORE | awk '{ print $5 }' >> parser.txt
# done

# echo "/*****  Solution 200MHZ  *****/" >> parser.txt
# for PAR in 4 8 16 32 64 
# do
	# echo "%%% [PAR $PAR] %%%" >> parser.txt

	# # Impl Result
	# cat /c/Polar_decoder_HLS/SC_Polar_decoder_v5_Full_Elag/Results_Implem_Latency/N_$NBIT-K_$KBIT-Q_$QUANT-NO_Elag/PAR_$PAR/Solution_200MHz/impl/my_module_export.rpt | grep SLICE: | awk '{ print $2 }' >> parser.txt
	# cat /c/Polar_decoder_HLS/SC_Polar_decoder_v5_Full_Elag/Results_Implem_Latency/N_$NBIT-K_$KBIT-Q_$QUANT-NO_Elag/PAR_$PAR/Solution_200MHz/impl/my_module_export.rpt | grep LUT: | awk '{ print $2 }' >> parser.txt
	# cat /c/Polar_decoder_HLS/SC_Polar_decoder_v5_Full_Elag/Results_Implem_Latency/N_$NBIT-K_$KBIT-Q_$QUANT-NO_Elag/PAR_$PAR/Solution_200MHz/impl/my_module_export.rpt | grep FF: | awk '{ print $2 }' >> parser.txt
	# cat /c/Polar_decoder_HLS/SC_Polar_decoder_v5_Full_Elag/Results_Implem_Latency/N_$NBIT-K_$KBIT-Q_$QUANT-NO_Elag/PAR_$PAR/Solution_200MHz/impl/my_module_export.rpt | grep DSP: | awk '{ print $2 }' >> parser.txt
	# cat /c/Polar_decoder_HLS/SC_Polar_decoder_v5_Full_Elag/Results_Implem_Latency/N_$NBIT-K_$KBIT-Q_$QUANT-NO_Elag/PAR_$PAR/Solution_200MHz/impl/my_module_export.rpt | grep BRAM: | awk '{ print $2 }' >> parser.txt
	# cat /c/Polar_decoder_HLS/SC_Polar_decoder_v5_Full_Elag/Results_Implem_Latency/N_$NBIT-K_$KBIT-Q_$QUANT-NO_Elag/PAR_$PAR/Solution_200MHz/impl/my_module_export.rpt | grep post-implementation: | awk '{ print $4 }' >> parser.txt
	# # Latence Mesurer
	# cat /c/Polar_decoder_HLS/SC_Polar_decoder_v5_Full_Elag/Results_Implem_Latency/N_$NBIT-K_$KBIT-Q_$QUANT-NO_Elag/PAR_$PAR/Solution_200MHz/cosim/my_module.log | grep PROCESS | awk '{ print $5 }' >> parser.txt
	# cat /c/Polar_decoder_HLS/SC_Polar_decoder_v5_Full_Elag/Results_Implem_Latency/N_$NBIT-K_$KBIT-Q_$QUANT-NO_Elag/PAR_$PAR/Solution_200MHz/cosim/my_module.log | grep STORE | awk '{ print $5 }' >> parser.txt
# done

# echo "/***** Solution 300MHZ  *****/" >> parser.txt
# for PAR in 4 8 16 32 64 
# do
	# echo "%%% [PAR $PAR] %%%" >> parser.txt

	# # Impl Result
	# cat /c/Polar_decoder_HLS/SC_Polar_decoder_v5_Full_Elag/Results_Implem_Latency/N_$NBIT-K_$KBIT-Q_$QUANT-NO_Elag/PAR_$PAR/Solution_300MHz/impl/my_module_export.rpt | grep SLICE: | awk '{ print $2 }' >> parser.txt
	# cat /c/Polar_decoder_HLS/SC_Polar_decoder_v5_Full_Elag/Results_Implem_Latency/N_$NBIT-K_$KBIT-Q_$QUANT-NO_Elag/PAR_$PAR/Solution_300MHz/impl/my_module_export.rpt | grep LUT: | awk '{ print $2 }' >> parser.txt
	# cat /c/Polar_decoder_HLS/SC_Polar_decoder_v5_Full_Elag/Results_Implem_Latency/N_$NBIT-K_$KBIT-Q_$QUANT-NO_Elag/PAR_$PAR/Solution_300MHz/impl/my_module_export.rpt | grep FF: | awk '{ print $2 }' >> parser.txt
	# cat /c/Polar_decoder_HLS/SC_Polar_decoder_v5_Full_Elag/Results_Implem_Latency/N_$NBIT-K_$KBIT-Q_$QUANT-NO_Elag/PAR_$PAR/Solution_300MHz/impl/my_module_export.rpt | grep DSP: | awk '{ print $2 }' >> parser.txt
	# cat /c/Polar_decoder_HLS/SC_Polar_decoder_v5_Full_Elag/Results_Implem_Latency/N_$NBIT-K_$KBIT-Q_$QUANT-NO_Elag/PAR_$PAR/Solution_300MHz/impl/my_module_export.rpt | grep BRAM: | awk '{ print $2 }' >> parser.txt
	# cat /c/Polar_decoder_HLS/SC_Polar_decoder_v5_Full_Elag/Results_Implem_Latency/N_$NBIT-K_$KBIT-Q_$QUANT-NO_Elag/PAR_$PAR/Solution_300MHz/impl/my_module_export.rpt | grep post-implementation: | awk '{ print $4 }' >> parser.txt
	# # Latence Mesuré
	# cat /c/Polar_decoder_HLS/SC_Polar_decoder_v5_Full_Elag/Results_Implem_Latency/N_$NBIT-K_$KBIT-Q_$QUANT-NO_Elag/PAR_$PAR/Solution_300MHz/cosim/my_module.log | grep PROCESS | awk '{ print $5 }' >> parser.txt
	# cat /c/Polar_decoder_HLS/SC_Polar_decoder_v5_Full_Elag/Results_Implem_Latency/N_$NBIT-K_$KBIT-Q_$QUANT-NO_Elag/PAR_$PAR/Solution_300MHz/cosim/my_module.log | grep STORE | awk '{ print $5 }' >> parser.txt

# done


exit 0