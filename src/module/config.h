 
#define LLR_BITS		6 
#define LLR_MAXV    	(+127)
#define LLR_MINV    	(-127)

#define LLR 			sc_bigint<LLR_BITS>
#define BIT 			sc_biguint<1>
#define TYPE_LLRS   	sc_bigint<PAR * LLR_BITS>
#define TYPE_BITS   	sc_biguint<PAR>

#define SIGMAG 				// CA2 or SIGMAG , define the representation format
						// CA2 : 2's complement , SIGMAG : sign and magnitude

#define EXTENDED		1	// 1: on, 0: off. Enable the extension of the quantization format for the internal signals of the operators

#define PRUNING_LEVEL 	2	// 0: No Pruning, 1: Leaf Pruning, 2: Full Pruning enable (R0 node pruned by default when enable)

// Choose the type of node to pruned. 1 : enable, 0 : off
#define ELAG_R1			1	// Node of Rate 1 						   , by default only for G_function
#define ELAG_REP		1	// Repetition Node 			(ex: 00000001) , by default only for F_function
#define ELAG_SPC		1	// Single Parity Check Node (ex: 01111111) , by default only for G_function
#define ELAG_REP2		0	// REP2 Node 				(ex: 00000011) , by default only for F_function
#define ELAG_SPC2		0	// SPC2 Node 				(ex: 00111111) , by default only for G_function

// These options enable pruning for the rariest cases  ( F_R1, G_R0, F_SPC, ...)
#define ELAG_RARE		0	// 1 enable, 0 off

#define ELAG_H0			1	// When pruning a node R0 from a F function, the F_R0 state can be bypass at the cost of a dedicated H function : H0

#define _MONITORING_ 
