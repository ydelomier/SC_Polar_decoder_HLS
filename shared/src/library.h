#ifndef LIBRARY
#define LIBRARY

#include "systemc.h"
#include "functions.h"
#include "sc_list_fct.h"

//#define DEBUG_LIB

//*************************************************************************//
//** 								 									 **//
//** 						FUNCTIONS LIBRARY							 **//
//** 								 									 **//
//*************************************************************************//



template <int Q>
sc_bigint<Q> Adapt_format(sc_bigint<Q> a){  // Convert CA2 to SIGMAG / SIGMAG to CA2
#pragma HLS INLINE

#if defined CA2
	sc_bigint<Q> res = a;
#elif defined SIGMAG
	sc_bigint <Q> res = (sc_bigint<Q>) qconv_format( (sc_biguint<Q>) a );
#endif
	return res;
}

//*************************************************************************//
//** 						PRUNNED NODE DEFINITION						 **//
//*************************************************************************//

#define NODE_R0		0x00  	// Node R0
#define NODE_R1 	0x0F  	// Node R1
#define NODE_REP 	0x02  	// Node REP
#define NODE_SPC 	0x04 	// Node SPC
#define NODE_REP2	0x03  	// Node REP2
#define NODE_SPC2 	0x05  	// Node SPC2
#define NODE_RN		0x08  	// Node ?


//*************************************************************************//
//** 						POLAR CODE FUNCTIONS						 **//
//*************************************************************************//

template <int P, int Q>  // P : Parallelism , Q : Quantification
sc_bigint<P*Q> PU_FUNCTION_F(sc_bigint<P*Q> la, sc_bigint<P*Q> lb){
#pragma HLS INLINE off
	return Function_F<P,Q> (la,lb);
}

template <int P, int Q>
sc_bigint<P*Q> PU_FUNCTION_G(sc_bigint<P*Q> la, sc_bigint<P*Q> lb, sc_biguint<P> sa){
#pragma HLS INLINE off
	return Function_G<P,Q> (la,lb,sa);
}

template <int P>
sc_biguint<P> PU_FUNCTION_H(sc_biguint<P> sa, sc_biguint<P> sb){
#pragma HLS INLINE off
	return Function_H<P> (sa,sb);
}


//*************************************************************************//
//** 					PRUNING SPECIALIZED FUNCTIONS					 **//
//*************************************************************************//


/////////////// Adder tree for REP Node ////////////////
template <int P, int Q>
sc_bigint< Q + LOG2_PAR + 1 > ADD_TREE_FUNCTION (sc_bigint< P * Q > llr, sc_bigint< Q + LOG2_PAR + 1 > old_sum)
{
#pragma HLS INLINE
#if PAR == 2
	return ADDER_TREE_2< Q > (llr, old_sum);
#elif PAR == 4
	return ADDER_TREE_4< Q > (llr, old_sum);
#elif PAR == 8
	return ADDER_TREE_8< Q > (llr, old_sum);
#elif PAR == 16
	return ADDER_TREE_16< Q > (llr, old_sum);
#elif PAR == 32
	return ADDER_TREE_32< Q > (llr, old_sum);
#elif PAR == 64
	return ADDER_TREE_64< Q > (llr, old_sum);
#elif PAR == 128
	return ADDER_TREE_128< Q > (llr, old_sum);
#elif PAR == 256
	return ADDER_TREE_256< Q > (llr, old_sum);
#endif
}

/////////////// XOR tree for SPC Node ////////////////

template <int P>
sc_biguint< 1 > PARITY_TREE_FUNCTION (sc_biguint< P > sign, sc_biguint< 1 > old_parity)
{
#pragma HLS INLINE
#if PAR == 2
	return Parity_TREE_2 (sign, old_parity);
#elif PAR == 4
	return Parity_TREE_4 (sign, old_parity);
#elif PAR == 8
	return Parity_TREE_8 (sign, old_parity);
#elif PAR == 16
	return Parity_TREE_16 (sign, old_parity);
#elif PAR == 32
	return Parity_TREE_32 (sign, old_parity);
#elif PAR == 64
	return Parity_TREE_64 (sign, old_parity);
#elif PAR == 128
	return Parity_TREE_128 (sign, old_parity);
#elif PAR == 256
	return Parity_TREE_256 (sign, old_parity);
#endif
}

/////////////// MIN tree and MASK for SPC Node ////////////////

template <int P, int Q>
sc_bigint< Q + P > MIN_MASK_TREE_FCT (sc_bigint< P * Q> llr)
{
#pragma HLS INLINE
#if PAR == 2
	return Min_Mask_TREE_2< Q > (llr);
#elif PAR == 4
	return Min_Mask_TREE_4< Q > (llr);
#elif PAR == 8
	return Min_Mask_TREE_8< Q > (llr);
#elif PAR == 16
	return Min_Mask_TREE_16< Q > (llr);
#elif PAR == 32
	return Min_Mask_TREE_32< Q > (llr);
#elif PAR == 64
	return Min_Mask_TREE_64< Q > (llr);
#elif PAR == 128
	return Min_Mask_TREE_128< Q > (llr);
#elif PAR == 256
	return Min_Mask_TREE_256< Q > (llr);
#endif
}

//*************************************************************************//
//** 					LEAF SPECIALIZED FUNCTIONS						 **//
//*************************************************************************//

template <int P, int Q>
sc_biguint< P > Spec_Polar_Decoder ( sc_bigint< P *Q> llr, sc_biguint< P > fb )
{
#pragma HLS INLINE
#if PAR == 1
	return Spec_PolarDec_1< Q > (llr, fb);
#elif PAR == 2
	return Spec_PolarDec_2< Q > (llr, fb);
#elif PAR == 4
	return Spec_PolarDec_4< Q > (llr, fb);
#elif PAR == 8
	return Spec_PolarDec_8< Q > (llr, fb);
#elif PAR == 16
	return Spec_PolarDec_16< Q > (llr, fb);
#elif PAR == 32
	return Spec_PolarDec_32< Q > (llr, fb);
#elif PAR == 64
	return Spec_PolarDec_64< Q > (llr, fb);
#elif PAR == 128
	return Spec_PolarDec_128< Q > (llr, fb);
#elif PAR == 256
	return Spec_PolarDec_256< Q > (llr, fb);
#endif
}

template <int P, int Q>
sc_biguint< P > Spec_Node_R0 (sc_bigint< P * Q> llr)
{
#pragma HLS INLINE
	return (sc_biguint< P >) 0;
}

template <int P, int Q>
sc_biguint< P > Spec_Node_R1 (sc_bigint< P * Q> llr)
{
#pragma HLS INLINE
	return VECTOR_SIGN<P, Q>( llr);
}

template <int P, int Q>
sc_biguint< P > Spec_REP_Node (sc_bigint< P * Q> llr)
{
#pragma HLS INLINE
#if PAR == 2
	return REP_2_Node< Q > (llr);
#elif PAR == 4
	return REP_4_Node< Q > (llr);
#elif PAR == 8
	return REP_8_Node< Q > (llr);
#elif PAR == 16
	return REP_16_Node< Q > (llr);
#elif PAR == 32
	return REP_32_Node< Q > (llr);
#elif PAR == 64
	return REP_64_Node< Q > (llr);
#elif PAR == 128
	return REP_128_Node< Q > (llr);
#elif PAR == 256
	return REP_256_Node< Q > (llr);
#endif
}

template <int P, int Q>
sc_biguint< P > Spec_REP_REP2_Node (sc_bigint< P * Q> llr, bool sel) // Sel 0 : REP, Sel 1 : REP2
{
#pragma HLS INLINE
#if PAR == 2
	return REP_REP2_2_Node< Q > (llr, sel);
#elif PAR == 4
	return REP_REP2_4_Node< Q > (llr, sel);
#elif PAR == 8
	return REP_REP2_8_Node< Q > (llr, sel);
#elif PAR == 16
	return REP_REP2_16_Node< Q > (llr, sel);
#elif PAR == 32
	return REP_REP2_32_Node< Q > (llr, sel);
#elif PAR == 64
	return REP_REP2_64_Node< Q > (llr, sel);
#elif PAR == 128
	return REP_REP2_128_Node< Q > (llr, sel);
#elif PAR == 256
	return REP_REP2_256_Node< Q > (llr, sel);
#endif
}

template <int P, int Q>
sc_biguint< P > Spec_SPC_Node (sc_bigint< P * Q> llr)
{
#pragma HLS INLINE
#if PAR == 2
	return SPC_Node_2< Q > (llr);
#elif PAR == 4
	return SPC_Node_4< Q > (llr);
#elif PAR == 8
	return SPC_Node_8< Q > (llr);
#elif PAR == 16
	return SPC_Node_16< Q > (llr);
#elif PAR == 32
	return SPC_Node_32< Q > (llr);
#elif PAR == 64
	return SPC_Node_64< Q > (llr);
#elif PAR == 128
	return SPC_Node_128< Q > (llr);
#elif PAR == 256
	return SPC_Node_256< Q > (llr);
#endif
}

template <int P, int Q>
sc_biguint< P > Spec_SPC_SPC2_Node (sc_bigint< P * Q> llr, bool sel)
{
#pragma HLS INLINE
#if PAR == 2
	return SPC_SPC2_Node_2< Q > (llr, sel);
#elif PAR == 4
	return SPC_SPC2_Node_4< Q > (llr, sel);
#elif PAR == 8
	return SPC_SPC2_Node_8< Q > (llr, sel);
#elif PAR == 16
	return SPC_SPC2_Node_16< Q > (llr, sel);
#elif PAR == 32
	return SPC_SPC2_Node_32< Q > (llr, sel);
#elif PAR == 64
	return SPC_SPC2_Node_64< Q > (llr, sel);
#elif PAR == 128
	return SPC_SPC2_Node_128< Q > (llr, sel);
#elif PAR == 256
	return SPC_SPC2_Node_256< Q > (llr, sel);
#endif
}


//*************************************************************************//
//** 						PRUNNED NODE DEFINITION						 **//
//*************************************************************************//

template <int P, int L, int Q, int LOG2L, int D>
void Spec_LIST_Decoder ( LLR_struct<P,Q,LOG2L> input[L], PS_struct<P,Q,LOG2L> output[L], sc_biguint<P> fb, sc_biguint<LOG2L*D> LIST_STACK[L] )
{
#pragma HLS INLINE off
#if PAR == 1
	return List_R_P1 <L, Q, LOG2L, D> (input, output, fb, LIST_STACK);
#elif PAR == 2
	return List_R_P2 <L, Q, LOG2L, D> (input, output, fb, LIST_STACK);
#elif PAR == 4
	return List_R_P4 <L, Q, LOG2L, D> (input, output, fb, LIST_STACK);
#elif PAR == 8
	return List_R_P8 <L, Q, LOG2L, D> (input, output, fb, LIST_STACK);
#elif PAR == 16
	return List_R_P16 <L, Q, LOG2L, D> (input, output, fb, LIST_STACK);
#elif PAR == 32
	return List_R_P32 <L, Q, LOG2L, D> (input, output, fb, LIST_STACK);
#elif PAR == 64
	return List_R_P64 <L, Q, LOG2L, D> (input, output, fb, LIST_STACK);
#elif PAR == 128
	return List_R_P128 <L, Q, LOG2L, D> (input, output, fb, LIST_STACK);
#elif PAR == 256
	return List_R_P256 <L, Q, LOG2L, D> (input, output, fb, LIST_STACK);
#endif
}

#endif
