#ifndef SC_LIST_FCT
#define SC_LIST_FCT

#include "systemc.h"
#include "vector.h"
#include "functions.h"

//*************************************************************************//
//** 					LIST POLAR DECODER FUNCTIONS					 **//
//*************************************************************************//

template <int P, int Q, int LOG2L>
struct LLR_struct {
	sc_bigint < P *Q > 	llr;    // LLR
	sc_biguint< Q > 	metric; // metric
	sc_biguint< LOG2L > path;	// stack path
};

template <int P, int Q, int LOG2L>
struct PS_struct {
	sc_biguint< P > 	bit;    // PS
	sc_biguint< Q > 	metric; // metric
	sc_biguint< LOG2L > path;	// stack path
};


//*************************************************************************//
//** 						LIST MULTIPLEXERS							 **//
//*************************************************************************//

template <int L, int Q, int LOG2L>
sc_bigint<Q> MUX_LIST ( sc_bigint<Q> tab[L], sc_uint<LOG2L> sel)
{
#if L_SIZE == 2
	return LIST_MUX2<L,Q,LOG2L> ( tab, sel );
#elif L_SIZE == 4
	return LIST_MUX4<L,Q,LOG2L> ( tab, sel );
#elif L_SIZE == 8
	return LIST_MUX8<L,Q,LOG2L> ( tab, sel );
#elif L_SIZE == 16
	return LIST_MUX16<L,Q,LOG2L> ( tab, sel );
#elif L_SIZE == 32
	return LIST_MUX32<L,Q,LOG2L> ( tab, sel );
#else
	return LIST_MUX64<L,Q,LOG2L> ( tab, sel );
#endif

}

template <class T, int L, int Q, int LOG2L>
T MUX_LIST_t ( T tab[L], sc_uint<LOG2L> sel)
{
#if L_SIZE == 2
	return LIST_MUX2<T,L,Q,LOG2L> ( tab, sel );
#elif L_SIZE == 4
	return LIST_MUX4<T,L,Q,LOG2L> ( tab, sel );
#elif L_SIZE == 8
	return LIST_MUX8<T,L,Q,LOG2L> ( tab, sel );
#elif L_SIZE == 16
	return LIST_MUX16<T,L,Q,LOG2L> ( tab, sel );
#elif L_SIZE == 32
	return LIST_MUX32<T,L,Q,LOG2L> ( tab, sel );
#else
	return LIST_MUX64<T,L,Q,LOG2L> ( tab, sel );
#endif

}

//*************************************************************************//
//** 						UPDATE METRIC								 **//
//*************************************************************************//

template <int L, int Q, int LOG2L>
void UPDATE_METRIC (LLR_struct<1,Q,LOG2L> input[L], PS_struct<1,Q,LOG2L> output[2*L], sc_biguint<1> fb )
{
#pragma HLS INLINE
#pragma HLS ARRAY_PARTITION variable=input complete dim=1
#pragma HLS ARRAY_PARTITION variable=output complete dim=1

	if(fb == 1)
	{
		for(int i = 0; i < L; i++)
		{
		#pragma HLS UNROLL
			sc_biguint<1> sign_llr = qsign<Q> (input[i].llr);
			sc_bigint<Q> abs_llr = qabs<Q> (input[i].llr);

			// first path is created for bit = sign(LLR)
			output[2*i].bit 	= sign_llr;
			output[2*i].metric 	= input[i].metric; 	// no penalty for this metric
			output[2*i].path 	= input[i].path; 	// path is the same

			// second path is created for bit = ~sign(LLR)
			output[2*i+1].bit 	= ~sign_llr;
			output[2*i+1].metric= qadd<Q> ( input[i].metric , abs_llr) ; // penalty added for this metric
			output[2*i+1].path 	= input[i].path; 	// path is the same
		}
	}
	else
	{
		for(int i = 0; i < L; i++)
		{
		#pragma HLS UNROLL
			sc_biguint<1> sign_llr = qsign<Q> (input[i].llr);
			sc_bigint<Q> abs_llr = qabs<Q> (input[i].llr);

			output[i].bit = (sc_biguint< 1 >) 0;

			if( sign_llr != (sc_biguint< 1 >) 0 )
				output[i].metric = qadd<Q> ( input[i].metric , abs_llr);
			else
				output[i].metric = input[i].metric;

			output[i].path = input[i].path; 	// path is the same

		}

		// The rest of output array is fixed with a metric of +inf for sorting
		for(int i = 0; i < L; i++)
		{
		#pragma HLS UNROLL
			output[L+i].bit 	= (sc_biguint< 1 >) 0;
			output[L+i].metric 	= (sc_bigint<Q>)((sc_uint<1>)0, (sc_uint<Q-1>)0xFFFF); 	// +infinity
			output[L+i].path 	= (sc_biguint< 1 >) 0;
		}
	}
}


//*************************************************************************//
//** 			 		 SORTER ELEMENTS								 **//
//*************************************************************************//

// swap unit
template <int Q, int LOG2L>
void swap( PS_struct<1,Q,LOG2L> *xp, PS_struct<1,Q,LOG2L> *yp)
{
#pragma HLS INLINE
	PS_struct<1,Q,LOG2L> temp = *xp;
	*xp = *yp;
	*yp = temp;
}

// Compare and Select unit
template <int Q>
sc_uint<1> COMP( sc_biguint< Q > mx, sc_biguint< Q > my)
{
#pragma HLS INLINE
	if ( mx > my )
		return (sc_uint<1>) 1;
	else
		return (sc_uint<1>) 0;
}

template <int Q, int LOG2L>
sc_uint<1> COMP( PS_struct<1,Q,LOG2L> *xp, PS_struct<1,Q,LOG2L> *yp)
{
#pragma HLS INLINE
	if ( (xp->metric) > (yp->metric) )
		return (sc_uint<1>) 1;
	else
		return (sc_uint<1>) 0;
}

// Compare and Select unit
template <int Q, int LOG2L>
void CAS( PS_struct<1,Q,LOG2L> *xp, PS_struct<1,Q,LOG2L> *yp)
{
#pragma HLS INLINE
	if ( (xp->metric) > (yp->metric) )
	{
		PS_struct<1,Q,LOG2L> temp = *xp;
		*xp = *yp;
		*yp = temp;
	}
}


//*************************************************************************//
//** 						BUBBLE SORTER								 **//
//*************************************************************************//

template <int L, int Q, int LOG2L>
void BUBBLE_SORT (PS_struct<1,Q,LOG2L> input[2*L], PS_struct<1,Q,LOG2L> output[L])
{
#pragma HLS INLINE
#pragma HLS ARRAY_PARTITION variable=input complete dim=1
#pragma HLS ARRAY_PARTITION variable=output complete dim=1

	// Sort the 2*L input
	for (int i = 0; i < 2*L-1; i++){
	#pragma HLS UNROLL

	   for (int j = 0; j < 2*L-i-1; j++){
		#pragma HLS UNROLL
		   CAS<Q, LOG2L> (&input[j], &input[j+1]);
	   }
	}

   // output equal the L smallest input
   for (int i = 0; i < L; i++){
	#pragma HLS UNROLL

	   output[i] = input[i];
   }
}


//*************************************************************************//
//** 						CUSTOM SORTER								 **//
//*************************************************************************//

template <int Q>
void CUSTOM_SORT_L2 (PS_struct<1,Q,1> input[4], PS_struct<1,Q,1> output[2], sc_biguint<1> fb)
{
#pragma HLS INLINE
#pragma HLS ARRAY_PARTITION variable=input complete dim=1
#pragma HLS ARRAY_PARTITION variable=output complete dim=1

	// Adjust input when FB=0
	if (fb == 0){
		input[2] = input[1];
		input[1] = input[0];
	}

	CAS<Q, 1> (&input[1], &input[2]);

   // output equal the L smallest input
	if (fb == 0){
		output[0] = input[1];
		output[1] = input[2];
	}
	else{
		output[0] = input[0];
		output[1] = input[1];
	}
}

template <int Q>
void CUSTOM_SORT_L4 (PS_struct<1,Q,2> input[8], PS_struct<1,Q,2> output[4], sc_biguint<1> fb)
{
#pragma HLS INLINE
#pragma HLS ARRAY_PARTITION variable=input complete dim=1
#pragma HLS ARRAY_PARTITION variable=output complete dim=1

	// Adjust input when FB=0
	if (fb == 0){
		input[6] = input[3];
		input[5] = input[2];
		input[4] = input[1];
		input[3] = input[0];
		input[2].metric = 0;
		input[1].metric = 0;
	}

	CAS<Q, 2> (&input[1], &input[2]);
	CAS<Q, 2> (&input[3], &input[4]);
	CAS<Q, 2> (&input[5], &input[6]);

	CAS<Q, 2> (&input[2], &input[3]);
	CAS<Q, 2> (&input[4], &input[5]);

	CAS<Q, 2> (&input[3], &input[4]);
	CAS<Q, 2> (&input[5], &input[6]);

	CAS<Q, 2> (&input[4], &input[5]);

   // output equal the L smallest input
	if (fb == 0){
		output[0] = input[3];
		output[1] = input[4];
		output[2] = input[5];
		output[3] = input[6];
	}
	else{
		output[0] = input[0];
		output[1] = input[1];
		output[2] = input[2];
		output[3] = input[3];
	}
}

template <int Q>
void CUSTOM_SORT_L8 (PS_struct<1,Q,3> input[16], PS_struct<1,Q,3> output[8], sc_biguint<1> fb)
{
#pragma HLS INLINE
#pragma HLS ARRAY_PARTITION variable=input complete dim=1
#pragma HLS ARRAY_PARTITION variable=output complete dim=1

	// Adjust input when FB=0
	if (fb == 0){
		input[14].metric = MAX_VAL;
		input[13].metric = MAX_VAL;
		input[12] = input[7];
		input[11] = input[6];
		input[10] = input[5];
		input[9] = input[4];
		input[8] = input[3];
		input[7] = input[2];
		input[6] = input[1];
		input[5] = input[0];
		input[4].metric = 0;
		input[3].metric = 0;
		input[2].metric = 0;
		input[1].metric = 0;
	}

	CAS<Q, 3> (&input[1], &input[2]);
	CAS<Q, 3> (&input[3], &input[4]);
	CAS<Q, 3> (&input[5], &input[6]);
	CAS<Q, 3> (&input[7], &input[8]);
	CAS<Q, 3> (&input[9], &input[10]);
	CAS<Q, 3> (&input[11], &input[12]);
	CAS<Q, 3> (&input[13], &input[14]);

	CAS<Q, 3> (&input[2], &input[3]);
	CAS<Q, 3> (&input[4], &input[5]);
	CAS<Q, 3> (&input[6], &input[7]);
	CAS<Q, 3> (&input[8], &input[9]);
	CAS<Q, 3> (&input[10], &input[11]);
	CAS<Q, 3> (&input[12], &input[13]);

	CAS<Q, 3> (&input[3], &input[4]);
	CAS<Q, 3> (&input[5], &input[6]);
	CAS<Q, 3> (&input[7], &input[8]);
	CAS<Q, 3> (&input[9], &input[10]);
	CAS<Q, 3> (&input[11], &input[12]);

	CAS<Q, 3> (&input[4], &input[5]);
	CAS<Q, 3> (&input[6], &input[7]);
	CAS<Q, 3> (&input[8], &input[9]);
	CAS<Q, 3> (&input[10], &input[11]);

	CAS<Q, 3> (&input[5], &input[6]);
	CAS<Q, 3> (&input[7], &input[8]);
	CAS<Q, 3> (&input[9], &input[10]);
	CAS<Q, 3> (&input[11], &input[12]);

	CAS<Q, 3> (&input[6], &input[7]);
	CAS<Q, 3> (&input[8], &input[9]);
	CAS<Q, 3> (&input[10], &input[11]);

	CAS<Q, 3> (&input[5], &input[6]);
	CAS<Q, 3> (&input[7], &input[8]);
	CAS<Q, 3> (&input[9], &input[10]);
	CAS<Q, 3> (&input[11], &input[12]);

	CAS<Q, 3> (&input[6], &input[7]);
	CAS<Q, 3> (&input[8], &input[9]);
	CAS<Q, 3> (&input[10], &input[11]);

   // output equal the L smallest input
	if (fb == 0){
		output[0] = input[5];
		output[1] = input[6];
		output[2] = input[7];
		output[3] = input[8];
		output[4] = input[9];
		output[5] = input[10];
		output[6] = input[11];
		output[7] = input[12];
	}
	else{
		output[0] = input[0];
		output[1] = input[1];
		output[2] = input[2];
		output[3] = input[3];
		output[4] = input[4];
		output[5] = input[5];
		output[6] = input[6];
		output[7] = input[7];
	}
}

template <int Q>
void CUSTOM_SORT_L16 (PS_struct<1,Q,4> input[32], PS_struct<1,Q,4> output[16], sc_biguint<1> fb)
{
#pragma HLS INLINE
#pragma HLS ARRAY_PARTITION variable=input complete dim=1
#pragma HLS ARRAY_PARTITION variable=output complete dim=1

	// Adjust input when FB=0
	if (fb == 0){
		input[30].metric = MAX_VAL;
		input[29].metric = MAX_VAL;
		input[28].metric = MAX_VAL;
		input[27].metric = MAX_VAL;
		input[26].metric = MAX_VAL;
		input[25].metric = MAX_VAL;
		input[24] = input[15];
		input[23] = input[14];
		input[22] = input[13];
		input[21] = input[12];
		input[20] = input[11];
		input[19] = input[10];
		input[18] = input[9];
		input[17] = input[8];
		input[16] = input[7];
		input[15] = input[6];
		input[14] = input[5];
		input[13] = input[4];
		input[12] = input[3];
		input[11] = input[2];
		input[10] = input[1];
		input[9] = input[0];
		input[8].metric = 0;
		input[7].metric = 0;
		input[6].metric = 0;
		input[5].metric = 0;
		input[4].metric = 0;
		input[3].metric = 0;
		input[2].metric = 0;
		input[1].metric = 0;
	}

	//1
	CAS<Q, 4> (&input[1], &input[2]);
	CAS<Q, 4> (&input[3], &input[4]);
	CAS<Q, 4> (&input[5], &input[6]);
	CAS<Q, 4> (&input[7], &input[8]);
	CAS<Q, 4> (&input[9], &input[10]);
	CAS<Q, 4> (&input[11], &input[12]);
	CAS<Q, 4> (&input[13], &input[14]);
	CAS<Q, 4> (&input[15], &input[16]);
	CAS<Q, 4> (&input[17], &input[18]);
	CAS<Q, 4> (&input[19], &input[20]);
	CAS<Q, 4> (&input[21], &input[22]);
	CAS<Q, 4> (&input[23], &input[24]);
	CAS<Q, 4> (&input[25], &input[26]);
	CAS<Q, 4> (&input[27], &input[28]);
	CAS<Q, 4> (&input[29], &input[30]);
	//2
	CAS<Q, 4> (&input[2], &input[3]);
	CAS<Q, 4> (&input[4], &input[5]);
	CAS<Q, 4> (&input[6], &input[7]);
	CAS<Q, 4> (&input[8], &input[9]);
	CAS<Q, 4> (&input[10], &input[11]);
	CAS<Q, 4> (&input[12], &input[13]);
	CAS<Q, 4> (&input[14], &input[15]);
	CAS<Q, 4> (&input[16], &input[17]);
	CAS<Q, 4> (&input[18], &input[19]);
	CAS<Q, 4> (&input[20], &input[21]);
	CAS<Q, 4> (&input[22], &input[23]);
	CAS<Q, 4> (&input[24], &input[25]);
	CAS<Q, 4> (&input[26], &input[27]);
	CAS<Q, 4> (&input[28], &input[29]);
	//3
	CAS<Q, 4> (&input[3], &input[4]);
	CAS<Q, 4> (&input[5], &input[6]);
	CAS<Q, 4> (&input[7], &input[8]);
	CAS<Q, 4> (&input[9], &input[10]);
	CAS<Q, 4> (&input[11], &input[12]);
	CAS<Q, 4> (&input[13], &input[14]);
	CAS<Q, 4> (&input[15], &input[16]);
	CAS<Q, 4> (&input[17], &input[18]);
	CAS<Q, 4> (&input[19], &input[20]);
	CAS<Q, 4> (&input[21], &input[22]);
	CAS<Q, 4> (&input[23], &input[24]);
	CAS<Q, 4> (&input[25], &input[26]);
	CAS<Q, 4> (&input[27], &input[28]);
	//4
	CAS<Q, 4> (&input[4], &input[5]);
	CAS<Q, 4> (&input[6], &input[7]);
	CAS<Q, 4> (&input[8], &input[9]);
	CAS<Q, 4> (&input[10], &input[11]);
	CAS<Q, 4> (&input[12], &input[13]);
	CAS<Q, 4> (&input[14], &input[15]);
	CAS<Q, 4> (&input[16], &input[17]);
	CAS<Q, 4> (&input[18], &input[19]);
	CAS<Q, 4> (&input[20], &input[21]);
	CAS<Q, 4> (&input[22], &input[23]);
	CAS<Q, 4> (&input[24], &input[25]);
	CAS<Q, 4> (&input[26], &input[27]);
	//5
	CAS<Q, 4> (&input[5], &input[6]);
	CAS<Q, 4> (&input[7], &input[8]);
	CAS<Q, 4> (&input[9], &input[10]);
	CAS<Q, 4> (&input[11], &input[12]);
	CAS<Q, 4> (&input[13], &input[14]);
	CAS<Q, 4> (&input[15], &input[16]);
	CAS<Q, 4> (&input[17], &input[18]);
	CAS<Q, 4> (&input[19], &input[20]);
	CAS<Q, 4> (&input[21], &input[22]);
	CAS<Q, 4> (&input[23], &input[24]);
	CAS<Q, 4> (&input[25], &input[26]);
	//6
	CAS<Q, 4> (&input[6], &input[7]);
	CAS<Q, 4> (&input[8], &input[9]);
	CAS<Q, 4> (&input[10], &input[11]);
	CAS<Q, 4> (&input[12], &input[13]);
	CAS<Q, 4> (&input[14], &input[15]);
	CAS<Q, 4> (&input[16], &input[17]);
	CAS<Q, 4> (&input[18], &input[19]);
	CAS<Q, 4> (&input[20], &input[21]);
	CAS<Q, 4> (&input[22], &input[23]);
	CAS<Q, 4> (&input[24], &input[25]);
	//7
	CAS<Q, 4> (&input[7], &input[8]);
	CAS<Q, 4> (&input[9], &input[10]);
	CAS<Q, 4> (&input[11], &input[12]);
	CAS<Q, 4> (&input[13], &input[14]);
	CAS<Q, 4> (&input[15], &input[16]);
	CAS<Q, 4> (&input[17], &input[18]);
	CAS<Q, 4> (&input[19], &input[20]);
	CAS<Q, 4> (&input[21], &input[22]);
	CAS<Q, 4> (&input[23], &input[24]);
	//8
	CAS<Q, 4> (&input[8], &input[9]);
	CAS<Q, 4> (&input[10], &input[11]);
	CAS<Q, 4> (&input[12], &input[13]);
	CAS<Q, 4> (&input[14], &input[15]);
	CAS<Q, 4> (&input[16], &input[17]);
	CAS<Q, 4> (&input[18], &input[19]);
	CAS<Q, 4> (&input[20], &input[21]);
	CAS<Q, 4> (&input[22], &input[23]);
	//9
	CAS<Q, 4> (&input[9], &input[10]);
	CAS<Q, 4> (&input[11], &input[12]);
	CAS<Q, 4> (&input[13], &input[14]);
	CAS<Q, 4> (&input[15], &input[16]);
	CAS<Q, 4> (&input[17], &input[18]);
	CAS<Q, 4> (&input[19], &input[20]);
	CAS<Q, 4> (&input[21], &input[22]);
	CAS<Q, 4> (&input[23], &input[24]);
	//10
	CAS<Q, 4> (&input[10], &input[11]);
	CAS<Q, 4> (&input[12], &input[13]);
	CAS<Q, 4> (&input[14], &input[15]);
	CAS<Q, 4> (&input[16], &input[17]);
	CAS<Q, 4> (&input[18], &input[19]);
	CAS<Q, 4> (&input[20], &input[21]);
	CAS<Q, 4> (&input[22], &input[23]);
	//11
	CAS<Q, 4> (&input[9], &input[10]);
	CAS<Q, 4> (&input[11], &input[12]);
	CAS<Q, 4> (&input[13], &input[14]);
	CAS<Q, 4> (&input[15], &input[16]);
	CAS<Q, 4> (&input[17], &input[18]);
	CAS<Q, 4> (&input[19], &input[20]);
	CAS<Q, 4> (&input[21], &input[22]);
	CAS<Q, 4> (&input[23], &input[24]);
	//12
	CAS<Q, 4> (&input[10], &input[11]);
	CAS<Q, 4> (&input[12], &input[13]);
	CAS<Q, 4> (&input[14], &input[15]);
	CAS<Q, 4> (&input[16], &input[17]);
	CAS<Q, 4> (&input[18], &input[19]);
	CAS<Q, 4> (&input[20], &input[21]);
	CAS<Q, 4> (&input[22], &input[23]);
	//13
	CAS<Q, 4> (&input[9], &input[10]);
	CAS<Q, 4> (&input[11], &input[12]);
	CAS<Q, 4> (&input[13], &input[14]);
	CAS<Q, 4> (&input[15], &input[16]);
	CAS<Q, 4> (&input[17], &input[18]);
	CAS<Q, 4> (&input[19], &input[20]);
	CAS<Q, 4> (&input[21], &input[22]);
	CAS<Q, 4> (&input[23], &input[24]);
	//14
	CAS<Q, 4> (&input[10], &input[11]);
	CAS<Q, 4> (&input[12], &input[13]);
	CAS<Q, 4> (&input[14], &input[15]);
	CAS<Q, 4> (&input[16], &input[17]);
	CAS<Q, 4> (&input[18], &input[19]);
	CAS<Q, 4> (&input[20], &input[21]);
	CAS<Q, 4> (&input[22], &input[23]);
	//15
	CAS<Q, 4> (&input[9], &input[10]);
	CAS<Q, 4> (&input[11], &input[12]);
	CAS<Q, 4> (&input[13], &input[14]);
	CAS<Q, 4> (&input[15], &input[16]);
	CAS<Q, 4> (&input[17], &input[18]);
	CAS<Q, 4> (&input[19], &input[20]);
	CAS<Q, 4> (&input[21], &input[22]);
	CAS<Q, 4> (&input[23], &input[24]);
	//16
	CAS<Q, 4> (&input[10], &input[11]);
	CAS<Q, 4> (&input[12], &input[13]);
	CAS<Q, 4> (&input[14], &input[15]);
	CAS<Q, 4> (&input[16], &input[17]);
	CAS<Q, 4> (&input[18], &input[19]);
	CAS<Q, 4> (&input[20], &input[21]);
	CAS<Q, 4> (&input[22], &input[23]);

   // output equal the L smallest input
	if (fb == 0){
		output[0] = input[9];
		output[1] = input[10];
		output[2] = input[11];
		output[3] = input[12];
		output[4] = input[13];
		output[5] = input[14];
		output[6] = input[15];
		output[7] = input[16];
		output[8] = input[17];
		output[9] = input[18];
		output[10] = input[19];
		output[11] = input[20];
		output[12] = input[21];
		output[13] = input[22];
		output[14] = input[23];
		output[15] = input[24];
	}
	else{
		output[0] = input[0];
		output[1] = input[1];
		output[2] = input[2];
		output[3] = input[3];
		output[4] = input[4];
		output[5] = input[5];
		output[6] = input[6];
		output[7] = input[7];
		output[8] = input[8];
		output[9] = input[9];
		output[10] = input[10];
		output[11] = input[11];
		output[12] = input[12];
		output[13] = input[13];
		output[14] = input[14];
		output[15] = input[15];
	}
}

template <int Q>
void CUSTOM_SORT_L32 (PS_struct<1,Q,5> input[64], PS_struct<1,Q,5> output[32], sc_biguint<1> fb)
{
#pragma HLS INLINE
#pragma HLS ARRAY_PARTITION variable=input complete dim=1
#pragma HLS ARRAY_PARTITION variable=output complete dim=1

	// Adjust input when FB=0
	if (fb == 0){
		input[63].metric = MAX_VAL;
		input[62].metric = MAX_VAL;
		input[61].metric = MAX_VAL;
		input[60].metric = MAX_VAL;
		input[59].metric = MAX_VAL;
		input[58].metric = MAX_VAL;
		input[57].metric = MAX_VAL;
		input[56].metric = MAX_VAL;
		input[55].metric = MAX_VAL;
		input[54].metric = MAX_VAL;
		input[53].metric = MAX_VAL;
		input[52].metric = MAX_VAL;
		input[51].metric = MAX_VAL;
		input[50].metric = MAX_VAL;
		input[49].metric = MAX_VAL;
		input[48] = input[31];
		input[47] = input[30];
		input[46] = input[29];
		input[45] = input[28];
		input[44] = input[27];
		input[43] = input[26];
		input[42] = input[25];
		input[41] = input[24];
		input[40] = input[23];
		input[39] = input[22];
		input[38] = input[21];
		input[37] = input[20];
		input[36] = input[19];
		input[35] = input[18];
		input[34] = input[17];
		input[33] = input[16];
		input[32] = input[15];
		input[31] = input[14];
		input[30] = input[13];
		input[29] = input[12];
		input[28] = input[11];
		input[27] = input[10];
		input[26] = input[9];
		input[25] = input[8];
		input[24] = input[7];
		input[23] = input[6];
		input[22] = input[5];
		input[21] = input[4];
		input[20] = input[3];
		input[19] = input[2];
		input[18] = input[1];
		input[17] = input[0];
		input[16].metric = 0;
		input[15].metric = 0;
		input[14].metric = 0;
		input[13].metric = 0;
		input[12].metric = 0;
		input[11].metric = 0;
		input[10].metric = 0;
		input[9].metric = 0;
		input[8].metric = 0;
		input[7].metric = 0;
		input[6].metric = 0;
		input[5].metric = 0;
		input[4].metric = 0;
		input[3].metric = 0;
		input[2].metric = 0;
		input[1].metric = 0;
	}

	//1
	CAS<Q, 5> (&input[1], &input[2]);
	CAS<Q, 5> (&input[3], &input[4]);
	CAS<Q, 5> (&input[5], &input[6]);
	CAS<Q, 5> (&input[7], &input[8]);
	CAS<Q, 5> (&input[9], &input[10]);
	CAS<Q, 5> (&input[11], &input[12]);
	CAS<Q, 5> (&input[13], &input[14]);
	CAS<Q, 5> (&input[15], &input[16]);
	CAS<Q, 5> (&input[17], &input[18]);
	CAS<Q, 5> (&input[19], &input[20]);
	CAS<Q, 5> (&input[21], &input[22]);
	CAS<Q, 5> (&input[23], &input[24]);
	CAS<Q, 5> (&input[25], &input[26]);
	CAS<Q, 5> (&input[27], &input[28]);
	CAS<Q, 5> (&input[29], &input[30]);
	CAS<Q, 5> (&input[31], &input[32]);
	CAS<Q, 5> (&input[33], &input[34]);
	CAS<Q, 5> (&input[35], &input[36]);
	CAS<Q, 5> (&input[37], &input[38]);
	CAS<Q, 5> (&input[39], &input[40]);
	CAS<Q, 5> (&input[41], &input[42]);
	CAS<Q, 5> (&input[43], &input[44]);
	CAS<Q, 5> (&input[45], &input[46]);
	CAS<Q, 5> (&input[47], &input[48]);
	CAS<Q, 5> (&input[49], &input[50]);
	CAS<Q, 5> (&input[51], &input[52]);
	CAS<Q, 5> (&input[53], &input[54]);
	CAS<Q, 5> (&input[55], &input[56]);
	CAS<Q, 5> (&input[57], &input[58]);
	CAS<Q, 5> (&input[59], &input[60]);
	CAS<Q, 5> (&input[61], &input[62]);
	//2
	CAS<Q, 5> (&input[2], &input[3]);
	CAS<Q, 5> (&input[4], &input[5]);
	CAS<Q, 5> (&input[6], &input[7]);
	CAS<Q, 5> (&input[8], &input[9]);
	CAS<Q, 5> (&input[10], &input[11]);
	CAS<Q, 5> (&input[12], &input[13]);
	CAS<Q, 5> (&input[14], &input[15]);
	CAS<Q, 5> (&input[16], &input[17]);
	CAS<Q, 5> (&input[18], &input[19]);
	CAS<Q, 5> (&input[20], &input[21]);
	CAS<Q, 5> (&input[22], &input[23]);
	CAS<Q, 5> (&input[24], &input[25]);
	CAS<Q, 5> (&input[26], &input[27]);
	CAS<Q, 5> (&input[28], &input[29]);
	CAS<Q, 5> (&input[30], &input[31]);
	CAS<Q, 5> (&input[32], &input[33]);
	CAS<Q, 5> (&input[34], &input[35]);
	CAS<Q, 5> (&input[36], &input[37]);
	CAS<Q, 5> (&input[38], &input[39]);
	CAS<Q, 5> (&input[40], &input[41]);
	CAS<Q, 5> (&input[42], &input[43]);
	CAS<Q, 5> (&input[44], &input[45]);
	CAS<Q, 5> (&input[46], &input[47]);
	CAS<Q, 5> (&input[48], &input[49]);
	CAS<Q, 5> (&input[50], &input[51]);
	CAS<Q, 5> (&input[52], &input[53]);
	CAS<Q, 5> (&input[54], &input[55]);
	CAS<Q, 5> (&input[56], &input[57]);
	CAS<Q, 5> (&input[58], &input[59]);
	CAS<Q, 5> (&input[60], &input[61]);
	//3
	CAS<Q, 5> (&input[3], &input[4]);
	CAS<Q, 5> (&input[5], &input[6]);
	CAS<Q, 5> (&input[7], &input[8]);
	CAS<Q, 5> (&input[9], &input[10]);
	CAS<Q, 5> (&input[11], &input[12]);
	CAS<Q, 5> (&input[13], &input[14]);
	CAS<Q, 5> (&input[15], &input[16]);
	CAS<Q, 5> (&input[17], &input[18]);
	CAS<Q, 5> (&input[19], &input[20]);
	CAS<Q, 5> (&input[21], &input[22]);
	CAS<Q, 5> (&input[23], &input[24]);
	CAS<Q, 5> (&input[25], &input[26]);
	CAS<Q, 5> (&input[27], &input[28]);
	CAS<Q, 5> (&input[29], &input[30]);
	CAS<Q, 5> (&input[31], &input[32]);
	CAS<Q, 5> (&input[33], &input[34]);
	CAS<Q, 5> (&input[35], &input[36]);
	CAS<Q, 5> (&input[37], &input[38]);
	CAS<Q, 5> (&input[39], &input[40]);
	CAS<Q, 5> (&input[41], &input[42]);
	CAS<Q, 5> (&input[43], &input[44]);
	CAS<Q, 5> (&input[45], &input[46]);
	CAS<Q, 5> (&input[47], &input[48]);
	CAS<Q, 5> (&input[49], &input[50]);
	CAS<Q, 5> (&input[51], &input[52]);
	CAS<Q, 5> (&input[53], &input[54]);
	CAS<Q, 5> (&input[55], &input[56]);
	CAS<Q, 5> (&input[57], &input[58]);
	CAS<Q, 5> (&input[59], &input[60]);
	//4
	CAS<Q, 5> (&input[4], &input[5]);
	CAS<Q, 5> (&input[6], &input[7]);
	CAS<Q, 5> (&input[8], &input[9]);
	CAS<Q, 5> (&input[10], &input[11]);
	CAS<Q, 5> (&input[12], &input[13]);
	CAS<Q, 5> (&input[14], &input[15]);
	CAS<Q, 5> (&input[16], &input[17]);
	CAS<Q, 5> (&input[18], &input[19]);
	CAS<Q, 5> (&input[20], &input[21]);
	CAS<Q, 5> (&input[22], &input[23]);
	CAS<Q, 5> (&input[24], &input[25]);
	CAS<Q, 5> (&input[26], &input[27]);
	CAS<Q, 5> (&input[28], &input[29]);
	CAS<Q, 5> (&input[30], &input[31]);
	CAS<Q, 5> (&input[32], &input[33]);
	CAS<Q, 5> (&input[34], &input[35]);
	CAS<Q, 5> (&input[36], &input[37]);
	CAS<Q, 5> (&input[38], &input[39]);
	CAS<Q, 5> (&input[40], &input[41]);
	CAS<Q, 5> (&input[42], &input[43]);
	CAS<Q, 5> (&input[44], &input[45]);
	CAS<Q, 5> (&input[46], &input[47]);
	CAS<Q, 5> (&input[48], &input[49]);
	CAS<Q, 5> (&input[50], &input[51]);
	CAS<Q, 5> (&input[52], &input[53]);
	CAS<Q, 5> (&input[54], &input[55]);
	CAS<Q, 5> (&input[56], &input[57]);
	CAS<Q, 5> (&input[58], &input[59]);
	//5
	CAS<Q, 5> (&input[5], &input[6]);
	CAS<Q, 5> (&input[7], &input[8]);
	CAS<Q, 5> (&input[9], &input[10]);
	CAS<Q, 5> (&input[11], &input[12]);
	CAS<Q, 5> (&input[13], &input[14]);
	CAS<Q, 5> (&input[15], &input[16]);
	CAS<Q, 5> (&input[17], &input[18]);
	CAS<Q, 5> (&input[19], &input[20]);
	CAS<Q, 5> (&input[21], &input[22]);
	CAS<Q, 5> (&input[23], &input[24]);
	CAS<Q, 5> (&input[25], &input[26]);
	CAS<Q, 5> (&input[27], &input[28]);
	CAS<Q, 5> (&input[29], &input[30]);
	CAS<Q, 5> (&input[31], &input[32]);
	CAS<Q, 5> (&input[33], &input[34]);
	CAS<Q, 5> (&input[35], &input[36]);
	CAS<Q, 5> (&input[37], &input[38]);
	CAS<Q, 5> (&input[39], &input[40]);
	CAS<Q, 5> (&input[41], &input[42]);
	CAS<Q, 5> (&input[43], &input[44]);
	CAS<Q, 5> (&input[45], &input[46]);
	CAS<Q, 5> (&input[47], &input[48]);
	CAS<Q, 5> (&input[49], &input[50]);
	CAS<Q, 5> (&input[51], &input[52]);
	CAS<Q, 5> (&input[53], &input[54]);
	CAS<Q, 5> (&input[55], &input[56]);
	CAS<Q, 5> (&input[57], &input[58]);
	//6
	CAS<Q, 5> (&input[6], &input[7]);
	CAS<Q, 5> (&input[8], &input[9]);
	CAS<Q, 5> (&input[10], &input[11]);
	CAS<Q, 5> (&input[12], &input[13]);
	CAS<Q, 5> (&input[14], &input[15]);
	CAS<Q, 5> (&input[16], &input[17]);
	CAS<Q, 5> (&input[18], &input[19]);
	CAS<Q, 5> (&input[20], &input[21]);
	CAS<Q, 5> (&input[22], &input[23]);
	CAS<Q, 5> (&input[24], &input[25]);
	CAS<Q, 5> (&input[26], &input[27]);
	CAS<Q, 5> (&input[28], &input[29]);
	CAS<Q, 5> (&input[30], &input[31]);
	CAS<Q, 5> (&input[32], &input[33]);
	CAS<Q, 5> (&input[34], &input[35]);
	CAS<Q, 5> (&input[36], &input[37]);
	CAS<Q, 5> (&input[38], &input[39]);
	CAS<Q, 5> (&input[40], &input[41]);
	CAS<Q, 5> (&input[42], &input[43]);
	CAS<Q, 5> (&input[44], &input[45]);
	CAS<Q, 5> (&input[46], &input[47]);
	CAS<Q, 5> (&input[48], &input[49]);
	CAS<Q, 5> (&input[50], &input[51]);
	CAS<Q, 5> (&input[52], &input[53]);
	CAS<Q, 5> (&input[54], &input[55]);
	CAS<Q, 5> (&input[56], &input[57]);
	//7
	CAS<Q, 5> (&input[7], &input[8]);
	CAS<Q, 5> (&input[9], &input[10]);
	CAS<Q, 5> (&input[11], &input[12]);
	CAS<Q, 5> (&input[13], &input[14]);
	CAS<Q, 5> (&input[15], &input[16]);
	CAS<Q, 5> (&input[17], &input[18]);
	CAS<Q, 5> (&input[19], &input[20]);
	CAS<Q, 5> (&input[21], &input[22]);
	CAS<Q, 5> (&input[23], &input[24]);
	CAS<Q, 5> (&input[25], &input[26]);
	CAS<Q, 5> (&input[27], &input[28]);
	CAS<Q, 5> (&input[29], &input[30]);
	CAS<Q, 5> (&input[31], &input[32]);
	CAS<Q, 5> (&input[33], &input[34]);
	CAS<Q, 5> (&input[35], &input[36]);
	CAS<Q, 5> (&input[37], &input[38]);
	CAS<Q, 5> (&input[39], &input[40]);
	CAS<Q, 5> (&input[41], &input[42]);
	CAS<Q, 5> (&input[43], &input[44]);
	CAS<Q, 5> (&input[45], &input[46]);
	CAS<Q, 5> (&input[47], &input[48]);
	CAS<Q, 5> (&input[49], &input[50]);
	CAS<Q, 5> (&input[51], &input[52]);
	CAS<Q, 5> (&input[53], &input[54]);
	CAS<Q, 5> (&input[55], &input[56]);
	//8
	CAS<Q, 5> (&input[8], &input[9]);
	CAS<Q, 5> (&input[10], &input[11]);
	CAS<Q, 5> (&input[12], &input[13]);
	CAS<Q, 5> (&input[14], &input[15]);
	CAS<Q, 5> (&input[16], &input[17]);
	CAS<Q, 5> (&input[18], &input[19]);
	CAS<Q, 5> (&input[20], &input[21]);
	CAS<Q, 5> (&input[22], &input[23]);
	CAS<Q, 5> (&input[24], &input[25]);
	CAS<Q, 5> (&input[26], &input[27]);
	CAS<Q, 5> (&input[28], &input[29]);
	CAS<Q, 5> (&input[30], &input[31]);
	CAS<Q, 5> (&input[32], &input[33]);
	CAS<Q, 5> (&input[34], &input[35]);
	CAS<Q, 5> (&input[36], &input[37]);
	CAS<Q, 5> (&input[38], &input[39]);
	CAS<Q, 5> (&input[40], &input[41]);
	CAS<Q, 5> (&input[42], &input[43]);
	CAS<Q, 5> (&input[44], &input[45]);
	CAS<Q, 5> (&input[46], &input[47]);
	CAS<Q, 5> (&input[48], &input[49]);
	CAS<Q, 5> (&input[50], &input[51]);
	CAS<Q, 5> (&input[52], &input[53]);
	CAS<Q, 5> (&input[54], &input[55]);
	//9
	CAS<Q, 5> (&input[9], &input[10]);
	CAS<Q, 5> (&input[11], &input[12]);
	CAS<Q, 5> (&input[13], &input[14]);
	CAS<Q, 5> (&input[15], &input[16]);
	CAS<Q, 5> (&input[17], &input[18]);
	CAS<Q, 5> (&input[19], &input[20]);
	CAS<Q, 5> (&input[21], &input[22]);
	CAS<Q, 5> (&input[23], &input[24]);
	CAS<Q, 5> (&input[25], &input[26]);
	CAS<Q, 5> (&input[27], &input[28]);
	CAS<Q, 5> (&input[29], &input[30]);
	CAS<Q, 5> (&input[31], &input[32]);
	CAS<Q, 5> (&input[33], &input[34]);
	CAS<Q, 5> (&input[35], &input[36]);
	CAS<Q, 5> (&input[37], &input[38]);
	CAS<Q, 5> (&input[39], &input[40]);
	CAS<Q, 5> (&input[41], &input[42]);
	CAS<Q, 5> (&input[43], &input[44]);
	CAS<Q, 5> (&input[45], &input[46]);
	CAS<Q, 5> (&input[47], &input[48]);
	CAS<Q, 5> (&input[49], &input[50]);
	CAS<Q, 5> (&input[51], &input[52]);
	CAS<Q, 5> (&input[53], &input[54]);
	//10
	CAS<Q, 5> (&input[10], &input[11]);
	CAS<Q, 5> (&input[12], &input[13]);
	CAS<Q, 5> (&input[14], &input[15]);
	CAS<Q, 5> (&input[16], &input[17]);
	CAS<Q, 5> (&input[18], &input[19]);
	CAS<Q, 5> (&input[20], &input[21]);
	CAS<Q, 5> (&input[22], &input[23]);
	CAS<Q, 5> (&input[24], &input[25]);
	CAS<Q, 5> (&input[26], &input[27]);
	CAS<Q, 5> (&input[28], &input[29]);
	CAS<Q, 5> (&input[30], &input[31]);
	CAS<Q, 5> (&input[32], &input[33]);
	CAS<Q, 5> (&input[34], &input[35]);
	CAS<Q, 5> (&input[36], &input[37]);
	CAS<Q, 5> (&input[38], &input[39]);
	CAS<Q, 5> (&input[40], &input[41]);
	CAS<Q, 5> (&input[42], &input[43]);
	CAS<Q, 5> (&input[44], &input[45]);
	CAS<Q, 5> (&input[46], &input[47]);
	CAS<Q, 5> (&input[48], &input[49]);
	CAS<Q, 5> (&input[50], &input[51]);
	CAS<Q, 5> (&input[52], &input[53]);
	//11
	CAS<Q, 5> (&input[11], &input[12]);
	CAS<Q, 5> (&input[13], &input[14]);
	CAS<Q, 5> (&input[15], &input[16]);
	CAS<Q, 5> (&input[17], &input[18]);
	CAS<Q, 5> (&input[19], &input[20]);
	CAS<Q, 5> (&input[21], &input[22]);
	CAS<Q, 5> (&input[23], &input[24]);
	CAS<Q, 5> (&input[25], &input[26]);
	CAS<Q, 5> (&input[27], &input[28]);
	CAS<Q, 5> (&input[29], &input[30]);
	CAS<Q, 5> (&input[31], &input[32]);
	CAS<Q, 5> (&input[33], &input[34]);
	CAS<Q, 5> (&input[35], &input[36]);
	CAS<Q, 5> (&input[37], &input[38]);
	CAS<Q, 5> (&input[39], &input[40]);
	CAS<Q, 5> (&input[41], &input[42]);
	CAS<Q, 5> (&input[43], &input[44]);
	CAS<Q, 5> (&input[45], &input[46]);
	CAS<Q, 5> (&input[47], &input[48]);
	CAS<Q, 5> (&input[49], &input[50]);
	CAS<Q, 5> (&input[51], &input[52]);
	//12
	CAS<Q, 5> (&input[12], &input[13]);
	CAS<Q, 5> (&input[14], &input[15]);
	CAS<Q, 5> (&input[16], &input[17]);
	CAS<Q, 5> (&input[18], &input[19]);
	CAS<Q, 5> (&input[20], &input[21]);
	CAS<Q, 5> (&input[22], &input[23]);
	CAS<Q, 5> (&input[24], &input[25]);
	CAS<Q, 5> (&input[26], &input[27]);
	CAS<Q, 5> (&input[28], &input[29]);
	CAS<Q, 5> (&input[30], &input[31]);
	CAS<Q, 5> (&input[32], &input[33]);
	CAS<Q, 5> (&input[34], &input[35]);
	CAS<Q, 5> (&input[36], &input[37]);
	CAS<Q, 5> (&input[38], &input[39]);
	CAS<Q, 5> (&input[40], &input[41]);
	CAS<Q, 5> (&input[42], &input[43]);
	CAS<Q, 5> (&input[44], &input[45]);
	CAS<Q, 5> (&input[46], &input[47]);
	CAS<Q, 5> (&input[48], &input[49]);
	CAS<Q, 5> (&input[50], &input[51]);
	//13
	CAS<Q, 5> (&input[13], &input[14]);
	CAS<Q, 5> (&input[15], &input[16]);
	CAS<Q, 5> (&input[17], &input[18]);
	CAS<Q, 5> (&input[19], &input[20]);
	CAS<Q, 5> (&input[21], &input[22]);
	CAS<Q, 5> (&input[23], &input[24]);
	CAS<Q, 5> (&input[25], &input[26]);
	CAS<Q, 5> (&input[27], &input[28]);
	CAS<Q, 5> (&input[29], &input[30]);
	CAS<Q, 5> (&input[31], &input[32]);
	CAS<Q, 5> (&input[33], &input[34]);
	CAS<Q, 5> (&input[35], &input[36]);
	CAS<Q, 5> (&input[37], &input[38]);
	CAS<Q, 5> (&input[39], &input[40]);
	CAS<Q, 5> (&input[41], &input[42]);
	CAS<Q, 5> (&input[43], &input[44]);
	CAS<Q, 5> (&input[45], &input[46]);
	CAS<Q, 5> (&input[47], &input[48]);
	CAS<Q, 5> (&input[49], &input[50]);
	//14
	CAS<Q, 5> (&input[14], &input[15]);
	CAS<Q, 5> (&input[16], &input[17]);
	CAS<Q, 5> (&input[18], &input[19]);
	CAS<Q, 5> (&input[20], &input[21]);
	CAS<Q, 5> (&input[22], &input[23]);
	CAS<Q, 5> (&input[24], &input[25]);
	CAS<Q, 5> (&input[26], &input[27]);
	CAS<Q, 5> (&input[28], &input[29]);
	CAS<Q, 5> (&input[30], &input[31]);
	CAS<Q, 5> (&input[32], &input[33]);
	CAS<Q, 5> (&input[34], &input[35]);
	CAS<Q, 5> (&input[36], &input[37]);
	CAS<Q, 5> (&input[38], &input[39]);
	CAS<Q, 5> (&input[40], &input[41]);
	CAS<Q, 5> (&input[42], &input[43]);
	CAS<Q, 5> (&input[44], &input[45]);
	CAS<Q, 5> (&input[46], &input[47]);
	CAS<Q, 5> (&input[48], &input[49]);
	//15
	CAS<Q, 5> (&input[15], &input[16]);
	CAS<Q, 5> (&input[17], &input[18]);
	CAS<Q, 5> (&input[19], &input[20]);
	CAS<Q, 5> (&input[21], &input[22]);
	CAS<Q, 5> (&input[23], &input[24]);
	CAS<Q, 5> (&input[25], &input[26]);
	CAS<Q, 5> (&input[27], &input[28]);
	CAS<Q, 5> (&input[29], &input[30]);
	CAS<Q, 5> (&input[31], &input[32]);
	CAS<Q, 5> (&input[33], &input[34]);
	CAS<Q, 5> (&input[35], &input[36]);
	CAS<Q, 5> (&input[37], &input[38]);
	CAS<Q, 5> (&input[39], &input[40]);
	CAS<Q, 5> (&input[41], &input[42]);
	CAS<Q, 5> (&input[43], &input[44]);
	CAS<Q, 5> (&input[45], &input[46]);
	CAS<Q, 5> (&input[47], &input[48]);
	//16
	CAS<Q, 5> (&input[16], &input[17]);
	CAS<Q, 5> (&input[18], &input[19]);
	CAS<Q, 5> (&input[20], &input[21]);
	CAS<Q, 5> (&input[22], &input[23]);
	CAS<Q, 5> (&input[24], &input[25]);
	CAS<Q, 5> (&input[26], &input[27]);
	CAS<Q, 5> (&input[28], &input[29]);
	CAS<Q, 5> (&input[30], &input[31]);
	CAS<Q, 5> (&input[32], &input[33]);
	CAS<Q, 5> (&input[34], &input[35]);
	CAS<Q, 5> (&input[36], &input[37]);
	CAS<Q, 5> (&input[38], &input[39]);
	CAS<Q, 5> (&input[40], &input[41]);
	CAS<Q, 5> (&input[42], &input[43]);
	CAS<Q, 5> (&input[44], &input[45]);
	CAS<Q, 5> (&input[46], &input[47]);
	//17
	CAS<Q, 5> (&input[17], &input[18]);
	CAS<Q, 5> (&input[19], &input[20]);
	CAS<Q, 5> (&input[21], &input[22]);
	CAS<Q, 5> (&input[23], &input[24]);
	CAS<Q, 5> (&input[25], &input[26]);
	CAS<Q, 5> (&input[27], &input[28]);
	CAS<Q, 5> (&input[29], &input[30]);
	CAS<Q, 5> (&input[31], &input[32]);
	CAS<Q, 5> (&input[33], &input[34]);
	CAS<Q, 5> (&input[35], &input[36]);
	CAS<Q, 5> (&input[37], &input[38]);
	CAS<Q, 5> (&input[39], &input[40]);
	CAS<Q, 5> (&input[41], &input[42]);
	CAS<Q, 5> (&input[43], &input[44]);
	CAS<Q, 5> (&input[45], &input[46]);
	CAS<Q, 5> (&input[47], &input[48]);
	//18
	CAS<Q, 5> (&input[18], &input[19]);
	CAS<Q, 5> (&input[20], &input[21]);
	CAS<Q, 5> (&input[22], &input[23]);
	CAS<Q, 5> (&input[24], &input[25]);
	CAS<Q, 5> (&input[26], &input[27]);
	CAS<Q, 5> (&input[28], &input[29]);
	CAS<Q, 5> (&input[30], &input[31]);
	CAS<Q, 5> (&input[32], &input[33]);
	CAS<Q, 5> (&input[34], &input[35]);
	CAS<Q, 5> (&input[36], &input[37]);
	CAS<Q, 5> (&input[38], &input[39]);
	CAS<Q, 5> (&input[40], &input[41]);
	CAS<Q, 5> (&input[42], &input[43]);
	CAS<Q, 5> (&input[44], &input[45]);
	CAS<Q, 5> (&input[46], &input[47]);
	//19
	CAS<Q, 5> (&input[17], &input[18]);
	CAS<Q, 5> (&input[19], &input[20]);
	CAS<Q, 5> (&input[21], &input[22]);
	CAS<Q, 5> (&input[23], &input[24]);
	CAS<Q, 5> (&input[25], &input[26]);
	CAS<Q, 5> (&input[27], &input[28]);
	CAS<Q, 5> (&input[29], &input[30]);
	CAS<Q, 5> (&input[31], &input[32]);
	CAS<Q, 5> (&input[33], &input[34]);
	CAS<Q, 5> (&input[35], &input[36]);
	CAS<Q, 5> (&input[37], &input[38]);
	CAS<Q, 5> (&input[39], &input[40]);
	CAS<Q, 5> (&input[41], &input[42]);
	CAS<Q, 5> (&input[43], &input[44]);
	CAS<Q, 5> (&input[45], &input[46]);
	CAS<Q, 5> (&input[47], &input[48]);
	//20
	CAS<Q, 5> (&input[18], &input[19]);
	CAS<Q, 5> (&input[20], &input[21]);
	CAS<Q, 5> (&input[22], &input[23]);
	CAS<Q, 5> (&input[24], &input[25]);
	CAS<Q, 5> (&input[26], &input[27]);
	CAS<Q, 5> (&input[28], &input[29]);
	CAS<Q, 5> (&input[30], &input[31]);
	CAS<Q, 5> (&input[32], &input[33]);
	CAS<Q, 5> (&input[34], &input[35]);
	CAS<Q, 5> (&input[36], &input[37]);
	CAS<Q, 5> (&input[38], &input[39]);
	CAS<Q, 5> (&input[40], &input[41]);
	CAS<Q, 5> (&input[42], &input[43]);
	CAS<Q, 5> (&input[44], &input[45]);
	CAS<Q, 5> (&input[46], &input[47]);
	//21
	CAS<Q, 5> (&input[17], &input[18]);
	CAS<Q, 5> (&input[19], &input[20]);
	CAS<Q, 5> (&input[21], &input[22]);
	CAS<Q, 5> (&input[23], &input[24]);
	CAS<Q, 5> (&input[25], &input[26]);
	CAS<Q, 5> (&input[27], &input[28]);
	CAS<Q, 5> (&input[29], &input[30]);
	CAS<Q, 5> (&input[31], &input[32]);
	CAS<Q, 5> (&input[33], &input[34]);
	CAS<Q, 5> (&input[35], &input[36]);
	CAS<Q, 5> (&input[37], &input[38]);
	CAS<Q, 5> (&input[39], &input[40]);
	CAS<Q, 5> (&input[41], &input[42]);
	CAS<Q, 5> (&input[43], &input[44]);
	CAS<Q, 5> (&input[45], &input[46]);
	CAS<Q, 5> (&input[47], &input[48]);
	//22
	CAS<Q, 5> (&input[18], &input[19]);
	CAS<Q, 5> (&input[20], &input[21]);
	CAS<Q, 5> (&input[22], &input[23]);
	CAS<Q, 5> (&input[24], &input[25]);
	CAS<Q, 5> (&input[26], &input[27]);
	CAS<Q, 5> (&input[28], &input[29]);
	CAS<Q, 5> (&input[30], &input[31]);
	CAS<Q, 5> (&input[32], &input[33]);
	CAS<Q, 5> (&input[34], &input[35]);
	CAS<Q, 5> (&input[36], &input[37]);
	CAS<Q, 5> (&input[38], &input[39]);
	CAS<Q, 5> (&input[40], &input[41]);
	CAS<Q, 5> (&input[42], &input[43]);
	CAS<Q, 5> (&input[44], &input[45]);
	CAS<Q, 5> (&input[46], &input[47]);
	//23
	CAS<Q, 5> (&input[17], &input[18]);
	CAS<Q, 5> (&input[19], &input[20]);
	CAS<Q, 5> (&input[21], &input[22]);
	CAS<Q, 5> (&input[23], &input[24]);
	CAS<Q, 5> (&input[25], &input[26]);
	CAS<Q, 5> (&input[27], &input[28]);
	CAS<Q, 5> (&input[29], &input[30]);
	CAS<Q, 5> (&input[31], &input[32]);
	CAS<Q, 5> (&input[33], &input[34]);
	CAS<Q, 5> (&input[35], &input[36]);
	CAS<Q, 5> (&input[37], &input[38]);
	CAS<Q, 5> (&input[39], &input[40]);
	CAS<Q, 5> (&input[41], &input[42]);
	CAS<Q, 5> (&input[43], &input[44]);
	CAS<Q, 5> (&input[45], &input[46]);
	CAS<Q, 5> (&input[47], &input[48]);
	//24
	CAS<Q, 5> (&input[18], &input[19]);
	CAS<Q, 5> (&input[20], &input[21]);
	CAS<Q, 5> (&input[22], &input[23]);
	CAS<Q, 5> (&input[24], &input[25]);
	CAS<Q, 5> (&input[26], &input[27]);
	CAS<Q, 5> (&input[28], &input[29]);
	CAS<Q, 5> (&input[30], &input[31]);
	CAS<Q, 5> (&input[32], &input[33]);
	CAS<Q, 5> (&input[34], &input[35]);
	CAS<Q, 5> (&input[36], &input[37]);
	CAS<Q, 5> (&input[38], &input[39]);
	CAS<Q, 5> (&input[40], &input[41]);
	CAS<Q, 5> (&input[42], &input[43]);
	CAS<Q, 5> (&input[44], &input[45]);
	CAS<Q, 5> (&input[46], &input[47]);
	//25
	CAS<Q, 5> (&input[17], &input[18]);
	CAS<Q, 5> (&input[19], &input[20]);
	CAS<Q, 5> (&input[21], &input[22]);
	CAS<Q, 5> (&input[23], &input[24]);
	CAS<Q, 5> (&input[25], &input[26]);
	CAS<Q, 5> (&input[27], &input[28]);
	CAS<Q, 5> (&input[29], &input[30]);
	CAS<Q, 5> (&input[31], &input[32]);
	CAS<Q, 5> (&input[33], &input[34]);
	CAS<Q, 5> (&input[35], &input[36]);
	CAS<Q, 5> (&input[37], &input[38]);
	CAS<Q, 5> (&input[39], &input[40]);
	CAS<Q, 5> (&input[41], &input[42]);
	CAS<Q, 5> (&input[43], &input[44]);
	CAS<Q, 5> (&input[45], &input[46]);
	CAS<Q, 5> (&input[47], &input[48]);
	//26
	CAS<Q, 5> (&input[18], &input[19]);
	CAS<Q, 5> (&input[20], &input[21]);
	CAS<Q, 5> (&input[22], &input[23]);
	CAS<Q, 5> (&input[24], &input[25]);
	CAS<Q, 5> (&input[26], &input[27]);
	CAS<Q, 5> (&input[28], &input[29]);
	CAS<Q, 5> (&input[30], &input[31]);
	CAS<Q, 5> (&input[32], &input[33]);
	CAS<Q, 5> (&input[34], &input[35]);
	CAS<Q, 5> (&input[36], &input[37]);
	CAS<Q, 5> (&input[38], &input[39]);
	CAS<Q, 5> (&input[40], &input[41]);
	CAS<Q, 5> (&input[42], &input[43]);
	CAS<Q, 5> (&input[44], &input[45]);
	CAS<Q, 5> (&input[46], &input[47]);
	//27
	CAS<Q, 5> (&input[17], &input[18]);
	CAS<Q, 5> (&input[19], &input[20]);
	CAS<Q, 5> (&input[21], &input[22]);
	CAS<Q, 5> (&input[23], &input[24]);
	CAS<Q, 5> (&input[25], &input[26]);
	CAS<Q, 5> (&input[27], &input[28]);
	CAS<Q, 5> (&input[29], &input[30]);
	CAS<Q, 5> (&input[31], &input[32]);
	CAS<Q, 5> (&input[33], &input[34]);
	CAS<Q, 5> (&input[35], &input[36]);
	CAS<Q, 5> (&input[37], &input[38]);
	CAS<Q, 5> (&input[39], &input[40]);
	CAS<Q, 5> (&input[41], &input[42]);
	CAS<Q, 5> (&input[43], &input[44]);
	CAS<Q, 5> (&input[45], &input[46]);
	CAS<Q, 5> (&input[47], &input[48]);
	//28
	CAS<Q, 5> (&input[18], &input[19]);
	CAS<Q, 5> (&input[20], &input[21]);
	CAS<Q, 5> (&input[22], &input[23]);
	CAS<Q, 5> (&input[24], &input[25]);
	CAS<Q, 5> (&input[26], &input[27]);
	CAS<Q, 5> (&input[28], &input[29]);
	CAS<Q, 5> (&input[30], &input[31]);
	CAS<Q, 5> (&input[32], &input[33]);
	CAS<Q, 5> (&input[34], &input[35]);
	CAS<Q, 5> (&input[36], &input[37]);
	CAS<Q, 5> (&input[38], &input[39]);
	CAS<Q, 5> (&input[40], &input[41]);
	CAS<Q, 5> (&input[42], &input[43]);
	CAS<Q, 5> (&input[44], &input[45]);
	CAS<Q, 5> (&input[46], &input[47]);
	//29
	CAS<Q, 5> (&input[17], &input[18]);
	CAS<Q, 5> (&input[19], &input[20]);
	CAS<Q, 5> (&input[21], &input[22]);
	CAS<Q, 5> (&input[23], &input[24]);
	CAS<Q, 5> (&input[25], &input[26]);
	CAS<Q, 5> (&input[27], &input[28]);
	CAS<Q, 5> (&input[29], &input[30]);
	CAS<Q, 5> (&input[31], &input[32]);
	CAS<Q, 5> (&input[33], &input[34]);
	CAS<Q, 5> (&input[35], &input[36]);
	CAS<Q, 5> (&input[37], &input[38]);
	CAS<Q, 5> (&input[39], &input[40]);
	CAS<Q, 5> (&input[41], &input[42]);
	CAS<Q, 5> (&input[43], &input[44]);
	CAS<Q, 5> (&input[45], &input[46]);
	CAS<Q, 5> (&input[47], &input[48]);
	//30
	CAS<Q, 5> (&input[18], &input[19]);
	CAS<Q, 5> (&input[20], &input[21]);
	CAS<Q, 5> (&input[22], &input[23]);
	CAS<Q, 5> (&input[24], &input[25]);
	CAS<Q, 5> (&input[26], &input[27]);
	CAS<Q, 5> (&input[28], &input[29]);
	CAS<Q, 5> (&input[30], &input[31]);
	CAS<Q, 5> (&input[32], &input[33]);
	CAS<Q, 5> (&input[34], &input[35]);
	CAS<Q, 5> (&input[36], &input[37]);
	CAS<Q, 5> (&input[38], &input[39]);
	CAS<Q, 5> (&input[40], &input[41]);
	CAS<Q, 5> (&input[42], &input[43]);
	CAS<Q, 5> (&input[44], &input[45]);
	CAS<Q, 5> (&input[46], &input[47]);
	//31
	CAS<Q, 5> (&input[17], &input[18]);
	CAS<Q, 5> (&input[19], &input[20]);
	CAS<Q, 5> (&input[21], &input[22]);
	CAS<Q, 5> (&input[23], &input[24]);
	CAS<Q, 5> (&input[25], &input[26]);
	CAS<Q, 5> (&input[27], &input[28]);
	CAS<Q, 5> (&input[29], &input[30]);
	CAS<Q, 5> (&input[31], &input[32]);
	CAS<Q, 5> (&input[33], &input[34]);
	CAS<Q, 5> (&input[35], &input[36]);
	CAS<Q, 5> (&input[37], &input[38]);
	CAS<Q, 5> (&input[39], &input[40]);
	CAS<Q, 5> (&input[41], &input[42]);
	CAS<Q, 5> (&input[43], &input[44]);
	CAS<Q, 5> (&input[45], &input[46]);
	CAS<Q, 5> (&input[47], &input[48]);
	//32
	CAS<Q, 5> (&input[18], &input[19]);
	CAS<Q, 5> (&input[20], &input[21]);
	CAS<Q, 5> (&input[22], &input[23]);
	CAS<Q, 5> (&input[24], &input[25]);
	CAS<Q, 5> (&input[26], &input[27]);
	CAS<Q, 5> (&input[28], &input[29]);
	CAS<Q, 5> (&input[30], &input[31]);
	CAS<Q, 5> (&input[32], &input[33]);
	CAS<Q, 5> (&input[34], &input[35]);
	CAS<Q, 5> (&input[36], &input[37]);
	CAS<Q, 5> (&input[38], &input[39]);
	CAS<Q, 5> (&input[40], &input[41]);
	CAS<Q, 5> (&input[42], &input[43]);
	CAS<Q, 5> (&input[44], &input[45]);
	CAS<Q, 5> (&input[46], &input[47]);

	/////////////////////////////////////////////////

   // output equal the L smallest input
	if (fb == 0){
		output[0] = input[17];
		output[1] = input[18];
		output[2] = input[19];
		output[3] = input[20];
		output[4] = input[21];
		output[5] = input[22];
		output[6] = input[23];
		output[7] = input[24];
		output[8] = input[25];
		output[9] = input[26];
		output[10] = input[27];
		output[11] = input[28];
		output[12] = input[29];
		output[13] = input[30];
		output[14] = input[31];
		output[15] = input[32];
		output[16] = input[33];
		output[17] = input[34];
		output[18] = input[35];
		output[19] = input[36];
		output[20] = input[37];
		output[21] = input[38];
		output[22] = input[39];
		output[23] = input[40];
		output[24] = input[41];
		output[25] = input[42];
		output[26] = input[43];
		output[27] = input[44];
		output[28] = input[45];
		output[29] = input[46];
		output[30] = input[47];
		output[31] = input[48];
	}
	else{
		output[0] = input[0];
		output[1] = input[1];
		output[2] = input[2];
		output[3] = input[3];
		output[4] = input[4];
		output[5] = input[5];
		output[6] = input[6];
		output[7] = input[7];
		output[8] = input[8];
		output[9] = input[9];
		output[10] = input[10];
		output[11] = input[11];
		output[12] = input[12];
		output[13] = input[13];
		output[14] = input[14];
		output[15] = input[15];
		output[16] = input[16];
		output[17] = input[17];
		output[18] = input[18];
		output[19] = input[19];
		output[20] = input[20];
		output[21] = input[21];
		output[22] = input[22];
		output[23] = input[23];
		output[24] = input[24];
		output[25] = input[25];
		output[26] = input[26];
		output[27] = input[27];
		output[28] = input[28];
		output[29] = input[29];
		output[30] = input[30];
		output[31] = input[31];
	}
}

template <int Q>
void CUSTOM_SORT_L64 (PS_struct<1,Q,6> input[128], PS_struct<1,Q,6> output[64], sc_biguint<1> fb)
{
#pragma HLS INLINE
#pragma HLS ARRAY_PARTITION variable=input complete dim=1
#pragma HLS ARRAY_PARTITION variable=output complete dim=1

	// Adjust input when FB=0
	if (fb == 0){
		input[126].metric = MAX_VAL;
		input[125].metric = MAX_VAL;
		input[124].metric = MAX_VAL;
		input[123].metric = MAX_VAL;
		input[122].metric = MAX_VAL;
		input[121].metric = MAX_VAL;
		input[120].metric = MAX_VAL;
		input[119].metric = MAX_VAL;
		input[118].metric = MAX_VAL;
		input[117].metric = MAX_VAL;
		input[116].metric = MAX_VAL;
		input[115].metric = MAX_VAL;
		input[114].metric = MAX_VAL;
		input[113].metric = MAX_VAL;
		input[112].metric = MAX_VAL;
		input[111].metric = MAX_VAL;
		input[110].metric = MAX_VAL;
		input[109].metric = MAX_VAL;
		input[108].metric = MAX_VAL;
		input[107].metric = MAX_VAL;
		input[106].metric = MAX_VAL;
		input[105].metric = MAX_VAL;
		input[104].metric = MAX_VAL;
		input[103].metric = MAX_VAL;
		input[102].metric = MAX_VAL;
		input[101].metric = MAX_VAL;
		input[100].metric = MAX_VAL;
		input[99].metric = MAX_VAL;
		input[98].metric = MAX_VAL;
		input[97].metric = MAX_VAL;
		input[96] = input[63];
		input[95] = input[62];
		input[94] = input[61];
		input[93] = input[60];
		input[92] = input[59];
		input[91] = input[58];
		input[90] = input[57];
		input[89] = input[56];
		input[88] = input[55];
		input[87] = input[54];
		input[86] = input[53];
		input[85] = input[52];
		input[84] = input[51];
		input[83] = input[50];
		input[82] = input[49];
		input[81] = input[48];
		input[80] = input[47];
		input[79] = input[46];
		input[78] = input[45];
		input[77] = input[44];
		input[76] = input[43];
		input[75] = input[42];
		input[74] = input[41];
		input[73] = input[40];
		input[72] = input[39];
		input[71] = input[38];
		input[70] = input[37];
		input[69] = input[36];
		input[68] = input[35];
		input[67] = input[34];
		input[66] = input[33];
		input[65] = input[32];
		input[64] = input[31];
		input[63] = input[30];
		input[62] = input[29];
		input[61] = input[28];
		input[60] = input[27];
		input[59] = input[26];
		input[58] = input[25];
		input[57] = input[24];
		input[56] = input[23];
		input[55] = input[22];
		input[54] = input[21];
		input[53] = input[20];
		input[52] = input[19];
		input[51] = input[18];
		input[50] = input[17];
		input[49] = input[16];
		input[48] = input[15];
		input[47] = input[14];
		input[46] = input[13];
		input[45] = input[12];
		input[44] = input[11];
		input[43] = input[10];
		input[42] = input[9];
		input[41] = input[8];
		input[40] = input[7];
		input[39] = input[6];
		input[38] = input[5];
		input[37] = input[4];
		input[36] = input[3];
		input[35] = input[2];
		input[34] = input[1];
		input[33] = input[0];
		input[32].metric = 0;
		input[31].metric = 0;
		input[30].metric = 0;
		input[29].metric = 0;
		input[28].metric = 0;
		input[27].metric = 0;
		input[26].metric = 0;
		input[25].metric = 0;
		input[24].metric = 0;
		input[23].metric = 0;
		input[22].metric = 0;
		input[21].metric = 0;
		input[20].metric = 0;
		input[19].metric = 0;
		input[18].metric = 0;
		input[17].metric = 0;
		input[16].metric = 0;
		input[15].metric = 0;
		input[14].metric = 0;
		input[13].metric = 0;
		input[12].metric = 0;
		input[11].metric = 0;
		input[10].metric = 0;
		input[9].metric = 0;
		input[8].metric = 0;
		input[7].metric = 0;
		input[6].metric = 0;
		input[5].metric = 0;
		input[4].metric = 0;
		input[3].metric = 0;
		input[2].metric = 0;
		input[1].metric = 0;
	}

	//1
	CAS<Q, 6> (&input[1], &input[2]);
	CAS<Q, 6> (&input[3], &input[4]);
	CAS<Q, 6> (&input[5], &input[6]);
	CAS<Q, 6> (&input[7], &input[8]);
	CAS<Q, 6> (&input[9], &input[10]);
	CAS<Q, 6> (&input[11], &input[12]);
	CAS<Q, 6> (&input[13], &input[14]);
	CAS<Q, 6> (&input[15], &input[16]);
	CAS<Q, 6> (&input[17], &input[18]);
	CAS<Q, 6> (&input[19], &input[20]);
	CAS<Q, 6> (&input[21], &input[22]);
	CAS<Q, 6> (&input[23], &input[24]);
	CAS<Q, 6> (&input[25], &input[26]);
	CAS<Q, 6> (&input[27], &input[28]);
	CAS<Q, 6> (&input[29], &input[30]);
	CAS<Q, 6> (&input[31], &input[32]);
	CAS<Q, 6> (&input[33], &input[34]);
	CAS<Q, 6> (&input[35], &input[36]);
	CAS<Q, 6> (&input[37], &input[38]);
	CAS<Q, 6> (&input[39], &input[40]);
	CAS<Q, 6> (&input[41], &input[42]);
	CAS<Q, 6> (&input[43], &input[44]);
	CAS<Q, 6> (&input[45], &input[46]);
	CAS<Q, 6> (&input[47], &input[48]);
	CAS<Q, 6> (&input[49], &input[50]);
	CAS<Q, 6> (&input[51], &input[52]);
	CAS<Q, 6> (&input[53], &input[54]);
	CAS<Q, 6> (&input[55], &input[56]);
	CAS<Q, 6> (&input[57], &input[58]);
	CAS<Q, 6> (&input[59], &input[60]);
	CAS<Q, 6> (&input[61], &input[62]);
	CAS<Q, 6> (&input[63], &input[64]);
	CAS<Q, 6> (&input[65], &input[66]);
	CAS<Q, 6> (&input[67], &input[68]);
	CAS<Q, 6> (&input[69], &input[70]);
	CAS<Q, 6> (&input[71], &input[72]);
	CAS<Q, 6> (&input[73], &input[74]);
	CAS<Q, 6> (&input[75], &input[76]);
	CAS<Q, 6> (&input[77], &input[78]);
	CAS<Q, 6> (&input[79], &input[80]);
	CAS<Q, 6> (&input[81], &input[82]);
	CAS<Q, 6> (&input[83], &input[84]);
	CAS<Q, 6> (&input[85], &input[86]);
	CAS<Q, 6> (&input[87], &input[88]);
	CAS<Q, 6> (&input[89], &input[90]);
	CAS<Q, 6> (&input[81], &input[92]);
	CAS<Q, 6> (&input[93], &input[94]);
	CAS<Q, 6> (&input[95], &input[96]);
	CAS<Q, 6> (&input[97], &input[98]);
	CAS<Q, 6> (&input[99], &input[100]);
	CAS<Q, 6> (&input[101], &input[102]);
	CAS<Q, 6> (&input[103], &input[104]);
	CAS<Q, 6> (&input[105], &input[106]);
	CAS<Q, 6> (&input[107], &input[108]);
	CAS<Q, 6> (&input[109], &input[110]);
	CAS<Q, 6> (&input[111], &input[112]);
	CAS<Q, 6> (&input[113], &input[114]);
	CAS<Q, 6> (&input[115], &input[116]);
	CAS<Q, 6> (&input[117], &input[118]);
	CAS<Q, 6> (&input[119], &input[120]);
	CAS<Q, 6> (&input[121], &input[122]);
	CAS<Q, 6> (&input[123], &input[124]);
	CAS<Q, 6> (&input[125], &input[126]);
	//2
	CAS<Q, 6> (&input[2], &input[3]);
	CAS<Q, 6> (&input[4], &input[5]);
	CAS<Q, 6> (&input[6], &input[7]);
	CAS<Q, 6> (&input[8], &input[9]);
	CAS<Q, 6> (&input[10], &input[11]);
	CAS<Q, 6> (&input[12], &input[13]);
	CAS<Q, 6> (&input[14], &input[15]);
	CAS<Q, 6> (&input[16], &input[17]);
	CAS<Q, 6> (&input[18], &input[19]);
	CAS<Q, 6> (&input[20], &input[21]);
	CAS<Q, 6> (&input[22], &input[23]);
	CAS<Q, 6> (&input[24], &input[25]);
	CAS<Q, 6> (&input[26], &input[27]);
	CAS<Q, 6> (&input[28], &input[29]);
	CAS<Q, 6> (&input[30], &input[31]);
	CAS<Q, 6> (&input[32], &input[33]);
	CAS<Q, 6> (&input[34], &input[35]);
	CAS<Q, 6> (&input[36], &input[37]);
	CAS<Q, 6> (&input[38], &input[39]);
	CAS<Q, 6> (&input[40], &input[41]);
	CAS<Q, 6> (&input[42], &input[43]);
	CAS<Q, 6> (&input[44], &input[45]);
	CAS<Q, 6> (&input[46], &input[47]);
	CAS<Q, 6> (&input[48], &input[49]);
	CAS<Q, 6> (&input[50], &input[51]);
	CAS<Q, 6> (&input[52], &input[53]);
	CAS<Q, 6> (&input[54], &input[55]);
	CAS<Q, 6> (&input[56], &input[57]);
	CAS<Q, 6> (&input[58], &input[59]);
	CAS<Q, 6> (&input[60], &input[61]);
	CAS<Q, 6> (&input[62], &input[63]);
	CAS<Q, 6> (&input[64], &input[65]);
	CAS<Q, 6> (&input[66], &input[67]);
	CAS<Q, 6> (&input[68], &input[69]);
	CAS<Q, 6> (&input[70], &input[71]);
	CAS<Q, 6> (&input[72], &input[73]);
	CAS<Q, 6> (&input[74], &input[75]);
	CAS<Q, 6> (&input[76], &input[77]);
	CAS<Q, 6> (&input[78], &input[79]);
	CAS<Q, 6> (&input[80], &input[81]);
	CAS<Q, 6> (&input[82], &input[83]);
	CAS<Q, 6> (&input[84], &input[85]);
	CAS<Q, 6> (&input[86], &input[87]);
	CAS<Q, 6> (&input[88], &input[89]);
	CAS<Q, 6> (&input[90], &input[91]);
	CAS<Q, 6> (&input[92], &input[93]);
	CAS<Q, 6> (&input[94], &input[95]);
	CAS<Q, 6> (&input[96], &input[97]);
	CAS<Q, 6> (&input[98], &input[99]);
	CAS<Q, 6> (&input[100], &input[101]);
	CAS<Q, 6> (&input[102], &input[103]);
	CAS<Q, 6> (&input[104], &input[105]);
	CAS<Q, 6> (&input[106], &input[107]);
	CAS<Q, 6> (&input[108], &input[109]);
	CAS<Q, 6> (&input[110], &input[111]);
	CAS<Q, 6> (&input[112], &input[113]);
	CAS<Q, 6> (&input[114], &input[115]);
	CAS<Q, 6> (&input[116], &input[117]);
	CAS<Q, 6> (&input[118], &input[119]);
	CAS<Q, 6> (&input[120], &input[121]);
	CAS<Q, 6> (&input[122], &input[123]);
	CAS<Q, 6> (&input[124], &input[125]);
	//3
	CAS<Q, 6> (&input[3], &input[4]);
	CAS<Q, 6> (&input[5], &input[6]);
	CAS<Q, 6> (&input[7], &input[8]);
	CAS<Q, 6> (&input[9], &input[10]);
	CAS<Q, 6> (&input[11], &input[12]);
	CAS<Q, 6> (&input[13], &input[14]);
	CAS<Q, 6> (&input[15], &input[16]);
	CAS<Q, 6> (&input[17], &input[18]);
	CAS<Q, 6> (&input[19], &input[20]);
	CAS<Q, 6> (&input[21], &input[22]);
	CAS<Q, 6> (&input[23], &input[24]);
	CAS<Q, 6> (&input[25], &input[26]);
	CAS<Q, 6> (&input[27], &input[28]);
	CAS<Q, 6> (&input[29], &input[30]);
	CAS<Q, 6> (&input[31], &input[32]);
	CAS<Q, 6> (&input[33], &input[34]);
	CAS<Q, 6> (&input[35], &input[36]);
	CAS<Q, 6> (&input[37], &input[38]);
	CAS<Q, 6> (&input[39], &input[40]);
	CAS<Q, 6> (&input[41], &input[42]);
	CAS<Q, 6> (&input[43], &input[44]);
	CAS<Q, 6> (&input[45], &input[46]);
	CAS<Q, 6> (&input[47], &input[48]);
	CAS<Q, 6> (&input[49], &input[50]);
	CAS<Q, 6> (&input[51], &input[52]);
	CAS<Q, 6> (&input[53], &input[54]);
	CAS<Q, 6> (&input[55], &input[56]);
	CAS<Q, 6> (&input[57], &input[58]);
	CAS<Q, 6> (&input[59], &input[60]);
	CAS<Q, 6> (&input[61], &input[62]);
	CAS<Q, 6> (&input[63], &input[64]);
	CAS<Q, 6> (&input[65], &input[66]);
	CAS<Q, 6> (&input[67], &input[68]);
	CAS<Q, 6> (&input[69], &input[70]);
	CAS<Q, 6> (&input[71], &input[72]);
	CAS<Q, 6> (&input[73], &input[74]);
	CAS<Q, 6> (&input[75], &input[76]);
	CAS<Q, 6> (&input[77], &input[78]);
	CAS<Q, 6> (&input[79], &input[80]);
	CAS<Q, 6> (&input[81], &input[82]);
	CAS<Q, 6> (&input[83], &input[84]);
	CAS<Q, 6> (&input[85], &input[86]);
	CAS<Q, 6> (&input[87], &input[88]);
	CAS<Q, 6> (&input[89], &input[90]);
	CAS<Q, 6> (&input[81], &input[92]);
	CAS<Q, 6> (&input[93], &input[94]);
	CAS<Q, 6> (&input[95], &input[96]);
	CAS<Q, 6> (&input[97], &input[98]);
	CAS<Q, 6> (&input[99], &input[100]);
	CAS<Q, 6> (&input[101], &input[102]);
	CAS<Q, 6> (&input[103], &input[104]);
	CAS<Q, 6> (&input[105], &input[106]);
	CAS<Q, 6> (&input[107], &input[108]);
	CAS<Q, 6> (&input[109], &input[110]);
	CAS<Q, 6> (&input[111], &input[112]);
	CAS<Q, 6> (&input[113], &input[114]);
	CAS<Q, 6> (&input[115], &input[116]);
	CAS<Q, 6> (&input[117], &input[118]);
	CAS<Q, 6> (&input[119], &input[120]);
	CAS<Q, 6> (&input[121], &input[122]);
	CAS<Q, 6> (&input[123], &input[124]);
	//4
	CAS<Q, 6> (&input[4], &input[5]);
	CAS<Q, 6> (&input[6], &input[7]);
	CAS<Q, 6> (&input[8], &input[9]);
	CAS<Q, 6> (&input[10], &input[11]);
	CAS<Q, 6> (&input[12], &input[13]);
	CAS<Q, 6> (&input[14], &input[15]);
	CAS<Q, 6> (&input[16], &input[17]);
	CAS<Q, 6> (&input[18], &input[19]);
	CAS<Q, 6> (&input[20], &input[21]);
	CAS<Q, 6> (&input[22], &input[23]);
	CAS<Q, 6> (&input[24], &input[25]);
	CAS<Q, 6> (&input[26], &input[27]);
	CAS<Q, 6> (&input[28], &input[29]);
	CAS<Q, 6> (&input[30], &input[31]);
	CAS<Q, 6> (&input[32], &input[33]);
	CAS<Q, 6> (&input[34], &input[35]);
	CAS<Q, 6> (&input[36], &input[37]);
	CAS<Q, 6> (&input[38], &input[39]);
	CAS<Q, 6> (&input[40], &input[41]);
	CAS<Q, 6> (&input[42], &input[43]);
	CAS<Q, 6> (&input[44], &input[45]);
	CAS<Q, 6> (&input[46], &input[47]);
	CAS<Q, 6> (&input[48], &input[49]);
	CAS<Q, 6> (&input[50], &input[51]);
	CAS<Q, 6> (&input[52], &input[53]);
	CAS<Q, 6> (&input[54], &input[55]);
	CAS<Q, 6> (&input[56], &input[57]);
	CAS<Q, 6> (&input[58], &input[59]);
	CAS<Q, 6> (&input[60], &input[61]);
	CAS<Q, 6> (&input[62], &input[63]);
	CAS<Q, 6> (&input[64], &input[65]);
	CAS<Q, 6> (&input[66], &input[67]);
	CAS<Q, 6> (&input[68], &input[69]);
	CAS<Q, 6> (&input[70], &input[71]);
	CAS<Q, 6> (&input[72], &input[73]);
	CAS<Q, 6> (&input[74], &input[75]);
	CAS<Q, 6> (&input[76], &input[77]);
	CAS<Q, 6> (&input[78], &input[79]);
	CAS<Q, 6> (&input[80], &input[81]);
	CAS<Q, 6> (&input[82], &input[83]);
	CAS<Q, 6> (&input[84], &input[85]);
	CAS<Q, 6> (&input[86], &input[87]);
	CAS<Q, 6> (&input[88], &input[89]);
	CAS<Q, 6> (&input[90], &input[91]);
	CAS<Q, 6> (&input[92], &input[93]);
	CAS<Q, 6> (&input[94], &input[95]);
	CAS<Q, 6> (&input[96], &input[97]);
	CAS<Q, 6> (&input[98], &input[99]);
	CAS<Q, 6> (&input[100], &input[101]);
	CAS<Q, 6> (&input[102], &input[103]);
	CAS<Q, 6> (&input[104], &input[105]);
	CAS<Q, 6> (&input[106], &input[107]);
	CAS<Q, 6> (&input[108], &input[109]);
	CAS<Q, 6> (&input[110], &input[111]);
	CAS<Q, 6> (&input[112], &input[113]);
	CAS<Q, 6> (&input[114], &input[115]);
	CAS<Q, 6> (&input[116], &input[117]);
	CAS<Q, 6> (&input[118], &input[119]);
	CAS<Q, 6> (&input[120], &input[121]);
	CAS<Q, 6> (&input[122], &input[123]);
	//5
	CAS<Q, 6> (&input[5], &input[6]);
	CAS<Q, 6> (&input[7], &input[8]);
	CAS<Q, 6> (&input[9], &input[10]);
	CAS<Q, 6> (&input[11], &input[12]);
	CAS<Q, 6> (&input[13], &input[14]);
	CAS<Q, 6> (&input[15], &input[16]);
	CAS<Q, 6> (&input[17], &input[18]);
	CAS<Q, 6> (&input[19], &input[20]);
	CAS<Q, 6> (&input[21], &input[22]);
	CAS<Q, 6> (&input[23], &input[24]);
	CAS<Q, 6> (&input[25], &input[26]);
	CAS<Q, 6> (&input[27], &input[28]);
	CAS<Q, 6> (&input[29], &input[30]);
	CAS<Q, 6> (&input[31], &input[32]);
	CAS<Q, 6> (&input[33], &input[34]);
	CAS<Q, 6> (&input[35], &input[36]);
	CAS<Q, 6> (&input[37], &input[38]);
	CAS<Q, 6> (&input[39], &input[40]);
	CAS<Q, 6> (&input[41], &input[42]);
	CAS<Q, 6> (&input[43], &input[44]);
	CAS<Q, 6> (&input[45], &input[46]);
	CAS<Q, 6> (&input[47], &input[48]);
	CAS<Q, 6> (&input[49], &input[50]);
	CAS<Q, 6> (&input[51], &input[52]);
	CAS<Q, 6> (&input[53], &input[54]);
	CAS<Q, 6> (&input[55], &input[56]);
	CAS<Q, 6> (&input[57], &input[58]);
	CAS<Q, 6> (&input[59], &input[60]);
	CAS<Q, 6> (&input[61], &input[62]);
	CAS<Q, 6> (&input[63], &input[64]);
	CAS<Q, 6> (&input[65], &input[66]);
	CAS<Q, 6> (&input[67], &input[68]);
	CAS<Q, 6> (&input[69], &input[70]);
	CAS<Q, 6> (&input[71], &input[72]);
	CAS<Q, 6> (&input[73], &input[74]);
	CAS<Q, 6> (&input[75], &input[76]);
	CAS<Q, 6> (&input[77], &input[78]);
	CAS<Q, 6> (&input[79], &input[80]);
	CAS<Q, 6> (&input[81], &input[82]);
	CAS<Q, 6> (&input[83], &input[84]);
	CAS<Q, 6> (&input[85], &input[86]);
	CAS<Q, 6> (&input[87], &input[88]);
	CAS<Q, 6> (&input[89], &input[90]);
	CAS<Q, 6> (&input[81], &input[92]);
	CAS<Q, 6> (&input[93], &input[94]);
	CAS<Q, 6> (&input[95], &input[96]);
	CAS<Q, 6> (&input[97], &input[98]);
	CAS<Q, 6> (&input[99], &input[100]);
	CAS<Q, 6> (&input[101], &input[102]);
	CAS<Q, 6> (&input[103], &input[104]);
	CAS<Q, 6> (&input[105], &input[106]);
	CAS<Q, 6> (&input[107], &input[108]);
	CAS<Q, 6> (&input[109], &input[110]);
	CAS<Q, 6> (&input[111], &input[112]);
	CAS<Q, 6> (&input[113], &input[114]);
	CAS<Q, 6> (&input[115], &input[116]);
	CAS<Q, 6> (&input[117], &input[118]);
	CAS<Q, 6> (&input[119], &input[120]);
	CAS<Q, 6> (&input[121], &input[122]);
	//6
	CAS<Q, 6> (&input[6], &input[7]);
	CAS<Q, 6> (&input[8], &input[9]);
	CAS<Q, 6> (&input[10], &input[11]);
	CAS<Q, 6> (&input[12], &input[13]);
	CAS<Q, 6> (&input[14], &input[15]);
	CAS<Q, 6> (&input[16], &input[17]);
	CAS<Q, 6> (&input[18], &input[19]);
	CAS<Q, 6> (&input[20], &input[21]);
	CAS<Q, 6> (&input[22], &input[23]);
	CAS<Q, 6> (&input[24], &input[25]);
	CAS<Q, 6> (&input[26], &input[27]);
	CAS<Q, 6> (&input[28], &input[29]);
	CAS<Q, 6> (&input[30], &input[31]);
	CAS<Q, 6> (&input[32], &input[33]);
	CAS<Q, 6> (&input[34], &input[35]);
	CAS<Q, 6> (&input[36], &input[37]);
	CAS<Q, 6> (&input[38], &input[39]);
	CAS<Q, 6> (&input[40], &input[41]);
	CAS<Q, 6> (&input[42], &input[43]);
	CAS<Q, 6> (&input[44], &input[45]);
	CAS<Q, 6> (&input[46], &input[47]);
	CAS<Q, 6> (&input[48], &input[49]);
	CAS<Q, 6> (&input[50], &input[51]);
	CAS<Q, 6> (&input[52], &input[53]);
	CAS<Q, 6> (&input[54], &input[55]);
	CAS<Q, 6> (&input[56], &input[57]);
	CAS<Q, 6> (&input[58], &input[59]);
	CAS<Q, 6> (&input[60], &input[61]);
	CAS<Q, 6> (&input[62], &input[63]);
	CAS<Q, 6> (&input[64], &input[65]);
	CAS<Q, 6> (&input[66], &input[67]);
	CAS<Q, 6> (&input[68], &input[69]);
	CAS<Q, 6> (&input[70], &input[71]);
	CAS<Q, 6> (&input[72], &input[73]);
	CAS<Q, 6> (&input[74], &input[75]);
	CAS<Q, 6> (&input[76], &input[77]);
	CAS<Q, 6> (&input[78], &input[79]);
	CAS<Q, 6> (&input[80], &input[81]);
	CAS<Q, 6> (&input[82], &input[83]);
	CAS<Q, 6> (&input[84], &input[85]);
	CAS<Q, 6> (&input[86], &input[87]);
	CAS<Q, 6> (&input[88], &input[89]);
	CAS<Q, 6> (&input[90], &input[91]);
	CAS<Q, 6> (&input[92], &input[93]);
	CAS<Q, 6> (&input[94], &input[95]);
	CAS<Q, 6> (&input[96], &input[97]);
	CAS<Q, 6> (&input[98], &input[99]);
	CAS<Q, 6> (&input[100], &input[101]);
	CAS<Q, 6> (&input[102], &input[103]);
	CAS<Q, 6> (&input[104], &input[105]);
	CAS<Q, 6> (&input[106], &input[107]);
	CAS<Q, 6> (&input[108], &input[109]);
	CAS<Q, 6> (&input[110], &input[111]);
	CAS<Q, 6> (&input[112], &input[113]);
	CAS<Q, 6> (&input[114], &input[115]);
	CAS<Q, 6> (&input[116], &input[117]);
	CAS<Q, 6> (&input[118], &input[119]);
	CAS<Q, 6> (&input[120], &input[121]);
	//7
	CAS<Q, 6> (&input[7], &input[8]);
	CAS<Q, 6> (&input[9], &input[10]);
	CAS<Q, 6> (&input[11], &input[12]);
	CAS<Q, 6> (&input[13], &input[14]);
	CAS<Q, 6> (&input[15], &input[16]);
	CAS<Q, 6> (&input[17], &input[18]);
	CAS<Q, 6> (&input[19], &input[20]);
	CAS<Q, 6> (&input[21], &input[22]);
	CAS<Q, 6> (&input[23], &input[24]);
	CAS<Q, 6> (&input[25], &input[26]);
	CAS<Q, 6> (&input[27], &input[28]);
	CAS<Q, 6> (&input[29], &input[30]);
	CAS<Q, 6> (&input[31], &input[32]);
	CAS<Q, 6> (&input[33], &input[34]);
	CAS<Q, 6> (&input[35], &input[36]);
	CAS<Q, 6> (&input[37], &input[38]);
	CAS<Q, 6> (&input[39], &input[40]);
	CAS<Q, 6> (&input[41], &input[42]);
	CAS<Q, 6> (&input[43], &input[44]);
	CAS<Q, 6> (&input[45], &input[46]);
	CAS<Q, 6> (&input[47], &input[48]);
	CAS<Q, 6> (&input[49], &input[50]);
	CAS<Q, 6> (&input[51], &input[52]);
	CAS<Q, 6> (&input[53], &input[54]);
	CAS<Q, 6> (&input[55], &input[56]);
	CAS<Q, 6> (&input[57], &input[58]);
	CAS<Q, 6> (&input[59], &input[60]);
	CAS<Q, 6> (&input[61], &input[62]);
	CAS<Q, 6> (&input[63], &input[64]);
	CAS<Q, 6> (&input[65], &input[66]);
	CAS<Q, 6> (&input[67], &input[68]);
	CAS<Q, 6> (&input[69], &input[70]);
	CAS<Q, 6> (&input[71], &input[72]);
	CAS<Q, 6> (&input[73], &input[74]);
	CAS<Q, 6> (&input[75], &input[76]);
	CAS<Q, 6> (&input[77], &input[78]);
	CAS<Q, 6> (&input[79], &input[80]);
	CAS<Q, 6> (&input[81], &input[82]);
	CAS<Q, 6> (&input[83], &input[84]);
	CAS<Q, 6> (&input[85], &input[86]);
	CAS<Q, 6> (&input[87], &input[88]);
	CAS<Q, 6> (&input[89], &input[90]);
	CAS<Q, 6> (&input[81], &input[92]);
	CAS<Q, 6> (&input[93], &input[94]);
	CAS<Q, 6> (&input[95], &input[96]);
	CAS<Q, 6> (&input[97], &input[98]);
	CAS<Q, 6> (&input[99], &input[100]);
	CAS<Q, 6> (&input[101], &input[102]);
	CAS<Q, 6> (&input[103], &input[104]);
	CAS<Q, 6> (&input[105], &input[106]);
	CAS<Q, 6> (&input[107], &input[108]);
	CAS<Q, 6> (&input[109], &input[110]);
	CAS<Q, 6> (&input[111], &input[112]);
	CAS<Q, 6> (&input[113], &input[114]);
	CAS<Q, 6> (&input[115], &input[116]);
	CAS<Q, 6> (&input[117], &input[118]);
	CAS<Q, 6> (&input[119], &input[120]);
	//8
	CAS<Q, 6> (&input[8], &input[9]);
	CAS<Q, 6> (&input[10], &input[11]);
	CAS<Q, 6> (&input[12], &input[13]);
	CAS<Q, 6> (&input[14], &input[15]);
	CAS<Q, 6> (&input[16], &input[17]);
	CAS<Q, 6> (&input[18], &input[19]);
	CAS<Q, 6> (&input[20], &input[21]);
	CAS<Q, 6> (&input[22], &input[23]);
	CAS<Q, 6> (&input[24], &input[25]);
	CAS<Q, 6> (&input[26], &input[27]);
	CAS<Q, 6> (&input[28], &input[29]);
	CAS<Q, 6> (&input[30], &input[31]);
	CAS<Q, 6> (&input[32], &input[33]);
	CAS<Q, 6> (&input[34], &input[35]);
	CAS<Q, 6> (&input[36], &input[37]);
	CAS<Q, 6> (&input[38], &input[39]);
	CAS<Q, 6> (&input[40], &input[41]);
	CAS<Q, 6> (&input[42], &input[43]);
	CAS<Q, 6> (&input[44], &input[45]);
	CAS<Q, 6> (&input[46], &input[47]);
	CAS<Q, 6> (&input[48], &input[49]);
	CAS<Q, 6> (&input[50], &input[51]);
	CAS<Q, 6> (&input[52], &input[53]);
	CAS<Q, 6> (&input[54], &input[55]);
	CAS<Q, 6> (&input[56], &input[57]);
	CAS<Q, 6> (&input[58], &input[59]);
	CAS<Q, 6> (&input[60], &input[61]);
	CAS<Q, 6> (&input[62], &input[63]);
	CAS<Q, 6> (&input[64], &input[65]);
	CAS<Q, 6> (&input[66], &input[67]);
	CAS<Q, 6> (&input[68], &input[69]);
	CAS<Q, 6> (&input[70], &input[71]);
	CAS<Q, 6> (&input[72], &input[73]);
	CAS<Q, 6> (&input[74], &input[75]);
	CAS<Q, 6> (&input[76], &input[77]);
	CAS<Q, 6> (&input[78], &input[79]);
	CAS<Q, 6> (&input[80], &input[81]);
	CAS<Q, 6> (&input[82], &input[83]);
	CAS<Q, 6> (&input[84], &input[85]);
	CAS<Q, 6> (&input[86], &input[87]);
	CAS<Q, 6> (&input[88], &input[89]);
	CAS<Q, 6> (&input[90], &input[91]);
	CAS<Q, 6> (&input[92], &input[93]);
	CAS<Q, 6> (&input[94], &input[95]);
	CAS<Q, 6> (&input[96], &input[97]);
	CAS<Q, 6> (&input[98], &input[99]);
	CAS<Q, 6> (&input[100], &input[101]);
	CAS<Q, 6> (&input[102], &input[103]);
	CAS<Q, 6> (&input[104], &input[105]);
	CAS<Q, 6> (&input[106], &input[107]);
	CAS<Q, 6> (&input[108], &input[109]);
	CAS<Q, 6> (&input[110], &input[111]);
	CAS<Q, 6> (&input[112], &input[113]);
	CAS<Q, 6> (&input[114], &input[115]);
	CAS<Q, 6> (&input[116], &input[117]);
	CAS<Q, 6> (&input[118], &input[119]);
	//9
	CAS<Q, 6> (&input[9], &input[10]);
	CAS<Q, 6> (&input[11], &input[12]);
	CAS<Q, 6> (&input[13], &input[14]);
	CAS<Q, 6> (&input[15], &input[16]);
	CAS<Q, 6> (&input[17], &input[18]);
	CAS<Q, 6> (&input[19], &input[20]);
	CAS<Q, 6> (&input[21], &input[22]);
	CAS<Q, 6> (&input[23], &input[24]);
	CAS<Q, 6> (&input[25], &input[26]);
	CAS<Q, 6> (&input[27], &input[28]);
	CAS<Q, 6> (&input[29], &input[30]);
	CAS<Q, 6> (&input[31], &input[32]);
	CAS<Q, 6> (&input[33], &input[34]);
	CAS<Q, 6> (&input[35], &input[36]);
	CAS<Q, 6> (&input[37], &input[38]);
	CAS<Q, 6> (&input[39], &input[40]);
	CAS<Q, 6> (&input[41], &input[42]);
	CAS<Q, 6> (&input[43], &input[44]);
	CAS<Q, 6> (&input[45], &input[46]);
	CAS<Q, 6> (&input[47], &input[48]);
	CAS<Q, 6> (&input[49], &input[50]);
	CAS<Q, 6> (&input[51], &input[52]);
	CAS<Q, 6> (&input[53], &input[54]);
	CAS<Q, 6> (&input[55], &input[56]);
	CAS<Q, 6> (&input[57], &input[58]);
	CAS<Q, 6> (&input[59], &input[60]);
	CAS<Q, 6> (&input[61], &input[62]);
	CAS<Q, 6> (&input[63], &input[64]);
	CAS<Q, 6> (&input[65], &input[66]);
	CAS<Q, 6> (&input[67], &input[68]);
	CAS<Q, 6> (&input[69], &input[70]);
	CAS<Q, 6> (&input[71], &input[72]);
	CAS<Q, 6> (&input[73], &input[74]);
	CAS<Q, 6> (&input[75], &input[76]);
	CAS<Q, 6> (&input[77], &input[78]);
	CAS<Q, 6> (&input[79], &input[80]);
	CAS<Q, 6> (&input[81], &input[82]);
	CAS<Q, 6> (&input[83], &input[84]);
	CAS<Q, 6> (&input[85], &input[86]);
	CAS<Q, 6> (&input[87], &input[88]);
	CAS<Q, 6> (&input[89], &input[90]);
	CAS<Q, 6> (&input[81], &input[92]);
	CAS<Q, 6> (&input[93], &input[94]);
	CAS<Q, 6> (&input[95], &input[96]);
	CAS<Q, 6> (&input[97], &input[98]);
	CAS<Q, 6> (&input[99], &input[100]);
	CAS<Q, 6> (&input[101], &input[102]);
	CAS<Q, 6> (&input[103], &input[104]);
	CAS<Q, 6> (&input[105], &input[106]);
	CAS<Q, 6> (&input[107], &input[108]);
	CAS<Q, 6> (&input[109], &input[110]);
	CAS<Q, 6> (&input[111], &input[112]);
	CAS<Q, 6> (&input[113], &input[114]);
	CAS<Q, 6> (&input[115], &input[116]);
	CAS<Q, 6> (&input[117], &input[118]);
	//10
	CAS<Q, 6> (&input[10], &input[11]);
	CAS<Q, 6> (&input[12], &input[13]);
	CAS<Q, 6> (&input[14], &input[15]);
	CAS<Q, 6> (&input[16], &input[17]);
	CAS<Q, 6> (&input[18], &input[19]);
	CAS<Q, 6> (&input[20], &input[21]);
	CAS<Q, 6> (&input[22], &input[23]);
	CAS<Q, 6> (&input[24], &input[25]);
	CAS<Q, 6> (&input[26], &input[27]);
	CAS<Q, 6> (&input[28], &input[29]);
	CAS<Q, 6> (&input[30], &input[31]);
	CAS<Q, 6> (&input[32], &input[33]);
	CAS<Q, 6> (&input[34], &input[35]);
	CAS<Q, 6> (&input[36], &input[37]);
	CAS<Q, 6> (&input[38], &input[39]);
	CAS<Q, 6> (&input[40], &input[41]);
	CAS<Q, 6> (&input[42], &input[43]);
	CAS<Q, 6> (&input[44], &input[45]);
	CAS<Q, 6> (&input[46], &input[47]);
	CAS<Q, 6> (&input[48], &input[49]);
	CAS<Q, 6> (&input[50], &input[51]);
	CAS<Q, 6> (&input[52], &input[53]);
	CAS<Q, 6> (&input[54], &input[55]);
	CAS<Q, 6> (&input[56], &input[57]);
	CAS<Q, 6> (&input[58], &input[59]);
	CAS<Q, 6> (&input[60], &input[61]);
	CAS<Q, 6> (&input[62], &input[63]);
	CAS<Q, 6> (&input[64], &input[65]);
	CAS<Q, 6> (&input[66], &input[67]);
	CAS<Q, 6> (&input[68], &input[69]);
	CAS<Q, 6> (&input[70], &input[71]);
	CAS<Q, 6> (&input[72], &input[73]);
	CAS<Q, 6> (&input[74], &input[75]);
	CAS<Q, 6> (&input[76], &input[77]);
	CAS<Q, 6> (&input[78], &input[79]);
	CAS<Q, 6> (&input[80], &input[81]);
	CAS<Q, 6> (&input[82], &input[83]);
	CAS<Q, 6> (&input[84], &input[85]);
	CAS<Q, 6> (&input[86], &input[87]);
	CAS<Q, 6> (&input[88], &input[89]);
	CAS<Q, 6> (&input[90], &input[91]);
	CAS<Q, 6> (&input[92], &input[93]);
	CAS<Q, 6> (&input[94], &input[95]);
	CAS<Q, 6> (&input[96], &input[97]);
	CAS<Q, 6> (&input[98], &input[99]);
	CAS<Q, 6> (&input[100], &input[101]);
	CAS<Q, 6> (&input[102], &input[103]);
	CAS<Q, 6> (&input[104], &input[105]);
	CAS<Q, 6> (&input[106], &input[107]);
	CAS<Q, 6> (&input[108], &input[109]);
	CAS<Q, 6> (&input[110], &input[111]);
	CAS<Q, 6> (&input[112], &input[113]);
	CAS<Q, 6> (&input[114], &input[115]);
	CAS<Q, 6> (&input[116], &input[117]);
	//11
	CAS<Q, 6> (&input[11], &input[12]);
	CAS<Q, 6> (&input[13], &input[14]);
	CAS<Q, 6> (&input[15], &input[16]);
	CAS<Q, 6> (&input[17], &input[18]);
	CAS<Q, 6> (&input[19], &input[20]);
	CAS<Q, 6> (&input[21], &input[22]);
	CAS<Q, 6> (&input[23], &input[24]);
	CAS<Q, 6> (&input[25], &input[26]);
	CAS<Q, 6> (&input[27], &input[28]);
	CAS<Q, 6> (&input[29], &input[30]);
	CAS<Q, 6> (&input[31], &input[32]);
	CAS<Q, 6> (&input[33], &input[34]);
	CAS<Q, 6> (&input[35], &input[36]);
	CAS<Q, 6> (&input[37], &input[38]);
	CAS<Q, 6> (&input[39], &input[40]);
	CAS<Q, 6> (&input[41], &input[42]);
	CAS<Q, 6> (&input[43], &input[44]);
	CAS<Q, 6> (&input[45], &input[46]);
	CAS<Q, 6> (&input[47], &input[48]);
	CAS<Q, 6> (&input[49], &input[50]);
	CAS<Q, 6> (&input[51], &input[52]);
	CAS<Q, 6> (&input[53], &input[54]);
	CAS<Q, 6> (&input[55], &input[56]);
	CAS<Q, 6> (&input[57], &input[58]);
	CAS<Q, 6> (&input[59], &input[60]);
	CAS<Q, 6> (&input[61], &input[62]);
	CAS<Q, 6> (&input[63], &input[64]);
	CAS<Q, 6> (&input[65], &input[66]);
	CAS<Q, 6> (&input[67], &input[68]);
	CAS<Q, 6> (&input[69], &input[70]);
	CAS<Q, 6> (&input[71], &input[72]);
	CAS<Q, 6> (&input[73], &input[74]);
	CAS<Q, 6> (&input[75], &input[76]);
	CAS<Q, 6> (&input[77], &input[78]);
	CAS<Q, 6> (&input[79], &input[80]);
	CAS<Q, 6> (&input[81], &input[82]);
	CAS<Q, 6> (&input[83], &input[84]);
	CAS<Q, 6> (&input[85], &input[86]);
	CAS<Q, 6> (&input[87], &input[88]);
	CAS<Q, 6> (&input[89], &input[90]);
	CAS<Q, 6> (&input[81], &input[92]);
	CAS<Q, 6> (&input[93], &input[94]);
	CAS<Q, 6> (&input[95], &input[96]);
	CAS<Q, 6> (&input[97], &input[98]);
	CAS<Q, 6> (&input[99], &input[100]);
	CAS<Q, 6> (&input[101], &input[102]);
	CAS<Q, 6> (&input[103], &input[104]);
	CAS<Q, 6> (&input[105], &input[106]);
	CAS<Q, 6> (&input[107], &input[108]);
	CAS<Q, 6> (&input[109], &input[110]);
	CAS<Q, 6> (&input[111], &input[112]);
	CAS<Q, 6> (&input[113], &input[114]);
	CAS<Q, 6> (&input[115], &input[116]);
	//12
	CAS<Q, 6> (&input[12], &input[13]);
	CAS<Q, 6> (&input[14], &input[15]);
	CAS<Q, 6> (&input[16], &input[17]);
	CAS<Q, 6> (&input[18], &input[19]);
	CAS<Q, 6> (&input[20], &input[21]);
	CAS<Q, 6> (&input[22], &input[23]);
	CAS<Q, 6> (&input[24], &input[25]);
	CAS<Q, 6> (&input[26], &input[27]);
	CAS<Q, 6> (&input[28], &input[29]);
	CAS<Q, 6> (&input[30], &input[31]);
	CAS<Q, 6> (&input[32], &input[33]);
	CAS<Q, 6> (&input[34], &input[35]);
	CAS<Q, 6> (&input[36], &input[37]);
	CAS<Q, 6> (&input[38], &input[39]);
	CAS<Q, 6> (&input[40], &input[41]);
	CAS<Q, 6> (&input[42], &input[43]);
	CAS<Q, 6> (&input[44], &input[45]);
	CAS<Q, 6> (&input[46], &input[47]);
	CAS<Q, 6> (&input[48], &input[49]);
	CAS<Q, 6> (&input[50], &input[51]);
	CAS<Q, 6> (&input[52], &input[53]);
	CAS<Q, 6> (&input[54], &input[55]);
	CAS<Q, 6> (&input[56], &input[57]);
	CAS<Q, 6> (&input[58], &input[59]);
	CAS<Q, 6> (&input[60], &input[61]);
	CAS<Q, 6> (&input[62], &input[63]);
	CAS<Q, 6> (&input[64], &input[65]);
	CAS<Q, 6> (&input[66], &input[67]);
	CAS<Q, 6> (&input[68], &input[69]);
	CAS<Q, 6> (&input[70], &input[71]);
	CAS<Q, 6> (&input[72], &input[73]);
	CAS<Q, 6> (&input[74], &input[75]);
	CAS<Q, 6> (&input[76], &input[77]);
	CAS<Q, 6> (&input[78], &input[79]);
	CAS<Q, 6> (&input[80], &input[81]);
	CAS<Q, 6> (&input[82], &input[83]);
	CAS<Q, 6> (&input[84], &input[85]);
	CAS<Q, 6> (&input[86], &input[87]);
	CAS<Q, 6> (&input[88], &input[89]);
	CAS<Q, 6> (&input[90], &input[91]);
	CAS<Q, 6> (&input[92], &input[93]);
	CAS<Q, 6> (&input[94], &input[95]);
	CAS<Q, 6> (&input[96], &input[97]);
	CAS<Q, 6> (&input[98], &input[99]);
	CAS<Q, 6> (&input[100], &input[101]);
	CAS<Q, 6> (&input[102], &input[103]);
	CAS<Q, 6> (&input[104], &input[105]);
	CAS<Q, 6> (&input[106], &input[107]);
	CAS<Q, 6> (&input[108], &input[109]);
	CAS<Q, 6> (&input[110], &input[111]);
	CAS<Q, 6> (&input[112], &input[113]);
	CAS<Q, 6> (&input[114], &input[115]);
	//13
	CAS<Q, 6> (&input[13], &input[14]);
	CAS<Q, 6> (&input[15], &input[16]);
	CAS<Q, 6> (&input[17], &input[18]);
	CAS<Q, 6> (&input[19], &input[20]);
	CAS<Q, 6> (&input[21], &input[22]);
	CAS<Q, 6> (&input[23], &input[24]);
	CAS<Q, 6> (&input[25], &input[26]);
	CAS<Q, 6> (&input[27], &input[28]);
	CAS<Q, 6> (&input[29], &input[30]);
	CAS<Q, 6> (&input[31], &input[32]);
	CAS<Q, 6> (&input[33], &input[34]);
	CAS<Q, 6> (&input[35], &input[36]);
	CAS<Q, 6> (&input[37], &input[38]);
	CAS<Q, 6> (&input[39], &input[40]);
	CAS<Q, 6> (&input[41], &input[42]);
	CAS<Q, 6> (&input[43], &input[44]);
	CAS<Q, 6> (&input[45], &input[46]);
	CAS<Q, 6> (&input[47], &input[48]);
	CAS<Q, 6> (&input[49], &input[50]);
	CAS<Q, 6> (&input[51], &input[52]);
	CAS<Q, 6> (&input[53], &input[54]);
	CAS<Q, 6> (&input[55], &input[56]);
	CAS<Q, 6> (&input[57], &input[58]);
	CAS<Q, 6> (&input[59], &input[60]);
	CAS<Q, 6> (&input[61], &input[62]);
	CAS<Q, 6> (&input[63], &input[64]);
	CAS<Q, 6> (&input[65], &input[66]);
	CAS<Q, 6> (&input[67], &input[68]);
	CAS<Q, 6> (&input[69], &input[70]);
	CAS<Q, 6> (&input[71], &input[72]);
	CAS<Q, 6> (&input[73], &input[74]);
	CAS<Q, 6> (&input[75], &input[76]);
	CAS<Q, 6> (&input[77], &input[78]);
	CAS<Q, 6> (&input[79], &input[80]);
	CAS<Q, 6> (&input[81], &input[82]);
	CAS<Q, 6> (&input[83], &input[84]);
	CAS<Q, 6> (&input[85], &input[86]);
	CAS<Q, 6> (&input[87], &input[88]);
	CAS<Q, 6> (&input[89], &input[90]);
	CAS<Q, 6> (&input[81], &input[92]);
	CAS<Q, 6> (&input[93], &input[94]);
	CAS<Q, 6> (&input[95], &input[96]);
	CAS<Q, 6> (&input[97], &input[98]);
	CAS<Q, 6> (&input[99], &input[100]);
	CAS<Q, 6> (&input[101], &input[102]);
	CAS<Q, 6> (&input[103], &input[104]);
	CAS<Q, 6> (&input[105], &input[106]);
	CAS<Q, 6> (&input[107], &input[108]);
	CAS<Q, 6> (&input[109], &input[110]);
	CAS<Q, 6> (&input[111], &input[112]);
	CAS<Q, 6> (&input[113], &input[114]);
	//14
	CAS<Q, 6> (&input[14], &input[15]);
	CAS<Q, 6> (&input[16], &input[17]);
	CAS<Q, 6> (&input[18], &input[19]);
	CAS<Q, 6> (&input[20], &input[21]);
	CAS<Q, 6> (&input[22], &input[23]);
	CAS<Q, 6> (&input[24], &input[25]);
	CAS<Q, 6> (&input[26], &input[27]);
	CAS<Q, 6> (&input[28], &input[29]);
	CAS<Q, 6> (&input[30], &input[31]);
	CAS<Q, 6> (&input[32], &input[33]);
	CAS<Q, 6> (&input[34], &input[35]);
	CAS<Q, 6> (&input[36], &input[37]);
	CAS<Q, 6> (&input[38], &input[39]);
	CAS<Q, 6> (&input[40], &input[41]);
	CAS<Q, 6> (&input[42], &input[43]);
	CAS<Q, 6> (&input[44], &input[45]);
	CAS<Q, 6> (&input[46], &input[47]);
	CAS<Q, 6> (&input[48], &input[49]);
	CAS<Q, 6> (&input[50], &input[51]);
	CAS<Q, 6> (&input[52], &input[53]);
	CAS<Q, 6> (&input[54], &input[55]);
	CAS<Q, 6> (&input[56], &input[57]);
	CAS<Q, 6> (&input[58], &input[59]);
	CAS<Q, 6> (&input[60], &input[61]);
	CAS<Q, 6> (&input[62], &input[63]);
	CAS<Q, 6> (&input[64], &input[65]);
	CAS<Q, 6> (&input[66], &input[67]);
	CAS<Q, 6> (&input[68], &input[69]);
	CAS<Q, 6> (&input[70], &input[71]);
	CAS<Q, 6> (&input[72], &input[73]);
	CAS<Q, 6> (&input[74], &input[75]);
	CAS<Q, 6> (&input[76], &input[77]);
	CAS<Q, 6> (&input[78], &input[79]);
	CAS<Q, 6> (&input[80], &input[81]);
	CAS<Q, 6> (&input[82], &input[83]);
	CAS<Q, 6> (&input[84], &input[85]);
	CAS<Q, 6> (&input[86], &input[87]);
	CAS<Q, 6> (&input[88], &input[89]);
	CAS<Q, 6> (&input[90], &input[91]);
	CAS<Q, 6> (&input[92], &input[93]);
	CAS<Q, 6> (&input[94], &input[95]);
	CAS<Q, 6> (&input[96], &input[97]);
	CAS<Q, 6> (&input[98], &input[99]);
	CAS<Q, 6> (&input[100], &input[101]);
	CAS<Q, 6> (&input[102], &input[103]);
	CAS<Q, 6> (&input[104], &input[105]);
	CAS<Q, 6> (&input[106], &input[107]);
	CAS<Q, 6> (&input[108], &input[109]);
	CAS<Q, 6> (&input[110], &input[111]);
	CAS<Q, 6> (&input[112], &input[113]);
	//15
	CAS<Q, 6> (&input[15], &input[16]);
	CAS<Q, 6> (&input[17], &input[18]);
	CAS<Q, 6> (&input[19], &input[20]);
	CAS<Q, 6> (&input[21], &input[22]);
	CAS<Q, 6> (&input[23], &input[24]);
	CAS<Q, 6> (&input[25], &input[26]);
	CAS<Q, 6> (&input[27], &input[28]);
	CAS<Q, 6> (&input[29], &input[30]);
	CAS<Q, 6> (&input[31], &input[32]);
	CAS<Q, 6> (&input[33], &input[34]);
	CAS<Q, 6> (&input[35], &input[36]);
	CAS<Q, 6> (&input[37], &input[38]);
	CAS<Q, 6> (&input[39], &input[40]);
	CAS<Q, 6> (&input[41], &input[42]);
	CAS<Q, 6> (&input[43], &input[44]);
	CAS<Q, 6> (&input[45], &input[46]);
	CAS<Q, 6> (&input[47], &input[48]);
	CAS<Q, 6> (&input[49], &input[50]);
	CAS<Q, 6> (&input[51], &input[52]);
	CAS<Q, 6> (&input[53], &input[54]);
	CAS<Q, 6> (&input[55], &input[56]);
	CAS<Q, 6> (&input[57], &input[58]);
	CAS<Q, 6> (&input[59], &input[60]);
	CAS<Q, 6> (&input[61], &input[62]);
	CAS<Q, 6> (&input[63], &input[64]);
	CAS<Q, 6> (&input[65], &input[66]);
	CAS<Q, 6> (&input[67], &input[68]);
	CAS<Q, 6> (&input[69], &input[70]);
	CAS<Q, 6> (&input[71], &input[72]);
	CAS<Q, 6> (&input[73], &input[74]);
	CAS<Q, 6> (&input[75], &input[76]);
	CAS<Q, 6> (&input[77], &input[78]);
	CAS<Q, 6> (&input[79], &input[80]);
	CAS<Q, 6> (&input[81], &input[82]);
	CAS<Q, 6> (&input[83], &input[84]);
	CAS<Q, 6> (&input[85], &input[86]);
	CAS<Q, 6> (&input[87], &input[88]);
	CAS<Q, 6> (&input[89], &input[90]);
	CAS<Q, 6> (&input[81], &input[92]);
	CAS<Q, 6> (&input[93], &input[94]);
	CAS<Q, 6> (&input[95], &input[96]);
	CAS<Q, 6> (&input[97], &input[98]);
	CAS<Q, 6> (&input[99], &input[100]);
	CAS<Q, 6> (&input[101], &input[102]);
	CAS<Q, 6> (&input[103], &input[104]);
	CAS<Q, 6> (&input[105], &input[106]);
	CAS<Q, 6> (&input[107], &input[108]);
	CAS<Q, 6> (&input[109], &input[110]);
	CAS<Q, 6> (&input[111], &input[112]);
	//16
	CAS<Q, 6> (&input[16], &input[17]);
	CAS<Q, 6> (&input[18], &input[19]);
	CAS<Q, 6> (&input[20], &input[21]);
	CAS<Q, 6> (&input[22], &input[23]);
	CAS<Q, 6> (&input[24], &input[25]);
	CAS<Q, 6> (&input[26], &input[27]);
	CAS<Q, 6> (&input[28], &input[29]);
	CAS<Q, 6> (&input[30], &input[31]);
	CAS<Q, 6> (&input[32], &input[33]);
	CAS<Q, 6> (&input[34], &input[35]);
	CAS<Q, 6> (&input[36], &input[37]);
	CAS<Q, 6> (&input[38], &input[39]);
	CAS<Q, 6> (&input[40], &input[41]);
	CAS<Q, 6> (&input[42], &input[43]);
	CAS<Q, 6> (&input[44], &input[45]);
	CAS<Q, 6> (&input[46], &input[47]);
	CAS<Q, 6> (&input[48], &input[49]);
	CAS<Q, 6> (&input[50], &input[51]);
	CAS<Q, 6> (&input[52], &input[53]);
	CAS<Q, 6> (&input[54], &input[55]);
	CAS<Q, 6> (&input[56], &input[57]);
	CAS<Q, 6> (&input[58], &input[59]);
	CAS<Q, 6> (&input[60], &input[61]);
	CAS<Q, 6> (&input[62], &input[63]);
	CAS<Q, 6> (&input[64], &input[65]);
	CAS<Q, 6> (&input[66], &input[67]);
	CAS<Q, 6> (&input[68], &input[69]);
	CAS<Q, 6> (&input[70], &input[71]);
	CAS<Q, 6> (&input[72], &input[73]);
	CAS<Q, 6> (&input[74], &input[75]);
	CAS<Q, 6> (&input[76], &input[77]);
	CAS<Q, 6> (&input[78], &input[79]);
	CAS<Q, 6> (&input[80], &input[81]);
	CAS<Q, 6> (&input[82], &input[83]);
	CAS<Q, 6> (&input[84], &input[85]);
	CAS<Q, 6> (&input[86], &input[87]);
	CAS<Q, 6> (&input[88], &input[89]);
	CAS<Q, 6> (&input[90], &input[91]);
	CAS<Q, 6> (&input[92], &input[93]);
	CAS<Q, 6> (&input[94], &input[95]);
	CAS<Q, 6> (&input[96], &input[97]);
	CAS<Q, 6> (&input[98], &input[99]);
	CAS<Q, 6> (&input[100], &input[101]);
	CAS<Q, 6> (&input[102], &input[103]);
	CAS<Q, 6> (&input[104], &input[105]);
	CAS<Q, 6> (&input[106], &input[107]);
	CAS<Q, 6> (&input[108], &input[109]);
	CAS<Q, 6> (&input[110], &input[111]);
	//17
	CAS<Q, 6> (&input[17], &input[18]);
	CAS<Q, 6> (&input[19], &input[20]);
	CAS<Q, 6> (&input[21], &input[22]);
	CAS<Q, 6> (&input[23], &input[24]);
	CAS<Q, 6> (&input[25], &input[26]);
	CAS<Q, 6> (&input[27], &input[28]);
	CAS<Q, 6> (&input[29], &input[30]);
	CAS<Q, 6> (&input[31], &input[32]);
	CAS<Q, 6> (&input[33], &input[34]);
	CAS<Q, 6> (&input[35], &input[36]);
	CAS<Q, 6> (&input[37], &input[38]);
	CAS<Q, 6> (&input[39], &input[40]);
	CAS<Q, 6> (&input[41], &input[42]);
	CAS<Q, 6> (&input[43], &input[44]);
	CAS<Q, 6> (&input[45], &input[46]);
	CAS<Q, 6> (&input[47], &input[48]);
	CAS<Q, 6> (&input[49], &input[50]);
	CAS<Q, 6> (&input[51], &input[52]);
	CAS<Q, 6> (&input[53], &input[54]);
	CAS<Q, 6> (&input[55], &input[56]);
	CAS<Q, 6> (&input[57], &input[58]);
	CAS<Q, 6> (&input[59], &input[60]);
	CAS<Q, 6> (&input[61], &input[62]);
	CAS<Q, 6> (&input[63], &input[64]);
	CAS<Q, 6> (&input[65], &input[66]);
	CAS<Q, 6> (&input[67], &input[68]);
	CAS<Q, 6> (&input[69], &input[70]);
	CAS<Q, 6> (&input[71], &input[72]);
	CAS<Q, 6> (&input[73], &input[74]);
	CAS<Q, 6> (&input[75], &input[76]);
	CAS<Q, 6> (&input[77], &input[78]);
	CAS<Q, 6> (&input[79], &input[80]);
	CAS<Q, 6> (&input[81], &input[82]);
	CAS<Q, 6> (&input[83], &input[84]);
	CAS<Q, 6> (&input[85], &input[86]);
	CAS<Q, 6> (&input[87], &input[88]);
	CAS<Q, 6> (&input[89], &input[90]);
	CAS<Q, 6> (&input[81], &input[92]);
	CAS<Q, 6> (&input[93], &input[94]);
	CAS<Q, 6> (&input[95], &input[96]);
	CAS<Q, 6> (&input[97], &input[98]);
	CAS<Q, 6> (&input[99], &input[100]);
	CAS<Q, 6> (&input[101], &input[102]);
	CAS<Q, 6> (&input[103], &input[104]);
	CAS<Q, 6> (&input[105], &input[106]);
	CAS<Q, 6> (&input[107], &input[108]);
	CAS<Q, 6> (&input[109], &input[110]);
	//18
	CAS<Q, 6> (&input[18], &input[19]);
	CAS<Q, 6> (&input[20], &input[21]);
	CAS<Q, 6> (&input[22], &input[23]);
	CAS<Q, 6> (&input[24], &input[25]);
	CAS<Q, 6> (&input[26], &input[27]);
	CAS<Q, 6> (&input[28], &input[29]);
	CAS<Q, 6> (&input[30], &input[31]);
	CAS<Q, 6> (&input[32], &input[33]);
	CAS<Q, 6> (&input[34], &input[35]);
	CAS<Q, 6> (&input[36], &input[37]);
	CAS<Q, 6> (&input[38], &input[39]);
	CAS<Q, 6> (&input[40], &input[41]);
	CAS<Q, 6> (&input[42], &input[43]);
	CAS<Q, 6> (&input[44], &input[45]);
	CAS<Q, 6> (&input[46], &input[47]);
	CAS<Q, 6> (&input[48], &input[49]);
	CAS<Q, 6> (&input[50], &input[51]);
	CAS<Q, 6> (&input[52], &input[53]);
	CAS<Q, 6> (&input[54], &input[55]);
	CAS<Q, 6> (&input[56], &input[57]);
	CAS<Q, 6> (&input[58], &input[59]);
	CAS<Q, 6> (&input[60], &input[61]);
	CAS<Q, 6> (&input[62], &input[63]);
	CAS<Q, 6> (&input[64], &input[65]);
	CAS<Q, 6> (&input[66], &input[67]);
	CAS<Q, 6> (&input[68], &input[69]);
	CAS<Q, 6> (&input[70], &input[71]);
	CAS<Q, 6> (&input[72], &input[73]);
	CAS<Q, 6> (&input[74], &input[75]);
	CAS<Q, 6> (&input[76], &input[77]);
	CAS<Q, 6> (&input[78], &input[79]);
	CAS<Q, 6> (&input[80], &input[81]);
	CAS<Q, 6> (&input[82], &input[83]);
	CAS<Q, 6> (&input[84], &input[85]);
	CAS<Q, 6> (&input[86], &input[87]);
	CAS<Q, 6> (&input[88], &input[89]);
	CAS<Q, 6> (&input[90], &input[91]);
	CAS<Q, 6> (&input[92], &input[93]);
	CAS<Q, 6> (&input[94], &input[95]);
	CAS<Q, 6> (&input[96], &input[97]);
	CAS<Q, 6> (&input[98], &input[99]);
	CAS<Q, 6> (&input[100], &input[101]);
	CAS<Q, 6> (&input[102], &input[103]);
	CAS<Q, 6> (&input[104], &input[105]);
	CAS<Q, 6> (&input[106], &input[107]);
	CAS<Q, 6> (&input[108], &input[109]);
	//19
	CAS<Q, 6> (&input[19], &input[20]);
	CAS<Q, 6> (&input[21], &input[22]);
	CAS<Q, 6> (&input[23], &input[24]);
	CAS<Q, 6> (&input[25], &input[26]);
	CAS<Q, 6> (&input[27], &input[28]);
	CAS<Q, 6> (&input[29], &input[30]);
	CAS<Q, 6> (&input[31], &input[32]);
	CAS<Q, 6> (&input[33], &input[34]);
	CAS<Q, 6> (&input[35], &input[36]);
	CAS<Q, 6> (&input[37], &input[38]);
	CAS<Q, 6> (&input[39], &input[40]);
	CAS<Q, 6> (&input[41], &input[42]);
	CAS<Q, 6> (&input[43], &input[44]);
	CAS<Q, 6> (&input[45], &input[46]);
	CAS<Q, 6> (&input[47], &input[48]);
	CAS<Q, 6> (&input[49], &input[50]);
	CAS<Q, 6> (&input[51], &input[52]);
	CAS<Q, 6> (&input[53], &input[54]);
	CAS<Q, 6> (&input[55], &input[56]);
	CAS<Q, 6> (&input[57], &input[58]);
	CAS<Q, 6> (&input[59], &input[60]);
	CAS<Q, 6> (&input[61], &input[62]);
	CAS<Q, 6> (&input[63], &input[64]);
	CAS<Q, 6> (&input[65], &input[66]);
	CAS<Q, 6> (&input[67], &input[68]);
	CAS<Q, 6> (&input[69], &input[70]);
	CAS<Q, 6> (&input[71], &input[72]);
	CAS<Q, 6> (&input[73], &input[74]);
	CAS<Q, 6> (&input[75], &input[76]);
	CAS<Q, 6> (&input[77], &input[78]);
	CAS<Q, 6> (&input[79], &input[80]);
	CAS<Q, 6> (&input[81], &input[82]);
	CAS<Q, 6> (&input[83], &input[84]);
	CAS<Q, 6> (&input[85], &input[86]);
	CAS<Q, 6> (&input[87], &input[88]);
	CAS<Q, 6> (&input[89], &input[90]);
	CAS<Q, 6> (&input[81], &input[92]);
	CAS<Q, 6> (&input[93], &input[94]);
	CAS<Q, 6> (&input[95], &input[96]);
	CAS<Q, 6> (&input[97], &input[98]);
	CAS<Q, 6> (&input[99], &input[100]);
	CAS<Q, 6> (&input[101], &input[102]);
	CAS<Q, 6> (&input[103], &input[104]);
	CAS<Q, 6> (&input[105], &input[106]);
	CAS<Q, 6> (&input[107], &input[108]);
	//20
	CAS<Q, 6> (&input[20], &input[21]);
	CAS<Q, 6> (&input[22], &input[23]);
	CAS<Q, 6> (&input[24], &input[25]);
	CAS<Q, 6> (&input[26], &input[27]);
	CAS<Q, 6> (&input[28], &input[29]);
	CAS<Q, 6> (&input[30], &input[31]);
	CAS<Q, 6> (&input[32], &input[33]);
	CAS<Q, 6> (&input[34], &input[35]);
	CAS<Q, 6> (&input[36], &input[37]);
	CAS<Q, 6> (&input[38], &input[39]);
	CAS<Q, 6> (&input[40], &input[41]);
	CAS<Q, 6> (&input[42], &input[43]);
	CAS<Q, 6> (&input[44], &input[45]);
	CAS<Q, 6> (&input[46], &input[47]);
	CAS<Q, 6> (&input[48], &input[49]);
	CAS<Q, 6> (&input[50], &input[51]);
	CAS<Q, 6> (&input[52], &input[53]);
	CAS<Q, 6> (&input[54], &input[55]);
	CAS<Q, 6> (&input[56], &input[57]);
	CAS<Q, 6> (&input[58], &input[59]);
	CAS<Q, 6> (&input[60], &input[61]);
	CAS<Q, 6> (&input[62], &input[63]);
	CAS<Q, 6> (&input[64], &input[65]);
	CAS<Q, 6> (&input[66], &input[67]);
	CAS<Q, 6> (&input[68], &input[69]);
	CAS<Q, 6> (&input[70], &input[71]);
	CAS<Q, 6> (&input[72], &input[73]);
	CAS<Q, 6> (&input[74], &input[75]);
	CAS<Q, 6> (&input[76], &input[77]);
	CAS<Q, 6> (&input[78], &input[79]);
	CAS<Q, 6> (&input[80], &input[81]);
	CAS<Q, 6> (&input[82], &input[83]);
	CAS<Q, 6> (&input[84], &input[85]);
	CAS<Q, 6> (&input[86], &input[87]);
	CAS<Q, 6> (&input[88], &input[89]);
	CAS<Q, 6> (&input[90], &input[91]);
	CAS<Q, 6> (&input[92], &input[93]);
	CAS<Q, 6> (&input[94], &input[95]);
	CAS<Q, 6> (&input[96], &input[97]);
	CAS<Q, 6> (&input[98], &input[99]);
	CAS<Q, 6> (&input[100], &input[101]);
	CAS<Q, 6> (&input[102], &input[103]);
	CAS<Q, 6> (&input[104], &input[105]);
	CAS<Q, 6> (&input[106], &input[107]);
	//21
	CAS<Q, 6> (&input[21], &input[22]);
	CAS<Q, 6> (&input[23], &input[24]);
	CAS<Q, 6> (&input[25], &input[26]);
	CAS<Q, 6> (&input[27], &input[28]);
	CAS<Q, 6> (&input[29], &input[30]);
	CAS<Q, 6> (&input[31], &input[32]);
	CAS<Q, 6> (&input[33], &input[34]);
	CAS<Q, 6> (&input[35], &input[36]);
	CAS<Q, 6> (&input[37], &input[38]);
	CAS<Q, 6> (&input[39], &input[40]);
	CAS<Q, 6> (&input[41], &input[42]);
	CAS<Q, 6> (&input[43], &input[44]);
	CAS<Q, 6> (&input[45], &input[46]);
	CAS<Q, 6> (&input[47], &input[48]);
	CAS<Q, 6> (&input[49], &input[50]);
	CAS<Q, 6> (&input[51], &input[52]);
	CAS<Q, 6> (&input[53], &input[54]);
	CAS<Q, 6> (&input[55], &input[56]);
	CAS<Q, 6> (&input[57], &input[58]);
	CAS<Q, 6> (&input[59], &input[60]);
	CAS<Q, 6> (&input[61], &input[62]);
	CAS<Q, 6> (&input[63], &input[64]);
	CAS<Q, 6> (&input[65], &input[66]);
	CAS<Q, 6> (&input[67], &input[68]);
	CAS<Q, 6> (&input[69], &input[70]);
	CAS<Q, 6> (&input[71], &input[72]);
	CAS<Q, 6> (&input[73], &input[74]);
	CAS<Q, 6> (&input[75], &input[76]);
	CAS<Q, 6> (&input[77], &input[78]);
	CAS<Q, 6> (&input[79], &input[80]);
	CAS<Q, 6> (&input[81], &input[82]);
	CAS<Q, 6> (&input[83], &input[84]);
	CAS<Q, 6> (&input[85], &input[86]);
	CAS<Q, 6> (&input[87], &input[88]);
	CAS<Q, 6> (&input[89], &input[90]);
	CAS<Q, 6> (&input[81], &input[92]);
	CAS<Q, 6> (&input[93], &input[94]);
	CAS<Q, 6> (&input[95], &input[96]);
	CAS<Q, 6> (&input[97], &input[98]);
	CAS<Q, 6> (&input[99], &input[100]);
	CAS<Q, 6> (&input[101], &input[102]);
	CAS<Q, 6> (&input[103], &input[104]);
	CAS<Q, 6> (&input[105], &input[106]);
	//22
	CAS<Q, 6> (&input[22], &input[23]);
	CAS<Q, 6> (&input[24], &input[25]);
	CAS<Q, 6> (&input[26], &input[27]);
	CAS<Q, 6> (&input[28], &input[29]);
	CAS<Q, 6> (&input[30], &input[31]);
	CAS<Q, 6> (&input[32], &input[33]);
	CAS<Q, 6> (&input[34], &input[35]);
	CAS<Q, 6> (&input[36], &input[37]);
	CAS<Q, 6> (&input[38], &input[39]);
	CAS<Q, 6> (&input[40], &input[41]);
	CAS<Q, 6> (&input[42], &input[43]);
	CAS<Q, 6> (&input[44], &input[45]);
	CAS<Q, 6> (&input[46], &input[47]);
	CAS<Q, 6> (&input[48], &input[49]);
	CAS<Q, 6> (&input[50], &input[51]);
	CAS<Q, 6> (&input[52], &input[53]);
	CAS<Q, 6> (&input[54], &input[55]);
	CAS<Q, 6> (&input[56], &input[57]);
	CAS<Q, 6> (&input[58], &input[59]);
	CAS<Q, 6> (&input[60], &input[61]);
	CAS<Q, 6> (&input[62], &input[63]);
	CAS<Q, 6> (&input[64], &input[65]);
	CAS<Q, 6> (&input[66], &input[67]);
	CAS<Q, 6> (&input[68], &input[69]);
	CAS<Q, 6> (&input[70], &input[71]);
	CAS<Q, 6> (&input[72], &input[73]);
	CAS<Q, 6> (&input[74], &input[75]);
	CAS<Q, 6> (&input[76], &input[77]);
	CAS<Q, 6> (&input[78], &input[79]);
	CAS<Q, 6> (&input[80], &input[81]);
	CAS<Q, 6> (&input[82], &input[83]);
	CAS<Q, 6> (&input[84], &input[85]);
	CAS<Q, 6> (&input[86], &input[87]);
	CAS<Q, 6> (&input[88], &input[89]);
	CAS<Q, 6> (&input[90], &input[91]);
	CAS<Q, 6> (&input[92], &input[93]);
	CAS<Q, 6> (&input[94], &input[95]);
	CAS<Q, 6> (&input[96], &input[97]);
	CAS<Q, 6> (&input[98], &input[99]);
	CAS<Q, 6> (&input[100], &input[101]);
	CAS<Q, 6> (&input[102], &input[103]);
	CAS<Q, 6> (&input[104], &input[105]);
	//23
	CAS<Q, 6> (&input[23], &input[24]);
	CAS<Q, 6> (&input[25], &input[26]);
	CAS<Q, 6> (&input[27], &input[28]);
	CAS<Q, 6> (&input[29], &input[30]);
	CAS<Q, 6> (&input[31], &input[32]);
	CAS<Q, 6> (&input[33], &input[34]);
	CAS<Q, 6> (&input[35], &input[36]);
	CAS<Q, 6> (&input[37], &input[38]);
	CAS<Q, 6> (&input[39], &input[40]);
	CAS<Q, 6> (&input[41], &input[42]);
	CAS<Q, 6> (&input[43], &input[44]);
	CAS<Q, 6> (&input[45], &input[46]);
	CAS<Q, 6> (&input[47], &input[48]);
	CAS<Q, 6> (&input[49], &input[50]);
	CAS<Q, 6> (&input[51], &input[52]);
	CAS<Q, 6> (&input[53], &input[54]);
	CAS<Q, 6> (&input[55], &input[56]);
	CAS<Q, 6> (&input[57], &input[58]);
	CAS<Q, 6> (&input[59], &input[60]);
	CAS<Q, 6> (&input[61], &input[62]);
	CAS<Q, 6> (&input[63], &input[64]);
	CAS<Q, 6> (&input[65], &input[66]);
	CAS<Q, 6> (&input[67], &input[68]);
	CAS<Q, 6> (&input[69], &input[70]);
	CAS<Q, 6> (&input[71], &input[72]);
	CAS<Q, 6> (&input[73], &input[74]);
	CAS<Q, 6> (&input[75], &input[76]);
	CAS<Q, 6> (&input[77], &input[78]);
	CAS<Q, 6> (&input[79], &input[80]);
	CAS<Q, 6> (&input[81], &input[82]);
	CAS<Q, 6> (&input[83], &input[84]);
	CAS<Q, 6> (&input[85], &input[86]);
	CAS<Q, 6> (&input[87], &input[88]);
	CAS<Q, 6> (&input[89], &input[90]);
	CAS<Q, 6> (&input[81], &input[92]);
	CAS<Q, 6> (&input[93], &input[94]);
	CAS<Q, 6> (&input[95], &input[96]);
	CAS<Q, 6> (&input[97], &input[98]);
	CAS<Q, 6> (&input[99], &input[100]);
	CAS<Q, 6> (&input[101], &input[102]);
	CAS<Q, 6> (&input[103], &input[104]);
	//24
	CAS<Q, 6> (&input[24], &input[25]);
	CAS<Q, 6> (&input[26], &input[27]);
	CAS<Q, 6> (&input[28], &input[29]);
	CAS<Q, 6> (&input[30], &input[31]);
	CAS<Q, 6> (&input[32], &input[33]);
	CAS<Q, 6> (&input[34], &input[35]);
	CAS<Q, 6> (&input[36], &input[37]);
	CAS<Q, 6> (&input[38], &input[39]);
	CAS<Q, 6> (&input[40], &input[41]);
	CAS<Q, 6> (&input[42], &input[43]);
	CAS<Q, 6> (&input[44], &input[45]);
	CAS<Q, 6> (&input[46], &input[47]);
	CAS<Q, 6> (&input[48], &input[49]);
	CAS<Q, 6> (&input[50], &input[51]);
	CAS<Q, 6> (&input[52], &input[53]);
	CAS<Q, 6> (&input[54], &input[55]);
	CAS<Q, 6> (&input[56], &input[57]);
	CAS<Q, 6> (&input[58], &input[59]);
	CAS<Q, 6> (&input[60], &input[61]);
	CAS<Q, 6> (&input[62], &input[63]);
	CAS<Q, 6> (&input[64], &input[65]);
	CAS<Q, 6> (&input[66], &input[67]);
	CAS<Q, 6> (&input[68], &input[69]);
	CAS<Q, 6> (&input[70], &input[71]);
	CAS<Q, 6> (&input[72], &input[73]);
	CAS<Q, 6> (&input[74], &input[75]);
	CAS<Q, 6> (&input[76], &input[77]);
	CAS<Q, 6> (&input[78], &input[79]);
	CAS<Q, 6> (&input[80], &input[81]);
	CAS<Q, 6> (&input[82], &input[83]);
	CAS<Q, 6> (&input[84], &input[85]);
	CAS<Q, 6> (&input[86], &input[87]);
	CAS<Q, 6> (&input[88], &input[89]);
	CAS<Q, 6> (&input[90], &input[91]);
	CAS<Q, 6> (&input[92], &input[93]);
	CAS<Q, 6> (&input[94], &input[95]);
	CAS<Q, 6> (&input[96], &input[97]);
	CAS<Q, 6> (&input[98], &input[99]);
	CAS<Q, 6> (&input[100], &input[101]);
	CAS<Q, 6> (&input[102], &input[103]);
	//25
	CAS<Q, 6> (&input[25], &input[26]);
	CAS<Q, 6> (&input[27], &input[28]);
	CAS<Q, 6> (&input[29], &input[30]);
	CAS<Q, 6> (&input[31], &input[32]);
	CAS<Q, 6> (&input[33], &input[34]);
	CAS<Q, 6> (&input[35], &input[36]);
	CAS<Q, 6> (&input[37], &input[38]);
	CAS<Q, 6> (&input[39], &input[40]);
	CAS<Q, 6> (&input[41], &input[42]);
	CAS<Q, 6> (&input[43], &input[44]);
	CAS<Q, 6> (&input[45], &input[46]);
	CAS<Q, 6> (&input[47], &input[48]);
	CAS<Q, 6> (&input[49], &input[50]);
	CAS<Q, 6> (&input[51], &input[52]);
	CAS<Q, 6> (&input[53], &input[54]);
	CAS<Q, 6> (&input[55], &input[56]);
	CAS<Q, 6> (&input[57], &input[58]);
	CAS<Q, 6> (&input[59], &input[60]);
	CAS<Q, 6> (&input[61], &input[62]);
	CAS<Q, 6> (&input[63], &input[64]);
	CAS<Q, 6> (&input[65], &input[66]);
	CAS<Q, 6> (&input[67], &input[68]);
	CAS<Q, 6> (&input[69], &input[70]);
	CAS<Q, 6> (&input[71], &input[72]);
	CAS<Q, 6> (&input[73], &input[74]);
	CAS<Q, 6> (&input[75], &input[76]);
	CAS<Q, 6> (&input[77], &input[78]);
	CAS<Q, 6> (&input[79], &input[80]);
	CAS<Q, 6> (&input[81], &input[82]);
	CAS<Q, 6> (&input[83], &input[84]);
	CAS<Q, 6> (&input[85], &input[86]);
	CAS<Q, 6> (&input[87], &input[88]);
	CAS<Q, 6> (&input[89], &input[90]);
	CAS<Q, 6> (&input[81], &input[92]);
	CAS<Q, 6> (&input[93], &input[94]);
	CAS<Q, 6> (&input[95], &input[96]);
	CAS<Q, 6> (&input[97], &input[98]);
	CAS<Q, 6> (&input[99], &input[100]);
	CAS<Q, 6> (&input[101], &input[102]);
	//26
	CAS<Q, 6> (&input[26], &input[27]);
	CAS<Q, 6> (&input[28], &input[29]);
	CAS<Q, 6> (&input[30], &input[31]);
	CAS<Q, 6> (&input[32], &input[33]);
	CAS<Q, 6> (&input[34], &input[35]);
	CAS<Q, 6> (&input[36], &input[37]);
	CAS<Q, 6> (&input[38], &input[39]);
	CAS<Q, 6> (&input[40], &input[41]);
	CAS<Q, 6> (&input[42], &input[43]);
	CAS<Q, 6> (&input[44], &input[45]);
	CAS<Q, 6> (&input[46], &input[47]);
	CAS<Q, 6> (&input[48], &input[49]);
	CAS<Q, 6> (&input[50], &input[51]);
	CAS<Q, 6> (&input[52], &input[53]);
	CAS<Q, 6> (&input[54], &input[55]);
	CAS<Q, 6> (&input[56], &input[57]);
	CAS<Q, 6> (&input[58], &input[59]);
	CAS<Q, 6> (&input[60], &input[61]);
	CAS<Q, 6> (&input[62], &input[63]);
	CAS<Q, 6> (&input[64], &input[65]);
	CAS<Q, 6> (&input[66], &input[67]);
	CAS<Q, 6> (&input[68], &input[69]);
	CAS<Q, 6> (&input[70], &input[71]);
	CAS<Q, 6> (&input[72], &input[73]);
	CAS<Q, 6> (&input[74], &input[75]);
	CAS<Q, 6> (&input[76], &input[77]);
	CAS<Q, 6> (&input[78], &input[79]);
	CAS<Q, 6> (&input[80], &input[81]);
	CAS<Q, 6> (&input[82], &input[83]);
	CAS<Q, 6> (&input[84], &input[85]);
	CAS<Q, 6> (&input[86], &input[87]);
	CAS<Q, 6> (&input[88], &input[89]);
	CAS<Q, 6> (&input[90], &input[91]);
	CAS<Q, 6> (&input[92], &input[93]);
	CAS<Q, 6> (&input[94], &input[95]);
	CAS<Q, 6> (&input[96], &input[97]);
	CAS<Q, 6> (&input[98], &input[99]);
	CAS<Q, 6> (&input[100], &input[101]);
	//27
	CAS<Q, 6> (&input[27], &input[28]);
	CAS<Q, 6> (&input[29], &input[30]);
	CAS<Q, 6> (&input[31], &input[32]);
	CAS<Q, 6> (&input[33], &input[34]);
	CAS<Q, 6> (&input[35], &input[36]);
	CAS<Q, 6> (&input[37], &input[38]);
	CAS<Q, 6> (&input[39], &input[40]);
	CAS<Q, 6> (&input[41], &input[42]);
	CAS<Q, 6> (&input[43], &input[44]);
	CAS<Q, 6> (&input[45], &input[46]);
	CAS<Q, 6> (&input[47], &input[48]);
	CAS<Q, 6> (&input[49], &input[50]);
	CAS<Q, 6> (&input[51], &input[52]);
	CAS<Q, 6> (&input[53], &input[54]);
	CAS<Q, 6> (&input[55], &input[56]);
	CAS<Q, 6> (&input[57], &input[58]);
	CAS<Q, 6> (&input[59], &input[60]);
	CAS<Q, 6> (&input[61], &input[62]);
	CAS<Q, 6> (&input[63], &input[64]);
	CAS<Q, 6> (&input[65], &input[66]);
	CAS<Q, 6> (&input[67], &input[68]);
	CAS<Q, 6> (&input[69], &input[70]);
	CAS<Q, 6> (&input[71], &input[72]);
	CAS<Q, 6> (&input[73], &input[74]);
	CAS<Q, 6> (&input[75], &input[76]);
	CAS<Q, 6> (&input[77], &input[78]);
	CAS<Q, 6> (&input[79], &input[80]);
	CAS<Q, 6> (&input[81], &input[82]);
	CAS<Q, 6> (&input[83], &input[84]);
	CAS<Q, 6> (&input[85], &input[86]);
	CAS<Q, 6> (&input[87], &input[88]);
	CAS<Q, 6> (&input[89], &input[90]);
	CAS<Q, 6> (&input[81], &input[92]);
	CAS<Q, 6> (&input[93], &input[94]);
	CAS<Q, 6> (&input[95], &input[96]);
	CAS<Q, 6> (&input[97], &input[98]);
	CAS<Q, 6> (&input[99], &input[100]);
	//28
	CAS<Q, 6> (&input[28], &input[29]);
	CAS<Q, 6> (&input[30], &input[31]);
	CAS<Q, 6> (&input[32], &input[33]);
	CAS<Q, 6> (&input[34], &input[35]);
	CAS<Q, 6> (&input[36], &input[37]);
	CAS<Q, 6> (&input[38], &input[39]);
	CAS<Q, 6> (&input[40], &input[41]);
	CAS<Q, 6> (&input[42], &input[43]);
	CAS<Q, 6> (&input[44], &input[45]);
	CAS<Q, 6> (&input[46], &input[47]);
	CAS<Q, 6> (&input[48], &input[49]);
	CAS<Q, 6> (&input[50], &input[51]);
	CAS<Q, 6> (&input[52], &input[53]);
	CAS<Q, 6> (&input[54], &input[55]);
	CAS<Q, 6> (&input[56], &input[57]);
	CAS<Q, 6> (&input[58], &input[59]);
	CAS<Q, 6> (&input[60], &input[61]);
	CAS<Q, 6> (&input[62], &input[63]);
	CAS<Q, 6> (&input[64], &input[65]);
	CAS<Q, 6> (&input[66], &input[67]);
	CAS<Q, 6> (&input[68], &input[69]);
	CAS<Q, 6> (&input[70], &input[71]);
	CAS<Q, 6> (&input[72], &input[73]);
	CAS<Q, 6> (&input[74], &input[75]);
	CAS<Q, 6> (&input[76], &input[77]);
	CAS<Q, 6> (&input[78], &input[79]);
	CAS<Q, 6> (&input[80], &input[81]);
	CAS<Q, 6> (&input[82], &input[83]);
	CAS<Q, 6> (&input[84], &input[85]);
	CAS<Q, 6> (&input[86], &input[87]);
	CAS<Q, 6> (&input[88], &input[89]);
	CAS<Q, 6> (&input[90], &input[91]);
	CAS<Q, 6> (&input[92], &input[93]);
	CAS<Q, 6> (&input[94], &input[95]);
	CAS<Q, 6> (&input[96], &input[97]);
	CAS<Q, 6> (&input[98], &input[99]);
	//29
	CAS<Q, 6> (&input[29], &input[30]);
	CAS<Q, 6> (&input[31], &input[32]);
	CAS<Q, 6> (&input[33], &input[34]);
	CAS<Q, 6> (&input[35], &input[36]);
	CAS<Q, 6> (&input[37], &input[38]);
	CAS<Q, 6> (&input[39], &input[40]);
	CAS<Q, 6> (&input[41], &input[42]);
	CAS<Q, 6> (&input[43], &input[44]);
	CAS<Q, 6> (&input[45], &input[46]);
	CAS<Q, 6> (&input[47], &input[48]);
	CAS<Q, 6> (&input[49], &input[50]);
	CAS<Q, 6> (&input[51], &input[52]);
	CAS<Q, 6> (&input[53], &input[54]);
	CAS<Q, 6> (&input[55], &input[56]);
	CAS<Q, 6> (&input[57], &input[58]);
	CAS<Q, 6> (&input[59], &input[60]);
	CAS<Q, 6> (&input[61], &input[62]);
	CAS<Q, 6> (&input[63], &input[64]);
	CAS<Q, 6> (&input[65], &input[66]);
	CAS<Q, 6> (&input[67], &input[68]);
	CAS<Q, 6> (&input[69], &input[70]);
	CAS<Q, 6> (&input[71], &input[72]);
	CAS<Q, 6> (&input[73], &input[74]);
	CAS<Q, 6> (&input[75], &input[76]);
	CAS<Q, 6> (&input[77], &input[78]);
	CAS<Q, 6> (&input[79], &input[80]);
	CAS<Q, 6> (&input[81], &input[82]);
	CAS<Q, 6> (&input[83], &input[84]);
	CAS<Q, 6> (&input[85], &input[86]);
	CAS<Q, 6> (&input[87], &input[88]);
	CAS<Q, 6> (&input[89], &input[90]);
	CAS<Q, 6> (&input[81], &input[92]);
	CAS<Q, 6> (&input[93], &input[94]);
	CAS<Q, 6> (&input[95], &input[96]);
	CAS<Q, 6> (&input[97], &input[98]);
	//30
	CAS<Q, 6> (&input[30], &input[31]);
	CAS<Q, 6> (&input[32], &input[33]);
	CAS<Q, 6> (&input[34], &input[35]);
	CAS<Q, 6> (&input[36], &input[37]);
	CAS<Q, 6> (&input[38], &input[39]);
	CAS<Q, 6> (&input[40], &input[41]);
	CAS<Q, 6> (&input[42], &input[43]);
	CAS<Q, 6> (&input[44], &input[45]);
	CAS<Q, 6> (&input[46], &input[47]);
	CAS<Q, 6> (&input[48], &input[49]);
	CAS<Q, 6> (&input[50], &input[51]);
	CAS<Q, 6> (&input[52], &input[53]);
	CAS<Q, 6> (&input[54], &input[55]);
	CAS<Q, 6> (&input[56], &input[57]);
	CAS<Q, 6> (&input[58], &input[59]);
	CAS<Q, 6> (&input[60], &input[61]);
	CAS<Q, 6> (&input[62], &input[63]);
	CAS<Q, 6> (&input[64], &input[65]);
	CAS<Q, 6> (&input[66], &input[67]);
	CAS<Q, 6> (&input[68], &input[69]);
	CAS<Q, 6> (&input[70], &input[71]);
	CAS<Q, 6> (&input[72], &input[73]);
	CAS<Q, 6> (&input[74], &input[75]);
	CAS<Q, 6> (&input[76], &input[77]);
	CAS<Q, 6> (&input[78], &input[79]);
	CAS<Q, 6> (&input[80], &input[81]);
	CAS<Q, 6> (&input[82], &input[83]);
	CAS<Q, 6> (&input[84], &input[85]);
	CAS<Q, 6> (&input[86], &input[87]);
	CAS<Q, 6> (&input[88], &input[89]);
	CAS<Q, 6> (&input[90], &input[91]);
	CAS<Q, 6> (&input[92], &input[93]);
	CAS<Q, 6> (&input[94], &input[95]);
	CAS<Q, 6> (&input[96], &input[97]);
	//31
	CAS<Q, 6> (&input[31], &input[32]);
	CAS<Q, 6> (&input[33], &input[34]);
	CAS<Q, 6> (&input[35], &input[36]);
	CAS<Q, 6> (&input[37], &input[38]);
	CAS<Q, 6> (&input[39], &input[40]);
	CAS<Q, 6> (&input[41], &input[42]);
	CAS<Q, 6> (&input[43], &input[44]);
	CAS<Q, 6> (&input[45], &input[46]);
	CAS<Q, 6> (&input[47], &input[48]);
	CAS<Q, 6> (&input[49], &input[50]);
	CAS<Q, 6> (&input[51], &input[52]);
	CAS<Q, 6> (&input[53], &input[54]);
	CAS<Q, 6> (&input[55], &input[56]);
	CAS<Q, 6> (&input[57], &input[58]);
	CAS<Q, 6> (&input[59], &input[60]);
	CAS<Q, 6> (&input[61], &input[62]);
	CAS<Q, 6> (&input[63], &input[64]);
	CAS<Q, 6> (&input[65], &input[66]);
	CAS<Q, 6> (&input[67], &input[68]);
	CAS<Q, 6> (&input[69], &input[70]);
	CAS<Q, 6> (&input[71], &input[72]);
	CAS<Q, 6> (&input[73], &input[74]);
	CAS<Q, 6> (&input[75], &input[76]);
	CAS<Q, 6> (&input[77], &input[78]);
	CAS<Q, 6> (&input[79], &input[80]);
	CAS<Q, 6> (&input[81], &input[82]);
	CAS<Q, 6> (&input[83], &input[84]);
	CAS<Q, 6> (&input[85], &input[86]);
	CAS<Q, 6> (&input[87], &input[88]);
	CAS<Q, 6> (&input[89], &input[90]);
	CAS<Q, 6> (&input[81], &input[92]);
	CAS<Q, 6> (&input[93], &input[94]);
	CAS<Q, 6> (&input[95], &input[96]);
	//32
	CAS<Q, 6> (&input[32], &input[33]);
	CAS<Q, 6> (&input[34], &input[35]);
	CAS<Q, 6> (&input[36], &input[37]);
	CAS<Q, 6> (&input[38], &input[39]);
	CAS<Q, 6> (&input[40], &input[41]);
	CAS<Q, 6> (&input[42], &input[43]);
	CAS<Q, 6> (&input[44], &input[45]);
	CAS<Q, 6> (&input[46], &input[47]);
	CAS<Q, 6> (&input[48], &input[49]);
	CAS<Q, 6> (&input[50], &input[51]);
	CAS<Q, 6> (&input[52], &input[53]);
	CAS<Q, 6> (&input[54], &input[55]);
	CAS<Q, 6> (&input[56], &input[57]);
	CAS<Q, 6> (&input[58], &input[59]);
	CAS<Q, 6> (&input[60], &input[61]);
	CAS<Q, 6> (&input[62], &input[63]);
	CAS<Q, 6> (&input[64], &input[65]);
	CAS<Q, 6> (&input[66], &input[67]);
	CAS<Q, 6> (&input[68], &input[69]);
	CAS<Q, 6> (&input[70], &input[71]);
	CAS<Q, 6> (&input[72], &input[73]);
	CAS<Q, 6> (&input[74], &input[75]);
	CAS<Q, 6> (&input[76], &input[77]);
	CAS<Q, 6> (&input[78], &input[79]);
	CAS<Q, 6> (&input[80], &input[81]);
	CAS<Q, 6> (&input[82], &input[83]);
	CAS<Q, 6> (&input[84], &input[85]);
	CAS<Q, 6> (&input[86], &input[87]);
	CAS<Q, 6> (&input[88], &input[89]);
	CAS<Q, 6> (&input[90], &input[91]);
	CAS<Q, 6> (&input[92], &input[93]);
	CAS<Q, 6> (&input[94], &input[95]);
	//33
	CAS<Q, 6> (&input[33], &input[34]);
	CAS<Q, 6> (&input[35], &input[36]);
	CAS<Q, 6> (&input[37], &input[38]);
	CAS<Q, 6> (&input[39], &input[40]);
	CAS<Q, 6> (&input[41], &input[42]);
	CAS<Q, 6> (&input[43], &input[44]);
	CAS<Q, 6> (&input[45], &input[46]);
	CAS<Q, 6> (&input[47], &input[48]);
	CAS<Q, 6> (&input[49], &input[50]);
	CAS<Q, 6> (&input[51], &input[52]);
	CAS<Q, 6> (&input[53], &input[54]);
	CAS<Q, 6> (&input[55], &input[56]);
	CAS<Q, 6> (&input[57], &input[58]);
	CAS<Q, 6> (&input[59], &input[60]);
	CAS<Q, 6> (&input[61], &input[62]);
	CAS<Q, 6> (&input[63], &input[64]);
	CAS<Q, 6> (&input[65], &input[66]);
	CAS<Q, 6> (&input[67], &input[68]);
	CAS<Q, 6> (&input[69], &input[70]);
	CAS<Q, 6> (&input[71], &input[72]);
	CAS<Q, 6> (&input[73], &input[74]);
	CAS<Q, 6> (&input[75], &input[76]);
	CAS<Q, 6> (&input[77], &input[78]);
	CAS<Q, 6> (&input[79], &input[80]);
	CAS<Q, 6> (&input[81], &input[82]);
	CAS<Q, 6> (&input[83], &input[84]);
	CAS<Q, 6> (&input[85], &input[86]);
	CAS<Q, 6> (&input[87], &input[88]);
	CAS<Q, 6> (&input[89], &input[90]);
	CAS<Q, 6> (&input[81], &input[92]);
	CAS<Q, 6> (&input[93], &input[94]);
	CAS<Q, 6> (&input[95], &input[96]);
	//34
	CAS<Q, 6> (&input[34], &input[35]);
	CAS<Q, 6> (&input[36], &input[37]);
	CAS<Q, 6> (&input[38], &input[39]);
	CAS<Q, 6> (&input[40], &input[41]);
	CAS<Q, 6> (&input[42], &input[43]);
	CAS<Q, 6> (&input[44], &input[45]);
	CAS<Q, 6> (&input[46], &input[47]);
	CAS<Q, 6> (&input[48], &input[49]);
	CAS<Q, 6> (&input[50], &input[51]);
	CAS<Q, 6> (&input[52], &input[53]);
	CAS<Q, 6> (&input[54], &input[55]);
	CAS<Q, 6> (&input[56], &input[57]);
	CAS<Q, 6> (&input[58], &input[59]);
	CAS<Q, 6> (&input[60], &input[61]);
	CAS<Q, 6> (&input[62], &input[63]);
	CAS<Q, 6> (&input[64], &input[65]);
	CAS<Q, 6> (&input[66], &input[67]);
	CAS<Q, 6> (&input[68], &input[69]);
	CAS<Q, 6> (&input[70], &input[71]);
	CAS<Q, 6> (&input[72], &input[73]);
	CAS<Q, 6> (&input[74], &input[75]);
	CAS<Q, 6> (&input[76], &input[77]);
	CAS<Q, 6> (&input[78], &input[79]);
	CAS<Q, 6> (&input[80], &input[81]);
	CAS<Q, 6> (&input[82], &input[83]);
	CAS<Q, 6> (&input[84], &input[85]);
	CAS<Q, 6> (&input[86], &input[87]);
	CAS<Q, 6> (&input[88], &input[89]);
	CAS<Q, 6> (&input[90], &input[91]);
	CAS<Q, 6> (&input[92], &input[93]);
	CAS<Q, 6> (&input[94], &input[95]);
	//35
	CAS<Q, 6> (&input[33], &input[34]);
	CAS<Q, 6> (&input[35], &input[36]);
	CAS<Q, 6> (&input[37], &input[38]);
	CAS<Q, 6> (&input[39], &input[40]);
	CAS<Q, 6> (&input[41], &input[42]);
	CAS<Q, 6> (&input[43], &input[44]);
	CAS<Q, 6> (&input[45], &input[46]);
	CAS<Q, 6> (&input[47], &input[48]);
	CAS<Q, 6> (&input[49], &input[50]);
	CAS<Q, 6> (&input[51], &input[52]);
	CAS<Q, 6> (&input[53], &input[54]);
	CAS<Q, 6> (&input[55], &input[56]);
	CAS<Q, 6> (&input[57], &input[58]);
	CAS<Q, 6> (&input[59], &input[60]);
	CAS<Q, 6> (&input[61], &input[62]);
	CAS<Q, 6> (&input[63], &input[64]);
	CAS<Q, 6> (&input[65], &input[66]);
	CAS<Q, 6> (&input[67], &input[68]);
	CAS<Q, 6> (&input[69], &input[70]);
	CAS<Q, 6> (&input[71], &input[72]);
	CAS<Q, 6> (&input[73], &input[74]);
	CAS<Q, 6> (&input[75], &input[76]);
	CAS<Q, 6> (&input[77], &input[78]);
	CAS<Q, 6> (&input[79], &input[80]);
	CAS<Q, 6> (&input[81], &input[82]);
	CAS<Q, 6> (&input[83], &input[84]);
	CAS<Q, 6> (&input[85], &input[86]);
	CAS<Q, 6> (&input[87], &input[88]);
	CAS<Q, 6> (&input[89], &input[90]);
	CAS<Q, 6> (&input[81], &input[92]);
	CAS<Q, 6> (&input[93], &input[94]);
	CAS<Q, 6> (&input[95], &input[96]);
	//36
	CAS<Q, 6> (&input[34], &input[35]);
	CAS<Q, 6> (&input[36], &input[37]);
	CAS<Q, 6> (&input[38], &input[39]);
	CAS<Q, 6> (&input[40], &input[41]);
	CAS<Q, 6> (&input[42], &input[43]);
	CAS<Q, 6> (&input[44], &input[45]);
	CAS<Q, 6> (&input[46], &input[47]);
	CAS<Q, 6> (&input[48], &input[49]);
	CAS<Q, 6> (&input[50], &input[51]);
	CAS<Q, 6> (&input[52], &input[53]);
	CAS<Q, 6> (&input[54], &input[55]);
	CAS<Q, 6> (&input[56], &input[57]);
	CAS<Q, 6> (&input[58], &input[59]);
	CAS<Q, 6> (&input[60], &input[61]);
	CAS<Q, 6> (&input[62], &input[63]);
	CAS<Q, 6> (&input[64], &input[65]);
	CAS<Q, 6> (&input[66], &input[67]);
	CAS<Q, 6> (&input[68], &input[69]);
	CAS<Q, 6> (&input[70], &input[71]);
	CAS<Q, 6> (&input[72], &input[73]);
	CAS<Q, 6> (&input[74], &input[75]);
	CAS<Q, 6> (&input[76], &input[77]);
	CAS<Q, 6> (&input[78], &input[79]);
	CAS<Q, 6> (&input[80], &input[81]);
	CAS<Q, 6> (&input[82], &input[83]);
	CAS<Q, 6> (&input[84], &input[85]);
	CAS<Q, 6> (&input[86], &input[87]);
	CAS<Q, 6> (&input[88], &input[89]);
	CAS<Q, 6> (&input[90], &input[91]);
	CAS<Q, 6> (&input[92], &input[93]);
	CAS<Q, 6> (&input[94], &input[95]);
	//37
	CAS<Q, 6> (&input[33], &input[34]);
	CAS<Q, 6> (&input[35], &input[36]);
	CAS<Q, 6> (&input[37], &input[38]);
	CAS<Q, 6> (&input[39], &input[40]);
	CAS<Q, 6> (&input[41], &input[42]);
	CAS<Q, 6> (&input[43], &input[44]);
	CAS<Q, 6> (&input[45], &input[46]);
	CAS<Q, 6> (&input[47], &input[48]);
	CAS<Q, 6> (&input[49], &input[50]);
	CAS<Q, 6> (&input[51], &input[52]);
	CAS<Q, 6> (&input[53], &input[54]);
	CAS<Q, 6> (&input[55], &input[56]);
	CAS<Q, 6> (&input[57], &input[58]);
	CAS<Q, 6> (&input[59], &input[60]);
	CAS<Q, 6> (&input[61], &input[62]);
	CAS<Q, 6> (&input[63], &input[64]);
	CAS<Q, 6> (&input[65], &input[66]);
	CAS<Q, 6> (&input[67], &input[68]);
	CAS<Q, 6> (&input[69], &input[70]);
	CAS<Q, 6> (&input[71], &input[72]);
	CAS<Q, 6> (&input[73], &input[74]);
	CAS<Q, 6> (&input[75], &input[76]);
	CAS<Q, 6> (&input[77], &input[78]);
	CAS<Q, 6> (&input[79], &input[80]);
	CAS<Q, 6> (&input[81], &input[82]);
	CAS<Q, 6> (&input[83], &input[84]);
	CAS<Q, 6> (&input[85], &input[86]);
	CAS<Q, 6> (&input[87], &input[88]);
	CAS<Q, 6> (&input[89], &input[90]);
	CAS<Q, 6> (&input[81], &input[92]);
	CAS<Q, 6> (&input[93], &input[94]);
	CAS<Q, 6> (&input[95], &input[96]);
	//38
	CAS<Q, 6> (&input[34], &input[35]);
	CAS<Q, 6> (&input[36], &input[37]);
	CAS<Q, 6> (&input[38], &input[39]);
	CAS<Q, 6> (&input[40], &input[41]);
	CAS<Q, 6> (&input[42], &input[43]);
	CAS<Q, 6> (&input[44], &input[45]);
	CAS<Q, 6> (&input[46], &input[47]);
	CAS<Q, 6> (&input[48], &input[49]);
	CAS<Q, 6> (&input[50], &input[51]);
	CAS<Q, 6> (&input[52], &input[53]);
	CAS<Q, 6> (&input[54], &input[55]);
	CAS<Q, 6> (&input[56], &input[57]);
	CAS<Q, 6> (&input[58], &input[59]);
	CAS<Q, 6> (&input[60], &input[61]);
	CAS<Q, 6> (&input[62], &input[63]);
	CAS<Q, 6> (&input[64], &input[65]);
	CAS<Q, 6> (&input[66], &input[67]);
	CAS<Q, 6> (&input[68], &input[69]);
	CAS<Q, 6> (&input[70], &input[71]);
	CAS<Q, 6> (&input[72], &input[73]);
	CAS<Q, 6> (&input[74], &input[75]);
	CAS<Q, 6> (&input[76], &input[77]);
	CAS<Q, 6> (&input[78], &input[79]);
	CAS<Q, 6> (&input[80], &input[81]);
	CAS<Q, 6> (&input[82], &input[83]);
	CAS<Q, 6> (&input[84], &input[85]);
	CAS<Q, 6> (&input[86], &input[87]);
	CAS<Q, 6> (&input[88], &input[89]);
	CAS<Q, 6> (&input[90], &input[91]);
	CAS<Q, 6> (&input[92], &input[93]);
	CAS<Q, 6> (&input[94], &input[95]);
	//39
	CAS<Q, 6> (&input[33], &input[34]);
	CAS<Q, 6> (&input[35], &input[36]);
	CAS<Q, 6> (&input[37], &input[38]);
	CAS<Q, 6> (&input[39], &input[40]);
	CAS<Q, 6> (&input[41], &input[42]);
	CAS<Q, 6> (&input[43], &input[44]);
	CAS<Q, 6> (&input[45], &input[46]);
	CAS<Q, 6> (&input[47], &input[48]);
	CAS<Q, 6> (&input[49], &input[50]);
	CAS<Q, 6> (&input[51], &input[52]);
	CAS<Q, 6> (&input[53], &input[54]);
	CAS<Q, 6> (&input[55], &input[56]);
	CAS<Q, 6> (&input[57], &input[58]);
	CAS<Q, 6> (&input[59], &input[60]);
	CAS<Q, 6> (&input[61], &input[62]);
	CAS<Q, 6> (&input[63], &input[64]);
	CAS<Q, 6> (&input[65], &input[66]);
	CAS<Q, 6> (&input[67], &input[68]);
	CAS<Q, 6> (&input[69], &input[70]);
	CAS<Q, 6> (&input[71], &input[72]);
	CAS<Q, 6> (&input[73], &input[74]);
	CAS<Q, 6> (&input[75], &input[76]);
	CAS<Q, 6> (&input[77], &input[78]);
	CAS<Q, 6> (&input[79], &input[80]);
	CAS<Q, 6> (&input[81], &input[82]);
	CAS<Q, 6> (&input[83], &input[84]);
	CAS<Q, 6> (&input[85], &input[86]);
	CAS<Q, 6> (&input[87], &input[88]);
	CAS<Q, 6> (&input[89], &input[90]);
	CAS<Q, 6> (&input[81], &input[92]);
	CAS<Q, 6> (&input[93], &input[94]);
	CAS<Q, 6> (&input[95], &input[96]);
	//40
	CAS<Q, 6> (&input[34], &input[35]);
	CAS<Q, 6> (&input[36], &input[37]);
	CAS<Q, 6> (&input[38], &input[39]);
	CAS<Q, 6> (&input[40], &input[41]);
	CAS<Q, 6> (&input[42], &input[43]);
	CAS<Q, 6> (&input[44], &input[45]);
	CAS<Q, 6> (&input[46], &input[47]);
	CAS<Q, 6> (&input[48], &input[49]);
	CAS<Q, 6> (&input[50], &input[51]);
	CAS<Q, 6> (&input[52], &input[53]);
	CAS<Q, 6> (&input[54], &input[55]);
	CAS<Q, 6> (&input[56], &input[57]);
	CAS<Q, 6> (&input[58], &input[59]);
	CAS<Q, 6> (&input[60], &input[61]);
	CAS<Q, 6> (&input[62], &input[63]);
	CAS<Q, 6> (&input[64], &input[65]);
	CAS<Q, 6> (&input[66], &input[67]);
	CAS<Q, 6> (&input[68], &input[69]);
	CAS<Q, 6> (&input[70], &input[71]);
	CAS<Q, 6> (&input[72], &input[73]);
	CAS<Q, 6> (&input[74], &input[75]);
	CAS<Q, 6> (&input[76], &input[77]);
	CAS<Q, 6> (&input[78], &input[79]);
	CAS<Q, 6> (&input[80], &input[81]);
	CAS<Q, 6> (&input[82], &input[83]);
	CAS<Q, 6> (&input[84], &input[85]);
	CAS<Q, 6> (&input[86], &input[87]);
	CAS<Q, 6> (&input[88], &input[89]);
	CAS<Q, 6> (&input[90], &input[91]);
	CAS<Q, 6> (&input[92], &input[93]);
	CAS<Q, 6> (&input[94], &input[95]);
	//41
	CAS<Q, 6> (&input[33], &input[34]);
	CAS<Q, 6> (&input[35], &input[36]);
	CAS<Q, 6> (&input[37], &input[38]);
	CAS<Q, 6> (&input[39], &input[40]);
	CAS<Q, 6> (&input[41], &input[42]);
	CAS<Q, 6> (&input[43], &input[44]);
	CAS<Q, 6> (&input[45], &input[46]);
	CAS<Q, 6> (&input[47], &input[48]);
	CAS<Q, 6> (&input[49], &input[50]);
	CAS<Q, 6> (&input[51], &input[52]);
	CAS<Q, 6> (&input[53], &input[54]);
	CAS<Q, 6> (&input[55], &input[56]);
	CAS<Q, 6> (&input[57], &input[58]);
	CAS<Q, 6> (&input[59], &input[60]);
	CAS<Q, 6> (&input[61], &input[62]);
	CAS<Q, 6> (&input[63], &input[64]);
	CAS<Q, 6> (&input[65], &input[66]);
	CAS<Q, 6> (&input[67], &input[68]);
	CAS<Q, 6> (&input[69], &input[70]);
	CAS<Q, 6> (&input[71], &input[72]);
	CAS<Q, 6> (&input[73], &input[74]);
	CAS<Q, 6> (&input[75], &input[76]);
	CAS<Q, 6> (&input[77], &input[78]);
	CAS<Q, 6> (&input[79], &input[80]);
	CAS<Q, 6> (&input[81], &input[82]);
	CAS<Q, 6> (&input[83], &input[84]);
	CAS<Q, 6> (&input[85], &input[86]);
	CAS<Q, 6> (&input[87], &input[88]);
	CAS<Q, 6> (&input[89], &input[90]);
	CAS<Q, 6> (&input[81], &input[92]);
	CAS<Q, 6> (&input[93], &input[94]);
	CAS<Q, 6> (&input[95], &input[96]);
	//42
	CAS<Q, 6> (&input[34], &input[35]);
	CAS<Q, 6> (&input[36], &input[37]);
	CAS<Q, 6> (&input[38], &input[39]);
	CAS<Q, 6> (&input[40], &input[41]);
	CAS<Q, 6> (&input[42], &input[43]);
	CAS<Q, 6> (&input[44], &input[45]);
	CAS<Q, 6> (&input[46], &input[47]);
	CAS<Q, 6> (&input[48], &input[49]);
	CAS<Q, 6> (&input[50], &input[51]);
	CAS<Q, 6> (&input[52], &input[53]);
	CAS<Q, 6> (&input[54], &input[55]);
	CAS<Q, 6> (&input[56], &input[57]);
	CAS<Q, 6> (&input[58], &input[59]);
	CAS<Q, 6> (&input[60], &input[61]);
	CAS<Q, 6> (&input[62], &input[63]);
	CAS<Q, 6> (&input[64], &input[65]);
	CAS<Q, 6> (&input[66], &input[67]);
	CAS<Q, 6> (&input[68], &input[69]);
	CAS<Q, 6> (&input[70], &input[71]);
	CAS<Q, 6> (&input[72], &input[73]);
	CAS<Q, 6> (&input[74], &input[75]);
	CAS<Q, 6> (&input[76], &input[77]);
	CAS<Q, 6> (&input[78], &input[79]);
	CAS<Q, 6> (&input[80], &input[81]);
	CAS<Q, 6> (&input[82], &input[83]);
	CAS<Q, 6> (&input[84], &input[85]);
	CAS<Q, 6> (&input[86], &input[87]);
	CAS<Q, 6> (&input[88], &input[89]);
	CAS<Q, 6> (&input[90], &input[91]);
	CAS<Q, 6> (&input[92], &input[93]);
	CAS<Q, 6> (&input[94], &input[95]);
	//43
	CAS<Q, 6> (&input[33], &input[34]);
	CAS<Q, 6> (&input[35], &input[36]);
	CAS<Q, 6> (&input[37], &input[38]);
	CAS<Q, 6> (&input[39], &input[40]);
	CAS<Q, 6> (&input[41], &input[42]);
	CAS<Q, 6> (&input[43], &input[44]);
	CAS<Q, 6> (&input[45], &input[46]);
	CAS<Q, 6> (&input[47], &input[48]);
	CAS<Q, 6> (&input[49], &input[50]);
	CAS<Q, 6> (&input[51], &input[52]);
	CAS<Q, 6> (&input[53], &input[54]);
	CAS<Q, 6> (&input[55], &input[56]);
	CAS<Q, 6> (&input[57], &input[58]);
	CAS<Q, 6> (&input[59], &input[60]);
	CAS<Q, 6> (&input[61], &input[62]);
	CAS<Q, 6> (&input[63], &input[64]);
	CAS<Q, 6> (&input[65], &input[66]);
	CAS<Q, 6> (&input[67], &input[68]);
	CAS<Q, 6> (&input[69], &input[70]);
	CAS<Q, 6> (&input[71], &input[72]);
	CAS<Q, 6> (&input[73], &input[74]);
	CAS<Q, 6> (&input[75], &input[76]);
	CAS<Q, 6> (&input[77], &input[78]);
	CAS<Q, 6> (&input[79], &input[80]);
	CAS<Q, 6> (&input[81], &input[82]);
	CAS<Q, 6> (&input[83], &input[84]);
	CAS<Q, 6> (&input[85], &input[86]);
	CAS<Q, 6> (&input[87], &input[88]);
	CAS<Q, 6> (&input[89], &input[90]);
	CAS<Q, 6> (&input[81], &input[92]);
	CAS<Q, 6> (&input[93], &input[94]);
	CAS<Q, 6> (&input[95], &input[96]);
	//44
	CAS<Q, 6> (&input[34], &input[35]);
	CAS<Q, 6> (&input[36], &input[37]);
	CAS<Q, 6> (&input[38], &input[39]);
	CAS<Q, 6> (&input[40], &input[41]);
	CAS<Q, 6> (&input[42], &input[43]);
	CAS<Q, 6> (&input[44], &input[45]);
	CAS<Q, 6> (&input[46], &input[47]);
	CAS<Q, 6> (&input[48], &input[49]);
	CAS<Q, 6> (&input[50], &input[51]);
	CAS<Q, 6> (&input[52], &input[53]);
	CAS<Q, 6> (&input[54], &input[55]);
	CAS<Q, 6> (&input[56], &input[57]);
	CAS<Q, 6> (&input[58], &input[59]);
	CAS<Q, 6> (&input[60], &input[61]);
	CAS<Q, 6> (&input[62], &input[63]);
	CAS<Q, 6> (&input[64], &input[65]);
	CAS<Q, 6> (&input[66], &input[67]);
	CAS<Q, 6> (&input[68], &input[69]);
	CAS<Q, 6> (&input[70], &input[71]);
	CAS<Q, 6> (&input[72], &input[73]);
	CAS<Q, 6> (&input[74], &input[75]);
	CAS<Q, 6> (&input[76], &input[77]);
	CAS<Q, 6> (&input[78], &input[79]);
	CAS<Q, 6> (&input[80], &input[81]);
	CAS<Q, 6> (&input[82], &input[83]);
	CAS<Q, 6> (&input[84], &input[85]);
	CAS<Q, 6> (&input[86], &input[87]);
	CAS<Q, 6> (&input[88], &input[89]);
	CAS<Q, 6> (&input[90], &input[91]);
	CAS<Q, 6> (&input[92], &input[93]);
	CAS<Q, 6> (&input[94], &input[95]);
	//45
	CAS<Q, 6> (&input[33], &input[34]);
	CAS<Q, 6> (&input[35], &input[36]);
	CAS<Q, 6> (&input[37], &input[38]);
	CAS<Q, 6> (&input[39], &input[40]);
	CAS<Q, 6> (&input[41], &input[42]);
	CAS<Q, 6> (&input[43], &input[44]);
	CAS<Q, 6> (&input[45], &input[46]);
	CAS<Q, 6> (&input[47], &input[48]);
	CAS<Q, 6> (&input[49], &input[50]);
	CAS<Q, 6> (&input[51], &input[52]);
	CAS<Q, 6> (&input[53], &input[54]);
	CAS<Q, 6> (&input[55], &input[56]);
	CAS<Q, 6> (&input[57], &input[58]);
	CAS<Q, 6> (&input[59], &input[60]);
	CAS<Q, 6> (&input[61], &input[62]);
	CAS<Q, 6> (&input[63], &input[64]);
	CAS<Q, 6> (&input[65], &input[66]);
	CAS<Q, 6> (&input[67], &input[68]);
	CAS<Q, 6> (&input[69], &input[70]);
	CAS<Q, 6> (&input[71], &input[72]);
	CAS<Q, 6> (&input[73], &input[74]);
	CAS<Q, 6> (&input[75], &input[76]);
	CAS<Q, 6> (&input[77], &input[78]);
	CAS<Q, 6> (&input[79], &input[80]);
	CAS<Q, 6> (&input[81], &input[82]);
	CAS<Q, 6> (&input[83], &input[84]);
	CAS<Q, 6> (&input[85], &input[86]);
	CAS<Q, 6> (&input[87], &input[88]);
	CAS<Q, 6> (&input[89], &input[90]);
	CAS<Q, 6> (&input[81], &input[92]);
	CAS<Q, 6> (&input[93], &input[94]);
	CAS<Q, 6> (&input[95], &input[96]);
	//46
	CAS<Q, 6> (&input[34], &input[35]);
	CAS<Q, 6> (&input[36], &input[37]);
	CAS<Q, 6> (&input[38], &input[39]);
	CAS<Q, 6> (&input[40], &input[41]);
	CAS<Q, 6> (&input[42], &input[43]);
	CAS<Q, 6> (&input[44], &input[45]);
	CAS<Q, 6> (&input[46], &input[47]);
	CAS<Q, 6> (&input[48], &input[49]);
	CAS<Q, 6> (&input[50], &input[51]);
	CAS<Q, 6> (&input[52], &input[53]);
	CAS<Q, 6> (&input[54], &input[55]);
	CAS<Q, 6> (&input[56], &input[57]);
	CAS<Q, 6> (&input[58], &input[59]);
	CAS<Q, 6> (&input[60], &input[61]);
	CAS<Q, 6> (&input[62], &input[63]);
	CAS<Q, 6> (&input[64], &input[65]);
	CAS<Q, 6> (&input[66], &input[67]);
	CAS<Q, 6> (&input[68], &input[69]);
	CAS<Q, 6> (&input[70], &input[71]);
	CAS<Q, 6> (&input[72], &input[73]);
	CAS<Q, 6> (&input[74], &input[75]);
	CAS<Q, 6> (&input[76], &input[77]);
	CAS<Q, 6> (&input[78], &input[79]);
	CAS<Q, 6> (&input[80], &input[81]);
	CAS<Q, 6> (&input[82], &input[83]);
	CAS<Q, 6> (&input[84], &input[85]);
	CAS<Q, 6> (&input[86], &input[87]);
	CAS<Q, 6> (&input[88], &input[89]);
	CAS<Q, 6> (&input[90], &input[91]);
	CAS<Q, 6> (&input[92], &input[93]);
	CAS<Q, 6> (&input[94], &input[95]);
	//47
	CAS<Q, 6> (&input[33], &input[34]);
	CAS<Q, 6> (&input[35], &input[36]);
	CAS<Q, 6> (&input[37], &input[38]);
	CAS<Q, 6> (&input[39], &input[40]);
	CAS<Q, 6> (&input[41], &input[42]);
	CAS<Q, 6> (&input[43], &input[44]);
	CAS<Q, 6> (&input[45], &input[46]);
	CAS<Q, 6> (&input[47], &input[48]);
	CAS<Q, 6> (&input[49], &input[50]);
	CAS<Q, 6> (&input[51], &input[52]);
	CAS<Q, 6> (&input[53], &input[54]);
	CAS<Q, 6> (&input[55], &input[56]);
	CAS<Q, 6> (&input[57], &input[58]);
	CAS<Q, 6> (&input[59], &input[60]);
	CAS<Q, 6> (&input[61], &input[62]);
	CAS<Q, 6> (&input[63], &input[64]);
	CAS<Q, 6> (&input[65], &input[66]);
	CAS<Q, 6> (&input[67], &input[68]);
	CAS<Q, 6> (&input[69], &input[70]);
	CAS<Q, 6> (&input[71], &input[72]);
	CAS<Q, 6> (&input[73], &input[74]);
	CAS<Q, 6> (&input[75], &input[76]);
	CAS<Q, 6> (&input[77], &input[78]);
	CAS<Q, 6> (&input[79], &input[80]);
	CAS<Q, 6> (&input[81], &input[82]);
	CAS<Q, 6> (&input[83], &input[84]);
	CAS<Q, 6> (&input[85], &input[86]);
	CAS<Q, 6> (&input[87], &input[88]);
	CAS<Q, 6> (&input[89], &input[90]);
	CAS<Q, 6> (&input[81], &input[92]);
	CAS<Q, 6> (&input[93], &input[94]);
	CAS<Q, 6> (&input[95], &input[96]);
	//48
	CAS<Q, 6> (&input[34], &input[35]);
	CAS<Q, 6> (&input[36], &input[37]);
	CAS<Q, 6> (&input[38], &input[39]);
	CAS<Q, 6> (&input[40], &input[41]);
	CAS<Q, 6> (&input[42], &input[43]);
	CAS<Q, 6> (&input[44], &input[45]);
	CAS<Q, 6> (&input[46], &input[47]);
	CAS<Q, 6> (&input[48], &input[49]);
	CAS<Q, 6> (&input[50], &input[51]);
	CAS<Q, 6> (&input[52], &input[53]);
	CAS<Q, 6> (&input[54], &input[55]);
	CAS<Q, 6> (&input[56], &input[57]);
	CAS<Q, 6> (&input[58], &input[59]);
	CAS<Q, 6> (&input[60], &input[61]);
	CAS<Q, 6> (&input[62], &input[63]);
	CAS<Q, 6> (&input[64], &input[65]);
	CAS<Q, 6> (&input[66], &input[67]);
	CAS<Q, 6> (&input[68], &input[69]);
	CAS<Q, 6> (&input[70], &input[71]);
	CAS<Q, 6> (&input[72], &input[73]);
	CAS<Q, 6> (&input[74], &input[75]);
	CAS<Q, 6> (&input[76], &input[77]);
	CAS<Q, 6> (&input[78], &input[79]);
	CAS<Q, 6> (&input[80], &input[81]);
	CAS<Q, 6> (&input[82], &input[83]);
	CAS<Q, 6> (&input[84], &input[85]);
	CAS<Q, 6> (&input[86], &input[87]);
	CAS<Q, 6> (&input[88], &input[89]);
	CAS<Q, 6> (&input[90], &input[91]);
	CAS<Q, 6> (&input[92], &input[93]);
	CAS<Q, 6> (&input[94], &input[95]);
	//49
	CAS<Q, 6> (&input[33], &input[34]);
	CAS<Q, 6> (&input[35], &input[36]);
	CAS<Q, 6> (&input[37], &input[38]);
	CAS<Q, 6> (&input[39], &input[40]);
	CAS<Q, 6> (&input[41], &input[42]);
	CAS<Q, 6> (&input[43], &input[44]);
	CAS<Q, 6> (&input[45], &input[46]);
	CAS<Q, 6> (&input[47], &input[48]);
	CAS<Q, 6> (&input[49], &input[50]);
	CAS<Q, 6> (&input[51], &input[52]);
	CAS<Q, 6> (&input[53], &input[54]);
	CAS<Q, 6> (&input[55], &input[56]);
	CAS<Q, 6> (&input[57], &input[58]);
	CAS<Q, 6> (&input[59], &input[60]);
	CAS<Q, 6> (&input[61], &input[62]);
	CAS<Q, 6> (&input[63], &input[64]);
	CAS<Q, 6> (&input[65], &input[66]);
	CAS<Q, 6> (&input[67], &input[68]);
	CAS<Q, 6> (&input[69], &input[70]);
	CAS<Q, 6> (&input[71], &input[72]);
	CAS<Q, 6> (&input[73], &input[74]);
	CAS<Q, 6> (&input[75], &input[76]);
	CAS<Q, 6> (&input[77], &input[78]);
	CAS<Q, 6> (&input[79], &input[80]);
	CAS<Q, 6> (&input[81], &input[82]);
	CAS<Q, 6> (&input[83], &input[84]);
	CAS<Q, 6> (&input[85], &input[86]);
	CAS<Q, 6> (&input[87], &input[88]);
	CAS<Q, 6> (&input[89], &input[90]);
	CAS<Q, 6> (&input[81], &input[92]);
	CAS<Q, 6> (&input[93], &input[94]);
	CAS<Q, 6> (&input[95], &input[96]);
	//50
	CAS<Q, 6> (&input[34], &input[35]);
	CAS<Q, 6> (&input[36], &input[37]);
	CAS<Q, 6> (&input[38], &input[39]);
	CAS<Q, 6> (&input[40], &input[41]);
	CAS<Q, 6> (&input[42], &input[43]);
	CAS<Q, 6> (&input[44], &input[45]);
	CAS<Q, 6> (&input[46], &input[47]);
	CAS<Q, 6> (&input[48], &input[49]);
	CAS<Q, 6> (&input[50], &input[51]);
	CAS<Q, 6> (&input[52], &input[53]);
	CAS<Q, 6> (&input[54], &input[55]);
	CAS<Q, 6> (&input[56], &input[57]);
	CAS<Q, 6> (&input[58], &input[59]);
	CAS<Q, 6> (&input[60], &input[61]);
	CAS<Q, 6> (&input[62], &input[63]);
	CAS<Q, 6> (&input[64], &input[65]);
	CAS<Q, 6> (&input[66], &input[67]);
	CAS<Q, 6> (&input[68], &input[69]);
	CAS<Q, 6> (&input[70], &input[71]);
	CAS<Q, 6> (&input[72], &input[73]);
	CAS<Q, 6> (&input[74], &input[75]);
	CAS<Q, 6> (&input[76], &input[77]);
	CAS<Q, 6> (&input[78], &input[79]);
	CAS<Q, 6> (&input[80], &input[81]);
	CAS<Q, 6> (&input[82], &input[83]);
	CAS<Q, 6> (&input[84], &input[85]);
	CAS<Q, 6> (&input[86], &input[87]);
	CAS<Q, 6> (&input[88], &input[89]);
	CAS<Q, 6> (&input[90], &input[91]);
	CAS<Q, 6> (&input[92], &input[93]);
	CAS<Q, 6> (&input[94], &input[95]);
	//51
	CAS<Q, 6> (&input[33], &input[34]);
	CAS<Q, 6> (&input[35], &input[36]);
	CAS<Q, 6> (&input[37], &input[38]);
	CAS<Q, 6> (&input[39], &input[40]);
	CAS<Q, 6> (&input[41], &input[42]);
	CAS<Q, 6> (&input[43], &input[44]);
	CAS<Q, 6> (&input[45], &input[46]);
	CAS<Q, 6> (&input[47], &input[48]);
	CAS<Q, 6> (&input[49], &input[50]);
	CAS<Q, 6> (&input[51], &input[52]);
	CAS<Q, 6> (&input[53], &input[54]);
	CAS<Q, 6> (&input[55], &input[56]);
	CAS<Q, 6> (&input[57], &input[58]);
	CAS<Q, 6> (&input[59], &input[60]);
	CAS<Q, 6> (&input[61], &input[62]);
	CAS<Q, 6> (&input[63], &input[64]);
	CAS<Q, 6> (&input[65], &input[66]);
	CAS<Q, 6> (&input[67], &input[68]);
	CAS<Q, 6> (&input[69], &input[70]);
	CAS<Q, 6> (&input[71], &input[72]);
	CAS<Q, 6> (&input[73], &input[74]);
	CAS<Q, 6> (&input[75], &input[76]);
	CAS<Q, 6> (&input[77], &input[78]);
	CAS<Q, 6> (&input[79], &input[80]);
	CAS<Q, 6> (&input[81], &input[82]);
	CAS<Q, 6> (&input[83], &input[84]);
	CAS<Q, 6> (&input[85], &input[86]);
	CAS<Q, 6> (&input[87], &input[88]);
	CAS<Q, 6> (&input[89], &input[90]);
	CAS<Q, 6> (&input[81], &input[92]);
	CAS<Q, 6> (&input[93], &input[94]);
	CAS<Q, 6> (&input[95], &input[96]);
	//52
	CAS<Q, 6> (&input[34], &input[35]);
	CAS<Q, 6> (&input[36], &input[37]);
	CAS<Q, 6> (&input[38], &input[39]);
	CAS<Q, 6> (&input[40], &input[41]);
	CAS<Q, 6> (&input[42], &input[43]);
	CAS<Q, 6> (&input[44], &input[45]);
	CAS<Q, 6> (&input[46], &input[47]);
	CAS<Q, 6> (&input[48], &input[49]);
	CAS<Q, 6> (&input[50], &input[51]);
	CAS<Q, 6> (&input[52], &input[53]);
	CAS<Q, 6> (&input[54], &input[55]);
	CAS<Q, 6> (&input[56], &input[57]);
	CAS<Q, 6> (&input[58], &input[59]);
	CAS<Q, 6> (&input[60], &input[61]);
	CAS<Q, 6> (&input[62], &input[63]);
	CAS<Q, 6> (&input[64], &input[65]);
	CAS<Q, 6> (&input[66], &input[67]);
	CAS<Q, 6> (&input[68], &input[69]);
	CAS<Q, 6> (&input[70], &input[71]);
	CAS<Q, 6> (&input[72], &input[73]);
	CAS<Q, 6> (&input[74], &input[75]);
	CAS<Q, 6> (&input[76], &input[77]);
	CAS<Q, 6> (&input[78], &input[79]);
	CAS<Q, 6> (&input[80], &input[81]);
	CAS<Q, 6> (&input[82], &input[83]);
	CAS<Q, 6> (&input[84], &input[85]);
	CAS<Q, 6> (&input[86], &input[87]);
	CAS<Q, 6> (&input[88], &input[89]);
	CAS<Q, 6> (&input[90], &input[91]);
	CAS<Q, 6> (&input[92], &input[93]);
	CAS<Q, 6> (&input[94], &input[95]);
	//53
	CAS<Q, 6> (&input[33], &input[34]);
	CAS<Q, 6> (&input[35], &input[36]);
	CAS<Q, 6> (&input[37], &input[38]);
	CAS<Q, 6> (&input[39], &input[40]);
	CAS<Q, 6> (&input[41], &input[42]);
	CAS<Q, 6> (&input[43], &input[44]);
	CAS<Q, 6> (&input[45], &input[46]);
	CAS<Q, 6> (&input[47], &input[48]);
	CAS<Q, 6> (&input[49], &input[50]);
	CAS<Q, 6> (&input[51], &input[52]);
	CAS<Q, 6> (&input[53], &input[54]);
	CAS<Q, 6> (&input[55], &input[56]);
	CAS<Q, 6> (&input[57], &input[58]);
	CAS<Q, 6> (&input[59], &input[60]);
	CAS<Q, 6> (&input[61], &input[62]);
	CAS<Q, 6> (&input[63], &input[64]);
	CAS<Q, 6> (&input[65], &input[66]);
	CAS<Q, 6> (&input[67], &input[68]);
	CAS<Q, 6> (&input[69], &input[70]);
	CAS<Q, 6> (&input[71], &input[72]);
	CAS<Q, 6> (&input[73], &input[74]);
	CAS<Q, 6> (&input[75], &input[76]);
	CAS<Q, 6> (&input[77], &input[78]);
	CAS<Q, 6> (&input[79], &input[80]);
	CAS<Q, 6> (&input[81], &input[82]);
	CAS<Q, 6> (&input[83], &input[84]);
	CAS<Q, 6> (&input[85], &input[86]);
	CAS<Q, 6> (&input[87], &input[88]);
	CAS<Q, 6> (&input[89], &input[90]);
	CAS<Q, 6> (&input[81], &input[92]);
	CAS<Q, 6> (&input[93], &input[94]);
	CAS<Q, 6> (&input[95], &input[96]);
	//54
	CAS<Q, 6> (&input[34], &input[35]);
	CAS<Q, 6> (&input[36], &input[37]);
	CAS<Q, 6> (&input[38], &input[39]);
	CAS<Q, 6> (&input[40], &input[41]);
	CAS<Q, 6> (&input[42], &input[43]);
	CAS<Q, 6> (&input[44], &input[45]);
	CAS<Q, 6> (&input[46], &input[47]);
	CAS<Q, 6> (&input[48], &input[49]);
	CAS<Q, 6> (&input[50], &input[51]);
	CAS<Q, 6> (&input[52], &input[53]);
	CAS<Q, 6> (&input[54], &input[55]);
	CAS<Q, 6> (&input[56], &input[57]);
	CAS<Q, 6> (&input[58], &input[59]);
	CAS<Q, 6> (&input[60], &input[61]);
	CAS<Q, 6> (&input[62], &input[63]);
	CAS<Q, 6> (&input[64], &input[65]);
	CAS<Q, 6> (&input[66], &input[67]);
	CAS<Q, 6> (&input[68], &input[69]);
	CAS<Q, 6> (&input[70], &input[71]);
	CAS<Q, 6> (&input[72], &input[73]);
	CAS<Q, 6> (&input[74], &input[75]);
	CAS<Q, 6> (&input[76], &input[77]);
	CAS<Q, 6> (&input[78], &input[79]);
	CAS<Q, 6> (&input[80], &input[81]);
	CAS<Q, 6> (&input[82], &input[83]);
	CAS<Q, 6> (&input[84], &input[85]);
	CAS<Q, 6> (&input[86], &input[87]);
	CAS<Q, 6> (&input[88], &input[89]);
	CAS<Q, 6> (&input[90], &input[91]);
	CAS<Q, 6> (&input[92], &input[93]);
	CAS<Q, 6> (&input[94], &input[95]);
	//55
	CAS<Q, 6> (&input[33], &input[34]);
	CAS<Q, 6> (&input[35], &input[36]);
	CAS<Q, 6> (&input[37], &input[38]);
	CAS<Q, 6> (&input[39], &input[40]);
	CAS<Q, 6> (&input[41], &input[42]);
	CAS<Q, 6> (&input[43], &input[44]);
	CAS<Q, 6> (&input[45], &input[46]);
	CAS<Q, 6> (&input[47], &input[48]);
	CAS<Q, 6> (&input[49], &input[50]);
	CAS<Q, 6> (&input[51], &input[52]);
	CAS<Q, 6> (&input[53], &input[54]);
	CAS<Q, 6> (&input[55], &input[56]);
	CAS<Q, 6> (&input[57], &input[58]);
	CAS<Q, 6> (&input[59], &input[60]);
	CAS<Q, 6> (&input[61], &input[62]);
	CAS<Q, 6> (&input[63], &input[64]);
	CAS<Q, 6> (&input[65], &input[66]);
	CAS<Q, 6> (&input[67], &input[68]);
	CAS<Q, 6> (&input[69], &input[70]);
	CAS<Q, 6> (&input[71], &input[72]);
	CAS<Q, 6> (&input[73], &input[74]);
	CAS<Q, 6> (&input[75], &input[76]);
	CAS<Q, 6> (&input[77], &input[78]);
	CAS<Q, 6> (&input[79], &input[80]);
	CAS<Q, 6> (&input[81], &input[82]);
	CAS<Q, 6> (&input[83], &input[84]);
	CAS<Q, 6> (&input[85], &input[86]);
	CAS<Q, 6> (&input[87], &input[88]);
	CAS<Q, 6> (&input[89], &input[90]);
	CAS<Q, 6> (&input[81], &input[92]);
	CAS<Q, 6> (&input[93], &input[94]);
	CAS<Q, 6> (&input[95], &input[96]);
	//56
	CAS<Q, 6> (&input[34], &input[35]);
	CAS<Q, 6> (&input[36], &input[37]);
	CAS<Q, 6> (&input[38], &input[39]);
	CAS<Q, 6> (&input[40], &input[41]);
	CAS<Q, 6> (&input[42], &input[43]);
	CAS<Q, 6> (&input[44], &input[45]);
	CAS<Q, 6> (&input[46], &input[47]);
	CAS<Q, 6> (&input[48], &input[49]);
	CAS<Q, 6> (&input[50], &input[51]);
	CAS<Q, 6> (&input[52], &input[53]);
	CAS<Q, 6> (&input[54], &input[55]);
	CAS<Q, 6> (&input[56], &input[57]);
	CAS<Q, 6> (&input[58], &input[59]);
	CAS<Q, 6> (&input[60], &input[61]);
	CAS<Q, 6> (&input[62], &input[63]);
	CAS<Q, 6> (&input[64], &input[65]);
	CAS<Q, 6> (&input[66], &input[67]);
	CAS<Q, 6> (&input[68], &input[69]);
	CAS<Q, 6> (&input[70], &input[71]);
	CAS<Q, 6> (&input[72], &input[73]);
	CAS<Q, 6> (&input[74], &input[75]);
	CAS<Q, 6> (&input[76], &input[77]);
	CAS<Q, 6> (&input[78], &input[79]);
	CAS<Q, 6> (&input[80], &input[81]);
	CAS<Q, 6> (&input[82], &input[83]);
	CAS<Q, 6> (&input[84], &input[85]);
	CAS<Q, 6> (&input[86], &input[87]);
	CAS<Q, 6> (&input[88], &input[89]);
	CAS<Q, 6> (&input[90], &input[91]);
	CAS<Q, 6> (&input[92], &input[93]);
	CAS<Q, 6> (&input[94], &input[95]);
	//57
	CAS<Q, 6> (&input[33], &input[34]);
	CAS<Q, 6> (&input[35], &input[36]);
	CAS<Q, 6> (&input[37], &input[38]);
	CAS<Q, 6> (&input[39], &input[40]);
	CAS<Q, 6> (&input[41], &input[42]);
	CAS<Q, 6> (&input[43], &input[44]);
	CAS<Q, 6> (&input[45], &input[46]);
	CAS<Q, 6> (&input[47], &input[48]);
	CAS<Q, 6> (&input[49], &input[50]);
	CAS<Q, 6> (&input[51], &input[52]);
	CAS<Q, 6> (&input[53], &input[54]);
	CAS<Q, 6> (&input[55], &input[56]);
	CAS<Q, 6> (&input[57], &input[58]);
	CAS<Q, 6> (&input[59], &input[60]);
	CAS<Q, 6> (&input[61], &input[62]);
	CAS<Q, 6> (&input[63], &input[64]);
	CAS<Q, 6> (&input[65], &input[66]);
	CAS<Q, 6> (&input[67], &input[68]);
	CAS<Q, 6> (&input[69], &input[70]);
	CAS<Q, 6> (&input[71], &input[72]);
	CAS<Q, 6> (&input[73], &input[74]);
	CAS<Q, 6> (&input[75], &input[76]);
	CAS<Q, 6> (&input[77], &input[78]);
	CAS<Q, 6> (&input[79], &input[80]);
	CAS<Q, 6> (&input[81], &input[82]);
	CAS<Q, 6> (&input[83], &input[84]);
	CAS<Q, 6> (&input[85], &input[86]);
	CAS<Q, 6> (&input[87], &input[88]);
	CAS<Q, 6> (&input[89], &input[90]);
	CAS<Q, 6> (&input[81], &input[92]);
	CAS<Q, 6> (&input[93], &input[94]);
	CAS<Q, 6> (&input[95], &input[96]);
	//58
	CAS<Q, 6> (&input[34], &input[35]);
	CAS<Q, 6> (&input[36], &input[37]);
	CAS<Q, 6> (&input[38], &input[39]);
	CAS<Q, 6> (&input[40], &input[41]);
	CAS<Q, 6> (&input[42], &input[43]);
	CAS<Q, 6> (&input[44], &input[45]);
	CAS<Q, 6> (&input[46], &input[47]);
	CAS<Q, 6> (&input[48], &input[49]);
	CAS<Q, 6> (&input[50], &input[51]);
	CAS<Q, 6> (&input[52], &input[53]);
	CAS<Q, 6> (&input[54], &input[55]);
	CAS<Q, 6> (&input[56], &input[57]);
	CAS<Q, 6> (&input[58], &input[59]);
	CAS<Q, 6> (&input[60], &input[61]);
	CAS<Q, 6> (&input[62], &input[63]);
	CAS<Q, 6> (&input[64], &input[65]);
	CAS<Q, 6> (&input[66], &input[67]);
	CAS<Q, 6> (&input[68], &input[69]);
	CAS<Q, 6> (&input[70], &input[71]);
	CAS<Q, 6> (&input[72], &input[73]);
	CAS<Q, 6> (&input[74], &input[75]);
	CAS<Q, 6> (&input[76], &input[77]);
	CAS<Q, 6> (&input[78], &input[79]);
	CAS<Q, 6> (&input[80], &input[81]);
	CAS<Q, 6> (&input[82], &input[83]);
	CAS<Q, 6> (&input[84], &input[85]);
	CAS<Q, 6> (&input[86], &input[87]);
	CAS<Q, 6> (&input[88], &input[89]);
	CAS<Q, 6> (&input[90], &input[91]);
	CAS<Q, 6> (&input[92], &input[93]);
	CAS<Q, 6> (&input[94], &input[95]);
	//59
	CAS<Q, 6> (&input[33], &input[34]);
	CAS<Q, 6> (&input[35], &input[36]);
	CAS<Q, 6> (&input[37], &input[38]);
	CAS<Q, 6> (&input[39], &input[40]);
	CAS<Q, 6> (&input[41], &input[42]);
	CAS<Q, 6> (&input[43], &input[44]);
	CAS<Q, 6> (&input[45], &input[46]);
	CAS<Q, 6> (&input[47], &input[48]);
	CAS<Q, 6> (&input[49], &input[50]);
	CAS<Q, 6> (&input[51], &input[52]);
	CAS<Q, 6> (&input[53], &input[54]);
	CAS<Q, 6> (&input[55], &input[56]);
	CAS<Q, 6> (&input[57], &input[58]);
	CAS<Q, 6> (&input[59], &input[60]);
	CAS<Q, 6> (&input[61], &input[62]);
	CAS<Q, 6> (&input[63], &input[64]);
	CAS<Q, 6> (&input[65], &input[66]);
	CAS<Q, 6> (&input[67], &input[68]);
	CAS<Q, 6> (&input[69], &input[70]);
	CAS<Q, 6> (&input[71], &input[72]);
	CAS<Q, 6> (&input[73], &input[74]);
	CAS<Q, 6> (&input[75], &input[76]);
	CAS<Q, 6> (&input[77], &input[78]);
	CAS<Q, 6> (&input[79], &input[80]);
	CAS<Q, 6> (&input[81], &input[82]);
	CAS<Q, 6> (&input[83], &input[84]);
	CAS<Q, 6> (&input[85], &input[86]);
	CAS<Q, 6> (&input[87], &input[88]);
	CAS<Q, 6> (&input[89], &input[90]);
	CAS<Q, 6> (&input[81], &input[92]);
	CAS<Q, 6> (&input[93], &input[94]);
	CAS<Q, 6> (&input[95], &input[96]);
	//60
	CAS<Q, 6> (&input[34], &input[35]);
	CAS<Q, 6> (&input[36], &input[37]);
	CAS<Q, 6> (&input[38], &input[39]);
	CAS<Q, 6> (&input[40], &input[41]);
	CAS<Q, 6> (&input[42], &input[43]);
	CAS<Q, 6> (&input[44], &input[45]);
	CAS<Q, 6> (&input[46], &input[47]);
	CAS<Q, 6> (&input[48], &input[49]);
	CAS<Q, 6> (&input[50], &input[51]);
	CAS<Q, 6> (&input[52], &input[53]);
	CAS<Q, 6> (&input[54], &input[55]);
	CAS<Q, 6> (&input[56], &input[57]);
	CAS<Q, 6> (&input[58], &input[59]);
	CAS<Q, 6> (&input[60], &input[61]);
	CAS<Q, 6> (&input[62], &input[63]);
	CAS<Q, 6> (&input[64], &input[65]);
	CAS<Q, 6> (&input[66], &input[67]);
	CAS<Q, 6> (&input[68], &input[69]);
	CAS<Q, 6> (&input[70], &input[71]);
	CAS<Q, 6> (&input[72], &input[73]);
	CAS<Q, 6> (&input[74], &input[75]);
	CAS<Q, 6> (&input[76], &input[77]);
	CAS<Q, 6> (&input[78], &input[79]);
	CAS<Q, 6> (&input[80], &input[81]);
	CAS<Q, 6> (&input[82], &input[83]);
	CAS<Q, 6> (&input[84], &input[85]);
	CAS<Q, 6> (&input[86], &input[87]);
	CAS<Q, 6> (&input[88], &input[89]);
	CAS<Q, 6> (&input[90], &input[91]);
	CAS<Q, 6> (&input[92], &input[93]);
	CAS<Q, 6> (&input[94], &input[95]);
	//61
	CAS<Q, 6> (&input[33], &input[34]);
	CAS<Q, 6> (&input[35], &input[36]);
	CAS<Q, 6> (&input[37], &input[38]);
	CAS<Q, 6> (&input[39], &input[40]);
	CAS<Q, 6> (&input[41], &input[42]);
	CAS<Q, 6> (&input[43], &input[44]);
	CAS<Q, 6> (&input[45], &input[46]);
	CAS<Q, 6> (&input[47], &input[48]);
	CAS<Q, 6> (&input[49], &input[50]);
	CAS<Q, 6> (&input[51], &input[52]);
	CAS<Q, 6> (&input[53], &input[54]);
	CAS<Q, 6> (&input[55], &input[56]);
	CAS<Q, 6> (&input[57], &input[58]);
	CAS<Q, 6> (&input[59], &input[60]);
	CAS<Q, 6> (&input[61], &input[62]);
	CAS<Q, 6> (&input[63], &input[64]);
	CAS<Q, 6> (&input[65], &input[66]);
	CAS<Q, 6> (&input[67], &input[68]);
	CAS<Q, 6> (&input[69], &input[70]);
	CAS<Q, 6> (&input[71], &input[72]);
	CAS<Q, 6> (&input[73], &input[74]);
	CAS<Q, 6> (&input[75], &input[76]);
	CAS<Q, 6> (&input[77], &input[78]);
	CAS<Q, 6> (&input[79], &input[80]);
	CAS<Q, 6> (&input[81], &input[82]);
	CAS<Q, 6> (&input[83], &input[84]);
	CAS<Q, 6> (&input[85], &input[86]);
	CAS<Q, 6> (&input[87], &input[88]);
	CAS<Q, 6> (&input[89], &input[90]);
	CAS<Q, 6> (&input[81], &input[92]);
	CAS<Q, 6> (&input[93], &input[94]);
	CAS<Q, 6> (&input[95], &input[96]);
	//62
	CAS<Q, 6> (&input[34], &input[35]);
	CAS<Q, 6> (&input[36], &input[37]);
	CAS<Q, 6> (&input[38], &input[39]);
	CAS<Q, 6> (&input[40], &input[41]);
	CAS<Q, 6> (&input[42], &input[43]);
	CAS<Q, 6> (&input[44], &input[45]);
	CAS<Q, 6> (&input[46], &input[47]);
	CAS<Q, 6> (&input[48], &input[49]);
	CAS<Q, 6> (&input[50], &input[51]);
	CAS<Q, 6> (&input[52], &input[53]);
	CAS<Q, 6> (&input[54], &input[55]);
	CAS<Q, 6> (&input[56], &input[57]);
	CAS<Q, 6> (&input[58], &input[59]);
	CAS<Q, 6> (&input[60], &input[61]);
	CAS<Q, 6> (&input[62], &input[63]);
	CAS<Q, 6> (&input[64], &input[65]);
	CAS<Q, 6> (&input[66], &input[67]);
	CAS<Q, 6> (&input[68], &input[69]);
	CAS<Q, 6> (&input[70], &input[71]);
	CAS<Q, 6> (&input[72], &input[73]);
	CAS<Q, 6> (&input[74], &input[75]);
	CAS<Q, 6> (&input[76], &input[77]);
	CAS<Q, 6> (&input[78], &input[79]);
	CAS<Q, 6> (&input[80], &input[81]);
	CAS<Q, 6> (&input[82], &input[83]);
	CAS<Q, 6> (&input[84], &input[85]);
	CAS<Q, 6> (&input[86], &input[87]);
	CAS<Q, 6> (&input[88], &input[89]);
	CAS<Q, 6> (&input[90], &input[91]);
	CAS<Q, 6> (&input[92], &input[93]);
	CAS<Q, 6> (&input[94], &input[95]);
	//63
	CAS<Q, 6> (&input[33], &input[34]);
	CAS<Q, 6> (&input[35], &input[36]);
	CAS<Q, 6> (&input[37], &input[38]);
	CAS<Q, 6> (&input[39], &input[40]);
	CAS<Q, 6> (&input[41], &input[42]);
	CAS<Q, 6> (&input[43], &input[44]);
	CAS<Q, 6> (&input[45], &input[46]);
	CAS<Q, 6> (&input[47], &input[48]);
	CAS<Q, 6> (&input[49], &input[50]);
	CAS<Q, 6> (&input[51], &input[52]);
	CAS<Q, 6> (&input[53], &input[54]);
	CAS<Q, 6> (&input[55], &input[56]);
	CAS<Q, 6> (&input[57], &input[58]);
	CAS<Q, 6> (&input[59], &input[60]);
	CAS<Q, 6> (&input[61], &input[62]);
	CAS<Q, 6> (&input[63], &input[64]);
	CAS<Q, 6> (&input[65], &input[66]);
	CAS<Q, 6> (&input[67], &input[68]);
	CAS<Q, 6> (&input[69], &input[70]);
	CAS<Q, 6> (&input[71], &input[72]);
	CAS<Q, 6> (&input[73], &input[74]);
	CAS<Q, 6> (&input[75], &input[76]);
	CAS<Q, 6> (&input[77], &input[78]);
	CAS<Q, 6> (&input[79], &input[80]);
	CAS<Q, 6> (&input[81], &input[82]);
	CAS<Q, 6> (&input[83], &input[84]);
	CAS<Q, 6> (&input[85], &input[86]);
	CAS<Q, 6> (&input[87], &input[88]);
	CAS<Q, 6> (&input[89], &input[90]);
	CAS<Q, 6> (&input[81], &input[92]);
	CAS<Q, 6> (&input[93], &input[94]);
	CAS<Q, 6> (&input[95], &input[96]);
	//64
	CAS<Q, 6> (&input[34], &input[35]);
	CAS<Q, 6> (&input[36], &input[37]);
	CAS<Q, 6> (&input[38], &input[39]);
	CAS<Q, 6> (&input[40], &input[41]);
	CAS<Q, 6> (&input[42], &input[43]);
	CAS<Q, 6> (&input[44], &input[45]);
	CAS<Q, 6> (&input[46], &input[47]);
	CAS<Q, 6> (&input[48], &input[49]);
	CAS<Q, 6> (&input[50], &input[51]);
	CAS<Q, 6> (&input[52], &input[53]);
	CAS<Q, 6> (&input[54], &input[55]);
	CAS<Q, 6> (&input[56], &input[57]);
	CAS<Q, 6> (&input[58], &input[59]);
	CAS<Q, 6> (&input[60], &input[61]);
	CAS<Q, 6> (&input[62], &input[63]);
	CAS<Q, 6> (&input[64], &input[65]);
	CAS<Q, 6> (&input[66], &input[67]);
	CAS<Q, 6> (&input[68], &input[69]);
	CAS<Q, 6> (&input[70], &input[71]);
	CAS<Q, 6> (&input[72], &input[73]);
	CAS<Q, 6> (&input[74], &input[75]);
	CAS<Q, 6> (&input[76], &input[77]);
	CAS<Q, 6> (&input[78], &input[79]);
	CAS<Q, 6> (&input[80], &input[81]);
	CAS<Q, 6> (&input[82], &input[83]);
	CAS<Q, 6> (&input[84], &input[85]);
	CAS<Q, 6> (&input[86], &input[87]);
	CAS<Q, 6> (&input[88], &input[89]);
	CAS<Q, 6> (&input[90], &input[91]);
	CAS<Q, 6> (&input[92], &input[93]);
	CAS<Q, 6> (&input[94], &input[95]);
	
	/////////////////////////////////////////////////

   // output equal the L smallest input
	if (fb == 0){
		output[0] = input[33];
		output[1] = input[34];
		output[2] = input[35];
		output[3] = input[36];
		output[4] = input[37];
		output[5] = input[38];
		output[6] = input[39];
		output[7] = input[40];
		output[8] = input[41];
		output[9] = input[42];
		output[10] = input[43];
		output[11] = input[44];
		output[12] = input[45];
		output[13] = input[46];
		output[14] = input[47];
		output[15] = input[48];
		output[16] = input[49];
		output[17] = input[50];
		output[18] = input[51];
		output[19] = input[52];
		output[20] = input[53];
		output[21] = input[54];
		output[22] = input[55];
		output[23] = input[56];
		output[24] = input[57];
		output[25] = input[58];
		output[26] = input[59];
		output[27] = input[60];
		output[28] = input[61];
		output[29] = input[62];
		output[30] = input[63];
		output[31] = input[64];
		output[32] = input[65];
		output[33] = input[66];
		output[34] = input[67];
		output[35] = input[68];
		output[36] = input[69];
		output[37] = input[70];
		output[38] = input[71];
		output[39] = input[72];
		output[40] = input[73];
		output[41] = input[74];
		output[42] = input[75];
		output[43] = input[76];
		output[44] = input[77];
		output[45] = input[78];
		output[46] = input[79];
		output[47] = input[80];
		output[48] = input[81];
		output[49] = input[82];
		output[50] = input[83];
		output[51] = input[84];
		output[52] = input[85];
		output[53] = input[86];
		output[54] = input[87];
		output[55] = input[88];
		output[56] = input[89];
		output[57] = input[90];
		output[58] = input[91];
		output[59] = input[92];
		output[60] = input[93];
		output[61] = input[94];
		output[62] = input[95];
		output[63] = input[96];
	}
	else{
		output[0] = input[0];
		output[1] = input[1];
		output[2] = input[2];
		output[3] = input[3];
		output[4] = input[4];
		output[5] = input[5];
		output[6] = input[6];
		output[7] = input[7];
		output[8] = input[8];
		output[9] = input[9];
		output[10] = input[10];
		output[11] = input[11];
		output[12] = input[12];
		output[13] = input[13];
		output[14] = input[14];
		output[15] = input[15];
		output[16] = input[16];
		output[17] = input[17];
		output[18] = input[18];
		output[19] = input[19];
		output[20] = input[20];
		output[21] = input[21];
		output[22] = input[22];
		output[23] = input[23];
		output[24] = input[24];
		output[25] = input[25];
		output[26] = input[26];
		output[27] = input[27];
		output[28] = input[28];
		output[29] = input[29];
		output[30] = input[30];
		output[31] = input[31];
		output[32] = input[32];
		output[33] = input[33];
		output[34] = input[34];
		output[35] = input[35];
		output[36] = input[36];
		output[37] = input[37];
		output[38] = input[38];
		output[39] = input[39];
		output[40] = input[40];
		output[41] = input[41];
		output[42] = input[42];
		output[43] = input[43];
		output[44] = input[44];
		output[45] = input[45];
		output[46] = input[46];
		output[47] = input[47];
		output[48] = input[48];
		output[49] = input[49];
		output[50] = input[50];
		output[51] = input[51];
		output[52] = input[52];
		output[53] = input[53];
		output[54] = input[54];
		output[55] = input[55];
		output[56] = input[56];
		output[57] = input[57];
		output[58] = input[58];
		output[59] = input[59];
		output[60] = input[60];
		output[61] = input[61];
		output[62] = input[62];
		output[63] = input[63];
	}
}


template < int L, int Q, int LOG2L>
void CUSTOM_SORT (PS_struct<1,Q,LOG2L> input[2*L], PS_struct<1,Q,LOG2L> output[L], sc_biguint<1> fb)
{
#pragma HLS INLINE
#if L_SIZE == 2
	CUSTOM_SORT_L2 < Q > ( input, output, fb);
#elif L_SIZE == 4
	CUSTOM_SORT_L4 < Q > ( input, output, fb);
#elif L_SIZE == 8
	CUSTOM_SORT_L8 < Q > ( input, output, fb);
#elif L_SIZE == 16
	CUSTOM_SORT_L16 < Q > ( input, output, fb);
#elif L_SIZE == 32
	CUSTOM_SORT_L32 < Q > ( input, output, fb);
#elif L_SIZE == 64
	CUSTOM_SORT_L64 < Q > ( input, output, fb);
#endif

}

//*************************************************************************//
//** 						RANK ORDER SORTER							 **//
//*************************************************************************//

/* template < int L, int Q, int LOG2L>
void RANKORDER_SORT (PS_struct<1,Q,LOG2L> input[2*L], PS_struct<1,Q,LOG2L> output[L], sc_biguint<1> fb)
{
#pragma HLS INLINE off

// Adjust input when FB=0

	if (fb == 0){
		input[2*L-2] = input[(L-1)];
		input[2*L-1].metric = MAX_VAL;
		loop0:for(sc_int<LOG2L+1> k = L-2; k >= 0; k--)
		{
		#pragma HLS UNROLL
			input[2*k+1] = input[k];
		}
		loop1:for(sc_uint<LOG2L+1> k = 0; k < (L-1); k++)
		{
		#pragma HLS UNROLL
			input[2*k].metric = 0;

		}
	}

	PS_struct<1,Q,LOG2L> temp[2*L-1];
	sc_uint<1> comp_r [2*L_2*L];
	sc_uint<LOG2L + 1 > position [2*L];

	loop2:for(sc_uint<LOG2L+2> i = 0; i < 2*L; i++)
	{
	#pragma HLS UNROLL
		position[i] = (sc_uint<LOG2L + 1 >) 0;

	}
#pragma HLS ARRAY_PARTITION variable=temp complete dim=1
#pragma HLS ARRAY_PARTITION variable=comp_r complete dim=1
#pragma HLS ARRAY_PARTITION variable=position complete dim=1

// COMPARE inputs f

	// comparison odd value
	loop3:for(sc_uint<LOG2L+1> i = 1; i <= (2*L-3); i+=2)
	{
	#pragma HLS UNROLL
		loop31:for(sc_uint<LOG2L+1> j = i+1 ; j <= (2*L-2); j++)
			comp_r_i_j] = COMP <Q, LOG2L> (&input[i], &input[j]);
	}

	// comparison even value
	loop4:for(sc_uint<LOG2L+1> i = 0; i <= (2*L-4); i+=2)
	{
	#pragma HLS UNROLL
		loop41:for(sc_uint<LOG2L+1> j = i+1 ; j <= (2*L-2); j++){
			comp_r_i_j] = 0;
		}
	}

// Compute positions
	position[0] = 0;
#ifdef _MONITORING_
	position[2*L-1] = 2*L-1;
#endif
	loop5:for(sc_uint<LOG2L+1> k = 1; k <= (2*L-2); k++)
	{
	#pragma HLS UNROLL

		loop51:for(sc_uint<LOG2L+3> i = 0; i <= (k-1); i ++)
			position[k] += (sc_uint<1>) (~(comp_r_i_k]));

		loop52:for(sc_uint<LOG2L+3> j = k+1; j <= (2*L-2); j ++)
			position[k] += comp_r_k_j];
	}

// Multiplex the L smallest outputs
	temp[0] = input[0];
	loop6:for (sc_uint<LOG2L+1> k = 1; k < 2*L-1; k++){
	#pragma HLS UNROLL

	   int indice;
	   loop61:for(sc_uint<LOG2L+1> odd = 1; odd <= k; odd+= 2 )  // can be from all the odd input above
	   {
		   if (position[odd] == k )
			   indice = odd;
	   }
	   int a = (2*k);
	   int b = (2*L-2);
	   sc_uint<LOG2L+1> uper_bond = ( a > b ) ? b : a;
	   loop62:for(sc_uint<LOG2L+1> ptr = k; ptr <= uper_bond; ptr++ )  // or from all input under until 2*k
	   {
		   if (position[ptr] == k )
			   indice = ptr;
	   }
#ifndef _MONITORING_  // For RTL simulation MUX must be explicitly declared
	   temp[k] = input[indice];
#else

	   temp[k] = MUX_LIST_t < PS_struct<1,Q,LOG2L>, L, Q, LOG2L > (input, indice);
#endif
	}

// output the L smallest input
	if (fb == 0){
		loop7:for(sc_uint<LOG2L+1> i = 0; i < L; i++)
		{
		#pragma HLS UNROLL
			output[i] = temp[i+(L-1)];
		}
	}
	else{
		loop8:for(sc_uint<LOG2L+1> i = 0; i < L; i++)
		{
		#pragma HLS UNROLL
			output[i] = temp[i];
		}
	}

}*/

template <int Q>
void RANKORDER_SORT_L2 (PS_struct<1,Q,1> input[4], PS_struct<1,Q,1> output[2], sc_biguint<1> fb)
{
#pragma HLS INLINE off
#pragma HLS ARRAY_PARTITION variable=input complete dim=1
#pragma HLS ARRAY_PARTITION variable=output complete dim=1

// Adjust input when FB=0

	if (fb == 0){
		input[2] = input[1];
		input[1] = input[0];
	}

	sc_uint<2> position[4];
#pragma HLS ARRAY_PARTITION variable=position complete dim=1

// COMPARE inputs

	sc_uint<1> comp_r_1_2 = COMP < Q > (input[1].metric, input[2].metric);

// Compute positions

	position[0] = 0;
	position[1] = 1 + comp_r_1_2;
	position[2] = 1 + (sc_uint<1>) ~(comp_r_1_2);
	position[3] = (sc_uint<2>) 3;

// Multiplex the inputs in the order
	PS_struct<1,Q,1> temp[3];
#pragma HLS ARRAY_PARTITION variable=temp complete dim=1
	temp[0] = input[0];
	temp[1] = RO_MUX4 < PS_struct<1,Q,1> > (input, position, (sc_uint<2>) 1);
	temp[2] = RO_MUX4 < PS_struct<1,Q,1> > (input, position, (sc_uint<2>) 2);

// output the L smallest input
	if (fb == 0){
		output[0] = temp[1];
		output[1] = temp[2];
	}
	else{
		output[0] = temp[0];
		output[1] = temp[1];
	}

}

template <int Q>
void RANKORDER_SORT_L4 (PS_struct<1,Q,2> input[8], PS_struct<1,Q,2> output[4], sc_biguint<1> fb)
{
#pragma HLS INLINE
#pragma HLS ARRAY_PARTITION variable=input complete dim=1
#pragma HLS ARRAY_PARTITION variable=output complete dim=1

// Adjust input when FB=0

	if (fb == 0){
		input[7].metric = MAX_VAL;
		input[6] = input[3];
		input[5] = input[2];
		input[4].metric = 0;
		input[3] = input[1];
		input[2].metric = 0;
		input[1] = input[0];
	}

	sc_uint<3> position[8];
#pragma HLS ARRAY_PARTITION variable=position complete dim=1

// COMPARE inputs

	sc_uint<1> comp_r_1_2 = COMP < Q > (input[1].metric, input[2].metric);
	sc_uint<1> comp_r_1_3 = COMP < Q > (input[1].metric, input[3].metric);
	sc_uint<1> comp_r_1_4 = COMP < Q > (input[1].metric, input[4].metric);
	sc_uint<1> comp_r_1_5 = COMP < Q > (input[1].metric, input[5].metric);
	sc_uint<1> comp_r_1_6 = COMP < Q > (input[1].metric, input[6].metric);

	sc_uint<1> comp_r_3_4 = COMP < Q > (input[3].metric, input[4].metric);
	sc_uint<1> comp_r_3_5 = COMP < Q > (input[3].metric, input[5].metric);
	sc_uint<1> comp_r_3_6 = COMP < Q > (input[3].metric, input[6].metric);
	
	sc_uint<1> comp_r_5_6 = COMP < Q > (input[5].metric, input[6].metric);
	
// Compute positions
 
	position[0] = 0;
	position[1] = 1 + comp_r_1_2 + comp_r_1_3 + comp_r_1_4 + comp_r_1_5 + comp_r_1_6;
	position[2] = 1 + (sc_uint<1>) ~(comp_r_1_2);
	position[3] = 2 + (sc_uint<1>) ~(comp_r_1_3) + comp_r_3_4 + comp_r_3_5 + comp_r_3_6;
	position[4] = 2 + (sc_uint<1>) ~(comp_r_1_4) + (sc_uint<1>) ~(comp_r_3_4);
	position[5] = 3 + (sc_uint<1>) ~(comp_r_1_5) + (sc_uint<1>) ~(comp_r_3_5) + comp_r_5_6;
	position[6] = 3 + (sc_uint<1>) ~(comp_r_1_6) + (sc_uint<1>) ~(comp_r_3_6) + (sc_uint<1>) ~(comp_r_5_6);
	position[7] = (sc_uint<3>) 7;

// Multiplex the inputs in the order
	PS_struct<1,Q,2> temp[7];
#pragma HLS ARRAY_PARTITION variable=temp complete dim=1
	temp[0] = input[0];
	temp[1] = RO_MUX8 < PS_struct<1,Q,2> > (input, position, (sc_uint<3>) 1);
	temp[2] = RO_MUX8 < PS_struct<1,Q,2> > (input, position, (sc_uint<3>) 2);
	temp[3] = RO_MUX8 < PS_struct<1,Q,2> > (input, position, (sc_uint<3>) 3);
	temp[4] = RO_MUX8 < PS_struct<1,Q,2> > (input, position, (sc_uint<3>) 4);
	temp[5] = RO_MUX8 < PS_struct<1,Q,2> > (input, position, (sc_uint<3>) 5);
	temp[6] = RO_MUX8 < PS_struct<1,Q,2> > (input, position, (sc_uint<3>) 6);


// output the L smallest input
	if (fb == 0){
		output[0] = temp[3];
		output[1] = temp[4];
		output[2] = temp[5];
		output[3] = temp[6];
	}
	else{
		output[0] = temp[0];
		output[1] = temp[1];
		output[2] = temp[2];
		output[3] = temp[3];
	}
}

template <int Q>
void RANKORDER_SORT_L8 (PS_struct<1,Q,3> input[16], PS_struct<1,Q,3> output[8], sc_biguint<1> fb)
{
#pragma HLS INLINE
#pragma HLS ARRAY_PARTITION variable=input complete dim=1
#pragma HLS ARRAY_PARTITION variable=output complete dim=1
// Adjust input when FB=0

	if (fb == 0){
		input[15].metric = MAX_VAL;
		input[14] = input[7];
		input[13] = input[6];
		input[12].metric = 0;
		input[11] = input[5];
		input[10].metric = 0;
		input[9] = input[4];
		input[8].metric = 0;
		input[7] = input[3]; 
		input[6].metric = 0;
		input[5] = input[2];
		input[4].metric = 0;
		input[3] = input[1];
		input[2].metric = 0;
		input[1] = input[0];
	}

	sc_uint<4> position[16];
#pragma HLS ARRAY_PARTITION variable=position complete dim=1

// COMPARE inputs

	sc_uint<1> comp_r_1_2 = COMP < Q > (input[1].metric, input[2].metric);
	sc_uint<1> comp_r_1_3 = COMP < Q > (input[1].metric, input[3].metric);
	sc_uint<1> comp_r_1_4 = COMP < Q > (input[1].metric, input[4].metric);
	sc_uint<1> comp_r_1_5 = COMP < Q > (input[1].metric, input[5].metric);
	sc_uint<1> comp_r_1_6 = COMP < Q > (input[1].metric, input[6].metric);
	sc_uint<1> comp_r_1_7 = COMP < Q > (input[1].metric, input[7].metric);
	sc_uint<1> comp_r_1_8 = COMP < Q > (input[1].metric, input[8].metric);
	sc_uint<1> comp_r_1_9 = COMP < Q > (input[1].metric, input[9].metric);
	sc_uint<1> comp_r_1_10 = COMP < Q > (input[1].metric, input[10].metric);
	sc_uint<1> comp_r_1_11 = COMP < Q > (input[1].metric, input[11].metric);
	sc_uint<1> comp_r_1_12 = COMP < Q > (input[1].metric, input[12].metric);
	sc_uint<1> comp_r_1_13 = COMP < Q > (input[1].metric, input[13].metric);
	sc_uint<1> comp_r_1_14 = COMP < Q > (input[1].metric, input[14].metric);

	sc_uint<1> comp_r_3_4 = COMP < Q > (input[3].metric, input[4].metric);
	sc_uint<1> comp_r_3_5 = COMP < Q > (input[3].metric, input[5].metric);
	sc_uint<1> comp_r_3_6 = COMP < Q > (input[3].metric, input[6].metric);
	sc_uint<1> comp_r_3_7 = COMP < Q > (input[3].metric, input[7].metric);
	sc_uint<1> comp_r_3_8 = COMP < Q > (input[3].metric, input[8].metric);
	sc_uint<1> comp_r_3_9 = COMP < Q > (input[3].metric, input[9].metric);
	sc_uint<1> comp_r_3_10 = COMP < Q > (input[3].metric, input[10].metric);
	sc_uint<1> comp_r_3_11 = COMP < Q > (input[3].metric, input[11].metric);
	sc_uint<1> comp_r_3_12 = COMP < Q > (input[3].metric, input[12].metric);
	sc_uint<1> comp_r_3_13 = COMP < Q > (input[3].metric, input[13].metric);
	sc_uint<1> comp_r_3_14 = COMP < Q > (input[3].metric, input[14].metric);
	
	sc_uint<1> comp_r_5_6 = COMP < Q > (input[5].metric, input[6].metric);
	sc_uint<1> comp_r_5_7 = COMP < Q > (input[5].metric, input[7].metric);
	sc_uint<1> comp_r_5_8 = COMP < Q > (input[5].metric, input[8].metric);
	sc_uint<1> comp_r_5_9 = COMP < Q > (input[5].metric, input[9].metric);
	sc_uint<1> comp_r_5_10 = COMP < Q > (input[5].metric, input[10].metric);
	sc_uint<1> comp_r_5_11 = COMP < Q > (input[5].metric, input[11].metric);
	sc_uint<1> comp_r_5_12 = COMP < Q > (input[5].metric, input[12].metric);
	sc_uint<1> comp_r_5_13 = COMP < Q > (input[5].metric, input[13].metric);
	sc_uint<1> comp_r_5_14 = COMP < Q > (input[5].metric, input[14].metric);
	
	sc_uint<1> comp_r_7_8 = COMP < Q > (input[7].metric, input[8].metric);
	sc_uint<1> comp_r_7_9 = COMP < Q > (input[7].metric, input[9].metric);
	sc_uint<1> comp_r_7_10 = COMP < Q > (input[7].metric, input[10].metric);
	sc_uint<1> comp_r_7_11 = COMP < Q > (input[7].metric, input[11].metric);
	sc_uint<1> comp_r_7_12 = COMP < Q > (input[7].metric, input[12].metric);
	sc_uint<1> comp_r_7_13 = COMP < Q > (input[7].metric, input[13].metric);
	sc_uint<1> comp_r_7_14 = COMP < Q > (input[7].metric, input[14].metric);
	
	sc_uint<1> comp_r_9_10 = COMP < Q > (input[9].metric, input[10].metric);
	sc_uint<1> comp_r_9_11 = COMP < Q > (input[9].metric, input[11].metric);
	sc_uint<1> comp_r_9_12 = COMP < Q > (input[9].metric, input[12].metric);
	sc_uint<1> comp_r_9_13 = COMP < Q > (input[9].metric, input[13].metric);
	sc_uint<1> comp_r_9_14 = COMP < Q > (input[9].metric, input[14].metric);
	
	sc_uint<1> comp_r_11_12 = COMP < Q > (input[11].metric, input[12].metric);
	sc_uint<1> comp_r_11_13 = COMP < Q > (input[11].metric, input[13].metric);
	sc_uint<1> comp_r_11_14 = COMP < Q > (input[11].metric, input[14].metric);
	
	sc_uint<1> comp_r_13_14 = COMP < Q > (input[13].metric, input[14].metric);
	
// Compute positions
 
	position[0] = 0;
	position[1] = 1 + comp_r_1_2 + comp_r_1_3 + comp_r_1_4 + comp_r_1_5 + comp_r_1_6 + comp_r_1_7 + comp_r_1_8 + comp_r_1_9 + comp_r_1_10 + comp_r_1_11 + comp_r_1_12 + comp_r_1_13 + comp_r_1_14;
	position[2] = 1 + (sc_uint<1>) ~(comp_r_1_2);
	position[3] = 2 + (sc_uint<1>) ~(comp_r_1_3) + comp_r_3_4 + comp_r_3_5 + comp_r_3_6 + comp_r_3_7 + comp_r_3_8 + comp_r_3_9 + comp_r_3_10 + comp_r_3_11 + comp_r_3_12 + comp_r_3_13 + comp_r_3_14;
	position[4] = 2 + (sc_uint<1>) ~(comp_r_1_4) + (sc_uint<1>) ~(comp_r_3_4);
	position[5] = 3 + (sc_uint<1>) ~(comp_r_1_5) + + (sc_uint<1>) ~(comp_r_3_5) + comp_r_5_6 + comp_r_5_7 + comp_r_5_8 + comp_r_5_9 + comp_r_5_10 + comp_r_5_11 + comp_r_5_12 + comp_r_5_13 + comp_r_5_14;
	position[6] = 3 + (sc_uint<1>) ~(comp_r_1_6) + (sc_uint<1>) ~(comp_r_3_6) + (sc_uint<1>) ~(comp_r_5_6);
	position[7] = 4 + (sc_uint<1>) ~(comp_r_1_7) + + (sc_uint<1>) ~(comp_r_3_7) + (sc_uint<1>) ~(comp_r_5_7) + comp_r_7_8 + comp_r_7_9 + comp_r_7_10 + comp_r_7_11 + comp_r_7_12 + comp_r_7_13 + comp_r_7_14;
	position[8] = 4 + (sc_uint<1>) ~(comp_r_1_8) + (sc_uint<1>) ~(comp_r_3_8) + (sc_uint<1>) ~(comp_r_5_8) + (sc_uint<1>) ~(comp_r_7_8);
	position[9] = 5 + (sc_uint<1>) ~(comp_r_1_9) + + (sc_uint<1>) ~(comp_r_3_9) + (sc_uint<1>) ~(comp_r_5_9) + (sc_uint<1>) ~(comp_r_7_9) + comp_r_9_10 + comp_r_9_11 + comp_r_9_12 + comp_r_9_13 + comp_r_9_14;
	position[10] = 5 + (sc_uint<1>) ~(comp_r_1_10) + (sc_uint<1>) ~(comp_r_3_10) + (sc_uint<1>) ~(comp_r_5_10) + (sc_uint<1>) ~(comp_r_7_10) + (sc_uint<1>) ~(comp_r_9_10);
	position[11] = 6 + (sc_uint<1>) ~(comp_r_1_11) + + (sc_uint<1>) ~(comp_r_3_11) + (sc_uint<1>) ~(comp_r_5_11) + (sc_uint<1>) ~(comp_r_7_11) + (sc_uint<1>) ~(comp_r_9_11) + comp_r_11_12 + comp_r_11_13 + comp_r_11_14;
	position[12] = 6 + (sc_uint<1>) ~(comp_r_1_12) + (sc_uint<1>) ~(comp_r_3_12) + (sc_uint<1>) ~(comp_r_5_12) + (sc_uint<1>) ~(comp_r_7_12) + (sc_uint<1>) ~(comp_r_9_12) + (sc_uint<1>) ~(comp_r_11_12);
	position[13] = 7 + (sc_uint<1>) ~(comp_r_1_13) + + (sc_uint<1>) ~(comp_r_3_13) + (sc_uint<1>) ~(comp_r_5_13) + (sc_uint<1>) ~(comp_r_7_13) + (sc_uint<1>) ~(comp_r_9_13) + (sc_uint<1>) ~(comp_r_11_13) + comp_r_13_14;
	position[14] = 7 + (sc_uint<1>) ~(comp_r_1_14) + (sc_uint<1>) ~(comp_r_3_14) + (sc_uint<1>) ~(comp_r_5_14) + (sc_uint<1>) ~(comp_r_7_14) + (sc_uint<1>) ~(comp_r_9_14) + (sc_uint<1>) ~(comp_r_11_14) + (sc_uint<1>) ~(comp_r_13_14);
	position[15] = (sc_uint<4>) 15;

// Multiplex the inputs in the order
	PS_struct<1,Q,3> temp[15];
#pragma HLS ARRAY_PARTITION variable=temp complete dim=1
	temp[0] = input[0];
	temp[1] = RO_MUX16 < PS_struct<1,Q,3> > (input, position, (sc_uint<4>) 1);
	temp[2] = RO_MUX16 < PS_struct<1,Q,3> > (input, position, (sc_uint<4>) 2);
	temp[3] = RO_MUX16 < PS_struct<1,Q,3> > (input, position, (sc_uint<4>) 3);
	temp[4] = RO_MUX16 < PS_struct<1,Q,3> > (input, position, (sc_uint<4>) 4);
	temp[5] = RO_MUX16 < PS_struct<1,Q,3> > (input, position, (sc_uint<4>) 5);
	temp[6] = RO_MUX16 < PS_struct<1,Q,3> > (input, position, (sc_uint<4>) 6);
	temp[7] = RO_MUX16 < PS_struct<1,Q,3> > (input, position, (sc_uint<4>) 7);
	temp[8] = RO_MUX16 < PS_struct<1,Q,3> > (input, position, (sc_uint<4>) 8);
	temp[9] = RO_MUX16 < PS_struct<1,Q,3> > (input, position, (sc_uint<4>) 9);
	temp[10] = RO_MUX16 < PS_struct<1,Q,3> > (input, position, (sc_uint<4>) 10);
	temp[11] = RO_MUX16 < PS_struct<1,Q,3> > (input, position, (sc_uint<4>) 11);
	temp[12] = RO_MUX16 < PS_struct<1,Q,3> > (input, position, (sc_uint<4>) 12);
	temp[13] = RO_MUX16 < PS_struct<1,Q,3> > (input, position, (sc_uint<4>) 13);
	temp[14] = RO_MUX16 < PS_struct<1,Q,3> > (input, position, (sc_uint<4>) 14);

// output the L smallest input
	if (fb == 0){
		output[0] = temp[7];
		output[1] = temp[8];
		output[2] = temp[9];
		output[3] = temp[10];
		output[4] = temp[11];
		output[5] = temp[12];
		output[6] = temp[13];
		output[7] = temp[14];
	}
	else{
		output[0] = temp[0];
		output[1] = temp[1];
		output[2] = temp[2];
		output[3] = temp[3];
		output[4] = temp[4];
		output[5] = temp[5];
		output[6] = temp[6];
		output[7] = temp[7];
	}
}

template <int Q>
void RANKORDER_SORT_L16 (PS_struct<1,Q,4> input[32], PS_struct<1,Q,4> output[16], sc_biguint<1> fb)
{
#pragma HLS INLINE off
#pragma HLS ARRAY_PARTITION variable=input complete dim=1
#pragma HLS ARRAY_PARTITION variable=output complete dim=1

// Adjust input when FB=0

    if (fb == 0){
        input[31].metric = MAX_VAL;
        input[30] = input[15];
        input[29] = input[14];
        input[28].metric = 0;
        input[27] = input[13];
        input[26].metric = 0;
        input[25] = input[12];
        input[24].metric = 0;
        input[23] = input[11];
        input[22].metric = 0;
        input[21] = input[10];
        input[20].metric = 0;
        input[19] = input[9];
        input[18].metric = 0;
        input[17] = input[8];
        input[16].metric = 0;
        input[15] = input[7];
        input[14].metric = 0;
        input[13] = input[6];
        input[12].metric = 0;
        input[11] = input[5];
        input[10].metric = 0;
        input[9] = input[4];
        input[8].metric = 0;
        input[7] = input[3];
        input[6].metric = 0;
        input[5] = input[2];
        input[4].metric = 0;
        input[3] = input[1];
        input[2].metric = 0;
        input[1] = input[0];
    }

    sc_uint<5> position[32];
#pragma HLS ARRAY_PARTITION variable=position complete dim=1

// COMPARE inputs

    sc_uint<1> comp_r_1_2 = COMP < Q > (input[1].metric, input[2].metric);
    sc_uint<1> comp_r_1_3 = COMP < Q > (input[1].metric, input[3].metric);
    sc_uint<1> comp_r_1_4 = COMP < Q > (input[1].metric, input[4].metric);
    sc_uint<1> comp_r_1_5 = COMP < Q > (input[1].metric, input[5].metric);
    sc_uint<1> comp_r_1_6 = COMP < Q > (input[1].metric, input[6].metric);
    sc_uint<1> comp_r_1_7 = COMP < Q > (input[1].metric, input[7].metric);
    sc_uint<1> comp_r_1_8 = COMP < Q > (input[1].metric, input[8].metric);
    sc_uint<1> comp_r_1_9 = COMP < Q > (input[1].metric, input[9].metric);
    sc_uint<1> comp_r_1_10 = COMP < Q > (input[1].metric, input[10].metric);
    sc_uint<1> comp_r_1_11 = COMP < Q > (input[1].metric, input[11].metric);
    sc_uint<1> comp_r_1_12 = COMP < Q > (input[1].metric, input[12].metric);
    sc_uint<1> comp_r_1_13 = COMP < Q > (input[1].metric, input[13].metric);
    sc_uint<1> comp_r_1_14 = COMP < Q > (input[1].metric, input[14].metric);
    sc_uint<1> comp_r_1_15 = COMP < Q > (input[1].metric, input[15].metric);
    sc_uint<1> comp_r_1_16 = COMP < Q > (input[1].metric, input[16].metric);
    sc_uint<1> comp_r_1_17 = COMP < Q > (input[1].metric, input[17].metric);
    sc_uint<1> comp_r_1_18 = COMP < Q > (input[1].metric, input[18].metric);
    sc_uint<1> comp_r_1_19 = COMP < Q > (input[1].metric, input[19].metric);
    sc_uint<1> comp_r_1_20 = COMP < Q > (input[1].metric, input[20].metric);
    sc_uint<1> comp_r_1_21 = COMP < Q > (input[1].metric, input[21].metric);
    sc_uint<1> comp_r_1_22 = COMP < Q > (input[1].metric, input[22].metric);
    sc_uint<1> comp_r_1_23 = COMP < Q > (input[1].metric, input[23].metric);
    sc_uint<1> comp_r_1_24 = COMP < Q > (input[1].metric, input[24].metric);
    sc_uint<1> comp_r_1_25 = COMP < Q > (input[1].metric, input[25].metric);
    sc_uint<1> comp_r_1_26 = COMP < Q > (input[1].metric, input[26].metric);
    sc_uint<1> comp_r_1_27 = COMP < Q > (input[1].metric, input[27].metric);
    sc_uint<1> comp_r_1_28 = COMP < Q > (input[1].metric, input[28].metric);
    sc_uint<1> comp_r_1_29 = COMP < Q > (input[1].metric, input[29].metric);
    sc_uint<1> comp_r_1_30 = COMP < Q > (input[1].metric, input[30].metric);

    sc_uint<1> comp_r_3_4 = COMP < Q > (input[3].metric, input[4].metric);
    sc_uint<1> comp_r_3_5 = COMP < Q > (input[3].metric, input[5].metric);
    sc_uint<1> comp_r_3_6 = COMP < Q > (input[3].metric, input[6].metric);
    sc_uint<1> comp_r_3_7 = COMP < Q > (input[3].metric, input[7].metric);
    sc_uint<1> comp_r_3_8 = COMP < Q > (input[3].metric, input[8].metric);
    sc_uint<1> comp_r_3_9 = COMP < Q > (input[3].metric, input[9].metric);
    sc_uint<1> comp_r_3_10 = COMP < Q > (input[3].metric, input[10].metric);
    sc_uint<1> comp_r_3_11 = COMP < Q > (input[3].metric, input[11].metric);
    sc_uint<1> comp_r_3_12 = COMP < Q > (input[3].metric, input[12].metric);
    sc_uint<1> comp_r_3_13 = COMP < Q > (input[3].metric, input[13].metric);
    sc_uint<1> comp_r_3_14 = COMP < Q > (input[3].metric, input[14].metric);
    sc_uint<1> comp_r_3_15 = COMP < Q > (input[3].metric, input[15].metric);
    sc_uint<1> comp_r_3_16 = COMP < Q > (input[3].metric, input[16].metric);
    sc_uint<1> comp_r_3_17 = COMP < Q > (input[3].metric, input[17].metric);
    sc_uint<1> comp_r_3_18 = COMP < Q > (input[3].metric, input[18].metric);
    sc_uint<1> comp_r_3_19 = COMP < Q > (input[3].metric, input[19].metric);
    sc_uint<1> comp_r_3_20 = COMP < Q > (input[3].metric, input[20].metric);
    sc_uint<1> comp_r_3_21 = COMP < Q > (input[3].metric, input[21].metric);
    sc_uint<1> comp_r_3_22 = COMP < Q > (input[3].metric, input[22].metric);
    sc_uint<1> comp_r_3_23 = COMP < Q > (input[3].metric, input[23].metric);
    sc_uint<1> comp_r_3_24 = COMP < Q > (input[3].metric, input[24].metric);
    sc_uint<1> comp_r_3_25 = COMP < Q > (input[3].metric, input[25].metric);
    sc_uint<1> comp_r_3_26 = COMP < Q > (input[3].metric, input[26].metric);
    sc_uint<1> comp_r_3_27 = COMP < Q > (input[3].metric, input[27].metric);
    sc_uint<1> comp_r_3_28 = COMP < Q > (input[3].metric, input[28].metric);
    sc_uint<1> comp_r_3_29 = COMP < Q > (input[3].metric, input[29].metric);
    sc_uint<1> comp_r_3_30 = COMP < Q > (input[3].metric, input[30].metric);

    sc_uint<1> comp_r_5_6 = COMP < Q > (input[5].metric, input[6].metric);
    sc_uint<1> comp_r_5_7 = COMP < Q > (input[5].metric, input[7].metric);
    sc_uint<1> comp_r_5_8 = COMP < Q > (input[5].metric, input[8].metric);
    sc_uint<1> comp_r_5_9 = COMP < Q > (input[5].metric, input[9].metric);
    sc_uint<1> comp_r_5_10 = COMP < Q > (input[5].metric, input[10].metric);
    sc_uint<1> comp_r_5_11 = COMP < Q > (input[5].metric, input[11].metric);
    sc_uint<1> comp_r_5_12 = COMP < Q > (input[5].metric, input[12].metric);
    sc_uint<1> comp_r_5_13 = COMP < Q > (input[5].metric, input[13].metric);
    sc_uint<1> comp_r_5_14 = COMP < Q > (input[5].metric, input[14].metric);
    sc_uint<1> comp_r_5_15 = COMP < Q > (input[5].metric, input[15].metric);
    sc_uint<1> comp_r_5_16 = COMP < Q > (input[5].metric, input[16].metric);
    sc_uint<1> comp_r_5_17 = COMP < Q > (input[5].metric, input[17].metric);
    sc_uint<1> comp_r_5_18 = COMP < Q > (input[5].metric, input[18].metric);
    sc_uint<1> comp_r_5_19 = COMP < Q > (input[5].metric, input[19].metric);
    sc_uint<1> comp_r_5_20 = COMP < Q > (input[5].metric, input[20].metric);
    sc_uint<1> comp_r_5_21 = COMP < Q > (input[5].metric, input[21].metric);
    sc_uint<1> comp_r_5_22 = COMP < Q > (input[5].metric, input[22].metric);
    sc_uint<1> comp_r_5_23 = COMP < Q > (input[5].metric, input[23].metric);
    sc_uint<1> comp_r_5_24 = COMP < Q > (input[5].metric, input[24].metric);
    sc_uint<1> comp_r_5_25 = COMP < Q > (input[5].metric, input[25].metric);
    sc_uint<1> comp_r_5_26 = COMP < Q > (input[5].metric, input[26].metric);
    sc_uint<1> comp_r_5_27 = COMP < Q > (input[5].metric, input[27].metric);
    sc_uint<1> comp_r_5_28 = COMP < Q > (input[5].metric, input[28].metric);
    sc_uint<1> comp_r_5_29 = COMP < Q > (input[5].metric, input[29].metric);
    sc_uint<1> comp_r_5_30 = COMP < Q > (input[5].metric, input[30].metric);

    sc_uint<1> comp_r_7_8 = COMP < Q > (input[7].metric, input[8].metric);
    sc_uint<1> comp_r_7_9 = COMP < Q > (input[7].metric, input[9].metric);
    sc_uint<1> comp_r_7_10 = COMP < Q > (input[7].metric, input[10].metric);
    sc_uint<1> comp_r_7_11 = COMP < Q > (input[7].metric, input[11].metric);
    sc_uint<1> comp_r_7_12 = COMP < Q > (input[7].metric, input[12].metric);
    sc_uint<1> comp_r_7_13 = COMP < Q > (input[7].metric, input[13].metric);
    sc_uint<1> comp_r_7_14 = COMP < Q > (input[7].metric, input[14].metric);
    sc_uint<1> comp_r_7_15 = COMP < Q > (input[7].metric, input[15].metric);
    sc_uint<1> comp_r_7_16 = COMP < Q > (input[7].metric, input[16].metric);
    sc_uint<1> comp_r_7_17 = COMP < Q > (input[7].metric, input[17].metric);
    sc_uint<1> comp_r_7_18 = COMP < Q > (input[7].metric, input[18].metric);
    sc_uint<1> comp_r_7_19 = COMP < Q > (input[7].metric, input[19].metric);
    sc_uint<1> comp_r_7_20 = COMP < Q > (input[7].metric, input[20].metric);
    sc_uint<1> comp_r_7_21 = COMP < Q > (input[7].metric, input[21].metric);
    sc_uint<1> comp_r_7_22 = COMP < Q > (input[7].metric, input[22].metric);
    sc_uint<1> comp_r_7_23 = COMP < Q > (input[7].metric, input[23].metric);
    sc_uint<1> comp_r_7_24 = COMP < Q > (input[7].metric, input[24].metric);
    sc_uint<1> comp_r_7_25 = COMP < Q > (input[7].metric, input[25].metric);
    sc_uint<1> comp_r_7_26 = COMP < Q > (input[7].metric, input[26].metric);
    sc_uint<1> comp_r_7_27 = COMP < Q > (input[7].metric, input[27].metric);
    sc_uint<1> comp_r_7_28 = COMP < Q > (input[7].metric, input[28].metric);
    sc_uint<1> comp_r_7_29 = COMP < Q > (input[7].metric, input[29].metric);
    sc_uint<1> comp_r_7_30 = COMP < Q > (input[7].metric, input[30].metric);

    sc_uint<1> comp_r_9_10 = COMP < Q > (input[9].metric, input[10].metric);
    sc_uint<1> comp_r_9_11 = COMP < Q > (input[9].metric, input[11].metric);
    sc_uint<1> comp_r_9_12 = COMP < Q > (input[9].metric, input[12].metric);
    sc_uint<1> comp_r_9_13 = COMP < Q > (input[9].metric, input[13].metric);
    sc_uint<1> comp_r_9_14 = COMP < Q > (input[9].metric, input[14].metric);
    sc_uint<1> comp_r_9_15 = COMP < Q > (input[9].metric, input[15].metric);
    sc_uint<1> comp_r_9_16 = COMP < Q > (input[9].metric, input[16].metric);
    sc_uint<1> comp_r_9_17 = COMP < Q > (input[9].metric, input[17].metric);
    sc_uint<1> comp_r_9_18 = COMP < Q > (input[9].metric, input[18].metric);
    sc_uint<1> comp_r_9_19 = COMP < Q > (input[9].metric, input[19].metric);
    sc_uint<1> comp_r_9_20 = COMP < Q > (input[9].metric, input[20].metric);
    sc_uint<1> comp_r_9_21 = COMP < Q > (input[9].metric, input[21].metric);
    sc_uint<1> comp_r_9_22 = COMP < Q > (input[9].metric, input[22].metric);
    sc_uint<1> comp_r_9_23 = COMP < Q > (input[9].metric, input[23].metric);
    sc_uint<1> comp_r_9_24 = COMP < Q > (input[9].metric, input[24].metric);
    sc_uint<1> comp_r_9_25 = COMP < Q > (input[9].metric, input[25].metric);
    sc_uint<1> comp_r_9_26 = COMP < Q > (input[9].metric, input[26].metric);
    sc_uint<1> comp_r_9_27 = COMP < Q > (input[9].metric, input[27].metric);
    sc_uint<1> comp_r_9_28 = COMP < Q > (input[9].metric, input[28].metric);
    sc_uint<1> comp_r_9_29 = COMP < Q > (input[9].metric, input[29].metric);
    sc_uint<1> comp_r_9_30 = COMP < Q > (input[9].metric, input[30].metric);

    sc_uint<1> comp_r_11_12 = COMP < Q > (input[11].metric, input[12].metric);
    sc_uint<1> comp_r_11_13 = COMP < Q > (input[11].metric, input[13].metric);
    sc_uint<1> comp_r_11_14 = COMP < Q > (input[11].metric, input[14].metric);
    sc_uint<1> comp_r_11_15 = COMP < Q > (input[11].metric, input[15].metric);
    sc_uint<1> comp_r_11_16 = COMP < Q > (input[11].metric, input[16].metric);
    sc_uint<1> comp_r_11_17 = COMP < Q > (input[11].metric, input[17].metric);
    sc_uint<1> comp_r_11_18 = COMP < Q > (input[11].metric, input[18].metric);
    sc_uint<1> comp_r_11_19 = COMP < Q > (input[11].metric, input[19].metric);
    sc_uint<1> comp_r_11_20 = COMP < Q > (input[11].metric, input[20].metric);
    sc_uint<1> comp_r_11_21 = COMP < Q > (input[11].metric, input[21].metric);
    sc_uint<1> comp_r_11_22 = COMP < Q > (input[11].metric, input[22].metric);
    sc_uint<1> comp_r_11_23 = COMP < Q > (input[11].metric, input[23].metric);
    sc_uint<1> comp_r_11_24 = COMP < Q > (input[11].metric, input[24].metric);
    sc_uint<1> comp_r_11_25 = COMP < Q > (input[11].metric, input[25].metric);
    sc_uint<1> comp_r_11_26 = COMP < Q > (input[11].metric, input[26].metric);
    sc_uint<1> comp_r_11_27 = COMP < Q > (input[11].metric, input[27].metric);
    sc_uint<1> comp_r_11_28 = COMP < Q > (input[11].metric, input[28].metric);
    sc_uint<1> comp_r_11_29 = COMP < Q > (input[11].metric, input[29].metric);
    sc_uint<1> comp_r_11_30 = COMP < Q > (input[11].metric, input[30].metric);

    sc_uint<1> comp_r_13_14 = COMP < Q > (input[13].metric, input[14].metric);
    sc_uint<1> comp_r_13_15 = COMP < Q > (input[13].metric, input[15].metric);
    sc_uint<1> comp_r_13_16 = COMP < Q > (input[13].metric, input[16].metric);
    sc_uint<1> comp_r_13_17 = COMP < Q > (input[13].metric, input[17].metric);
    sc_uint<1> comp_r_13_18 = COMP < Q > (input[13].metric, input[18].metric);
    sc_uint<1> comp_r_13_19 = COMP < Q > (input[13].metric, input[19].metric);
    sc_uint<1> comp_r_13_20 = COMP < Q > (input[13].metric, input[20].metric);
    sc_uint<1> comp_r_13_21 = COMP < Q > (input[13].metric, input[21].metric);
    sc_uint<1> comp_r_13_22 = COMP < Q > (input[13].metric, input[22].metric);
    sc_uint<1> comp_r_13_23 = COMP < Q > (input[13].metric, input[23].metric);
    sc_uint<1> comp_r_13_24 = COMP < Q > (input[13].metric, input[24].metric);
    sc_uint<1> comp_r_13_25 = COMP < Q > (input[13].metric, input[25].metric);
    sc_uint<1> comp_r_13_26 = COMP < Q > (input[13].metric, input[26].metric);
    sc_uint<1> comp_r_13_27 = COMP < Q > (input[13].metric, input[27].metric);
    sc_uint<1> comp_r_13_28 = COMP < Q > (input[13].metric, input[28].metric);
    sc_uint<1> comp_r_13_29 = COMP < Q > (input[13].metric, input[29].metric);
    sc_uint<1> comp_r_13_30 = COMP < Q > (input[13].metric, input[30].metric);

    sc_uint<1> comp_r_15_16 = COMP < Q > (input[15].metric, input[16].metric);
    sc_uint<1> comp_r_15_17 = COMP < Q > (input[15].metric, input[17].metric);
    sc_uint<1> comp_r_15_18 = COMP < Q > (input[15].metric, input[18].metric);
    sc_uint<1> comp_r_15_19 = COMP < Q > (input[15].metric, input[19].metric);
    sc_uint<1> comp_r_15_20 = COMP < Q > (input[15].metric, input[20].metric);
    sc_uint<1> comp_r_15_21 = COMP < Q > (input[15].metric, input[21].metric);
    sc_uint<1> comp_r_15_22 = COMP < Q > (input[15].metric, input[22].metric);
    sc_uint<1> comp_r_15_23 = COMP < Q > (input[15].metric, input[23].metric);
    sc_uint<1> comp_r_15_24 = COMP < Q > (input[15].metric, input[24].metric);
    sc_uint<1> comp_r_15_25 = COMP < Q > (input[15].metric, input[25].metric);
    sc_uint<1> comp_r_15_26 = COMP < Q > (input[15].metric, input[26].metric);
    sc_uint<1> comp_r_15_27 = COMP < Q > (input[15].metric, input[27].metric);
    sc_uint<1> comp_r_15_28 = COMP < Q > (input[15].metric, input[28].metric);
    sc_uint<1> comp_r_15_29 = COMP < Q > (input[15].metric, input[29].metric);
    sc_uint<1> comp_r_15_30 = COMP < Q > (input[15].metric, input[30].metric);

    sc_uint<1> comp_r_17_18 = COMP < Q > (input[17].metric, input[18].metric);
    sc_uint<1> comp_r_17_19 = COMP < Q > (input[17].metric, input[19].metric);
    sc_uint<1> comp_r_17_20 = COMP < Q > (input[17].metric, input[20].metric);
    sc_uint<1> comp_r_17_21 = COMP < Q > (input[17].metric, input[21].metric);
    sc_uint<1> comp_r_17_22 = COMP < Q > (input[17].metric, input[22].metric);
    sc_uint<1> comp_r_17_23 = COMP < Q > (input[17].metric, input[23].metric);
    sc_uint<1> comp_r_17_24 = COMP < Q > (input[17].metric, input[24].metric);
    sc_uint<1> comp_r_17_25 = COMP < Q > (input[17].metric, input[25].metric);
    sc_uint<1> comp_r_17_26 = COMP < Q > (input[17].metric, input[26].metric);
    sc_uint<1> comp_r_17_27 = COMP < Q > (input[17].metric, input[27].metric);
    sc_uint<1> comp_r_17_28 = COMP < Q > (input[17].metric, input[28].metric);
    sc_uint<1> comp_r_17_29 = COMP < Q > (input[17].metric, input[29].metric);
    sc_uint<1> comp_r_17_30 = COMP < Q > (input[17].metric, input[30].metric);

    sc_uint<1> comp_r_19_20 = COMP < Q > (input[19].metric, input[20].metric);
    sc_uint<1> comp_r_19_21 = COMP < Q > (input[19].metric, input[21].metric);
    sc_uint<1> comp_r_19_22 = COMP < Q > (input[19].metric, input[22].metric);
    sc_uint<1> comp_r_19_23 = COMP < Q > (input[19].metric, input[23].metric);
    sc_uint<1> comp_r_19_24 = COMP < Q > (input[19].metric, input[24].metric);
    sc_uint<1> comp_r_19_25 = COMP < Q > (input[19].metric, input[25].metric);
    sc_uint<1> comp_r_19_26 = COMP < Q > (input[19].metric, input[26].metric);
    sc_uint<1> comp_r_19_27 = COMP < Q > (input[19].metric, input[27].metric);
    sc_uint<1> comp_r_19_28 = COMP < Q > (input[19].metric, input[28].metric);
    sc_uint<1> comp_r_19_29 = COMP < Q > (input[19].metric, input[29].metric);
    sc_uint<1> comp_r_19_30 = COMP < Q > (input[19].metric, input[30].metric);

    sc_uint<1> comp_r_21_22 = COMP < Q > (input[21].metric, input[22].metric);
    sc_uint<1> comp_r_21_23 = COMP < Q > (input[21].metric, input[23].metric);
    sc_uint<1> comp_r_21_24 = COMP < Q > (input[21].metric, input[24].metric);
    sc_uint<1> comp_r_21_25 = COMP < Q > (input[21].metric, input[25].metric);
    sc_uint<1> comp_r_21_26 = COMP < Q > (input[21].metric, input[26].metric);
    sc_uint<1> comp_r_21_27 = COMP < Q > (input[21].metric, input[27].metric);
    sc_uint<1> comp_r_21_28 = COMP < Q > (input[21].metric, input[28].metric);
    sc_uint<1> comp_r_21_29 = COMP < Q > (input[21].metric, input[29].metric);
    sc_uint<1> comp_r_21_30 = COMP < Q > (input[21].metric, input[30].metric);

    sc_uint<1> comp_r_23_24 = COMP < Q > (input[23].metric, input[24].metric);
    sc_uint<1> comp_r_23_25 = COMP < Q > (input[23].metric, input[25].metric);
    sc_uint<1> comp_r_23_26 = COMP < Q > (input[23].metric, input[26].metric);
    sc_uint<1> comp_r_23_27 = COMP < Q > (input[23].metric, input[27].metric);
    sc_uint<1> comp_r_23_28 = COMP < Q > (input[23].metric, input[28].metric);
    sc_uint<1> comp_r_23_29 = COMP < Q > (input[23].metric, input[29].metric);
    sc_uint<1> comp_r_23_30 = COMP < Q > (input[23].metric, input[30].metric);

    sc_uint<1> comp_r_25_26 = COMP < Q > (input[25].metric, input[26].metric);
    sc_uint<1> comp_r_25_27 = COMP < Q > (input[25].metric, input[27].metric);
    sc_uint<1> comp_r_25_28 = COMP < Q > (input[25].metric, input[28].metric);
    sc_uint<1> comp_r_25_29 = COMP < Q > (input[25].metric, input[29].metric);
    sc_uint<1> comp_r_25_30 = COMP < Q > (input[25].metric, input[30].metric);

    sc_uint<1> comp_r_27_28 = COMP < Q > (input[27].metric, input[28].metric);
    sc_uint<1> comp_r_27_29 = COMP < Q > (input[27].metric, input[29].metric);
    sc_uint<1> comp_r_27_30 = COMP < Q > (input[27].metric, input[30].metric);

    sc_uint<1> comp_r_29_30 = COMP < Q > (input[29].metric, input[30].metric);

// COMPUTE positions

    position[0] = 0;
    position[1] = 1 + comp_r_1_2 + comp_r_1_3 + comp_r_1_4 + comp_r_1_5 + comp_r_1_6 + comp_r_1_7 + comp_r_1_8 + comp_r_1_9 + comp_r_1_10 + comp_r_1_11 + comp_r_1_12 + comp_r_1_13 + comp_r_1_14 + comp_r_1_15 + comp_r_1_16 + comp_r_1_17 + comp_r_1_18 + comp_r_1_19 + comp_r_1_20 + comp_r_1_21 + comp_r_1_22 + comp_r_1_23 + comp_r_1_24 + comp_r_1_25 + comp_r_1_26 + comp_r_1_27 + comp_r_1_28 + comp_r_1_29 + comp_r_1_30;
    position[2] = 1 + (sc_uint<1>) ~(comp_r_1_2);
    position[3] = 2 + (sc_uint<1>) ~(comp_r_1_3) + comp_r_3_4 + comp_r_3_5 + comp_r_3_6 + comp_r_3_7 + comp_r_3_8 + comp_r_3_9 + comp_r_3_10 + comp_r_3_11 + comp_r_3_12 + comp_r_3_13 + comp_r_3_14 + comp_r_3_15 + comp_r_3_16 + comp_r_3_17 + comp_r_3_18 + comp_r_3_19 + comp_r_3_20 + comp_r_3_21 + comp_r_3_22 + comp_r_3_23 + comp_r_3_24 + comp_r_3_25 + comp_r_3_26 + comp_r_3_27 + comp_r_3_28 + comp_r_3_29 + comp_r_3_30;
    position[4] = 2 + (sc_uint<1>) ~(comp_r_1_4) + (sc_uint<1>) ~(comp_r_3_4);
    position[5] = 3 + (sc_uint<1>) ~(comp_r_1_5) + (sc_uint<1>) ~(comp_r_3_5) + comp_r_5_6 + comp_r_5_7 + comp_r_5_8 + comp_r_5_9 + comp_r_5_10 + comp_r_5_11 + comp_r_5_12 + comp_r_5_13 + comp_r_5_14 + comp_r_5_15 + comp_r_5_16 + comp_r_5_17 + comp_r_5_18 + comp_r_5_19 + comp_r_5_20 + comp_r_5_21 + comp_r_5_22 + comp_r_5_23 + comp_r_5_24 + comp_r_5_25 + comp_r_5_26 + comp_r_5_27 + comp_r_5_28 + comp_r_5_29 + comp_r_5_30;
    position[6] = 3 + (sc_uint<1>) ~(comp_r_1_6) + (sc_uint<1>) ~(comp_r_3_6) + (sc_uint<1>) ~(comp_r_5_6);
    position[7] = 4 + (sc_uint<1>) ~(comp_r_1_7) + (sc_uint<1>) ~(comp_r_3_7) + (sc_uint<1>) ~(comp_r_5_7) + comp_r_7_8 + comp_r_7_9 + comp_r_7_10 + comp_r_7_11 + comp_r_7_12 + comp_r_7_13 + comp_r_7_14 + comp_r_7_15 + comp_r_7_16 + comp_r_7_17 + comp_r_7_18 + comp_r_7_19 + comp_r_7_20 + comp_r_7_21 + comp_r_7_22 + comp_r_7_23 + comp_r_7_24 + comp_r_7_25 + comp_r_7_26 + comp_r_7_27 + comp_r_7_28 + comp_r_7_29 + comp_r_7_30;
    position[8] = 4 + (sc_uint<1>) ~(comp_r_1_8) + (sc_uint<1>) ~(comp_r_3_8) + (sc_uint<1>) ~(comp_r_5_8) + (sc_uint<1>) ~(comp_r_7_8);
    position[9] = 5 + (sc_uint<1>) ~(comp_r_1_9) + (sc_uint<1>) ~(comp_r_3_9) + (sc_uint<1>) ~(comp_r_5_9) + (sc_uint<1>) ~(comp_r_7_9) + comp_r_9_10 + comp_r_9_11 + comp_r_9_12 + comp_r_9_13 + comp_r_9_14 + comp_r_9_15 + comp_r_9_16 + comp_r_9_17 + comp_r_9_18 + comp_r_9_19 + comp_r_9_20 + comp_r_9_21 + comp_r_9_22 + comp_r_9_23 + comp_r_9_24 + comp_r_9_25 + comp_r_9_26 + comp_r_9_27 + comp_r_9_28 + comp_r_9_29 + comp_r_9_30;
    position[10] = 5 + (sc_uint<1>) ~(comp_r_1_10) + (sc_uint<1>) ~(comp_r_3_10) + (sc_uint<1>) ~(comp_r_5_10) + (sc_uint<1>) ~(comp_r_7_10) + (sc_uint<1>) ~(comp_r_9_10);
    position[11] = 6 + (sc_uint<1>) ~(comp_r_1_11) + (sc_uint<1>) ~(comp_r_3_11) + (sc_uint<1>) ~(comp_r_5_11) + (sc_uint<1>) ~(comp_r_7_11) + (sc_uint<1>) ~(comp_r_9_11) + comp_r_11_12 + comp_r_11_13 + comp_r_11_14 + comp_r_11_15 + comp_r_11_16 + comp_r_11_17 + comp_r_11_18 + comp_r_11_19 + comp_r_11_20 + comp_r_11_21 + comp_r_11_22 + comp_r_11_23 + comp_r_11_24 + comp_r_11_25 + comp_r_11_26 + comp_r_11_27 + comp_r_11_28 + comp_r_11_29 + comp_r_11_30;
    position[12] = 6 + (sc_uint<1>) ~(comp_r_1_12) + (sc_uint<1>) ~(comp_r_3_12) + (sc_uint<1>) ~(comp_r_5_12) + (sc_uint<1>) ~(comp_r_7_12) + (sc_uint<1>) ~(comp_r_9_12) + (sc_uint<1>) ~(comp_r_11_12);
    position[13] = 7 + (sc_uint<1>) ~(comp_r_1_13) + (sc_uint<1>) ~(comp_r_3_13) + (sc_uint<1>) ~(comp_r_5_13) + (sc_uint<1>) ~(comp_r_7_13) + (sc_uint<1>) ~(comp_r_9_13) + (sc_uint<1>) ~(comp_r_11_13) + comp_r_13_14 + comp_r_13_15 + comp_r_13_16 + comp_r_13_17 + comp_r_13_18 + comp_r_13_19 + comp_r_13_20 + comp_r_13_21 + comp_r_13_22 + comp_r_13_23 + comp_r_13_24 + comp_r_13_25 + comp_r_13_26 + comp_r_13_27 + comp_r_13_28 + comp_r_13_29 + comp_r_13_30;
    position[14] = 7 + (sc_uint<1>) ~(comp_r_1_14) + (sc_uint<1>) ~(comp_r_3_14) + (sc_uint<1>) ~(comp_r_5_14) + (sc_uint<1>) ~(comp_r_7_14) + (sc_uint<1>) ~(comp_r_9_14) + (sc_uint<1>) ~(comp_r_11_14) + (sc_uint<1>) ~(comp_r_13_14);
    position[15] = 8 + (sc_uint<1>) ~(comp_r_1_15) + (sc_uint<1>) ~(comp_r_3_15) + (sc_uint<1>) ~(comp_r_5_15) + (sc_uint<1>) ~(comp_r_7_15) + (sc_uint<1>) ~(comp_r_9_15) + (sc_uint<1>) ~(comp_r_11_15) + (sc_uint<1>) ~(comp_r_13_15) + comp_r_15_16 + comp_r_15_17 + comp_r_15_18 + comp_r_15_19 + comp_r_15_20 + comp_r_15_21 + comp_r_15_22 + comp_r_15_23 + comp_r_15_24 + comp_r_15_25 + comp_r_15_26 + comp_r_15_27 + comp_r_15_28 + comp_r_15_29 + comp_r_15_30;
    position[16] = 8 + (sc_uint<1>) ~(comp_r_1_16) + (sc_uint<1>) ~(comp_r_3_16) + (sc_uint<1>) ~(comp_r_5_16) + (sc_uint<1>) ~(comp_r_7_16) + (sc_uint<1>) ~(comp_r_9_16) + (sc_uint<1>) ~(comp_r_11_16) + (sc_uint<1>) ~(comp_r_13_16) + (sc_uint<1>) ~(comp_r_15_16);
    position[17] = 9 + (sc_uint<1>) ~(comp_r_1_17) + (sc_uint<1>) ~(comp_r_3_17) + (sc_uint<1>) ~(comp_r_5_17) + (sc_uint<1>) ~(comp_r_7_17) + (sc_uint<1>) ~(comp_r_9_17) + (sc_uint<1>) ~(comp_r_11_17) + (sc_uint<1>) ~(comp_r_13_17) + (sc_uint<1>) ~(comp_r_15_17) + comp_r_17_18 + comp_r_17_19 + comp_r_17_20 + comp_r_17_21 + comp_r_17_22 + comp_r_17_23 + comp_r_17_24 + comp_r_17_25 + comp_r_17_26 + comp_r_17_27 + comp_r_17_28 + comp_r_17_29 + comp_r_17_30;
    position[18] = 9 + (sc_uint<1>) ~(comp_r_1_18) + (sc_uint<1>) ~(comp_r_3_18) + (sc_uint<1>) ~(comp_r_5_18) + (sc_uint<1>) ~(comp_r_7_18) + (sc_uint<1>) ~(comp_r_9_18) + (sc_uint<1>) ~(comp_r_11_18) + (sc_uint<1>) ~(comp_r_13_18) + (sc_uint<1>) ~(comp_r_15_18) + (sc_uint<1>) ~(comp_r_17_18);
    position[19] = 10 + (sc_uint<1>) ~(comp_r_1_19) + (sc_uint<1>) ~(comp_r_3_19) + (sc_uint<1>) ~(comp_r_5_19) + (sc_uint<1>) ~(comp_r_7_19) + (sc_uint<1>) ~(comp_r_9_19) + (sc_uint<1>) ~(comp_r_11_19) + (sc_uint<1>) ~(comp_r_13_19) + (sc_uint<1>) ~(comp_r_15_19) + (sc_uint<1>) ~(comp_r_17_19) + comp_r_19_20 + comp_r_19_21 + comp_r_19_22 + comp_r_19_23 + comp_r_19_24 + comp_r_19_25 + comp_r_19_26 + comp_r_19_27 + comp_r_19_28 + comp_r_19_29 + comp_r_19_30;
    position[20] = 10 + (sc_uint<1>) ~(comp_r_1_20) + (sc_uint<1>) ~(comp_r_3_20) + (sc_uint<1>) ~(comp_r_5_20) + (sc_uint<1>) ~(comp_r_7_20) + (sc_uint<1>) ~(comp_r_9_20) + (sc_uint<1>) ~(comp_r_11_20) + (sc_uint<1>) ~(comp_r_13_20) + (sc_uint<1>) ~(comp_r_15_20) + (sc_uint<1>) ~(comp_r_17_20) + (sc_uint<1>) ~(comp_r_19_20);
    position[21] = 11 + (sc_uint<1>) ~(comp_r_1_21) + (sc_uint<1>) ~(comp_r_3_21) + (sc_uint<1>) ~(comp_r_5_21) + (sc_uint<1>) ~(comp_r_7_21) + (sc_uint<1>) ~(comp_r_9_21) + (sc_uint<1>) ~(comp_r_11_21) + (sc_uint<1>) ~(comp_r_13_21) + (sc_uint<1>) ~(comp_r_15_21) + (sc_uint<1>) ~(comp_r_17_21) + (sc_uint<1>) ~(comp_r_19_21) + comp_r_21_22 + comp_r_21_23 + comp_r_21_24 + comp_r_21_25 + comp_r_21_26 + comp_r_21_27 + comp_r_21_28 + comp_r_21_29 + comp_r_21_30;
    position[22] = 11 + (sc_uint<1>) ~(comp_r_1_22) + (sc_uint<1>) ~(comp_r_3_22) + (sc_uint<1>) ~(comp_r_5_22) + (sc_uint<1>) ~(comp_r_7_22) + (sc_uint<1>) ~(comp_r_9_22) + (sc_uint<1>) ~(comp_r_11_22) + (sc_uint<1>) ~(comp_r_13_22) + (sc_uint<1>) ~(comp_r_15_22) + (sc_uint<1>) ~(comp_r_17_22) + (sc_uint<1>) ~(comp_r_19_22) + (sc_uint<1>) ~(comp_r_21_22);
    position[23] = 12 + (sc_uint<1>) ~(comp_r_1_23) + (sc_uint<1>) ~(comp_r_3_23) + (sc_uint<1>) ~(comp_r_5_23) + (sc_uint<1>) ~(comp_r_7_23) + (sc_uint<1>) ~(comp_r_9_23) + (sc_uint<1>) ~(comp_r_11_23) + (sc_uint<1>) ~(comp_r_13_23) + (sc_uint<1>) ~(comp_r_15_23) + (sc_uint<1>) ~(comp_r_17_23) + (sc_uint<1>) ~(comp_r_19_23) + (sc_uint<1>) ~(comp_r_21_23) + comp_r_23_24 + comp_r_23_25 + comp_r_23_26 + comp_r_23_27 + comp_r_23_28 + comp_r_23_29 + comp_r_23_30;
    position[24] = 12 + (sc_uint<1>) ~(comp_r_1_24) + (sc_uint<1>) ~(comp_r_3_24) + (sc_uint<1>) ~(comp_r_5_24) + (sc_uint<1>) ~(comp_r_7_24) + (sc_uint<1>) ~(comp_r_9_24) + (sc_uint<1>) ~(comp_r_11_24) + (sc_uint<1>) ~(comp_r_13_24) + (sc_uint<1>) ~(comp_r_15_24) + (sc_uint<1>) ~(comp_r_17_24) + (sc_uint<1>) ~(comp_r_19_24) + (sc_uint<1>) ~(comp_r_21_24) + (sc_uint<1>) ~(comp_r_23_24);
    position[25] = 13 + (sc_uint<1>) ~(comp_r_1_25) + (sc_uint<1>) ~(comp_r_3_25) + (sc_uint<1>) ~(comp_r_5_25) + (sc_uint<1>) ~(comp_r_7_25) + (sc_uint<1>) ~(comp_r_9_25) + (sc_uint<1>) ~(comp_r_11_25) + (sc_uint<1>) ~(comp_r_13_25) + (sc_uint<1>) ~(comp_r_15_25) + (sc_uint<1>) ~(comp_r_17_25) + (sc_uint<1>) ~(comp_r_19_25) + (sc_uint<1>) ~(comp_r_21_25) + (sc_uint<1>) ~(comp_r_23_25) + comp_r_25_26 + comp_r_25_27 + comp_r_25_28 + comp_r_25_29 + comp_r_25_30;
    position[26] = 13 + (sc_uint<1>) ~(comp_r_1_26) + (sc_uint<1>) ~(comp_r_3_26) + (sc_uint<1>) ~(comp_r_5_26) + (sc_uint<1>) ~(comp_r_7_26) + (sc_uint<1>) ~(comp_r_9_26) + (sc_uint<1>) ~(comp_r_11_26) + (sc_uint<1>) ~(comp_r_13_26) + (sc_uint<1>) ~(comp_r_15_26) + (sc_uint<1>) ~(comp_r_17_26) + (sc_uint<1>) ~(comp_r_19_26) + (sc_uint<1>) ~(comp_r_21_26) + (sc_uint<1>) ~(comp_r_23_26) + (sc_uint<1>) ~(comp_r_25_26);
    position[27] = 14 + (sc_uint<1>) ~(comp_r_1_27) + (sc_uint<1>) ~(comp_r_3_27) + (sc_uint<1>) ~(comp_r_5_27) + (sc_uint<1>) ~(comp_r_7_27) + (sc_uint<1>) ~(comp_r_9_27) + (sc_uint<1>) ~(comp_r_11_27) + (sc_uint<1>) ~(comp_r_13_27) + (sc_uint<1>) ~(comp_r_15_27) + (sc_uint<1>) ~(comp_r_17_27) + (sc_uint<1>) ~(comp_r_19_27) + (sc_uint<1>) ~(comp_r_21_27) + (sc_uint<1>) ~(comp_r_23_27) + (sc_uint<1>) ~(comp_r_25_27) + comp_r_27_28 + comp_r_27_29 + comp_r_27_30;
    position[28] = 14 + (sc_uint<1>) ~(comp_r_1_28) + (sc_uint<1>) ~(comp_r_3_28) + (sc_uint<1>) ~(comp_r_5_28) + (sc_uint<1>) ~(comp_r_7_28) + (sc_uint<1>) ~(comp_r_9_28) + (sc_uint<1>) ~(comp_r_11_28) + (sc_uint<1>) ~(comp_r_13_28) + (sc_uint<1>) ~(comp_r_15_28) + (sc_uint<1>) ~(comp_r_17_28) + (sc_uint<1>) ~(comp_r_19_28) + (sc_uint<1>) ~(comp_r_21_28) + (sc_uint<1>) ~(comp_r_23_28) + (sc_uint<1>) ~(comp_r_25_28) + (sc_uint<1>) ~(comp_r_27_28);
    position[29] = 15 + (sc_uint<1>) ~(comp_r_1_29) + (sc_uint<1>) ~(comp_r_3_29) + (sc_uint<1>) ~(comp_r_5_29) + (sc_uint<1>) ~(comp_r_7_29) + (sc_uint<1>) ~(comp_r_9_29) + (sc_uint<1>) ~(comp_r_11_29) + (sc_uint<1>) ~(comp_r_13_29) + (sc_uint<1>) ~(comp_r_15_29) + (sc_uint<1>) ~(comp_r_17_29) + (sc_uint<1>) ~(comp_r_19_29) + (sc_uint<1>) ~(comp_r_21_29) + (sc_uint<1>) ~(comp_r_23_29) + (sc_uint<1>) ~(comp_r_25_29) + (sc_uint<1>) ~(comp_r_27_29) + comp_r_29_30;
    position[30] = 15 + (sc_uint<1>) ~(comp_r_1_30) + (sc_uint<1>) ~(comp_r_3_30) + (sc_uint<1>) ~(comp_r_5_30) + (sc_uint<1>) ~(comp_r_7_30) + (sc_uint<1>) ~(comp_r_9_30) + (sc_uint<1>) ~(comp_r_11_30) + (sc_uint<1>) ~(comp_r_13_30) + (sc_uint<1>) ~(comp_r_15_30) + (sc_uint<1>) ~(comp_r_17_30) + (sc_uint<1>) ~(comp_r_19_30) + (sc_uint<1>) ~(comp_r_21_30) + (sc_uint<1>) ~(comp_r_23_30) + (sc_uint<1>) ~(comp_r_25_30) + (sc_uint<1>) ~(comp_r_27_30) + (sc_uint<1>) ~(comp_r_29_30);
    position[31] = (sc_uint<5>) 31;

 // Multiplex the inputs in the order
    PS_struct<1,Q,4> temp[31];
#pragma HLS ARRAY_PARTITION variable=temp complete dim=1
    temp[0] = input[0];
    temp[1] = RO_MUX32 < PS_struct<1,Q,4> > (input, position, (sc_uint<5>) 1);
    temp[2] = RO_MUX32 < PS_struct<1,Q,4> > (input, position, (sc_uint<5>) 2);
    temp[3] = RO_MUX32 < PS_struct<1,Q,4> > (input, position, (sc_uint<5>) 3);
    temp[4] = RO_MUX32 < PS_struct<1,Q,4> > (input, position, (sc_uint<5>) 4);
    temp[5] = RO_MUX32 < PS_struct<1,Q,4> > (input, position, (sc_uint<5>) 5);
    temp[6] = RO_MUX32 < PS_struct<1,Q,4> > (input, position, (sc_uint<5>) 6);
    temp[7] = RO_MUX32 < PS_struct<1,Q,4> > (input, position, (sc_uint<5>) 7);
    temp[8] = RO_MUX32 < PS_struct<1,Q,4> > (input, position, (sc_uint<5>) 8);
    temp[9] = RO_MUX32 < PS_struct<1,Q,4> > (input, position, (sc_uint<5>) 9);
    temp[10] = RO_MUX32 < PS_struct<1,Q,4> > (input, position, (sc_uint<5>) 10);
    temp[11] = RO_MUX32 < PS_struct<1,Q,4> > (input, position, (sc_uint<5>) 11);
    temp[12] = RO_MUX32 < PS_struct<1,Q,4> > (input, position, (sc_uint<5>) 12);
    temp[13] = RO_MUX32 < PS_struct<1,Q,4> > (input, position, (sc_uint<5>) 13);
    temp[14] = RO_MUX32 < PS_struct<1,Q,4> > (input, position, (sc_uint<5>) 14);
    temp[15] = RO_MUX32 < PS_struct<1,Q,4> > (input, position, (sc_uint<5>) 15);
    temp[16] = RO_MUX32 < PS_struct<1,Q,4> > (input, position, (sc_uint<5>) 16);
    temp[17] = RO_MUX32 < PS_struct<1,Q,4> > (input, position, (sc_uint<5>) 17);
    temp[18] = RO_MUX32 < PS_struct<1,Q,4> > (input, position, (sc_uint<5>) 18);
    temp[19] = RO_MUX32 < PS_struct<1,Q,4> > (input, position, (sc_uint<5>) 19);
    temp[20] = RO_MUX32 < PS_struct<1,Q,4> > (input, position, (sc_uint<5>) 20);
    temp[21] = RO_MUX32 < PS_struct<1,Q,4> > (input, position, (sc_uint<5>) 21);
    temp[22] = RO_MUX32 < PS_struct<1,Q,4> > (input, position, (sc_uint<5>) 22);
    temp[23] = RO_MUX32 < PS_struct<1,Q,4> > (input, position, (sc_uint<5>) 23);
    temp[24] = RO_MUX32 < PS_struct<1,Q,4> > (input, position, (sc_uint<5>) 24);
    temp[25] = RO_MUX32 < PS_struct<1,Q,4> > (input, position, (sc_uint<5>) 25);
    temp[26] = RO_MUX32 < PS_struct<1,Q,4> > (input, position, (sc_uint<5>) 26);
    temp[27] = RO_MUX32 < PS_struct<1,Q,4> > (input, position, (sc_uint<5>) 27);
    temp[28] = RO_MUX32 < PS_struct<1,Q,4> > (input, position, (sc_uint<5>) 28);
    temp[29] = RO_MUX32 < PS_struct<1,Q,4> > (input, position, (sc_uint<5>) 29);
    temp[30] = RO_MUX32 < PS_struct<1,Q,4> > (input, position, (sc_uint<5>) 30);

// output the L smallest input
    if (fb == 0){
    output[0] = temp[15];
    output[1] = temp[16];
    output[2] = temp[17];
    output[3] = temp[18];
    output[4] = temp[19];
    output[5] = temp[20];
    output[6] = temp[21];
    output[7] = temp[22];
    output[8] = temp[23];
    output[9] = temp[24];
    output[10] = temp[25];
    output[11] = temp[26];
    output[12] = temp[27];
    output[13] = temp[28];
    output[14] = temp[29];
    output[15] = temp[30];
    }
    else{
    output[0] = temp[0];
    output[1] = temp[1];
    output[2] = temp[2];
    output[3] = temp[3];
    output[4] = temp[4];
    output[5] = temp[5];
    output[6] = temp[6];
    output[7] = temp[7];
    output[8] = temp[8];
    output[9] = temp[9];
    output[10] = temp[10];
    output[11] = temp[11];
    output[12] = temp[12];
    output[13] = temp[13];
    output[14] = temp[14];
    output[15] = temp[15];
    }
}

template <int Q> 
void RANKORDER_SORT_L32 (PS_struct<1,Q,5> input[64], PS_struct<1,Q,5> output[32], sc_biguint<1> fb)
{
#pragma HLS INLINE off
#pragma HLS ARRAY_PARTITION variable=input complete dim=1
#pragma HLS ARRAY_PARTITION variable=output complete dim=1

// Adjust input when FB=0

    if (fb == 0){
        input[63].metric = MAX_VAL;
        input[62] = input[31];
        input[61] = input[30];
        input[60].metric = 0;
        input[59] = input[29];
        input[58].metric = 0;
        input[57] = input[28];
        input[56].metric = 0;
        input[55] = input[27];
        input[54].metric = 0;
        input[53] = input[26];
        input[52].metric = 0;
        input[51] = input[25];
        input[50].metric = 0;
        input[49] = input[24];
        input[48].metric = 0;
        input[47] = input[23];
        input[46].metric = 0;
        input[45] = input[22];
        input[44].metric = 0;
        input[43] = input[21];
        input[42].metric = 0;
        input[41] = input[20];
        input[40].metric = 0;
        input[39] = input[19];
        input[38].metric = 0;
        input[37] = input[18];
        input[36].metric = 0;
        input[35] = input[17];
        input[34].metric = 0;
        input[33] = input[16];
        input[32].metric = 0;
        input[31] = input[15];
        input[30].metric = 0;
        input[29] = input[14];
        input[28].metric = 0;
        input[27] = input[13];
        input[26].metric = 0;
        input[25] = input[12];
        input[24].metric = 0;
        input[23] = input[11];
        input[22].metric = 0;
        input[21] = input[10];
        input[20].metric = 0;
        input[19] = input[9];
        input[18].metric = 0;
        input[17] = input[8];
        input[16].metric = 0;
        input[15] = input[7];
        input[14].metric = 0;
        input[13] = input[6];
        input[12].metric = 0;
        input[11] = input[5];
        input[10].metric = 0;
        input[9] = input[4];
        input[8].metric = 0;
        input[7] = input[3];
        input[6].metric = 0;
        input[5] = input[2];
        input[4].metric = 0;
        input[3] = input[1];
        input[2].metric = 0;
        input[1] = input[0];
    }

    sc_uint<6> position[64];
#pragma HLS ARRAY_PARTITION variable=position complete dim=1

// COMPARE inputs

    sc_uint<1> comp_r_1_2 = COMP < Q > (input[1].metric, input[2].metric);
    sc_uint<1> comp_r_1_3 = COMP < Q > (input[1].metric, input[3].metric);
    sc_uint<1> comp_r_1_4 = COMP < Q > (input[1].metric, input[4].metric);
    sc_uint<1> comp_r_1_5 = COMP < Q > (input[1].metric, input[5].metric);
    sc_uint<1> comp_r_1_6 = COMP < Q > (input[1].metric, input[6].metric);
    sc_uint<1> comp_r_1_7 = COMP < Q > (input[1].metric, input[7].metric);
    sc_uint<1> comp_r_1_8 = COMP < Q > (input[1].metric, input[8].metric);
    sc_uint<1> comp_r_1_9 = COMP < Q > (input[1].metric, input[9].metric);
    sc_uint<1> comp_r_1_10 = COMP < Q > (input[1].metric, input[10].metric);
    sc_uint<1> comp_r_1_11 = COMP < Q > (input[1].metric, input[11].metric);
    sc_uint<1> comp_r_1_12 = COMP < Q > (input[1].metric, input[12].metric);
    sc_uint<1> comp_r_1_13 = COMP < Q > (input[1].metric, input[13].metric);
    sc_uint<1> comp_r_1_14 = COMP < Q > (input[1].metric, input[14].metric);
    sc_uint<1> comp_r_1_15 = COMP < Q > (input[1].metric, input[15].metric);
    sc_uint<1> comp_r_1_16 = COMP < Q > (input[1].metric, input[16].metric);
    sc_uint<1> comp_r_1_17 = COMP < Q > (input[1].metric, input[17].metric);
    sc_uint<1> comp_r_1_18 = COMP < Q > (input[1].metric, input[18].metric);
    sc_uint<1> comp_r_1_19 = COMP < Q > (input[1].metric, input[19].metric);
    sc_uint<1> comp_r_1_20 = COMP < Q > (input[1].metric, input[20].metric);
    sc_uint<1> comp_r_1_21 = COMP < Q > (input[1].metric, input[21].metric);
    sc_uint<1> comp_r_1_22 = COMP < Q > (input[1].metric, input[22].metric);
    sc_uint<1> comp_r_1_23 = COMP < Q > (input[1].metric, input[23].metric);
    sc_uint<1> comp_r_1_24 = COMP < Q > (input[1].metric, input[24].metric);
    sc_uint<1> comp_r_1_25 = COMP < Q > (input[1].metric, input[25].metric);
    sc_uint<1> comp_r_1_26 = COMP < Q > (input[1].metric, input[26].metric);
    sc_uint<1> comp_r_1_27 = COMP < Q > (input[1].metric, input[27].metric);
    sc_uint<1> comp_r_1_28 = COMP < Q > (input[1].metric, input[28].metric);
    sc_uint<1> comp_r_1_29 = COMP < Q > (input[1].metric, input[29].metric);
    sc_uint<1> comp_r_1_30 = COMP < Q > (input[1].metric, input[30].metric);
    sc_uint<1> comp_r_1_31 = COMP < Q > (input[1].metric, input[31].metric);
    sc_uint<1> comp_r_1_32 = COMP < Q > (input[1].metric, input[32].metric);
    sc_uint<1> comp_r_1_33 = COMP < Q > (input[1].metric, input[33].metric);
    sc_uint<1> comp_r_1_34 = COMP < Q > (input[1].metric, input[34].metric);
    sc_uint<1> comp_r_1_35 = COMP < Q > (input[1].metric, input[35].metric);
    sc_uint<1> comp_r_1_36 = COMP < Q > (input[1].metric, input[36].metric);
    sc_uint<1> comp_r_1_37 = COMP < Q > (input[1].metric, input[37].metric);
    sc_uint<1> comp_r_1_38 = COMP < Q > (input[1].metric, input[38].metric);
    sc_uint<1> comp_r_1_39 = COMP < Q > (input[1].metric, input[39].metric);
    sc_uint<1> comp_r_1_40 = COMP < Q > (input[1].metric, input[40].metric);
    sc_uint<1> comp_r_1_41 = COMP < Q > (input[1].metric, input[41].metric);
    sc_uint<1> comp_r_1_42 = COMP < Q > (input[1].metric, input[42].metric);
    sc_uint<1> comp_r_1_43 = COMP < Q > (input[1].metric, input[43].metric);
    sc_uint<1> comp_r_1_44 = COMP < Q > (input[1].metric, input[44].metric);
    sc_uint<1> comp_r_1_45 = COMP < Q > (input[1].metric, input[45].metric);
    sc_uint<1> comp_r_1_46 = COMP < Q > (input[1].metric, input[46].metric);
    sc_uint<1> comp_r_1_47 = COMP < Q > (input[1].metric, input[47].metric);
    sc_uint<1> comp_r_1_48 = COMP < Q > (input[1].metric, input[48].metric);
    sc_uint<1> comp_r_1_49 = COMP < Q > (input[1].metric, input[49].metric);
    sc_uint<1> comp_r_1_50 = COMP < Q > (input[1].metric, input[50].metric);
    sc_uint<1> comp_r_1_51 = COMP < Q > (input[1].metric, input[51].metric);
    sc_uint<1> comp_r_1_52 = COMP < Q > (input[1].metric, input[52].metric);
    sc_uint<1> comp_r_1_53 = COMP < Q > (input[1].metric, input[53].metric);
    sc_uint<1> comp_r_1_54 = COMP < Q > (input[1].metric, input[54].metric);
    sc_uint<1> comp_r_1_55 = COMP < Q > (input[1].metric, input[55].metric);
    sc_uint<1> comp_r_1_56 = COMP < Q > (input[1].metric, input[56].metric);
    sc_uint<1> comp_r_1_57 = COMP < Q > (input[1].metric, input[57].metric);
    sc_uint<1> comp_r_1_58 = COMP < Q > (input[1].metric, input[58].metric);
    sc_uint<1> comp_r_1_59 = COMP < Q > (input[1].metric, input[59].metric);
    sc_uint<1> comp_r_1_60 = COMP < Q > (input[1].metric, input[60].metric);
    sc_uint<1> comp_r_1_61 = COMP < Q > (input[1].metric, input[61].metric);
    sc_uint<1> comp_r_1_62 = COMP < Q > (input[1].metric, input[62].metric);

    sc_uint<1> comp_r_3_4 = COMP < Q > (input[3].metric, input[4].metric);
    sc_uint<1> comp_r_3_5 = COMP < Q > (input[3].metric, input[5].metric);
    sc_uint<1> comp_r_3_6 = COMP < Q > (input[3].metric, input[6].metric);
    sc_uint<1> comp_r_3_7 = COMP < Q > (input[3].metric, input[7].metric);
    sc_uint<1> comp_r_3_8 = COMP < Q > (input[3].metric, input[8].metric);
    sc_uint<1> comp_r_3_9 = COMP < Q > (input[3].metric, input[9].metric);
    sc_uint<1> comp_r_3_10 = COMP < Q > (input[3].metric, input[10].metric);
    sc_uint<1> comp_r_3_11 = COMP < Q > (input[3].metric, input[11].metric);
    sc_uint<1> comp_r_3_12 = COMP < Q > (input[3].metric, input[12].metric);
    sc_uint<1> comp_r_3_13 = COMP < Q > (input[3].metric, input[13].metric);
    sc_uint<1> comp_r_3_14 = COMP < Q > (input[3].metric, input[14].metric);
    sc_uint<1> comp_r_3_15 = COMP < Q > (input[3].metric, input[15].metric);
    sc_uint<1> comp_r_3_16 = COMP < Q > (input[3].metric, input[16].metric);
    sc_uint<1> comp_r_3_17 = COMP < Q > (input[3].metric, input[17].metric);
    sc_uint<1> comp_r_3_18 = COMP < Q > (input[3].metric, input[18].metric);
    sc_uint<1> comp_r_3_19 = COMP < Q > (input[3].metric, input[19].metric);
    sc_uint<1> comp_r_3_20 = COMP < Q > (input[3].metric, input[20].metric);
    sc_uint<1> comp_r_3_21 = COMP < Q > (input[3].metric, input[21].metric);
    sc_uint<1> comp_r_3_22 = COMP < Q > (input[3].metric, input[22].metric);
    sc_uint<1> comp_r_3_23 = COMP < Q > (input[3].metric, input[23].metric);
    sc_uint<1> comp_r_3_24 = COMP < Q > (input[3].metric, input[24].metric);
    sc_uint<1> comp_r_3_25 = COMP < Q > (input[3].metric, input[25].metric);
    sc_uint<1> comp_r_3_26 = COMP < Q > (input[3].metric, input[26].metric);
    sc_uint<1> comp_r_3_27 = COMP < Q > (input[3].metric, input[27].metric);
    sc_uint<1> comp_r_3_28 = COMP < Q > (input[3].metric, input[28].metric);
    sc_uint<1> comp_r_3_29 = COMP < Q > (input[3].metric, input[29].metric);
    sc_uint<1> comp_r_3_30 = COMP < Q > (input[3].metric, input[30].metric);
    sc_uint<1> comp_r_3_31 = COMP < Q > (input[3].metric, input[31].metric);
    sc_uint<1> comp_r_3_32 = COMP < Q > (input[3].metric, input[32].metric);
    sc_uint<1> comp_r_3_33 = COMP < Q > (input[3].metric, input[33].metric);
    sc_uint<1> comp_r_3_34 = COMP < Q > (input[3].metric, input[34].metric);
    sc_uint<1> comp_r_3_35 = COMP < Q > (input[3].metric, input[35].metric);
    sc_uint<1> comp_r_3_36 = COMP < Q > (input[3].metric, input[36].metric);
    sc_uint<1> comp_r_3_37 = COMP < Q > (input[3].metric, input[37].metric);
    sc_uint<1> comp_r_3_38 = COMP < Q > (input[3].metric, input[38].metric);
    sc_uint<1> comp_r_3_39 = COMP < Q > (input[3].metric, input[39].metric);
    sc_uint<1> comp_r_3_40 = COMP < Q > (input[3].metric, input[40].metric);
    sc_uint<1> comp_r_3_41 = COMP < Q > (input[3].metric, input[41].metric);
    sc_uint<1> comp_r_3_42 = COMP < Q > (input[3].metric, input[42].metric);
    sc_uint<1> comp_r_3_43 = COMP < Q > (input[3].metric, input[43].metric);
    sc_uint<1> comp_r_3_44 = COMP < Q > (input[3].metric, input[44].metric);
    sc_uint<1> comp_r_3_45 = COMP < Q > (input[3].metric, input[45].metric);
    sc_uint<1> comp_r_3_46 = COMP < Q > (input[3].metric, input[46].metric);
    sc_uint<1> comp_r_3_47 = COMP < Q > (input[3].metric, input[47].metric);
    sc_uint<1> comp_r_3_48 = COMP < Q > (input[3].metric, input[48].metric);
    sc_uint<1> comp_r_3_49 = COMP < Q > (input[3].metric, input[49].metric);
    sc_uint<1> comp_r_3_50 = COMP < Q > (input[3].metric, input[50].metric);
    sc_uint<1> comp_r_3_51 = COMP < Q > (input[3].metric, input[51].metric);
    sc_uint<1> comp_r_3_52 = COMP < Q > (input[3].metric, input[52].metric);
    sc_uint<1> comp_r_3_53 = COMP < Q > (input[3].metric, input[53].metric);
    sc_uint<1> comp_r_3_54 = COMP < Q > (input[3].metric, input[54].metric);
    sc_uint<1> comp_r_3_55 = COMP < Q > (input[3].metric, input[55].metric);
    sc_uint<1> comp_r_3_56 = COMP < Q > (input[3].metric, input[56].metric);
    sc_uint<1> comp_r_3_57 = COMP < Q > (input[3].metric, input[57].metric);
    sc_uint<1> comp_r_3_58 = COMP < Q > (input[3].metric, input[58].metric);
    sc_uint<1> comp_r_3_59 = COMP < Q > (input[3].metric, input[59].metric);
    sc_uint<1> comp_r_3_60 = COMP < Q > (input[3].metric, input[60].metric);
    sc_uint<1> comp_r_3_61 = COMP < Q > (input[3].metric, input[61].metric);
    sc_uint<1> comp_r_3_62 = COMP < Q > (input[3].metric, input[62].metric);

    sc_uint<1> comp_r_5_6 = COMP < Q > (input[5].metric, input[6].metric);
    sc_uint<1> comp_r_5_7 = COMP < Q > (input[5].metric, input[7].metric);
    sc_uint<1> comp_r_5_8 = COMP < Q > (input[5].metric, input[8].metric);
    sc_uint<1> comp_r_5_9 = COMP < Q > (input[5].metric, input[9].metric);
    sc_uint<1> comp_r_5_10 = COMP < Q > (input[5].metric, input[10].metric);
    sc_uint<1> comp_r_5_11 = COMP < Q > (input[5].metric, input[11].metric);
    sc_uint<1> comp_r_5_12 = COMP < Q > (input[5].metric, input[12].metric);
    sc_uint<1> comp_r_5_13 = COMP < Q > (input[5].metric, input[13].metric);
    sc_uint<1> comp_r_5_14 = COMP < Q > (input[5].metric, input[14].metric);
    sc_uint<1> comp_r_5_15 = COMP < Q > (input[5].metric, input[15].metric);
    sc_uint<1> comp_r_5_16 = COMP < Q > (input[5].metric, input[16].metric);
    sc_uint<1> comp_r_5_17 = COMP < Q > (input[5].metric, input[17].metric);
    sc_uint<1> comp_r_5_18 = COMP < Q > (input[5].metric, input[18].metric);
    sc_uint<1> comp_r_5_19 = COMP < Q > (input[5].metric, input[19].metric);
    sc_uint<1> comp_r_5_20 = COMP < Q > (input[5].metric, input[20].metric);
    sc_uint<1> comp_r_5_21 = COMP < Q > (input[5].metric, input[21].metric);
    sc_uint<1> comp_r_5_22 = COMP < Q > (input[5].metric, input[22].metric);
    sc_uint<1> comp_r_5_23 = COMP < Q > (input[5].metric, input[23].metric);
    sc_uint<1> comp_r_5_24 = COMP < Q > (input[5].metric, input[24].metric);
    sc_uint<1> comp_r_5_25 = COMP < Q > (input[5].metric, input[25].metric);
    sc_uint<1> comp_r_5_26 = COMP < Q > (input[5].metric, input[26].metric);
    sc_uint<1> comp_r_5_27 = COMP < Q > (input[5].metric, input[27].metric);
    sc_uint<1> comp_r_5_28 = COMP < Q > (input[5].metric, input[28].metric);
    sc_uint<1> comp_r_5_29 = COMP < Q > (input[5].metric, input[29].metric);
    sc_uint<1> comp_r_5_30 = COMP < Q > (input[5].metric, input[30].metric);
    sc_uint<1> comp_r_5_31 = COMP < Q > (input[5].metric, input[31].metric);
    sc_uint<1> comp_r_5_32 = COMP < Q > (input[5].metric, input[32].metric);
    sc_uint<1> comp_r_5_33 = COMP < Q > (input[5].metric, input[33].metric);
    sc_uint<1> comp_r_5_34 = COMP < Q > (input[5].metric, input[34].metric);
    sc_uint<1> comp_r_5_35 = COMP < Q > (input[5].metric, input[35].metric);
    sc_uint<1> comp_r_5_36 = COMP < Q > (input[5].metric, input[36].metric);
    sc_uint<1> comp_r_5_37 = COMP < Q > (input[5].metric, input[37].metric);
    sc_uint<1> comp_r_5_38 = COMP < Q > (input[5].metric, input[38].metric);
    sc_uint<1> comp_r_5_39 = COMP < Q > (input[5].metric, input[39].metric);
    sc_uint<1> comp_r_5_40 = COMP < Q > (input[5].metric, input[40].metric);
    sc_uint<1> comp_r_5_41 = COMP < Q > (input[5].metric, input[41].metric);
    sc_uint<1> comp_r_5_42 = COMP < Q > (input[5].metric, input[42].metric);
    sc_uint<1> comp_r_5_43 = COMP < Q > (input[5].metric, input[43].metric);
    sc_uint<1> comp_r_5_44 = COMP < Q > (input[5].metric, input[44].metric);
    sc_uint<1> comp_r_5_45 = COMP < Q > (input[5].metric, input[45].metric);
    sc_uint<1> comp_r_5_46 = COMP < Q > (input[5].metric, input[46].metric);
    sc_uint<1> comp_r_5_47 = COMP < Q > (input[5].metric, input[47].metric);
    sc_uint<1> comp_r_5_48 = COMP < Q > (input[5].metric, input[48].metric);
    sc_uint<1> comp_r_5_49 = COMP < Q > (input[5].metric, input[49].metric);
    sc_uint<1> comp_r_5_50 = COMP < Q > (input[5].metric, input[50].metric);
    sc_uint<1> comp_r_5_51 = COMP < Q > (input[5].metric, input[51].metric);
    sc_uint<1> comp_r_5_52 = COMP < Q > (input[5].metric, input[52].metric);
    sc_uint<1> comp_r_5_53 = COMP < Q > (input[5].metric, input[53].metric);
    sc_uint<1> comp_r_5_54 = COMP < Q > (input[5].metric, input[54].metric);
    sc_uint<1> comp_r_5_55 = COMP < Q > (input[5].metric, input[55].metric);
    sc_uint<1> comp_r_5_56 = COMP < Q > (input[5].metric, input[56].metric);
    sc_uint<1> comp_r_5_57 = COMP < Q > (input[5].metric, input[57].metric);
    sc_uint<1> comp_r_5_58 = COMP < Q > (input[5].metric, input[58].metric);
    sc_uint<1> comp_r_5_59 = COMP < Q > (input[5].metric, input[59].metric);
    sc_uint<1> comp_r_5_60 = COMP < Q > (input[5].metric, input[60].metric);
    sc_uint<1> comp_r_5_61 = COMP < Q > (input[5].metric, input[61].metric);
    sc_uint<1> comp_r_5_62 = COMP < Q > (input[5].metric, input[62].metric);

    sc_uint<1> comp_r_7_8 = COMP < Q > (input[7].metric, input[8].metric);
    sc_uint<1> comp_r_7_9 = COMP < Q > (input[7].metric, input[9].metric);
    sc_uint<1> comp_r_7_10 = COMP < Q > (input[7].metric, input[10].metric);
    sc_uint<1> comp_r_7_11 = COMP < Q > (input[7].metric, input[11].metric);
    sc_uint<1> comp_r_7_12 = COMP < Q > (input[7].metric, input[12].metric);
    sc_uint<1> comp_r_7_13 = COMP < Q > (input[7].metric, input[13].metric);
    sc_uint<1> comp_r_7_14 = COMP < Q > (input[7].metric, input[14].metric);
    sc_uint<1> comp_r_7_15 = COMP < Q > (input[7].metric, input[15].metric);
    sc_uint<1> comp_r_7_16 = COMP < Q > (input[7].metric, input[16].metric);
    sc_uint<1> comp_r_7_17 = COMP < Q > (input[7].metric, input[17].metric);
    sc_uint<1> comp_r_7_18 = COMP < Q > (input[7].metric, input[18].metric);
    sc_uint<1> comp_r_7_19 = COMP < Q > (input[7].metric, input[19].metric);
    sc_uint<1> comp_r_7_20 = COMP < Q > (input[7].metric, input[20].metric);
    sc_uint<1> comp_r_7_21 = COMP < Q > (input[7].metric, input[21].metric);
    sc_uint<1> comp_r_7_22 = COMP < Q > (input[7].metric, input[22].metric);
    sc_uint<1> comp_r_7_23 = COMP < Q > (input[7].metric, input[23].metric);
    sc_uint<1> comp_r_7_24 = COMP < Q > (input[7].metric, input[24].metric);
    sc_uint<1> comp_r_7_25 = COMP < Q > (input[7].metric, input[25].metric);
    sc_uint<1> comp_r_7_26 = COMP < Q > (input[7].metric, input[26].metric);
    sc_uint<1> comp_r_7_27 = COMP < Q > (input[7].metric, input[27].metric);
    sc_uint<1> comp_r_7_28 = COMP < Q > (input[7].metric, input[28].metric);
    sc_uint<1> comp_r_7_29 = COMP < Q > (input[7].metric, input[29].metric);
    sc_uint<1> comp_r_7_30 = COMP < Q > (input[7].metric, input[30].metric);
    sc_uint<1> comp_r_7_31 = COMP < Q > (input[7].metric, input[31].metric);
    sc_uint<1> comp_r_7_32 = COMP < Q > (input[7].metric, input[32].metric);
    sc_uint<1> comp_r_7_33 = COMP < Q > (input[7].metric, input[33].metric);
    sc_uint<1> comp_r_7_34 = COMP < Q > (input[7].metric, input[34].metric);
    sc_uint<1> comp_r_7_35 = COMP < Q > (input[7].metric, input[35].metric);
    sc_uint<1> comp_r_7_36 = COMP < Q > (input[7].metric, input[36].metric);
    sc_uint<1> comp_r_7_37 = COMP < Q > (input[7].metric, input[37].metric);
    sc_uint<1> comp_r_7_38 = COMP < Q > (input[7].metric, input[38].metric);
    sc_uint<1> comp_r_7_39 = COMP < Q > (input[7].metric, input[39].metric);
    sc_uint<1> comp_r_7_40 = COMP < Q > (input[7].metric, input[40].metric);
    sc_uint<1> comp_r_7_41 = COMP < Q > (input[7].metric, input[41].metric);
    sc_uint<1> comp_r_7_42 = COMP < Q > (input[7].metric, input[42].metric);
    sc_uint<1> comp_r_7_43 = COMP < Q > (input[7].metric, input[43].metric);
    sc_uint<1> comp_r_7_44 = COMP < Q > (input[7].metric, input[44].metric);
    sc_uint<1> comp_r_7_45 = COMP < Q > (input[7].metric, input[45].metric);
    sc_uint<1> comp_r_7_46 = COMP < Q > (input[7].metric, input[46].metric);
    sc_uint<1> comp_r_7_47 = COMP < Q > (input[7].metric, input[47].metric);
    sc_uint<1> comp_r_7_48 = COMP < Q > (input[7].metric, input[48].metric);
    sc_uint<1> comp_r_7_49 = COMP < Q > (input[7].metric, input[49].metric);
    sc_uint<1> comp_r_7_50 = COMP < Q > (input[7].metric, input[50].metric);
    sc_uint<1> comp_r_7_51 = COMP < Q > (input[7].metric, input[51].metric);
    sc_uint<1> comp_r_7_52 = COMP < Q > (input[7].metric, input[52].metric);
    sc_uint<1> comp_r_7_53 = COMP < Q > (input[7].metric, input[53].metric);
    sc_uint<1> comp_r_7_54 = COMP < Q > (input[7].metric, input[54].metric);
    sc_uint<1> comp_r_7_55 = COMP < Q > (input[7].metric, input[55].metric);
    sc_uint<1> comp_r_7_56 = COMP < Q > (input[7].metric, input[56].metric);
    sc_uint<1> comp_r_7_57 = COMP < Q > (input[7].metric, input[57].metric);
    sc_uint<1> comp_r_7_58 = COMP < Q > (input[7].metric, input[58].metric);
    sc_uint<1> comp_r_7_59 = COMP < Q > (input[7].metric, input[59].metric);
    sc_uint<1> comp_r_7_60 = COMP < Q > (input[7].metric, input[60].metric);
    sc_uint<1> comp_r_7_61 = COMP < Q > (input[7].metric, input[61].metric);
    sc_uint<1> comp_r_7_62 = COMP < Q > (input[7].metric, input[62].metric);

    sc_uint<1> comp_r_9_10 = COMP < Q > (input[9].metric, input[10].metric);
    sc_uint<1> comp_r_9_11 = COMP < Q > (input[9].metric, input[11].metric);
    sc_uint<1> comp_r_9_12 = COMP < Q > (input[9].metric, input[12].metric);
    sc_uint<1> comp_r_9_13 = COMP < Q > (input[9].metric, input[13].metric);
    sc_uint<1> comp_r_9_14 = COMP < Q > (input[9].metric, input[14].metric);
    sc_uint<1> comp_r_9_15 = COMP < Q > (input[9].metric, input[15].metric);
    sc_uint<1> comp_r_9_16 = COMP < Q > (input[9].metric, input[16].metric);
    sc_uint<1> comp_r_9_17 = COMP < Q > (input[9].metric, input[17].metric);
    sc_uint<1> comp_r_9_18 = COMP < Q > (input[9].metric, input[18].metric);
    sc_uint<1> comp_r_9_19 = COMP < Q > (input[9].metric, input[19].metric);
    sc_uint<1> comp_r_9_20 = COMP < Q > (input[9].metric, input[20].metric);
    sc_uint<1> comp_r_9_21 = COMP < Q > (input[9].metric, input[21].metric);
    sc_uint<1> comp_r_9_22 = COMP < Q > (input[9].metric, input[22].metric);
    sc_uint<1> comp_r_9_23 = COMP < Q > (input[9].metric, input[23].metric);
    sc_uint<1> comp_r_9_24 = COMP < Q > (input[9].metric, input[24].metric);
    sc_uint<1> comp_r_9_25 = COMP < Q > (input[9].metric, input[25].metric);
    sc_uint<1> comp_r_9_26 = COMP < Q > (input[9].metric, input[26].metric);
    sc_uint<1> comp_r_9_27 = COMP < Q > (input[9].metric, input[27].metric);
    sc_uint<1> comp_r_9_28 = COMP < Q > (input[9].metric, input[28].metric);
    sc_uint<1> comp_r_9_29 = COMP < Q > (input[9].metric, input[29].metric);
    sc_uint<1> comp_r_9_30 = COMP < Q > (input[9].metric, input[30].metric);
    sc_uint<1> comp_r_9_31 = COMP < Q > (input[9].metric, input[31].metric);
    sc_uint<1> comp_r_9_32 = COMP < Q > (input[9].metric, input[32].metric);
    sc_uint<1> comp_r_9_33 = COMP < Q > (input[9].metric, input[33].metric);
    sc_uint<1> comp_r_9_34 = COMP < Q > (input[9].metric, input[34].metric);
    sc_uint<1> comp_r_9_35 = COMP < Q > (input[9].metric, input[35].metric);
    sc_uint<1> comp_r_9_36 = COMP < Q > (input[9].metric, input[36].metric);
    sc_uint<1> comp_r_9_37 = COMP < Q > (input[9].metric, input[37].metric);
    sc_uint<1> comp_r_9_38 = COMP < Q > (input[9].metric, input[38].metric);
    sc_uint<1> comp_r_9_39 = COMP < Q > (input[9].metric, input[39].metric);
    sc_uint<1> comp_r_9_40 = COMP < Q > (input[9].metric, input[40].metric);
    sc_uint<1> comp_r_9_41 = COMP < Q > (input[9].metric, input[41].metric);
    sc_uint<1> comp_r_9_42 = COMP < Q > (input[9].metric, input[42].metric);
    sc_uint<1> comp_r_9_43 = COMP < Q > (input[9].metric, input[43].metric);
    sc_uint<1> comp_r_9_44 = COMP < Q > (input[9].metric, input[44].metric);
    sc_uint<1> comp_r_9_45 = COMP < Q > (input[9].metric, input[45].metric);
    sc_uint<1> comp_r_9_46 = COMP < Q > (input[9].metric, input[46].metric);
    sc_uint<1> comp_r_9_47 = COMP < Q > (input[9].metric, input[47].metric);
    sc_uint<1> comp_r_9_48 = COMP < Q > (input[9].metric, input[48].metric);
    sc_uint<1> comp_r_9_49 = COMP < Q > (input[9].metric, input[49].metric);
    sc_uint<1> comp_r_9_50 = COMP < Q > (input[9].metric, input[50].metric);
    sc_uint<1> comp_r_9_51 = COMP < Q > (input[9].metric, input[51].metric);
    sc_uint<1> comp_r_9_52 = COMP < Q > (input[9].metric, input[52].metric);
    sc_uint<1> comp_r_9_53 = COMP < Q > (input[9].metric, input[53].metric);
    sc_uint<1> comp_r_9_54 = COMP < Q > (input[9].metric, input[54].metric);
    sc_uint<1> comp_r_9_55 = COMP < Q > (input[9].metric, input[55].metric);
    sc_uint<1> comp_r_9_56 = COMP < Q > (input[9].metric, input[56].metric);
    sc_uint<1> comp_r_9_57 = COMP < Q > (input[9].metric, input[57].metric);
    sc_uint<1> comp_r_9_58 = COMP < Q > (input[9].metric, input[58].metric);
    sc_uint<1> comp_r_9_59 = COMP < Q > (input[9].metric, input[59].metric);
    sc_uint<1> comp_r_9_60 = COMP < Q > (input[9].metric, input[60].metric);
    sc_uint<1> comp_r_9_61 = COMP < Q > (input[9].metric, input[61].metric);
    sc_uint<1> comp_r_9_62 = COMP < Q > (input[9].metric, input[62].metric);

    sc_uint<1> comp_r_11_12 = COMP < Q > (input[11].metric, input[12].metric);
    sc_uint<1> comp_r_11_13 = COMP < Q > (input[11].metric, input[13].metric);
    sc_uint<1> comp_r_11_14 = COMP < Q > (input[11].metric, input[14].metric);
    sc_uint<1> comp_r_11_15 = COMP < Q > (input[11].metric, input[15].metric);
    sc_uint<1> comp_r_11_16 = COMP < Q > (input[11].metric, input[16].metric);
    sc_uint<1> comp_r_11_17 = COMP < Q > (input[11].metric, input[17].metric);
    sc_uint<1> comp_r_11_18 = COMP < Q > (input[11].metric, input[18].metric);
    sc_uint<1> comp_r_11_19 = COMP < Q > (input[11].metric, input[19].metric);
    sc_uint<1> comp_r_11_20 = COMP < Q > (input[11].metric, input[20].metric);
    sc_uint<1> comp_r_11_21 = COMP < Q > (input[11].metric, input[21].metric);
    sc_uint<1> comp_r_11_22 = COMP < Q > (input[11].metric, input[22].metric);
    sc_uint<1> comp_r_11_23 = COMP < Q > (input[11].metric, input[23].metric);
    sc_uint<1> comp_r_11_24 = COMP < Q > (input[11].metric, input[24].metric);
    sc_uint<1> comp_r_11_25 = COMP < Q > (input[11].metric, input[25].metric);
    sc_uint<1> comp_r_11_26 = COMP < Q > (input[11].metric, input[26].metric);
    sc_uint<1> comp_r_11_27 = COMP < Q > (input[11].metric, input[27].metric);
    sc_uint<1> comp_r_11_28 = COMP < Q > (input[11].metric, input[28].metric);
    sc_uint<1> comp_r_11_29 = COMP < Q > (input[11].metric, input[29].metric);
    sc_uint<1> comp_r_11_30 = COMP < Q > (input[11].metric, input[30].metric);
    sc_uint<1> comp_r_11_31 = COMP < Q > (input[11].metric, input[31].metric);
    sc_uint<1> comp_r_11_32 = COMP < Q > (input[11].metric, input[32].metric);
    sc_uint<1> comp_r_11_33 = COMP < Q > (input[11].metric, input[33].metric);
    sc_uint<1> comp_r_11_34 = COMP < Q > (input[11].metric, input[34].metric);
    sc_uint<1> comp_r_11_35 = COMP < Q > (input[11].metric, input[35].metric);
    sc_uint<1> comp_r_11_36 = COMP < Q > (input[11].metric, input[36].metric);
    sc_uint<1> comp_r_11_37 = COMP < Q > (input[11].metric, input[37].metric);
    sc_uint<1> comp_r_11_38 = COMP < Q > (input[11].metric, input[38].metric);
    sc_uint<1> comp_r_11_39 = COMP < Q > (input[11].metric, input[39].metric);
    sc_uint<1> comp_r_11_40 = COMP < Q > (input[11].metric, input[40].metric);
    sc_uint<1> comp_r_11_41 = COMP < Q > (input[11].metric, input[41].metric);
    sc_uint<1> comp_r_11_42 = COMP < Q > (input[11].metric, input[42].metric);
    sc_uint<1> comp_r_11_43 = COMP < Q > (input[11].metric, input[43].metric);
    sc_uint<1> comp_r_11_44 = COMP < Q > (input[11].metric, input[44].metric);
    sc_uint<1> comp_r_11_45 = COMP < Q > (input[11].metric, input[45].metric);
    sc_uint<1> comp_r_11_46 = COMP < Q > (input[11].metric, input[46].metric);
    sc_uint<1> comp_r_11_47 = COMP < Q > (input[11].metric, input[47].metric);
    sc_uint<1> comp_r_11_48 = COMP < Q > (input[11].metric, input[48].metric);
    sc_uint<1> comp_r_11_49 = COMP < Q > (input[11].metric, input[49].metric);
    sc_uint<1> comp_r_11_50 = COMP < Q > (input[11].metric, input[50].metric);
    sc_uint<1> comp_r_11_51 = COMP < Q > (input[11].metric, input[51].metric);
    sc_uint<1> comp_r_11_52 = COMP < Q > (input[11].metric, input[52].metric);
    sc_uint<1> comp_r_11_53 = COMP < Q > (input[11].metric, input[53].metric);
    sc_uint<1> comp_r_11_54 = COMP < Q > (input[11].metric, input[54].metric);
    sc_uint<1> comp_r_11_55 = COMP < Q > (input[11].metric, input[55].metric);
    sc_uint<1> comp_r_11_56 = COMP < Q > (input[11].metric, input[56].metric);
    sc_uint<1> comp_r_11_57 = COMP < Q > (input[11].metric, input[57].metric);
    sc_uint<1> comp_r_11_58 = COMP < Q > (input[11].metric, input[58].metric);
    sc_uint<1> comp_r_11_59 = COMP < Q > (input[11].metric, input[59].metric);
    sc_uint<1> comp_r_11_60 = COMP < Q > (input[11].metric, input[60].metric);
    sc_uint<1> comp_r_11_61 = COMP < Q > (input[11].metric, input[61].metric);
    sc_uint<1> comp_r_11_62 = COMP < Q > (input[11].metric, input[62].metric);

    sc_uint<1> comp_r_13_14 = COMP < Q > (input[13].metric, input[14].metric);
    sc_uint<1> comp_r_13_15 = COMP < Q > (input[13].metric, input[15].metric);
    sc_uint<1> comp_r_13_16 = COMP < Q > (input[13].metric, input[16].metric);
    sc_uint<1> comp_r_13_17 = COMP < Q > (input[13].metric, input[17].metric);
    sc_uint<1> comp_r_13_18 = COMP < Q > (input[13].metric, input[18].metric);
    sc_uint<1> comp_r_13_19 = COMP < Q > (input[13].metric, input[19].metric);
    sc_uint<1> comp_r_13_20 = COMP < Q > (input[13].metric, input[20].metric);
    sc_uint<1> comp_r_13_21 = COMP < Q > (input[13].metric, input[21].metric);
    sc_uint<1> comp_r_13_22 = COMP < Q > (input[13].metric, input[22].metric);
    sc_uint<1> comp_r_13_23 = COMP < Q > (input[13].metric, input[23].metric);
    sc_uint<1> comp_r_13_24 = COMP < Q > (input[13].metric, input[24].metric);
    sc_uint<1> comp_r_13_25 = COMP < Q > (input[13].metric, input[25].metric);
    sc_uint<1> comp_r_13_26 = COMP < Q > (input[13].metric, input[26].metric);
    sc_uint<1> comp_r_13_27 = COMP < Q > (input[13].metric, input[27].metric);
    sc_uint<1> comp_r_13_28 = COMP < Q > (input[13].metric, input[28].metric);
    sc_uint<1> comp_r_13_29 = COMP < Q > (input[13].metric, input[29].metric);
    sc_uint<1> comp_r_13_30 = COMP < Q > (input[13].metric, input[30].metric);
    sc_uint<1> comp_r_13_31 = COMP < Q > (input[13].metric, input[31].metric);
    sc_uint<1> comp_r_13_32 = COMP < Q > (input[13].metric, input[32].metric);
    sc_uint<1> comp_r_13_33 = COMP < Q > (input[13].metric, input[33].metric);
    sc_uint<1> comp_r_13_34 = COMP < Q > (input[13].metric, input[34].metric);
    sc_uint<1> comp_r_13_35 = COMP < Q > (input[13].metric, input[35].metric);
    sc_uint<1> comp_r_13_36 = COMP < Q > (input[13].metric, input[36].metric);
    sc_uint<1> comp_r_13_37 = COMP < Q > (input[13].metric, input[37].metric);
    sc_uint<1> comp_r_13_38 = COMP < Q > (input[13].metric, input[38].metric);
    sc_uint<1> comp_r_13_39 = COMP < Q > (input[13].metric, input[39].metric);
    sc_uint<1> comp_r_13_40 = COMP < Q > (input[13].metric, input[40].metric);
    sc_uint<1> comp_r_13_41 = COMP < Q > (input[13].metric, input[41].metric);
    sc_uint<1> comp_r_13_42 = COMP < Q > (input[13].metric, input[42].metric);
    sc_uint<1> comp_r_13_43 = COMP < Q > (input[13].metric, input[43].metric);
    sc_uint<1> comp_r_13_44 = COMP < Q > (input[13].metric, input[44].metric);
    sc_uint<1> comp_r_13_45 = COMP < Q > (input[13].metric, input[45].metric);
    sc_uint<1> comp_r_13_46 = COMP < Q > (input[13].metric, input[46].metric);
    sc_uint<1> comp_r_13_47 = COMP < Q > (input[13].metric, input[47].metric);
    sc_uint<1> comp_r_13_48 = COMP < Q > (input[13].metric, input[48].metric);
    sc_uint<1> comp_r_13_49 = COMP < Q > (input[13].metric, input[49].metric);
    sc_uint<1> comp_r_13_50 = COMP < Q > (input[13].metric, input[50].metric);
    sc_uint<1> comp_r_13_51 = COMP < Q > (input[13].metric, input[51].metric);
    sc_uint<1> comp_r_13_52 = COMP < Q > (input[13].metric, input[52].metric);
    sc_uint<1> comp_r_13_53 = COMP < Q > (input[13].metric, input[53].metric);
    sc_uint<1> comp_r_13_54 = COMP < Q > (input[13].metric, input[54].metric);
    sc_uint<1> comp_r_13_55 = COMP < Q > (input[13].metric, input[55].metric);
    sc_uint<1> comp_r_13_56 = COMP < Q > (input[13].metric, input[56].metric);
    sc_uint<1> comp_r_13_57 = COMP < Q > (input[13].metric, input[57].metric);
    sc_uint<1> comp_r_13_58 = COMP < Q > (input[13].metric, input[58].metric);
    sc_uint<1> comp_r_13_59 = COMP < Q > (input[13].metric, input[59].metric);
    sc_uint<1> comp_r_13_60 = COMP < Q > (input[13].metric, input[60].metric);
    sc_uint<1> comp_r_13_61 = COMP < Q > (input[13].metric, input[61].metric);
    sc_uint<1> comp_r_13_62 = COMP < Q > (input[13].metric, input[62].metric);

    sc_uint<1> comp_r_15_16 = COMP < Q > (input[15].metric, input[16].metric);
    sc_uint<1> comp_r_15_17 = COMP < Q > (input[15].metric, input[17].metric);
    sc_uint<1> comp_r_15_18 = COMP < Q > (input[15].metric, input[18].metric);
    sc_uint<1> comp_r_15_19 = COMP < Q > (input[15].metric, input[19].metric);
    sc_uint<1> comp_r_15_20 = COMP < Q > (input[15].metric, input[20].metric);
    sc_uint<1> comp_r_15_21 = COMP < Q > (input[15].metric, input[21].metric);
    sc_uint<1> comp_r_15_22 = COMP < Q > (input[15].metric, input[22].metric);
    sc_uint<1> comp_r_15_23 = COMP < Q > (input[15].metric, input[23].metric);
    sc_uint<1> comp_r_15_24 = COMP < Q > (input[15].metric, input[24].metric);
    sc_uint<1> comp_r_15_25 = COMP < Q > (input[15].metric, input[25].metric);
    sc_uint<1> comp_r_15_26 = COMP < Q > (input[15].metric, input[26].metric);
    sc_uint<1> comp_r_15_27 = COMP < Q > (input[15].metric, input[27].metric);
    sc_uint<1> comp_r_15_28 = COMP < Q > (input[15].metric, input[28].metric);
    sc_uint<1> comp_r_15_29 = COMP < Q > (input[15].metric, input[29].metric);
    sc_uint<1> comp_r_15_30 = COMP < Q > (input[15].metric, input[30].metric);
    sc_uint<1> comp_r_15_31 = COMP < Q > (input[15].metric, input[31].metric);
    sc_uint<1> comp_r_15_32 = COMP < Q > (input[15].metric, input[32].metric);
    sc_uint<1> comp_r_15_33 = COMP < Q > (input[15].metric, input[33].metric);
    sc_uint<1> comp_r_15_34 = COMP < Q > (input[15].metric, input[34].metric);
    sc_uint<1> comp_r_15_35 = COMP < Q > (input[15].metric, input[35].metric);
    sc_uint<1> comp_r_15_36 = COMP < Q > (input[15].metric, input[36].metric);
    sc_uint<1> comp_r_15_37 = COMP < Q > (input[15].metric, input[37].metric);
    sc_uint<1> comp_r_15_38 = COMP < Q > (input[15].metric, input[38].metric);
    sc_uint<1> comp_r_15_39 = COMP < Q > (input[15].metric, input[39].metric);
    sc_uint<1> comp_r_15_40 = COMP < Q > (input[15].metric, input[40].metric);
    sc_uint<1> comp_r_15_41 = COMP < Q > (input[15].metric, input[41].metric);
    sc_uint<1> comp_r_15_42 = COMP < Q > (input[15].metric, input[42].metric);
    sc_uint<1> comp_r_15_43 = COMP < Q > (input[15].metric, input[43].metric);
    sc_uint<1> comp_r_15_44 = COMP < Q > (input[15].metric, input[44].metric);
    sc_uint<1> comp_r_15_45 = COMP < Q > (input[15].metric, input[45].metric);
    sc_uint<1> comp_r_15_46 = COMP < Q > (input[15].metric, input[46].metric);
    sc_uint<1> comp_r_15_47 = COMP < Q > (input[15].metric, input[47].metric);
    sc_uint<1> comp_r_15_48 = COMP < Q > (input[15].metric, input[48].metric);
    sc_uint<1> comp_r_15_49 = COMP < Q > (input[15].metric, input[49].metric);
    sc_uint<1> comp_r_15_50 = COMP < Q > (input[15].metric, input[50].metric);
    sc_uint<1> comp_r_15_51 = COMP < Q > (input[15].metric, input[51].metric);
    sc_uint<1> comp_r_15_52 = COMP < Q > (input[15].metric, input[52].metric);
    sc_uint<1> comp_r_15_53 = COMP < Q > (input[15].metric, input[53].metric);
    sc_uint<1> comp_r_15_54 = COMP < Q > (input[15].metric, input[54].metric);
    sc_uint<1> comp_r_15_55 = COMP < Q > (input[15].metric, input[55].metric);
    sc_uint<1> comp_r_15_56 = COMP < Q > (input[15].metric, input[56].metric);
    sc_uint<1> comp_r_15_57 = COMP < Q > (input[15].metric, input[57].metric);
    sc_uint<1> comp_r_15_58 = COMP < Q > (input[15].metric, input[58].metric);
    sc_uint<1> comp_r_15_59 = COMP < Q > (input[15].metric, input[59].metric);
    sc_uint<1> comp_r_15_60 = COMP < Q > (input[15].metric, input[60].metric);
    sc_uint<1> comp_r_15_61 = COMP < Q > (input[15].metric, input[61].metric);
    sc_uint<1> comp_r_15_62 = COMP < Q > (input[15].metric, input[62].metric);

    sc_uint<1> comp_r_17_18 = COMP < Q > (input[17].metric, input[18].metric);
    sc_uint<1> comp_r_17_19 = COMP < Q > (input[17].metric, input[19].metric);
    sc_uint<1> comp_r_17_20 = COMP < Q > (input[17].metric, input[20].metric);
    sc_uint<1> comp_r_17_21 = COMP < Q > (input[17].metric, input[21].metric);
    sc_uint<1> comp_r_17_22 = COMP < Q > (input[17].metric, input[22].metric);
    sc_uint<1> comp_r_17_23 = COMP < Q > (input[17].metric, input[23].metric);
    sc_uint<1> comp_r_17_24 = COMP < Q > (input[17].metric, input[24].metric);
    sc_uint<1> comp_r_17_25 = COMP < Q > (input[17].metric, input[25].metric);
    sc_uint<1> comp_r_17_26 = COMP < Q > (input[17].metric, input[26].metric);
    sc_uint<1> comp_r_17_27 = COMP < Q > (input[17].metric, input[27].metric);
    sc_uint<1> comp_r_17_28 = COMP < Q > (input[17].metric, input[28].metric);
    sc_uint<1> comp_r_17_29 = COMP < Q > (input[17].metric, input[29].metric);
    sc_uint<1> comp_r_17_30 = COMP < Q > (input[17].metric, input[30].metric);
    sc_uint<1> comp_r_17_31 = COMP < Q > (input[17].metric, input[31].metric);
    sc_uint<1> comp_r_17_32 = COMP < Q > (input[17].metric, input[32].metric);
    sc_uint<1> comp_r_17_33 = COMP < Q > (input[17].metric, input[33].metric);
    sc_uint<1> comp_r_17_34 = COMP < Q > (input[17].metric, input[34].metric);
    sc_uint<1> comp_r_17_35 = COMP < Q > (input[17].metric, input[35].metric);
    sc_uint<1> comp_r_17_36 = COMP < Q > (input[17].metric, input[36].metric);
    sc_uint<1> comp_r_17_37 = COMP < Q > (input[17].metric, input[37].metric);
    sc_uint<1> comp_r_17_38 = COMP < Q > (input[17].metric, input[38].metric);
    sc_uint<1> comp_r_17_39 = COMP < Q > (input[17].metric, input[39].metric);
    sc_uint<1> comp_r_17_40 = COMP < Q > (input[17].metric, input[40].metric);
    sc_uint<1> comp_r_17_41 = COMP < Q > (input[17].metric, input[41].metric);
    sc_uint<1> comp_r_17_42 = COMP < Q > (input[17].metric, input[42].metric);
    sc_uint<1> comp_r_17_43 = COMP < Q > (input[17].metric, input[43].metric);
    sc_uint<1> comp_r_17_44 = COMP < Q > (input[17].metric, input[44].metric);
    sc_uint<1> comp_r_17_45 = COMP < Q > (input[17].metric, input[45].metric);
    sc_uint<1> comp_r_17_46 = COMP < Q > (input[17].metric, input[46].metric);
    sc_uint<1> comp_r_17_47 = COMP < Q > (input[17].metric, input[47].metric);
    sc_uint<1> comp_r_17_48 = COMP < Q > (input[17].metric, input[48].metric);
    sc_uint<1> comp_r_17_49 = COMP < Q > (input[17].metric, input[49].metric);
    sc_uint<1> comp_r_17_50 = COMP < Q > (input[17].metric, input[50].metric);
    sc_uint<1> comp_r_17_51 = COMP < Q > (input[17].metric, input[51].metric);
    sc_uint<1> comp_r_17_52 = COMP < Q > (input[17].metric, input[52].metric);
    sc_uint<1> comp_r_17_53 = COMP < Q > (input[17].metric, input[53].metric);
    sc_uint<1> comp_r_17_54 = COMP < Q > (input[17].metric, input[54].metric);
    sc_uint<1> comp_r_17_55 = COMP < Q > (input[17].metric, input[55].metric);
    sc_uint<1> comp_r_17_56 = COMP < Q > (input[17].metric, input[56].metric);
    sc_uint<1> comp_r_17_57 = COMP < Q > (input[17].metric, input[57].metric);
    sc_uint<1> comp_r_17_58 = COMP < Q > (input[17].metric, input[58].metric);
    sc_uint<1> comp_r_17_59 = COMP < Q > (input[17].metric, input[59].metric);
    sc_uint<1> comp_r_17_60 = COMP < Q > (input[17].metric, input[60].metric);
    sc_uint<1> comp_r_17_61 = COMP < Q > (input[17].metric, input[61].metric);
    sc_uint<1> comp_r_17_62 = COMP < Q > (input[17].metric, input[62].metric);

    sc_uint<1> comp_r_19_20 = COMP < Q > (input[19].metric, input[20].metric);
    sc_uint<1> comp_r_19_21 = COMP < Q > (input[19].metric, input[21].metric);
    sc_uint<1> comp_r_19_22 = COMP < Q > (input[19].metric, input[22].metric);
    sc_uint<1> comp_r_19_23 = COMP < Q > (input[19].metric, input[23].metric);
    sc_uint<1> comp_r_19_24 = COMP < Q > (input[19].metric, input[24].metric);
    sc_uint<1> comp_r_19_25 = COMP < Q > (input[19].metric, input[25].metric);
    sc_uint<1> comp_r_19_26 = COMP < Q > (input[19].metric, input[26].metric);
    sc_uint<1> comp_r_19_27 = COMP < Q > (input[19].metric, input[27].metric);
    sc_uint<1> comp_r_19_28 = COMP < Q > (input[19].metric, input[28].metric);
    sc_uint<1> comp_r_19_29 = COMP < Q > (input[19].metric, input[29].metric);
    sc_uint<1> comp_r_19_30 = COMP < Q > (input[19].metric, input[30].metric);
    sc_uint<1> comp_r_19_31 = COMP < Q > (input[19].metric, input[31].metric);
    sc_uint<1> comp_r_19_32 = COMP < Q > (input[19].metric, input[32].metric);
    sc_uint<1> comp_r_19_33 = COMP < Q > (input[19].metric, input[33].metric);
    sc_uint<1> comp_r_19_34 = COMP < Q > (input[19].metric, input[34].metric);
    sc_uint<1> comp_r_19_35 = COMP < Q > (input[19].metric, input[35].metric);
    sc_uint<1> comp_r_19_36 = COMP < Q > (input[19].metric, input[36].metric);
    sc_uint<1> comp_r_19_37 = COMP < Q > (input[19].metric, input[37].metric);
    sc_uint<1> comp_r_19_38 = COMP < Q > (input[19].metric, input[38].metric);
    sc_uint<1> comp_r_19_39 = COMP < Q > (input[19].metric, input[39].metric);
    sc_uint<1> comp_r_19_40 = COMP < Q > (input[19].metric, input[40].metric);
    sc_uint<1> comp_r_19_41 = COMP < Q > (input[19].metric, input[41].metric);
    sc_uint<1> comp_r_19_42 = COMP < Q > (input[19].metric, input[42].metric);
    sc_uint<1> comp_r_19_43 = COMP < Q > (input[19].metric, input[43].metric);
    sc_uint<1> comp_r_19_44 = COMP < Q > (input[19].metric, input[44].metric);
    sc_uint<1> comp_r_19_45 = COMP < Q > (input[19].metric, input[45].metric);
    sc_uint<1> comp_r_19_46 = COMP < Q > (input[19].metric, input[46].metric);
    sc_uint<1> comp_r_19_47 = COMP < Q > (input[19].metric, input[47].metric);
    sc_uint<1> comp_r_19_48 = COMP < Q > (input[19].metric, input[48].metric);
    sc_uint<1> comp_r_19_49 = COMP < Q > (input[19].metric, input[49].metric);
    sc_uint<1> comp_r_19_50 = COMP < Q > (input[19].metric, input[50].metric);
    sc_uint<1> comp_r_19_51 = COMP < Q > (input[19].metric, input[51].metric);
    sc_uint<1> comp_r_19_52 = COMP < Q > (input[19].metric, input[52].metric);
    sc_uint<1> comp_r_19_53 = COMP < Q > (input[19].metric, input[53].metric);
    sc_uint<1> comp_r_19_54 = COMP < Q > (input[19].metric, input[54].metric);
    sc_uint<1> comp_r_19_55 = COMP < Q > (input[19].metric, input[55].metric);
    sc_uint<1> comp_r_19_56 = COMP < Q > (input[19].metric, input[56].metric);
    sc_uint<1> comp_r_19_57 = COMP < Q > (input[19].metric, input[57].metric);
    sc_uint<1> comp_r_19_58 = COMP < Q > (input[19].metric, input[58].metric);
    sc_uint<1> comp_r_19_59 = COMP < Q > (input[19].metric, input[59].metric);
    sc_uint<1> comp_r_19_60 = COMP < Q > (input[19].metric, input[60].metric);
    sc_uint<1> comp_r_19_61 = COMP < Q > (input[19].metric, input[61].metric);
    sc_uint<1> comp_r_19_62 = COMP < Q > (input[19].metric, input[62].metric);

    sc_uint<1> comp_r_21_22 = COMP < Q > (input[21].metric, input[22].metric);
    sc_uint<1> comp_r_21_23 = COMP < Q > (input[21].metric, input[23].metric);
    sc_uint<1> comp_r_21_24 = COMP < Q > (input[21].metric, input[24].metric);
    sc_uint<1> comp_r_21_25 = COMP < Q > (input[21].metric, input[25].metric);
    sc_uint<1> comp_r_21_26 = COMP < Q > (input[21].metric, input[26].metric);
    sc_uint<1> comp_r_21_27 = COMP < Q > (input[21].metric, input[27].metric);
    sc_uint<1> comp_r_21_28 = COMP < Q > (input[21].metric, input[28].metric);
    sc_uint<1> comp_r_21_29 = COMP < Q > (input[21].metric, input[29].metric);
    sc_uint<1> comp_r_21_30 = COMP < Q > (input[21].metric, input[30].metric);
    sc_uint<1> comp_r_21_31 = COMP < Q > (input[21].metric, input[31].metric);
    sc_uint<1> comp_r_21_32 = COMP < Q > (input[21].metric, input[32].metric);
    sc_uint<1> comp_r_21_33 = COMP < Q > (input[21].metric, input[33].metric);
    sc_uint<1> comp_r_21_34 = COMP < Q > (input[21].metric, input[34].metric);
    sc_uint<1> comp_r_21_35 = COMP < Q > (input[21].metric, input[35].metric);
    sc_uint<1> comp_r_21_36 = COMP < Q > (input[21].metric, input[36].metric);
    sc_uint<1> comp_r_21_37 = COMP < Q > (input[21].metric, input[37].metric);
    sc_uint<1> comp_r_21_38 = COMP < Q > (input[21].metric, input[38].metric);
    sc_uint<1> comp_r_21_39 = COMP < Q > (input[21].metric, input[39].metric);
    sc_uint<1> comp_r_21_40 = COMP < Q > (input[21].metric, input[40].metric);
    sc_uint<1> comp_r_21_41 = COMP < Q > (input[21].metric, input[41].metric);
    sc_uint<1> comp_r_21_42 = COMP < Q > (input[21].metric, input[42].metric);
    sc_uint<1> comp_r_21_43 = COMP < Q > (input[21].metric, input[43].metric);
    sc_uint<1> comp_r_21_44 = COMP < Q > (input[21].metric, input[44].metric);
    sc_uint<1> comp_r_21_45 = COMP < Q > (input[21].metric, input[45].metric);
    sc_uint<1> comp_r_21_46 = COMP < Q > (input[21].metric, input[46].metric);
    sc_uint<1> comp_r_21_47 = COMP < Q > (input[21].metric, input[47].metric);
    sc_uint<1> comp_r_21_48 = COMP < Q > (input[21].metric, input[48].metric);
    sc_uint<1> comp_r_21_49 = COMP < Q > (input[21].metric, input[49].metric);
    sc_uint<1> comp_r_21_50 = COMP < Q > (input[21].metric, input[50].metric);
    sc_uint<1> comp_r_21_51 = COMP < Q > (input[21].metric, input[51].metric);
    sc_uint<1> comp_r_21_52 = COMP < Q > (input[21].metric, input[52].metric);
    sc_uint<1> comp_r_21_53 = COMP < Q > (input[21].metric, input[53].metric);
    sc_uint<1> comp_r_21_54 = COMP < Q > (input[21].metric, input[54].metric);
    sc_uint<1> comp_r_21_55 = COMP < Q > (input[21].metric, input[55].metric);
    sc_uint<1> comp_r_21_56 = COMP < Q > (input[21].metric, input[56].metric);
    sc_uint<1> comp_r_21_57 = COMP < Q > (input[21].metric, input[57].metric);
    sc_uint<1> comp_r_21_58 = COMP < Q > (input[21].metric, input[58].metric);
    sc_uint<1> comp_r_21_59 = COMP < Q > (input[21].metric, input[59].metric);
    sc_uint<1> comp_r_21_60 = COMP < Q > (input[21].metric, input[60].metric);
    sc_uint<1> comp_r_21_61 = COMP < Q > (input[21].metric, input[61].metric);
    sc_uint<1> comp_r_21_62 = COMP < Q > (input[21].metric, input[62].metric);

    sc_uint<1> comp_r_23_24 = COMP < Q > (input[23].metric, input[24].metric);
    sc_uint<1> comp_r_23_25 = COMP < Q > (input[23].metric, input[25].metric);
    sc_uint<1> comp_r_23_26 = COMP < Q > (input[23].metric, input[26].metric);
    sc_uint<1> comp_r_23_27 = COMP < Q > (input[23].metric, input[27].metric);
    sc_uint<1> comp_r_23_28 = COMP < Q > (input[23].metric, input[28].metric);
    sc_uint<1> comp_r_23_29 = COMP < Q > (input[23].metric, input[29].metric);
    sc_uint<1> comp_r_23_30 = COMP < Q > (input[23].metric, input[30].metric);
    sc_uint<1> comp_r_23_31 = COMP < Q > (input[23].metric, input[31].metric);
    sc_uint<1> comp_r_23_32 = COMP < Q > (input[23].metric, input[32].metric);
    sc_uint<1> comp_r_23_33 = COMP < Q > (input[23].metric, input[33].metric);
    sc_uint<1> comp_r_23_34 = COMP < Q > (input[23].metric, input[34].metric);
    sc_uint<1> comp_r_23_35 = COMP < Q > (input[23].metric, input[35].metric);
    sc_uint<1> comp_r_23_36 = COMP < Q > (input[23].metric, input[36].metric);
    sc_uint<1> comp_r_23_37 = COMP < Q > (input[23].metric, input[37].metric);
    sc_uint<1> comp_r_23_38 = COMP < Q > (input[23].metric, input[38].metric);
    sc_uint<1> comp_r_23_39 = COMP < Q > (input[23].metric, input[39].metric);
    sc_uint<1> comp_r_23_40 = COMP < Q > (input[23].metric, input[40].metric);
    sc_uint<1> comp_r_23_41 = COMP < Q > (input[23].metric, input[41].metric);
    sc_uint<1> comp_r_23_42 = COMP < Q > (input[23].metric, input[42].metric);
    sc_uint<1> comp_r_23_43 = COMP < Q > (input[23].metric, input[43].metric);
    sc_uint<1> comp_r_23_44 = COMP < Q > (input[23].metric, input[44].metric);
    sc_uint<1> comp_r_23_45 = COMP < Q > (input[23].metric, input[45].metric);
    sc_uint<1> comp_r_23_46 = COMP < Q > (input[23].metric, input[46].metric);
    sc_uint<1> comp_r_23_47 = COMP < Q > (input[23].metric, input[47].metric);
    sc_uint<1> comp_r_23_48 = COMP < Q > (input[23].metric, input[48].metric);
    sc_uint<1> comp_r_23_49 = COMP < Q > (input[23].metric, input[49].metric);
    sc_uint<1> comp_r_23_50 = COMP < Q > (input[23].metric, input[50].metric);
    sc_uint<1> comp_r_23_51 = COMP < Q > (input[23].metric, input[51].metric);
    sc_uint<1> comp_r_23_52 = COMP < Q > (input[23].metric, input[52].metric);
    sc_uint<1> comp_r_23_53 = COMP < Q > (input[23].metric, input[53].metric);
    sc_uint<1> comp_r_23_54 = COMP < Q > (input[23].metric, input[54].metric);
    sc_uint<1> comp_r_23_55 = COMP < Q > (input[23].metric, input[55].metric);
    sc_uint<1> comp_r_23_56 = COMP < Q > (input[23].metric, input[56].metric);
    sc_uint<1> comp_r_23_57 = COMP < Q > (input[23].metric, input[57].metric);
    sc_uint<1> comp_r_23_58 = COMP < Q > (input[23].metric, input[58].metric);
    sc_uint<1> comp_r_23_59 = COMP < Q > (input[23].metric, input[59].metric);
    sc_uint<1> comp_r_23_60 = COMP < Q > (input[23].metric, input[60].metric);
    sc_uint<1> comp_r_23_61 = COMP < Q > (input[23].metric, input[61].metric);
    sc_uint<1> comp_r_23_62 = COMP < Q > (input[23].metric, input[62].metric);

    sc_uint<1> comp_r_25_26 = COMP < Q > (input[25].metric, input[26].metric);
    sc_uint<1> comp_r_25_27 = COMP < Q > (input[25].metric, input[27].metric);
    sc_uint<1> comp_r_25_28 = COMP < Q > (input[25].metric, input[28].metric);
    sc_uint<1> comp_r_25_29 = COMP < Q > (input[25].metric, input[29].metric);
    sc_uint<1> comp_r_25_30 = COMP < Q > (input[25].metric, input[30].metric);
    sc_uint<1> comp_r_25_31 = COMP < Q > (input[25].metric, input[31].metric);
    sc_uint<1> comp_r_25_32 = COMP < Q > (input[25].metric, input[32].metric);
    sc_uint<1> comp_r_25_33 = COMP < Q > (input[25].metric, input[33].metric);
    sc_uint<1> comp_r_25_34 = COMP < Q > (input[25].metric, input[34].metric);
    sc_uint<1> comp_r_25_35 = COMP < Q > (input[25].metric, input[35].metric);
    sc_uint<1> comp_r_25_36 = COMP < Q > (input[25].metric, input[36].metric);
    sc_uint<1> comp_r_25_37 = COMP < Q > (input[25].metric, input[37].metric);
    sc_uint<1> comp_r_25_38 = COMP < Q > (input[25].metric, input[38].metric);
    sc_uint<1> comp_r_25_39 = COMP < Q > (input[25].metric, input[39].metric);
    sc_uint<1> comp_r_25_40 = COMP < Q > (input[25].metric, input[40].metric);
    sc_uint<1> comp_r_25_41 = COMP < Q > (input[25].metric, input[41].metric);
    sc_uint<1> comp_r_25_42 = COMP < Q > (input[25].metric, input[42].metric);
    sc_uint<1> comp_r_25_43 = COMP < Q > (input[25].metric, input[43].metric);
    sc_uint<1> comp_r_25_44 = COMP < Q > (input[25].metric, input[44].metric);
    sc_uint<1> comp_r_25_45 = COMP < Q > (input[25].metric, input[45].metric);
    sc_uint<1> comp_r_25_46 = COMP < Q > (input[25].metric, input[46].metric);
    sc_uint<1> comp_r_25_47 = COMP < Q > (input[25].metric, input[47].metric);
    sc_uint<1> comp_r_25_48 = COMP < Q > (input[25].metric, input[48].metric);
    sc_uint<1> comp_r_25_49 = COMP < Q > (input[25].metric, input[49].metric);
    sc_uint<1> comp_r_25_50 = COMP < Q > (input[25].metric, input[50].metric);
    sc_uint<1> comp_r_25_51 = COMP < Q > (input[25].metric, input[51].metric);
    sc_uint<1> comp_r_25_52 = COMP < Q > (input[25].metric, input[52].metric);
    sc_uint<1> comp_r_25_53 = COMP < Q > (input[25].metric, input[53].metric);
    sc_uint<1> comp_r_25_54 = COMP < Q > (input[25].metric, input[54].metric);
    sc_uint<1> comp_r_25_55 = COMP < Q > (input[25].metric, input[55].metric);
    sc_uint<1> comp_r_25_56 = COMP < Q > (input[25].metric, input[56].metric);
    sc_uint<1> comp_r_25_57 = COMP < Q > (input[25].metric, input[57].metric);
    sc_uint<1> comp_r_25_58 = COMP < Q > (input[25].metric, input[58].metric);
    sc_uint<1> comp_r_25_59 = COMP < Q > (input[25].metric, input[59].metric);
    sc_uint<1> comp_r_25_60 = COMP < Q > (input[25].metric, input[60].metric);
    sc_uint<1> comp_r_25_61 = COMP < Q > (input[25].metric, input[61].metric);
    sc_uint<1> comp_r_25_62 = COMP < Q > (input[25].metric, input[62].metric);

    sc_uint<1> comp_r_27_28 = COMP < Q > (input[27].metric, input[28].metric);
    sc_uint<1> comp_r_27_29 = COMP < Q > (input[27].metric, input[29].metric);
    sc_uint<1> comp_r_27_30 = COMP < Q > (input[27].metric, input[30].metric);
    sc_uint<1> comp_r_27_31 = COMP < Q > (input[27].metric, input[31].metric);
    sc_uint<1> comp_r_27_32 = COMP < Q > (input[27].metric, input[32].metric);
    sc_uint<1> comp_r_27_33 = COMP < Q > (input[27].metric, input[33].metric);
    sc_uint<1> comp_r_27_34 = COMP < Q > (input[27].metric, input[34].metric);
    sc_uint<1> comp_r_27_35 = COMP < Q > (input[27].metric, input[35].metric);
    sc_uint<1> comp_r_27_36 = COMP < Q > (input[27].metric, input[36].metric);
    sc_uint<1> comp_r_27_37 = COMP < Q > (input[27].metric, input[37].metric);
    sc_uint<1> comp_r_27_38 = COMP < Q > (input[27].metric, input[38].metric);
    sc_uint<1> comp_r_27_39 = COMP < Q > (input[27].metric, input[39].metric);
    sc_uint<1> comp_r_27_40 = COMP < Q > (input[27].metric, input[40].metric);
    sc_uint<1> comp_r_27_41 = COMP < Q > (input[27].metric, input[41].metric);
    sc_uint<1> comp_r_27_42 = COMP < Q > (input[27].metric, input[42].metric);
    sc_uint<1> comp_r_27_43 = COMP < Q > (input[27].metric, input[43].metric);
    sc_uint<1> comp_r_27_44 = COMP < Q > (input[27].metric, input[44].metric);
    sc_uint<1> comp_r_27_45 = COMP < Q > (input[27].metric, input[45].metric);
    sc_uint<1> comp_r_27_46 = COMP < Q > (input[27].metric, input[46].metric);
    sc_uint<1> comp_r_27_47 = COMP < Q > (input[27].metric, input[47].metric);
    sc_uint<1> comp_r_27_48 = COMP < Q > (input[27].metric, input[48].metric);
    sc_uint<1> comp_r_27_49 = COMP < Q > (input[27].metric, input[49].metric);
    sc_uint<1> comp_r_27_50 = COMP < Q > (input[27].metric, input[50].metric);
    sc_uint<1> comp_r_27_51 = COMP < Q > (input[27].metric, input[51].metric);
    sc_uint<1> comp_r_27_52 = COMP < Q > (input[27].metric, input[52].metric);
    sc_uint<1> comp_r_27_53 = COMP < Q > (input[27].metric, input[53].metric);
    sc_uint<1> comp_r_27_54 = COMP < Q > (input[27].metric, input[54].metric);
    sc_uint<1> comp_r_27_55 = COMP < Q > (input[27].metric, input[55].metric);
    sc_uint<1> comp_r_27_56 = COMP < Q > (input[27].metric, input[56].metric);
    sc_uint<1> comp_r_27_57 = COMP < Q > (input[27].metric, input[57].metric);
    sc_uint<1> comp_r_27_58 = COMP < Q > (input[27].metric, input[58].metric);
    sc_uint<1> comp_r_27_59 = COMP < Q > (input[27].metric, input[59].metric);
    sc_uint<1> comp_r_27_60 = COMP < Q > (input[27].metric, input[60].metric);
    sc_uint<1> comp_r_27_61 = COMP < Q > (input[27].metric, input[61].metric);
    sc_uint<1> comp_r_27_62 = COMP < Q > (input[27].metric, input[62].metric);

    sc_uint<1> comp_r_29_30 = COMP < Q > (input[29].metric, input[30].metric);
    sc_uint<1> comp_r_29_31 = COMP < Q > (input[29].metric, input[31].metric);
    sc_uint<1> comp_r_29_32 = COMP < Q > (input[29].metric, input[32].metric);
    sc_uint<1> comp_r_29_33 = COMP < Q > (input[29].metric, input[33].metric);
    sc_uint<1> comp_r_29_34 = COMP < Q > (input[29].metric, input[34].metric);
    sc_uint<1> comp_r_29_35 = COMP < Q > (input[29].metric, input[35].metric);
    sc_uint<1> comp_r_29_36 = COMP < Q > (input[29].metric, input[36].metric);
    sc_uint<1> comp_r_29_37 = COMP < Q > (input[29].metric, input[37].metric);
    sc_uint<1> comp_r_29_38 = COMP < Q > (input[29].metric, input[38].metric);
    sc_uint<1> comp_r_29_39 = COMP < Q > (input[29].metric, input[39].metric);
    sc_uint<1> comp_r_29_40 = COMP < Q > (input[29].metric, input[40].metric);
    sc_uint<1> comp_r_29_41 = COMP < Q > (input[29].metric, input[41].metric);
    sc_uint<1> comp_r_29_42 = COMP < Q > (input[29].metric, input[42].metric);
    sc_uint<1> comp_r_29_43 = COMP < Q > (input[29].metric, input[43].metric);
    sc_uint<1> comp_r_29_44 = COMP < Q > (input[29].metric, input[44].metric);
    sc_uint<1> comp_r_29_45 = COMP < Q > (input[29].metric, input[45].metric);
    sc_uint<1> comp_r_29_46 = COMP < Q > (input[29].metric, input[46].metric);
    sc_uint<1> comp_r_29_47 = COMP < Q > (input[29].metric, input[47].metric);
    sc_uint<1> comp_r_29_48 = COMP < Q > (input[29].metric, input[48].metric);
    sc_uint<1> comp_r_29_49 = COMP < Q > (input[29].metric, input[49].metric);
    sc_uint<1> comp_r_29_50 = COMP < Q > (input[29].metric, input[50].metric);
    sc_uint<1> comp_r_29_51 = COMP < Q > (input[29].metric, input[51].metric);
    sc_uint<1> comp_r_29_52 = COMP < Q > (input[29].metric, input[52].metric);
    sc_uint<1> comp_r_29_53 = COMP < Q > (input[29].metric, input[53].metric);
    sc_uint<1> comp_r_29_54 = COMP < Q > (input[29].metric, input[54].metric);
    sc_uint<1> comp_r_29_55 = COMP < Q > (input[29].metric, input[55].metric);
    sc_uint<1> comp_r_29_56 = COMP < Q > (input[29].metric, input[56].metric);
    sc_uint<1> comp_r_29_57 = COMP < Q > (input[29].metric, input[57].metric);
    sc_uint<1> comp_r_29_58 = COMP < Q > (input[29].metric, input[58].metric);
    sc_uint<1> comp_r_29_59 = COMP < Q > (input[29].metric, input[59].metric);
    sc_uint<1> comp_r_29_60 = COMP < Q > (input[29].metric, input[60].metric);
    sc_uint<1> comp_r_29_61 = COMP < Q > (input[29].metric, input[61].metric);
    sc_uint<1> comp_r_29_62 = COMP < Q > (input[29].metric, input[62].metric);

    sc_uint<1> comp_r_31_32 = COMP < Q > (input[31].metric, input[32].metric);
    sc_uint<1> comp_r_31_33 = COMP < Q > (input[31].metric, input[33].metric);
    sc_uint<1> comp_r_31_34 = COMP < Q > (input[31].metric, input[34].metric);
    sc_uint<1> comp_r_31_35 = COMP < Q > (input[31].metric, input[35].metric);
    sc_uint<1> comp_r_31_36 = COMP < Q > (input[31].metric, input[36].metric);
    sc_uint<1> comp_r_31_37 = COMP < Q > (input[31].metric, input[37].metric);
    sc_uint<1> comp_r_31_38 = COMP < Q > (input[31].metric, input[38].metric);
    sc_uint<1> comp_r_31_39 = COMP < Q > (input[31].metric, input[39].metric);
    sc_uint<1> comp_r_31_40 = COMP < Q > (input[31].metric, input[40].metric);
    sc_uint<1> comp_r_31_41 = COMP < Q > (input[31].metric, input[41].metric);
    sc_uint<1> comp_r_31_42 = COMP < Q > (input[31].metric, input[42].metric);
    sc_uint<1> comp_r_31_43 = COMP < Q > (input[31].metric, input[43].metric);
    sc_uint<1> comp_r_31_44 = COMP < Q > (input[31].metric, input[44].metric);
    sc_uint<1> comp_r_31_45 = COMP < Q > (input[31].metric, input[45].metric);
    sc_uint<1> comp_r_31_46 = COMP < Q > (input[31].metric, input[46].metric);
    sc_uint<1> comp_r_31_47 = COMP < Q > (input[31].metric, input[47].metric);
    sc_uint<1> comp_r_31_48 = COMP < Q > (input[31].metric, input[48].metric);
    sc_uint<1> comp_r_31_49 = COMP < Q > (input[31].metric, input[49].metric);
    sc_uint<1> comp_r_31_50 = COMP < Q > (input[31].metric, input[50].metric);
    sc_uint<1> comp_r_31_51 = COMP < Q > (input[31].metric, input[51].metric);
    sc_uint<1> comp_r_31_52 = COMP < Q > (input[31].metric, input[52].metric);
    sc_uint<1> comp_r_31_53 = COMP < Q > (input[31].metric, input[53].metric);
    sc_uint<1> comp_r_31_54 = COMP < Q > (input[31].metric, input[54].metric);
    sc_uint<1> comp_r_31_55 = COMP < Q > (input[31].metric, input[55].metric);
    sc_uint<1> comp_r_31_56 = COMP < Q > (input[31].metric, input[56].metric);
    sc_uint<1> comp_r_31_57 = COMP < Q > (input[31].metric, input[57].metric);
    sc_uint<1> comp_r_31_58 = COMP < Q > (input[31].metric, input[58].metric);
    sc_uint<1> comp_r_31_59 = COMP < Q > (input[31].metric, input[59].metric);
    sc_uint<1> comp_r_31_60 = COMP < Q > (input[31].metric, input[60].metric);
    sc_uint<1> comp_r_31_61 = COMP < Q > (input[31].metric, input[61].metric);
    sc_uint<1> comp_r_31_62 = COMP < Q > (input[31].metric, input[62].metric);

    sc_uint<1> comp_r_33_34 = COMP < Q > (input[33].metric, input[34].metric);
    sc_uint<1> comp_r_33_35 = COMP < Q > (input[33].metric, input[35].metric);
    sc_uint<1> comp_r_33_36 = COMP < Q > (input[33].metric, input[36].metric);
    sc_uint<1> comp_r_33_37 = COMP < Q > (input[33].metric, input[37].metric);
    sc_uint<1> comp_r_33_38 = COMP < Q > (input[33].metric, input[38].metric);
    sc_uint<1> comp_r_33_39 = COMP < Q > (input[33].metric, input[39].metric);
    sc_uint<1> comp_r_33_40 = COMP < Q > (input[33].metric, input[40].metric);
    sc_uint<1> comp_r_33_41 = COMP < Q > (input[33].metric, input[41].metric);
    sc_uint<1> comp_r_33_42 = COMP < Q > (input[33].metric, input[42].metric);
    sc_uint<1> comp_r_33_43 = COMP < Q > (input[33].metric, input[43].metric);
    sc_uint<1> comp_r_33_44 = COMP < Q > (input[33].metric, input[44].metric);
    sc_uint<1> comp_r_33_45 = COMP < Q > (input[33].metric, input[45].metric);
    sc_uint<1> comp_r_33_46 = COMP < Q > (input[33].metric, input[46].metric);
    sc_uint<1> comp_r_33_47 = COMP < Q > (input[33].metric, input[47].metric);
    sc_uint<1> comp_r_33_48 = COMP < Q > (input[33].metric, input[48].metric);
    sc_uint<1> comp_r_33_49 = COMP < Q > (input[33].metric, input[49].metric);
    sc_uint<1> comp_r_33_50 = COMP < Q > (input[33].metric, input[50].metric);
    sc_uint<1> comp_r_33_51 = COMP < Q > (input[33].metric, input[51].metric);
    sc_uint<1> comp_r_33_52 = COMP < Q > (input[33].metric, input[52].metric);
    sc_uint<1> comp_r_33_53 = COMP < Q > (input[33].metric, input[53].metric);
    sc_uint<1> comp_r_33_54 = COMP < Q > (input[33].metric, input[54].metric);
    sc_uint<1> comp_r_33_55 = COMP < Q > (input[33].metric, input[55].metric);
    sc_uint<1> comp_r_33_56 = COMP < Q > (input[33].metric, input[56].metric);
    sc_uint<1> comp_r_33_57 = COMP < Q > (input[33].metric, input[57].metric);
    sc_uint<1> comp_r_33_58 = COMP < Q > (input[33].metric, input[58].metric);
    sc_uint<1> comp_r_33_59 = COMP < Q > (input[33].metric, input[59].metric);
    sc_uint<1> comp_r_33_60 = COMP < Q > (input[33].metric, input[60].metric);
    sc_uint<1> comp_r_33_61 = COMP < Q > (input[33].metric, input[61].metric);
    sc_uint<1> comp_r_33_62 = COMP < Q > (input[33].metric, input[62].metric);

    sc_uint<1> comp_r_35_36 = COMP < Q > (input[35].metric, input[36].metric);
    sc_uint<1> comp_r_35_37 = COMP < Q > (input[35].metric, input[37].metric);
    sc_uint<1> comp_r_35_38 = COMP < Q > (input[35].metric, input[38].metric);
    sc_uint<1> comp_r_35_39 = COMP < Q > (input[35].metric, input[39].metric);
    sc_uint<1> comp_r_35_40 = COMP < Q > (input[35].metric, input[40].metric);
    sc_uint<1> comp_r_35_41 = COMP < Q > (input[35].metric, input[41].metric);
    sc_uint<1> comp_r_35_42 = COMP < Q > (input[35].metric, input[42].metric);
    sc_uint<1> comp_r_35_43 = COMP < Q > (input[35].metric, input[43].metric);
    sc_uint<1> comp_r_35_44 = COMP < Q > (input[35].metric, input[44].metric);
    sc_uint<1> comp_r_35_45 = COMP < Q > (input[35].metric, input[45].metric);
    sc_uint<1> comp_r_35_46 = COMP < Q > (input[35].metric, input[46].metric);
    sc_uint<1> comp_r_35_47 = COMP < Q > (input[35].metric, input[47].metric);
    sc_uint<1> comp_r_35_48 = COMP < Q > (input[35].metric, input[48].metric);
    sc_uint<1> comp_r_35_49 = COMP < Q > (input[35].metric, input[49].metric);
    sc_uint<1> comp_r_35_50 = COMP < Q > (input[35].metric, input[50].metric);
    sc_uint<1> comp_r_35_51 = COMP < Q > (input[35].metric, input[51].metric);
    sc_uint<1> comp_r_35_52 = COMP < Q > (input[35].metric, input[52].metric);
    sc_uint<1> comp_r_35_53 = COMP < Q > (input[35].metric, input[53].metric);
    sc_uint<1> comp_r_35_54 = COMP < Q > (input[35].metric, input[54].metric);
    sc_uint<1> comp_r_35_55 = COMP < Q > (input[35].metric, input[55].metric);
    sc_uint<1> comp_r_35_56 = COMP < Q > (input[35].metric, input[56].metric);
    sc_uint<1> comp_r_35_57 = COMP < Q > (input[35].metric, input[57].metric);
    sc_uint<1> comp_r_35_58 = COMP < Q > (input[35].metric, input[58].metric);
    sc_uint<1> comp_r_35_59 = COMP < Q > (input[35].metric, input[59].metric);
    sc_uint<1> comp_r_35_60 = COMP < Q > (input[35].metric, input[60].metric);
    sc_uint<1> comp_r_35_61 = COMP < Q > (input[35].metric, input[61].metric);
    sc_uint<1> comp_r_35_62 = COMP < Q > (input[35].metric, input[62].metric);

    sc_uint<1> comp_r_37_38 = COMP < Q > (input[37].metric, input[38].metric);
    sc_uint<1> comp_r_37_39 = COMP < Q > (input[37].metric, input[39].metric);
    sc_uint<1> comp_r_37_40 = COMP < Q > (input[37].metric, input[40].metric);
    sc_uint<1> comp_r_37_41 = COMP < Q > (input[37].metric, input[41].metric);
    sc_uint<1> comp_r_37_42 = COMP < Q > (input[37].metric, input[42].metric);
    sc_uint<1> comp_r_37_43 = COMP < Q > (input[37].metric, input[43].metric);
    sc_uint<1> comp_r_37_44 = COMP < Q > (input[37].metric, input[44].metric);
    sc_uint<1> comp_r_37_45 = COMP < Q > (input[37].metric, input[45].metric);
    sc_uint<1> comp_r_37_46 = COMP < Q > (input[37].metric, input[46].metric);
    sc_uint<1> comp_r_37_47 = COMP < Q > (input[37].metric, input[47].metric);
    sc_uint<1> comp_r_37_48 = COMP < Q > (input[37].metric, input[48].metric);
    sc_uint<1> comp_r_37_49 = COMP < Q > (input[37].metric, input[49].metric);
    sc_uint<1> comp_r_37_50 = COMP < Q > (input[37].metric, input[50].metric);
    sc_uint<1> comp_r_37_51 = COMP < Q > (input[37].metric, input[51].metric);
    sc_uint<1> comp_r_37_52 = COMP < Q > (input[37].metric, input[52].metric);
    sc_uint<1> comp_r_37_53 = COMP < Q > (input[37].metric, input[53].metric);
    sc_uint<1> comp_r_37_54 = COMP < Q > (input[37].metric, input[54].metric);
    sc_uint<1> comp_r_37_55 = COMP < Q > (input[37].metric, input[55].metric);
    sc_uint<1> comp_r_37_56 = COMP < Q > (input[37].metric, input[56].metric);
    sc_uint<1> comp_r_37_57 = COMP < Q > (input[37].metric, input[57].metric);
    sc_uint<1> comp_r_37_58 = COMP < Q > (input[37].metric, input[58].metric);
    sc_uint<1> comp_r_37_59 = COMP < Q > (input[37].metric, input[59].metric);
    sc_uint<1> comp_r_37_60 = COMP < Q > (input[37].metric, input[60].metric);
    sc_uint<1> comp_r_37_61 = COMP < Q > (input[37].metric, input[61].metric);
    sc_uint<1> comp_r_37_62 = COMP < Q > (input[37].metric, input[62].metric);

    sc_uint<1> comp_r_39_40 = COMP < Q > (input[39].metric, input[40].metric);
    sc_uint<1> comp_r_39_41 = COMP < Q > (input[39].metric, input[41].metric);
    sc_uint<1> comp_r_39_42 = COMP < Q > (input[39].metric, input[42].metric);
    sc_uint<1> comp_r_39_43 = COMP < Q > (input[39].metric, input[43].metric);
    sc_uint<1> comp_r_39_44 = COMP < Q > (input[39].metric, input[44].metric);
    sc_uint<1> comp_r_39_45 = COMP < Q > (input[39].metric, input[45].metric);
    sc_uint<1> comp_r_39_46 = COMP < Q > (input[39].metric, input[46].metric);
    sc_uint<1> comp_r_39_47 = COMP < Q > (input[39].metric, input[47].metric);
    sc_uint<1> comp_r_39_48 = COMP < Q > (input[39].metric, input[48].metric);
    sc_uint<1> comp_r_39_49 = COMP < Q > (input[39].metric, input[49].metric);
    sc_uint<1> comp_r_39_50 = COMP < Q > (input[39].metric, input[50].metric);
    sc_uint<1> comp_r_39_51 = COMP < Q > (input[39].metric, input[51].metric);
    sc_uint<1> comp_r_39_52 = COMP < Q > (input[39].metric, input[52].metric);
    sc_uint<1> comp_r_39_53 = COMP < Q > (input[39].metric, input[53].metric);
    sc_uint<1> comp_r_39_54 = COMP < Q > (input[39].metric, input[54].metric);
    sc_uint<1> comp_r_39_55 = COMP < Q > (input[39].metric, input[55].metric);
    sc_uint<1> comp_r_39_56 = COMP < Q > (input[39].metric, input[56].metric);
    sc_uint<1> comp_r_39_57 = COMP < Q > (input[39].metric, input[57].metric);
    sc_uint<1> comp_r_39_58 = COMP < Q > (input[39].metric, input[58].metric);
    sc_uint<1> comp_r_39_59 = COMP < Q > (input[39].metric, input[59].metric);
    sc_uint<1> comp_r_39_60 = COMP < Q > (input[39].metric, input[60].metric);
    sc_uint<1> comp_r_39_61 = COMP < Q > (input[39].metric, input[61].metric);
    sc_uint<1> comp_r_39_62 = COMP < Q > (input[39].metric, input[62].metric);

    sc_uint<1> comp_r_41_42 = COMP < Q > (input[41].metric, input[42].metric);
    sc_uint<1> comp_r_41_43 = COMP < Q > (input[41].metric, input[43].metric);
    sc_uint<1> comp_r_41_44 = COMP < Q > (input[41].metric, input[44].metric);
    sc_uint<1> comp_r_41_45 = COMP < Q > (input[41].metric, input[45].metric);
    sc_uint<1> comp_r_41_46 = COMP < Q > (input[41].metric, input[46].metric);
    sc_uint<1> comp_r_41_47 = COMP < Q > (input[41].metric, input[47].metric);
    sc_uint<1> comp_r_41_48 = COMP < Q > (input[41].metric, input[48].metric);
    sc_uint<1> comp_r_41_49 = COMP < Q > (input[41].metric, input[49].metric);
    sc_uint<1> comp_r_41_50 = COMP < Q > (input[41].metric, input[50].metric);
    sc_uint<1> comp_r_41_51 = COMP < Q > (input[41].metric, input[51].metric);
    sc_uint<1> comp_r_41_52 = COMP < Q > (input[41].metric, input[52].metric);
    sc_uint<1> comp_r_41_53 = COMP < Q > (input[41].metric, input[53].metric);
    sc_uint<1> comp_r_41_54 = COMP < Q > (input[41].metric, input[54].metric);
    sc_uint<1> comp_r_41_55 = COMP < Q > (input[41].metric, input[55].metric);
    sc_uint<1> comp_r_41_56 = COMP < Q > (input[41].metric, input[56].metric);
    sc_uint<1> comp_r_41_57 = COMP < Q > (input[41].metric, input[57].metric);
    sc_uint<1> comp_r_41_58 = COMP < Q > (input[41].metric, input[58].metric);
    sc_uint<1> comp_r_41_59 = COMP < Q > (input[41].metric, input[59].metric);
    sc_uint<1> comp_r_41_60 = COMP < Q > (input[41].metric, input[60].metric);
    sc_uint<1> comp_r_41_61 = COMP < Q > (input[41].metric, input[61].metric);
    sc_uint<1> comp_r_41_62 = COMP < Q > (input[41].metric, input[62].metric);

    sc_uint<1> comp_r_43_44 = COMP < Q > (input[43].metric, input[44].metric);
    sc_uint<1> comp_r_43_45 = COMP < Q > (input[43].metric, input[45].metric);
    sc_uint<1> comp_r_43_46 = COMP < Q > (input[43].metric, input[46].metric);
    sc_uint<1> comp_r_43_47 = COMP < Q > (input[43].metric, input[47].metric);
    sc_uint<1> comp_r_43_48 = COMP < Q > (input[43].metric, input[48].metric);
    sc_uint<1> comp_r_43_49 = COMP < Q > (input[43].metric, input[49].metric);
    sc_uint<1> comp_r_43_50 = COMP < Q > (input[43].metric, input[50].metric);
    sc_uint<1> comp_r_43_51 = COMP < Q > (input[43].metric, input[51].metric);
    sc_uint<1> comp_r_43_52 = COMP < Q > (input[43].metric, input[52].metric);
    sc_uint<1> comp_r_43_53 = COMP < Q > (input[43].metric, input[53].metric);
    sc_uint<1> comp_r_43_54 = COMP < Q > (input[43].metric, input[54].metric);
    sc_uint<1> comp_r_43_55 = COMP < Q > (input[43].metric, input[55].metric);
    sc_uint<1> comp_r_43_56 = COMP < Q > (input[43].metric, input[56].metric);
    sc_uint<1> comp_r_43_57 = COMP < Q > (input[43].metric, input[57].metric);
    sc_uint<1> comp_r_43_58 = COMP < Q > (input[43].metric, input[58].metric);
    sc_uint<1> comp_r_43_59 = COMP < Q > (input[43].metric, input[59].metric);
    sc_uint<1> comp_r_43_60 = COMP < Q > (input[43].metric, input[60].metric);
    sc_uint<1> comp_r_43_61 = COMP < Q > (input[43].metric, input[61].metric);
    sc_uint<1> comp_r_43_62 = COMP < Q > (input[43].metric, input[62].metric);

    sc_uint<1> comp_r_45_46 = COMP < Q > (input[45].metric, input[46].metric);
    sc_uint<1> comp_r_45_47 = COMP < Q > (input[45].metric, input[47].metric);
    sc_uint<1> comp_r_45_48 = COMP < Q > (input[45].metric, input[48].metric);
    sc_uint<1> comp_r_45_49 = COMP < Q > (input[45].metric, input[49].metric);
    sc_uint<1> comp_r_45_50 = COMP < Q > (input[45].metric, input[50].metric);
    sc_uint<1> comp_r_45_51 = COMP < Q > (input[45].metric, input[51].metric);
    sc_uint<1> comp_r_45_52 = COMP < Q > (input[45].metric, input[52].metric);
    sc_uint<1> comp_r_45_53 = COMP < Q > (input[45].metric, input[53].metric);
    sc_uint<1> comp_r_45_54 = COMP < Q > (input[45].metric, input[54].metric);
    sc_uint<1> comp_r_45_55 = COMP < Q > (input[45].metric, input[55].metric);
    sc_uint<1> comp_r_45_56 = COMP < Q > (input[45].metric, input[56].metric);
    sc_uint<1> comp_r_45_57 = COMP < Q > (input[45].metric, input[57].metric);
    sc_uint<1> comp_r_45_58 = COMP < Q > (input[45].metric, input[58].metric);
    sc_uint<1> comp_r_45_59 = COMP < Q > (input[45].metric, input[59].metric);
    sc_uint<1> comp_r_45_60 = COMP < Q > (input[45].metric, input[60].metric);
    sc_uint<1> comp_r_45_61 = COMP < Q > (input[45].metric, input[61].metric);
    sc_uint<1> comp_r_45_62 = COMP < Q > (input[45].metric, input[62].metric);

    sc_uint<1> comp_r_47_48 = COMP < Q > (input[47].metric, input[48].metric);
    sc_uint<1> comp_r_47_49 = COMP < Q > (input[47].metric, input[49].metric);
    sc_uint<1> comp_r_47_50 = COMP < Q > (input[47].metric, input[50].metric);
    sc_uint<1> comp_r_47_51 = COMP < Q > (input[47].metric, input[51].metric);
    sc_uint<1> comp_r_47_52 = COMP < Q > (input[47].metric, input[52].metric);
    sc_uint<1> comp_r_47_53 = COMP < Q > (input[47].metric, input[53].metric);
    sc_uint<1> comp_r_47_54 = COMP < Q > (input[47].metric, input[54].metric);
    sc_uint<1> comp_r_47_55 = COMP < Q > (input[47].metric, input[55].metric);
    sc_uint<1> comp_r_47_56 = COMP < Q > (input[47].metric, input[56].metric);
    sc_uint<1> comp_r_47_57 = COMP < Q > (input[47].metric, input[57].metric);
    sc_uint<1> comp_r_47_58 = COMP < Q > (input[47].metric, input[58].metric);
    sc_uint<1> comp_r_47_59 = COMP < Q > (input[47].metric, input[59].metric);
    sc_uint<1> comp_r_47_60 = COMP < Q > (input[47].metric, input[60].metric);
    sc_uint<1> comp_r_47_61 = COMP < Q > (input[47].metric, input[61].metric);
    sc_uint<1> comp_r_47_62 = COMP < Q > (input[47].metric, input[62].metric);

    sc_uint<1> comp_r_49_50 = COMP < Q > (input[49].metric, input[50].metric);
    sc_uint<1> comp_r_49_51 = COMP < Q > (input[49].metric, input[51].metric);
    sc_uint<1> comp_r_49_52 = COMP < Q > (input[49].metric, input[52].metric);
    sc_uint<1> comp_r_49_53 = COMP < Q > (input[49].metric, input[53].metric);
    sc_uint<1> comp_r_49_54 = COMP < Q > (input[49].metric, input[54].metric);
    sc_uint<1> comp_r_49_55 = COMP < Q > (input[49].metric, input[55].metric);
    sc_uint<1> comp_r_49_56 = COMP < Q > (input[49].metric, input[56].metric);
    sc_uint<1> comp_r_49_57 = COMP < Q > (input[49].metric, input[57].metric);
    sc_uint<1> comp_r_49_58 = COMP < Q > (input[49].metric, input[58].metric);
    sc_uint<1> comp_r_49_59 = COMP < Q > (input[49].metric, input[59].metric);
    sc_uint<1> comp_r_49_60 = COMP < Q > (input[49].metric, input[60].metric);
    sc_uint<1> comp_r_49_61 = COMP < Q > (input[49].metric, input[61].metric);
    sc_uint<1> comp_r_49_62 = COMP < Q > (input[49].metric, input[62].metric);

    sc_uint<1> comp_r_51_52 = COMP < Q > (input[51].metric, input[52].metric);
    sc_uint<1> comp_r_51_53 = COMP < Q > (input[51].metric, input[53].metric);
    sc_uint<1> comp_r_51_54 = COMP < Q > (input[51].metric, input[54].metric);
    sc_uint<1> comp_r_51_55 = COMP < Q > (input[51].metric, input[55].metric);
    sc_uint<1> comp_r_51_56 = COMP < Q > (input[51].metric, input[56].metric);
    sc_uint<1> comp_r_51_57 = COMP < Q > (input[51].metric, input[57].metric);
    sc_uint<1> comp_r_51_58 = COMP < Q > (input[51].metric, input[58].metric);
    sc_uint<1> comp_r_51_59 = COMP < Q > (input[51].metric, input[59].metric);
    sc_uint<1> comp_r_51_60 = COMP < Q > (input[51].metric, input[60].metric);
    sc_uint<1> comp_r_51_61 = COMP < Q > (input[51].metric, input[61].metric);
    sc_uint<1> comp_r_51_62 = COMP < Q > (input[51].metric, input[62].metric);

    sc_uint<1> comp_r_53_54 = COMP < Q > (input[53].metric, input[54].metric);
    sc_uint<1> comp_r_53_55 = COMP < Q > (input[53].metric, input[55].metric);
    sc_uint<1> comp_r_53_56 = COMP < Q > (input[53].metric, input[56].metric);
    sc_uint<1> comp_r_53_57 = COMP < Q > (input[53].metric, input[57].metric);
    sc_uint<1> comp_r_53_58 = COMP < Q > (input[53].metric, input[58].metric);
    sc_uint<1> comp_r_53_59 = COMP < Q > (input[53].metric, input[59].metric);
    sc_uint<1> comp_r_53_60 = COMP < Q > (input[53].metric, input[60].metric);
    sc_uint<1> comp_r_53_61 = COMP < Q > (input[53].metric, input[61].metric);
    sc_uint<1> comp_r_53_62 = COMP < Q > (input[53].metric, input[62].metric);

    sc_uint<1> comp_r_55_56 = COMP < Q > (input[55].metric, input[56].metric);
    sc_uint<1> comp_r_55_57 = COMP < Q > (input[55].metric, input[57].metric);
    sc_uint<1> comp_r_55_58 = COMP < Q > (input[55].metric, input[58].metric);
    sc_uint<1> comp_r_55_59 = COMP < Q > (input[55].metric, input[59].metric);
    sc_uint<1> comp_r_55_60 = COMP < Q > (input[55].metric, input[60].metric);
    sc_uint<1> comp_r_55_61 = COMP < Q > (input[55].metric, input[61].metric);
    sc_uint<1> comp_r_55_62 = COMP < Q > (input[55].metric, input[62].metric);

    sc_uint<1> comp_r_57_58 = COMP < Q > (input[57].metric, input[58].metric);
    sc_uint<1> comp_r_57_59 = COMP < Q > (input[57].metric, input[59].metric);
    sc_uint<1> comp_r_57_60 = COMP < Q > (input[57].metric, input[60].metric);
    sc_uint<1> comp_r_57_61 = COMP < Q > (input[57].metric, input[61].metric);
    sc_uint<1> comp_r_57_62 = COMP < Q > (input[57].metric, input[62].metric);

    sc_uint<1> comp_r_59_60 = COMP < Q > (input[59].metric, input[60].metric);
    sc_uint<1> comp_r_59_61 = COMP < Q > (input[59].metric, input[61].metric);
    sc_uint<1> comp_r_59_62 = COMP < Q > (input[59].metric, input[62].metric);

    sc_uint<1> comp_r_61_62 = COMP < Q > (input[61].metric, input[62].metric);

// COMPUTE positions

    position[0] = 0;
    position[1] = 1 + comp_r_1_2 + comp_r_1_3 + comp_r_1_4 + comp_r_1_5 + comp_r_1_6 + comp_r_1_7 + comp_r_1_8 + comp_r_1_9 + comp_r_1_10 + comp_r_1_11 + comp_r_1_12 + comp_r_1_13 + comp_r_1_14 + comp_r_1_15 + comp_r_1_16 + comp_r_1_17 + comp_r_1_18 + comp_r_1_19 + comp_r_1_20 + comp_r_1_21 + comp_r_1_22 + comp_r_1_23 + comp_r_1_24 + comp_r_1_25 + comp_r_1_26 + comp_r_1_27 + comp_r_1_28 + comp_r_1_29 + comp_r_1_30 + comp_r_1_31 + comp_r_1_32 + comp_r_1_33 + comp_r_1_34 + comp_r_1_35 + comp_r_1_36 + comp_r_1_37 + comp_r_1_38 + comp_r_1_39 + comp_r_1_40 + comp_r_1_41 + comp_r_1_42 + comp_r_1_43 + comp_r_1_44 + comp_r_1_45 + comp_r_1_46 + comp_r_1_47 + comp_r_1_48 + comp_r_1_49 + comp_r_1_50 + comp_r_1_51 + comp_r_1_52 + comp_r_1_53 + comp_r_1_54 + comp_r_1_55 + comp_r_1_56 + comp_r_1_57 + comp_r_1_58 + comp_r_1_59 + comp_r_1_60 + comp_r_1_61 + comp_r_1_62;
    position[2] = 1 + (sc_uint<1>) ~(comp_r_1_2);
    position[3] = 2 + (sc_uint<1>) ~(comp_r_1_3) + comp_r_3_4 + comp_r_3_5 + comp_r_3_6 + comp_r_3_7 + comp_r_3_8 + comp_r_3_9 + comp_r_3_10 + comp_r_3_11 + comp_r_3_12 + comp_r_3_13 + comp_r_3_14 + comp_r_3_15 + comp_r_3_16 + comp_r_3_17 + comp_r_3_18 + comp_r_3_19 + comp_r_3_20 + comp_r_3_21 + comp_r_3_22 + comp_r_3_23 + comp_r_3_24 + comp_r_3_25 + comp_r_3_26 + comp_r_3_27 + comp_r_3_28 + comp_r_3_29 + comp_r_3_30 + comp_r_3_31 + comp_r_3_32 + comp_r_3_33 + comp_r_3_34 + comp_r_3_35 + comp_r_3_36 + comp_r_3_37 + comp_r_3_38 + comp_r_3_39 + comp_r_3_40 + comp_r_3_41 + comp_r_3_42 + comp_r_3_43 + comp_r_3_44 + comp_r_3_45 + comp_r_3_46 + comp_r_3_47 + comp_r_3_48 + comp_r_3_49 + comp_r_3_50 + comp_r_3_51 + comp_r_3_52 + comp_r_3_53 + comp_r_3_54 + comp_r_3_55 + comp_r_3_56 + comp_r_3_57 + comp_r_3_58 + comp_r_3_59 + comp_r_3_60 + comp_r_3_61 + comp_r_3_62;
    position[4] = 2 + (sc_uint<1>) ~(comp_r_1_4) + (sc_uint<1>) ~(comp_r_3_4);
    position[5] = 3 + (sc_uint<1>) ~(comp_r_1_5) + (sc_uint<1>) ~(comp_r_3_5) + comp_r_5_6 + comp_r_5_7 + comp_r_5_8 + comp_r_5_9 + comp_r_5_10 + comp_r_5_11 + comp_r_5_12 + comp_r_5_13 + comp_r_5_14 + comp_r_5_15 + comp_r_5_16 + comp_r_5_17 + comp_r_5_18 + comp_r_5_19 + comp_r_5_20 + comp_r_5_21 + comp_r_5_22 + comp_r_5_23 + comp_r_5_24 + comp_r_5_25 + comp_r_5_26 + comp_r_5_27 + comp_r_5_28 + comp_r_5_29 + comp_r_5_30 + comp_r_5_31 + comp_r_5_32 + comp_r_5_33 + comp_r_5_34 + comp_r_5_35 + comp_r_5_36 + comp_r_5_37 + comp_r_5_38 + comp_r_5_39 + comp_r_5_40 + comp_r_5_41 + comp_r_5_42 + comp_r_5_43 + comp_r_5_44 + comp_r_5_45 + comp_r_5_46 + comp_r_5_47 + comp_r_5_48 + comp_r_5_49 + comp_r_5_50 + comp_r_5_51 + comp_r_5_52 + comp_r_5_53 + comp_r_5_54 + comp_r_5_55 + comp_r_5_56 + comp_r_5_57 + comp_r_5_58 + comp_r_5_59 + comp_r_5_60 + comp_r_5_61 + comp_r_5_62;
    position[6] = 3 + (sc_uint<1>) ~(comp_r_1_6) + (sc_uint<1>) ~(comp_r_3_6) + (sc_uint<1>) ~(comp_r_5_6);
    position[7] = 4 + (sc_uint<1>) ~(comp_r_1_7) + (sc_uint<1>) ~(comp_r_3_7) + (sc_uint<1>) ~(comp_r_5_7) + comp_r_7_8 + comp_r_7_9 + comp_r_7_10 + comp_r_7_11 + comp_r_7_12 + comp_r_7_13 + comp_r_7_14 + comp_r_7_15 + comp_r_7_16 + comp_r_7_17 + comp_r_7_18 + comp_r_7_19 + comp_r_7_20 + comp_r_7_21 + comp_r_7_22 + comp_r_7_23 + comp_r_7_24 + comp_r_7_25 + comp_r_7_26 + comp_r_7_27 + comp_r_7_28 + comp_r_7_29 + comp_r_7_30 + comp_r_7_31 + comp_r_7_32 + comp_r_7_33 + comp_r_7_34 + comp_r_7_35 + comp_r_7_36 + comp_r_7_37 + comp_r_7_38 + comp_r_7_39 + comp_r_7_40 + comp_r_7_41 + comp_r_7_42 + comp_r_7_43 + comp_r_7_44 + comp_r_7_45 + comp_r_7_46 + comp_r_7_47 + comp_r_7_48 + comp_r_7_49 + comp_r_7_50 + comp_r_7_51 + comp_r_7_52 + comp_r_7_53 + comp_r_7_54 + comp_r_7_55 + comp_r_7_56 + comp_r_7_57 + comp_r_7_58 + comp_r_7_59 + comp_r_7_60 + comp_r_7_61 + comp_r_7_62;
    position[8] = 4 + (sc_uint<1>) ~(comp_r_1_8) + (sc_uint<1>) ~(comp_r_3_8) + (sc_uint<1>) ~(comp_r_5_8) + (sc_uint<1>) ~(comp_r_7_8);
    position[9] = 5 + (sc_uint<1>) ~(comp_r_1_9) + (sc_uint<1>) ~(comp_r_3_9) + (sc_uint<1>) ~(comp_r_5_9) + (sc_uint<1>) ~(comp_r_7_9) + comp_r_9_10 + comp_r_9_11 + comp_r_9_12 + comp_r_9_13 + comp_r_9_14 + comp_r_9_15 + comp_r_9_16 + comp_r_9_17 + comp_r_9_18 + comp_r_9_19 + comp_r_9_20 + comp_r_9_21 + comp_r_9_22 + comp_r_9_23 + comp_r_9_24 + comp_r_9_25 + comp_r_9_26 + comp_r_9_27 + comp_r_9_28 + comp_r_9_29 + comp_r_9_30 + comp_r_9_31 + comp_r_9_32 + comp_r_9_33 + comp_r_9_34 + comp_r_9_35 + comp_r_9_36 + comp_r_9_37 + comp_r_9_38 + comp_r_9_39 + comp_r_9_40 + comp_r_9_41 + comp_r_9_42 + comp_r_9_43 + comp_r_9_44 + comp_r_9_45 + comp_r_9_46 + comp_r_9_47 + comp_r_9_48 + comp_r_9_49 + comp_r_9_50 + comp_r_9_51 + comp_r_9_52 + comp_r_9_53 + comp_r_9_54 + comp_r_9_55 + comp_r_9_56 + comp_r_9_57 + comp_r_9_58 + comp_r_9_59 + comp_r_9_60 + comp_r_9_61 + comp_r_9_62;
    position[10] = 5 + (sc_uint<1>) ~(comp_r_1_10) + (sc_uint<1>) ~(comp_r_3_10) + (sc_uint<1>) ~(comp_r_5_10) + (sc_uint<1>) ~(comp_r_7_10) + (sc_uint<1>) ~(comp_r_9_10);
    position[11] = 6 + (sc_uint<1>) ~(comp_r_1_11) + (sc_uint<1>) ~(comp_r_3_11) + (sc_uint<1>) ~(comp_r_5_11) + (sc_uint<1>) ~(comp_r_7_11) + (sc_uint<1>) ~(comp_r_9_11) + comp_r_11_12 + comp_r_11_13 + comp_r_11_14 + comp_r_11_15 + comp_r_11_16 + comp_r_11_17 + comp_r_11_18 + comp_r_11_19 + comp_r_11_20 + comp_r_11_21 + comp_r_11_22 + comp_r_11_23 + comp_r_11_24 + comp_r_11_25 + comp_r_11_26 + comp_r_11_27 + comp_r_11_28 + comp_r_11_29 + comp_r_11_30 + comp_r_11_31 + comp_r_11_32 + comp_r_11_33 + comp_r_11_34 + comp_r_11_35 + comp_r_11_36 + comp_r_11_37 + comp_r_11_38 + comp_r_11_39 + comp_r_11_40 + comp_r_11_41 + comp_r_11_42 + comp_r_11_43 + comp_r_11_44 + comp_r_11_45 + comp_r_11_46 + comp_r_11_47 + comp_r_11_48 + comp_r_11_49 + comp_r_11_50 + comp_r_11_51 + comp_r_11_52 + comp_r_11_53 + comp_r_11_54 + comp_r_11_55 + comp_r_11_56 + comp_r_11_57 + comp_r_11_58 + comp_r_11_59 + comp_r_11_60 + comp_r_11_61 + comp_r_11_62;
    position[12] = 6 + (sc_uint<1>) ~(comp_r_1_12) + (sc_uint<1>) ~(comp_r_3_12) + (sc_uint<1>) ~(comp_r_5_12) + (sc_uint<1>) ~(comp_r_7_12) + (sc_uint<1>) ~(comp_r_9_12) + (sc_uint<1>) ~(comp_r_11_12);
    position[13] = 7 + (sc_uint<1>) ~(comp_r_1_13) + (sc_uint<1>) ~(comp_r_3_13) + (sc_uint<1>) ~(comp_r_5_13) + (sc_uint<1>) ~(comp_r_7_13) + (sc_uint<1>) ~(comp_r_9_13) + (sc_uint<1>) ~(comp_r_11_13) + comp_r_13_14 + comp_r_13_15 + comp_r_13_16 + comp_r_13_17 + comp_r_13_18 + comp_r_13_19 + comp_r_13_20 + comp_r_13_21 + comp_r_13_22 + comp_r_13_23 + comp_r_13_24 + comp_r_13_25 + comp_r_13_26 + comp_r_13_27 + comp_r_13_28 + comp_r_13_29 + comp_r_13_30 + comp_r_13_31 + comp_r_13_32 + comp_r_13_33 + comp_r_13_34 + comp_r_13_35 + comp_r_13_36 + comp_r_13_37 + comp_r_13_38 + comp_r_13_39 + comp_r_13_40 + comp_r_13_41 + comp_r_13_42 + comp_r_13_43 + comp_r_13_44 + comp_r_13_45 + comp_r_13_46 + comp_r_13_47 + comp_r_13_48 + comp_r_13_49 + comp_r_13_50 + comp_r_13_51 + comp_r_13_52 + comp_r_13_53 + comp_r_13_54 + comp_r_13_55 + comp_r_13_56 + comp_r_13_57 + comp_r_13_58 + comp_r_13_59 + comp_r_13_60 + comp_r_13_61 + comp_r_13_62;
    position[14] = 7 + (sc_uint<1>) ~(comp_r_1_14) + (sc_uint<1>) ~(comp_r_3_14) + (sc_uint<1>) ~(comp_r_5_14) + (sc_uint<1>) ~(comp_r_7_14) + (sc_uint<1>) ~(comp_r_9_14) + (sc_uint<1>) ~(comp_r_11_14) + (sc_uint<1>) ~(comp_r_13_14);
    position[15] = 8 + (sc_uint<1>) ~(comp_r_1_15) + (sc_uint<1>) ~(comp_r_3_15) + (sc_uint<1>) ~(comp_r_5_15) + (sc_uint<1>) ~(comp_r_7_15) + (sc_uint<1>) ~(comp_r_9_15) + (sc_uint<1>) ~(comp_r_11_15) + (sc_uint<1>) ~(comp_r_13_15) + comp_r_15_16 + comp_r_15_17 + comp_r_15_18 + comp_r_15_19 + comp_r_15_20 + comp_r_15_21 + comp_r_15_22 + comp_r_15_23 + comp_r_15_24 + comp_r_15_25 + comp_r_15_26 + comp_r_15_27 + comp_r_15_28 + comp_r_15_29 + comp_r_15_30 + comp_r_15_31 + comp_r_15_32 + comp_r_15_33 + comp_r_15_34 + comp_r_15_35 + comp_r_15_36 + comp_r_15_37 + comp_r_15_38 + comp_r_15_39 + comp_r_15_40 + comp_r_15_41 + comp_r_15_42 + comp_r_15_43 + comp_r_15_44 + comp_r_15_45 + comp_r_15_46 + comp_r_15_47 + comp_r_15_48 + comp_r_15_49 + comp_r_15_50 + comp_r_15_51 + comp_r_15_52 + comp_r_15_53 + comp_r_15_54 + comp_r_15_55 + comp_r_15_56 + comp_r_15_57 + comp_r_15_58 + comp_r_15_59 + comp_r_15_60 + comp_r_15_61 + comp_r_15_62;
    position[16] = 8 + (sc_uint<1>) ~(comp_r_1_16) + (sc_uint<1>) ~(comp_r_3_16) + (sc_uint<1>) ~(comp_r_5_16) + (sc_uint<1>) ~(comp_r_7_16) + (sc_uint<1>) ~(comp_r_9_16) + (sc_uint<1>) ~(comp_r_11_16) + (sc_uint<1>) ~(comp_r_13_16) + (sc_uint<1>) ~(comp_r_15_16);
    position[17] = 9 + (sc_uint<1>) ~(comp_r_1_17) + (sc_uint<1>) ~(comp_r_3_17) + (sc_uint<1>) ~(comp_r_5_17) + (sc_uint<1>) ~(comp_r_7_17) + (sc_uint<1>) ~(comp_r_9_17) + (sc_uint<1>) ~(comp_r_11_17) + (sc_uint<1>) ~(comp_r_13_17) + (sc_uint<1>) ~(comp_r_15_17) + comp_r_17_18 + comp_r_17_19 + comp_r_17_20 + comp_r_17_21 + comp_r_17_22 + comp_r_17_23 + comp_r_17_24 + comp_r_17_25 + comp_r_17_26 + comp_r_17_27 + comp_r_17_28 + comp_r_17_29 + comp_r_17_30 + comp_r_17_31 + comp_r_17_32 + comp_r_17_33 + comp_r_17_34 + comp_r_17_35 + comp_r_17_36 + comp_r_17_37 + comp_r_17_38 + comp_r_17_39 + comp_r_17_40 + comp_r_17_41 + comp_r_17_42 + comp_r_17_43 + comp_r_17_44 + comp_r_17_45 + comp_r_17_46 + comp_r_17_47 + comp_r_17_48 + comp_r_17_49 + comp_r_17_50 + comp_r_17_51 + comp_r_17_52 + comp_r_17_53 + comp_r_17_54 + comp_r_17_55 + comp_r_17_56 + comp_r_17_57 + comp_r_17_58 + comp_r_17_59 + comp_r_17_60 + comp_r_17_61 + comp_r_17_62;
    position[18] = 9 + (sc_uint<1>) ~(comp_r_1_18) + (sc_uint<1>) ~(comp_r_3_18) + (sc_uint<1>) ~(comp_r_5_18) + (sc_uint<1>) ~(comp_r_7_18) + (sc_uint<1>) ~(comp_r_9_18) + (sc_uint<1>) ~(comp_r_11_18) + (sc_uint<1>) ~(comp_r_13_18) + (sc_uint<1>) ~(comp_r_15_18) + (sc_uint<1>) ~(comp_r_17_18);
    position[19] = 10 + (sc_uint<1>) ~(comp_r_1_19) + (sc_uint<1>) ~(comp_r_3_19) + (sc_uint<1>) ~(comp_r_5_19) + (sc_uint<1>) ~(comp_r_7_19) + (sc_uint<1>) ~(comp_r_9_19) + (sc_uint<1>) ~(comp_r_11_19) + (sc_uint<1>) ~(comp_r_13_19) + (sc_uint<1>) ~(comp_r_15_19) + (sc_uint<1>) ~(comp_r_17_19) + comp_r_19_20 + comp_r_19_21 + comp_r_19_22 + comp_r_19_23 + comp_r_19_24 + comp_r_19_25 + comp_r_19_26 + comp_r_19_27 + comp_r_19_28 + comp_r_19_29 + comp_r_19_30 + comp_r_19_31 + comp_r_19_32 + comp_r_19_33 + comp_r_19_34 + comp_r_19_35 + comp_r_19_36 + comp_r_19_37 + comp_r_19_38 + comp_r_19_39 + comp_r_19_40 + comp_r_19_41 + comp_r_19_42 + comp_r_19_43 + comp_r_19_44 + comp_r_19_45 + comp_r_19_46 + comp_r_19_47 + comp_r_19_48 + comp_r_19_49 + comp_r_19_50 + comp_r_19_51 + comp_r_19_52 + comp_r_19_53 + comp_r_19_54 + comp_r_19_55 + comp_r_19_56 + comp_r_19_57 + comp_r_19_58 + comp_r_19_59 + comp_r_19_60 + comp_r_19_61 + comp_r_19_62;
    position[20] = 10 + (sc_uint<1>) ~(comp_r_1_20) + (sc_uint<1>) ~(comp_r_3_20) + (sc_uint<1>) ~(comp_r_5_20) + (sc_uint<1>) ~(comp_r_7_20) + (sc_uint<1>) ~(comp_r_9_20) + (sc_uint<1>) ~(comp_r_11_20) + (sc_uint<1>) ~(comp_r_13_20) + (sc_uint<1>) ~(comp_r_15_20) + (sc_uint<1>) ~(comp_r_17_20) + (sc_uint<1>) ~(comp_r_19_20);
    position[21] = 11 + (sc_uint<1>) ~(comp_r_1_21) + (sc_uint<1>) ~(comp_r_3_21) + (sc_uint<1>) ~(comp_r_5_21) + (sc_uint<1>) ~(comp_r_7_21) + (sc_uint<1>) ~(comp_r_9_21) + (sc_uint<1>) ~(comp_r_11_21) + (sc_uint<1>) ~(comp_r_13_21) + (sc_uint<1>) ~(comp_r_15_21) + (sc_uint<1>) ~(comp_r_17_21) + (sc_uint<1>) ~(comp_r_19_21) + comp_r_21_22 + comp_r_21_23 + comp_r_21_24 + comp_r_21_25 + comp_r_21_26 + comp_r_21_27 + comp_r_21_28 + comp_r_21_29 + comp_r_21_30 + comp_r_21_31 + comp_r_21_32 + comp_r_21_33 + comp_r_21_34 + comp_r_21_35 + comp_r_21_36 + comp_r_21_37 + comp_r_21_38 + comp_r_21_39 + comp_r_21_40 + comp_r_21_41 + comp_r_21_42 + comp_r_21_43 + comp_r_21_44 + comp_r_21_45 + comp_r_21_46 + comp_r_21_47 + comp_r_21_48 + comp_r_21_49 + comp_r_21_50 + comp_r_21_51 + comp_r_21_52 + comp_r_21_53 + comp_r_21_54 + comp_r_21_55 + comp_r_21_56 + comp_r_21_57 + comp_r_21_58 + comp_r_21_59 + comp_r_21_60 + comp_r_21_61 + comp_r_21_62;
    position[22] = 11 + (sc_uint<1>) ~(comp_r_1_22) + (sc_uint<1>) ~(comp_r_3_22) + (sc_uint<1>) ~(comp_r_5_22) + (sc_uint<1>) ~(comp_r_7_22) + (sc_uint<1>) ~(comp_r_9_22) + (sc_uint<1>) ~(comp_r_11_22) + (sc_uint<1>) ~(comp_r_13_22) + (sc_uint<1>) ~(comp_r_15_22) + (sc_uint<1>) ~(comp_r_17_22) + (sc_uint<1>) ~(comp_r_19_22) + (sc_uint<1>) ~(comp_r_21_22);
    position[23] = 12 + (sc_uint<1>) ~(comp_r_1_23) + (sc_uint<1>) ~(comp_r_3_23) + (sc_uint<1>) ~(comp_r_5_23) + (sc_uint<1>) ~(comp_r_7_23) + (sc_uint<1>) ~(comp_r_9_23) + (sc_uint<1>) ~(comp_r_11_23) + (sc_uint<1>) ~(comp_r_13_23) + (sc_uint<1>) ~(comp_r_15_23) + (sc_uint<1>) ~(comp_r_17_23) + (sc_uint<1>) ~(comp_r_19_23) + (sc_uint<1>) ~(comp_r_21_23) + comp_r_23_24 + comp_r_23_25 + comp_r_23_26 + comp_r_23_27 + comp_r_23_28 + comp_r_23_29 + comp_r_23_30 + comp_r_23_31 + comp_r_23_32 + comp_r_23_33 + comp_r_23_34 + comp_r_23_35 + comp_r_23_36 + comp_r_23_37 + comp_r_23_38 + comp_r_23_39 + comp_r_23_40 + comp_r_23_41 + comp_r_23_42 + comp_r_23_43 + comp_r_23_44 + comp_r_23_45 + comp_r_23_46 + comp_r_23_47 + comp_r_23_48 + comp_r_23_49 + comp_r_23_50 + comp_r_23_51 + comp_r_23_52 + comp_r_23_53 + comp_r_23_54 + comp_r_23_55 + comp_r_23_56 + comp_r_23_57 + comp_r_23_58 + comp_r_23_59 + comp_r_23_60 + comp_r_23_61 + comp_r_23_62;
    position[24] = 12 + (sc_uint<1>) ~(comp_r_1_24) + (sc_uint<1>) ~(comp_r_3_24) + (sc_uint<1>) ~(comp_r_5_24) + (sc_uint<1>) ~(comp_r_7_24) + (sc_uint<1>) ~(comp_r_9_24) + (sc_uint<1>) ~(comp_r_11_24) + (sc_uint<1>) ~(comp_r_13_24) + (sc_uint<1>) ~(comp_r_15_24) + (sc_uint<1>) ~(comp_r_17_24) + (sc_uint<1>) ~(comp_r_19_24) + (sc_uint<1>) ~(comp_r_21_24) + (sc_uint<1>) ~(comp_r_23_24);
    position[25] = 13 + (sc_uint<1>) ~(comp_r_1_25) + (sc_uint<1>) ~(comp_r_3_25) + (sc_uint<1>) ~(comp_r_5_25) + (sc_uint<1>) ~(comp_r_7_25) + (sc_uint<1>) ~(comp_r_9_25) + (sc_uint<1>) ~(comp_r_11_25) + (sc_uint<1>) ~(comp_r_13_25) + (sc_uint<1>) ~(comp_r_15_25) + (sc_uint<1>) ~(comp_r_17_25) + (sc_uint<1>) ~(comp_r_19_25) + (sc_uint<1>) ~(comp_r_21_25) + (sc_uint<1>) ~(comp_r_23_25) + comp_r_25_26 + comp_r_25_27 + comp_r_25_28 + comp_r_25_29 + comp_r_25_30 + comp_r_25_31 + comp_r_25_32 + comp_r_25_33 + comp_r_25_34 + comp_r_25_35 + comp_r_25_36 + comp_r_25_37 + comp_r_25_38 + comp_r_25_39 + comp_r_25_40 + comp_r_25_41 + comp_r_25_42 + comp_r_25_43 + comp_r_25_44 + comp_r_25_45 + comp_r_25_46 + comp_r_25_47 + comp_r_25_48 + comp_r_25_49 + comp_r_25_50 + comp_r_25_51 + comp_r_25_52 + comp_r_25_53 + comp_r_25_54 + comp_r_25_55 + comp_r_25_56 + comp_r_25_57 + comp_r_25_58 + comp_r_25_59 + comp_r_25_60 + comp_r_25_61 + comp_r_25_62;
    position[26] = 13 + (sc_uint<1>) ~(comp_r_1_26) + (sc_uint<1>) ~(comp_r_3_26) + (sc_uint<1>) ~(comp_r_5_26) + (sc_uint<1>) ~(comp_r_7_26) + (sc_uint<1>) ~(comp_r_9_26) + (sc_uint<1>) ~(comp_r_11_26) + (sc_uint<1>) ~(comp_r_13_26) + (sc_uint<1>) ~(comp_r_15_26) + (sc_uint<1>) ~(comp_r_17_26) + (sc_uint<1>) ~(comp_r_19_26) + (sc_uint<1>) ~(comp_r_21_26) + (sc_uint<1>) ~(comp_r_23_26) + (sc_uint<1>) ~(comp_r_25_26);
    position[27] = 14 + (sc_uint<1>) ~(comp_r_1_27) + (sc_uint<1>) ~(comp_r_3_27) + (sc_uint<1>) ~(comp_r_5_27) + (sc_uint<1>) ~(comp_r_7_27) + (sc_uint<1>) ~(comp_r_9_27) + (sc_uint<1>) ~(comp_r_11_27) + (sc_uint<1>) ~(comp_r_13_27) + (sc_uint<1>) ~(comp_r_15_27) + (sc_uint<1>) ~(comp_r_17_27) + (sc_uint<1>) ~(comp_r_19_27) + (sc_uint<1>) ~(comp_r_21_27) + (sc_uint<1>) ~(comp_r_23_27) + (sc_uint<1>) ~(comp_r_25_27) + comp_r_27_28 + comp_r_27_29 + comp_r_27_30 + comp_r_27_31 + comp_r_27_32 + comp_r_27_33 + comp_r_27_34 + comp_r_27_35 + comp_r_27_36 + comp_r_27_37 + comp_r_27_38 + comp_r_27_39 + comp_r_27_40 + comp_r_27_41 + comp_r_27_42 + comp_r_27_43 + comp_r_27_44 + comp_r_27_45 + comp_r_27_46 + comp_r_27_47 + comp_r_27_48 + comp_r_27_49 + comp_r_27_50 + comp_r_27_51 + comp_r_27_52 + comp_r_27_53 + comp_r_27_54 + comp_r_27_55 + comp_r_27_56 + comp_r_27_57 + comp_r_27_58 + comp_r_27_59 + comp_r_27_60 + comp_r_27_61 + comp_r_27_62;
    position[28] = 14 + (sc_uint<1>) ~(comp_r_1_28) + (sc_uint<1>) ~(comp_r_3_28) + (sc_uint<1>) ~(comp_r_5_28) + (sc_uint<1>) ~(comp_r_7_28) + (sc_uint<1>) ~(comp_r_9_28) + (sc_uint<1>) ~(comp_r_11_28) + (sc_uint<1>) ~(comp_r_13_28) + (sc_uint<1>) ~(comp_r_15_28) + (sc_uint<1>) ~(comp_r_17_28) + (sc_uint<1>) ~(comp_r_19_28) + (sc_uint<1>) ~(comp_r_21_28) + (sc_uint<1>) ~(comp_r_23_28) + (sc_uint<1>) ~(comp_r_25_28) + (sc_uint<1>) ~(comp_r_27_28);
    position[29] = 15 + (sc_uint<1>) ~(comp_r_1_29) + (sc_uint<1>) ~(comp_r_3_29) + (sc_uint<1>) ~(comp_r_5_29) + (sc_uint<1>) ~(comp_r_7_29) + (sc_uint<1>) ~(comp_r_9_29) + (sc_uint<1>) ~(comp_r_11_29) + (sc_uint<1>) ~(comp_r_13_29) + (sc_uint<1>) ~(comp_r_15_29) + (sc_uint<1>) ~(comp_r_17_29) + (sc_uint<1>) ~(comp_r_19_29) + (sc_uint<1>) ~(comp_r_21_29) + (sc_uint<1>) ~(comp_r_23_29) + (sc_uint<1>) ~(comp_r_25_29) + (sc_uint<1>) ~(comp_r_27_29) + comp_r_29_30 + comp_r_29_31 + comp_r_29_32 + comp_r_29_33 + comp_r_29_34 + comp_r_29_35 + comp_r_29_36 + comp_r_29_37 + comp_r_29_38 + comp_r_29_39 + comp_r_29_40 + comp_r_29_41 + comp_r_29_42 + comp_r_29_43 + comp_r_29_44 + comp_r_29_45 + comp_r_29_46 + comp_r_29_47 + comp_r_29_48 + comp_r_29_49 + comp_r_29_50 + comp_r_29_51 + comp_r_29_52 + comp_r_29_53 + comp_r_29_54 + comp_r_29_55 + comp_r_29_56 + comp_r_29_57 + comp_r_29_58 + comp_r_29_59 + comp_r_29_60 + comp_r_29_61 + comp_r_29_62;
    position[30] = 15 + (sc_uint<1>) ~(comp_r_1_30) + (sc_uint<1>) ~(comp_r_3_30) + (sc_uint<1>) ~(comp_r_5_30) + (sc_uint<1>) ~(comp_r_7_30) + (sc_uint<1>) ~(comp_r_9_30) + (sc_uint<1>) ~(comp_r_11_30) + (sc_uint<1>) ~(comp_r_13_30) + (sc_uint<1>) ~(comp_r_15_30) + (sc_uint<1>) ~(comp_r_17_30) + (sc_uint<1>) ~(comp_r_19_30) + (sc_uint<1>) ~(comp_r_21_30) + (sc_uint<1>) ~(comp_r_23_30) + (sc_uint<1>) ~(comp_r_25_30) + (sc_uint<1>) ~(comp_r_27_30) + (sc_uint<1>) ~(comp_r_29_30);
    position[31] = 16 + (sc_uint<1>) ~(comp_r_1_31) + (sc_uint<1>) ~(comp_r_3_31) + (sc_uint<1>) ~(comp_r_5_31) + (sc_uint<1>) ~(comp_r_7_31) + (sc_uint<1>) ~(comp_r_9_31) + (sc_uint<1>) ~(comp_r_11_31) + (sc_uint<1>) ~(comp_r_13_31) + (sc_uint<1>) ~(comp_r_15_31) + (sc_uint<1>) ~(comp_r_17_31) + (sc_uint<1>) ~(comp_r_19_31) + (sc_uint<1>) ~(comp_r_21_31) + (sc_uint<1>) ~(comp_r_23_31) + (sc_uint<1>) ~(comp_r_25_31) + (sc_uint<1>) ~(comp_r_27_31) + (sc_uint<1>) ~(comp_r_29_31) + comp_r_31_32 + comp_r_31_33 + comp_r_31_34 + comp_r_31_35 + comp_r_31_36 + comp_r_31_37 + comp_r_31_38 + comp_r_31_39 + comp_r_31_40 + comp_r_31_41 + comp_r_31_42 + comp_r_31_43 + comp_r_31_44 + comp_r_31_45 + comp_r_31_46 + comp_r_31_47 + comp_r_31_48 + comp_r_31_49 + comp_r_31_50 + comp_r_31_51 + comp_r_31_52 + comp_r_31_53 + comp_r_31_54 + comp_r_31_55 + comp_r_31_56 + comp_r_31_57 + comp_r_31_58 + comp_r_31_59 + comp_r_31_60 + comp_r_31_61 + comp_r_31_62;
    position[32] = 16 + (sc_uint<1>) ~(comp_r_1_32) + (sc_uint<1>) ~(comp_r_3_32) + (sc_uint<1>) ~(comp_r_5_32) + (sc_uint<1>) ~(comp_r_7_32) + (sc_uint<1>) ~(comp_r_9_32) + (sc_uint<1>) ~(comp_r_11_32) + (sc_uint<1>) ~(comp_r_13_32) + (sc_uint<1>) ~(comp_r_15_32) + (sc_uint<1>) ~(comp_r_17_32) + (sc_uint<1>) ~(comp_r_19_32) + (sc_uint<1>) ~(comp_r_21_32) + (sc_uint<1>) ~(comp_r_23_32) + (sc_uint<1>) ~(comp_r_25_32) + (sc_uint<1>) ~(comp_r_27_32) + (sc_uint<1>) ~(comp_r_29_32) + (sc_uint<1>) ~(comp_r_31_32);
    position[33] = 17 + (sc_uint<1>) ~(comp_r_1_33) + (sc_uint<1>) ~(comp_r_3_33) + (sc_uint<1>) ~(comp_r_5_33) + (sc_uint<1>) ~(comp_r_7_33) + (sc_uint<1>) ~(comp_r_9_33) + (sc_uint<1>) ~(comp_r_11_33) + (sc_uint<1>) ~(comp_r_13_33) + (sc_uint<1>) ~(comp_r_15_33) + (sc_uint<1>) ~(comp_r_17_33) + (sc_uint<1>) ~(comp_r_19_33) + (sc_uint<1>) ~(comp_r_21_33) + (sc_uint<1>) ~(comp_r_23_33) + (sc_uint<1>) ~(comp_r_25_33) + (sc_uint<1>) ~(comp_r_27_33) + (sc_uint<1>) ~(comp_r_29_33) + (sc_uint<1>) ~(comp_r_31_33) + comp_r_33_34 + comp_r_33_35 + comp_r_33_36 + comp_r_33_37 + comp_r_33_38 + comp_r_33_39 + comp_r_33_40 + comp_r_33_41 + comp_r_33_42 + comp_r_33_43 + comp_r_33_44 + comp_r_33_45 + comp_r_33_46 + comp_r_33_47 + comp_r_33_48 + comp_r_33_49 + comp_r_33_50 + comp_r_33_51 + comp_r_33_52 + comp_r_33_53 + comp_r_33_54 + comp_r_33_55 + comp_r_33_56 + comp_r_33_57 + comp_r_33_58 + comp_r_33_59 + comp_r_33_60 + comp_r_33_61 + comp_r_33_62;
    position[34] = 17 + (sc_uint<1>) ~(comp_r_1_34) + (sc_uint<1>) ~(comp_r_3_34) + (sc_uint<1>) ~(comp_r_5_34) + (sc_uint<1>) ~(comp_r_7_34) + (sc_uint<1>) ~(comp_r_9_34) + (sc_uint<1>) ~(comp_r_11_34) + (sc_uint<1>) ~(comp_r_13_34) + (sc_uint<1>) ~(comp_r_15_34) + (sc_uint<1>) ~(comp_r_17_34) + (sc_uint<1>) ~(comp_r_19_34) + (sc_uint<1>) ~(comp_r_21_34) + (sc_uint<1>) ~(comp_r_23_34) + (sc_uint<1>) ~(comp_r_25_34) + (sc_uint<1>) ~(comp_r_27_34) + (sc_uint<1>) ~(comp_r_29_34) + (sc_uint<1>) ~(comp_r_31_34) + (sc_uint<1>) ~(comp_r_33_34);
    position[35] = 18 + (sc_uint<1>) ~(comp_r_1_35) + (sc_uint<1>) ~(comp_r_3_35) + (sc_uint<1>) ~(comp_r_5_35) + (sc_uint<1>) ~(comp_r_7_35) + (sc_uint<1>) ~(comp_r_9_35) + (sc_uint<1>) ~(comp_r_11_35) + (sc_uint<1>) ~(comp_r_13_35) + (sc_uint<1>) ~(comp_r_15_35) + (sc_uint<1>) ~(comp_r_17_35) + (sc_uint<1>) ~(comp_r_19_35) + (sc_uint<1>) ~(comp_r_21_35) + (sc_uint<1>) ~(comp_r_23_35) + (sc_uint<1>) ~(comp_r_25_35) + (sc_uint<1>) ~(comp_r_27_35) + (sc_uint<1>) ~(comp_r_29_35) + (sc_uint<1>) ~(comp_r_31_35) + (sc_uint<1>) ~(comp_r_33_35) + comp_r_35_36 + comp_r_35_37 + comp_r_35_38 + comp_r_35_39 + comp_r_35_40 + comp_r_35_41 + comp_r_35_42 + comp_r_35_43 + comp_r_35_44 + comp_r_35_45 + comp_r_35_46 + comp_r_35_47 + comp_r_35_48 + comp_r_35_49 + comp_r_35_50 + comp_r_35_51 + comp_r_35_52 + comp_r_35_53 + comp_r_35_54 + comp_r_35_55 + comp_r_35_56 + comp_r_35_57 + comp_r_35_58 + comp_r_35_59 + comp_r_35_60 + comp_r_35_61 + comp_r_35_62;
    position[36] = 18 + (sc_uint<1>) ~(comp_r_1_36) + (sc_uint<1>) ~(comp_r_3_36) + (sc_uint<1>) ~(comp_r_5_36) + (sc_uint<1>) ~(comp_r_7_36) + (sc_uint<1>) ~(comp_r_9_36) + (sc_uint<1>) ~(comp_r_11_36) + (sc_uint<1>) ~(comp_r_13_36) + (sc_uint<1>) ~(comp_r_15_36) + (sc_uint<1>) ~(comp_r_17_36) + (sc_uint<1>) ~(comp_r_19_36) + (sc_uint<1>) ~(comp_r_21_36) + (sc_uint<1>) ~(comp_r_23_36) + (sc_uint<1>) ~(comp_r_25_36) + (sc_uint<1>) ~(comp_r_27_36) + (sc_uint<1>) ~(comp_r_29_36) + (sc_uint<1>) ~(comp_r_31_36) + (sc_uint<1>) ~(comp_r_33_36) + (sc_uint<1>) ~(comp_r_35_36);
    position[37] = 19 + (sc_uint<1>) ~(comp_r_1_37) + (sc_uint<1>) ~(comp_r_3_37) + (sc_uint<1>) ~(comp_r_5_37) + (sc_uint<1>) ~(comp_r_7_37) + (sc_uint<1>) ~(comp_r_9_37) + (sc_uint<1>) ~(comp_r_11_37) + (sc_uint<1>) ~(comp_r_13_37) + (sc_uint<1>) ~(comp_r_15_37) + (sc_uint<1>) ~(comp_r_17_37) + (sc_uint<1>) ~(comp_r_19_37) + (sc_uint<1>) ~(comp_r_21_37) + (sc_uint<1>) ~(comp_r_23_37) + (sc_uint<1>) ~(comp_r_25_37) + (sc_uint<1>) ~(comp_r_27_37) + (sc_uint<1>) ~(comp_r_29_37) + (sc_uint<1>) ~(comp_r_31_37) + (sc_uint<1>) ~(comp_r_33_37) + (sc_uint<1>) ~(comp_r_35_37) + comp_r_37_38 + comp_r_37_39 + comp_r_37_40 + comp_r_37_41 + comp_r_37_42 + comp_r_37_43 + comp_r_37_44 + comp_r_37_45 + comp_r_37_46 + comp_r_37_47 + comp_r_37_48 + comp_r_37_49 + comp_r_37_50 + comp_r_37_51 + comp_r_37_52 + comp_r_37_53 + comp_r_37_54 + comp_r_37_55 + comp_r_37_56 + comp_r_37_57 + comp_r_37_58 + comp_r_37_59 + comp_r_37_60 + comp_r_37_61 + comp_r_37_62;
    position[38] = 19 + (sc_uint<1>) ~(comp_r_1_38) + (sc_uint<1>) ~(comp_r_3_38) + (sc_uint<1>) ~(comp_r_5_38) + (sc_uint<1>) ~(comp_r_7_38) + (sc_uint<1>) ~(comp_r_9_38) + (sc_uint<1>) ~(comp_r_11_38) + (sc_uint<1>) ~(comp_r_13_38) + (sc_uint<1>) ~(comp_r_15_38) + (sc_uint<1>) ~(comp_r_17_38) + (sc_uint<1>) ~(comp_r_19_38) + (sc_uint<1>) ~(comp_r_21_38) + (sc_uint<1>) ~(comp_r_23_38) + (sc_uint<1>) ~(comp_r_25_38) + (sc_uint<1>) ~(comp_r_27_38) + (sc_uint<1>) ~(comp_r_29_38) + (sc_uint<1>) ~(comp_r_31_38) + (sc_uint<1>) ~(comp_r_33_38) + (sc_uint<1>) ~(comp_r_35_38) + (sc_uint<1>) ~(comp_r_37_38);
    position[39] = 20 + (sc_uint<1>) ~(comp_r_1_39) + (sc_uint<1>) ~(comp_r_3_39) + (sc_uint<1>) ~(comp_r_5_39) + (sc_uint<1>) ~(comp_r_7_39) + (sc_uint<1>) ~(comp_r_9_39) + (sc_uint<1>) ~(comp_r_11_39) + (sc_uint<1>) ~(comp_r_13_39) + (sc_uint<1>) ~(comp_r_15_39) + (sc_uint<1>) ~(comp_r_17_39) + (sc_uint<1>) ~(comp_r_19_39) + (sc_uint<1>) ~(comp_r_21_39) + (sc_uint<1>) ~(comp_r_23_39) + (sc_uint<1>) ~(comp_r_25_39) + (sc_uint<1>) ~(comp_r_27_39) + (sc_uint<1>) ~(comp_r_29_39) + (sc_uint<1>) ~(comp_r_31_39) + (sc_uint<1>) ~(comp_r_33_39) + (sc_uint<1>) ~(comp_r_35_39) + (sc_uint<1>) ~(comp_r_37_39) + comp_r_39_40 + comp_r_39_41 + comp_r_39_42 + comp_r_39_43 + comp_r_39_44 + comp_r_39_45 + comp_r_39_46 + comp_r_39_47 + comp_r_39_48 + comp_r_39_49 + comp_r_39_50 + comp_r_39_51 + comp_r_39_52 + comp_r_39_53 + comp_r_39_54 + comp_r_39_55 + comp_r_39_56 + comp_r_39_57 + comp_r_39_58 + comp_r_39_59 + comp_r_39_60 + comp_r_39_61 + comp_r_39_62;
    position[40] = 20 + (sc_uint<1>) ~(comp_r_1_40) + (sc_uint<1>) ~(comp_r_3_40) + (sc_uint<1>) ~(comp_r_5_40) + (sc_uint<1>) ~(comp_r_7_40) + (sc_uint<1>) ~(comp_r_9_40) + (sc_uint<1>) ~(comp_r_11_40) + (sc_uint<1>) ~(comp_r_13_40) + (sc_uint<1>) ~(comp_r_15_40) + (sc_uint<1>) ~(comp_r_17_40) + (sc_uint<1>) ~(comp_r_19_40) + (sc_uint<1>) ~(comp_r_21_40) + (sc_uint<1>) ~(comp_r_23_40) + (sc_uint<1>) ~(comp_r_25_40) + (sc_uint<1>) ~(comp_r_27_40) + (sc_uint<1>) ~(comp_r_29_40) + (sc_uint<1>) ~(comp_r_31_40) + (sc_uint<1>) ~(comp_r_33_40) + (sc_uint<1>) ~(comp_r_35_40) + (sc_uint<1>) ~(comp_r_37_40) + (sc_uint<1>) ~(comp_r_39_40);
    position[41] = 21 + (sc_uint<1>) ~(comp_r_1_41) + (sc_uint<1>) ~(comp_r_3_41) + (sc_uint<1>) ~(comp_r_5_41) + (sc_uint<1>) ~(comp_r_7_41) + (sc_uint<1>) ~(comp_r_9_41) + (sc_uint<1>) ~(comp_r_11_41) + (sc_uint<1>) ~(comp_r_13_41) + (sc_uint<1>) ~(comp_r_15_41) + (sc_uint<1>) ~(comp_r_17_41) + (sc_uint<1>) ~(comp_r_19_41) + (sc_uint<1>) ~(comp_r_21_41) + (sc_uint<1>) ~(comp_r_23_41) + (sc_uint<1>) ~(comp_r_25_41) + (sc_uint<1>) ~(comp_r_27_41) + (sc_uint<1>) ~(comp_r_29_41) + (sc_uint<1>) ~(comp_r_31_41) + (sc_uint<1>) ~(comp_r_33_41) + (sc_uint<1>) ~(comp_r_35_41) + (sc_uint<1>) ~(comp_r_37_41) + (sc_uint<1>) ~(comp_r_39_41) + comp_r_41_42 + comp_r_41_43 + comp_r_41_44 + comp_r_41_45 + comp_r_41_46 + comp_r_41_47 + comp_r_41_48 + comp_r_41_49 + comp_r_41_50 + comp_r_41_51 + comp_r_41_52 + comp_r_41_53 + comp_r_41_54 + comp_r_41_55 + comp_r_41_56 + comp_r_41_57 + comp_r_41_58 + comp_r_41_59 + comp_r_41_60 + comp_r_41_61 + comp_r_41_62;
    position[42] = 21 + (sc_uint<1>) ~(comp_r_1_42) + (sc_uint<1>) ~(comp_r_3_42) + (sc_uint<1>) ~(comp_r_5_42) + (sc_uint<1>) ~(comp_r_7_42) + (sc_uint<1>) ~(comp_r_9_42) + (sc_uint<1>) ~(comp_r_11_42) + (sc_uint<1>) ~(comp_r_13_42) + (sc_uint<1>) ~(comp_r_15_42) + (sc_uint<1>) ~(comp_r_17_42) + (sc_uint<1>) ~(comp_r_19_42) + (sc_uint<1>) ~(comp_r_21_42) + (sc_uint<1>) ~(comp_r_23_42) + (sc_uint<1>) ~(comp_r_25_42) + (sc_uint<1>) ~(comp_r_27_42) + (sc_uint<1>) ~(comp_r_29_42) + (sc_uint<1>) ~(comp_r_31_42) + (sc_uint<1>) ~(comp_r_33_42) + (sc_uint<1>) ~(comp_r_35_42) + (sc_uint<1>) ~(comp_r_37_42) + (sc_uint<1>) ~(comp_r_39_42) + (sc_uint<1>) ~(comp_r_41_42);
    position[43] = 22 + (sc_uint<1>) ~(comp_r_1_43) + (sc_uint<1>) ~(comp_r_3_43) + (sc_uint<1>) ~(comp_r_5_43) + (sc_uint<1>) ~(comp_r_7_43) + (sc_uint<1>) ~(comp_r_9_43) + (sc_uint<1>) ~(comp_r_11_43) + (sc_uint<1>) ~(comp_r_13_43) + (sc_uint<1>) ~(comp_r_15_43) + (sc_uint<1>) ~(comp_r_17_43) + (sc_uint<1>) ~(comp_r_19_43) + (sc_uint<1>) ~(comp_r_21_43) + (sc_uint<1>) ~(comp_r_23_43) + (sc_uint<1>) ~(comp_r_25_43) + (sc_uint<1>) ~(comp_r_27_43) + (sc_uint<1>) ~(comp_r_29_43) + (sc_uint<1>) ~(comp_r_31_43) + (sc_uint<1>) ~(comp_r_33_43) + (sc_uint<1>) ~(comp_r_35_43) + (sc_uint<1>) ~(comp_r_37_43) + (sc_uint<1>) ~(comp_r_39_43) + (sc_uint<1>) ~(comp_r_41_43) + comp_r_43_44 + comp_r_43_45 + comp_r_43_46 + comp_r_43_47 + comp_r_43_48 + comp_r_43_49 + comp_r_43_50 + comp_r_43_51 + comp_r_43_52 + comp_r_43_53 + comp_r_43_54 + comp_r_43_55 + comp_r_43_56 + comp_r_43_57 + comp_r_43_58 + comp_r_43_59 + comp_r_43_60 + comp_r_43_61 + comp_r_43_62;
    position[44] = 22 + (sc_uint<1>) ~(comp_r_1_44) + (sc_uint<1>) ~(comp_r_3_44) + (sc_uint<1>) ~(comp_r_5_44) + (sc_uint<1>) ~(comp_r_7_44) + (sc_uint<1>) ~(comp_r_9_44) + (sc_uint<1>) ~(comp_r_11_44) + (sc_uint<1>) ~(comp_r_13_44) + (sc_uint<1>) ~(comp_r_15_44) + (sc_uint<1>) ~(comp_r_17_44) + (sc_uint<1>) ~(comp_r_19_44) + (sc_uint<1>) ~(comp_r_21_44) + (sc_uint<1>) ~(comp_r_23_44) + (sc_uint<1>) ~(comp_r_25_44) + (sc_uint<1>) ~(comp_r_27_44) + (sc_uint<1>) ~(comp_r_29_44) + (sc_uint<1>) ~(comp_r_31_44) + (sc_uint<1>) ~(comp_r_33_44) + (sc_uint<1>) ~(comp_r_35_44) + (sc_uint<1>) ~(comp_r_37_44) + (sc_uint<1>) ~(comp_r_39_44) + (sc_uint<1>) ~(comp_r_41_44) + (sc_uint<1>) ~(comp_r_43_44);
    position[45] = 23 + (sc_uint<1>) ~(comp_r_1_45) + (sc_uint<1>) ~(comp_r_3_45) + (sc_uint<1>) ~(comp_r_5_45) + (sc_uint<1>) ~(comp_r_7_45) + (sc_uint<1>) ~(comp_r_9_45) + (sc_uint<1>) ~(comp_r_11_45) + (sc_uint<1>) ~(comp_r_13_45) + (sc_uint<1>) ~(comp_r_15_45) + (sc_uint<1>) ~(comp_r_17_45) + (sc_uint<1>) ~(comp_r_19_45) + (sc_uint<1>) ~(comp_r_21_45) + (sc_uint<1>) ~(comp_r_23_45) + (sc_uint<1>) ~(comp_r_25_45) + (sc_uint<1>) ~(comp_r_27_45) + (sc_uint<1>) ~(comp_r_29_45) + (sc_uint<1>) ~(comp_r_31_45) + (sc_uint<1>) ~(comp_r_33_45) + (sc_uint<1>) ~(comp_r_35_45) + (sc_uint<1>) ~(comp_r_37_45) + (sc_uint<1>) ~(comp_r_39_45) + (sc_uint<1>) ~(comp_r_41_45) + (sc_uint<1>) ~(comp_r_43_45) + comp_r_45_46 + comp_r_45_47 + comp_r_45_48 + comp_r_45_49 + comp_r_45_50 + comp_r_45_51 + comp_r_45_52 + comp_r_45_53 + comp_r_45_54 + comp_r_45_55 + comp_r_45_56 + comp_r_45_57 + comp_r_45_58 + comp_r_45_59 + comp_r_45_60 + comp_r_45_61 + comp_r_45_62;
    position[46] = 23 + (sc_uint<1>) ~(comp_r_1_46) + (sc_uint<1>) ~(comp_r_3_46) + (sc_uint<1>) ~(comp_r_5_46) + (sc_uint<1>) ~(comp_r_7_46) + (sc_uint<1>) ~(comp_r_9_46) + (sc_uint<1>) ~(comp_r_11_46) + (sc_uint<1>) ~(comp_r_13_46) + (sc_uint<1>) ~(comp_r_15_46) + (sc_uint<1>) ~(comp_r_17_46) + (sc_uint<1>) ~(comp_r_19_46) + (sc_uint<1>) ~(comp_r_21_46) + (sc_uint<1>) ~(comp_r_23_46) + (sc_uint<1>) ~(comp_r_25_46) + (sc_uint<1>) ~(comp_r_27_46) + (sc_uint<1>) ~(comp_r_29_46) + (sc_uint<1>) ~(comp_r_31_46) + (sc_uint<1>) ~(comp_r_33_46) + (sc_uint<1>) ~(comp_r_35_46) + (sc_uint<1>) ~(comp_r_37_46) + (sc_uint<1>) ~(comp_r_39_46) + (sc_uint<1>) ~(comp_r_41_46) + (sc_uint<1>) ~(comp_r_43_46) + (sc_uint<1>) ~(comp_r_45_46);
    position[47] = 24 + (sc_uint<1>) ~(comp_r_1_47) + (sc_uint<1>) ~(comp_r_3_47) + (sc_uint<1>) ~(comp_r_5_47) + (sc_uint<1>) ~(comp_r_7_47) + (sc_uint<1>) ~(comp_r_9_47) + (sc_uint<1>) ~(comp_r_11_47) + (sc_uint<1>) ~(comp_r_13_47) + (sc_uint<1>) ~(comp_r_15_47) + (sc_uint<1>) ~(comp_r_17_47) + (sc_uint<1>) ~(comp_r_19_47) + (sc_uint<1>) ~(comp_r_21_47) + (sc_uint<1>) ~(comp_r_23_47) + (sc_uint<1>) ~(comp_r_25_47) + (sc_uint<1>) ~(comp_r_27_47) + (sc_uint<1>) ~(comp_r_29_47) + (sc_uint<1>) ~(comp_r_31_47) + (sc_uint<1>) ~(comp_r_33_47) + (sc_uint<1>) ~(comp_r_35_47) + (sc_uint<1>) ~(comp_r_37_47) + (sc_uint<1>) ~(comp_r_39_47) + (sc_uint<1>) ~(comp_r_41_47) + (sc_uint<1>) ~(comp_r_43_47) + (sc_uint<1>) ~(comp_r_45_47) + comp_r_47_48 + comp_r_47_49 + comp_r_47_50 + comp_r_47_51 + comp_r_47_52 + comp_r_47_53 + comp_r_47_54 + comp_r_47_55 + comp_r_47_56 + comp_r_47_57 + comp_r_47_58 + comp_r_47_59 + comp_r_47_60 + comp_r_47_61 + comp_r_47_62;
    position[48] = 24 + (sc_uint<1>) ~(comp_r_1_48) + (sc_uint<1>) ~(comp_r_3_48) + (sc_uint<1>) ~(comp_r_5_48) + (sc_uint<1>) ~(comp_r_7_48) + (sc_uint<1>) ~(comp_r_9_48) + (sc_uint<1>) ~(comp_r_11_48) + (sc_uint<1>) ~(comp_r_13_48) + (sc_uint<1>) ~(comp_r_15_48) + (sc_uint<1>) ~(comp_r_17_48) + (sc_uint<1>) ~(comp_r_19_48) + (sc_uint<1>) ~(comp_r_21_48) + (sc_uint<1>) ~(comp_r_23_48) + (sc_uint<1>) ~(comp_r_25_48) + (sc_uint<1>) ~(comp_r_27_48) + (sc_uint<1>) ~(comp_r_29_48) + (sc_uint<1>) ~(comp_r_31_48) + (sc_uint<1>) ~(comp_r_33_48) + (sc_uint<1>) ~(comp_r_35_48) + (sc_uint<1>) ~(comp_r_37_48) + (sc_uint<1>) ~(comp_r_39_48) + (sc_uint<1>) ~(comp_r_41_48) + (sc_uint<1>) ~(comp_r_43_48) + (sc_uint<1>) ~(comp_r_45_48) + (sc_uint<1>) ~(comp_r_47_48);
    position[49] = 25 + (sc_uint<1>) ~(comp_r_1_49) + (sc_uint<1>) ~(comp_r_3_49) + (sc_uint<1>) ~(comp_r_5_49) + (sc_uint<1>) ~(comp_r_7_49) + (sc_uint<1>) ~(comp_r_9_49) + (sc_uint<1>) ~(comp_r_11_49) + (sc_uint<1>) ~(comp_r_13_49) + (sc_uint<1>) ~(comp_r_15_49) + (sc_uint<1>) ~(comp_r_17_49) + (sc_uint<1>) ~(comp_r_19_49) + (sc_uint<1>) ~(comp_r_21_49) + (sc_uint<1>) ~(comp_r_23_49) + (sc_uint<1>) ~(comp_r_25_49) + (sc_uint<1>) ~(comp_r_27_49) + (sc_uint<1>) ~(comp_r_29_49) + (sc_uint<1>) ~(comp_r_31_49) + (sc_uint<1>) ~(comp_r_33_49) + (sc_uint<1>) ~(comp_r_35_49) + (sc_uint<1>) ~(comp_r_37_49) + (sc_uint<1>) ~(comp_r_39_49) + (sc_uint<1>) ~(comp_r_41_49) + (sc_uint<1>) ~(comp_r_43_49) + (sc_uint<1>) ~(comp_r_45_49) + (sc_uint<1>) ~(comp_r_47_49) + comp_r_49_50 + comp_r_49_51 + comp_r_49_52 + comp_r_49_53 + comp_r_49_54 + comp_r_49_55 + comp_r_49_56 + comp_r_49_57 + comp_r_49_58 + comp_r_49_59 + comp_r_49_60 + comp_r_49_61 + comp_r_49_62;
    position[50] = 25 + (sc_uint<1>) ~(comp_r_1_50) + (sc_uint<1>) ~(comp_r_3_50) + (sc_uint<1>) ~(comp_r_5_50) + (sc_uint<1>) ~(comp_r_7_50) + (sc_uint<1>) ~(comp_r_9_50) + (sc_uint<1>) ~(comp_r_11_50) + (sc_uint<1>) ~(comp_r_13_50) + (sc_uint<1>) ~(comp_r_15_50) + (sc_uint<1>) ~(comp_r_17_50) + (sc_uint<1>) ~(comp_r_19_50) + (sc_uint<1>) ~(comp_r_21_50) + (sc_uint<1>) ~(comp_r_23_50) + (sc_uint<1>) ~(comp_r_25_50) + (sc_uint<1>) ~(comp_r_27_50) + (sc_uint<1>) ~(comp_r_29_50) + (sc_uint<1>) ~(comp_r_31_50) + (sc_uint<1>) ~(comp_r_33_50) + (sc_uint<1>) ~(comp_r_35_50) + (sc_uint<1>) ~(comp_r_37_50) + (sc_uint<1>) ~(comp_r_39_50) + (sc_uint<1>) ~(comp_r_41_50) + (sc_uint<1>) ~(comp_r_43_50) + (sc_uint<1>) ~(comp_r_45_50) + (sc_uint<1>) ~(comp_r_47_50) + (sc_uint<1>) ~(comp_r_49_50);
    position[51] = 26 + (sc_uint<1>) ~(comp_r_1_51) + (sc_uint<1>) ~(comp_r_3_51) + (sc_uint<1>) ~(comp_r_5_51) + (sc_uint<1>) ~(comp_r_7_51) + (sc_uint<1>) ~(comp_r_9_51) + (sc_uint<1>) ~(comp_r_11_51) + (sc_uint<1>) ~(comp_r_13_51) + (sc_uint<1>) ~(comp_r_15_51) + (sc_uint<1>) ~(comp_r_17_51) + (sc_uint<1>) ~(comp_r_19_51) + (sc_uint<1>) ~(comp_r_21_51) + (sc_uint<1>) ~(comp_r_23_51) + (sc_uint<1>) ~(comp_r_25_51) + (sc_uint<1>) ~(comp_r_27_51) + (sc_uint<1>) ~(comp_r_29_51) + (sc_uint<1>) ~(comp_r_31_51) + (sc_uint<1>) ~(comp_r_33_51) + (sc_uint<1>) ~(comp_r_35_51) + (sc_uint<1>) ~(comp_r_37_51) + (sc_uint<1>) ~(comp_r_39_51) + (sc_uint<1>) ~(comp_r_41_51) + (sc_uint<1>) ~(comp_r_43_51) + (sc_uint<1>) ~(comp_r_45_51) + (sc_uint<1>) ~(comp_r_47_51) + (sc_uint<1>) ~(comp_r_49_51) + comp_r_51_52 + comp_r_51_53 + comp_r_51_54 + comp_r_51_55 + comp_r_51_56 + comp_r_51_57 + comp_r_51_58 + comp_r_51_59 + comp_r_51_60 + comp_r_51_61 + comp_r_51_62;
    position[52] = 26 + (sc_uint<1>) ~(comp_r_1_52) + (sc_uint<1>) ~(comp_r_3_52) + (sc_uint<1>) ~(comp_r_5_52) + (sc_uint<1>) ~(comp_r_7_52) + (sc_uint<1>) ~(comp_r_9_52) + (sc_uint<1>) ~(comp_r_11_52) + (sc_uint<1>) ~(comp_r_13_52) + (sc_uint<1>) ~(comp_r_15_52) + (sc_uint<1>) ~(comp_r_17_52) + (sc_uint<1>) ~(comp_r_19_52) + (sc_uint<1>) ~(comp_r_21_52) + (sc_uint<1>) ~(comp_r_23_52) + (sc_uint<1>) ~(comp_r_25_52) + (sc_uint<1>) ~(comp_r_27_52) + (sc_uint<1>) ~(comp_r_29_52) + (sc_uint<1>) ~(comp_r_31_52) + (sc_uint<1>) ~(comp_r_33_52) + (sc_uint<1>) ~(comp_r_35_52) + (sc_uint<1>) ~(comp_r_37_52) + (sc_uint<1>) ~(comp_r_39_52) + (sc_uint<1>) ~(comp_r_41_52) + (sc_uint<1>) ~(comp_r_43_52) + (sc_uint<1>) ~(comp_r_45_52) + (sc_uint<1>) ~(comp_r_47_52) + (sc_uint<1>) ~(comp_r_49_52) + (sc_uint<1>) ~(comp_r_51_52);
    position[53] = 27 + (sc_uint<1>) ~(comp_r_1_53) + (sc_uint<1>) ~(comp_r_3_53) + (sc_uint<1>) ~(comp_r_5_53) + (sc_uint<1>) ~(comp_r_7_53) + (sc_uint<1>) ~(comp_r_9_53) + (sc_uint<1>) ~(comp_r_11_53) + (sc_uint<1>) ~(comp_r_13_53) + (sc_uint<1>) ~(comp_r_15_53) + (sc_uint<1>) ~(comp_r_17_53) + (sc_uint<1>) ~(comp_r_19_53) + (sc_uint<1>) ~(comp_r_21_53) + (sc_uint<1>) ~(comp_r_23_53) + (sc_uint<1>) ~(comp_r_25_53) + (sc_uint<1>) ~(comp_r_27_53) + (sc_uint<1>) ~(comp_r_29_53) + (sc_uint<1>) ~(comp_r_31_53) + (sc_uint<1>) ~(comp_r_33_53) + (sc_uint<1>) ~(comp_r_35_53) + (sc_uint<1>) ~(comp_r_37_53) + (sc_uint<1>) ~(comp_r_39_53) + (sc_uint<1>) ~(comp_r_41_53) + (sc_uint<1>) ~(comp_r_43_53) + (sc_uint<1>) ~(comp_r_45_53) + (sc_uint<1>) ~(comp_r_47_53) + (sc_uint<1>) ~(comp_r_49_53) + (sc_uint<1>) ~(comp_r_51_53) + comp_r_53_54 + comp_r_53_55 + comp_r_53_56 + comp_r_53_57 + comp_r_53_58 + comp_r_53_59 + comp_r_53_60 + comp_r_53_61 + comp_r_53_62;
    position[54] = 27 + (sc_uint<1>) ~(comp_r_1_54) + (sc_uint<1>) ~(comp_r_3_54) + (sc_uint<1>) ~(comp_r_5_54) + (sc_uint<1>) ~(comp_r_7_54) + (sc_uint<1>) ~(comp_r_9_54) + (sc_uint<1>) ~(comp_r_11_54) + (sc_uint<1>) ~(comp_r_13_54) + (sc_uint<1>) ~(comp_r_15_54) + (sc_uint<1>) ~(comp_r_17_54) + (sc_uint<1>) ~(comp_r_19_54) + (sc_uint<1>) ~(comp_r_21_54) + (sc_uint<1>) ~(comp_r_23_54) + (sc_uint<1>) ~(comp_r_25_54) + (sc_uint<1>) ~(comp_r_27_54) + (sc_uint<1>) ~(comp_r_29_54) + (sc_uint<1>) ~(comp_r_31_54) + (sc_uint<1>) ~(comp_r_33_54) + (sc_uint<1>) ~(comp_r_35_54) + (sc_uint<1>) ~(comp_r_37_54) + (sc_uint<1>) ~(comp_r_39_54) + (sc_uint<1>) ~(comp_r_41_54) + (sc_uint<1>) ~(comp_r_43_54) + (sc_uint<1>) ~(comp_r_45_54) + (sc_uint<1>) ~(comp_r_47_54) + (sc_uint<1>) ~(comp_r_49_54) + (sc_uint<1>) ~(comp_r_51_54) + (sc_uint<1>) ~(comp_r_53_54);
    position[55] = 28 + (sc_uint<1>) ~(comp_r_1_55) + (sc_uint<1>) ~(comp_r_3_55) + (sc_uint<1>) ~(comp_r_5_55) + (sc_uint<1>) ~(comp_r_7_55) + (sc_uint<1>) ~(comp_r_9_55) + (sc_uint<1>) ~(comp_r_11_55) + (sc_uint<1>) ~(comp_r_13_55) + (sc_uint<1>) ~(comp_r_15_55) + (sc_uint<1>) ~(comp_r_17_55) + (sc_uint<1>) ~(comp_r_19_55) + (sc_uint<1>) ~(comp_r_21_55) + (sc_uint<1>) ~(comp_r_23_55) + (sc_uint<1>) ~(comp_r_25_55) + (sc_uint<1>) ~(comp_r_27_55) + (sc_uint<1>) ~(comp_r_29_55) + (sc_uint<1>) ~(comp_r_31_55) + (sc_uint<1>) ~(comp_r_33_55) + (sc_uint<1>) ~(comp_r_35_55) + (sc_uint<1>) ~(comp_r_37_55) + (sc_uint<1>) ~(comp_r_39_55) + (sc_uint<1>) ~(comp_r_41_55) + (sc_uint<1>) ~(comp_r_43_55) + (sc_uint<1>) ~(comp_r_45_55) + (sc_uint<1>) ~(comp_r_47_55) + (sc_uint<1>) ~(comp_r_49_55) + (sc_uint<1>) ~(comp_r_51_55) + (sc_uint<1>) ~(comp_r_53_55) + comp_r_55_56 + comp_r_55_57 + comp_r_55_58 + comp_r_55_59 + comp_r_55_60 + comp_r_55_61 + comp_r_55_62;
    position[56] = 28 + (sc_uint<1>) ~(comp_r_1_56) + (sc_uint<1>) ~(comp_r_3_56) + (sc_uint<1>) ~(comp_r_5_56) + (sc_uint<1>) ~(comp_r_7_56) + (sc_uint<1>) ~(comp_r_9_56) + (sc_uint<1>) ~(comp_r_11_56) + (sc_uint<1>) ~(comp_r_13_56) + (sc_uint<1>) ~(comp_r_15_56) + (sc_uint<1>) ~(comp_r_17_56) + (sc_uint<1>) ~(comp_r_19_56) + (sc_uint<1>) ~(comp_r_21_56) + (sc_uint<1>) ~(comp_r_23_56) + (sc_uint<1>) ~(comp_r_25_56) + (sc_uint<1>) ~(comp_r_27_56) + (sc_uint<1>) ~(comp_r_29_56) + (sc_uint<1>) ~(comp_r_31_56) + (sc_uint<1>) ~(comp_r_33_56) + (sc_uint<1>) ~(comp_r_35_56) + (sc_uint<1>) ~(comp_r_37_56) + (sc_uint<1>) ~(comp_r_39_56) + (sc_uint<1>) ~(comp_r_41_56) + (sc_uint<1>) ~(comp_r_43_56) + (sc_uint<1>) ~(comp_r_45_56) + (sc_uint<1>) ~(comp_r_47_56) + (sc_uint<1>) ~(comp_r_49_56) + (sc_uint<1>) ~(comp_r_51_56) + (sc_uint<1>) ~(comp_r_53_56) + (sc_uint<1>) ~(comp_r_55_56);
    position[57] = 29 + (sc_uint<1>) ~(comp_r_1_57) + (sc_uint<1>) ~(comp_r_3_57) + (sc_uint<1>) ~(comp_r_5_57) + (sc_uint<1>) ~(comp_r_7_57) + (sc_uint<1>) ~(comp_r_9_57) + (sc_uint<1>) ~(comp_r_11_57) + (sc_uint<1>) ~(comp_r_13_57) + (sc_uint<1>) ~(comp_r_15_57) + (sc_uint<1>) ~(comp_r_17_57) + (sc_uint<1>) ~(comp_r_19_57) + (sc_uint<1>) ~(comp_r_21_57) + (sc_uint<1>) ~(comp_r_23_57) + (sc_uint<1>) ~(comp_r_25_57) + (sc_uint<1>) ~(comp_r_27_57) + (sc_uint<1>) ~(comp_r_29_57) + (sc_uint<1>) ~(comp_r_31_57) + (sc_uint<1>) ~(comp_r_33_57) + (sc_uint<1>) ~(comp_r_35_57) + (sc_uint<1>) ~(comp_r_37_57) + (sc_uint<1>) ~(comp_r_39_57) + (sc_uint<1>) ~(comp_r_41_57) + (sc_uint<1>) ~(comp_r_43_57) + (sc_uint<1>) ~(comp_r_45_57) + (sc_uint<1>) ~(comp_r_47_57) + (sc_uint<1>) ~(comp_r_49_57) + (sc_uint<1>) ~(comp_r_51_57) + (sc_uint<1>) ~(comp_r_53_57) + (sc_uint<1>) ~(comp_r_55_57) + comp_r_57_58 + comp_r_57_59 + comp_r_57_60 + comp_r_57_61 + comp_r_57_62;
    position[58] = 29 + (sc_uint<1>) ~(comp_r_1_58) + (sc_uint<1>) ~(comp_r_3_58) + (sc_uint<1>) ~(comp_r_5_58) + (sc_uint<1>) ~(comp_r_7_58) + (sc_uint<1>) ~(comp_r_9_58) + (sc_uint<1>) ~(comp_r_11_58) + (sc_uint<1>) ~(comp_r_13_58) + (sc_uint<1>) ~(comp_r_15_58) + (sc_uint<1>) ~(comp_r_17_58) + (sc_uint<1>) ~(comp_r_19_58) + (sc_uint<1>) ~(comp_r_21_58) + (sc_uint<1>) ~(comp_r_23_58) + (sc_uint<1>) ~(comp_r_25_58) + (sc_uint<1>) ~(comp_r_27_58) + (sc_uint<1>) ~(comp_r_29_58) + (sc_uint<1>) ~(comp_r_31_58) + (sc_uint<1>) ~(comp_r_33_58) + (sc_uint<1>) ~(comp_r_35_58) + (sc_uint<1>) ~(comp_r_37_58) + (sc_uint<1>) ~(comp_r_39_58) + (sc_uint<1>) ~(comp_r_41_58) + (sc_uint<1>) ~(comp_r_43_58) + (sc_uint<1>) ~(comp_r_45_58) + (sc_uint<1>) ~(comp_r_47_58) + (sc_uint<1>) ~(comp_r_49_58) + (sc_uint<1>) ~(comp_r_51_58) + (sc_uint<1>) ~(comp_r_53_58) + (sc_uint<1>) ~(comp_r_55_58) + (sc_uint<1>) ~(comp_r_57_58);
    position[59] = 30 + (sc_uint<1>) ~(comp_r_1_59) + (sc_uint<1>) ~(comp_r_3_59) + (sc_uint<1>) ~(comp_r_5_59) + (sc_uint<1>) ~(comp_r_7_59) + (sc_uint<1>) ~(comp_r_9_59) + (sc_uint<1>) ~(comp_r_11_59) + (sc_uint<1>) ~(comp_r_13_59) + (sc_uint<1>) ~(comp_r_15_59) + (sc_uint<1>) ~(comp_r_17_59) + (sc_uint<1>) ~(comp_r_19_59) + (sc_uint<1>) ~(comp_r_21_59) + (sc_uint<1>) ~(comp_r_23_59) + (sc_uint<1>) ~(comp_r_25_59) + (sc_uint<1>) ~(comp_r_27_59) + (sc_uint<1>) ~(comp_r_29_59) + (sc_uint<1>) ~(comp_r_31_59) + (sc_uint<1>) ~(comp_r_33_59) + (sc_uint<1>) ~(comp_r_35_59) + (sc_uint<1>) ~(comp_r_37_59) + (sc_uint<1>) ~(comp_r_39_59) + (sc_uint<1>) ~(comp_r_41_59) + (sc_uint<1>) ~(comp_r_43_59) + (sc_uint<1>) ~(comp_r_45_59) + (sc_uint<1>) ~(comp_r_47_59) + (sc_uint<1>) ~(comp_r_49_59) + (sc_uint<1>) ~(comp_r_51_59) + (sc_uint<1>) ~(comp_r_53_59) + (sc_uint<1>) ~(comp_r_55_59) + (sc_uint<1>) ~(comp_r_57_59) + comp_r_59_60 + comp_r_59_61 + comp_r_59_62;
    position[60] = 30 + (sc_uint<1>) ~(comp_r_1_60) + (sc_uint<1>) ~(comp_r_3_60) + (sc_uint<1>) ~(comp_r_5_60) + (sc_uint<1>) ~(comp_r_7_60) + (sc_uint<1>) ~(comp_r_9_60) + (sc_uint<1>) ~(comp_r_11_60) + (sc_uint<1>) ~(comp_r_13_60) + (sc_uint<1>) ~(comp_r_15_60) + (sc_uint<1>) ~(comp_r_17_60) + (sc_uint<1>) ~(comp_r_19_60) + (sc_uint<1>) ~(comp_r_21_60) + (sc_uint<1>) ~(comp_r_23_60) + (sc_uint<1>) ~(comp_r_25_60) + (sc_uint<1>) ~(comp_r_27_60) + (sc_uint<1>) ~(comp_r_29_60) + (sc_uint<1>) ~(comp_r_31_60) + (sc_uint<1>) ~(comp_r_33_60) + (sc_uint<1>) ~(comp_r_35_60) + (sc_uint<1>) ~(comp_r_37_60) + (sc_uint<1>) ~(comp_r_39_60) + (sc_uint<1>) ~(comp_r_41_60) + (sc_uint<1>) ~(comp_r_43_60) + (sc_uint<1>) ~(comp_r_45_60) + (sc_uint<1>) ~(comp_r_47_60) + (sc_uint<1>) ~(comp_r_49_60) + (sc_uint<1>) ~(comp_r_51_60) + (sc_uint<1>) ~(comp_r_53_60) + (sc_uint<1>) ~(comp_r_55_60) + (sc_uint<1>) ~(comp_r_57_60) + (sc_uint<1>) ~(comp_r_59_60);
    position[61] = 31 + (sc_uint<1>) ~(comp_r_1_61) + (sc_uint<1>) ~(comp_r_3_61) + (sc_uint<1>) ~(comp_r_5_61) + (sc_uint<1>) ~(comp_r_7_61) + (sc_uint<1>) ~(comp_r_9_61) + (sc_uint<1>) ~(comp_r_11_61) + (sc_uint<1>) ~(comp_r_13_61) + (sc_uint<1>) ~(comp_r_15_61) + (sc_uint<1>) ~(comp_r_17_61) + (sc_uint<1>) ~(comp_r_19_61) + (sc_uint<1>) ~(comp_r_21_61) + (sc_uint<1>) ~(comp_r_23_61) + (sc_uint<1>) ~(comp_r_25_61) + (sc_uint<1>) ~(comp_r_27_61) + (sc_uint<1>) ~(comp_r_29_61) + (sc_uint<1>) ~(comp_r_31_61) + (sc_uint<1>) ~(comp_r_33_61) + (sc_uint<1>) ~(comp_r_35_61) + (sc_uint<1>) ~(comp_r_37_61) + (sc_uint<1>) ~(comp_r_39_61) + (sc_uint<1>) ~(comp_r_41_61) + (sc_uint<1>) ~(comp_r_43_61) + (sc_uint<1>) ~(comp_r_45_61) + (sc_uint<1>) ~(comp_r_47_61) + (sc_uint<1>) ~(comp_r_49_61) + (sc_uint<1>) ~(comp_r_51_61) + (sc_uint<1>) ~(comp_r_53_61) + (sc_uint<1>) ~(comp_r_55_61) + (sc_uint<1>) ~(comp_r_57_61) + (sc_uint<1>) ~(comp_r_59_61) + comp_r_61_62;
    position[62] = 31 + (sc_uint<1>) ~(comp_r_1_62) + (sc_uint<1>) ~(comp_r_3_62) + (sc_uint<1>) ~(comp_r_5_62) + (sc_uint<1>) ~(comp_r_7_62) + (sc_uint<1>) ~(comp_r_9_62) + (sc_uint<1>) ~(comp_r_11_62) + (sc_uint<1>) ~(comp_r_13_62) + (sc_uint<1>) ~(comp_r_15_62) + (sc_uint<1>) ~(comp_r_17_62) + (sc_uint<1>) ~(comp_r_19_62) + (sc_uint<1>) ~(comp_r_21_62) + (sc_uint<1>) ~(comp_r_23_62) + (sc_uint<1>) ~(comp_r_25_62) + (sc_uint<1>) ~(comp_r_27_62) + (sc_uint<1>) ~(comp_r_29_62) + (sc_uint<1>) ~(comp_r_31_62) + (sc_uint<1>) ~(comp_r_33_62) + (sc_uint<1>) ~(comp_r_35_62) + (sc_uint<1>) ~(comp_r_37_62) + (sc_uint<1>) ~(comp_r_39_62) + (sc_uint<1>) ~(comp_r_41_62) + (sc_uint<1>) ~(comp_r_43_62) + (sc_uint<1>) ~(comp_r_45_62) + (sc_uint<1>) ~(comp_r_47_62) + (sc_uint<1>) ~(comp_r_49_62) + (sc_uint<1>) ~(comp_r_51_62) + (sc_uint<1>) ~(comp_r_53_62) + (sc_uint<1>) ~(comp_r_55_62) + (sc_uint<1>) ~(comp_r_57_62) + (sc_uint<1>) ~(comp_r_59_62) + (sc_uint<1>) ~(comp_r_61_62);
    position[63] = (sc_uint<6>) 63;

 // Multiplex the inputs in the order
    PS_struct<1,Q,5> temp[63];
#pragma HLS ARRAY_PARTITION variable=temp complete dim=1
    temp[0] = input[0];
    temp[1] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 1);
    temp[2] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 2);
    temp[3] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 3);
    temp[4] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 4);
    temp[5] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 5);
    temp[6] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 6);
    temp[7] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 7);
    temp[8] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 8);
    temp[9] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 9);
    temp[10] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 10);
    temp[11] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 11);
    temp[12] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 12);
    temp[13] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 13);
    temp[14] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 14);
    temp[15] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 15);
    temp[16] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 16);
    temp[17] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 17);
    temp[18] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 18);
    temp[19] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 19);
    temp[20] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 20);
    temp[21] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 21);
    temp[22] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 22);
    temp[23] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 23);
    temp[24] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 24);
    temp[25] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 25);
    temp[26] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 26);
    temp[27] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 27);
    temp[28] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 28);
    temp[29] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 29);
    temp[30] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 30);
    temp[31] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 31);
    temp[32] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 32);
    temp[33] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 33);
    temp[34] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 34);
    temp[35] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 35);
    temp[36] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 36);
    temp[37] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 37);
    temp[38] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 38);
    temp[39] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 39);
    temp[40] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 40);
    temp[41] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 41);
    temp[42] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 42);
    temp[43] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 43);
    temp[44] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 44);
    temp[45] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 45);
    temp[46] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 46);
    temp[47] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 47);
    temp[48] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 48);
    temp[49] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 49);
    temp[50] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 50);
    temp[51] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 51);
    temp[52] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 52);
    temp[53] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 53);
    temp[54] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 54);
    temp[55] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 55);
    temp[56] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 56);
    temp[57] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 57);
    temp[58] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 58);
    temp[59] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 59);
    temp[60] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 60);
    temp[61] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 61);
    temp[62] = RO_MUX64 < PS_struct<1,Q,5> > (input, position, (sc_uint<6>) 62);

// output the L smallest input
    if (fb == 0){
    output[0] = temp[31]; 
    output[1] = temp[32]; 
    output[2] = temp[33]; 
    output[3] = temp[34]; 
    output[4] = temp[35]; 
    output[5] = temp[36]; 
    output[6] = temp[37]; 
    output[7] = temp[38]; 
    output[8] = temp[39]; 
    output[9] = temp[40]; 
    output[10] = temp[41]; 
    output[11] = temp[42]; 
    output[12] = temp[43]; 
    output[13] = temp[44]; 
    output[14] = temp[45]; 
    output[15] = temp[46]; 
    output[16] = temp[47]; 
    output[17] = temp[48]; 
    output[18] = temp[49]; 
    output[19] = temp[50]; 
    output[20] = temp[51]; 
    output[21] = temp[52]; 
    output[22] = temp[53]; 
    output[23] = temp[54]; 
    output[24] = temp[55]; 
    output[25] = temp[56]; 
    output[26] = temp[57]; 
    output[27] = temp[58]; 
    output[28] = temp[59]; 
    output[29] = temp[60]; 
    output[30] = temp[61]; 
    output[31] = temp[62]; 
    }
    else{
    output[0] = temp[0]; 
    output[1] = temp[1]; 
    output[2] = temp[2]; 
    output[3] = temp[3]; 
    output[4] = temp[4]; 
    output[5] = temp[5]; 
    output[6] = temp[6]; 
    output[7] = temp[7]; 
    output[8] = temp[8]; 
    output[9] = temp[9]; 
    output[10] = temp[10]; 
    output[11] = temp[11]; 
    output[12] = temp[12]; 
    output[13] = temp[13]; 
    output[14] = temp[14]; 
    output[15] = temp[15]; 
    output[16] = temp[16]; 
    output[17] = temp[17]; 
    output[18] = temp[18]; 
    output[19] = temp[19]; 
    output[20] = temp[20]; 
    output[21] = temp[21]; 
    output[22] = temp[22]; 
    output[23] = temp[23]; 
    output[24] = temp[24]; 
    output[25] = temp[25]; 
    output[26] = temp[26]; 
    output[27] = temp[27]; 
    output[28] = temp[28]; 
    output[29] = temp[29]; 
    output[30] = temp[30]; 
    output[31] = temp[31]; 
    }
}

template <int Q> 
void RANKORDER_SORT_L64 (PS_struct<1,Q,6> input[128], PS_struct<1,Q,6> output[64], sc_biguint<1> fb)
{
#pragma HLS INLINE off
#pragma HLS ARRAY_PARTITION variable=input complete dim=1
#pragma HLS ARRAY_PARTITION variable=output complete dim=1

// Adjust input when FB=0

    if (fb == 0){
        input[127].metric = MAX_VAL;
        input[126] = input[63];
        input[125] = input[62];
        input[124].metric = 0;
        input[123] = input[61];
        input[122].metric = 0;
        input[121] = input[60];
        input[120].metric = 0;
        input[119] = input[59];
        input[118].metric = 0;
        input[117] = input[58];
        input[116].metric = 0;
        input[115] = input[57];
        input[114].metric = 0;
        input[113] = input[56];
        input[112].metric = 0;
        input[111] = input[55];
        input[110].metric = 0;
        input[109] = input[54];
        input[108].metric = 0;
        input[107] = input[53];
        input[106].metric = 0;
        input[105] = input[52];
        input[104].metric = 0;
        input[103] = input[51];
        input[102].metric = 0;
        input[101] = input[50];
        input[100].metric = 0;
        input[99] = input[49];
        input[98].metric = 0;
        input[97] = input[48];
        input[96].metric = 0;
        input[95] = input[47];
        input[94].metric = 0;
        input[93] = input[46];
        input[92].metric = 0;
        input[91] = input[45];
        input[90].metric = 0;
        input[89] = input[44];
        input[88].metric = 0;
        input[87] = input[43];
        input[86].metric = 0;
        input[85] = input[42];
        input[84].metric = 0;
        input[83] = input[41];
        input[82].metric = 0;
        input[81] = input[40];
        input[80].metric = 0;
        input[79] = input[39];
        input[78].metric = 0;
        input[77] = input[38];
        input[76].metric = 0;
        input[75] = input[37];
        input[74].metric = 0;
        input[73] = input[36];
        input[72].metric = 0;
        input[71] = input[35];
        input[70].metric = 0;
        input[69] = input[34];
        input[68].metric = 0;
        input[67] = input[33];
        input[66].metric = 0;
        input[65] = input[32];
        input[64].metric = 0;
        input[63] = input[31];
        input[62].metric = 0;
        input[61] = input[30];
        input[60].metric = 0;
        input[59] = input[29];
        input[58].metric = 0;
        input[57] = input[28];
        input[56].metric = 0;
        input[55] = input[27];
        input[54].metric = 0;
        input[53] = input[26];
        input[52].metric = 0;
        input[51] = input[25];
        input[50].metric = 0;
        input[49] = input[24];
        input[48].metric = 0;
        input[47] = input[23];
        input[46].metric = 0;
        input[45] = input[22];
        input[44].metric = 0;
        input[43] = input[21];
        input[42].metric = 0;
        input[41] = input[20];
        input[40].metric = 0;
        input[39] = input[19];
        input[38].metric = 0;
        input[37] = input[18];
        input[36].metric = 0;
        input[35] = input[17];
        input[34].metric = 0;
        input[33] = input[16];
        input[32].metric = 0;
        input[31] = input[15];
        input[30].metric = 0;
        input[29] = input[14];
        input[28].metric = 0;
        input[27] = input[13];
        input[26].metric = 0;
        input[25] = input[12];
        input[24].metric = 0;
        input[23] = input[11];
        input[22].metric = 0;
        input[21] = input[10];
        input[20].metric = 0;
        input[19] = input[9];
        input[18].metric = 0;
        input[17] = input[8];
        input[16].metric = 0;
        input[15] = input[7];
        input[14].metric = 0;
        input[13] = input[6];
        input[12].metric = 0;
        input[11] = input[5];
        input[10].metric = 0;
        input[9] = input[4];
        input[8].metric = 0;
        input[7] = input[3];
        input[6].metric = 0;
        input[5] = input[2];
        input[4].metric = 0;
        input[3] = input[1];
        input[2].metric = 0;
        input[1] = input[0];
    }

    sc_uint<7> position[128];
#pragma HLS ARRAY_PARTITION variable=position complete dim=1

// COMPARE inputs

    sc_uint<1> comp_r_1_2 = COMP < Q > (input[1].metric, input[2].metric);
    sc_uint<1> comp_r_1_3 = COMP < Q > (input[1].metric, input[3].metric);
    sc_uint<1> comp_r_1_4 = COMP < Q > (input[1].metric, input[4].metric);
    sc_uint<1> comp_r_1_5 = COMP < Q > (input[1].metric, input[5].metric);
    sc_uint<1> comp_r_1_6 = COMP < Q > (input[1].metric, input[6].metric);
    sc_uint<1> comp_r_1_7 = COMP < Q > (input[1].metric, input[7].metric);
    sc_uint<1> comp_r_1_8 = COMP < Q > (input[1].metric, input[8].metric);
    sc_uint<1> comp_r_1_9 = COMP < Q > (input[1].metric, input[9].metric);
    sc_uint<1> comp_r_1_10 = COMP < Q > (input[1].metric, input[10].metric);
    sc_uint<1> comp_r_1_11 = COMP < Q > (input[1].metric, input[11].metric);
    sc_uint<1> comp_r_1_12 = COMP < Q > (input[1].metric, input[12].metric);
    sc_uint<1> comp_r_1_13 = COMP < Q > (input[1].metric, input[13].metric);
    sc_uint<1> comp_r_1_14 = COMP < Q > (input[1].metric, input[14].metric);
    sc_uint<1> comp_r_1_15 = COMP < Q > (input[1].metric, input[15].metric);
    sc_uint<1> comp_r_1_16 = COMP < Q > (input[1].metric, input[16].metric);
    sc_uint<1> comp_r_1_17 = COMP < Q > (input[1].metric, input[17].metric);
    sc_uint<1> comp_r_1_18 = COMP < Q > (input[1].metric, input[18].metric);
    sc_uint<1> comp_r_1_19 = COMP < Q > (input[1].metric, input[19].metric);
    sc_uint<1> comp_r_1_20 = COMP < Q > (input[1].metric, input[20].metric);
    sc_uint<1> comp_r_1_21 = COMP < Q > (input[1].metric, input[21].metric);
    sc_uint<1> comp_r_1_22 = COMP < Q > (input[1].metric, input[22].metric);
    sc_uint<1> comp_r_1_23 = COMP < Q > (input[1].metric, input[23].metric);
    sc_uint<1> comp_r_1_24 = COMP < Q > (input[1].metric, input[24].metric);
    sc_uint<1> comp_r_1_25 = COMP < Q > (input[1].metric, input[25].metric);
    sc_uint<1> comp_r_1_26 = COMP < Q > (input[1].metric, input[26].metric);
    sc_uint<1> comp_r_1_27 = COMP < Q > (input[1].metric, input[27].metric);
    sc_uint<1> comp_r_1_28 = COMP < Q > (input[1].metric, input[28].metric);
    sc_uint<1> comp_r_1_29 = COMP < Q > (input[1].metric, input[29].metric);
    sc_uint<1> comp_r_1_30 = COMP < Q > (input[1].metric, input[30].metric);
    sc_uint<1> comp_r_1_31 = COMP < Q > (input[1].metric, input[31].metric);
    sc_uint<1> comp_r_1_32 = COMP < Q > (input[1].metric, input[32].metric);
    sc_uint<1> comp_r_1_33 = COMP < Q > (input[1].metric, input[33].metric);
    sc_uint<1> comp_r_1_34 = COMP < Q > (input[1].metric, input[34].metric);
    sc_uint<1> comp_r_1_35 = COMP < Q > (input[1].metric, input[35].metric);
    sc_uint<1> comp_r_1_36 = COMP < Q > (input[1].metric, input[36].metric);
    sc_uint<1> comp_r_1_37 = COMP < Q > (input[1].metric, input[37].metric);
    sc_uint<1> comp_r_1_38 = COMP < Q > (input[1].metric, input[38].metric);
    sc_uint<1> comp_r_1_39 = COMP < Q > (input[1].metric, input[39].metric);
    sc_uint<1> comp_r_1_40 = COMP < Q > (input[1].metric, input[40].metric);
    sc_uint<1> comp_r_1_41 = COMP < Q > (input[1].metric, input[41].metric);
    sc_uint<1> comp_r_1_42 = COMP < Q > (input[1].metric, input[42].metric);
    sc_uint<1> comp_r_1_43 = COMP < Q > (input[1].metric, input[43].metric);
    sc_uint<1> comp_r_1_44 = COMP < Q > (input[1].metric, input[44].metric);
    sc_uint<1> comp_r_1_45 = COMP < Q > (input[1].metric, input[45].metric);
    sc_uint<1> comp_r_1_46 = COMP < Q > (input[1].metric, input[46].metric);
    sc_uint<1> comp_r_1_47 = COMP < Q > (input[1].metric, input[47].metric);
    sc_uint<1> comp_r_1_48 = COMP < Q > (input[1].metric, input[48].metric);
    sc_uint<1> comp_r_1_49 = COMP < Q > (input[1].metric, input[49].metric);
    sc_uint<1> comp_r_1_50 = COMP < Q > (input[1].metric, input[50].metric);
    sc_uint<1> comp_r_1_51 = COMP < Q > (input[1].metric, input[51].metric);
    sc_uint<1> comp_r_1_52 = COMP < Q > (input[1].metric, input[52].metric);
    sc_uint<1> comp_r_1_53 = COMP < Q > (input[1].metric, input[53].metric);
    sc_uint<1> comp_r_1_54 = COMP < Q > (input[1].metric, input[54].metric);
    sc_uint<1> comp_r_1_55 = COMP < Q > (input[1].metric, input[55].metric);
    sc_uint<1> comp_r_1_56 = COMP < Q > (input[1].metric, input[56].metric);
    sc_uint<1> comp_r_1_57 = COMP < Q > (input[1].metric, input[57].metric);
    sc_uint<1> comp_r_1_58 = COMP < Q > (input[1].metric, input[58].metric);
    sc_uint<1> comp_r_1_59 = COMP < Q > (input[1].metric, input[59].metric);
    sc_uint<1> comp_r_1_60 = COMP < Q > (input[1].metric, input[60].metric);
    sc_uint<1> comp_r_1_61 = COMP < Q > (input[1].metric, input[61].metric);
    sc_uint<1> comp_r_1_62 = COMP < Q > (input[1].metric, input[62].metric);
    sc_uint<1> comp_r_1_63 = COMP < Q > (input[1].metric, input[63].metric);
    sc_uint<1> comp_r_1_64 = COMP < Q > (input[1].metric, input[64].metric);
    sc_uint<1> comp_r_1_65 = COMP < Q > (input[1].metric, input[65].metric);
    sc_uint<1> comp_r_1_66 = COMP < Q > (input[1].metric, input[66].metric);
    sc_uint<1> comp_r_1_67 = COMP < Q > (input[1].metric, input[67].metric);
    sc_uint<1> comp_r_1_68 = COMP < Q > (input[1].metric, input[68].metric);
    sc_uint<1> comp_r_1_69 = COMP < Q > (input[1].metric, input[69].metric);
    sc_uint<1> comp_r_1_70 = COMP < Q > (input[1].metric, input[70].metric);
    sc_uint<1> comp_r_1_71 = COMP < Q > (input[1].metric, input[71].metric);
    sc_uint<1> comp_r_1_72 = COMP < Q > (input[1].metric, input[72].metric);
    sc_uint<1> comp_r_1_73 = COMP < Q > (input[1].metric, input[73].metric);
    sc_uint<1> comp_r_1_74 = COMP < Q > (input[1].metric, input[74].metric);
    sc_uint<1> comp_r_1_75 = COMP < Q > (input[1].metric, input[75].metric);
    sc_uint<1> comp_r_1_76 = COMP < Q > (input[1].metric, input[76].metric);
    sc_uint<1> comp_r_1_77 = COMP < Q > (input[1].metric, input[77].metric);
    sc_uint<1> comp_r_1_78 = COMP < Q > (input[1].metric, input[78].metric);
    sc_uint<1> comp_r_1_79 = COMP < Q > (input[1].metric, input[79].metric);
    sc_uint<1> comp_r_1_80 = COMP < Q > (input[1].metric, input[80].metric);
    sc_uint<1> comp_r_1_81 = COMP < Q > (input[1].metric, input[81].metric);
    sc_uint<1> comp_r_1_82 = COMP < Q > (input[1].metric, input[82].metric);
    sc_uint<1> comp_r_1_83 = COMP < Q > (input[1].metric, input[83].metric);
    sc_uint<1> comp_r_1_84 = COMP < Q > (input[1].metric, input[84].metric);
    sc_uint<1> comp_r_1_85 = COMP < Q > (input[1].metric, input[85].metric);
    sc_uint<1> comp_r_1_86 = COMP < Q > (input[1].metric, input[86].metric);
    sc_uint<1> comp_r_1_87 = COMP < Q > (input[1].metric, input[87].metric);
    sc_uint<1> comp_r_1_88 = COMP < Q > (input[1].metric, input[88].metric);
    sc_uint<1> comp_r_1_89 = COMP < Q > (input[1].metric, input[89].metric);
    sc_uint<1> comp_r_1_90 = COMP < Q > (input[1].metric, input[90].metric);
    sc_uint<1> comp_r_1_91 = COMP < Q > (input[1].metric, input[91].metric);
    sc_uint<1> comp_r_1_92 = COMP < Q > (input[1].metric, input[92].metric);
    sc_uint<1> comp_r_1_93 = COMP < Q > (input[1].metric, input[93].metric);
    sc_uint<1> comp_r_1_94 = COMP < Q > (input[1].metric, input[94].metric);
    sc_uint<1> comp_r_1_95 = COMP < Q > (input[1].metric, input[95].metric);
    sc_uint<1> comp_r_1_96 = COMP < Q > (input[1].metric, input[96].metric);
    sc_uint<1> comp_r_1_97 = COMP < Q > (input[1].metric, input[97].metric);
    sc_uint<1> comp_r_1_98 = COMP < Q > (input[1].metric, input[98].metric);
    sc_uint<1> comp_r_1_99 = COMP < Q > (input[1].metric, input[99].metric);
    sc_uint<1> comp_r_1_100 = COMP < Q > (input[1].metric, input[100].metric);
    sc_uint<1> comp_r_1_101 = COMP < Q > (input[1].metric, input[101].metric);
    sc_uint<1> comp_r_1_102 = COMP < Q > (input[1].metric, input[102].metric);
    sc_uint<1> comp_r_1_103 = COMP < Q > (input[1].metric, input[103].metric);
    sc_uint<1> comp_r_1_104 = COMP < Q > (input[1].metric, input[104].metric);
    sc_uint<1> comp_r_1_105 = COMP < Q > (input[1].metric, input[105].metric);
    sc_uint<1> comp_r_1_106 = COMP < Q > (input[1].metric, input[106].metric);
    sc_uint<1> comp_r_1_107 = COMP < Q > (input[1].metric, input[107].metric);
    sc_uint<1> comp_r_1_108 = COMP < Q > (input[1].metric, input[108].metric);
    sc_uint<1> comp_r_1_109 = COMP < Q > (input[1].metric, input[109].metric);
    sc_uint<1> comp_r_1_110 = COMP < Q > (input[1].metric, input[110].metric);
    sc_uint<1> comp_r_1_111 = COMP < Q > (input[1].metric, input[111].metric);
    sc_uint<1> comp_r_1_112 = COMP < Q > (input[1].metric, input[112].metric);
    sc_uint<1> comp_r_1_113 = COMP < Q > (input[1].metric, input[113].metric);
    sc_uint<1> comp_r_1_114 = COMP < Q > (input[1].metric, input[114].metric);
    sc_uint<1> comp_r_1_115 = COMP < Q > (input[1].metric, input[115].metric);
    sc_uint<1> comp_r_1_116 = COMP < Q > (input[1].metric, input[116].metric);
    sc_uint<1> comp_r_1_117 = COMP < Q > (input[1].metric, input[117].metric);
    sc_uint<1> comp_r_1_118 = COMP < Q > (input[1].metric, input[118].metric);
    sc_uint<1> comp_r_1_119 = COMP < Q > (input[1].metric, input[119].metric);
    sc_uint<1> comp_r_1_120 = COMP < Q > (input[1].metric, input[120].metric);
    sc_uint<1> comp_r_1_121 = COMP < Q > (input[1].metric, input[121].metric);
    sc_uint<1> comp_r_1_122 = COMP < Q > (input[1].metric, input[122].metric);
    sc_uint<1> comp_r_1_123 = COMP < Q > (input[1].metric, input[123].metric);
    sc_uint<1> comp_r_1_124 = COMP < Q > (input[1].metric, input[124].metric);
    sc_uint<1> comp_r_1_125 = COMP < Q > (input[1].metric, input[125].metric);
    sc_uint<1> comp_r_1_126 = COMP < Q > (input[1].metric, input[126].metric);

    sc_uint<1> comp_r_3_4 = COMP < Q > (input[3].metric, input[4].metric);
    sc_uint<1> comp_r_3_5 = COMP < Q > (input[3].metric, input[5].metric);
    sc_uint<1> comp_r_3_6 = COMP < Q > (input[3].metric, input[6].metric);
    sc_uint<1> comp_r_3_7 = COMP < Q > (input[3].metric, input[7].metric);
    sc_uint<1> comp_r_3_8 = COMP < Q > (input[3].metric, input[8].metric);
    sc_uint<1> comp_r_3_9 = COMP < Q > (input[3].metric, input[9].metric);
    sc_uint<1> comp_r_3_10 = COMP < Q > (input[3].metric, input[10].metric);
    sc_uint<1> comp_r_3_11 = COMP < Q > (input[3].metric, input[11].metric);
    sc_uint<1> comp_r_3_12 = COMP < Q > (input[3].metric, input[12].metric);
    sc_uint<1> comp_r_3_13 = COMP < Q > (input[3].metric, input[13].metric);
    sc_uint<1> comp_r_3_14 = COMP < Q > (input[3].metric, input[14].metric);
    sc_uint<1> comp_r_3_15 = COMP < Q > (input[3].metric, input[15].metric);
    sc_uint<1> comp_r_3_16 = COMP < Q > (input[3].metric, input[16].metric);
    sc_uint<1> comp_r_3_17 = COMP < Q > (input[3].metric, input[17].metric);
    sc_uint<1> comp_r_3_18 = COMP < Q > (input[3].metric, input[18].metric);
    sc_uint<1> comp_r_3_19 = COMP < Q > (input[3].metric, input[19].metric);
    sc_uint<1> comp_r_3_20 = COMP < Q > (input[3].metric, input[20].metric);
    sc_uint<1> comp_r_3_21 = COMP < Q > (input[3].metric, input[21].metric);
    sc_uint<1> comp_r_3_22 = COMP < Q > (input[3].metric, input[22].metric);
    sc_uint<1> comp_r_3_23 = COMP < Q > (input[3].metric, input[23].metric);
    sc_uint<1> comp_r_3_24 = COMP < Q > (input[3].metric, input[24].metric);
    sc_uint<1> comp_r_3_25 = COMP < Q > (input[3].metric, input[25].metric);
    sc_uint<1> comp_r_3_26 = COMP < Q > (input[3].metric, input[26].metric);
    sc_uint<1> comp_r_3_27 = COMP < Q > (input[3].metric, input[27].metric);
    sc_uint<1> comp_r_3_28 = COMP < Q > (input[3].metric, input[28].metric);
    sc_uint<1> comp_r_3_29 = COMP < Q > (input[3].metric, input[29].metric);
    sc_uint<1> comp_r_3_30 = COMP < Q > (input[3].metric, input[30].metric);
    sc_uint<1> comp_r_3_31 = COMP < Q > (input[3].metric, input[31].metric);
    sc_uint<1> comp_r_3_32 = COMP < Q > (input[3].metric, input[32].metric);
    sc_uint<1> comp_r_3_33 = COMP < Q > (input[3].metric, input[33].metric);
    sc_uint<1> comp_r_3_34 = COMP < Q > (input[3].metric, input[34].metric);
    sc_uint<1> comp_r_3_35 = COMP < Q > (input[3].metric, input[35].metric);
    sc_uint<1> comp_r_3_36 = COMP < Q > (input[3].metric, input[36].metric);
    sc_uint<1> comp_r_3_37 = COMP < Q > (input[3].metric, input[37].metric);
    sc_uint<1> comp_r_3_38 = COMP < Q > (input[3].metric, input[38].metric);
    sc_uint<1> comp_r_3_39 = COMP < Q > (input[3].metric, input[39].metric);
    sc_uint<1> comp_r_3_40 = COMP < Q > (input[3].metric, input[40].metric);
    sc_uint<1> comp_r_3_41 = COMP < Q > (input[3].metric, input[41].metric);
    sc_uint<1> comp_r_3_42 = COMP < Q > (input[3].metric, input[42].metric);
    sc_uint<1> comp_r_3_43 = COMP < Q > (input[3].metric, input[43].metric);
    sc_uint<1> comp_r_3_44 = COMP < Q > (input[3].metric, input[44].metric);
    sc_uint<1> comp_r_3_45 = COMP < Q > (input[3].metric, input[45].metric);
    sc_uint<1> comp_r_3_46 = COMP < Q > (input[3].metric, input[46].metric);
    sc_uint<1> comp_r_3_47 = COMP < Q > (input[3].metric, input[47].metric);
    sc_uint<1> comp_r_3_48 = COMP < Q > (input[3].metric, input[48].metric);
    sc_uint<1> comp_r_3_49 = COMP < Q > (input[3].metric, input[49].metric);
    sc_uint<1> comp_r_3_50 = COMP < Q > (input[3].metric, input[50].metric);
    sc_uint<1> comp_r_3_51 = COMP < Q > (input[3].metric, input[51].metric);
    sc_uint<1> comp_r_3_52 = COMP < Q > (input[3].metric, input[52].metric);
    sc_uint<1> comp_r_3_53 = COMP < Q > (input[3].metric, input[53].metric);
    sc_uint<1> comp_r_3_54 = COMP < Q > (input[3].metric, input[54].metric);
    sc_uint<1> comp_r_3_55 = COMP < Q > (input[3].metric, input[55].metric);
    sc_uint<1> comp_r_3_56 = COMP < Q > (input[3].metric, input[56].metric);
    sc_uint<1> comp_r_3_57 = COMP < Q > (input[3].metric, input[57].metric);
    sc_uint<1> comp_r_3_58 = COMP < Q > (input[3].metric, input[58].metric);
    sc_uint<1> comp_r_3_59 = COMP < Q > (input[3].metric, input[59].metric);
    sc_uint<1> comp_r_3_60 = COMP < Q > (input[3].metric, input[60].metric);
    sc_uint<1> comp_r_3_61 = COMP < Q > (input[3].metric, input[61].metric);
    sc_uint<1> comp_r_3_62 = COMP < Q > (input[3].metric, input[62].metric);
    sc_uint<1> comp_r_3_63 = COMP < Q > (input[3].metric, input[63].metric);
    sc_uint<1> comp_r_3_64 = COMP < Q > (input[3].metric, input[64].metric);
    sc_uint<1> comp_r_3_65 = COMP < Q > (input[3].metric, input[65].metric);
    sc_uint<1> comp_r_3_66 = COMP < Q > (input[3].metric, input[66].metric);
    sc_uint<1> comp_r_3_67 = COMP < Q > (input[3].metric, input[67].metric);
    sc_uint<1> comp_r_3_68 = COMP < Q > (input[3].metric, input[68].metric);
    sc_uint<1> comp_r_3_69 = COMP < Q > (input[3].metric, input[69].metric);
    sc_uint<1> comp_r_3_70 = COMP < Q > (input[3].metric, input[70].metric);
    sc_uint<1> comp_r_3_71 = COMP < Q > (input[3].metric, input[71].metric);
    sc_uint<1> comp_r_3_72 = COMP < Q > (input[3].metric, input[72].metric);
    sc_uint<1> comp_r_3_73 = COMP < Q > (input[3].metric, input[73].metric);
    sc_uint<1> comp_r_3_74 = COMP < Q > (input[3].metric, input[74].metric);
    sc_uint<1> comp_r_3_75 = COMP < Q > (input[3].metric, input[75].metric);
    sc_uint<1> comp_r_3_76 = COMP < Q > (input[3].metric, input[76].metric);
    sc_uint<1> comp_r_3_77 = COMP < Q > (input[3].metric, input[77].metric);
    sc_uint<1> comp_r_3_78 = COMP < Q > (input[3].metric, input[78].metric);
    sc_uint<1> comp_r_3_79 = COMP < Q > (input[3].metric, input[79].metric);
    sc_uint<1> comp_r_3_80 = COMP < Q > (input[3].metric, input[80].metric);
    sc_uint<1> comp_r_3_81 = COMP < Q > (input[3].metric, input[81].metric);
    sc_uint<1> comp_r_3_82 = COMP < Q > (input[3].metric, input[82].metric);
    sc_uint<1> comp_r_3_83 = COMP < Q > (input[3].metric, input[83].metric);
    sc_uint<1> comp_r_3_84 = COMP < Q > (input[3].metric, input[84].metric);
    sc_uint<1> comp_r_3_85 = COMP < Q > (input[3].metric, input[85].metric);
    sc_uint<1> comp_r_3_86 = COMP < Q > (input[3].metric, input[86].metric);
    sc_uint<1> comp_r_3_87 = COMP < Q > (input[3].metric, input[87].metric);
    sc_uint<1> comp_r_3_88 = COMP < Q > (input[3].metric, input[88].metric);
    sc_uint<1> comp_r_3_89 = COMP < Q > (input[3].metric, input[89].metric);
    sc_uint<1> comp_r_3_90 = COMP < Q > (input[3].metric, input[90].metric);
    sc_uint<1> comp_r_3_91 = COMP < Q > (input[3].metric, input[91].metric);
    sc_uint<1> comp_r_3_92 = COMP < Q > (input[3].metric, input[92].metric);
    sc_uint<1> comp_r_3_93 = COMP < Q > (input[3].metric, input[93].metric);
    sc_uint<1> comp_r_3_94 = COMP < Q > (input[3].metric, input[94].metric);
    sc_uint<1> comp_r_3_95 = COMP < Q > (input[3].metric, input[95].metric);
    sc_uint<1> comp_r_3_96 = COMP < Q > (input[3].metric, input[96].metric);
    sc_uint<1> comp_r_3_97 = COMP < Q > (input[3].metric, input[97].metric);
    sc_uint<1> comp_r_3_98 = COMP < Q > (input[3].metric, input[98].metric);
    sc_uint<1> comp_r_3_99 = COMP < Q > (input[3].metric, input[99].metric);
    sc_uint<1> comp_r_3_100 = COMP < Q > (input[3].metric, input[100].metric);
    sc_uint<1> comp_r_3_101 = COMP < Q > (input[3].metric, input[101].metric);
    sc_uint<1> comp_r_3_102 = COMP < Q > (input[3].metric, input[102].metric);
    sc_uint<1> comp_r_3_103 = COMP < Q > (input[3].metric, input[103].metric);
    sc_uint<1> comp_r_3_104 = COMP < Q > (input[3].metric, input[104].metric);
    sc_uint<1> comp_r_3_105 = COMP < Q > (input[3].metric, input[105].metric);
    sc_uint<1> comp_r_3_106 = COMP < Q > (input[3].metric, input[106].metric);
    sc_uint<1> comp_r_3_107 = COMP < Q > (input[3].metric, input[107].metric);
    sc_uint<1> comp_r_3_108 = COMP < Q > (input[3].metric, input[108].metric);
    sc_uint<1> comp_r_3_109 = COMP < Q > (input[3].metric, input[109].metric);
    sc_uint<1> comp_r_3_110 = COMP < Q > (input[3].metric, input[110].metric);
    sc_uint<1> comp_r_3_111 = COMP < Q > (input[3].metric, input[111].metric);
    sc_uint<1> comp_r_3_112 = COMP < Q > (input[3].metric, input[112].metric);
    sc_uint<1> comp_r_3_113 = COMP < Q > (input[3].metric, input[113].metric);
    sc_uint<1> comp_r_3_114 = COMP < Q > (input[3].metric, input[114].metric);
    sc_uint<1> comp_r_3_115 = COMP < Q > (input[3].metric, input[115].metric);
    sc_uint<1> comp_r_3_116 = COMP < Q > (input[3].metric, input[116].metric);
    sc_uint<1> comp_r_3_117 = COMP < Q > (input[3].metric, input[117].metric);
    sc_uint<1> comp_r_3_118 = COMP < Q > (input[3].metric, input[118].metric);
    sc_uint<1> comp_r_3_119 = COMP < Q > (input[3].metric, input[119].metric);
    sc_uint<1> comp_r_3_120 = COMP < Q > (input[3].metric, input[120].metric);
    sc_uint<1> comp_r_3_121 = COMP < Q > (input[3].metric, input[121].metric);
    sc_uint<1> comp_r_3_122 = COMP < Q > (input[3].metric, input[122].metric);
    sc_uint<1> comp_r_3_123 = COMP < Q > (input[3].metric, input[123].metric);
    sc_uint<1> comp_r_3_124 = COMP < Q > (input[3].metric, input[124].metric);
    sc_uint<1> comp_r_3_125 = COMP < Q > (input[3].metric, input[125].metric);
    sc_uint<1> comp_r_3_126 = COMP < Q > (input[3].metric, input[126].metric);

    sc_uint<1> comp_r_5_6 = COMP < Q > (input[5].metric, input[6].metric);
    sc_uint<1> comp_r_5_7 = COMP < Q > (input[5].metric, input[7].metric);
    sc_uint<1> comp_r_5_8 = COMP < Q > (input[5].metric, input[8].metric);
    sc_uint<1> comp_r_5_9 = COMP < Q > (input[5].metric, input[9].metric);
    sc_uint<1> comp_r_5_10 = COMP < Q > (input[5].metric, input[10].metric);
    sc_uint<1> comp_r_5_11 = COMP < Q > (input[5].metric, input[11].metric);
    sc_uint<1> comp_r_5_12 = COMP < Q > (input[5].metric, input[12].metric);
    sc_uint<1> comp_r_5_13 = COMP < Q > (input[5].metric, input[13].metric);
    sc_uint<1> comp_r_5_14 = COMP < Q > (input[5].metric, input[14].metric);
    sc_uint<1> comp_r_5_15 = COMP < Q > (input[5].metric, input[15].metric);
    sc_uint<1> comp_r_5_16 = COMP < Q > (input[5].metric, input[16].metric);
    sc_uint<1> comp_r_5_17 = COMP < Q > (input[5].metric, input[17].metric);
    sc_uint<1> comp_r_5_18 = COMP < Q > (input[5].metric, input[18].metric);
    sc_uint<1> comp_r_5_19 = COMP < Q > (input[5].metric, input[19].metric);
    sc_uint<1> comp_r_5_20 = COMP < Q > (input[5].metric, input[20].metric);
    sc_uint<1> comp_r_5_21 = COMP < Q > (input[5].metric, input[21].metric);
    sc_uint<1> comp_r_5_22 = COMP < Q > (input[5].metric, input[22].metric);
    sc_uint<1> comp_r_5_23 = COMP < Q > (input[5].metric, input[23].metric);
    sc_uint<1> comp_r_5_24 = COMP < Q > (input[5].metric, input[24].metric);
    sc_uint<1> comp_r_5_25 = COMP < Q > (input[5].metric, input[25].metric);
    sc_uint<1> comp_r_5_26 = COMP < Q > (input[5].metric, input[26].metric);
    sc_uint<1> comp_r_5_27 = COMP < Q > (input[5].metric, input[27].metric);
    sc_uint<1> comp_r_5_28 = COMP < Q > (input[5].metric, input[28].metric);
    sc_uint<1> comp_r_5_29 = COMP < Q > (input[5].metric, input[29].metric);
    sc_uint<1> comp_r_5_30 = COMP < Q > (input[5].metric, input[30].metric);
    sc_uint<1> comp_r_5_31 = COMP < Q > (input[5].metric, input[31].metric);
    sc_uint<1> comp_r_5_32 = COMP < Q > (input[5].metric, input[32].metric);
    sc_uint<1> comp_r_5_33 = COMP < Q > (input[5].metric, input[33].metric);
    sc_uint<1> comp_r_5_34 = COMP < Q > (input[5].metric, input[34].metric);
    sc_uint<1> comp_r_5_35 = COMP < Q > (input[5].metric, input[35].metric);
    sc_uint<1> comp_r_5_36 = COMP < Q > (input[5].metric, input[36].metric);
    sc_uint<1> comp_r_5_37 = COMP < Q > (input[5].metric, input[37].metric);
    sc_uint<1> comp_r_5_38 = COMP < Q > (input[5].metric, input[38].metric);
    sc_uint<1> comp_r_5_39 = COMP < Q > (input[5].metric, input[39].metric);
    sc_uint<1> comp_r_5_40 = COMP < Q > (input[5].metric, input[40].metric);
    sc_uint<1> comp_r_5_41 = COMP < Q > (input[5].metric, input[41].metric);
    sc_uint<1> comp_r_5_42 = COMP < Q > (input[5].metric, input[42].metric);
    sc_uint<1> comp_r_5_43 = COMP < Q > (input[5].metric, input[43].metric);
    sc_uint<1> comp_r_5_44 = COMP < Q > (input[5].metric, input[44].metric);
    sc_uint<1> comp_r_5_45 = COMP < Q > (input[5].metric, input[45].metric);
    sc_uint<1> comp_r_5_46 = COMP < Q > (input[5].metric, input[46].metric);
    sc_uint<1> comp_r_5_47 = COMP < Q > (input[5].metric, input[47].metric);
    sc_uint<1> comp_r_5_48 = COMP < Q > (input[5].metric, input[48].metric);
    sc_uint<1> comp_r_5_49 = COMP < Q > (input[5].metric, input[49].metric);
    sc_uint<1> comp_r_5_50 = COMP < Q > (input[5].metric, input[50].metric);
    sc_uint<1> comp_r_5_51 = COMP < Q > (input[5].metric, input[51].metric);
    sc_uint<1> comp_r_5_52 = COMP < Q > (input[5].metric, input[52].metric);
    sc_uint<1> comp_r_5_53 = COMP < Q > (input[5].metric, input[53].metric);
    sc_uint<1> comp_r_5_54 = COMP < Q > (input[5].metric, input[54].metric);
    sc_uint<1> comp_r_5_55 = COMP < Q > (input[5].metric, input[55].metric);
    sc_uint<1> comp_r_5_56 = COMP < Q > (input[5].metric, input[56].metric);
    sc_uint<1> comp_r_5_57 = COMP < Q > (input[5].metric, input[57].metric);
    sc_uint<1> comp_r_5_58 = COMP < Q > (input[5].metric, input[58].metric);
    sc_uint<1> comp_r_5_59 = COMP < Q > (input[5].metric, input[59].metric);
    sc_uint<1> comp_r_5_60 = COMP < Q > (input[5].metric, input[60].metric);
    sc_uint<1> comp_r_5_61 = COMP < Q > (input[5].metric, input[61].metric);
    sc_uint<1> comp_r_5_62 = COMP < Q > (input[5].metric, input[62].metric);
    sc_uint<1> comp_r_5_63 = COMP < Q > (input[5].metric, input[63].metric);
    sc_uint<1> comp_r_5_64 = COMP < Q > (input[5].metric, input[64].metric);
    sc_uint<1> comp_r_5_65 = COMP < Q > (input[5].metric, input[65].metric);
    sc_uint<1> comp_r_5_66 = COMP < Q > (input[5].metric, input[66].metric);
    sc_uint<1> comp_r_5_67 = COMP < Q > (input[5].metric, input[67].metric);
    sc_uint<1> comp_r_5_68 = COMP < Q > (input[5].metric, input[68].metric);
    sc_uint<1> comp_r_5_69 = COMP < Q > (input[5].metric, input[69].metric);
    sc_uint<1> comp_r_5_70 = COMP < Q > (input[5].metric, input[70].metric);
    sc_uint<1> comp_r_5_71 = COMP < Q > (input[5].metric, input[71].metric);
    sc_uint<1> comp_r_5_72 = COMP < Q > (input[5].metric, input[72].metric);
    sc_uint<1> comp_r_5_73 = COMP < Q > (input[5].metric, input[73].metric);
    sc_uint<1> comp_r_5_74 = COMP < Q > (input[5].metric, input[74].metric);
    sc_uint<1> comp_r_5_75 = COMP < Q > (input[5].metric, input[75].metric);
    sc_uint<1> comp_r_5_76 = COMP < Q > (input[5].metric, input[76].metric);
    sc_uint<1> comp_r_5_77 = COMP < Q > (input[5].metric, input[77].metric);
    sc_uint<1> comp_r_5_78 = COMP < Q > (input[5].metric, input[78].metric);
    sc_uint<1> comp_r_5_79 = COMP < Q > (input[5].metric, input[79].metric);
    sc_uint<1> comp_r_5_80 = COMP < Q > (input[5].metric, input[80].metric);
    sc_uint<1> comp_r_5_81 = COMP < Q > (input[5].metric, input[81].metric);
    sc_uint<1> comp_r_5_82 = COMP < Q > (input[5].metric, input[82].metric);
    sc_uint<1> comp_r_5_83 = COMP < Q > (input[5].metric, input[83].metric);
    sc_uint<1> comp_r_5_84 = COMP < Q > (input[5].metric, input[84].metric);
    sc_uint<1> comp_r_5_85 = COMP < Q > (input[5].metric, input[85].metric);
    sc_uint<1> comp_r_5_86 = COMP < Q > (input[5].metric, input[86].metric);
    sc_uint<1> comp_r_5_87 = COMP < Q > (input[5].metric, input[87].metric);
    sc_uint<1> comp_r_5_88 = COMP < Q > (input[5].metric, input[88].metric);
    sc_uint<1> comp_r_5_89 = COMP < Q > (input[5].metric, input[89].metric);
    sc_uint<1> comp_r_5_90 = COMP < Q > (input[5].metric, input[90].metric);
    sc_uint<1> comp_r_5_91 = COMP < Q > (input[5].metric, input[91].metric);
    sc_uint<1> comp_r_5_92 = COMP < Q > (input[5].metric, input[92].metric);
    sc_uint<1> comp_r_5_93 = COMP < Q > (input[5].metric, input[93].metric);
    sc_uint<1> comp_r_5_94 = COMP < Q > (input[5].metric, input[94].metric);
    sc_uint<1> comp_r_5_95 = COMP < Q > (input[5].metric, input[95].metric);
    sc_uint<1> comp_r_5_96 = COMP < Q > (input[5].metric, input[96].metric);
    sc_uint<1> comp_r_5_97 = COMP < Q > (input[5].metric, input[97].metric);
    sc_uint<1> comp_r_5_98 = COMP < Q > (input[5].metric, input[98].metric);
    sc_uint<1> comp_r_5_99 = COMP < Q > (input[5].metric, input[99].metric);
    sc_uint<1> comp_r_5_100 = COMP < Q > (input[5].metric, input[100].metric);
    sc_uint<1> comp_r_5_101 = COMP < Q > (input[5].metric, input[101].metric);
    sc_uint<1> comp_r_5_102 = COMP < Q > (input[5].metric, input[102].metric);
    sc_uint<1> comp_r_5_103 = COMP < Q > (input[5].metric, input[103].metric);
    sc_uint<1> comp_r_5_104 = COMP < Q > (input[5].metric, input[104].metric);
    sc_uint<1> comp_r_5_105 = COMP < Q > (input[5].metric, input[105].metric);
    sc_uint<1> comp_r_5_106 = COMP < Q > (input[5].metric, input[106].metric);
    sc_uint<1> comp_r_5_107 = COMP < Q > (input[5].metric, input[107].metric);
    sc_uint<1> comp_r_5_108 = COMP < Q > (input[5].metric, input[108].metric);
    sc_uint<1> comp_r_5_109 = COMP < Q > (input[5].metric, input[109].metric);
    sc_uint<1> comp_r_5_110 = COMP < Q > (input[5].metric, input[110].metric);
    sc_uint<1> comp_r_5_111 = COMP < Q > (input[5].metric, input[111].metric);
    sc_uint<1> comp_r_5_112 = COMP < Q > (input[5].metric, input[112].metric);
    sc_uint<1> comp_r_5_113 = COMP < Q > (input[5].metric, input[113].metric);
    sc_uint<1> comp_r_5_114 = COMP < Q > (input[5].metric, input[114].metric);
    sc_uint<1> comp_r_5_115 = COMP < Q > (input[5].metric, input[115].metric);
    sc_uint<1> comp_r_5_116 = COMP < Q > (input[5].metric, input[116].metric);
    sc_uint<1> comp_r_5_117 = COMP < Q > (input[5].metric, input[117].metric);
    sc_uint<1> comp_r_5_118 = COMP < Q > (input[5].metric, input[118].metric);
    sc_uint<1> comp_r_5_119 = COMP < Q > (input[5].metric, input[119].metric);
    sc_uint<1> comp_r_5_120 = COMP < Q > (input[5].metric, input[120].metric);
    sc_uint<1> comp_r_5_121 = COMP < Q > (input[5].metric, input[121].metric);
    sc_uint<1> comp_r_5_122 = COMP < Q > (input[5].metric, input[122].metric);
    sc_uint<1> comp_r_5_123 = COMP < Q > (input[5].metric, input[123].metric);
    sc_uint<1> comp_r_5_124 = COMP < Q > (input[5].metric, input[124].metric);
    sc_uint<1> comp_r_5_125 = COMP < Q > (input[5].metric, input[125].metric);
    sc_uint<1> comp_r_5_126 = COMP < Q > (input[5].metric, input[126].metric);

    sc_uint<1> comp_r_7_8 = COMP < Q > (input[7].metric, input[8].metric);
    sc_uint<1> comp_r_7_9 = COMP < Q > (input[7].metric, input[9].metric);
    sc_uint<1> comp_r_7_10 = COMP < Q > (input[7].metric, input[10].metric);
    sc_uint<1> comp_r_7_11 = COMP < Q > (input[7].metric, input[11].metric);
    sc_uint<1> comp_r_7_12 = COMP < Q > (input[7].metric, input[12].metric);
    sc_uint<1> comp_r_7_13 = COMP < Q > (input[7].metric, input[13].metric);
    sc_uint<1> comp_r_7_14 = COMP < Q > (input[7].metric, input[14].metric);
    sc_uint<1> comp_r_7_15 = COMP < Q > (input[7].metric, input[15].metric);
    sc_uint<1> comp_r_7_16 = COMP < Q > (input[7].metric, input[16].metric);
    sc_uint<1> comp_r_7_17 = COMP < Q > (input[7].metric, input[17].metric);
    sc_uint<1> comp_r_7_18 = COMP < Q > (input[7].metric, input[18].metric);
    sc_uint<1> comp_r_7_19 = COMP < Q > (input[7].metric, input[19].metric);
    sc_uint<1> comp_r_7_20 = COMP < Q > (input[7].metric, input[20].metric);
    sc_uint<1> comp_r_7_21 = COMP < Q > (input[7].metric, input[21].metric);
    sc_uint<1> comp_r_7_22 = COMP < Q > (input[7].metric, input[22].metric);
    sc_uint<1> comp_r_7_23 = COMP < Q > (input[7].metric, input[23].metric);
    sc_uint<1> comp_r_7_24 = COMP < Q > (input[7].metric, input[24].metric);
    sc_uint<1> comp_r_7_25 = COMP < Q > (input[7].metric, input[25].metric);
    sc_uint<1> comp_r_7_26 = COMP < Q > (input[7].metric, input[26].metric);
    sc_uint<1> comp_r_7_27 = COMP < Q > (input[7].metric, input[27].metric);
    sc_uint<1> comp_r_7_28 = COMP < Q > (input[7].metric, input[28].metric);
    sc_uint<1> comp_r_7_29 = COMP < Q > (input[7].metric, input[29].metric);
    sc_uint<1> comp_r_7_30 = COMP < Q > (input[7].metric, input[30].metric);
    sc_uint<1> comp_r_7_31 = COMP < Q > (input[7].metric, input[31].metric);
    sc_uint<1> comp_r_7_32 = COMP < Q > (input[7].metric, input[32].metric);
    sc_uint<1> comp_r_7_33 = COMP < Q > (input[7].metric, input[33].metric);
    sc_uint<1> comp_r_7_34 = COMP < Q > (input[7].metric, input[34].metric);
    sc_uint<1> comp_r_7_35 = COMP < Q > (input[7].metric, input[35].metric);
    sc_uint<1> comp_r_7_36 = COMP < Q > (input[7].metric, input[36].metric);
    sc_uint<1> comp_r_7_37 = COMP < Q > (input[7].metric, input[37].metric);
    sc_uint<1> comp_r_7_38 = COMP < Q > (input[7].metric, input[38].metric);
    sc_uint<1> comp_r_7_39 = COMP < Q > (input[7].metric, input[39].metric);
    sc_uint<1> comp_r_7_40 = COMP < Q > (input[7].metric, input[40].metric);
    sc_uint<1> comp_r_7_41 = COMP < Q > (input[7].metric, input[41].metric);
    sc_uint<1> comp_r_7_42 = COMP < Q > (input[7].metric, input[42].metric);
    sc_uint<1> comp_r_7_43 = COMP < Q > (input[7].metric, input[43].metric);
    sc_uint<1> comp_r_7_44 = COMP < Q > (input[7].metric, input[44].metric);
    sc_uint<1> comp_r_7_45 = COMP < Q > (input[7].metric, input[45].metric);
    sc_uint<1> comp_r_7_46 = COMP < Q > (input[7].metric, input[46].metric);
    sc_uint<1> comp_r_7_47 = COMP < Q > (input[7].metric, input[47].metric);
    sc_uint<1> comp_r_7_48 = COMP < Q > (input[7].metric, input[48].metric);
    sc_uint<1> comp_r_7_49 = COMP < Q > (input[7].metric, input[49].metric);
    sc_uint<1> comp_r_7_50 = COMP < Q > (input[7].metric, input[50].metric);
    sc_uint<1> comp_r_7_51 = COMP < Q > (input[7].metric, input[51].metric);
    sc_uint<1> comp_r_7_52 = COMP < Q > (input[7].metric, input[52].metric);
    sc_uint<1> comp_r_7_53 = COMP < Q > (input[7].metric, input[53].metric);
    sc_uint<1> comp_r_7_54 = COMP < Q > (input[7].metric, input[54].metric);
    sc_uint<1> comp_r_7_55 = COMP < Q > (input[7].metric, input[55].metric);
    sc_uint<1> comp_r_7_56 = COMP < Q > (input[7].metric, input[56].metric);
    sc_uint<1> comp_r_7_57 = COMP < Q > (input[7].metric, input[57].metric);
    sc_uint<1> comp_r_7_58 = COMP < Q > (input[7].metric, input[58].metric);
    sc_uint<1> comp_r_7_59 = COMP < Q > (input[7].metric, input[59].metric);
    sc_uint<1> comp_r_7_60 = COMP < Q > (input[7].metric, input[60].metric);
    sc_uint<1> comp_r_7_61 = COMP < Q > (input[7].metric, input[61].metric);
    sc_uint<1> comp_r_7_62 = COMP < Q > (input[7].metric, input[62].metric);
    sc_uint<1> comp_r_7_63 = COMP < Q > (input[7].metric, input[63].metric);
    sc_uint<1> comp_r_7_64 = COMP < Q > (input[7].metric, input[64].metric);
    sc_uint<1> comp_r_7_65 = COMP < Q > (input[7].metric, input[65].metric);
    sc_uint<1> comp_r_7_66 = COMP < Q > (input[7].metric, input[66].metric);
    sc_uint<1> comp_r_7_67 = COMP < Q > (input[7].metric, input[67].metric);
    sc_uint<1> comp_r_7_68 = COMP < Q > (input[7].metric, input[68].metric);
    sc_uint<1> comp_r_7_69 = COMP < Q > (input[7].metric, input[69].metric);
    sc_uint<1> comp_r_7_70 = COMP < Q > (input[7].metric, input[70].metric);
    sc_uint<1> comp_r_7_71 = COMP < Q > (input[7].metric, input[71].metric);
    sc_uint<1> comp_r_7_72 = COMP < Q > (input[7].metric, input[72].metric);
    sc_uint<1> comp_r_7_73 = COMP < Q > (input[7].metric, input[73].metric);
    sc_uint<1> comp_r_7_74 = COMP < Q > (input[7].metric, input[74].metric);
    sc_uint<1> comp_r_7_75 = COMP < Q > (input[7].metric, input[75].metric);
    sc_uint<1> comp_r_7_76 = COMP < Q > (input[7].metric, input[76].metric);
    sc_uint<1> comp_r_7_77 = COMP < Q > (input[7].metric, input[77].metric);
    sc_uint<1> comp_r_7_78 = COMP < Q > (input[7].metric, input[78].metric);
    sc_uint<1> comp_r_7_79 = COMP < Q > (input[7].metric, input[79].metric);
    sc_uint<1> comp_r_7_80 = COMP < Q > (input[7].metric, input[80].metric);
    sc_uint<1> comp_r_7_81 = COMP < Q > (input[7].metric, input[81].metric);
    sc_uint<1> comp_r_7_82 = COMP < Q > (input[7].metric, input[82].metric);
    sc_uint<1> comp_r_7_83 = COMP < Q > (input[7].metric, input[83].metric);
    sc_uint<1> comp_r_7_84 = COMP < Q > (input[7].metric, input[84].metric);
    sc_uint<1> comp_r_7_85 = COMP < Q > (input[7].metric, input[85].metric);
    sc_uint<1> comp_r_7_86 = COMP < Q > (input[7].metric, input[86].metric);
    sc_uint<1> comp_r_7_87 = COMP < Q > (input[7].metric, input[87].metric);
    sc_uint<1> comp_r_7_88 = COMP < Q > (input[7].metric, input[88].metric);
    sc_uint<1> comp_r_7_89 = COMP < Q > (input[7].metric, input[89].metric);
    sc_uint<1> comp_r_7_90 = COMP < Q > (input[7].metric, input[90].metric);
    sc_uint<1> comp_r_7_91 = COMP < Q > (input[7].metric, input[91].metric);
    sc_uint<1> comp_r_7_92 = COMP < Q > (input[7].metric, input[92].metric);
    sc_uint<1> comp_r_7_93 = COMP < Q > (input[7].metric, input[93].metric);
    sc_uint<1> comp_r_7_94 = COMP < Q > (input[7].metric, input[94].metric);
    sc_uint<1> comp_r_7_95 = COMP < Q > (input[7].metric, input[95].metric);
    sc_uint<1> comp_r_7_96 = COMP < Q > (input[7].metric, input[96].metric);
    sc_uint<1> comp_r_7_97 = COMP < Q > (input[7].metric, input[97].metric);
    sc_uint<1> comp_r_7_98 = COMP < Q > (input[7].metric, input[98].metric);
    sc_uint<1> comp_r_7_99 = COMP < Q > (input[7].metric, input[99].metric);
    sc_uint<1> comp_r_7_100 = COMP < Q > (input[7].metric, input[100].metric);
    sc_uint<1> comp_r_7_101 = COMP < Q > (input[7].metric, input[101].metric);
    sc_uint<1> comp_r_7_102 = COMP < Q > (input[7].metric, input[102].metric);
    sc_uint<1> comp_r_7_103 = COMP < Q > (input[7].metric, input[103].metric);
    sc_uint<1> comp_r_7_104 = COMP < Q > (input[7].metric, input[104].metric);
    sc_uint<1> comp_r_7_105 = COMP < Q > (input[7].metric, input[105].metric);
    sc_uint<1> comp_r_7_106 = COMP < Q > (input[7].metric, input[106].metric);
    sc_uint<1> comp_r_7_107 = COMP < Q > (input[7].metric, input[107].metric);
    sc_uint<1> comp_r_7_108 = COMP < Q > (input[7].metric, input[108].metric);
    sc_uint<1> comp_r_7_109 = COMP < Q > (input[7].metric, input[109].metric);
    sc_uint<1> comp_r_7_110 = COMP < Q > (input[7].metric, input[110].metric);
    sc_uint<1> comp_r_7_111 = COMP < Q > (input[7].metric, input[111].metric);
    sc_uint<1> comp_r_7_112 = COMP < Q > (input[7].metric, input[112].metric);
    sc_uint<1> comp_r_7_113 = COMP < Q > (input[7].metric, input[113].metric);
    sc_uint<1> comp_r_7_114 = COMP < Q > (input[7].metric, input[114].metric);
    sc_uint<1> comp_r_7_115 = COMP < Q > (input[7].metric, input[115].metric);
    sc_uint<1> comp_r_7_116 = COMP < Q > (input[7].metric, input[116].metric);
    sc_uint<1> comp_r_7_117 = COMP < Q > (input[7].metric, input[117].metric);
    sc_uint<1> comp_r_7_118 = COMP < Q > (input[7].metric, input[118].metric);
    sc_uint<1> comp_r_7_119 = COMP < Q > (input[7].metric, input[119].metric);
    sc_uint<1> comp_r_7_120 = COMP < Q > (input[7].metric, input[120].metric);
    sc_uint<1> comp_r_7_121 = COMP < Q > (input[7].metric, input[121].metric);
    sc_uint<1> comp_r_7_122 = COMP < Q > (input[7].metric, input[122].metric);
    sc_uint<1> comp_r_7_123 = COMP < Q > (input[7].metric, input[123].metric);
    sc_uint<1> comp_r_7_124 = COMP < Q > (input[7].metric, input[124].metric);
    sc_uint<1> comp_r_7_125 = COMP < Q > (input[7].metric, input[125].metric);
    sc_uint<1> comp_r_7_126 = COMP < Q > (input[7].metric, input[126].metric);

    sc_uint<1> comp_r_9_10 = COMP < Q > (input[9].metric, input[10].metric);
    sc_uint<1> comp_r_9_11 = COMP < Q > (input[9].metric, input[11].metric);
    sc_uint<1> comp_r_9_12 = COMP < Q > (input[9].metric, input[12].metric);
    sc_uint<1> comp_r_9_13 = COMP < Q > (input[9].metric, input[13].metric);
    sc_uint<1> comp_r_9_14 = COMP < Q > (input[9].metric, input[14].metric);
    sc_uint<1> comp_r_9_15 = COMP < Q > (input[9].metric, input[15].metric);
    sc_uint<1> comp_r_9_16 = COMP < Q > (input[9].metric, input[16].metric);
    sc_uint<1> comp_r_9_17 = COMP < Q > (input[9].metric, input[17].metric);
    sc_uint<1> comp_r_9_18 = COMP < Q > (input[9].metric, input[18].metric);
    sc_uint<1> comp_r_9_19 = COMP < Q > (input[9].metric, input[19].metric);
    sc_uint<1> comp_r_9_20 = COMP < Q > (input[9].metric, input[20].metric);
    sc_uint<1> comp_r_9_21 = COMP < Q > (input[9].metric, input[21].metric);
    sc_uint<1> comp_r_9_22 = COMP < Q > (input[9].metric, input[22].metric);
    sc_uint<1> comp_r_9_23 = COMP < Q > (input[9].metric, input[23].metric);
    sc_uint<1> comp_r_9_24 = COMP < Q > (input[9].metric, input[24].metric);
    sc_uint<1> comp_r_9_25 = COMP < Q > (input[9].metric, input[25].metric);
    sc_uint<1> comp_r_9_26 = COMP < Q > (input[9].metric, input[26].metric);
    sc_uint<1> comp_r_9_27 = COMP < Q > (input[9].metric, input[27].metric);
    sc_uint<1> comp_r_9_28 = COMP < Q > (input[9].metric, input[28].metric);
    sc_uint<1> comp_r_9_29 = COMP < Q > (input[9].metric, input[29].metric);
    sc_uint<1> comp_r_9_30 = COMP < Q > (input[9].metric, input[30].metric);
    sc_uint<1> comp_r_9_31 = COMP < Q > (input[9].metric, input[31].metric);
    sc_uint<1> comp_r_9_32 = COMP < Q > (input[9].metric, input[32].metric);
    sc_uint<1> comp_r_9_33 = COMP < Q > (input[9].metric, input[33].metric);
    sc_uint<1> comp_r_9_34 = COMP < Q > (input[9].metric, input[34].metric);
    sc_uint<1> comp_r_9_35 = COMP < Q > (input[9].metric, input[35].metric);
    sc_uint<1> comp_r_9_36 = COMP < Q > (input[9].metric, input[36].metric);
    sc_uint<1> comp_r_9_37 = COMP < Q > (input[9].metric, input[37].metric);
    sc_uint<1> comp_r_9_38 = COMP < Q > (input[9].metric, input[38].metric);
    sc_uint<1> comp_r_9_39 = COMP < Q > (input[9].metric, input[39].metric);
    sc_uint<1> comp_r_9_40 = COMP < Q > (input[9].metric, input[40].metric);
    sc_uint<1> comp_r_9_41 = COMP < Q > (input[9].metric, input[41].metric);
    sc_uint<1> comp_r_9_42 = COMP < Q > (input[9].metric, input[42].metric);
    sc_uint<1> comp_r_9_43 = COMP < Q > (input[9].metric, input[43].metric);
    sc_uint<1> comp_r_9_44 = COMP < Q > (input[9].metric, input[44].metric);
    sc_uint<1> comp_r_9_45 = COMP < Q > (input[9].metric, input[45].metric);
    sc_uint<1> comp_r_9_46 = COMP < Q > (input[9].metric, input[46].metric);
    sc_uint<1> comp_r_9_47 = COMP < Q > (input[9].metric, input[47].metric);
    sc_uint<1> comp_r_9_48 = COMP < Q > (input[9].metric, input[48].metric);
    sc_uint<1> comp_r_9_49 = COMP < Q > (input[9].metric, input[49].metric);
    sc_uint<1> comp_r_9_50 = COMP < Q > (input[9].metric, input[50].metric);
    sc_uint<1> comp_r_9_51 = COMP < Q > (input[9].metric, input[51].metric);
    sc_uint<1> comp_r_9_52 = COMP < Q > (input[9].metric, input[52].metric);
    sc_uint<1> comp_r_9_53 = COMP < Q > (input[9].metric, input[53].metric);
    sc_uint<1> comp_r_9_54 = COMP < Q > (input[9].metric, input[54].metric);
    sc_uint<1> comp_r_9_55 = COMP < Q > (input[9].metric, input[55].metric);
    sc_uint<1> comp_r_9_56 = COMP < Q > (input[9].metric, input[56].metric);
    sc_uint<1> comp_r_9_57 = COMP < Q > (input[9].metric, input[57].metric);
    sc_uint<1> comp_r_9_58 = COMP < Q > (input[9].metric, input[58].metric);
    sc_uint<1> comp_r_9_59 = COMP < Q > (input[9].metric, input[59].metric);
    sc_uint<1> comp_r_9_60 = COMP < Q > (input[9].metric, input[60].metric);
    sc_uint<1> comp_r_9_61 = COMP < Q > (input[9].metric, input[61].metric);
    sc_uint<1> comp_r_9_62 = COMP < Q > (input[9].metric, input[62].metric);
    sc_uint<1> comp_r_9_63 = COMP < Q > (input[9].metric, input[63].metric);
    sc_uint<1> comp_r_9_64 = COMP < Q > (input[9].metric, input[64].metric);
    sc_uint<1> comp_r_9_65 = COMP < Q > (input[9].metric, input[65].metric);
    sc_uint<1> comp_r_9_66 = COMP < Q > (input[9].metric, input[66].metric);
    sc_uint<1> comp_r_9_67 = COMP < Q > (input[9].metric, input[67].metric);
    sc_uint<1> comp_r_9_68 = COMP < Q > (input[9].metric, input[68].metric);
    sc_uint<1> comp_r_9_69 = COMP < Q > (input[9].metric, input[69].metric);
    sc_uint<1> comp_r_9_70 = COMP < Q > (input[9].metric, input[70].metric);
    sc_uint<1> comp_r_9_71 = COMP < Q > (input[9].metric, input[71].metric);
    sc_uint<1> comp_r_9_72 = COMP < Q > (input[9].metric, input[72].metric);
    sc_uint<1> comp_r_9_73 = COMP < Q > (input[9].metric, input[73].metric);
    sc_uint<1> comp_r_9_74 = COMP < Q > (input[9].metric, input[74].metric);
    sc_uint<1> comp_r_9_75 = COMP < Q > (input[9].metric, input[75].metric);
    sc_uint<1> comp_r_9_76 = COMP < Q > (input[9].metric, input[76].metric);
    sc_uint<1> comp_r_9_77 = COMP < Q > (input[9].metric, input[77].metric);
    sc_uint<1> comp_r_9_78 = COMP < Q > (input[9].metric, input[78].metric);
    sc_uint<1> comp_r_9_79 = COMP < Q > (input[9].metric, input[79].metric);
    sc_uint<1> comp_r_9_80 = COMP < Q > (input[9].metric, input[80].metric);
    sc_uint<1> comp_r_9_81 = COMP < Q > (input[9].metric, input[81].metric);
    sc_uint<1> comp_r_9_82 = COMP < Q > (input[9].metric, input[82].metric);
    sc_uint<1> comp_r_9_83 = COMP < Q > (input[9].metric, input[83].metric);
    sc_uint<1> comp_r_9_84 = COMP < Q > (input[9].metric, input[84].metric);
    sc_uint<1> comp_r_9_85 = COMP < Q > (input[9].metric, input[85].metric);
    sc_uint<1> comp_r_9_86 = COMP < Q > (input[9].metric, input[86].metric);
    sc_uint<1> comp_r_9_87 = COMP < Q > (input[9].metric, input[87].metric);
    sc_uint<1> comp_r_9_88 = COMP < Q > (input[9].metric, input[88].metric);
    sc_uint<1> comp_r_9_89 = COMP < Q > (input[9].metric, input[89].metric);
    sc_uint<1> comp_r_9_90 = COMP < Q > (input[9].metric, input[90].metric);
    sc_uint<1> comp_r_9_91 = COMP < Q > (input[9].metric, input[91].metric);
    sc_uint<1> comp_r_9_92 = COMP < Q > (input[9].metric, input[92].metric);
    sc_uint<1> comp_r_9_93 = COMP < Q > (input[9].metric, input[93].metric);
    sc_uint<1> comp_r_9_94 = COMP < Q > (input[9].metric, input[94].metric);
    sc_uint<1> comp_r_9_95 = COMP < Q > (input[9].metric, input[95].metric);
    sc_uint<1> comp_r_9_96 = COMP < Q > (input[9].metric, input[96].metric);
    sc_uint<1> comp_r_9_97 = COMP < Q > (input[9].metric, input[97].metric);
    sc_uint<1> comp_r_9_98 = COMP < Q > (input[9].metric, input[98].metric);
    sc_uint<1> comp_r_9_99 = COMP < Q > (input[9].metric, input[99].metric);
    sc_uint<1> comp_r_9_100 = COMP < Q > (input[9].metric, input[100].metric);
    sc_uint<1> comp_r_9_101 = COMP < Q > (input[9].metric, input[101].metric);
    sc_uint<1> comp_r_9_102 = COMP < Q > (input[9].metric, input[102].metric);
    sc_uint<1> comp_r_9_103 = COMP < Q > (input[9].metric, input[103].metric);
    sc_uint<1> comp_r_9_104 = COMP < Q > (input[9].metric, input[104].metric);
    sc_uint<1> comp_r_9_105 = COMP < Q > (input[9].metric, input[105].metric);
    sc_uint<1> comp_r_9_106 = COMP < Q > (input[9].metric, input[106].metric);
    sc_uint<1> comp_r_9_107 = COMP < Q > (input[9].metric, input[107].metric);
    sc_uint<1> comp_r_9_108 = COMP < Q > (input[9].metric, input[108].metric);
    sc_uint<1> comp_r_9_109 = COMP < Q > (input[9].metric, input[109].metric);
    sc_uint<1> comp_r_9_110 = COMP < Q > (input[9].metric, input[110].metric);
    sc_uint<1> comp_r_9_111 = COMP < Q > (input[9].metric, input[111].metric);
    sc_uint<1> comp_r_9_112 = COMP < Q > (input[9].metric, input[112].metric);
    sc_uint<1> comp_r_9_113 = COMP < Q > (input[9].metric, input[113].metric);
    sc_uint<1> comp_r_9_114 = COMP < Q > (input[9].metric, input[114].metric);
    sc_uint<1> comp_r_9_115 = COMP < Q > (input[9].metric, input[115].metric);
    sc_uint<1> comp_r_9_116 = COMP < Q > (input[9].metric, input[116].metric);
    sc_uint<1> comp_r_9_117 = COMP < Q > (input[9].metric, input[117].metric);
    sc_uint<1> comp_r_9_118 = COMP < Q > (input[9].metric, input[118].metric);
    sc_uint<1> comp_r_9_119 = COMP < Q > (input[9].metric, input[119].metric);
    sc_uint<1> comp_r_9_120 = COMP < Q > (input[9].metric, input[120].metric);
    sc_uint<1> comp_r_9_121 = COMP < Q > (input[9].metric, input[121].metric);
    sc_uint<1> comp_r_9_122 = COMP < Q > (input[9].metric, input[122].metric);
    sc_uint<1> comp_r_9_123 = COMP < Q > (input[9].metric, input[123].metric);
    sc_uint<1> comp_r_9_124 = COMP < Q > (input[9].metric, input[124].metric);
    sc_uint<1> comp_r_9_125 = COMP < Q > (input[9].metric, input[125].metric);
    sc_uint<1> comp_r_9_126 = COMP < Q > (input[9].metric, input[126].metric);

    sc_uint<1> comp_r_11_12 = COMP < Q > (input[11].metric, input[12].metric);
    sc_uint<1> comp_r_11_13 = COMP < Q > (input[11].metric, input[13].metric);
    sc_uint<1> comp_r_11_14 = COMP < Q > (input[11].metric, input[14].metric);
    sc_uint<1> comp_r_11_15 = COMP < Q > (input[11].metric, input[15].metric);
    sc_uint<1> comp_r_11_16 = COMP < Q > (input[11].metric, input[16].metric);
    sc_uint<1> comp_r_11_17 = COMP < Q > (input[11].metric, input[17].metric);
    sc_uint<1> comp_r_11_18 = COMP < Q > (input[11].metric, input[18].metric);
    sc_uint<1> comp_r_11_19 = COMP < Q > (input[11].metric, input[19].metric);
    sc_uint<1> comp_r_11_20 = COMP < Q > (input[11].metric, input[20].metric);
    sc_uint<1> comp_r_11_21 = COMP < Q > (input[11].metric, input[21].metric);
    sc_uint<1> comp_r_11_22 = COMP < Q > (input[11].metric, input[22].metric);
    sc_uint<1> comp_r_11_23 = COMP < Q > (input[11].metric, input[23].metric);
    sc_uint<1> comp_r_11_24 = COMP < Q > (input[11].metric, input[24].metric);
    sc_uint<1> comp_r_11_25 = COMP < Q > (input[11].metric, input[25].metric);
    sc_uint<1> comp_r_11_26 = COMP < Q > (input[11].metric, input[26].metric);
    sc_uint<1> comp_r_11_27 = COMP < Q > (input[11].metric, input[27].metric);
    sc_uint<1> comp_r_11_28 = COMP < Q > (input[11].metric, input[28].metric);
    sc_uint<1> comp_r_11_29 = COMP < Q > (input[11].metric, input[29].metric);
    sc_uint<1> comp_r_11_30 = COMP < Q > (input[11].metric, input[30].metric);
    sc_uint<1> comp_r_11_31 = COMP < Q > (input[11].metric, input[31].metric);
    sc_uint<1> comp_r_11_32 = COMP < Q > (input[11].metric, input[32].metric);
    sc_uint<1> comp_r_11_33 = COMP < Q > (input[11].metric, input[33].metric);
    sc_uint<1> comp_r_11_34 = COMP < Q > (input[11].metric, input[34].metric);
    sc_uint<1> comp_r_11_35 = COMP < Q > (input[11].metric, input[35].metric);
    sc_uint<1> comp_r_11_36 = COMP < Q > (input[11].metric, input[36].metric);
    sc_uint<1> comp_r_11_37 = COMP < Q > (input[11].metric, input[37].metric);
    sc_uint<1> comp_r_11_38 = COMP < Q > (input[11].metric, input[38].metric);
    sc_uint<1> comp_r_11_39 = COMP < Q > (input[11].metric, input[39].metric);
    sc_uint<1> comp_r_11_40 = COMP < Q > (input[11].metric, input[40].metric);
    sc_uint<1> comp_r_11_41 = COMP < Q > (input[11].metric, input[41].metric);
    sc_uint<1> comp_r_11_42 = COMP < Q > (input[11].metric, input[42].metric);
    sc_uint<1> comp_r_11_43 = COMP < Q > (input[11].metric, input[43].metric);
    sc_uint<1> comp_r_11_44 = COMP < Q > (input[11].metric, input[44].metric);
    sc_uint<1> comp_r_11_45 = COMP < Q > (input[11].metric, input[45].metric);
    sc_uint<1> comp_r_11_46 = COMP < Q > (input[11].metric, input[46].metric);
    sc_uint<1> comp_r_11_47 = COMP < Q > (input[11].metric, input[47].metric);
    sc_uint<1> comp_r_11_48 = COMP < Q > (input[11].metric, input[48].metric);
    sc_uint<1> comp_r_11_49 = COMP < Q > (input[11].metric, input[49].metric);
    sc_uint<1> comp_r_11_50 = COMP < Q > (input[11].metric, input[50].metric);
    sc_uint<1> comp_r_11_51 = COMP < Q > (input[11].metric, input[51].metric);
    sc_uint<1> comp_r_11_52 = COMP < Q > (input[11].metric, input[52].metric);
    sc_uint<1> comp_r_11_53 = COMP < Q > (input[11].metric, input[53].metric);
    sc_uint<1> comp_r_11_54 = COMP < Q > (input[11].metric, input[54].metric);
    sc_uint<1> comp_r_11_55 = COMP < Q > (input[11].metric, input[55].metric);
    sc_uint<1> comp_r_11_56 = COMP < Q > (input[11].metric, input[56].metric);
    sc_uint<1> comp_r_11_57 = COMP < Q > (input[11].metric, input[57].metric);
    sc_uint<1> comp_r_11_58 = COMP < Q > (input[11].metric, input[58].metric);
    sc_uint<1> comp_r_11_59 = COMP < Q > (input[11].metric, input[59].metric);
    sc_uint<1> comp_r_11_60 = COMP < Q > (input[11].metric, input[60].metric);
    sc_uint<1> comp_r_11_61 = COMP < Q > (input[11].metric, input[61].metric);
    sc_uint<1> comp_r_11_62 = COMP < Q > (input[11].metric, input[62].metric);
    sc_uint<1> comp_r_11_63 = COMP < Q > (input[11].metric, input[63].metric);
    sc_uint<1> comp_r_11_64 = COMP < Q > (input[11].metric, input[64].metric);
    sc_uint<1> comp_r_11_65 = COMP < Q > (input[11].metric, input[65].metric);
    sc_uint<1> comp_r_11_66 = COMP < Q > (input[11].metric, input[66].metric);
    sc_uint<1> comp_r_11_67 = COMP < Q > (input[11].metric, input[67].metric);
    sc_uint<1> comp_r_11_68 = COMP < Q > (input[11].metric, input[68].metric);
    sc_uint<1> comp_r_11_69 = COMP < Q > (input[11].metric, input[69].metric);
    sc_uint<1> comp_r_11_70 = COMP < Q > (input[11].metric, input[70].metric);
    sc_uint<1> comp_r_11_71 = COMP < Q > (input[11].metric, input[71].metric);
    sc_uint<1> comp_r_11_72 = COMP < Q > (input[11].metric, input[72].metric);
    sc_uint<1> comp_r_11_73 = COMP < Q > (input[11].metric, input[73].metric);
    sc_uint<1> comp_r_11_74 = COMP < Q > (input[11].metric, input[74].metric);
    sc_uint<1> comp_r_11_75 = COMP < Q > (input[11].metric, input[75].metric);
    sc_uint<1> comp_r_11_76 = COMP < Q > (input[11].metric, input[76].metric);
    sc_uint<1> comp_r_11_77 = COMP < Q > (input[11].metric, input[77].metric);
    sc_uint<1> comp_r_11_78 = COMP < Q > (input[11].metric, input[78].metric);
    sc_uint<1> comp_r_11_79 = COMP < Q > (input[11].metric, input[79].metric);
    sc_uint<1> comp_r_11_80 = COMP < Q > (input[11].metric, input[80].metric);
    sc_uint<1> comp_r_11_81 = COMP < Q > (input[11].metric, input[81].metric);
    sc_uint<1> comp_r_11_82 = COMP < Q > (input[11].metric, input[82].metric);
    sc_uint<1> comp_r_11_83 = COMP < Q > (input[11].metric, input[83].metric);
    sc_uint<1> comp_r_11_84 = COMP < Q > (input[11].metric, input[84].metric);
    sc_uint<1> comp_r_11_85 = COMP < Q > (input[11].metric, input[85].metric);
    sc_uint<1> comp_r_11_86 = COMP < Q > (input[11].metric, input[86].metric);
    sc_uint<1> comp_r_11_87 = COMP < Q > (input[11].metric, input[87].metric);
    sc_uint<1> comp_r_11_88 = COMP < Q > (input[11].metric, input[88].metric);
    sc_uint<1> comp_r_11_89 = COMP < Q > (input[11].metric, input[89].metric);
    sc_uint<1> comp_r_11_90 = COMP < Q > (input[11].metric, input[90].metric);
    sc_uint<1> comp_r_11_91 = COMP < Q > (input[11].metric, input[91].metric);
    sc_uint<1> comp_r_11_92 = COMP < Q > (input[11].metric, input[92].metric);
    sc_uint<1> comp_r_11_93 = COMP < Q > (input[11].metric, input[93].metric);
    sc_uint<1> comp_r_11_94 = COMP < Q > (input[11].metric, input[94].metric);
    sc_uint<1> comp_r_11_95 = COMP < Q > (input[11].metric, input[95].metric);
    sc_uint<1> comp_r_11_96 = COMP < Q > (input[11].metric, input[96].metric);
    sc_uint<1> comp_r_11_97 = COMP < Q > (input[11].metric, input[97].metric);
    sc_uint<1> comp_r_11_98 = COMP < Q > (input[11].metric, input[98].metric);
    sc_uint<1> comp_r_11_99 = COMP < Q > (input[11].metric, input[99].metric);
    sc_uint<1> comp_r_11_100 = COMP < Q > (input[11].metric, input[100].metric);
    sc_uint<1> comp_r_11_101 = COMP < Q > (input[11].metric, input[101].metric);
    sc_uint<1> comp_r_11_102 = COMP < Q > (input[11].metric, input[102].metric);
    sc_uint<1> comp_r_11_103 = COMP < Q > (input[11].metric, input[103].metric);
    sc_uint<1> comp_r_11_104 = COMP < Q > (input[11].metric, input[104].metric);
    sc_uint<1> comp_r_11_105 = COMP < Q > (input[11].metric, input[105].metric);
    sc_uint<1> comp_r_11_106 = COMP < Q > (input[11].metric, input[106].metric);
    sc_uint<1> comp_r_11_107 = COMP < Q > (input[11].metric, input[107].metric);
    sc_uint<1> comp_r_11_108 = COMP < Q > (input[11].metric, input[108].metric);
    sc_uint<1> comp_r_11_109 = COMP < Q > (input[11].metric, input[109].metric);
    sc_uint<1> comp_r_11_110 = COMP < Q > (input[11].metric, input[110].metric);
    sc_uint<1> comp_r_11_111 = COMP < Q > (input[11].metric, input[111].metric);
    sc_uint<1> comp_r_11_112 = COMP < Q > (input[11].metric, input[112].metric);
    sc_uint<1> comp_r_11_113 = COMP < Q > (input[11].metric, input[113].metric);
    sc_uint<1> comp_r_11_114 = COMP < Q > (input[11].metric, input[114].metric);
    sc_uint<1> comp_r_11_115 = COMP < Q > (input[11].metric, input[115].metric);
    sc_uint<1> comp_r_11_116 = COMP < Q > (input[11].metric, input[116].metric);
    sc_uint<1> comp_r_11_117 = COMP < Q > (input[11].metric, input[117].metric);
    sc_uint<1> comp_r_11_118 = COMP < Q > (input[11].metric, input[118].metric);
    sc_uint<1> comp_r_11_119 = COMP < Q > (input[11].metric, input[119].metric);
    sc_uint<1> comp_r_11_120 = COMP < Q > (input[11].metric, input[120].metric);
    sc_uint<1> comp_r_11_121 = COMP < Q > (input[11].metric, input[121].metric);
    sc_uint<1> comp_r_11_122 = COMP < Q > (input[11].metric, input[122].metric);
    sc_uint<1> comp_r_11_123 = COMP < Q > (input[11].metric, input[123].metric);
    sc_uint<1> comp_r_11_124 = COMP < Q > (input[11].metric, input[124].metric);
    sc_uint<1> comp_r_11_125 = COMP < Q > (input[11].metric, input[125].metric);
    sc_uint<1> comp_r_11_126 = COMP < Q > (input[11].metric, input[126].metric);

    sc_uint<1> comp_r_13_14 = COMP < Q > (input[13].metric, input[14].metric);
    sc_uint<1> comp_r_13_15 = COMP < Q > (input[13].metric, input[15].metric);
    sc_uint<1> comp_r_13_16 = COMP < Q > (input[13].metric, input[16].metric);
    sc_uint<1> comp_r_13_17 = COMP < Q > (input[13].metric, input[17].metric);
    sc_uint<1> comp_r_13_18 = COMP < Q > (input[13].metric, input[18].metric);
    sc_uint<1> comp_r_13_19 = COMP < Q > (input[13].metric, input[19].metric);
    sc_uint<1> comp_r_13_20 = COMP < Q > (input[13].metric, input[20].metric);
    sc_uint<1> comp_r_13_21 = COMP < Q > (input[13].metric, input[21].metric);
    sc_uint<1> comp_r_13_22 = COMP < Q > (input[13].metric, input[22].metric);
    sc_uint<1> comp_r_13_23 = COMP < Q > (input[13].metric, input[23].metric);
    sc_uint<1> comp_r_13_24 = COMP < Q > (input[13].metric, input[24].metric);
    sc_uint<1> comp_r_13_25 = COMP < Q > (input[13].metric, input[25].metric);
    sc_uint<1> comp_r_13_26 = COMP < Q > (input[13].metric, input[26].metric);
    sc_uint<1> comp_r_13_27 = COMP < Q > (input[13].metric, input[27].metric);
    sc_uint<1> comp_r_13_28 = COMP < Q > (input[13].metric, input[28].metric);
    sc_uint<1> comp_r_13_29 = COMP < Q > (input[13].metric, input[29].metric);
    sc_uint<1> comp_r_13_30 = COMP < Q > (input[13].metric, input[30].metric);
    sc_uint<1> comp_r_13_31 = COMP < Q > (input[13].metric, input[31].metric);
    sc_uint<1> comp_r_13_32 = COMP < Q > (input[13].metric, input[32].metric);
    sc_uint<1> comp_r_13_33 = COMP < Q > (input[13].metric, input[33].metric);
    sc_uint<1> comp_r_13_34 = COMP < Q > (input[13].metric, input[34].metric);
    sc_uint<1> comp_r_13_35 = COMP < Q > (input[13].metric, input[35].metric);
    sc_uint<1> comp_r_13_36 = COMP < Q > (input[13].metric, input[36].metric);
    sc_uint<1> comp_r_13_37 = COMP < Q > (input[13].metric, input[37].metric);
    sc_uint<1> comp_r_13_38 = COMP < Q > (input[13].metric, input[38].metric);
    sc_uint<1> comp_r_13_39 = COMP < Q > (input[13].metric, input[39].metric);
    sc_uint<1> comp_r_13_40 = COMP < Q > (input[13].metric, input[40].metric);
    sc_uint<1> comp_r_13_41 = COMP < Q > (input[13].metric, input[41].metric);
    sc_uint<1> comp_r_13_42 = COMP < Q > (input[13].metric, input[42].metric);
    sc_uint<1> comp_r_13_43 = COMP < Q > (input[13].metric, input[43].metric);
    sc_uint<1> comp_r_13_44 = COMP < Q > (input[13].metric, input[44].metric);
    sc_uint<1> comp_r_13_45 = COMP < Q > (input[13].metric, input[45].metric);
    sc_uint<1> comp_r_13_46 = COMP < Q > (input[13].metric, input[46].metric);
    sc_uint<1> comp_r_13_47 = COMP < Q > (input[13].metric, input[47].metric);
    sc_uint<1> comp_r_13_48 = COMP < Q > (input[13].metric, input[48].metric);
    sc_uint<1> comp_r_13_49 = COMP < Q > (input[13].metric, input[49].metric);
    sc_uint<1> comp_r_13_50 = COMP < Q > (input[13].metric, input[50].metric);
    sc_uint<1> comp_r_13_51 = COMP < Q > (input[13].metric, input[51].metric);
    sc_uint<1> comp_r_13_52 = COMP < Q > (input[13].metric, input[52].metric);
    sc_uint<1> comp_r_13_53 = COMP < Q > (input[13].metric, input[53].metric);
    sc_uint<1> comp_r_13_54 = COMP < Q > (input[13].metric, input[54].metric);
    sc_uint<1> comp_r_13_55 = COMP < Q > (input[13].metric, input[55].metric);
    sc_uint<1> comp_r_13_56 = COMP < Q > (input[13].metric, input[56].metric);
    sc_uint<1> comp_r_13_57 = COMP < Q > (input[13].metric, input[57].metric);
    sc_uint<1> comp_r_13_58 = COMP < Q > (input[13].metric, input[58].metric);
    sc_uint<1> comp_r_13_59 = COMP < Q > (input[13].metric, input[59].metric);
    sc_uint<1> comp_r_13_60 = COMP < Q > (input[13].metric, input[60].metric);
    sc_uint<1> comp_r_13_61 = COMP < Q > (input[13].metric, input[61].metric);
    sc_uint<1> comp_r_13_62 = COMP < Q > (input[13].metric, input[62].metric);
    sc_uint<1> comp_r_13_63 = COMP < Q > (input[13].metric, input[63].metric);
    sc_uint<1> comp_r_13_64 = COMP < Q > (input[13].metric, input[64].metric);
    sc_uint<1> comp_r_13_65 = COMP < Q > (input[13].metric, input[65].metric);
    sc_uint<1> comp_r_13_66 = COMP < Q > (input[13].metric, input[66].metric);
    sc_uint<1> comp_r_13_67 = COMP < Q > (input[13].metric, input[67].metric);
    sc_uint<1> comp_r_13_68 = COMP < Q > (input[13].metric, input[68].metric);
    sc_uint<1> comp_r_13_69 = COMP < Q > (input[13].metric, input[69].metric);
    sc_uint<1> comp_r_13_70 = COMP < Q > (input[13].metric, input[70].metric);
    sc_uint<1> comp_r_13_71 = COMP < Q > (input[13].metric, input[71].metric);
    sc_uint<1> comp_r_13_72 = COMP < Q > (input[13].metric, input[72].metric);
    sc_uint<1> comp_r_13_73 = COMP < Q > (input[13].metric, input[73].metric);
    sc_uint<1> comp_r_13_74 = COMP < Q > (input[13].metric, input[74].metric);
    sc_uint<1> comp_r_13_75 = COMP < Q > (input[13].metric, input[75].metric);
    sc_uint<1> comp_r_13_76 = COMP < Q > (input[13].metric, input[76].metric);
    sc_uint<1> comp_r_13_77 = COMP < Q > (input[13].metric, input[77].metric);
    sc_uint<1> comp_r_13_78 = COMP < Q > (input[13].metric, input[78].metric);
    sc_uint<1> comp_r_13_79 = COMP < Q > (input[13].metric, input[79].metric);
    sc_uint<1> comp_r_13_80 = COMP < Q > (input[13].metric, input[80].metric);
    sc_uint<1> comp_r_13_81 = COMP < Q > (input[13].metric, input[81].metric);
    sc_uint<1> comp_r_13_82 = COMP < Q > (input[13].metric, input[82].metric);
    sc_uint<1> comp_r_13_83 = COMP < Q > (input[13].metric, input[83].metric);
    sc_uint<1> comp_r_13_84 = COMP < Q > (input[13].metric, input[84].metric);
    sc_uint<1> comp_r_13_85 = COMP < Q > (input[13].metric, input[85].metric);
    sc_uint<1> comp_r_13_86 = COMP < Q > (input[13].metric, input[86].metric);
    sc_uint<1> comp_r_13_87 = COMP < Q > (input[13].metric, input[87].metric);
    sc_uint<1> comp_r_13_88 = COMP < Q > (input[13].metric, input[88].metric);
    sc_uint<1> comp_r_13_89 = COMP < Q > (input[13].metric, input[89].metric);
    sc_uint<1> comp_r_13_90 = COMP < Q > (input[13].metric, input[90].metric);
    sc_uint<1> comp_r_13_91 = COMP < Q > (input[13].metric, input[91].metric);
    sc_uint<1> comp_r_13_92 = COMP < Q > (input[13].metric, input[92].metric);
    sc_uint<1> comp_r_13_93 = COMP < Q > (input[13].metric, input[93].metric);
    sc_uint<1> comp_r_13_94 = COMP < Q > (input[13].metric, input[94].metric);
    sc_uint<1> comp_r_13_95 = COMP < Q > (input[13].metric, input[95].metric);
    sc_uint<1> comp_r_13_96 = COMP < Q > (input[13].metric, input[96].metric);
    sc_uint<1> comp_r_13_97 = COMP < Q > (input[13].metric, input[97].metric);
    sc_uint<1> comp_r_13_98 = COMP < Q > (input[13].metric, input[98].metric);
    sc_uint<1> comp_r_13_99 = COMP < Q > (input[13].metric, input[99].metric);
    sc_uint<1> comp_r_13_100 = COMP < Q > (input[13].metric, input[100].metric);
    sc_uint<1> comp_r_13_101 = COMP < Q > (input[13].metric, input[101].metric);
    sc_uint<1> comp_r_13_102 = COMP < Q > (input[13].metric, input[102].metric);
    sc_uint<1> comp_r_13_103 = COMP < Q > (input[13].metric, input[103].metric);
    sc_uint<1> comp_r_13_104 = COMP < Q > (input[13].metric, input[104].metric);
    sc_uint<1> comp_r_13_105 = COMP < Q > (input[13].metric, input[105].metric);
    sc_uint<1> comp_r_13_106 = COMP < Q > (input[13].metric, input[106].metric);
    sc_uint<1> comp_r_13_107 = COMP < Q > (input[13].metric, input[107].metric);
    sc_uint<1> comp_r_13_108 = COMP < Q > (input[13].metric, input[108].metric);
    sc_uint<1> comp_r_13_109 = COMP < Q > (input[13].metric, input[109].metric);
    sc_uint<1> comp_r_13_110 = COMP < Q > (input[13].metric, input[110].metric);
    sc_uint<1> comp_r_13_111 = COMP < Q > (input[13].metric, input[111].metric);
    sc_uint<1> comp_r_13_112 = COMP < Q > (input[13].metric, input[112].metric);
    sc_uint<1> comp_r_13_113 = COMP < Q > (input[13].metric, input[113].metric);
    sc_uint<1> comp_r_13_114 = COMP < Q > (input[13].metric, input[114].metric);
    sc_uint<1> comp_r_13_115 = COMP < Q > (input[13].metric, input[115].metric);
    sc_uint<1> comp_r_13_116 = COMP < Q > (input[13].metric, input[116].metric);
    sc_uint<1> comp_r_13_117 = COMP < Q > (input[13].metric, input[117].metric);
    sc_uint<1> comp_r_13_118 = COMP < Q > (input[13].metric, input[118].metric);
    sc_uint<1> comp_r_13_119 = COMP < Q > (input[13].metric, input[119].metric);
    sc_uint<1> comp_r_13_120 = COMP < Q > (input[13].metric, input[120].metric);
    sc_uint<1> comp_r_13_121 = COMP < Q > (input[13].metric, input[121].metric);
    sc_uint<1> comp_r_13_122 = COMP < Q > (input[13].metric, input[122].metric);
    sc_uint<1> comp_r_13_123 = COMP < Q > (input[13].metric, input[123].metric);
    sc_uint<1> comp_r_13_124 = COMP < Q > (input[13].metric, input[124].metric);
    sc_uint<1> comp_r_13_125 = COMP < Q > (input[13].metric, input[125].metric);
    sc_uint<1> comp_r_13_126 = COMP < Q > (input[13].metric, input[126].metric);

    sc_uint<1> comp_r_15_16 = COMP < Q > (input[15].metric, input[16].metric);
    sc_uint<1> comp_r_15_17 = COMP < Q > (input[15].metric, input[17].metric);
    sc_uint<1> comp_r_15_18 = COMP < Q > (input[15].metric, input[18].metric);
    sc_uint<1> comp_r_15_19 = COMP < Q > (input[15].metric, input[19].metric);
    sc_uint<1> comp_r_15_20 = COMP < Q > (input[15].metric, input[20].metric);
    sc_uint<1> comp_r_15_21 = COMP < Q > (input[15].metric, input[21].metric);
    sc_uint<1> comp_r_15_22 = COMP < Q > (input[15].metric, input[22].metric);
    sc_uint<1> comp_r_15_23 = COMP < Q > (input[15].metric, input[23].metric);
    sc_uint<1> comp_r_15_24 = COMP < Q > (input[15].metric, input[24].metric);
    sc_uint<1> comp_r_15_25 = COMP < Q > (input[15].metric, input[25].metric);
    sc_uint<1> comp_r_15_26 = COMP < Q > (input[15].metric, input[26].metric);
    sc_uint<1> comp_r_15_27 = COMP < Q > (input[15].metric, input[27].metric);
    sc_uint<1> comp_r_15_28 = COMP < Q > (input[15].metric, input[28].metric);
    sc_uint<1> comp_r_15_29 = COMP < Q > (input[15].metric, input[29].metric);
    sc_uint<1> comp_r_15_30 = COMP < Q > (input[15].metric, input[30].metric);
    sc_uint<1> comp_r_15_31 = COMP < Q > (input[15].metric, input[31].metric);
    sc_uint<1> comp_r_15_32 = COMP < Q > (input[15].metric, input[32].metric);
    sc_uint<1> comp_r_15_33 = COMP < Q > (input[15].metric, input[33].metric);
    sc_uint<1> comp_r_15_34 = COMP < Q > (input[15].metric, input[34].metric);
    sc_uint<1> comp_r_15_35 = COMP < Q > (input[15].metric, input[35].metric);
    sc_uint<1> comp_r_15_36 = COMP < Q > (input[15].metric, input[36].metric);
    sc_uint<1> comp_r_15_37 = COMP < Q > (input[15].metric, input[37].metric);
    sc_uint<1> comp_r_15_38 = COMP < Q > (input[15].metric, input[38].metric);
    sc_uint<1> comp_r_15_39 = COMP < Q > (input[15].metric, input[39].metric);
    sc_uint<1> comp_r_15_40 = COMP < Q > (input[15].metric, input[40].metric);
    sc_uint<1> comp_r_15_41 = COMP < Q > (input[15].metric, input[41].metric);
    sc_uint<1> comp_r_15_42 = COMP < Q > (input[15].metric, input[42].metric);
    sc_uint<1> comp_r_15_43 = COMP < Q > (input[15].metric, input[43].metric);
    sc_uint<1> comp_r_15_44 = COMP < Q > (input[15].metric, input[44].metric);
    sc_uint<1> comp_r_15_45 = COMP < Q > (input[15].metric, input[45].metric);
    sc_uint<1> comp_r_15_46 = COMP < Q > (input[15].metric, input[46].metric);
    sc_uint<1> comp_r_15_47 = COMP < Q > (input[15].metric, input[47].metric);
    sc_uint<1> comp_r_15_48 = COMP < Q > (input[15].metric, input[48].metric);
    sc_uint<1> comp_r_15_49 = COMP < Q > (input[15].metric, input[49].metric);
    sc_uint<1> comp_r_15_50 = COMP < Q > (input[15].metric, input[50].metric);
    sc_uint<1> comp_r_15_51 = COMP < Q > (input[15].metric, input[51].metric);
    sc_uint<1> comp_r_15_52 = COMP < Q > (input[15].metric, input[52].metric);
    sc_uint<1> comp_r_15_53 = COMP < Q > (input[15].metric, input[53].metric);
    sc_uint<1> comp_r_15_54 = COMP < Q > (input[15].metric, input[54].metric);
    sc_uint<1> comp_r_15_55 = COMP < Q > (input[15].metric, input[55].metric);
    sc_uint<1> comp_r_15_56 = COMP < Q > (input[15].metric, input[56].metric);
    sc_uint<1> comp_r_15_57 = COMP < Q > (input[15].metric, input[57].metric);
    sc_uint<1> comp_r_15_58 = COMP < Q > (input[15].metric, input[58].metric);
    sc_uint<1> comp_r_15_59 = COMP < Q > (input[15].metric, input[59].metric);
    sc_uint<1> comp_r_15_60 = COMP < Q > (input[15].metric, input[60].metric);
    sc_uint<1> comp_r_15_61 = COMP < Q > (input[15].metric, input[61].metric);
    sc_uint<1> comp_r_15_62 = COMP < Q > (input[15].metric, input[62].metric);
    sc_uint<1> comp_r_15_63 = COMP < Q > (input[15].metric, input[63].metric);
    sc_uint<1> comp_r_15_64 = COMP < Q > (input[15].metric, input[64].metric);
    sc_uint<1> comp_r_15_65 = COMP < Q > (input[15].metric, input[65].metric);
    sc_uint<1> comp_r_15_66 = COMP < Q > (input[15].metric, input[66].metric);
    sc_uint<1> comp_r_15_67 = COMP < Q > (input[15].metric, input[67].metric);
    sc_uint<1> comp_r_15_68 = COMP < Q > (input[15].metric, input[68].metric);
    sc_uint<1> comp_r_15_69 = COMP < Q > (input[15].metric, input[69].metric);
    sc_uint<1> comp_r_15_70 = COMP < Q > (input[15].metric, input[70].metric);
    sc_uint<1> comp_r_15_71 = COMP < Q > (input[15].metric, input[71].metric);
    sc_uint<1> comp_r_15_72 = COMP < Q > (input[15].metric, input[72].metric);
    sc_uint<1> comp_r_15_73 = COMP < Q > (input[15].metric, input[73].metric);
    sc_uint<1> comp_r_15_74 = COMP < Q > (input[15].metric, input[74].metric);
    sc_uint<1> comp_r_15_75 = COMP < Q > (input[15].metric, input[75].metric);
    sc_uint<1> comp_r_15_76 = COMP < Q > (input[15].metric, input[76].metric);
    sc_uint<1> comp_r_15_77 = COMP < Q > (input[15].metric, input[77].metric);
    sc_uint<1> comp_r_15_78 = COMP < Q > (input[15].metric, input[78].metric);
    sc_uint<1> comp_r_15_79 = COMP < Q > (input[15].metric, input[79].metric);
    sc_uint<1> comp_r_15_80 = COMP < Q > (input[15].metric, input[80].metric);
    sc_uint<1> comp_r_15_81 = COMP < Q > (input[15].metric, input[81].metric);
    sc_uint<1> comp_r_15_82 = COMP < Q > (input[15].metric, input[82].metric);
    sc_uint<1> comp_r_15_83 = COMP < Q > (input[15].metric, input[83].metric);
    sc_uint<1> comp_r_15_84 = COMP < Q > (input[15].metric, input[84].metric);
    sc_uint<1> comp_r_15_85 = COMP < Q > (input[15].metric, input[85].metric);
    sc_uint<1> comp_r_15_86 = COMP < Q > (input[15].metric, input[86].metric);
    sc_uint<1> comp_r_15_87 = COMP < Q > (input[15].metric, input[87].metric);
    sc_uint<1> comp_r_15_88 = COMP < Q > (input[15].metric, input[88].metric);
    sc_uint<1> comp_r_15_89 = COMP < Q > (input[15].metric, input[89].metric);
    sc_uint<1> comp_r_15_90 = COMP < Q > (input[15].metric, input[90].metric);
    sc_uint<1> comp_r_15_91 = COMP < Q > (input[15].metric, input[91].metric);
    sc_uint<1> comp_r_15_92 = COMP < Q > (input[15].metric, input[92].metric);
    sc_uint<1> comp_r_15_93 = COMP < Q > (input[15].metric, input[93].metric);
    sc_uint<1> comp_r_15_94 = COMP < Q > (input[15].metric, input[94].metric);
    sc_uint<1> comp_r_15_95 = COMP < Q > (input[15].metric, input[95].metric);
    sc_uint<1> comp_r_15_96 = COMP < Q > (input[15].metric, input[96].metric);
    sc_uint<1> comp_r_15_97 = COMP < Q > (input[15].metric, input[97].metric);
    sc_uint<1> comp_r_15_98 = COMP < Q > (input[15].metric, input[98].metric);
    sc_uint<1> comp_r_15_99 = COMP < Q > (input[15].metric, input[99].metric);
    sc_uint<1> comp_r_15_100 = COMP < Q > (input[15].metric, input[100].metric);
    sc_uint<1> comp_r_15_101 = COMP < Q > (input[15].metric, input[101].metric);
    sc_uint<1> comp_r_15_102 = COMP < Q > (input[15].metric, input[102].metric);
    sc_uint<1> comp_r_15_103 = COMP < Q > (input[15].metric, input[103].metric);
    sc_uint<1> comp_r_15_104 = COMP < Q > (input[15].metric, input[104].metric);
    sc_uint<1> comp_r_15_105 = COMP < Q > (input[15].metric, input[105].metric);
    sc_uint<1> comp_r_15_106 = COMP < Q > (input[15].metric, input[106].metric);
    sc_uint<1> comp_r_15_107 = COMP < Q > (input[15].metric, input[107].metric);
    sc_uint<1> comp_r_15_108 = COMP < Q > (input[15].metric, input[108].metric);
    sc_uint<1> comp_r_15_109 = COMP < Q > (input[15].metric, input[109].metric);
    sc_uint<1> comp_r_15_110 = COMP < Q > (input[15].metric, input[110].metric);
    sc_uint<1> comp_r_15_111 = COMP < Q > (input[15].metric, input[111].metric);
    sc_uint<1> comp_r_15_112 = COMP < Q > (input[15].metric, input[112].metric);
    sc_uint<1> comp_r_15_113 = COMP < Q > (input[15].metric, input[113].metric);
    sc_uint<1> comp_r_15_114 = COMP < Q > (input[15].metric, input[114].metric);
    sc_uint<1> comp_r_15_115 = COMP < Q > (input[15].metric, input[115].metric);
    sc_uint<1> comp_r_15_116 = COMP < Q > (input[15].metric, input[116].metric);
    sc_uint<1> comp_r_15_117 = COMP < Q > (input[15].metric, input[117].metric);
    sc_uint<1> comp_r_15_118 = COMP < Q > (input[15].metric, input[118].metric);
    sc_uint<1> comp_r_15_119 = COMP < Q > (input[15].metric, input[119].metric);
    sc_uint<1> comp_r_15_120 = COMP < Q > (input[15].metric, input[120].metric);
    sc_uint<1> comp_r_15_121 = COMP < Q > (input[15].metric, input[121].metric);
    sc_uint<1> comp_r_15_122 = COMP < Q > (input[15].metric, input[122].metric);
    sc_uint<1> comp_r_15_123 = COMP < Q > (input[15].metric, input[123].metric);
    sc_uint<1> comp_r_15_124 = COMP < Q > (input[15].metric, input[124].metric);
    sc_uint<1> comp_r_15_125 = COMP < Q > (input[15].metric, input[125].metric);
    sc_uint<1> comp_r_15_126 = COMP < Q > (input[15].metric, input[126].metric);

    sc_uint<1> comp_r_17_18 = COMP < Q > (input[17].metric, input[18].metric);
    sc_uint<1> comp_r_17_19 = COMP < Q > (input[17].metric, input[19].metric);
    sc_uint<1> comp_r_17_20 = COMP < Q > (input[17].metric, input[20].metric);
    sc_uint<1> comp_r_17_21 = COMP < Q > (input[17].metric, input[21].metric);
    sc_uint<1> comp_r_17_22 = COMP < Q > (input[17].metric, input[22].metric);
    sc_uint<1> comp_r_17_23 = COMP < Q > (input[17].metric, input[23].metric);
    sc_uint<1> comp_r_17_24 = COMP < Q > (input[17].metric, input[24].metric);
    sc_uint<1> comp_r_17_25 = COMP < Q > (input[17].metric, input[25].metric);
    sc_uint<1> comp_r_17_26 = COMP < Q > (input[17].metric, input[26].metric);
    sc_uint<1> comp_r_17_27 = COMP < Q > (input[17].metric, input[27].metric);
    sc_uint<1> comp_r_17_28 = COMP < Q > (input[17].metric, input[28].metric);
    sc_uint<1> comp_r_17_29 = COMP < Q > (input[17].metric, input[29].metric);
    sc_uint<1> comp_r_17_30 = COMP < Q > (input[17].metric, input[30].metric);
    sc_uint<1> comp_r_17_31 = COMP < Q > (input[17].metric, input[31].metric);
    sc_uint<1> comp_r_17_32 = COMP < Q > (input[17].metric, input[32].metric);
    sc_uint<1> comp_r_17_33 = COMP < Q > (input[17].metric, input[33].metric);
    sc_uint<1> comp_r_17_34 = COMP < Q > (input[17].metric, input[34].metric);
    sc_uint<1> comp_r_17_35 = COMP < Q > (input[17].metric, input[35].metric);
    sc_uint<1> comp_r_17_36 = COMP < Q > (input[17].metric, input[36].metric);
    sc_uint<1> comp_r_17_37 = COMP < Q > (input[17].metric, input[37].metric);
    sc_uint<1> comp_r_17_38 = COMP < Q > (input[17].metric, input[38].metric);
    sc_uint<1> comp_r_17_39 = COMP < Q > (input[17].metric, input[39].metric);
    sc_uint<1> comp_r_17_40 = COMP < Q > (input[17].metric, input[40].metric);
    sc_uint<1> comp_r_17_41 = COMP < Q > (input[17].metric, input[41].metric);
    sc_uint<1> comp_r_17_42 = COMP < Q > (input[17].metric, input[42].metric);
    sc_uint<1> comp_r_17_43 = COMP < Q > (input[17].metric, input[43].metric);
    sc_uint<1> comp_r_17_44 = COMP < Q > (input[17].metric, input[44].metric);
    sc_uint<1> comp_r_17_45 = COMP < Q > (input[17].metric, input[45].metric);
    sc_uint<1> comp_r_17_46 = COMP < Q > (input[17].metric, input[46].metric);
    sc_uint<1> comp_r_17_47 = COMP < Q > (input[17].metric, input[47].metric);
    sc_uint<1> comp_r_17_48 = COMP < Q > (input[17].metric, input[48].metric);
    sc_uint<1> comp_r_17_49 = COMP < Q > (input[17].metric, input[49].metric);
    sc_uint<1> comp_r_17_50 = COMP < Q > (input[17].metric, input[50].metric);
    sc_uint<1> comp_r_17_51 = COMP < Q > (input[17].metric, input[51].metric);
    sc_uint<1> comp_r_17_52 = COMP < Q > (input[17].metric, input[52].metric);
    sc_uint<1> comp_r_17_53 = COMP < Q > (input[17].metric, input[53].metric);
    sc_uint<1> comp_r_17_54 = COMP < Q > (input[17].metric, input[54].metric);
    sc_uint<1> comp_r_17_55 = COMP < Q > (input[17].metric, input[55].metric);
    sc_uint<1> comp_r_17_56 = COMP < Q > (input[17].metric, input[56].metric);
    sc_uint<1> comp_r_17_57 = COMP < Q > (input[17].metric, input[57].metric);
    sc_uint<1> comp_r_17_58 = COMP < Q > (input[17].metric, input[58].metric);
    sc_uint<1> comp_r_17_59 = COMP < Q > (input[17].metric, input[59].metric);
    sc_uint<1> comp_r_17_60 = COMP < Q > (input[17].metric, input[60].metric);
    sc_uint<1> comp_r_17_61 = COMP < Q > (input[17].metric, input[61].metric);
    sc_uint<1> comp_r_17_62 = COMP < Q > (input[17].metric, input[62].metric);
    sc_uint<1> comp_r_17_63 = COMP < Q > (input[17].metric, input[63].metric);
    sc_uint<1> comp_r_17_64 = COMP < Q > (input[17].metric, input[64].metric);
    sc_uint<1> comp_r_17_65 = COMP < Q > (input[17].metric, input[65].metric);
    sc_uint<1> comp_r_17_66 = COMP < Q > (input[17].metric, input[66].metric);
    sc_uint<1> comp_r_17_67 = COMP < Q > (input[17].metric, input[67].metric);
    sc_uint<1> comp_r_17_68 = COMP < Q > (input[17].metric, input[68].metric);
    sc_uint<1> comp_r_17_69 = COMP < Q > (input[17].metric, input[69].metric);
    sc_uint<1> comp_r_17_70 = COMP < Q > (input[17].metric, input[70].metric);
    sc_uint<1> comp_r_17_71 = COMP < Q > (input[17].metric, input[71].metric);
    sc_uint<1> comp_r_17_72 = COMP < Q > (input[17].metric, input[72].metric);
    sc_uint<1> comp_r_17_73 = COMP < Q > (input[17].metric, input[73].metric);
    sc_uint<1> comp_r_17_74 = COMP < Q > (input[17].metric, input[74].metric);
    sc_uint<1> comp_r_17_75 = COMP < Q > (input[17].metric, input[75].metric);
    sc_uint<1> comp_r_17_76 = COMP < Q > (input[17].metric, input[76].metric);
    sc_uint<1> comp_r_17_77 = COMP < Q > (input[17].metric, input[77].metric);
    sc_uint<1> comp_r_17_78 = COMP < Q > (input[17].metric, input[78].metric);
    sc_uint<1> comp_r_17_79 = COMP < Q > (input[17].metric, input[79].metric);
    sc_uint<1> comp_r_17_80 = COMP < Q > (input[17].metric, input[80].metric);
    sc_uint<1> comp_r_17_81 = COMP < Q > (input[17].metric, input[81].metric);
    sc_uint<1> comp_r_17_82 = COMP < Q > (input[17].metric, input[82].metric);
    sc_uint<1> comp_r_17_83 = COMP < Q > (input[17].metric, input[83].metric);
    sc_uint<1> comp_r_17_84 = COMP < Q > (input[17].metric, input[84].metric);
    sc_uint<1> comp_r_17_85 = COMP < Q > (input[17].metric, input[85].metric);
    sc_uint<1> comp_r_17_86 = COMP < Q > (input[17].metric, input[86].metric);
    sc_uint<1> comp_r_17_87 = COMP < Q > (input[17].metric, input[87].metric);
    sc_uint<1> comp_r_17_88 = COMP < Q > (input[17].metric, input[88].metric);
    sc_uint<1> comp_r_17_89 = COMP < Q > (input[17].metric, input[89].metric);
    sc_uint<1> comp_r_17_90 = COMP < Q > (input[17].metric, input[90].metric);
    sc_uint<1> comp_r_17_91 = COMP < Q > (input[17].metric, input[91].metric);
    sc_uint<1> comp_r_17_92 = COMP < Q > (input[17].metric, input[92].metric);
    sc_uint<1> comp_r_17_93 = COMP < Q > (input[17].metric, input[93].metric);
    sc_uint<1> comp_r_17_94 = COMP < Q > (input[17].metric, input[94].metric);
    sc_uint<1> comp_r_17_95 = COMP < Q > (input[17].metric, input[95].metric);
    sc_uint<1> comp_r_17_96 = COMP < Q > (input[17].metric, input[96].metric);
    sc_uint<1> comp_r_17_97 = COMP < Q > (input[17].metric, input[97].metric);
    sc_uint<1> comp_r_17_98 = COMP < Q > (input[17].metric, input[98].metric);
    sc_uint<1> comp_r_17_99 = COMP < Q > (input[17].metric, input[99].metric);
    sc_uint<1> comp_r_17_100 = COMP < Q > (input[17].metric, input[100].metric);
    sc_uint<1> comp_r_17_101 = COMP < Q > (input[17].metric, input[101].metric);
    sc_uint<1> comp_r_17_102 = COMP < Q > (input[17].metric, input[102].metric);
    sc_uint<1> comp_r_17_103 = COMP < Q > (input[17].metric, input[103].metric);
    sc_uint<1> comp_r_17_104 = COMP < Q > (input[17].metric, input[104].metric);
    sc_uint<1> comp_r_17_105 = COMP < Q > (input[17].metric, input[105].metric);
    sc_uint<1> comp_r_17_106 = COMP < Q > (input[17].metric, input[106].metric);
    sc_uint<1> comp_r_17_107 = COMP < Q > (input[17].metric, input[107].metric);
    sc_uint<1> comp_r_17_108 = COMP < Q > (input[17].metric, input[108].metric);
    sc_uint<1> comp_r_17_109 = COMP < Q > (input[17].metric, input[109].metric);
    sc_uint<1> comp_r_17_110 = COMP < Q > (input[17].metric, input[110].metric);
    sc_uint<1> comp_r_17_111 = COMP < Q > (input[17].metric, input[111].metric);
    sc_uint<1> comp_r_17_112 = COMP < Q > (input[17].metric, input[112].metric);
    sc_uint<1> comp_r_17_113 = COMP < Q > (input[17].metric, input[113].metric);
    sc_uint<1> comp_r_17_114 = COMP < Q > (input[17].metric, input[114].metric);
    sc_uint<1> comp_r_17_115 = COMP < Q > (input[17].metric, input[115].metric);
    sc_uint<1> comp_r_17_116 = COMP < Q > (input[17].metric, input[116].metric);
    sc_uint<1> comp_r_17_117 = COMP < Q > (input[17].metric, input[117].metric);
    sc_uint<1> comp_r_17_118 = COMP < Q > (input[17].metric, input[118].metric);
    sc_uint<1> comp_r_17_119 = COMP < Q > (input[17].metric, input[119].metric);
    sc_uint<1> comp_r_17_120 = COMP < Q > (input[17].metric, input[120].metric);
    sc_uint<1> comp_r_17_121 = COMP < Q > (input[17].metric, input[121].metric);
    sc_uint<1> comp_r_17_122 = COMP < Q > (input[17].metric, input[122].metric);
    sc_uint<1> comp_r_17_123 = COMP < Q > (input[17].metric, input[123].metric);
    sc_uint<1> comp_r_17_124 = COMP < Q > (input[17].metric, input[124].metric);
    sc_uint<1> comp_r_17_125 = COMP < Q > (input[17].metric, input[125].metric);
    sc_uint<1> comp_r_17_126 = COMP < Q > (input[17].metric, input[126].metric);

    sc_uint<1> comp_r_19_20 = COMP < Q > (input[19].metric, input[20].metric);
    sc_uint<1> comp_r_19_21 = COMP < Q > (input[19].metric, input[21].metric);
    sc_uint<1> comp_r_19_22 = COMP < Q > (input[19].metric, input[22].metric);
    sc_uint<1> comp_r_19_23 = COMP < Q > (input[19].metric, input[23].metric);
    sc_uint<1> comp_r_19_24 = COMP < Q > (input[19].metric, input[24].metric);
    sc_uint<1> comp_r_19_25 = COMP < Q > (input[19].metric, input[25].metric);
    sc_uint<1> comp_r_19_26 = COMP < Q > (input[19].metric, input[26].metric);
    sc_uint<1> comp_r_19_27 = COMP < Q > (input[19].metric, input[27].metric);
    sc_uint<1> comp_r_19_28 = COMP < Q > (input[19].metric, input[28].metric);
    sc_uint<1> comp_r_19_29 = COMP < Q > (input[19].metric, input[29].metric);
    sc_uint<1> comp_r_19_30 = COMP < Q > (input[19].metric, input[30].metric);
    sc_uint<1> comp_r_19_31 = COMP < Q > (input[19].metric, input[31].metric);
    sc_uint<1> comp_r_19_32 = COMP < Q > (input[19].metric, input[32].metric);
    sc_uint<1> comp_r_19_33 = COMP < Q > (input[19].metric, input[33].metric);
    sc_uint<1> comp_r_19_34 = COMP < Q > (input[19].metric, input[34].metric);
    sc_uint<1> comp_r_19_35 = COMP < Q > (input[19].metric, input[35].metric);
    sc_uint<1> comp_r_19_36 = COMP < Q > (input[19].metric, input[36].metric);
    sc_uint<1> comp_r_19_37 = COMP < Q > (input[19].metric, input[37].metric);
    sc_uint<1> comp_r_19_38 = COMP < Q > (input[19].metric, input[38].metric);
    sc_uint<1> comp_r_19_39 = COMP < Q > (input[19].metric, input[39].metric);
    sc_uint<1> comp_r_19_40 = COMP < Q > (input[19].metric, input[40].metric);
    sc_uint<1> comp_r_19_41 = COMP < Q > (input[19].metric, input[41].metric);
    sc_uint<1> comp_r_19_42 = COMP < Q > (input[19].metric, input[42].metric);
    sc_uint<1> comp_r_19_43 = COMP < Q > (input[19].metric, input[43].metric);
    sc_uint<1> comp_r_19_44 = COMP < Q > (input[19].metric, input[44].metric);
    sc_uint<1> comp_r_19_45 = COMP < Q > (input[19].metric, input[45].metric);
    sc_uint<1> comp_r_19_46 = COMP < Q > (input[19].metric, input[46].metric);
    sc_uint<1> comp_r_19_47 = COMP < Q > (input[19].metric, input[47].metric);
    sc_uint<1> comp_r_19_48 = COMP < Q > (input[19].metric, input[48].metric);
    sc_uint<1> comp_r_19_49 = COMP < Q > (input[19].metric, input[49].metric);
    sc_uint<1> comp_r_19_50 = COMP < Q > (input[19].metric, input[50].metric);
    sc_uint<1> comp_r_19_51 = COMP < Q > (input[19].metric, input[51].metric);
    sc_uint<1> comp_r_19_52 = COMP < Q > (input[19].metric, input[52].metric);
    sc_uint<1> comp_r_19_53 = COMP < Q > (input[19].metric, input[53].metric);
    sc_uint<1> comp_r_19_54 = COMP < Q > (input[19].metric, input[54].metric);
    sc_uint<1> comp_r_19_55 = COMP < Q > (input[19].metric, input[55].metric);
    sc_uint<1> comp_r_19_56 = COMP < Q > (input[19].metric, input[56].metric);
    sc_uint<1> comp_r_19_57 = COMP < Q > (input[19].metric, input[57].metric);
    sc_uint<1> comp_r_19_58 = COMP < Q > (input[19].metric, input[58].metric);
    sc_uint<1> comp_r_19_59 = COMP < Q > (input[19].metric, input[59].metric);
    sc_uint<1> comp_r_19_60 = COMP < Q > (input[19].metric, input[60].metric);
    sc_uint<1> comp_r_19_61 = COMP < Q > (input[19].metric, input[61].metric);
    sc_uint<1> comp_r_19_62 = COMP < Q > (input[19].metric, input[62].metric);
    sc_uint<1> comp_r_19_63 = COMP < Q > (input[19].metric, input[63].metric);
    sc_uint<1> comp_r_19_64 = COMP < Q > (input[19].metric, input[64].metric);
    sc_uint<1> comp_r_19_65 = COMP < Q > (input[19].metric, input[65].metric);
    sc_uint<1> comp_r_19_66 = COMP < Q > (input[19].metric, input[66].metric);
    sc_uint<1> comp_r_19_67 = COMP < Q > (input[19].metric, input[67].metric);
    sc_uint<1> comp_r_19_68 = COMP < Q > (input[19].metric, input[68].metric);
    sc_uint<1> comp_r_19_69 = COMP < Q > (input[19].metric, input[69].metric);
    sc_uint<1> comp_r_19_70 = COMP < Q > (input[19].metric, input[70].metric);
    sc_uint<1> comp_r_19_71 = COMP < Q > (input[19].metric, input[71].metric);
    sc_uint<1> comp_r_19_72 = COMP < Q > (input[19].metric, input[72].metric);
    sc_uint<1> comp_r_19_73 = COMP < Q > (input[19].metric, input[73].metric);
    sc_uint<1> comp_r_19_74 = COMP < Q > (input[19].metric, input[74].metric);
    sc_uint<1> comp_r_19_75 = COMP < Q > (input[19].metric, input[75].metric);
    sc_uint<1> comp_r_19_76 = COMP < Q > (input[19].metric, input[76].metric);
    sc_uint<1> comp_r_19_77 = COMP < Q > (input[19].metric, input[77].metric);
    sc_uint<1> comp_r_19_78 = COMP < Q > (input[19].metric, input[78].metric);
    sc_uint<1> comp_r_19_79 = COMP < Q > (input[19].metric, input[79].metric);
    sc_uint<1> comp_r_19_80 = COMP < Q > (input[19].metric, input[80].metric);
    sc_uint<1> comp_r_19_81 = COMP < Q > (input[19].metric, input[81].metric);
    sc_uint<1> comp_r_19_82 = COMP < Q > (input[19].metric, input[82].metric);
    sc_uint<1> comp_r_19_83 = COMP < Q > (input[19].metric, input[83].metric);
    sc_uint<1> comp_r_19_84 = COMP < Q > (input[19].metric, input[84].metric);
    sc_uint<1> comp_r_19_85 = COMP < Q > (input[19].metric, input[85].metric);
    sc_uint<1> comp_r_19_86 = COMP < Q > (input[19].metric, input[86].metric);
    sc_uint<1> comp_r_19_87 = COMP < Q > (input[19].metric, input[87].metric);
    sc_uint<1> comp_r_19_88 = COMP < Q > (input[19].metric, input[88].metric);
    sc_uint<1> comp_r_19_89 = COMP < Q > (input[19].metric, input[89].metric);
    sc_uint<1> comp_r_19_90 = COMP < Q > (input[19].metric, input[90].metric);
    sc_uint<1> comp_r_19_91 = COMP < Q > (input[19].metric, input[91].metric);
    sc_uint<1> comp_r_19_92 = COMP < Q > (input[19].metric, input[92].metric);
    sc_uint<1> comp_r_19_93 = COMP < Q > (input[19].metric, input[93].metric);
    sc_uint<1> comp_r_19_94 = COMP < Q > (input[19].metric, input[94].metric);
    sc_uint<1> comp_r_19_95 = COMP < Q > (input[19].metric, input[95].metric);
    sc_uint<1> comp_r_19_96 = COMP < Q > (input[19].metric, input[96].metric);
    sc_uint<1> comp_r_19_97 = COMP < Q > (input[19].metric, input[97].metric);
    sc_uint<1> comp_r_19_98 = COMP < Q > (input[19].metric, input[98].metric);
    sc_uint<1> comp_r_19_99 = COMP < Q > (input[19].metric, input[99].metric);
    sc_uint<1> comp_r_19_100 = COMP < Q > (input[19].metric, input[100].metric);
    sc_uint<1> comp_r_19_101 = COMP < Q > (input[19].metric, input[101].metric);
    sc_uint<1> comp_r_19_102 = COMP < Q > (input[19].metric, input[102].metric);
    sc_uint<1> comp_r_19_103 = COMP < Q > (input[19].metric, input[103].metric);
    sc_uint<1> comp_r_19_104 = COMP < Q > (input[19].metric, input[104].metric);
    sc_uint<1> comp_r_19_105 = COMP < Q > (input[19].metric, input[105].metric);
    sc_uint<1> comp_r_19_106 = COMP < Q > (input[19].metric, input[106].metric);
    sc_uint<1> comp_r_19_107 = COMP < Q > (input[19].metric, input[107].metric);
    sc_uint<1> comp_r_19_108 = COMP < Q > (input[19].metric, input[108].metric);
    sc_uint<1> comp_r_19_109 = COMP < Q > (input[19].metric, input[109].metric);
    sc_uint<1> comp_r_19_110 = COMP < Q > (input[19].metric, input[110].metric);
    sc_uint<1> comp_r_19_111 = COMP < Q > (input[19].metric, input[111].metric);
    sc_uint<1> comp_r_19_112 = COMP < Q > (input[19].metric, input[112].metric);
    sc_uint<1> comp_r_19_113 = COMP < Q > (input[19].metric, input[113].metric);
    sc_uint<1> comp_r_19_114 = COMP < Q > (input[19].metric, input[114].metric);
    sc_uint<1> comp_r_19_115 = COMP < Q > (input[19].metric, input[115].metric);
    sc_uint<1> comp_r_19_116 = COMP < Q > (input[19].metric, input[116].metric);
    sc_uint<1> comp_r_19_117 = COMP < Q > (input[19].metric, input[117].metric);
    sc_uint<1> comp_r_19_118 = COMP < Q > (input[19].metric, input[118].metric);
    sc_uint<1> comp_r_19_119 = COMP < Q > (input[19].metric, input[119].metric);
    sc_uint<1> comp_r_19_120 = COMP < Q > (input[19].metric, input[120].metric);
    sc_uint<1> comp_r_19_121 = COMP < Q > (input[19].metric, input[121].metric);
    sc_uint<1> comp_r_19_122 = COMP < Q > (input[19].metric, input[122].metric);
    sc_uint<1> comp_r_19_123 = COMP < Q > (input[19].metric, input[123].metric);
    sc_uint<1> comp_r_19_124 = COMP < Q > (input[19].metric, input[124].metric);
    sc_uint<1> comp_r_19_125 = COMP < Q > (input[19].metric, input[125].metric);
    sc_uint<1> comp_r_19_126 = COMP < Q > (input[19].metric, input[126].metric);

    sc_uint<1> comp_r_21_22 = COMP < Q > (input[21].metric, input[22].metric);
    sc_uint<1> comp_r_21_23 = COMP < Q > (input[21].metric, input[23].metric);
    sc_uint<1> comp_r_21_24 = COMP < Q > (input[21].metric, input[24].metric);
    sc_uint<1> comp_r_21_25 = COMP < Q > (input[21].metric, input[25].metric);
    sc_uint<1> comp_r_21_26 = COMP < Q > (input[21].metric, input[26].metric);
    sc_uint<1> comp_r_21_27 = COMP < Q > (input[21].metric, input[27].metric);
    sc_uint<1> comp_r_21_28 = COMP < Q > (input[21].metric, input[28].metric);
    sc_uint<1> comp_r_21_29 = COMP < Q > (input[21].metric, input[29].metric);
    sc_uint<1> comp_r_21_30 = COMP < Q > (input[21].metric, input[30].metric);
    sc_uint<1> comp_r_21_31 = COMP < Q > (input[21].metric, input[31].metric);
    sc_uint<1> comp_r_21_32 = COMP < Q > (input[21].metric, input[32].metric);
    sc_uint<1> comp_r_21_33 = COMP < Q > (input[21].metric, input[33].metric);
    sc_uint<1> comp_r_21_34 = COMP < Q > (input[21].metric, input[34].metric);
    sc_uint<1> comp_r_21_35 = COMP < Q > (input[21].metric, input[35].metric);
    sc_uint<1> comp_r_21_36 = COMP < Q > (input[21].metric, input[36].metric);
    sc_uint<1> comp_r_21_37 = COMP < Q > (input[21].metric, input[37].metric);
    sc_uint<1> comp_r_21_38 = COMP < Q > (input[21].metric, input[38].metric);
    sc_uint<1> comp_r_21_39 = COMP < Q > (input[21].metric, input[39].metric);
    sc_uint<1> comp_r_21_40 = COMP < Q > (input[21].metric, input[40].metric);
    sc_uint<1> comp_r_21_41 = COMP < Q > (input[21].metric, input[41].metric);
    sc_uint<1> comp_r_21_42 = COMP < Q > (input[21].metric, input[42].metric);
    sc_uint<1> comp_r_21_43 = COMP < Q > (input[21].metric, input[43].metric);
    sc_uint<1> comp_r_21_44 = COMP < Q > (input[21].metric, input[44].metric);
    sc_uint<1> comp_r_21_45 = COMP < Q > (input[21].metric, input[45].metric);
    sc_uint<1> comp_r_21_46 = COMP < Q > (input[21].metric, input[46].metric);
    sc_uint<1> comp_r_21_47 = COMP < Q > (input[21].metric, input[47].metric);
    sc_uint<1> comp_r_21_48 = COMP < Q > (input[21].metric, input[48].metric);
    sc_uint<1> comp_r_21_49 = COMP < Q > (input[21].metric, input[49].metric);
    sc_uint<1> comp_r_21_50 = COMP < Q > (input[21].metric, input[50].metric);
    sc_uint<1> comp_r_21_51 = COMP < Q > (input[21].metric, input[51].metric);
    sc_uint<1> comp_r_21_52 = COMP < Q > (input[21].metric, input[52].metric);
    sc_uint<1> comp_r_21_53 = COMP < Q > (input[21].metric, input[53].metric);
    sc_uint<1> comp_r_21_54 = COMP < Q > (input[21].metric, input[54].metric);
    sc_uint<1> comp_r_21_55 = COMP < Q > (input[21].metric, input[55].metric);
    sc_uint<1> comp_r_21_56 = COMP < Q > (input[21].metric, input[56].metric);
    sc_uint<1> comp_r_21_57 = COMP < Q > (input[21].metric, input[57].metric);
    sc_uint<1> comp_r_21_58 = COMP < Q > (input[21].metric, input[58].metric);
    sc_uint<1> comp_r_21_59 = COMP < Q > (input[21].metric, input[59].metric);
    sc_uint<1> comp_r_21_60 = COMP < Q > (input[21].metric, input[60].metric);
    sc_uint<1> comp_r_21_61 = COMP < Q > (input[21].metric, input[61].metric);
    sc_uint<1> comp_r_21_62 = COMP < Q > (input[21].metric, input[62].metric);
    sc_uint<1> comp_r_21_63 = COMP < Q > (input[21].metric, input[63].metric);
    sc_uint<1> comp_r_21_64 = COMP < Q > (input[21].metric, input[64].metric);
    sc_uint<1> comp_r_21_65 = COMP < Q > (input[21].metric, input[65].metric);
    sc_uint<1> comp_r_21_66 = COMP < Q > (input[21].metric, input[66].metric);
    sc_uint<1> comp_r_21_67 = COMP < Q > (input[21].metric, input[67].metric);
    sc_uint<1> comp_r_21_68 = COMP < Q > (input[21].metric, input[68].metric);
    sc_uint<1> comp_r_21_69 = COMP < Q > (input[21].metric, input[69].metric);
    sc_uint<1> comp_r_21_70 = COMP < Q > (input[21].metric, input[70].metric);
    sc_uint<1> comp_r_21_71 = COMP < Q > (input[21].metric, input[71].metric);
    sc_uint<1> comp_r_21_72 = COMP < Q > (input[21].metric, input[72].metric);
    sc_uint<1> comp_r_21_73 = COMP < Q > (input[21].metric, input[73].metric);
    sc_uint<1> comp_r_21_74 = COMP < Q > (input[21].metric, input[74].metric);
    sc_uint<1> comp_r_21_75 = COMP < Q > (input[21].metric, input[75].metric);
    sc_uint<1> comp_r_21_76 = COMP < Q > (input[21].metric, input[76].metric);
    sc_uint<1> comp_r_21_77 = COMP < Q > (input[21].metric, input[77].metric);
    sc_uint<1> comp_r_21_78 = COMP < Q > (input[21].metric, input[78].metric);
    sc_uint<1> comp_r_21_79 = COMP < Q > (input[21].metric, input[79].metric);
    sc_uint<1> comp_r_21_80 = COMP < Q > (input[21].metric, input[80].metric);
    sc_uint<1> comp_r_21_81 = COMP < Q > (input[21].metric, input[81].metric);
    sc_uint<1> comp_r_21_82 = COMP < Q > (input[21].metric, input[82].metric);
    sc_uint<1> comp_r_21_83 = COMP < Q > (input[21].metric, input[83].metric);
    sc_uint<1> comp_r_21_84 = COMP < Q > (input[21].metric, input[84].metric);
    sc_uint<1> comp_r_21_85 = COMP < Q > (input[21].metric, input[85].metric);
    sc_uint<1> comp_r_21_86 = COMP < Q > (input[21].metric, input[86].metric);
    sc_uint<1> comp_r_21_87 = COMP < Q > (input[21].metric, input[87].metric);
    sc_uint<1> comp_r_21_88 = COMP < Q > (input[21].metric, input[88].metric);
    sc_uint<1> comp_r_21_89 = COMP < Q > (input[21].metric, input[89].metric);
    sc_uint<1> comp_r_21_90 = COMP < Q > (input[21].metric, input[90].metric);
    sc_uint<1> comp_r_21_91 = COMP < Q > (input[21].metric, input[91].metric);
    sc_uint<1> comp_r_21_92 = COMP < Q > (input[21].metric, input[92].metric);
    sc_uint<1> comp_r_21_93 = COMP < Q > (input[21].metric, input[93].metric);
    sc_uint<1> comp_r_21_94 = COMP < Q > (input[21].metric, input[94].metric);
    sc_uint<1> comp_r_21_95 = COMP < Q > (input[21].metric, input[95].metric);
    sc_uint<1> comp_r_21_96 = COMP < Q > (input[21].metric, input[96].metric);
    sc_uint<1> comp_r_21_97 = COMP < Q > (input[21].metric, input[97].metric);
    sc_uint<1> comp_r_21_98 = COMP < Q > (input[21].metric, input[98].metric);
    sc_uint<1> comp_r_21_99 = COMP < Q > (input[21].metric, input[99].metric);
    sc_uint<1> comp_r_21_100 = COMP < Q > (input[21].metric, input[100].metric);
    sc_uint<1> comp_r_21_101 = COMP < Q > (input[21].metric, input[101].metric);
    sc_uint<1> comp_r_21_102 = COMP < Q > (input[21].metric, input[102].metric);
    sc_uint<1> comp_r_21_103 = COMP < Q > (input[21].metric, input[103].metric);
    sc_uint<1> comp_r_21_104 = COMP < Q > (input[21].metric, input[104].metric);
    sc_uint<1> comp_r_21_105 = COMP < Q > (input[21].metric, input[105].metric);
    sc_uint<1> comp_r_21_106 = COMP < Q > (input[21].metric, input[106].metric);
    sc_uint<1> comp_r_21_107 = COMP < Q > (input[21].metric, input[107].metric);
    sc_uint<1> comp_r_21_108 = COMP < Q > (input[21].metric, input[108].metric);
    sc_uint<1> comp_r_21_109 = COMP < Q > (input[21].metric, input[109].metric);
    sc_uint<1> comp_r_21_110 = COMP < Q > (input[21].metric, input[110].metric);
    sc_uint<1> comp_r_21_111 = COMP < Q > (input[21].metric, input[111].metric);
    sc_uint<1> comp_r_21_112 = COMP < Q > (input[21].metric, input[112].metric);
    sc_uint<1> comp_r_21_113 = COMP < Q > (input[21].metric, input[113].metric);
    sc_uint<1> comp_r_21_114 = COMP < Q > (input[21].metric, input[114].metric);
    sc_uint<1> comp_r_21_115 = COMP < Q > (input[21].metric, input[115].metric);
    sc_uint<1> comp_r_21_116 = COMP < Q > (input[21].metric, input[116].metric);
    sc_uint<1> comp_r_21_117 = COMP < Q > (input[21].metric, input[117].metric);
    sc_uint<1> comp_r_21_118 = COMP < Q > (input[21].metric, input[118].metric);
    sc_uint<1> comp_r_21_119 = COMP < Q > (input[21].metric, input[119].metric);
    sc_uint<1> comp_r_21_120 = COMP < Q > (input[21].metric, input[120].metric);
    sc_uint<1> comp_r_21_121 = COMP < Q > (input[21].metric, input[121].metric);
    sc_uint<1> comp_r_21_122 = COMP < Q > (input[21].metric, input[122].metric);
    sc_uint<1> comp_r_21_123 = COMP < Q > (input[21].metric, input[123].metric);
    sc_uint<1> comp_r_21_124 = COMP < Q > (input[21].metric, input[124].metric);
    sc_uint<1> comp_r_21_125 = COMP < Q > (input[21].metric, input[125].metric);
    sc_uint<1> comp_r_21_126 = COMP < Q > (input[21].metric, input[126].metric);

    sc_uint<1> comp_r_23_24 = COMP < Q > (input[23].metric, input[24].metric);
    sc_uint<1> comp_r_23_25 = COMP < Q > (input[23].metric, input[25].metric);
    sc_uint<1> comp_r_23_26 = COMP < Q > (input[23].metric, input[26].metric);
    sc_uint<1> comp_r_23_27 = COMP < Q > (input[23].metric, input[27].metric);
    sc_uint<1> comp_r_23_28 = COMP < Q > (input[23].metric, input[28].metric);
    sc_uint<1> comp_r_23_29 = COMP < Q > (input[23].metric, input[29].metric);
    sc_uint<1> comp_r_23_30 = COMP < Q > (input[23].metric, input[30].metric);
    sc_uint<1> comp_r_23_31 = COMP < Q > (input[23].metric, input[31].metric);
    sc_uint<1> comp_r_23_32 = COMP < Q > (input[23].metric, input[32].metric);
    sc_uint<1> comp_r_23_33 = COMP < Q > (input[23].metric, input[33].metric);
    sc_uint<1> comp_r_23_34 = COMP < Q > (input[23].metric, input[34].metric);
    sc_uint<1> comp_r_23_35 = COMP < Q > (input[23].metric, input[35].metric);
    sc_uint<1> comp_r_23_36 = COMP < Q > (input[23].metric, input[36].metric);
    sc_uint<1> comp_r_23_37 = COMP < Q > (input[23].metric, input[37].metric);
    sc_uint<1> comp_r_23_38 = COMP < Q > (input[23].metric, input[38].metric);
    sc_uint<1> comp_r_23_39 = COMP < Q > (input[23].metric, input[39].metric);
    sc_uint<1> comp_r_23_40 = COMP < Q > (input[23].metric, input[40].metric);
    sc_uint<1> comp_r_23_41 = COMP < Q > (input[23].metric, input[41].metric);
    sc_uint<1> comp_r_23_42 = COMP < Q > (input[23].metric, input[42].metric);
    sc_uint<1> comp_r_23_43 = COMP < Q > (input[23].metric, input[43].metric);
    sc_uint<1> comp_r_23_44 = COMP < Q > (input[23].metric, input[44].metric);
    sc_uint<1> comp_r_23_45 = COMP < Q > (input[23].metric, input[45].metric);
    sc_uint<1> comp_r_23_46 = COMP < Q > (input[23].metric, input[46].metric);
    sc_uint<1> comp_r_23_47 = COMP < Q > (input[23].metric, input[47].metric);
    sc_uint<1> comp_r_23_48 = COMP < Q > (input[23].metric, input[48].metric);
    sc_uint<1> comp_r_23_49 = COMP < Q > (input[23].metric, input[49].metric);
    sc_uint<1> comp_r_23_50 = COMP < Q > (input[23].metric, input[50].metric);
    sc_uint<1> comp_r_23_51 = COMP < Q > (input[23].metric, input[51].metric);
    sc_uint<1> comp_r_23_52 = COMP < Q > (input[23].metric, input[52].metric);
    sc_uint<1> comp_r_23_53 = COMP < Q > (input[23].metric, input[53].metric);
    sc_uint<1> comp_r_23_54 = COMP < Q > (input[23].metric, input[54].metric);
    sc_uint<1> comp_r_23_55 = COMP < Q > (input[23].metric, input[55].metric);
    sc_uint<1> comp_r_23_56 = COMP < Q > (input[23].metric, input[56].metric);
    sc_uint<1> comp_r_23_57 = COMP < Q > (input[23].metric, input[57].metric);
    sc_uint<1> comp_r_23_58 = COMP < Q > (input[23].metric, input[58].metric);
    sc_uint<1> comp_r_23_59 = COMP < Q > (input[23].metric, input[59].metric);
    sc_uint<1> comp_r_23_60 = COMP < Q > (input[23].metric, input[60].metric);
    sc_uint<1> comp_r_23_61 = COMP < Q > (input[23].metric, input[61].metric);
    sc_uint<1> comp_r_23_62 = COMP < Q > (input[23].metric, input[62].metric);
    sc_uint<1> comp_r_23_63 = COMP < Q > (input[23].metric, input[63].metric);
    sc_uint<1> comp_r_23_64 = COMP < Q > (input[23].metric, input[64].metric);
    sc_uint<1> comp_r_23_65 = COMP < Q > (input[23].metric, input[65].metric);
    sc_uint<1> comp_r_23_66 = COMP < Q > (input[23].metric, input[66].metric);
    sc_uint<1> comp_r_23_67 = COMP < Q > (input[23].metric, input[67].metric);
    sc_uint<1> comp_r_23_68 = COMP < Q > (input[23].metric, input[68].metric);
    sc_uint<1> comp_r_23_69 = COMP < Q > (input[23].metric, input[69].metric);
    sc_uint<1> comp_r_23_70 = COMP < Q > (input[23].metric, input[70].metric);
    sc_uint<1> comp_r_23_71 = COMP < Q > (input[23].metric, input[71].metric);
    sc_uint<1> comp_r_23_72 = COMP < Q > (input[23].metric, input[72].metric);
    sc_uint<1> comp_r_23_73 = COMP < Q > (input[23].metric, input[73].metric);
    sc_uint<1> comp_r_23_74 = COMP < Q > (input[23].metric, input[74].metric);
    sc_uint<1> comp_r_23_75 = COMP < Q > (input[23].metric, input[75].metric);
    sc_uint<1> comp_r_23_76 = COMP < Q > (input[23].metric, input[76].metric);
    sc_uint<1> comp_r_23_77 = COMP < Q > (input[23].metric, input[77].metric);
    sc_uint<1> comp_r_23_78 = COMP < Q > (input[23].metric, input[78].metric);
    sc_uint<1> comp_r_23_79 = COMP < Q > (input[23].metric, input[79].metric);
    sc_uint<1> comp_r_23_80 = COMP < Q > (input[23].metric, input[80].metric);
    sc_uint<1> comp_r_23_81 = COMP < Q > (input[23].metric, input[81].metric);
    sc_uint<1> comp_r_23_82 = COMP < Q > (input[23].metric, input[82].metric);
    sc_uint<1> comp_r_23_83 = COMP < Q > (input[23].metric, input[83].metric);
    sc_uint<1> comp_r_23_84 = COMP < Q > (input[23].metric, input[84].metric);
    sc_uint<1> comp_r_23_85 = COMP < Q > (input[23].metric, input[85].metric);
    sc_uint<1> comp_r_23_86 = COMP < Q > (input[23].metric, input[86].metric);
    sc_uint<1> comp_r_23_87 = COMP < Q > (input[23].metric, input[87].metric);
    sc_uint<1> comp_r_23_88 = COMP < Q > (input[23].metric, input[88].metric);
    sc_uint<1> comp_r_23_89 = COMP < Q > (input[23].metric, input[89].metric);
    sc_uint<1> comp_r_23_90 = COMP < Q > (input[23].metric, input[90].metric);
    sc_uint<1> comp_r_23_91 = COMP < Q > (input[23].metric, input[91].metric);
    sc_uint<1> comp_r_23_92 = COMP < Q > (input[23].metric, input[92].metric);
    sc_uint<1> comp_r_23_93 = COMP < Q > (input[23].metric, input[93].metric);
    sc_uint<1> comp_r_23_94 = COMP < Q > (input[23].metric, input[94].metric);
    sc_uint<1> comp_r_23_95 = COMP < Q > (input[23].metric, input[95].metric);
    sc_uint<1> comp_r_23_96 = COMP < Q > (input[23].metric, input[96].metric);
    sc_uint<1> comp_r_23_97 = COMP < Q > (input[23].metric, input[97].metric);
    sc_uint<1> comp_r_23_98 = COMP < Q > (input[23].metric, input[98].metric);
    sc_uint<1> comp_r_23_99 = COMP < Q > (input[23].metric, input[99].metric);
    sc_uint<1> comp_r_23_100 = COMP < Q > (input[23].metric, input[100].metric);
    sc_uint<1> comp_r_23_101 = COMP < Q > (input[23].metric, input[101].metric);
    sc_uint<1> comp_r_23_102 = COMP < Q > (input[23].metric, input[102].metric);
    sc_uint<1> comp_r_23_103 = COMP < Q > (input[23].metric, input[103].metric);
    sc_uint<1> comp_r_23_104 = COMP < Q > (input[23].metric, input[104].metric);
    sc_uint<1> comp_r_23_105 = COMP < Q > (input[23].metric, input[105].metric);
    sc_uint<1> comp_r_23_106 = COMP < Q > (input[23].metric, input[106].metric);
    sc_uint<1> comp_r_23_107 = COMP < Q > (input[23].metric, input[107].metric);
    sc_uint<1> comp_r_23_108 = COMP < Q > (input[23].metric, input[108].metric);
    sc_uint<1> comp_r_23_109 = COMP < Q > (input[23].metric, input[109].metric);
    sc_uint<1> comp_r_23_110 = COMP < Q > (input[23].metric, input[110].metric);
    sc_uint<1> comp_r_23_111 = COMP < Q > (input[23].metric, input[111].metric);
    sc_uint<1> comp_r_23_112 = COMP < Q > (input[23].metric, input[112].metric);
    sc_uint<1> comp_r_23_113 = COMP < Q > (input[23].metric, input[113].metric);
    sc_uint<1> comp_r_23_114 = COMP < Q > (input[23].metric, input[114].metric);
    sc_uint<1> comp_r_23_115 = COMP < Q > (input[23].metric, input[115].metric);
    sc_uint<1> comp_r_23_116 = COMP < Q > (input[23].metric, input[116].metric);
    sc_uint<1> comp_r_23_117 = COMP < Q > (input[23].metric, input[117].metric);
    sc_uint<1> comp_r_23_118 = COMP < Q > (input[23].metric, input[118].metric);
    sc_uint<1> comp_r_23_119 = COMP < Q > (input[23].metric, input[119].metric);
    sc_uint<1> comp_r_23_120 = COMP < Q > (input[23].metric, input[120].metric);
    sc_uint<1> comp_r_23_121 = COMP < Q > (input[23].metric, input[121].metric);
    sc_uint<1> comp_r_23_122 = COMP < Q > (input[23].metric, input[122].metric);
    sc_uint<1> comp_r_23_123 = COMP < Q > (input[23].metric, input[123].metric);
    sc_uint<1> comp_r_23_124 = COMP < Q > (input[23].metric, input[124].metric);
    sc_uint<1> comp_r_23_125 = COMP < Q > (input[23].metric, input[125].metric);
    sc_uint<1> comp_r_23_126 = COMP < Q > (input[23].metric, input[126].metric);

    sc_uint<1> comp_r_25_26 = COMP < Q > (input[25].metric, input[26].metric);
    sc_uint<1> comp_r_25_27 = COMP < Q > (input[25].metric, input[27].metric);
    sc_uint<1> comp_r_25_28 = COMP < Q > (input[25].metric, input[28].metric);
    sc_uint<1> comp_r_25_29 = COMP < Q > (input[25].metric, input[29].metric);
    sc_uint<1> comp_r_25_30 = COMP < Q > (input[25].metric, input[30].metric);
    sc_uint<1> comp_r_25_31 = COMP < Q > (input[25].metric, input[31].metric);
    sc_uint<1> comp_r_25_32 = COMP < Q > (input[25].metric, input[32].metric);
    sc_uint<1> comp_r_25_33 = COMP < Q > (input[25].metric, input[33].metric);
    sc_uint<1> comp_r_25_34 = COMP < Q > (input[25].metric, input[34].metric);
    sc_uint<1> comp_r_25_35 = COMP < Q > (input[25].metric, input[35].metric);
    sc_uint<1> comp_r_25_36 = COMP < Q > (input[25].metric, input[36].metric);
    sc_uint<1> comp_r_25_37 = COMP < Q > (input[25].metric, input[37].metric);
    sc_uint<1> comp_r_25_38 = COMP < Q > (input[25].metric, input[38].metric);
    sc_uint<1> comp_r_25_39 = COMP < Q > (input[25].metric, input[39].metric);
    sc_uint<1> comp_r_25_40 = COMP < Q > (input[25].metric, input[40].metric);
    sc_uint<1> comp_r_25_41 = COMP < Q > (input[25].metric, input[41].metric);
    sc_uint<1> comp_r_25_42 = COMP < Q > (input[25].metric, input[42].metric);
    sc_uint<1> comp_r_25_43 = COMP < Q > (input[25].metric, input[43].metric);
    sc_uint<1> comp_r_25_44 = COMP < Q > (input[25].metric, input[44].metric);
    sc_uint<1> comp_r_25_45 = COMP < Q > (input[25].metric, input[45].metric);
    sc_uint<1> comp_r_25_46 = COMP < Q > (input[25].metric, input[46].metric);
    sc_uint<1> comp_r_25_47 = COMP < Q > (input[25].metric, input[47].metric);
    sc_uint<1> comp_r_25_48 = COMP < Q > (input[25].metric, input[48].metric);
    sc_uint<1> comp_r_25_49 = COMP < Q > (input[25].metric, input[49].metric);
    sc_uint<1> comp_r_25_50 = COMP < Q > (input[25].metric, input[50].metric);
    sc_uint<1> comp_r_25_51 = COMP < Q > (input[25].metric, input[51].metric);
    sc_uint<1> comp_r_25_52 = COMP < Q > (input[25].metric, input[52].metric);
    sc_uint<1> comp_r_25_53 = COMP < Q > (input[25].metric, input[53].metric);
    sc_uint<1> comp_r_25_54 = COMP < Q > (input[25].metric, input[54].metric);
    sc_uint<1> comp_r_25_55 = COMP < Q > (input[25].metric, input[55].metric);
    sc_uint<1> comp_r_25_56 = COMP < Q > (input[25].metric, input[56].metric);
    sc_uint<1> comp_r_25_57 = COMP < Q > (input[25].metric, input[57].metric);
    sc_uint<1> comp_r_25_58 = COMP < Q > (input[25].metric, input[58].metric);
    sc_uint<1> comp_r_25_59 = COMP < Q > (input[25].metric, input[59].metric);
    sc_uint<1> comp_r_25_60 = COMP < Q > (input[25].metric, input[60].metric);
    sc_uint<1> comp_r_25_61 = COMP < Q > (input[25].metric, input[61].metric);
    sc_uint<1> comp_r_25_62 = COMP < Q > (input[25].metric, input[62].metric);
    sc_uint<1> comp_r_25_63 = COMP < Q > (input[25].metric, input[63].metric);
    sc_uint<1> comp_r_25_64 = COMP < Q > (input[25].metric, input[64].metric);
    sc_uint<1> comp_r_25_65 = COMP < Q > (input[25].metric, input[65].metric);
    sc_uint<1> comp_r_25_66 = COMP < Q > (input[25].metric, input[66].metric);
    sc_uint<1> comp_r_25_67 = COMP < Q > (input[25].metric, input[67].metric);
    sc_uint<1> comp_r_25_68 = COMP < Q > (input[25].metric, input[68].metric);
    sc_uint<1> comp_r_25_69 = COMP < Q > (input[25].metric, input[69].metric);
    sc_uint<1> comp_r_25_70 = COMP < Q > (input[25].metric, input[70].metric);
    sc_uint<1> comp_r_25_71 = COMP < Q > (input[25].metric, input[71].metric);
    sc_uint<1> comp_r_25_72 = COMP < Q > (input[25].metric, input[72].metric);
    sc_uint<1> comp_r_25_73 = COMP < Q > (input[25].metric, input[73].metric);
    sc_uint<1> comp_r_25_74 = COMP < Q > (input[25].metric, input[74].metric);
    sc_uint<1> comp_r_25_75 = COMP < Q > (input[25].metric, input[75].metric);
    sc_uint<1> comp_r_25_76 = COMP < Q > (input[25].metric, input[76].metric);
    sc_uint<1> comp_r_25_77 = COMP < Q > (input[25].metric, input[77].metric);
    sc_uint<1> comp_r_25_78 = COMP < Q > (input[25].metric, input[78].metric);
    sc_uint<1> comp_r_25_79 = COMP < Q > (input[25].metric, input[79].metric);
    sc_uint<1> comp_r_25_80 = COMP < Q > (input[25].metric, input[80].metric);
    sc_uint<1> comp_r_25_81 = COMP < Q > (input[25].metric, input[81].metric);
    sc_uint<1> comp_r_25_82 = COMP < Q > (input[25].metric, input[82].metric);
    sc_uint<1> comp_r_25_83 = COMP < Q > (input[25].metric, input[83].metric);
    sc_uint<1> comp_r_25_84 = COMP < Q > (input[25].metric, input[84].metric);
    sc_uint<1> comp_r_25_85 = COMP < Q > (input[25].metric, input[85].metric);
    sc_uint<1> comp_r_25_86 = COMP < Q > (input[25].metric, input[86].metric);
    sc_uint<1> comp_r_25_87 = COMP < Q > (input[25].metric, input[87].metric);
    sc_uint<1> comp_r_25_88 = COMP < Q > (input[25].metric, input[88].metric);
    sc_uint<1> comp_r_25_89 = COMP < Q > (input[25].metric, input[89].metric);
    sc_uint<1> comp_r_25_90 = COMP < Q > (input[25].metric, input[90].metric);
    sc_uint<1> comp_r_25_91 = COMP < Q > (input[25].metric, input[91].metric);
    sc_uint<1> comp_r_25_92 = COMP < Q > (input[25].metric, input[92].metric);
    sc_uint<1> comp_r_25_93 = COMP < Q > (input[25].metric, input[93].metric);
    sc_uint<1> comp_r_25_94 = COMP < Q > (input[25].metric, input[94].metric);
    sc_uint<1> comp_r_25_95 = COMP < Q > (input[25].metric, input[95].metric);
    sc_uint<1> comp_r_25_96 = COMP < Q > (input[25].metric, input[96].metric);
    sc_uint<1> comp_r_25_97 = COMP < Q > (input[25].metric, input[97].metric);
    sc_uint<1> comp_r_25_98 = COMP < Q > (input[25].metric, input[98].metric);
    sc_uint<1> comp_r_25_99 = COMP < Q > (input[25].metric, input[99].metric);
    sc_uint<1> comp_r_25_100 = COMP < Q > (input[25].metric, input[100].metric);
    sc_uint<1> comp_r_25_101 = COMP < Q > (input[25].metric, input[101].metric);
    sc_uint<1> comp_r_25_102 = COMP < Q > (input[25].metric, input[102].metric);
    sc_uint<1> comp_r_25_103 = COMP < Q > (input[25].metric, input[103].metric);
    sc_uint<1> comp_r_25_104 = COMP < Q > (input[25].metric, input[104].metric);
    sc_uint<1> comp_r_25_105 = COMP < Q > (input[25].metric, input[105].metric);
    sc_uint<1> comp_r_25_106 = COMP < Q > (input[25].metric, input[106].metric);
    sc_uint<1> comp_r_25_107 = COMP < Q > (input[25].metric, input[107].metric);
    sc_uint<1> comp_r_25_108 = COMP < Q > (input[25].metric, input[108].metric);
    sc_uint<1> comp_r_25_109 = COMP < Q > (input[25].metric, input[109].metric);
    sc_uint<1> comp_r_25_110 = COMP < Q > (input[25].metric, input[110].metric);
    sc_uint<1> comp_r_25_111 = COMP < Q > (input[25].metric, input[111].metric);
    sc_uint<1> comp_r_25_112 = COMP < Q > (input[25].metric, input[112].metric);
    sc_uint<1> comp_r_25_113 = COMP < Q > (input[25].metric, input[113].metric);
    sc_uint<1> comp_r_25_114 = COMP < Q > (input[25].metric, input[114].metric);
    sc_uint<1> comp_r_25_115 = COMP < Q > (input[25].metric, input[115].metric);
    sc_uint<1> comp_r_25_116 = COMP < Q > (input[25].metric, input[116].metric);
    sc_uint<1> comp_r_25_117 = COMP < Q > (input[25].metric, input[117].metric);
    sc_uint<1> comp_r_25_118 = COMP < Q > (input[25].metric, input[118].metric);
    sc_uint<1> comp_r_25_119 = COMP < Q > (input[25].metric, input[119].metric);
    sc_uint<1> comp_r_25_120 = COMP < Q > (input[25].metric, input[120].metric);
    sc_uint<1> comp_r_25_121 = COMP < Q > (input[25].metric, input[121].metric);
    sc_uint<1> comp_r_25_122 = COMP < Q > (input[25].metric, input[122].metric);
    sc_uint<1> comp_r_25_123 = COMP < Q > (input[25].metric, input[123].metric);
    sc_uint<1> comp_r_25_124 = COMP < Q > (input[25].metric, input[124].metric);
    sc_uint<1> comp_r_25_125 = COMP < Q > (input[25].metric, input[125].metric);
    sc_uint<1> comp_r_25_126 = COMP < Q > (input[25].metric, input[126].metric);

    sc_uint<1> comp_r_27_28 = COMP < Q > (input[27].metric, input[28].metric);
    sc_uint<1> comp_r_27_29 = COMP < Q > (input[27].metric, input[29].metric);
    sc_uint<1> comp_r_27_30 = COMP < Q > (input[27].metric, input[30].metric);
    sc_uint<1> comp_r_27_31 = COMP < Q > (input[27].metric, input[31].metric);
    sc_uint<1> comp_r_27_32 = COMP < Q > (input[27].metric, input[32].metric);
    sc_uint<1> comp_r_27_33 = COMP < Q > (input[27].metric, input[33].metric);
    sc_uint<1> comp_r_27_34 = COMP < Q > (input[27].metric, input[34].metric);
    sc_uint<1> comp_r_27_35 = COMP < Q > (input[27].metric, input[35].metric);
    sc_uint<1> comp_r_27_36 = COMP < Q > (input[27].metric, input[36].metric);
    sc_uint<1> comp_r_27_37 = COMP < Q > (input[27].metric, input[37].metric);
    sc_uint<1> comp_r_27_38 = COMP < Q > (input[27].metric, input[38].metric);
    sc_uint<1> comp_r_27_39 = COMP < Q > (input[27].metric, input[39].metric);
    sc_uint<1> comp_r_27_40 = COMP < Q > (input[27].metric, input[40].metric);
    sc_uint<1> comp_r_27_41 = COMP < Q > (input[27].metric, input[41].metric);
    sc_uint<1> comp_r_27_42 = COMP < Q > (input[27].metric, input[42].metric);
    sc_uint<1> comp_r_27_43 = COMP < Q > (input[27].metric, input[43].metric);
    sc_uint<1> comp_r_27_44 = COMP < Q > (input[27].metric, input[44].metric);
    sc_uint<1> comp_r_27_45 = COMP < Q > (input[27].metric, input[45].metric);
    sc_uint<1> comp_r_27_46 = COMP < Q > (input[27].metric, input[46].metric);
    sc_uint<1> comp_r_27_47 = COMP < Q > (input[27].metric, input[47].metric);
    sc_uint<1> comp_r_27_48 = COMP < Q > (input[27].metric, input[48].metric);
    sc_uint<1> comp_r_27_49 = COMP < Q > (input[27].metric, input[49].metric);
    sc_uint<1> comp_r_27_50 = COMP < Q > (input[27].metric, input[50].metric);
    sc_uint<1> comp_r_27_51 = COMP < Q > (input[27].metric, input[51].metric);
    sc_uint<1> comp_r_27_52 = COMP < Q > (input[27].metric, input[52].metric);
    sc_uint<1> comp_r_27_53 = COMP < Q > (input[27].metric, input[53].metric);
    sc_uint<1> comp_r_27_54 = COMP < Q > (input[27].metric, input[54].metric);
    sc_uint<1> comp_r_27_55 = COMP < Q > (input[27].metric, input[55].metric);
    sc_uint<1> comp_r_27_56 = COMP < Q > (input[27].metric, input[56].metric);
    sc_uint<1> comp_r_27_57 = COMP < Q > (input[27].metric, input[57].metric);
    sc_uint<1> comp_r_27_58 = COMP < Q > (input[27].metric, input[58].metric);
    sc_uint<1> comp_r_27_59 = COMP < Q > (input[27].metric, input[59].metric);
    sc_uint<1> comp_r_27_60 = COMP < Q > (input[27].metric, input[60].metric);
    sc_uint<1> comp_r_27_61 = COMP < Q > (input[27].metric, input[61].metric);
    sc_uint<1> comp_r_27_62 = COMP < Q > (input[27].metric, input[62].metric);
    sc_uint<1> comp_r_27_63 = COMP < Q > (input[27].metric, input[63].metric);
    sc_uint<1> comp_r_27_64 = COMP < Q > (input[27].metric, input[64].metric);
    sc_uint<1> comp_r_27_65 = COMP < Q > (input[27].metric, input[65].metric);
    sc_uint<1> comp_r_27_66 = COMP < Q > (input[27].metric, input[66].metric);
    sc_uint<1> comp_r_27_67 = COMP < Q > (input[27].metric, input[67].metric);
    sc_uint<1> comp_r_27_68 = COMP < Q > (input[27].metric, input[68].metric);
    sc_uint<1> comp_r_27_69 = COMP < Q > (input[27].metric, input[69].metric);
    sc_uint<1> comp_r_27_70 = COMP < Q > (input[27].metric, input[70].metric);
    sc_uint<1> comp_r_27_71 = COMP < Q > (input[27].metric, input[71].metric);
    sc_uint<1> comp_r_27_72 = COMP < Q > (input[27].metric, input[72].metric);
    sc_uint<1> comp_r_27_73 = COMP < Q > (input[27].metric, input[73].metric);
    sc_uint<1> comp_r_27_74 = COMP < Q > (input[27].metric, input[74].metric);
    sc_uint<1> comp_r_27_75 = COMP < Q > (input[27].metric, input[75].metric);
    sc_uint<1> comp_r_27_76 = COMP < Q > (input[27].metric, input[76].metric);
    sc_uint<1> comp_r_27_77 = COMP < Q > (input[27].metric, input[77].metric);
    sc_uint<1> comp_r_27_78 = COMP < Q > (input[27].metric, input[78].metric);
    sc_uint<1> comp_r_27_79 = COMP < Q > (input[27].metric, input[79].metric);
    sc_uint<1> comp_r_27_80 = COMP < Q > (input[27].metric, input[80].metric);
    sc_uint<1> comp_r_27_81 = COMP < Q > (input[27].metric, input[81].metric);
    sc_uint<1> comp_r_27_82 = COMP < Q > (input[27].metric, input[82].metric);
    sc_uint<1> comp_r_27_83 = COMP < Q > (input[27].metric, input[83].metric);
    sc_uint<1> comp_r_27_84 = COMP < Q > (input[27].metric, input[84].metric);
    sc_uint<1> comp_r_27_85 = COMP < Q > (input[27].metric, input[85].metric);
    sc_uint<1> comp_r_27_86 = COMP < Q > (input[27].metric, input[86].metric);
    sc_uint<1> comp_r_27_87 = COMP < Q > (input[27].metric, input[87].metric);
    sc_uint<1> comp_r_27_88 = COMP < Q > (input[27].metric, input[88].metric);
    sc_uint<1> comp_r_27_89 = COMP < Q > (input[27].metric, input[89].metric);
    sc_uint<1> comp_r_27_90 = COMP < Q > (input[27].metric, input[90].metric);
    sc_uint<1> comp_r_27_91 = COMP < Q > (input[27].metric, input[91].metric);
    sc_uint<1> comp_r_27_92 = COMP < Q > (input[27].metric, input[92].metric);
    sc_uint<1> comp_r_27_93 = COMP < Q > (input[27].metric, input[93].metric);
    sc_uint<1> comp_r_27_94 = COMP < Q > (input[27].metric, input[94].metric);
    sc_uint<1> comp_r_27_95 = COMP < Q > (input[27].metric, input[95].metric);
    sc_uint<1> comp_r_27_96 = COMP < Q > (input[27].metric, input[96].metric);
    sc_uint<1> comp_r_27_97 = COMP < Q > (input[27].metric, input[97].metric);
    sc_uint<1> comp_r_27_98 = COMP < Q > (input[27].metric, input[98].metric);
    sc_uint<1> comp_r_27_99 = COMP < Q > (input[27].metric, input[99].metric);
    sc_uint<1> comp_r_27_100 = COMP < Q > (input[27].metric, input[100].metric);
    sc_uint<1> comp_r_27_101 = COMP < Q > (input[27].metric, input[101].metric);
    sc_uint<1> comp_r_27_102 = COMP < Q > (input[27].metric, input[102].metric);
    sc_uint<1> comp_r_27_103 = COMP < Q > (input[27].metric, input[103].metric);
    sc_uint<1> comp_r_27_104 = COMP < Q > (input[27].metric, input[104].metric);
    sc_uint<1> comp_r_27_105 = COMP < Q > (input[27].metric, input[105].metric);
    sc_uint<1> comp_r_27_106 = COMP < Q > (input[27].metric, input[106].metric);
    sc_uint<1> comp_r_27_107 = COMP < Q > (input[27].metric, input[107].metric);
    sc_uint<1> comp_r_27_108 = COMP < Q > (input[27].metric, input[108].metric);
    sc_uint<1> comp_r_27_109 = COMP < Q > (input[27].metric, input[109].metric);
    sc_uint<1> comp_r_27_110 = COMP < Q > (input[27].metric, input[110].metric);
    sc_uint<1> comp_r_27_111 = COMP < Q > (input[27].metric, input[111].metric);
    sc_uint<1> comp_r_27_112 = COMP < Q > (input[27].metric, input[112].metric);
    sc_uint<1> comp_r_27_113 = COMP < Q > (input[27].metric, input[113].metric);
    sc_uint<1> comp_r_27_114 = COMP < Q > (input[27].metric, input[114].metric);
    sc_uint<1> comp_r_27_115 = COMP < Q > (input[27].metric, input[115].metric);
    sc_uint<1> comp_r_27_116 = COMP < Q > (input[27].metric, input[116].metric);
    sc_uint<1> comp_r_27_117 = COMP < Q > (input[27].metric, input[117].metric);
    sc_uint<1> comp_r_27_118 = COMP < Q > (input[27].metric, input[118].metric);
    sc_uint<1> comp_r_27_119 = COMP < Q > (input[27].metric, input[119].metric);
    sc_uint<1> comp_r_27_120 = COMP < Q > (input[27].metric, input[120].metric);
    sc_uint<1> comp_r_27_121 = COMP < Q > (input[27].metric, input[121].metric);
    sc_uint<1> comp_r_27_122 = COMP < Q > (input[27].metric, input[122].metric);
    sc_uint<1> comp_r_27_123 = COMP < Q > (input[27].metric, input[123].metric);
    sc_uint<1> comp_r_27_124 = COMP < Q > (input[27].metric, input[124].metric);
    sc_uint<1> comp_r_27_125 = COMP < Q > (input[27].metric, input[125].metric);
    sc_uint<1> comp_r_27_126 = COMP < Q > (input[27].metric, input[126].metric);

    sc_uint<1> comp_r_29_30 = COMP < Q > (input[29].metric, input[30].metric);
    sc_uint<1> comp_r_29_31 = COMP < Q > (input[29].metric, input[31].metric);
    sc_uint<1> comp_r_29_32 = COMP < Q > (input[29].metric, input[32].metric);
    sc_uint<1> comp_r_29_33 = COMP < Q > (input[29].metric, input[33].metric);
    sc_uint<1> comp_r_29_34 = COMP < Q > (input[29].metric, input[34].metric);
    sc_uint<1> comp_r_29_35 = COMP < Q > (input[29].metric, input[35].metric);
    sc_uint<1> comp_r_29_36 = COMP < Q > (input[29].metric, input[36].metric);
    sc_uint<1> comp_r_29_37 = COMP < Q > (input[29].metric, input[37].metric);
    sc_uint<1> comp_r_29_38 = COMP < Q > (input[29].metric, input[38].metric);
    sc_uint<1> comp_r_29_39 = COMP < Q > (input[29].metric, input[39].metric);
    sc_uint<1> comp_r_29_40 = COMP < Q > (input[29].metric, input[40].metric);
    sc_uint<1> comp_r_29_41 = COMP < Q > (input[29].metric, input[41].metric);
    sc_uint<1> comp_r_29_42 = COMP < Q > (input[29].metric, input[42].metric);
    sc_uint<1> comp_r_29_43 = COMP < Q > (input[29].metric, input[43].metric);
    sc_uint<1> comp_r_29_44 = COMP < Q > (input[29].metric, input[44].metric);
    sc_uint<1> comp_r_29_45 = COMP < Q > (input[29].metric, input[45].metric);
    sc_uint<1> comp_r_29_46 = COMP < Q > (input[29].metric, input[46].metric);
    sc_uint<1> comp_r_29_47 = COMP < Q > (input[29].metric, input[47].metric);
    sc_uint<1> comp_r_29_48 = COMP < Q > (input[29].metric, input[48].metric);
    sc_uint<1> comp_r_29_49 = COMP < Q > (input[29].metric, input[49].metric);
    sc_uint<1> comp_r_29_50 = COMP < Q > (input[29].metric, input[50].metric);
    sc_uint<1> comp_r_29_51 = COMP < Q > (input[29].metric, input[51].metric);
    sc_uint<1> comp_r_29_52 = COMP < Q > (input[29].metric, input[52].metric);
    sc_uint<1> comp_r_29_53 = COMP < Q > (input[29].metric, input[53].metric);
    sc_uint<1> comp_r_29_54 = COMP < Q > (input[29].metric, input[54].metric);
    sc_uint<1> comp_r_29_55 = COMP < Q > (input[29].metric, input[55].metric);
    sc_uint<1> comp_r_29_56 = COMP < Q > (input[29].metric, input[56].metric);
    sc_uint<1> comp_r_29_57 = COMP < Q > (input[29].metric, input[57].metric);
    sc_uint<1> comp_r_29_58 = COMP < Q > (input[29].metric, input[58].metric);
    sc_uint<1> comp_r_29_59 = COMP < Q > (input[29].metric, input[59].metric);
    sc_uint<1> comp_r_29_60 = COMP < Q > (input[29].metric, input[60].metric);
    sc_uint<1> comp_r_29_61 = COMP < Q > (input[29].metric, input[61].metric);
    sc_uint<1> comp_r_29_62 = COMP < Q > (input[29].metric, input[62].metric);
    sc_uint<1> comp_r_29_63 = COMP < Q > (input[29].metric, input[63].metric);
    sc_uint<1> comp_r_29_64 = COMP < Q > (input[29].metric, input[64].metric);
    sc_uint<1> comp_r_29_65 = COMP < Q > (input[29].metric, input[65].metric);
    sc_uint<1> comp_r_29_66 = COMP < Q > (input[29].metric, input[66].metric);
    sc_uint<1> comp_r_29_67 = COMP < Q > (input[29].metric, input[67].metric);
    sc_uint<1> comp_r_29_68 = COMP < Q > (input[29].metric, input[68].metric);
    sc_uint<1> comp_r_29_69 = COMP < Q > (input[29].metric, input[69].metric);
    sc_uint<1> comp_r_29_70 = COMP < Q > (input[29].metric, input[70].metric);
    sc_uint<1> comp_r_29_71 = COMP < Q > (input[29].metric, input[71].metric);
    sc_uint<1> comp_r_29_72 = COMP < Q > (input[29].metric, input[72].metric);
    sc_uint<1> comp_r_29_73 = COMP < Q > (input[29].metric, input[73].metric);
    sc_uint<1> comp_r_29_74 = COMP < Q > (input[29].metric, input[74].metric);
    sc_uint<1> comp_r_29_75 = COMP < Q > (input[29].metric, input[75].metric);
    sc_uint<1> comp_r_29_76 = COMP < Q > (input[29].metric, input[76].metric);
    sc_uint<1> comp_r_29_77 = COMP < Q > (input[29].metric, input[77].metric);
    sc_uint<1> comp_r_29_78 = COMP < Q > (input[29].metric, input[78].metric);
    sc_uint<1> comp_r_29_79 = COMP < Q > (input[29].metric, input[79].metric);
    sc_uint<1> comp_r_29_80 = COMP < Q > (input[29].metric, input[80].metric);
    sc_uint<1> comp_r_29_81 = COMP < Q > (input[29].metric, input[81].metric);
    sc_uint<1> comp_r_29_82 = COMP < Q > (input[29].metric, input[82].metric);
    sc_uint<1> comp_r_29_83 = COMP < Q > (input[29].metric, input[83].metric);
    sc_uint<1> comp_r_29_84 = COMP < Q > (input[29].metric, input[84].metric);
    sc_uint<1> comp_r_29_85 = COMP < Q > (input[29].metric, input[85].metric);
    sc_uint<1> comp_r_29_86 = COMP < Q > (input[29].metric, input[86].metric);
    sc_uint<1> comp_r_29_87 = COMP < Q > (input[29].metric, input[87].metric);
    sc_uint<1> comp_r_29_88 = COMP < Q > (input[29].metric, input[88].metric);
    sc_uint<1> comp_r_29_89 = COMP < Q > (input[29].metric, input[89].metric);
    sc_uint<1> comp_r_29_90 = COMP < Q > (input[29].metric, input[90].metric);
    sc_uint<1> comp_r_29_91 = COMP < Q > (input[29].metric, input[91].metric);
    sc_uint<1> comp_r_29_92 = COMP < Q > (input[29].metric, input[92].metric);
    sc_uint<1> comp_r_29_93 = COMP < Q > (input[29].metric, input[93].metric);
    sc_uint<1> comp_r_29_94 = COMP < Q > (input[29].metric, input[94].metric);
    sc_uint<1> comp_r_29_95 = COMP < Q > (input[29].metric, input[95].metric);
    sc_uint<1> comp_r_29_96 = COMP < Q > (input[29].metric, input[96].metric);
    sc_uint<1> comp_r_29_97 = COMP < Q > (input[29].metric, input[97].metric);
    sc_uint<1> comp_r_29_98 = COMP < Q > (input[29].metric, input[98].metric);
    sc_uint<1> comp_r_29_99 = COMP < Q > (input[29].metric, input[99].metric);
    sc_uint<1> comp_r_29_100 = COMP < Q > (input[29].metric, input[100].metric);
    sc_uint<1> comp_r_29_101 = COMP < Q > (input[29].metric, input[101].metric);
    sc_uint<1> comp_r_29_102 = COMP < Q > (input[29].metric, input[102].metric);
    sc_uint<1> comp_r_29_103 = COMP < Q > (input[29].metric, input[103].metric);
    sc_uint<1> comp_r_29_104 = COMP < Q > (input[29].metric, input[104].metric);
    sc_uint<1> comp_r_29_105 = COMP < Q > (input[29].metric, input[105].metric);
    sc_uint<1> comp_r_29_106 = COMP < Q > (input[29].metric, input[106].metric);
    sc_uint<1> comp_r_29_107 = COMP < Q > (input[29].metric, input[107].metric);
    sc_uint<1> comp_r_29_108 = COMP < Q > (input[29].metric, input[108].metric);
    sc_uint<1> comp_r_29_109 = COMP < Q > (input[29].metric, input[109].metric);
    sc_uint<1> comp_r_29_110 = COMP < Q > (input[29].metric, input[110].metric);
    sc_uint<1> comp_r_29_111 = COMP < Q > (input[29].metric, input[111].metric);
    sc_uint<1> comp_r_29_112 = COMP < Q > (input[29].metric, input[112].metric);
    sc_uint<1> comp_r_29_113 = COMP < Q > (input[29].metric, input[113].metric);
    sc_uint<1> comp_r_29_114 = COMP < Q > (input[29].metric, input[114].metric);
    sc_uint<1> comp_r_29_115 = COMP < Q > (input[29].metric, input[115].metric);
    sc_uint<1> comp_r_29_116 = COMP < Q > (input[29].metric, input[116].metric);
    sc_uint<1> comp_r_29_117 = COMP < Q > (input[29].metric, input[117].metric);
    sc_uint<1> comp_r_29_118 = COMP < Q > (input[29].metric, input[118].metric);
    sc_uint<1> comp_r_29_119 = COMP < Q > (input[29].metric, input[119].metric);
    sc_uint<1> comp_r_29_120 = COMP < Q > (input[29].metric, input[120].metric);
    sc_uint<1> comp_r_29_121 = COMP < Q > (input[29].metric, input[121].metric);
    sc_uint<1> comp_r_29_122 = COMP < Q > (input[29].metric, input[122].metric);
    sc_uint<1> comp_r_29_123 = COMP < Q > (input[29].metric, input[123].metric);
    sc_uint<1> comp_r_29_124 = COMP < Q > (input[29].metric, input[124].metric);
    sc_uint<1> comp_r_29_125 = COMP < Q > (input[29].metric, input[125].metric);
    sc_uint<1> comp_r_29_126 = COMP < Q > (input[29].metric, input[126].metric);

    sc_uint<1> comp_r_31_32 = COMP < Q > (input[31].metric, input[32].metric);
    sc_uint<1> comp_r_31_33 = COMP < Q > (input[31].metric, input[33].metric);
    sc_uint<1> comp_r_31_34 = COMP < Q > (input[31].metric, input[34].metric);
    sc_uint<1> comp_r_31_35 = COMP < Q > (input[31].metric, input[35].metric);
    sc_uint<1> comp_r_31_36 = COMP < Q > (input[31].metric, input[36].metric);
    sc_uint<1> comp_r_31_37 = COMP < Q > (input[31].metric, input[37].metric);
    sc_uint<1> comp_r_31_38 = COMP < Q > (input[31].metric, input[38].metric);
    sc_uint<1> comp_r_31_39 = COMP < Q > (input[31].metric, input[39].metric);
    sc_uint<1> comp_r_31_40 = COMP < Q > (input[31].metric, input[40].metric);
    sc_uint<1> comp_r_31_41 = COMP < Q > (input[31].metric, input[41].metric);
    sc_uint<1> comp_r_31_42 = COMP < Q > (input[31].metric, input[42].metric);
    sc_uint<1> comp_r_31_43 = COMP < Q > (input[31].metric, input[43].metric);
    sc_uint<1> comp_r_31_44 = COMP < Q > (input[31].metric, input[44].metric);
    sc_uint<1> comp_r_31_45 = COMP < Q > (input[31].metric, input[45].metric);
    sc_uint<1> comp_r_31_46 = COMP < Q > (input[31].metric, input[46].metric);
    sc_uint<1> comp_r_31_47 = COMP < Q > (input[31].metric, input[47].metric);
    sc_uint<1> comp_r_31_48 = COMP < Q > (input[31].metric, input[48].metric);
    sc_uint<1> comp_r_31_49 = COMP < Q > (input[31].metric, input[49].metric);
    sc_uint<1> comp_r_31_50 = COMP < Q > (input[31].metric, input[50].metric);
    sc_uint<1> comp_r_31_51 = COMP < Q > (input[31].metric, input[51].metric);
    sc_uint<1> comp_r_31_52 = COMP < Q > (input[31].metric, input[52].metric);
    sc_uint<1> comp_r_31_53 = COMP < Q > (input[31].metric, input[53].metric);
    sc_uint<1> comp_r_31_54 = COMP < Q > (input[31].metric, input[54].metric);
    sc_uint<1> comp_r_31_55 = COMP < Q > (input[31].metric, input[55].metric);
    sc_uint<1> comp_r_31_56 = COMP < Q > (input[31].metric, input[56].metric);
    sc_uint<1> comp_r_31_57 = COMP < Q > (input[31].metric, input[57].metric);
    sc_uint<1> comp_r_31_58 = COMP < Q > (input[31].metric, input[58].metric);
    sc_uint<1> comp_r_31_59 = COMP < Q > (input[31].metric, input[59].metric);
    sc_uint<1> comp_r_31_60 = COMP < Q > (input[31].metric, input[60].metric);
    sc_uint<1> comp_r_31_61 = COMP < Q > (input[31].metric, input[61].metric);
    sc_uint<1> comp_r_31_62 = COMP < Q > (input[31].metric, input[62].metric);
    sc_uint<1> comp_r_31_63 = COMP < Q > (input[31].metric, input[63].metric);
    sc_uint<1> comp_r_31_64 = COMP < Q > (input[31].metric, input[64].metric);
    sc_uint<1> comp_r_31_65 = COMP < Q > (input[31].metric, input[65].metric);
    sc_uint<1> comp_r_31_66 = COMP < Q > (input[31].metric, input[66].metric);
    sc_uint<1> comp_r_31_67 = COMP < Q > (input[31].metric, input[67].metric);
    sc_uint<1> comp_r_31_68 = COMP < Q > (input[31].metric, input[68].metric);
    sc_uint<1> comp_r_31_69 = COMP < Q > (input[31].metric, input[69].metric);
    sc_uint<1> comp_r_31_70 = COMP < Q > (input[31].metric, input[70].metric);
    sc_uint<1> comp_r_31_71 = COMP < Q > (input[31].metric, input[71].metric);
    sc_uint<1> comp_r_31_72 = COMP < Q > (input[31].metric, input[72].metric);
    sc_uint<1> comp_r_31_73 = COMP < Q > (input[31].metric, input[73].metric);
    sc_uint<1> comp_r_31_74 = COMP < Q > (input[31].metric, input[74].metric);
    sc_uint<1> comp_r_31_75 = COMP < Q > (input[31].metric, input[75].metric);
    sc_uint<1> comp_r_31_76 = COMP < Q > (input[31].metric, input[76].metric);
    sc_uint<1> comp_r_31_77 = COMP < Q > (input[31].metric, input[77].metric);
    sc_uint<1> comp_r_31_78 = COMP < Q > (input[31].metric, input[78].metric);
    sc_uint<1> comp_r_31_79 = COMP < Q > (input[31].metric, input[79].metric);
    sc_uint<1> comp_r_31_80 = COMP < Q > (input[31].metric, input[80].metric);
    sc_uint<1> comp_r_31_81 = COMP < Q > (input[31].metric, input[81].metric);
    sc_uint<1> comp_r_31_82 = COMP < Q > (input[31].metric, input[82].metric);
    sc_uint<1> comp_r_31_83 = COMP < Q > (input[31].metric, input[83].metric);
    sc_uint<1> comp_r_31_84 = COMP < Q > (input[31].metric, input[84].metric);
    sc_uint<1> comp_r_31_85 = COMP < Q > (input[31].metric, input[85].metric);
    sc_uint<1> comp_r_31_86 = COMP < Q > (input[31].metric, input[86].metric);
    sc_uint<1> comp_r_31_87 = COMP < Q > (input[31].metric, input[87].metric);
    sc_uint<1> comp_r_31_88 = COMP < Q > (input[31].metric, input[88].metric);
    sc_uint<1> comp_r_31_89 = COMP < Q > (input[31].metric, input[89].metric);
    sc_uint<1> comp_r_31_90 = COMP < Q > (input[31].metric, input[90].metric);
    sc_uint<1> comp_r_31_91 = COMP < Q > (input[31].metric, input[91].metric);
    sc_uint<1> comp_r_31_92 = COMP < Q > (input[31].metric, input[92].metric);
    sc_uint<1> comp_r_31_93 = COMP < Q > (input[31].metric, input[93].metric);
    sc_uint<1> comp_r_31_94 = COMP < Q > (input[31].metric, input[94].metric);
    sc_uint<1> comp_r_31_95 = COMP < Q > (input[31].metric, input[95].metric);
    sc_uint<1> comp_r_31_96 = COMP < Q > (input[31].metric, input[96].metric);
    sc_uint<1> comp_r_31_97 = COMP < Q > (input[31].metric, input[97].metric);
    sc_uint<1> comp_r_31_98 = COMP < Q > (input[31].metric, input[98].metric);
    sc_uint<1> comp_r_31_99 = COMP < Q > (input[31].metric, input[99].metric);
    sc_uint<1> comp_r_31_100 = COMP < Q > (input[31].metric, input[100].metric);
    sc_uint<1> comp_r_31_101 = COMP < Q > (input[31].metric, input[101].metric);
    sc_uint<1> comp_r_31_102 = COMP < Q > (input[31].metric, input[102].metric);
    sc_uint<1> comp_r_31_103 = COMP < Q > (input[31].metric, input[103].metric);
    sc_uint<1> comp_r_31_104 = COMP < Q > (input[31].metric, input[104].metric);
    sc_uint<1> comp_r_31_105 = COMP < Q > (input[31].metric, input[105].metric);
    sc_uint<1> comp_r_31_106 = COMP < Q > (input[31].metric, input[106].metric);
    sc_uint<1> comp_r_31_107 = COMP < Q > (input[31].metric, input[107].metric);
    sc_uint<1> comp_r_31_108 = COMP < Q > (input[31].metric, input[108].metric);
    sc_uint<1> comp_r_31_109 = COMP < Q > (input[31].metric, input[109].metric);
    sc_uint<1> comp_r_31_110 = COMP < Q > (input[31].metric, input[110].metric);
    sc_uint<1> comp_r_31_111 = COMP < Q > (input[31].metric, input[111].metric);
    sc_uint<1> comp_r_31_112 = COMP < Q > (input[31].metric, input[112].metric);
    sc_uint<1> comp_r_31_113 = COMP < Q > (input[31].metric, input[113].metric);
    sc_uint<1> comp_r_31_114 = COMP < Q > (input[31].metric, input[114].metric);
    sc_uint<1> comp_r_31_115 = COMP < Q > (input[31].metric, input[115].metric);
    sc_uint<1> comp_r_31_116 = COMP < Q > (input[31].metric, input[116].metric);
    sc_uint<1> comp_r_31_117 = COMP < Q > (input[31].metric, input[117].metric);
    sc_uint<1> comp_r_31_118 = COMP < Q > (input[31].metric, input[118].metric);
    sc_uint<1> comp_r_31_119 = COMP < Q > (input[31].metric, input[119].metric);
    sc_uint<1> comp_r_31_120 = COMP < Q > (input[31].metric, input[120].metric);
    sc_uint<1> comp_r_31_121 = COMP < Q > (input[31].metric, input[121].metric);
    sc_uint<1> comp_r_31_122 = COMP < Q > (input[31].metric, input[122].metric);
    sc_uint<1> comp_r_31_123 = COMP < Q > (input[31].metric, input[123].metric);
    sc_uint<1> comp_r_31_124 = COMP < Q > (input[31].metric, input[124].metric);
    sc_uint<1> comp_r_31_125 = COMP < Q > (input[31].metric, input[125].metric);
    sc_uint<1> comp_r_31_126 = COMP < Q > (input[31].metric, input[126].metric);

    sc_uint<1> comp_r_33_34 = COMP < Q > (input[33].metric, input[34].metric);
    sc_uint<1> comp_r_33_35 = COMP < Q > (input[33].metric, input[35].metric);
    sc_uint<1> comp_r_33_36 = COMP < Q > (input[33].metric, input[36].metric);
    sc_uint<1> comp_r_33_37 = COMP < Q > (input[33].metric, input[37].metric);
    sc_uint<1> comp_r_33_38 = COMP < Q > (input[33].metric, input[38].metric);
    sc_uint<1> comp_r_33_39 = COMP < Q > (input[33].metric, input[39].metric);
    sc_uint<1> comp_r_33_40 = COMP < Q > (input[33].metric, input[40].metric);
    sc_uint<1> comp_r_33_41 = COMP < Q > (input[33].metric, input[41].metric);
    sc_uint<1> comp_r_33_42 = COMP < Q > (input[33].metric, input[42].metric);
    sc_uint<1> comp_r_33_43 = COMP < Q > (input[33].metric, input[43].metric);
    sc_uint<1> comp_r_33_44 = COMP < Q > (input[33].metric, input[44].metric);
    sc_uint<1> comp_r_33_45 = COMP < Q > (input[33].metric, input[45].metric);
    sc_uint<1> comp_r_33_46 = COMP < Q > (input[33].metric, input[46].metric);
    sc_uint<1> comp_r_33_47 = COMP < Q > (input[33].metric, input[47].metric);
    sc_uint<1> comp_r_33_48 = COMP < Q > (input[33].metric, input[48].metric);
    sc_uint<1> comp_r_33_49 = COMP < Q > (input[33].metric, input[49].metric);
    sc_uint<1> comp_r_33_50 = COMP < Q > (input[33].metric, input[50].metric);
    sc_uint<1> comp_r_33_51 = COMP < Q > (input[33].metric, input[51].metric);
    sc_uint<1> comp_r_33_52 = COMP < Q > (input[33].metric, input[52].metric);
    sc_uint<1> comp_r_33_53 = COMP < Q > (input[33].metric, input[53].metric);
    sc_uint<1> comp_r_33_54 = COMP < Q > (input[33].metric, input[54].metric);
    sc_uint<1> comp_r_33_55 = COMP < Q > (input[33].metric, input[55].metric);
    sc_uint<1> comp_r_33_56 = COMP < Q > (input[33].metric, input[56].metric);
    sc_uint<1> comp_r_33_57 = COMP < Q > (input[33].metric, input[57].metric);
    sc_uint<1> comp_r_33_58 = COMP < Q > (input[33].metric, input[58].metric);
    sc_uint<1> comp_r_33_59 = COMP < Q > (input[33].metric, input[59].metric);
    sc_uint<1> comp_r_33_60 = COMP < Q > (input[33].metric, input[60].metric);
    sc_uint<1> comp_r_33_61 = COMP < Q > (input[33].metric, input[61].metric);
    sc_uint<1> comp_r_33_62 = COMP < Q > (input[33].metric, input[62].metric);
    sc_uint<1> comp_r_33_63 = COMP < Q > (input[33].metric, input[63].metric);
    sc_uint<1> comp_r_33_64 = COMP < Q > (input[33].metric, input[64].metric);
    sc_uint<1> comp_r_33_65 = COMP < Q > (input[33].metric, input[65].metric);
    sc_uint<1> comp_r_33_66 = COMP < Q > (input[33].metric, input[66].metric);
    sc_uint<1> comp_r_33_67 = COMP < Q > (input[33].metric, input[67].metric);
    sc_uint<1> comp_r_33_68 = COMP < Q > (input[33].metric, input[68].metric);
    sc_uint<1> comp_r_33_69 = COMP < Q > (input[33].metric, input[69].metric);
    sc_uint<1> comp_r_33_70 = COMP < Q > (input[33].metric, input[70].metric);
    sc_uint<1> comp_r_33_71 = COMP < Q > (input[33].metric, input[71].metric);
    sc_uint<1> comp_r_33_72 = COMP < Q > (input[33].metric, input[72].metric);
    sc_uint<1> comp_r_33_73 = COMP < Q > (input[33].metric, input[73].metric);
    sc_uint<1> comp_r_33_74 = COMP < Q > (input[33].metric, input[74].metric);
    sc_uint<1> comp_r_33_75 = COMP < Q > (input[33].metric, input[75].metric);
    sc_uint<1> comp_r_33_76 = COMP < Q > (input[33].metric, input[76].metric);
    sc_uint<1> comp_r_33_77 = COMP < Q > (input[33].metric, input[77].metric);
    sc_uint<1> comp_r_33_78 = COMP < Q > (input[33].metric, input[78].metric);
    sc_uint<1> comp_r_33_79 = COMP < Q > (input[33].metric, input[79].metric);
    sc_uint<1> comp_r_33_80 = COMP < Q > (input[33].metric, input[80].metric);
    sc_uint<1> comp_r_33_81 = COMP < Q > (input[33].metric, input[81].metric);
    sc_uint<1> comp_r_33_82 = COMP < Q > (input[33].metric, input[82].metric);
    sc_uint<1> comp_r_33_83 = COMP < Q > (input[33].metric, input[83].metric);
    sc_uint<1> comp_r_33_84 = COMP < Q > (input[33].metric, input[84].metric);
    sc_uint<1> comp_r_33_85 = COMP < Q > (input[33].metric, input[85].metric);
    sc_uint<1> comp_r_33_86 = COMP < Q > (input[33].metric, input[86].metric);
    sc_uint<1> comp_r_33_87 = COMP < Q > (input[33].metric, input[87].metric);
    sc_uint<1> comp_r_33_88 = COMP < Q > (input[33].metric, input[88].metric);
    sc_uint<1> comp_r_33_89 = COMP < Q > (input[33].metric, input[89].metric);
    sc_uint<1> comp_r_33_90 = COMP < Q > (input[33].metric, input[90].metric);
    sc_uint<1> comp_r_33_91 = COMP < Q > (input[33].metric, input[91].metric);
    sc_uint<1> comp_r_33_92 = COMP < Q > (input[33].metric, input[92].metric);
    sc_uint<1> comp_r_33_93 = COMP < Q > (input[33].metric, input[93].metric);
    sc_uint<1> comp_r_33_94 = COMP < Q > (input[33].metric, input[94].metric);
    sc_uint<1> comp_r_33_95 = COMP < Q > (input[33].metric, input[95].metric);
    sc_uint<1> comp_r_33_96 = COMP < Q > (input[33].metric, input[96].metric);
    sc_uint<1> comp_r_33_97 = COMP < Q > (input[33].metric, input[97].metric);
    sc_uint<1> comp_r_33_98 = COMP < Q > (input[33].metric, input[98].metric);
    sc_uint<1> comp_r_33_99 = COMP < Q > (input[33].metric, input[99].metric);
    sc_uint<1> comp_r_33_100 = COMP < Q > (input[33].metric, input[100].metric);
    sc_uint<1> comp_r_33_101 = COMP < Q > (input[33].metric, input[101].metric);
    sc_uint<1> comp_r_33_102 = COMP < Q > (input[33].metric, input[102].metric);
    sc_uint<1> comp_r_33_103 = COMP < Q > (input[33].metric, input[103].metric);
    sc_uint<1> comp_r_33_104 = COMP < Q > (input[33].metric, input[104].metric);
    sc_uint<1> comp_r_33_105 = COMP < Q > (input[33].metric, input[105].metric);
    sc_uint<1> comp_r_33_106 = COMP < Q > (input[33].metric, input[106].metric);
    sc_uint<1> comp_r_33_107 = COMP < Q > (input[33].metric, input[107].metric);
    sc_uint<1> comp_r_33_108 = COMP < Q > (input[33].metric, input[108].metric);
    sc_uint<1> comp_r_33_109 = COMP < Q > (input[33].metric, input[109].metric);
    sc_uint<1> comp_r_33_110 = COMP < Q > (input[33].metric, input[110].metric);
    sc_uint<1> comp_r_33_111 = COMP < Q > (input[33].metric, input[111].metric);
    sc_uint<1> comp_r_33_112 = COMP < Q > (input[33].metric, input[112].metric);
    sc_uint<1> comp_r_33_113 = COMP < Q > (input[33].metric, input[113].metric);
    sc_uint<1> comp_r_33_114 = COMP < Q > (input[33].metric, input[114].metric);
    sc_uint<1> comp_r_33_115 = COMP < Q > (input[33].metric, input[115].metric);
    sc_uint<1> comp_r_33_116 = COMP < Q > (input[33].metric, input[116].metric);
    sc_uint<1> comp_r_33_117 = COMP < Q > (input[33].metric, input[117].metric);
    sc_uint<1> comp_r_33_118 = COMP < Q > (input[33].metric, input[118].metric);
    sc_uint<1> comp_r_33_119 = COMP < Q > (input[33].metric, input[119].metric);
    sc_uint<1> comp_r_33_120 = COMP < Q > (input[33].metric, input[120].metric);
    sc_uint<1> comp_r_33_121 = COMP < Q > (input[33].metric, input[121].metric);
    sc_uint<1> comp_r_33_122 = COMP < Q > (input[33].metric, input[122].metric);
    sc_uint<1> comp_r_33_123 = COMP < Q > (input[33].metric, input[123].metric);
    sc_uint<1> comp_r_33_124 = COMP < Q > (input[33].metric, input[124].metric);
    sc_uint<1> comp_r_33_125 = COMP < Q > (input[33].metric, input[125].metric);
    sc_uint<1> comp_r_33_126 = COMP < Q > (input[33].metric, input[126].metric);

    sc_uint<1> comp_r_35_36 = COMP < Q > (input[35].metric, input[36].metric);
    sc_uint<1> comp_r_35_37 = COMP < Q > (input[35].metric, input[37].metric);
    sc_uint<1> comp_r_35_38 = COMP < Q > (input[35].metric, input[38].metric);
    sc_uint<1> comp_r_35_39 = COMP < Q > (input[35].metric, input[39].metric);
    sc_uint<1> comp_r_35_40 = COMP < Q > (input[35].metric, input[40].metric);
    sc_uint<1> comp_r_35_41 = COMP < Q > (input[35].metric, input[41].metric);
    sc_uint<1> comp_r_35_42 = COMP < Q > (input[35].metric, input[42].metric);
    sc_uint<1> comp_r_35_43 = COMP < Q > (input[35].metric, input[43].metric);
    sc_uint<1> comp_r_35_44 = COMP < Q > (input[35].metric, input[44].metric);
    sc_uint<1> comp_r_35_45 = COMP < Q > (input[35].metric, input[45].metric);
    sc_uint<1> comp_r_35_46 = COMP < Q > (input[35].metric, input[46].metric);
    sc_uint<1> comp_r_35_47 = COMP < Q > (input[35].metric, input[47].metric);
    sc_uint<1> comp_r_35_48 = COMP < Q > (input[35].metric, input[48].metric);
    sc_uint<1> comp_r_35_49 = COMP < Q > (input[35].metric, input[49].metric);
    sc_uint<1> comp_r_35_50 = COMP < Q > (input[35].metric, input[50].metric);
    sc_uint<1> comp_r_35_51 = COMP < Q > (input[35].metric, input[51].metric);
    sc_uint<1> comp_r_35_52 = COMP < Q > (input[35].metric, input[52].metric);
    sc_uint<1> comp_r_35_53 = COMP < Q > (input[35].metric, input[53].metric);
    sc_uint<1> comp_r_35_54 = COMP < Q > (input[35].metric, input[54].metric);
    sc_uint<1> comp_r_35_55 = COMP < Q > (input[35].metric, input[55].metric);
    sc_uint<1> comp_r_35_56 = COMP < Q > (input[35].metric, input[56].metric);
    sc_uint<1> comp_r_35_57 = COMP < Q > (input[35].metric, input[57].metric);
    sc_uint<1> comp_r_35_58 = COMP < Q > (input[35].metric, input[58].metric);
    sc_uint<1> comp_r_35_59 = COMP < Q > (input[35].metric, input[59].metric);
    sc_uint<1> comp_r_35_60 = COMP < Q > (input[35].metric, input[60].metric);
    sc_uint<1> comp_r_35_61 = COMP < Q > (input[35].metric, input[61].metric);
    sc_uint<1> comp_r_35_62 = COMP < Q > (input[35].metric, input[62].metric);
    sc_uint<1> comp_r_35_63 = COMP < Q > (input[35].metric, input[63].metric);
    sc_uint<1> comp_r_35_64 = COMP < Q > (input[35].metric, input[64].metric);
    sc_uint<1> comp_r_35_65 = COMP < Q > (input[35].metric, input[65].metric);
    sc_uint<1> comp_r_35_66 = COMP < Q > (input[35].metric, input[66].metric);
    sc_uint<1> comp_r_35_67 = COMP < Q > (input[35].metric, input[67].metric);
    sc_uint<1> comp_r_35_68 = COMP < Q > (input[35].metric, input[68].metric);
    sc_uint<1> comp_r_35_69 = COMP < Q > (input[35].metric, input[69].metric);
    sc_uint<1> comp_r_35_70 = COMP < Q > (input[35].metric, input[70].metric);
    sc_uint<1> comp_r_35_71 = COMP < Q > (input[35].metric, input[71].metric);
    sc_uint<1> comp_r_35_72 = COMP < Q > (input[35].metric, input[72].metric);
    sc_uint<1> comp_r_35_73 = COMP < Q > (input[35].metric, input[73].metric);
    sc_uint<1> comp_r_35_74 = COMP < Q > (input[35].metric, input[74].metric);
    sc_uint<1> comp_r_35_75 = COMP < Q > (input[35].metric, input[75].metric);
    sc_uint<1> comp_r_35_76 = COMP < Q > (input[35].metric, input[76].metric);
    sc_uint<1> comp_r_35_77 = COMP < Q > (input[35].metric, input[77].metric);
    sc_uint<1> comp_r_35_78 = COMP < Q > (input[35].metric, input[78].metric);
    sc_uint<1> comp_r_35_79 = COMP < Q > (input[35].metric, input[79].metric);
    sc_uint<1> comp_r_35_80 = COMP < Q > (input[35].metric, input[80].metric);
    sc_uint<1> comp_r_35_81 = COMP < Q > (input[35].metric, input[81].metric);
    sc_uint<1> comp_r_35_82 = COMP < Q > (input[35].metric, input[82].metric);
    sc_uint<1> comp_r_35_83 = COMP < Q > (input[35].metric, input[83].metric);
    sc_uint<1> comp_r_35_84 = COMP < Q > (input[35].metric, input[84].metric);
    sc_uint<1> comp_r_35_85 = COMP < Q > (input[35].metric, input[85].metric);
    sc_uint<1> comp_r_35_86 = COMP < Q > (input[35].metric, input[86].metric);
    sc_uint<1> comp_r_35_87 = COMP < Q > (input[35].metric, input[87].metric);
    sc_uint<1> comp_r_35_88 = COMP < Q > (input[35].metric, input[88].metric);
    sc_uint<1> comp_r_35_89 = COMP < Q > (input[35].metric, input[89].metric);
    sc_uint<1> comp_r_35_90 = COMP < Q > (input[35].metric, input[90].metric);
    sc_uint<1> comp_r_35_91 = COMP < Q > (input[35].metric, input[91].metric);
    sc_uint<1> comp_r_35_92 = COMP < Q > (input[35].metric, input[92].metric);
    sc_uint<1> comp_r_35_93 = COMP < Q > (input[35].metric, input[93].metric);
    sc_uint<1> comp_r_35_94 = COMP < Q > (input[35].metric, input[94].metric);
    sc_uint<1> comp_r_35_95 = COMP < Q > (input[35].metric, input[95].metric);
    sc_uint<1> comp_r_35_96 = COMP < Q > (input[35].metric, input[96].metric);
    sc_uint<1> comp_r_35_97 = COMP < Q > (input[35].metric, input[97].metric);
    sc_uint<1> comp_r_35_98 = COMP < Q > (input[35].metric, input[98].metric);
    sc_uint<1> comp_r_35_99 = COMP < Q > (input[35].metric, input[99].metric);
    sc_uint<1> comp_r_35_100 = COMP < Q > (input[35].metric, input[100].metric);
    sc_uint<1> comp_r_35_101 = COMP < Q > (input[35].metric, input[101].metric);
    sc_uint<1> comp_r_35_102 = COMP < Q > (input[35].metric, input[102].metric);
    sc_uint<1> comp_r_35_103 = COMP < Q > (input[35].metric, input[103].metric);
    sc_uint<1> comp_r_35_104 = COMP < Q > (input[35].metric, input[104].metric);
    sc_uint<1> comp_r_35_105 = COMP < Q > (input[35].metric, input[105].metric);
    sc_uint<1> comp_r_35_106 = COMP < Q > (input[35].metric, input[106].metric);
    sc_uint<1> comp_r_35_107 = COMP < Q > (input[35].metric, input[107].metric);
    sc_uint<1> comp_r_35_108 = COMP < Q > (input[35].metric, input[108].metric);
    sc_uint<1> comp_r_35_109 = COMP < Q > (input[35].metric, input[109].metric);
    sc_uint<1> comp_r_35_110 = COMP < Q > (input[35].metric, input[110].metric);
    sc_uint<1> comp_r_35_111 = COMP < Q > (input[35].metric, input[111].metric);
    sc_uint<1> comp_r_35_112 = COMP < Q > (input[35].metric, input[112].metric);
    sc_uint<1> comp_r_35_113 = COMP < Q > (input[35].metric, input[113].metric);
    sc_uint<1> comp_r_35_114 = COMP < Q > (input[35].metric, input[114].metric);
    sc_uint<1> comp_r_35_115 = COMP < Q > (input[35].metric, input[115].metric);
    sc_uint<1> comp_r_35_116 = COMP < Q > (input[35].metric, input[116].metric);
    sc_uint<1> comp_r_35_117 = COMP < Q > (input[35].metric, input[117].metric);
    sc_uint<1> comp_r_35_118 = COMP < Q > (input[35].metric, input[118].metric);
    sc_uint<1> comp_r_35_119 = COMP < Q > (input[35].metric, input[119].metric);
    sc_uint<1> comp_r_35_120 = COMP < Q > (input[35].metric, input[120].metric);
    sc_uint<1> comp_r_35_121 = COMP < Q > (input[35].metric, input[121].metric);
    sc_uint<1> comp_r_35_122 = COMP < Q > (input[35].metric, input[122].metric);
    sc_uint<1> comp_r_35_123 = COMP < Q > (input[35].metric, input[123].metric);
    sc_uint<1> comp_r_35_124 = COMP < Q > (input[35].metric, input[124].metric);
    sc_uint<1> comp_r_35_125 = COMP < Q > (input[35].metric, input[125].metric);
    sc_uint<1> comp_r_35_126 = COMP < Q > (input[35].metric, input[126].metric);

    sc_uint<1> comp_r_37_38 = COMP < Q > (input[37].metric, input[38].metric);
    sc_uint<1> comp_r_37_39 = COMP < Q > (input[37].metric, input[39].metric);
    sc_uint<1> comp_r_37_40 = COMP < Q > (input[37].metric, input[40].metric);
    sc_uint<1> comp_r_37_41 = COMP < Q > (input[37].metric, input[41].metric);
    sc_uint<1> comp_r_37_42 = COMP < Q > (input[37].metric, input[42].metric);
    sc_uint<1> comp_r_37_43 = COMP < Q > (input[37].metric, input[43].metric);
    sc_uint<1> comp_r_37_44 = COMP < Q > (input[37].metric, input[44].metric);
    sc_uint<1> comp_r_37_45 = COMP < Q > (input[37].metric, input[45].metric);
    sc_uint<1> comp_r_37_46 = COMP < Q > (input[37].metric, input[46].metric);
    sc_uint<1> comp_r_37_47 = COMP < Q > (input[37].metric, input[47].metric);
    sc_uint<1> comp_r_37_48 = COMP < Q > (input[37].metric, input[48].metric);
    sc_uint<1> comp_r_37_49 = COMP < Q > (input[37].metric, input[49].metric);
    sc_uint<1> comp_r_37_50 = COMP < Q > (input[37].metric, input[50].metric);
    sc_uint<1> comp_r_37_51 = COMP < Q > (input[37].metric, input[51].metric);
    sc_uint<1> comp_r_37_52 = COMP < Q > (input[37].metric, input[52].metric);
    sc_uint<1> comp_r_37_53 = COMP < Q > (input[37].metric, input[53].metric);
    sc_uint<1> comp_r_37_54 = COMP < Q > (input[37].metric, input[54].metric);
    sc_uint<1> comp_r_37_55 = COMP < Q > (input[37].metric, input[55].metric);
    sc_uint<1> comp_r_37_56 = COMP < Q > (input[37].metric, input[56].metric);
    sc_uint<1> comp_r_37_57 = COMP < Q > (input[37].metric, input[57].metric);
    sc_uint<1> comp_r_37_58 = COMP < Q > (input[37].metric, input[58].metric);
    sc_uint<1> comp_r_37_59 = COMP < Q > (input[37].metric, input[59].metric);
    sc_uint<1> comp_r_37_60 = COMP < Q > (input[37].metric, input[60].metric);
    sc_uint<1> comp_r_37_61 = COMP < Q > (input[37].metric, input[61].metric);
    sc_uint<1> comp_r_37_62 = COMP < Q > (input[37].metric, input[62].metric);
    sc_uint<1> comp_r_37_63 = COMP < Q > (input[37].metric, input[63].metric);
    sc_uint<1> comp_r_37_64 = COMP < Q > (input[37].metric, input[64].metric);
    sc_uint<1> comp_r_37_65 = COMP < Q > (input[37].metric, input[65].metric);
    sc_uint<1> comp_r_37_66 = COMP < Q > (input[37].metric, input[66].metric);
    sc_uint<1> comp_r_37_67 = COMP < Q > (input[37].metric, input[67].metric);
    sc_uint<1> comp_r_37_68 = COMP < Q > (input[37].metric, input[68].metric);
    sc_uint<1> comp_r_37_69 = COMP < Q > (input[37].metric, input[69].metric);
    sc_uint<1> comp_r_37_70 = COMP < Q > (input[37].metric, input[70].metric);
    sc_uint<1> comp_r_37_71 = COMP < Q > (input[37].metric, input[71].metric);
    sc_uint<1> comp_r_37_72 = COMP < Q > (input[37].metric, input[72].metric);
    sc_uint<1> comp_r_37_73 = COMP < Q > (input[37].metric, input[73].metric);
    sc_uint<1> comp_r_37_74 = COMP < Q > (input[37].metric, input[74].metric);
    sc_uint<1> comp_r_37_75 = COMP < Q > (input[37].metric, input[75].metric);
    sc_uint<1> comp_r_37_76 = COMP < Q > (input[37].metric, input[76].metric);
    sc_uint<1> comp_r_37_77 = COMP < Q > (input[37].metric, input[77].metric);
    sc_uint<1> comp_r_37_78 = COMP < Q > (input[37].metric, input[78].metric);
    sc_uint<1> comp_r_37_79 = COMP < Q > (input[37].metric, input[79].metric);
    sc_uint<1> comp_r_37_80 = COMP < Q > (input[37].metric, input[80].metric);
    sc_uint<1> comp_r_37_81 = COMP < Q > (input[37].metric, input[81].metric);
    sc_uint<1> comp_r_37_82 = COMP < Q > (input[37].metric, input[82].metric);
    sc_uint<1> comp_r_37_83 = COMP < Q > (input[37].metric, input[83].metric);
    sc_uint<1> comp_r_37_84 = COMP < Q > (input[37].metric, input[84].metric);
    sc_uint<1> comp_r_37_85 = COMP < Q > (input[37].metric, input[85].metric);
    sc_uint<1> comp_r_37_86 = COMP < Q > (input[37].metric, input[86].metric);
    sc_uint<1> comp_r_37_87 = COMP < Q > (input[37].metric, input[87].metric);
    sc_uint<1> comp_r_37_88 = COMP < Q > (input[37].metric, input[88].metric);
    sc_uint<1> comp_r_37_89 = COMP < Q > (input[37].metric, input[89].metric);
    sc_uint<1> comp_r_37_90 = COMP < Q > (input[37].metric, input[90].metric);
    sc_uint<1> comp_r_37_91 = COMP < Q > (input[37].metric, input[91].metric);
    sc_uint<1> comp_r_37_92 = COMP < Q > (input[37].metric, input[92].metric);
    sc_uint<1> comp_r_37_93 = COMP < Q > (input[37].metric, input[93].metric);
    sc_uint<1> comp_r_37_94 = COMP < Q > (input[37].metric, input[94].metric);
    sc_uint<1> comp_r_37_95 = COMP < Q > (input[37].metric, input[95].metric);
    sc_uint<1> comp_r_37_96 = COMP < Q > (input[37].metric, input[96].metric);
    sc_uint<1> comp_r_37_97 = COMP < Q > (input[37].metric, input[97].metric);
    sc_uint<1> comp_r_37_98 = COMP < Q > (input[37].metric, input[98].metric);
    sc_uint<1> comp_r_37_99 = COMP < Q > (input[37].metric, input[99].metric);
    sc_uint<1> comp_r_37_100 = COMP < Q > (input[37].metric, input[100].metric);
    sc_uint<1> comp_r_37_101 = COMP < Q > (input[37].metric, input[101].metric);
    sc_uint<1> comp_r_37_102 = COMP < Q > (input[37].metric, input[102].metric);
    sc_uint<1> comp_r_37_103 = COMP < Q > (input[37].metric, input[103].metric);
    sc_uint<1> comp_r_37_104 = COMP < Q > (input[37].metric, input[104].metric);
    sc_uint<1> comp_r_37_105 = COMP < Q > (input[37].metric, input[105].metric);
    sc_uint<1> comp_r_37_106 = COMP < Q > (input[37].metric, input[106].metric);
    sc_uint<1> comp_r_37_107 = COMP < Q > (input[37].metric, input[107].metric);
    sc_uint<1> comp_r_37_108 = COMP < Q > (input[37].metric, input[108].metric);
    sc_uint<1> comp_r_37_109 = COMP < Q > (input[37].metric, input[109].metric);
    sc_uint<1> comp_r_37_110 = COMP < Q > (input[37].metric, input[110].metric);
    sc_uint<1> comp_r_37_111 = COMP < Q > (input[37].metric, input[111].metric);
    sc_uint<1> comp_r_37_112 = COMP < Q > (input[37].metric, input[112].metric);
    sc_uint<1> comp_r_37_113 = COMP < Q > (input[37].metric, input[113].metric);
    sc_uint<1> comp_r_37_114 = COMP < Q > (input[37].metric, input[114].metric);
    sc_uint<1> comp_r_37_115 = COMP < Q > (input[37].metric, input[115].metric);
    sc_uint<1> comp_r_37_116 = COMP < Q > (input[37].metric, input[116].metric);
    sc_uint<1> comp_r_37_117 = COMP < Q > (input[37].metric, input[117].metric);
    sc_uint<1> comp_r_37_118 = COMP < Q > (input[37].metric, input[118].metric);
    sc_uint<1> comp_r_37_119 = COMP < Q > (input[37].metric, input[119].metric);
    sc_uint<1> comp_r_37_120 = COMP < Q > (input[37].metric, input[120].metric);
    sc_uint<1> comp_r_37_121 = COMP < Q > (input[37].metric, input[121].metric);
    sc_uint<1> comp_r_37_122 = COMP < Q > (input[37].metric, input[122].metric);
    sc_uint<1> comp_r_37_123 = COMP < Q > (input[37].metric, input[123].metric);
    sc_uint<1> comp_r_37_124 = COMP < Q > (input[37].metric, input[124].metric);
    sc_uint<1> comp_r_37_125 = COMP < Q > (input[37].metric, input[125].metric);
    sc_uint<1> comp_r_37_126 = COMP < Q > (input[37].metric, input[126].metric);

    sc_uint<1> comp_r_39_40 = COMP < Q > (input[39].metric, input[40].metric);
    sc_uint<1> comp_r_39_41 = COMP < Q > (input[39].metric, input[41].metric);
    sc_uint<1> comp_r_39_42 = COMP < Q > (input[39].metric, input[42].metric);
    sc_uint<1> comp_r_39_43 = COMP < Q > (input[39].metric, input[43].metric);
    sc_uint<1> comp_r_39_44 = COMP < Q > (input[39].metric, input[44].metric);
    sc_uint<1> comp_r_39_45 = COMP < Q > (input[39].metric, input[45].metric);
    sc_uint<1> comp_r_39_46 = COMP < Q > (input[39].metric, input[46].metric);
    sc_uint<1> comp_r_39_47 = COMP < Q > (input[39].metric, input[47].metric);
    sc_uint<1> comp_r_39_48 = COMP < Q > (input[39].metric, input[48].metric);
    sc_uint<1> comp_r_39_49 = COMP < Q > (input[39].metric, input[49].metric);
    sc_uint<1> comp_r_39_50 = COMP < Q > (input[39].metric, input[50].metric);
    sc_uint<1> comp_r_39_51 = COMP < Q > (input[39].metric, input[51].metric);
    sc_uint<1> comp_r_39_52 = COMP < Q > (input[39].metric, input[52].metric);
    sc_uint<1> comp_r_39_53 = COMP < Q > (input[39].metric, input[53].metric);
    sc_uint<1> comp_r_39_54 = COMP < Q > (input[39].metric, input[54].metric);
    sc_uint<1> comp_r_39_55 = COMP < Q > (input[39].metric, input[55].metric);
    sc_uint<1> comp_r_39_56 = COMP < Q > (input[39].metric, input[56].metric);
    sc_uint<1> comp_r_39_57 = COMP < Q > (input[39].metric, input[57].metric);
    sc_uint<1> comp_r_39_58 = COMP < Q > (input[39].metric, input[58].metric);
    sc_uint<1> comp_r_39_59 = COMP < Q > (input[39].metric, input[59].metric);
    sc_uint<1> comp_r_39_60 = COMP < Q > (input[39].metric, input[60].metric);
    sc_uint<1> comp_r_39_61 = COMP < Q > (input[39].metric, input[61].metric);
    sc_uint<1> comp_r_39_62 = COMP < Q > (input[39].metric, input[62].metric);
    sc_uint<1> comp_r_39_63 = COMP < Q > (input[39].metric, input[63].metric);
    sc_uint<1> comp_r_39_64 = COMP < Q > (input[39].metric, input[64].metric);
    sc_uint<1> comp_r_39_65 = COMP < Q > (input[39].metric, input[65].metric);
    sc_uint<1> comp_r_39_66 = COMP < Q > (input[39].metric, input[66].metric);
    sc_uint<1> comp_r_39_67 = COMP < Q > (input[39].metric, input[67].metric);
    sc_uint<1> comp_r_39_68 = COMP < Q > (input[39].metric, input[68].metric);
    sc_uint<1> comp_r_39_69 = COMP < Q > (input[39].metric, input[69].metric);
    sc_uint<1> comp_r_39_70 = COMP < Q > (input[39].metric, input[70].metric);
    sc_uint<1> comp_r_39_71 = COMP < Q > (input[39].metric, input[71].metric);
    sc_uint<1> comp_r_39_72 = COMP < Q > (input[39].metric, input[72].metric);
    sc_uint<1> comp_r_39_73 = COMP < Q > (input[39].metric, input[73].metric);
    sc_uint<1> comp_r_39_74 = COMP < Q > (input[39].metric, input[74].metric);
    sc_uint<1> comp_r_39_75 = COMP < Q > (input[39].metric, input[75].metric);
    sc_uint<1> comp_r_39_76 = COMP < Q > (input[39].metric, input[76].metric);
    sc_uint<1> comp_r_39_77 = COMP < Q > (input[39].metric, input[77].metric);
    sc_uint<1> comp_r_39_78 = COMP < Q > (input[39].metric, input[78].metric);
    sc_uint<1> comp_r_39_79 = COMP < Q > (input[39].metric, input[79].metric);
    sc_uint<1> comp_r_39_80 = COMP < Q > (input[39].metric, input[80].metric);
    sc_uint<1> comp_r_39_81 = COMP < Q > (input[39].metric, input[81].metric);
    sc_uint<1> comp_r_39_82 = COMP < Q > (input[39].metric, input[82].metric);
    sc_uint<1> comp_r_39_83 = COMP < Q > (input[39].metric, input[83].metric);
    sc_uint<1> comp_r_39_84 = COMP < Q > (input[39].metric, input[84].metric);
    sc_uint<1> comp_r_39_85 = COMP < Q > (input[39].metric, input[85].metric);
    sc_uint<1> comp_r_39_86 = COMP < Q > (input[39].metric, input[86].metric);
    sc_uint<1> comp_r_39_87 = COMP < Q > (input[39].metric, input[87].metric);
    sc_uint<1> comp_r_39_88 = COMP < Q > (input[39].metric, input[88].metric);
    sc_uint<1> comp_r_39_89 = COMP < Q > (input[39].metric, input[89].metric);
    sc_uint<1> comp_r_39_90 = COMP < Q > (input[39].metric, input[90].metric);
    sc_uint<1> comp_r_39_91 = COMP < Q > (input[39].metric, input[91].metric);
    sc_uint<1> comp_r_39_92 = COMP < Q > (input[39].metric, input[92].metric);
    sc_uint<1> comp_r_39_93 = COMP < Q > (input[39].metric, input[93].metric);
    sc_uint<1> comp_r_39_94 = COMP < Q > (input[39].metric, input[94].metric);
    sc_uint<1> comp_r_39_95 = COMP < Q > (input[39].metric, input[95].metric);
    sc_uint<1> comp_r_39_96 = COMP < Q > (input[39].metric, input[96].metric);
    sc_uint<1> comp_r_39_97 = COMP < Q > (input[39].metric, input[97].metric);
    sc_uint<1> comp_r_39_98 = COMP < Q > (input[39].metric, input[98].metric);
    sc_uint<1> comp_r_39_99 = COMP < Q > (input[39].metric, input[99].metric);
    sc_uint<1> comp_r_39_100 = COMP < Q > (input[39].metric, input[100].metric);
    sc_uint<1> comp_r_39_101 = COMP < Q > (input[39].metric, input[101].metric);
    sc_uint<1> comp_r_39_102 = COMP < Q > (input[39].metric, input[102].metric);
    sc_uint<1> comp_r_39_103 = COMP < Q > (input[39].metric, input[103].metric);
    sc_uint<1> comp_r_39_104 = COMP < Q > (input[39].metric, input[104].metric);
    sc_uint<1> comp_r_39_105 = COMP < Q > (input[39].metric, input[105].metric);
    sc_uint<1> comp_r_39_106 = COMP < Q > (input[39].metric, input[106].metric);
    sc_uint<1> comp_r_39_107 = COMP < Q > (input[39].metric, input[107].metric);
    sc_uint<1> comp_r_39_108 = COMP < Q > (input[39].metric, input[108].metric);
    sc_uint<1> comp_r_39_109 = COMP < Q > (input[39].metric, input[109].metric);
    sc_uint<1> comp_r_39_110 = COMP < Q > (input[39].metric, input[110].metric);
    sc_uint<1> comp_r_39_111 = COMP < Q > (input[39].metric, input[111].metric);
    sc_uint<1> comp_r_39_112 = COMP < Q > (input[39].metric, input[112].metric);
    sc_uint<1> comp_r_39_113 = COMP < Q > (input[39].metric, input[113].metric);
    sc_uint<1> comp_r_39_114 = COMP < Q > (input[39].metric, input[114].metric);
    sc_uint<1> comp_r_39_115 = COMP < Q > (input[39].metric, input[115].metric);
    sc_uint<1> comp_r_39_116 = COMP < Q > (input[39].metric, input[116].metric);
    sc_uint<1> comp_r_39_117 = COMP < Q > (input[39].metric, input[117].metric);
    sc_uint<1> comp_r_39_118 = COMP < Q > (input[39].metric, input[118].metric);
    sc_uint<1> comp_r_39_119 = COMP < Q > (input[39].metric, input[119].metric);
    sc_uint<1> comp_r_39_120 = COMP < Q > (input[39].metric, input[120].metric);
    sc_uint<1> comp_r_39_121 = COMP < Q > (input[39].metric, input[121].metric);
    sc_uint<1> comp_r_39_122 = COMP < Q > (input[39].metric, input[122].metric);
    sc_uint<1> comp_r_39_123 = COMP < Q > (input[39].metric, input[123].metric);
    sc_uint<1> comp_r_39_124 = COMP < Q > (input[39].metric, input[124].metric);
    sc_uint<1> comp_r_39_125 = COMP < Q > (input[39].metric, input[125].metric);
    sc_uint<1> comp_r_39_126 = COMP < Q > (input[39].metric, input[126].metric);

    sc_uint<1> comp_r_41_42 = COMP < Q > (input[41].metric, input[42].metric);
    sc_uint<1> comp_r_41_43 = COMP < Q > (input[41].metric, input[43].metric);
    sc_uint<1> comp_r_41_44 = COMP < Q > (input[41].metric, input[44].metric);
    sc_uint<1> comp_r_41_45 = COMP < Q > (input[41].metric, input[45].metric);
    sc_uint<1> comp_r_41_46 = COMP < Q > (input[41].metric, input[46].metric);
    sc_uint<1> comp_r_41_47 = COMP < Q > (input[41].metric, input[47].metric);
    sc_uint<1> comp_r_41_48 = COMP < Q > (input[41].metric, input[48].metric);
    sc_uint<1> comp_r_41_49 = COMP < Q > (input[41].metric, input[49].metric);
    sc_uint<1> comp_r_41_50 = COMP < Q > (input[41].metric, input[50].metric);
    sc_uint<1> comp_r_41_51 = COMP < Q > (input[41].metric, input[51].metric);
    sc_uint<1> comp_r_41_52 = COMP < Q > (input[41].metric, input[52].metric);
    sc_uint<1> comp_r_41_53 = COMP < Q > (input[41].metric, input[53].metric);
    sc_uint<1> comp_r_41_54 = COMP < Q > (input[41].metric, input[54].metric);
    sc_uint<1> comp_r_41_55 = COMP < Q > (input[41].metric, input[55].metric);
    sc_uint<1> comp_r_41_56 = COMP < Q > (input[41].metric, input[56].metric);
    sc_uint<1> comp_r_41_57 = COMP < Q > (input[41].metric, input[57].metric);
    sc_uint<1> comp_r_41_58 = COMP < Q > (input[41].metric, input[58].metric);
    sc_uint<1> comp_r_41_59 = COMP < Q > (input[41].metric, input[59].metric);
    sc_uint<1> comp_r_41_60 = COMP < Q > (input[41].metric, input[60].metric);
    sc_uint<1> comp_r_41_61 = COMP < Q > (input[41].metric, input[61].metric);
    sc_uint<1> comp_r_41_62 = COMP < Q > (input[41].metric, input[62].metric);
    sc_uint<1> comp_r_41_63 = COMP < Q > (input[41].metric, input[63].metric);
    sc_uint<1> comp_r_41_64 = COMP < Q > (input[41].metric, input[64].metric);
    sc_uint<1> comp_r_41_65 = COMP < Q > (input[41].metric, input[65].metric);
    sc_uint<1> comp_r_41_66 = COMP < Q > (input[41].metric, input[66].metric);
    sc_uint<1> comp_r_41_67 = COMP < Q > (input[41].metric, input[67].metric);
    sc_uint<1> comp_r_41_68 = COMP < Q > (input[41].metric, input[68].metric);
    sc_uint<1> comp_r_41_69 = COMP < Q > (input[41].metric, input[69].metric);
    sc_uint<1> comp_r_41_70 = COMP < Q > (input[41].metric, input[70].metric);
    sc_uint<1> comp_r_41_71 = COMP < Q > (input[41].metric, input[71].metric);
    sc_uint<1> comp_r_41_72 = COMP < Q > (input[41].metric, input[72].metric);
    sc_uint<1> comp_r_41_73 = COMP < Q > (input[41].metric, input[73].metric);
    sc_uint<1> comp_r_41_74 = COMP < Q > (input[41].metric, input[74].metric);
    sc_uint<1> comp_r_41_75 = COMP < Q > (input[41].metric, input[75].metric);
    sc_uint<1> comp_r_41_76 = COMP < Q > (input[41].metric, input[76].metric);
    sc_uint<1> comp_r_41_77 = COMP < Q > (input[41].metric, input[77].metric);
    sc_uint<1> comp_r_41_78 = COMP < Q > (input[41].metric, input[78].metric);
    sc_uint<1> comp_r_41_79 = COMP < Q > (input[41].metric, input[79].metric);
    sc_uint<1> comp_r_41_80 = COMP < Q > (input[41].metric, input[80].metric);
    sc_uint<1> comp_r_41_81 = COMP < Q > (input[41].metric, input[81].metric);
    sc_uint<1> comp_r_41_82 = COMP < Q > (input[41].metric, input[82].metric);
    sc_uint<1> comp_r_41_83 = COMP < Q > (input[41].metric, input[83].metric);
    sc_uint<1> comp_r_41_84 = COMP < Q > (input[41].metric, input[84].metric);
    sc_uint<1> comp_r_41_85 = COMP < Q > (input[41].metric, input[85].metric);
    sc_uint<1> comp_r_41_86 = COMP < Q > (input[41].metric, input[86].metric);
    sc_uint<1> comp_r_41_87 = COMP < Q > (input[41].metric, input[87].metric);
    sc_uint<1> comp_r_41_88 = COMP < Q > (input[41].metric, input[88].metric);
    sc_uint<1> comp_r_41_89 = COMP < Q > (input[41].metric, input[89].metric);
    sc_uint<1> comp_r_41_90 = COMP < Q > (input[41].metric, input[90].metric);
    sc_uint<1> comp_r_41_91 = COMP < Q > (input[41].metric, input[91].metric);
    sc_uint<1> comp_r_41_92 = COMP < Q > (input[41].metric, input[92].metric);
    sc_uint<1> comp_r_41_93 = COMP < Q > (input[41].metric, input[93].metric);
    sc_uint<1> comp_r_41_94 = COMP < Q > (input[41].metric, input[94].metric);
    sc_uint<1> comp_r_41_95 = COMP < Q > (input[41].metric, input[95].metric);
    sc_uint<1> comp_r_41_96 = COMP < Q > (input[41].metric, input[96].metric);
    sc_uint<1> comp_r_41_97 = COMP < Q > (input[41].metric, input[97].metric);
    sc_uint<1> comp_r_41_98 = COMP < Q > (input[41].metric, input[98].metric);
    sc_uint<1> comp_r_41_99 = COMP < Q > (input[41].metric, input[99].metric);
    sc_uint<1> comp_r_41_100 = COMP < Q > (input[41].metric, input[100].metric);
    sc_uint<1> comp_r_41_101 = COMP < Q > (input[41].metric, input[101].metric);
    sc_uint<1> comp_r_41_102 = COMP < Q > (input[41].metric, input[102].metric);
    sc_uint<1> comp_r_41_103 = COMP < Q > (input[41].metric, input[103].metric);
    sc_uint<1> comp_r_41_104 = COMP < Q > (input[41].metric, input[104].metric);
    sc_uint<1> comp_r_41_105 = COMP < Q > (input[41].metric, input[105].metric);
    sc_uint<1> comp_r_41_106 = COMP < Q > (input[41].metric, input[106].metric);
    sc_uint<1> comp_r_41_107 = COMP < Q > (input[41].metric, input[107].metric);
    sc_uint<1> comp_r_41_108 = COMP < Q > (input[41].metric, input[108].metric);
    sc_uint<1> comp_r_41_109 = COMP < Q > (input[41].metric, input[109].metric);
    sc_uint<1> comp_r_41_110 = COMP < Q > (input[41].metric, input[110].metric);
    sc_uint<1> comp_r_41_111 = COMP < Q > (input[41].metric, input[111].metric);
    sc_uint<1> comp_r_41_112 = COMP < Q > (input[41].metric, input[112].metric);
    sc_uint<1> comp_r_41_113 = COMP < Q > (input[41].metric, input[113].metric);
    sc_uint<1> comp_r_41_114 = COMP < Q > (input[41].metric, input[114].metric);
    sc_uint<1> comp_r_41_115 = COMP < Q > (input[41].metric, input[115].metric);
    sc_uint<1> comp_r_41_116 = COMP < Q > (input[41].metric, input[116].metric);
    sc_uint<1> comp_r_41_117 = COMP < Q > (input[41].metric, input[117].metric);
    sc_uint<1> comp_r_41_118 = COMP < Q > (input[41].metric, input[118].metric);
    sc_uint<1> comp_r_41_119 = COMP < Q > (input[41].metric, input[119].metric);
    sc_uint<1> comp_r_41_120 = COMP < Q > (input[41].metric, input[120].metric);
    sc_uint<1> comp_r_41_121 = COMP < Q > (input[41].metric, input[121].metric);
    sc_uint<1> comp_r_41_122 = COMP < Q > (input[41].metric, input[122].metric);
    sc_uint<1> comp_r_41_123 = COMP < Q > (input[41].metric, input[123].metric);
    sc_uint<1> comp_r_41_124 = COMP < Q > (input[41].metric, input[124].metric);
    sc_uint<1> comp_r_41_125 = COMP < Q > (input[41].metric, input[125].metric);
    sc_uint<1> comp_r_41_126 = COMP < Q > (input[41].metric, input[126].metric);

    sc_uint<1> comp_r_43_44 = COMP < Q > (input[43].metric, input[44].metric);
    sc_uint<1> comp_r_43_45 = COMP < Q > (input[43].metric, input[45].metric);
    sc_uint<1> comp_r_43_46 = COMP < Q > (input[43].metric, input[46].metric);
    sc_uint<1> comp_r_43_47 = COMP < Q > (input[43].metric, input[47].metric);
    sc_uint<1> comp_r_43_48 = COMP < Q > (input[43].metric, input[48].metric);
    sc_uint<1> comp_r_43_49 = COMP < Q > (input[43].metric, input[49].metric);
    sc_uint<1> comp_r_43_50 = COMP < Q > (input[43].metric, input[50].metric);
    sc_uint<1> comp_r_43_51 = COMP < Q > (input[43].metric, input[51].metric);
    sc_uint<1> comp_r_43_52 = COMP < Q > (input[43].metric, input[52].metric);
    sc_uint<1> comp_r_43_53 = COMP < Q > (input[43].metric, input[53].metric);
    sc_uint<1> comp_r_43_54 = COMP < Q > (input[43].metric, input[54].metric);
    sc_uint<1> comp_r_43_55 = COMP < Q > (input[43].metric, input[55].metric);
    sc_uint<1> comp_r_43_56 = COMP < Q > (input[43].metric, input[56].metric);
    sc_uint<1> comp_r_43_57 = COMP < Q > (input[43].metric, input[57].metric);
    sc_uint<1> comp_r_43_58 = COMP < Q > (input[43].metric, input[58].metric);
    sc_uint<1> comp_r_43_59 = COMP < Q > (input[43].metric, input[59].metric);
    sc_uint<1> comp_r_43_60 = COMP < Q > (input[43].metric, input[60].metric);
    sc_uint<1> comp_r_43_61 = COMP < Q > (input[43].metric, input[61].metric);
    sc_uint<1> comp_r_43_62 = COMP < Q > (input[43].metric, input[62].metric);
    sc_uint<1> comp_r_43_63 = COMP < Q > (input[43].metric, input[63].metric);
    sc_uint<1> comp_r_43_64 = COMP < Q > (input[43].metric, input[64].metric);
    sc_uint<1> comp_r_43_65 = COMP < Q > (input[43].metric, input[65].metric);
    sc_uint<1> comp_r_43_66 = COMP < Q > (input[43].metric, input[66].metric);
    sc_uint<1> comp_r_43_67 = COMP < Q > (input[43].metric, input[67].metric);
    sc_uint<1> comp_r_43_68 = COMP < Q > (input[43].metric, input[68].metric);
    sc_uint<1> comp_r_43_69 = COMP < Q > (input[43].metric, input[69].metric);
    sc_uint<1> comp_r_43_70 = COMP < Q > (input[43].metric, input[70].metric);
    sc_uint<1> comp_r_43_71 = COMP < Q > (input[43].metric, input[71].metric);
    sc_uint<1> comp_r_43_72 = COMP < Q > (input[43].metric, input[72].metric);
    sc_uint<1> comp_r_43_73 = COMP < Q > (input[43].metric, input[73].metric);
    sc_uint<1> comp_r_43_74 = COMP < Q > (input[43].metric, input[74].metric);
    sc_uint<1> comp_r_43_75 = COMP < Q > (input[43].metric, input[75].metric);
    sc_uint<1> comp_r_43_76 = COMP < Q > (input[43].metric, input[76].metric);
    sc_uint<1> comp_r_43_77 = COMP < Q > (input[43].metric, input[77].metric);
    sc_uint<1> comp_r_43_78 = COMP < Q > (input[43].metric, input[78].metric);
    sc_uint<1> comp_r_43_79 = COMP < Q > (input[43].metric, input[79].metric);
    sc_uint<1> comp_r_43_80 = COMP < Q > (input[43].metric, input[80].metric);
    sc_uint<1> comp_r_43_81 = COMP < Q > (input[43].metric, input[81].metric);
    sc_uint<1> comp_r_43_82 = COMP < Q > (input[43].metric, input[82].metric);
    sc_uint<1> comp_r_43_83 = COMP < Q > (input[43].metric, input[83].metric);
    sc_uint<1> comp_r_43_84 = COMP < Q > (input[43].metric, input[84].metric);
    sc_uint<1> comp_r_43_85 = COMP < Q > (input[43].metric, input[85].metric);
    sc_uint<1> comp_r_43_86 = COMP < Q > (input[43].metric, input[86].metric);
    sc_uint<1> comp_r_43_87 = COMP < Q > (input[43].metric, input[87].metric);
    sc_uint<1> comp_r_43_88 = COMP < Q > (input[43].metric, input[88].metric);
    sc_uint<1> comp_r_43_89 = COMP < Q > (input[43].metric, input[89].metric);
    sc_uint<1> comp_r_43_90 = COMP < Q > (input[43].metric, input[90].metric);
    sc_uint<1> comp_r_43_91 = COMP < Q > (input[43].metric, input[91].metric);
    sc_uint<1> comp_r_43_92 = COMP < Q > (input[43].metric, input[92].metric);
    sc_uint<1> comp_r_43_93 = COMP < Q > (input[43].metric, input[93].metric);
    sc_uint<1> comp_r_43_94 = COMP < Q > (input[43].metric, input[94].metric);
    sc_uint<1> comp_r_43_95 = COMP < Q > (input[43].metric, input[95].metric);
    sc_uint<1> comp_r_43_96 = COMP < Q > (input[43].metric, input[96].metric);
    sc_uint<1> comp_r_43_97 = COMP < Q > (input[43].metric, input[97].metric);
    sc_uint<1> comp_r_43_98 = COMP < Q > (input[43].metric, input[98].metric);
    sc_uint<1> comp_r_43_99 = COMP < Q > (input[43].metric, input[99].metric);
    sc_uint<1> comp_r_43_100 = COMP < Q > (input[43].metric, input[100].metric);
    sc_uint<1> comp_r_43_101 = COMP < Q > (input[43].metric, input[101].metric);
    sc_uint<1> comp_r_43_102 = COMP < Q > (input[43].metric, input[102].metric);
    sc_uint<1> comp_r_43_103 = COMP < Q > (input[43].metric, input[103].metric);
    sc_uint<1> comp_r_43_104 = COMP < Q > (input[43].metric, input[104].metric);
    sc_uint<1> comp_r_43_105 = COMP < Q > (input[43].metric, input[105].metric);
    sc_uint<1> comp_r_43_106 = COMP < Q > (input[43].metric, input[106].metric);
    sc_uint<1> comp_r_43_107 = COMP < Q > (input[43].metric, input[107].metric);
    sc_uint<1> comp_r_43_108 = COMP < Q > (input[43].metric, input[108].metric);
    sc_uint<1> comp_r_43_109 = COMP < Q > (input[43].metric, input[109].metric);
    sc_uint<1> comp_r_43_110 = COMP < Q > (input[43].metric, input[110].metric);
    sc_uint<1> comp_r_43_111 = COMP < Q > (input[43].metric, input[111].metric);
    sc_uint<1> comp_r_43_112 = COMP < Q > (input[43].metric, input[112].metric);
    sc_uint<1> comp_r_43_113 = COMP < Q > (input[43].metric, input[113].metric);
    sc_uint<1> comp_r_43_114 = COMP < Q > (input[43].metric, input[114].metric);
    sc_uint<1> comp_r_43_115 = COMP < Q > (input[43].metric, input[115].metric);
    sc_uint<1> comp_r_43_116 = COMP < Q > (input[43].metric, input[116].metric);
    sc_uint<1> comp_r_43_117 = COMP < Q > (input[43].metric, input[117].metric);
    sc_uint<1> comp_r_43_118 = COMP < Q > (input[43].metric, input[118].metric);
    sc_uint<1> comp_r_43_119 = COMP < Q > (input[43].metric, input[119].metric);
    sc_uint<1> comp_r_43_120 = COMP < Q > (input[43].metric, input[120].metric);
    sc_uint<1> comp_r_43_121 = COMP < Q > (input[43].metric, input[121].metric);
    sc_uint<1> comp_r_43_122 = COMP < Q > (input[43].metric, input[122].metric);
    sc_uint<1> comp_r_43_123 = COMP < Q > (input[43].metric, input[123].metric);
    sc_uint<1> comp_r_43_124 = COMP < Q > (input[43].metric, input[124].metric);
    sc_uint<1> comp_r_43_125 = COMP < Q > (input[43].metric, input[125].metric);
    sc_uint<1> comp_r_43_126 = COMP < Q > (input[43].metric, input[126].metric);

    sc_uint<1> comp_r_45_46 = COMP < Q > (input[45].metric, input[46].metric);
    sc_uint<1> comp_r_45_47 = COMP < Q > (input[45].metric, input[47].metric);
    sc_uint<1> comp_r_45_48 = COMP < Q > (input[45].metric, input[48].metric);
    sc_uint<1> comp_r_45_49 = COMP < Q > (input[45].metric, input[49].metric);
    sc_uint<1> comp_r_45_50 = COMP < Q > (input[45].metric, input[50].metric);
    sc_uint<1> comp_r_45_51 = COMP < Q > (input[45].metric, input[51].metric);
    sc_uint<1> comp_r_45_52 = COMP < Q > (input[45].metric, input[52].metric);
    sc_uint<1> comp_r_45_53 = COMP < Q > (input[45].metric, input[53].metric);
    sc_uint<1> comp_r_45_54 = COMP < Q > (input[45].metric, input[54].metric);
    sc_uint<1> comp_r_45_55 = COMP < Q > (input[45].metric, input[55].metric);
    sc_uint<1> comp_r_45_56 = COMP < Q > (input[45].metric, input[56].metric);
    sc_uint<1> comp_r_45_57 = COMP < Q > (input[45].metric, input[57].metric);
    sc_uint<1> comp_r_45_58 = COMP < Q > (input[45].metric, input[58].metric);
    sc_uint<1> comp_r_45_59 = COMP < Q > (input[45].metric, input[59].metric);
    sc_uint<1> comp_r_45_60 = COMP < Q > (input[45].metric, input[60].metric);
    sc_uint<1> comp_r_45_61 = COMP < Q > (input[45].metric, input[61].metric);
    sc_uint<1> comp_r_45_62 = COMP < Q > (input[45].metric, input[62].metric);
    sc_uint<1> comp_r_45_63 = COMP < Q > (input[45].metric, input[63].metric);
    sc_uint<1> comp_r_45_64 = COMP < Q > (input[45].metric, input[64].metric);
    sc_uint<1> comp_r_45_65 = COMP < Q > (input[45].metric, input[65].metric);
    sc_uint<1> comp_r_45_66 = COMP < Q > (input[45].metric, input[66].metric);
    sc_uint<1> comp_r_45_67 = COMP < Q > (input[45].metric, input[67].metric);
    sc_uint<1> comp_r_45_68 = COMP < Q > (input[45].metric, input[68].metric);
    sc_uint<1> comp_r_45_69 = COMP < Q > (input[45].metric, input[69].metric);
    sc_uint<1> comp_r_45_70 = COMP < Q > (input[45].metric, input[70].metric);
    sc_uint<1> comp_r_45_71 = COMP < Q > (input[45].metric, input[71].metric);
    sc_uint<1> comp_r_45_72 = COMP < Q > (input[45].metric, input[72].metric);
    sc_uint<1> comp_r_45_73 = COMP < Q > (input[45].metric, input[73].metric);
    sc_uint<1> comp_r_45_74 = COMP < Q > (input[45].metric, input[74].metric);
    sc_uint<1> comp_r_45_75 = COMP < Q > (input[45].metric, input[75].metric);
    sc_uint<1> comp_r_45_76 = COMP < Q > (input[45].metric, input[76].metric);
    sc_uint<1> comp_r_45_77 = COMP < Q > (input[45].metric, input[77].metric);
    sc_uint<1> comp_r_45_78 = COMP < Q > (input[45].metric, input[78].metric);
    sc_uint<1> comp_r_45_79 = COMP < Q > (input[45].metric, input[79].metric);
    sc_uint<1> comp_r_45_80 = COMP < Q > (input[45].metric, input[80].metric);
    sc_uint<1> comp_r_45_81 = COMP < Q > (input[45].metric, input[81].metric);
    sc_uint<1> comp_r_45_82 = COMP < Q > (input[45].metric, input[82].metric);
    sc_uint<1> comp_r_45_83 = COMP < Q > (input[45].metric, input[83].metric);
    sc_uint<1> comp_r_45_84 = COMP < Q > (input[45].metric, input[84].metric);
    sc_uint<1> comp_r_45_85 = COMP < Q > (input[45].metric, input[85].metric);
    sc_uint<1> comp_r_45_86 = COMP < Q > (input[45].metric, input[86].metric);
    sc_uint<1> comp_r_45_87 = COMP < Q > (input[45].metric, input[87].metric);
    sc_uint<1> comp_r_45_88 = COMP < Q > (input[45].metric, input[88].metric);
    sc_uint<1> comp_r_45_89 = COMP < Q > (input[45].metric, input[89].metric);
    sc_uint<1> comp_r_45_90 = COMP < Q > (input[45].metric, input[90].metric);
    sc_uint<1> comp_r_45_91 = COMP < Q > (input[45].metric, input[91].metric);
    sc_uint<1> comp_r_45_92 = COMP < Q > (input[45].metric, input[92].metric);
    sc_uint<1> comp_r_45_93 = COMP < Q > (input[45].metric, input[93].metric);
    sc_uint<1> comp_r_45_94 = COMP < Q > (input[45].metric, input[94].metric);
    sc_uint<1> comp_r_45_95 = COMP < Q > (input[45].metric, input[95].metric);
    sc_uint<1> comp_r_45_96 = COMP < Q > (input[45].metric, input[96].metric);
    sc_uint<1> comp_r_45_97 = COMP < Q > (input[45].metric, input[97].metric);
    sc_uint<1> comp_r_45_98 = COMP < Q > (input[45].metric, input[98].metric);
    sc_uint<1> comp_r_45_99 = COMP < Q > (input[45].metric, input[99].metric);
    sc_uint<1> comp_r_45_100 = COMP < Q > (input[45].metric, input[100].metric);
    sc_uint<1> comp_r_45_101 = COMP < Q > (input[45].metric, input[101].metric);
    sc_uint<1> comp_r_45_102 = COMP < Q > (input[45].metric, input[102].metric);
    sc_uint<1> comp_r_45_103 = COMP < Q > (input[45].metric, input[103].metric);
    sc_uint<1> comp_r_45_104 = COMP < Q > (input[45].metric, input[104].metric);
    sc_uint<1> comp_r_45_105 = COMP < Q > (input[45].metric, input[105].metric);
    sc_uint<1> comp_r_45_106 = COMP < Q > (input[45].metric, input[106].metric);
    sc_uint<1> comp_r_45_107 = COMP < Q > (input[45].metric, input[107].metric);
    sc_uint<1> comp_r_45_108 = COMP < Q > (input[45].metric, input[108].metric);
    sc_uint<1> comp_r_45_109 = COMP < Q > (input[45].metric, input[109].metric);
    sc_uint<1> comp_r_45_110 = COMP < Q > (input[45].metric, input[110].metric);
    sc_uint<1> comp_r_45_111 = COMP < Q > (input[45].metric, input[111].metric);
    sc_uint<1> comp_r_45_112 = COMP < Q > (input[45].metric, input[112].metric);
    sc_uint<1> comp_r_45_113 = COMP < Q > (input[45].metric, input[113].metric);
    sc_uint<1> comp_r_45_114 = COMP < Q > (input[45].metric, input[114].metric);
    sc_uint<1> comp_r_45_115 = COMP < Q > (input[45].metric, input[115].metric);
    sc_uint<1> comp_r_45_116 = COMP < Q > (input[45].metric, input[116].metric);
    sc_uint<1> comp_r_45_117 = COMP < Q > (input[45].metric, input[117].metric);
    sc_uint<1> comp_r_45_118 = COMP < Q > (input[45].metric, input[118].metric);
    sc_uint<1> comp_r_45_119 = COMP < Q > (input[45].metric, input[119].metric);
    sc_uint<1> comp_r_45_120 = COMP < Q > (input[45].metric, input[120].metric);
    sc_uint<1> comp_r_45_121 = COMP < Q > (input[45].metric, input[121].metric);
    sc_uint<1> comp_r_45_122 = COMP < Q > (input[45].metric, input[122].metric);
    sc_uint<1> comp_r_45_123 = COMP < Q > (input[45].metric, input[123].metric);
    sc_uint<1> comp_r_45_124 = COMP < Q > (input[45].metric, input[124].metric);
    sc_uint<1> comp_r_45_125 = COMP < Q > (input[45].metric, input[125].metric);
    sc_uint<1> comp_r_45_126 = COMP < Q > (input[45].metric, input[126].metric);

    sc_uint<1> comp_r_47_48 = COMP < Q > (input[47].metric, input[48].metric);
    sc_uint<1> comp_r_47_49 = COMP < Q > (input[47].metric, input[49].metric);
    sc_uint<1> comp_r_47_50 = COMP < Q > (input[47].metric, input[50].metric);
    sc_uint<1> comp_r_47_51 = COMP < Q > (input[47].metric, input[51].metric);
    sc_uint<1> comp_r_47_52 = COMP < Q > (input[47].metric, input[52].metric);
    sc_uint<1> comp_r_47_53 = COMP < Q > (input[47].metric, input[53].metric);
    sc_uint<1> comp_r_47_54 = COMP < Q > (input[47].metric, input[54].metric);
    sc_uint<1> comp_r_47_55 = COMP < Q > (input[47].metric, input[55].metric);
    sc_uint<1> comp_r_47_56 = COMP < Q > (input[47].metric, input[56].metric);
    sc_uint<1> comp_r_47_57 = COMP < Q > (input[47].metric, input[57].metric);
    sc_uint<1> comp_r_47_58 = COMP < Q > (input[47].metric, input[58].metric);
    sc_uint<1> comp_r_47_59 = COMP < Q > (input[47].metric, input[59].metric);
    sc_uint<1> comp_r_47_60 = COMP < Q > (input[47].metric, input[60].metric);
    sc_uint<1> comp_r_47_61 = COMP < Q > (input[47].metric, input[61].metric);
    sc_uint<1> comp_r_47_62 = COMP < Q > (input[47].metric, input[62].metric);
    sc_uint<1> comp_r_47_63 = COMP < Q > (input[47].metric, input[63].metric);
    sc_uint<1> comp_r_47_64 = COMP < Q > (input[47].metric, input[64].metric);
    sc_uint<1> comp_r_47_65 = COMP < Q > (input[47].metric, input[65].metric);
    sc_uint<1> comp_r_47_66 = COMP < Q > (input[47].metric, input[66].metric);
    sc_uint<1> comp_r_47_67 = COMP < Q > (input[47].metric, input[67].metric);
    sc_uint<1> comp_r_47_68 = COMP < Q > (input[47].metric, input[68].metric);
    sc_uint<1> comp_r_47_69 = COMP < Q > (input[47].metric, input[69].metric);
    sc_uint<1> comp_r_47_70 = COMP < Q > (input[47].metric, input[70].metric);
    sc_uint<1> comp_r_47_71 = COMP < Q > (input[47].metric, input[71].metric);
    sc_uint<1> comp_r_47_72 = COMP < Q > (input[47].metric, input[72].metric);
    sc_uint<1> comp_r_47_73 = COMP < Q > (input[47].metric, input[73].metric);
    sc_uint<1> comp_r_47_74 = COMP < Q > (input[47].metric, input[74].metric);
    sc_uint<1> comp_r_47_75 = COMP < Q > (input[47].metric, input[75].metric);
    sc_uint<1> comp_r_47_76 = COMP < Q > (input[47].metric, input[76].metric);
    sc_uint<1> comp_r_47_77 = COMP < Q > (input[47].metric, input[77].metric);
    sc_uint<1> comp_r_47_78 = COMP < Q > (input[47].metric, input[78].metric);
    sc_uint<1> comp_r_47_79 = COMP < Q > (input[47].metric, input[79].metric);
    sc_uint<1> comp_r_47_80 = COMP < Q > (input[47].metric, input[80].metric);
    sc_uint<1> comp_r_47_81 = COMP < Q > (input[47].metric, input[81].metric);
    sc_uint<1> comp_r_47_82 = COMP < Q > (input[47].metric, input[82].metric);
    sc_uint<1> comp_r_47_83 = COMP < Q > (input[47].metric, input[83].metric);
    sc_uint<1> comp_r_47_84 = COMP < Q > (input[47].metric, input[84].metric);
    sc_uint<1> comp_r_47_85 = COMP < Q > (input[47].metric, input[85].metric);
    sc_uint<1> comp_r_47_86 = COMP < Q > (input[47].metric, input[86].metric);
    sc_uint<1> comp_r_47_87 = COMP < Q > (input[47].metric, input[87].metric);
    sc_uint<1> comp_r_47_88 = COMP < Q > (input[47].metric, input[88].metric);
    sc_uint<1> comp_r_47_89 = COMP < Q > (input[47].metric, input[89].metric);
    sc_uint<1> comp_r_47_90 = COMP < Q > (input[47].metric, input[90].metric);
    sc_uint<1> comp_r_47_91 = COMP < Q > (input[47].metric, input[91].metric);
    sc_uint<1> comp_r_47_92 = COMP < Q > (input[47].metric, input[92].metric);
    sc_uint<1> comp_r_47_93 = COMP < Q > (input[47].metric, input[93].metric);
    sc_uint<1> comp_r_47_94 = COMP < Q > (input[47].metric, input[94].metric);
    sc_uint<1> comp_r_47_95 = COMP < Q > (input[47].metric, input[95].metric);
    sc_uint<1> comp_r_47_96 = COMP < Q > (input[47].metric, input[96].metric);
    sc_uint<1> comp_r_47_97 = COMP < Q > (input[47].metric, input[97].metric);
    sc_uint<1> comp_r_47_98 = COMP < Q > (input[47].metric, input[98].metric);
    sc_uint<1> comp_r_47_99 = COMP < Q > (input[47].metric, input[99].metric);
    sc_uint<1> comp_r_47_100 = COMP < Q > (input[47].metric, input[100].metric);
    sc_uint<1> comp_r_47_101 = COMP < Q > (input[47].metric, input[101].metric);
    sc_uint<1> comp_r_47_102 = COMP < Q > (input[47].metric, input[102].metric);
    sc_uint<1> comp_r_47_103 = COMP < Q > (input[47].metric, input[103].metric);
    sc_uint<1> comp_r_47_104 = COMP < Q > (input[47].metric, input[104].metric);
    sc_uint<1> comp_r_47_105 = COMP < Q > (input[47].metric, input[105].metric);
    sc_uint<1> comp_r_47_106 = COMP < Q > (input[47].metric, input[106].metric);
    sc_uint<1> comp_r_47_107 = COMP < Q > (input[47].metric, input[107].metric);
    sc_uint<1> comp_r_47_108 = COMP < Q > (input[47].metric, input[108].metric);
    sc_uint<1> comp_r_47_109 = COMP < Q > (input[47].metric, input[109].metric);
    sc_uint<1> comp_r_47_110 = COMP < Q > (input[47].metric, input[110].metric);
    sc_uint<1> comp_r_47_111 = COMP < Q > (input[47].metric, input[111].metric);
    sc_uint<1> comp_r_47_112 = COMP < Q > (input[47].metric, input[112].metric);
    sc_uint<1> comp_r_47_113 = COMP < Q > (input[47].metric, input[113].metric);
    sc_uint<1> comp_r_47_114 = COMP < Q > (input[47].metric, input[114].metric);
    sc_uint<1> comp_r_47_115 = COMP < Q > (input[47].metric, input[115].metric);
    sc_uint<1> comp_r_47_116 = COMP < Q > (input[47].metric, input[116].metric);
    sc_uint<1> comp_r_47_117 = COMP < Q > (input[47].metric, input[117].metric);
    sc_uint<1> comp_r_47_118 = COMP < Q > (input[47].metric, input[118].metric);
    sc_uint<1> comp_r_47_119 = COMP < Q > (input[47].metric, input[119].metric);
    sc_uint<1> comp_r_47_120 = COMP < Q > (input[47].metric, input[120].metric);
    sc_uint<1> comp_r_47_121 = COMP < Q > (input[47].metric, input[121].metric);
    sc_uint<1> comp_r_47_122 = COMP < Q > (input[47].metric, input[122].metric);
    sc_uint<1> comp_r_47_123 = COMP < Q > (input[47].metric, input[123].metric);
    sc_uint<1> comp_r_47_124 = COMP < Q > (input[47].metric, input[124].metric);
    sc_uint<1> comp_r_47_125 = COMP < Q > (input[47].metric, input[125].metric);
    sc_uint<1> comp_r_47_126 = COMP < Q > (input[47].metric, input[126].metric);

    sc_uint<1> comp_r_49_50 = COMP < Q > (input[49].metric, input[50].metric);
    sc_uint<1> comp_r_49_51 = COMP < Q > (input[49].metric, input[51].metric);
    sc_uint<1> comp_r_49_52 = COMP < Q > (input[49].metric, input[52].metric);
    sc_uint<1> comp_r_49_53 = COMP < Q > (input[49].metric, input[53].metric);
    sc_uint<1> comp_r_49_54 = COMP < Q > (input[49].metric, input[54].metric);
    sc_uint<1> comp_r_49_55 = COMP < Q > (input[49].metric, input[55].metric);
    sc_uint<1> comp_r_49_56 = COMP < Q > (input[49].metric, input[56].metric);
    sc_uint<1> comp_r_49_57 = COMP < Q > (input[49].metric, input[57].metric);
    sc_uint<1> comp_r_49_58 = COMP < Q > (input[49].metric, input[58].metric);
    sc_uint<1> comp_r_49_59 = COMP < Q > (input[49].metric, input[59].metric);
    sc_uint<1> comp_r_49_60 = COMP < Q > (input[49].metric, input[60].metric);
    sc_uint<1> comp_r_49_61 = COMP < Q > (input[49].metric, input[61].metric);
    sc_uint<1> comp_r_49_62 = COMP < Q > (input[49].metric, input[62].metric);
    sc_uint<1> comp_r_49_63 = COMP < Q > (input[49].metric, input[63].metric);
    sc_uint<1> comp_r_49_64 = COMP < Q > (input[49].metric, input[64].metric);
    sc_uint<1> comp_r_49_65 = COMP < Q > (input[49].metric, input[65].metric);
    sc_uint<1> comp_r_49_66 = COMP < Q > (input[49].metric, input[66].metric);
    sc_uint<1> comp_r_49_67 = COMP < Q > (input[49].metric, input[67].metric);
    sc_uint<1> comp_r_49_68 = COMP < Q > (input[49].metric, input[68].metric);
    sc_uint<1> comp_r_49_69 = COMP < Q > (input[49].metric, input[69].metric);
    sc_uint<1> comp_r_49_70 = COMP < Q > (input[49].metric, input[70].metric);
    sc_uint<1> comp_r_49_71 = COMP < Q > (input[49].metric, input[71].metric);
    sc_uint<1> comp_r_49_72 = COMP < Q > (input[49].metric, input[72].metric);
    sc_uint<1> comp_r_49_73 = COMP < Q > (input[49].metric, input[73].metric);
    sc_uint<1> comp_r_49_74 = COMP < Q > (input[49].metric, input[74].metric);
    sc_uint<1> comp_r_49_75 = COMP < Q > (input[49].metric, input[75].metric);
    sc_uint<1> comp_r_49_76 = COMP < Q > (input[49].metric, input[76].metric);
    sc_uint<1> comp_r_49_77 = COMP < Q > (input[49].metric, input[77].metric);
    sc_uint<1> comp_r_49_78 = COMP < Q > (input[49].metric, input[78].metric);
    sc_uint<1> comp_r_49_79 = COMP < Q > (input[49].metric, input[79].metric);
    sc_uint<1> comp_r_49_80 = COMP < Q > (input[49].metric, input[80].metric);
    sc_uint<1> comp_r_49_81 = COMP < Q > (input[49].metric, input[81].metric);
    sc_uint<1> comp_r_49_82 = COMP < Q > (input[49].metric, input[82].metric);
    sc_uint<1> comp_r_49_83 = COMP < Q > (input[49].metric, input[83].metric);
    sc_uint<1> comp_r_49_84 = COMP < Q > (input[49].metric, input[84].metric);
    sc_uint<1> comp_r_49_85 = COMP < Q > (input[49].metric, input[85].metric);
    sc_uint<1> comp_r_49_86 = COMP < Q > (input[49].metric, input[86].metric);
    sc_uint<1> comp_r_49_87 = COMP < Q > (input[49].metric, input[87].metric);
    sc_uint<1> comp_r_49_88 = COMP < Q > (input[49].metric, input[88].metric);
    sc_uint<1> comp_r_49_89 = COMP < Q > (input[49].metric, input[89].metric);
    sc_uint<1> comp_r_49_90 = COMP < Q > (input[49].metric, input[90].metric);
    sc_uint<1> comp_r_49_91 = COMP < Q > (input[49].metric, input[91].metric);
    sc_uint<1> comp_r_49_92 = COMP < Q > (input[49].metric, input[92].metric);
    sc_uint<1> comp_r_49_93 = COMP < Q > (input[49].metric, input[93].metric);
    sc_uint<1> comp_r_49_94 = COMP < Q > (input[49].metric, input[94].metric);
    sc_uint<1> comp_r_49_95 = COMP < Q > (input[49].metric, input[95].metric);
    sc_uint<1> comp_r_49_96 = COMP < Q > (input[49].metric, input[96].metric);
    sc_uint<1> comp_r_49_97 = COMP < Q > (input[49].metric, input[97].metric);
    sc_uint<1> comp_r_49_98 = COMP < Q > (input[49].metric, input[98].metric);
    sc_uint<1> comp_r_49_99 = COMP < Q > (input[49].metric, input[99].metric);
    sc_uint<1> comp_r_49_100 = COMP < Q > (input[49].metric, input[100].metric);
    sc_uint<1> comp_r_49_101 = COMP < Q > (input[49].metric, input[101].metric);
    sc_uint<1> comp_r_49_102 = COMP < Q > (input[49].metric, input[102].metric);
    sc_uint<1> comp_r_49_103 = COMP < Q > (input[49].metric, input[103].metric);
    sc_uint<1> comp_r_49_104 = COMP < Q > (input[49].metric, input[104].metric);
    sc_uint<1> comp_r_49_105 = COMP < Q > (input[49].metric, input[105].metric);
    sc_uint<1> comp_r_49_106 = COMP < Q > (input[49].metric, input[106].metric);
    sc_uint<1> comp_r_49_107 = COMP < Q > (input[49].metric, input[107].metric);
    sc_uint<1> comp_r_49_108 = COMP < Q > (input[49].metric, input[108].metric);
    sc_uint<1> comp_r_49_109 = COMP < Q > (input[49].metric, input[109].metric);
    sc_uint<1> comp_r_49_110 = COMP < Q > (input[49].metric, input[110].metric);
    sc_uint<1> comp_r_49_111 = COMP < Q > (input[49].metric, input[111].metric);
    sc_uint<1> comp_r_49_112 = COMP < Q > (input[49].metric, input[112].metric);
    sc_uint<1> comp_r_49_113 = COMP < Q > (input[49].metric, input[113].metric);
    sc_uint<1> comp_r_49_114 = COMP < Q > (input[49].metric, input[114].metric);
    sc_uint<1> comp_r_49_115 = COMP < Q > (input[49].metric, input[115].metric);
    sc_uint<1> comp_r_49_116 = COMP < Q > (input[49].metric, input[116].metric);
    sc_uint<1> comp_r_49_117 = COMP < Q > (input[49].metric, input[117].metric);
    sc_uint<1> comp_r_49_118 = COMP < Q > (input[49].metric, input[118].metric);
    sc_uint<1> comp_r_49_119 = COMP < Q > (input[49].metric, input[119].metric);
    sc_uint<1> comp_r_49_120 = COMP < Q > (input[49].metric, input[120].metric);
    sc_uint<1> comp_r_49_121 = COMP < Q > (input[49].metric, input[121].metric);
    sc_uint<1> comp_r_49_122 = COMP < Q > (input[49].metric, input[122].metric);
    sc_uint<1> comp_r_49_123 = COMP < Q > (input[49].metric, input[123].metric);
    sc_uint<1> comp_r_49_124 = COMP < Q > (input[49].metric, input[124].metric);
    sc_uint<1> comp_r_49_125 = COMP < Q > (input[49].metric, input[125].metric);
    sc_uint<1> comp_r_49_126 = COMP < Q > (input[49].metric, input[126].metric);

    sc_uint<1> comp_r_51_52 = COMP < Q > (input[51].metric, input[52].metric);
    sc_uint<1> comp_r_51_53 = COMP < Q > (input[51].metric, input[53].metric);
    sc_uint<1> comp_r_51_54 = COMP < Q > (input[51].metric, input[54].metric);
    sc_uint<1> comp_r_51_55 = COMP < Q > (input[51].metric, input[55].metric);
    sc_uint<1> comp_r_51_56 = COMP < Q > (input[51].metric, input[56].metric);
    sc_uint<1> comp_r_51_57 = COMP < Q > (input[51].metric, input[57].metric);
    sc_uint<1> comp_r_51_58 = COMP < Q > (input[51].metric, input[58].metric);
    sc_uint<1> comp_r_51_59 = COMP < Q > (input[51].metric, input[59].metric);
    sc_uint<1> comp_r_51_60 = COMP < Q > (input[51].metric, input[60].metric);
    sc_uint<1> comp_r_51_61 = COMP < Q > (input[51].metric, input[61].metric);
    sc_uint<1> comp_r_51_62 = COMP < Q > (input[51].metric, input[62].metric);
    sc_uint<1> comp_r_51_63 = COMP < Q > (input[51].metric, input[63].metric);
    sc_uint<1> comp_r_51_64 = COMP < Q > (input[51].metric, input[64].metric);
    sc_uint<1> comp_r_51_65 = COMP < Q > (input[51].metric, input[65].metric);
    sc_uint<1> comp_r_51_66 = COMP < Q > (input[51].metric, input[66].metric);
    sc_uint<1> comp_r_51_67 = COMP < Q > (input[51].metric, input[67].metric);
    sc_uint<1> comp_r_51_68 = COMP < Q > (input[51].metric, input[68].metric);
    sc_uint<1> comp_r_51_69 = COMP < Q > (input[51].metric, input[69].metric);
    sc_uint<1> comp_r_51_70 = COMP < Q > (input[51].metric, input[70].metric);
    sc_uint<1> comp_r_51_71 = COMP < Q > (input[51].metric, input[71].metric);
    sc_uint<1> comp_r_51_72 = COMP < Q > (input[51].metric, input[72].metric);
    sc_uint<1> comp_r_51_73 = COMP < Q > (input[51].metric, input[73].metric);
    sc_uint<1> comp_r_51_74 = COMP < Q > (input[51].metric, input[74].metric);
    sc_uint<1> comp_r_51_75 = COMP < Q > (input[51].metric, input[75].metric);
    sc_uint<1> comp_r_51_76 = COMP < Q > (input[51].metric, input[76].metric);
    sc_uint<1> comp_r_51_77 = COMP < Q > (input[51].metric, input[77].metric);
    sc_uint<1> comp_r_51_78 = COMP < Q > (input[51].metric, input[78].metric);
    sc_uint<1> comp_r_51_79 = COMP < Q > (input[51].metric, input[79].metric);
    sc_uint<1> comp_r_51_80 = COMP < Q > (input[51].metric, input[80].metric);
    sc_uint<1> comp_r_51_81 = COMP < Q > (input[51].metric, input[81].metric);
    sc_uint<1> comp_r_51_82 = COMP < Q > (input[51].metric, input[82].metric);
    sc_uint<1> comp_r_51_83 = COMP < Q > (input[51].metric, input[83].metric);
    sc_uint<1> comp_r_51_84 = COMP < Q > (input[51].metric, input[84].metric);
    sc_uint<1> comp_r_51_85 = COMP < Q > (input[51].metric, input[85].metric);
    sc_uint<1> comp_r_51_86 = COMP < Q > (input[51].metric, input[86].metric);
    sc_uint<1> comp_r_51_87 = COMP < Q > (input[51].metric, input[87].metric);
    sc_uint<1> comp_r_51_88 = COMP < Q > (input[51].metric, input[88].metric);
    sc_uint<1> comp_r_51_89 = COMP < Q > (input[51].metric, input[89].metric);
    sc_uint<1> comp_r_51_90 = COMP < Q > (input[51].metric, input[90].metric);
    sc_uint<1> comp_r_51_91 = COMP < Q > (input[51].metric, input[91].metric);
    sc_uint<1> comp_r_51_92 = COMP < Q > (input[51].metric, input[92].metric);
    sc_uint<1> comp_r_51_93 = COMP < Q > (input[51].metric, input[93].metric);
    sc_uint<1> comp_r_51_94 = COMP < Q > (input[51].metric, input[94].metric);
    sc_uint<1> comp_r_51_95 = COMP < Q > (input[51].metric, input[95].metric);
    sc_uint<1> comp_r_51_96 = COMP < Q > (input[51].metric, input[96].metric);
    sc_uint<1> comp_r_51_97 = COMP < Q > (input[51].metric, input[97].metric);
    sc_uint<1> comp_r_51_98 = COMP < Q > (input[51].metric, input[98].metric);
    sc_uint<1> comp_r_51_99 = COMP < Q > (input[51].metric, input[99].metric);
    sc_uint<1> comp_r_51_100 = COMP < Q > (input[51].metric, input[100].metric);
    sc_uint<1> comp_r_51_101 = COMP < Q > (input[51].metric, input[101].metric);
    sc_uint<1> comp_r_51_102 = COMP < Q > (input[51].metric, input[102].metric);
    sc_uint<1> comp_r_51_103 = COMP < Q > (input[51].metric, input[103].metric);
    sc_uint<1> comp_r_51_104 = COMP < Q > (input[51].metric, input[104].metric);
    sc_uint<1> comp_r_51_105 = COMP < Q > (input[51].metric, input[105].metric);
    sc_uint<1> comp_r_51_106 = COMP < Q > (input[51].metric, input[106].metric);
    sc_uint<1> comp_r_51_107 = COMP < Q > (input[51].metric, input[107].metric);
    sc_uint<1> comp_r_51_108 = COMP < Q > (input[51].metric, input[108].metric);
    sc_uint<1> comp_r_51_109 = COMP < Q > (input[51].metric, input[109].metric);
    sc_uint<1> comp_r_51_110 = COMP < Q > (input[51].metric, input[110].metric);
    sc_uint<1> comp_r_51_111 = COMP < Q > (input[51].metric, input[111].metric);
    sc_uint<1> comp_r_51_112 = COMP < Q > (input[51].metric, input[112].metric);
    sc_uint<1> comp_r_51_113 = COMP < Q > (input[51].metric, input[113].metric);
    sc_uint<1> comp_r_51_114 = COMP < Q > (input[51].metric, input[114].metric);
    sc_uint<1> comp_r_51_115 = COMP < Q > (input[51].metric, input[115].metric);
    sc_uint<1> comp_r_51_116 = COMP < Q > (input[51].metric, input[116].metric);
    sc_uint<1> comp_r_51_117 = COMP < Q > (input[51].metric, input[117].metric);
    sc_uint<1> comp_r_51_118 = COMP < Q > (input[51].metric, input[118].metric);
    sc_uint<1> comp_r_51_119 = COMP < Q > (input[51].metric, input[119].metric);
    sc_uint<1> comp_r_51_120 = COMP < Q > (input[51].metric, input[120].metric);
    sc_uint<1> comp_r_51_121 = COMP < Q > (input[51].metric, input[121].metric);
    sc_uint<1> comp_r_51_122 = COMP < Q > (input[51].metric, input[122].metric);
    sc_uint<1> comp_r_51_123 = COMP < Q > (input[51].metric, input[123].metric);
    sc_uint<1> comp_r_51_124 = COMP < Q > (input[51].metric, input[124].metric);
    sc_uint<1> comp_r_51_125 = COMP < Q > (input[51].metric, input[125].metric);
    sc_uint<1> comp_r_51_126 = COMP < Q > (input[51].metric, input[126].metric);

    sc_uint<1> comp_r_53_54 = COMP < Q > (input[53].metric, input[54].metric);
    sc_uint<1> comp_r_53_55 = COMP < Q > (input[53].metric, input[55].metric);
    sc_uint<1> comp_r_53_56 = COMP < Q > (input[53].metric, input[56].metric);
    sc_uint<1> comp_r_53_57 = COMP < Q > (input[53].metric, input[57].metric);
    sc_uint<1> comp_r_53_58 = COMP < Q > (input[53].metric, input[58].metric);
    sc_uint<1> comp_r_53_59 = COMP < Q > (input[53].metric, input[59].metric);
    sc_uint<1> comp_r_53_60 = COMP < Q > (input[53].metric, input[60].metric);
    sc_uint<1> comp_r_53_61 = COMP < Q > (input[53].metric, input[61].metric);
    sc_uint<1> comp_r_53_62 = COMP < Q > (input[53].metric, input[62].metric);
    sc_uint<1> comp_r_53_63 = COMP < Q > (input[53].metric, input[63].metric);
    sc_uint<1> comp_r_53_64 = COMP < Q > (input[53].metric, input[64].metric);
    sc_uint<1> comp_r_53_65 = COMP < Q > (input[53].metric, input[65].metric);
    sc_uint<1> comp_r_53_66 = COMP < Q > (input[53].metric, input[66].metric);
    sc_uint<1> comp_r_53_67 = COMP < Q > (input[53].metric, input[67].metric);
    sc_uint<1> comp_r_53_68 = COMP < Q > (input[53].metric, input[68].metric);
    sc_uint<1> comp_r_53_69 = COMP < Q > (input[53].metric, input[69].metric);
    sc_uint<1> comp_r_53_70 = COMP < Q > (input[53].metric, input[70].metric);
    sc_uint<1> comp_r_53_71 = COMP < Q > (input[53].metric, input[71].metric);
    sc_uint<1> comp_r_53_72 = COMP < Q > (input[53].metric, input[72].metric);
    sc_uint<1> comp_r_53_73 = COMP < Q > (input[53].metric, input[73].metric);
    sc_uint<1> comp_r_53_74 = COMP < Q > (input[53].metric, input[74].metric);
    sc_uint<1> comp_r_53_75 = COMP < Q > (input[53].metric, input[75].metric);
    sc_uint<1> comp_r_53_76 = COMP < Q > (input[53].metric, input[76].metric);
    sc_uint<1> comp_r_53_77 = COMP < Q > (input[53].metric, input[77].metric);
    sc_uint<1> comp_r_53_78 = COMP < Q > (input[53].metric, input[78].metric);
    sc_uint<1> comp_r_53_79 = COMP < Q > (input[53].metric, input[79].metric);
    sc_uint<1> comp_r_53_80 = COMP < Q > (input[53].metric, input[80].metric);
    sc_uint<1> comp_r_53_81 = COMP < Q > (input[53].metric, input[81].metric);
    sc_uint<1> comp_r_53_82 = COMP < Q > (input[53].metric, input[82].metric);
    sc_uint<1> comp_r_53_83 = COMP < Q > (input[53].metric, input[83].metric);
    sc_uint<1> comp_r_53_84 = COMP < Q > (input[53].metric, input[84].metric);
    sc_uint<1> comp_r_53_85 = COMP < Q > (input[53].metric, input[85].metric);
    sc_uint<1> comp_r_53_86 = COMP < Q > (input[53].metric, input[86].metric);
    sc_uint<1> comp_r_53_87 = COMP < Q > (input[53].metric, input[87].metric);
    sc_uint<1> comp_r_53_88 = COMP < Q > (input[53].metric, input[88].metric);
    sc_uint<1> comp_r_53_89 = COMP < Q > (input[53].metric, input[89].metric);
    sc_uint<1> comp_r_53_90 = COMP < Q > (input[53].metric, input[90].metric);
    sc_uint<1> comp_r_53_91 = COMP < Q > (input[53].metric, input[91].metric);
    sc_uint<1> comp_r_53_92 = COMP < Q > (input[53].metric, input[92].metric);
    sc_uint<1> comp_r_53_93 = COMP < Q > (input[53].metric, input[93].metric);
    sc_uint<1> comp_r_53_94 = COMP < Q > (input[53].metric, input[94].metric);
    sc_uint<1> comp_r_53_95 = COMP < Q > (input[53].metric, input[95].metric);
    sc_uint<1> comp_r_53_96 = COMP < Q > (input[53].metric, input[96].metric);
    sc_uint<1> comp_r_53_97 = COMP < Q > (input[53].metric, input[97].metric);
    sc_uint<1> comp_r_53_98 = COMP < Q > (input[53].metric, input[98].metric);
    sc_uint<1> comp_r_53_99 = COMP < Q > (input[53].metric, input[99].metric);
    sc_uint<1> comp_r_53_100 = COMP < Q > (input[53].metric, input[100].metric);
    sc_uint<1> comp_r_53_101 = COMP < Q > (input[53].metric, input[101].metric);
    sc_uint<1> comp_r_53_102 = COMP < Q > (input[53].metric, input[102].metric);
    sc_uint<1> comp_r_53_103 = COMP < Q > (input[53].metric, input[103].metric);
    sc_uint<1> comp_r_53_104 = COMP < Q > (input[53].metric, input[104].metric);
    sc_uint<1> comp_r_53_105 = COMP < Q > (input[53].metric, input[105].metric);
    sc_uint<1> comp_r_53_106 = COMP < Q > (input[53].metric, input[106].metric);
    sc_uint<1> comp_r_53_107 = COMP < Q > (input[53].metric, input[107].metric);
    sc_uint<1> comp_r_53_108 = COMP < Q > (input[53].metric, input[108].metric);
    sc_uint<1> comp_r_53_109 = COMP < Q > (input[53].metric, input[109].metric);
    sc_uint<1> comp_r_53_110 = COMP < Q > (input[53].metric, input[110].metric);
    sc_uint<1> comp_r_53_111 = COMP < Q > (input[53].metric, input[111].metric);
    sc_uint<1> comp_r_53_112 = COMP < Q > (input[53].metric, input[112].metric);
    sc_uint<1> comp_r_53_113 = COMP < Q > (input[53].metric, input[113].metric);
    sc_uint<1> comp_r_53_114 = COMP < Q > (input[53].metric, input[114].metric);
    sc_uint<1> comp_r_53_115 = COMP < Q > (input[53].metric, input[115].metric);
    sc_uint<1> comp_r_53_116 = COMP < Q > (input[53].metric, input[116].metric);
    sc_uint<1> comp_r_53_117 = COMP < Q > (input[53].metric, input[117].metric);
    sc_uint<1> comp_r_53_118 = COMP < Q > (input[53].metric, input[118].metric);
    sc_uint<1> comp_r_53_119 = COMP < Q > (input[53].metric, input[119].metric);
    sc_uint<1> comp_r_53_120 = COMP < Q > (input[53].metric, input[120].metric);
    sc_uint<1> comp_r_53_121 = COMP < Q > (input[53].metric, input[121].metric);
    sc_uint<1> comp_r_53_122 = COMP < Q > (input[53].metric, input[122].metric);
    sc_uint<1> comp_r_53_123 = COMP < Q > (input[53].metric, input[123].metric);
    sc_uint<1> comp_r_53_124 = COMP < Q > (input[53].metric, input[124].metric);
    sc_uint<1> comp_r_53_125 = COMP < Q > (input[53].metric, input[125].metric);
    sc_uint<1> comp_r_53_126 = COMP < Q > (input[53].metric, input[126].metric);

    sc_uint<1> comp_r_55_56 = COMP < Q > (input[55].metric, input[56].metric);
    sc_uint<1> comp_r_55_57 = COMP < Q > (input[55].metric, input[57].metric);
    sc_uint<1> comp_r_55_58 = COMP < Q > (input[55].metric, input[58].metric);
    sc_uint<1> comp_r_55_59 = COMP < Q > (input[55].metric, input[59].metric);
    sc_uint<1> comp_r_55_60 = COMP < Q > (input[55].metric, input[60].metric);
    sc_uint<1> comp_r_55_61 = COMP < Q > (input[55].metric, input[61].metric);
    sc_uint<1> comp_r_55_62 = COMP < Q > (input[55].metric, input[62].metric);
    sc_uint<1> comp_r_55_63 = COMP < Q > (input[55].metric, input[63].metric);
    sc_uint<1> comp_r_55_64 = COMP < Q > (input[55].metric, input[64].metric);
    sc_uint<1> comp_r_55_65 = COMP < Q > (input[55].metric, input[65].metric);
    sc_uint<1> comp_r_55_66 = COMP < Q > (input[55].metric, input[66].metric);
    sc_uint<1> comp_r_55_67 = COMP < Q > (input[55].metric, input[67].metric);
    sc_uint<1> comp_r_55_68 = COMP < Q > (input[55].metric, input[68].metric);
    sc_uint<1> comp_r_55_69 = COMP < Q > (input[55].metric, input[69].metric);
    sc_uint<1> comp_r_55_70 = COMP < Q > (input[55].metric, input[70].metric);
    sc_uint<1> comp_r_55_71 = COMP < Q > (input[55].metric, input[71].metric);
    sc_uint<1> comp_r_55_72 = COMP < Q > (input[55].metric, input[72].metric);
    sc_uint<1> comp_r_55_73 = COMP < Q > (input[55].metric, input[73].metric);
    sc_uint<1> comp_r_55_74 = COMP < Q > (input[55].metric, input[74].metric);
    sc_uint<1> comp_r_55_75 = COMP < Q > (input[55].metric, input[75].metric);
    sc_uint<1> comp_r_55_76 = COMP < Q > (input[55].metric, input[76].metric);
    sc_uint<1> comp_r_55_77 = COMP < Q > (input[55].metric, input[77].metric);
    sc_uint<1> comp_r_55_78 = COMP < Q > (input[55].metric, input[78].metric);
    sc_uint<1> comp_r_55_79 = COMP < Q > (input[55].metric, input[79].metric);
    sc_uint<1> comp_r_55_80 = COMP < Q > (input[55].metric, input[80].metric);
    sc_uint<1> comp_r_55_81 = COMP < Q > (input[55].metric, input[81].metric);
    sc_uint<1> comp_r_55_82 = COMP < Q > (input[55].metric, input[82].metric);
    sc_uint<1> comp_r_55_83 = COMP < Q > (input[55].metric, input[83].metric);
    sc_uint<1> comp_r_55_84 = COMP < Q > (input[55].metric, input[84].metric);
    sc_uint<1> comp_r_55_85 = COMP < Q > (input[55].metric, input[85].metric);
    sc_uint<1> comp_r_55_86 = COMP < Q > (input[55].metric, input[86].metric);
    sc_uint<1> comp_r_55_87 = COMP < Q > (input[55].metric, input[87].metric);
    sc_uint<1> comp_r_55_88 = COMP < Q > (input[55].metric, input[88].metric);
    sc_uint<1> comp_r_55_89 = COMP < Q > (input[55].metric, input[89].metric);
    sc_uint<1> comp_r_55_90 = COMP < Q > (input[55].metric, input[90].metric);
    sc_uint<1> comp_r_55_91 = COMP < Q > (input[55].metric, input[91].metric);
    sc_uint<1> comp_r_55_92 = COMP < Q > (input[55].metric, input[92].metric);
    sc_uint<1> comp_r_55_93 = COMP < Q > (input[55].metric, input[93].metric);
    sc_uint<1> comp_r_55_94 = COMP < Q > (input[55].metric, input[94].metric);
    sc_uint<1> comp_r_55_95 = COMP < Q > (input[55].metric, input[95].metric);
    sc_uint<1> comp_r_55_96 = COMP < Q > (input[55].metric, input[96].metric);
    sc_uint<1> comp_r_55_97 = COMP < Q > (input[55].metric, input[97].metric);
    sc_uint<1> comp_r_55_98 = COMP < Q > (input[55].metric, input[98].metric);
    sc_uint<1> comp_r_55_99 = COMP < Q > (input[55].metric, input[99].metric);
    sc_uint<1> comp_r_55_100 = COMP < Q > (input[55].metric, input[100].metric);
    sc_uint<1> comp_r_55_101 = COMP < Q > (input[55].metric, input[101].metric);
    sc_uint<1> comp_r_55_102 = COMP < Q > (input[55].metric, input[102].metric);
    sc_uint<1> comp_r_55_103 = COMP < Q > (input[55].metric, input[103].metric);
    sc_uint<1> comp_r_55_104 = COMP < Q > (input[55].metric, input[104].metric);
    sc_uint<1> comp_r_55_105 = COMP < Q > (input[55].metric, input[105].metric);
    sc_uint<1> comp_r_55_106 = COMP < Q > (input[55].metric, input[106].metric);
    sc_uint<1> comp_r_55_107 = COMP < Q > (input[55].metric, input[107].metric);
    sc_uint<1> comp_r_55_108 = COMP < Q > (input[55].metric, input[108].metric);
    sc_uint<1> comp_r_55_109 = COMP < Q > (input[55].metric, input[109].metric);
    sc_uint<1> comp_r_55_110 = COMP < Q > (input[55].metric, input[110].metric);
    sc_uint<1> comp_r_55_111 = COMP < Q > (input[55].metric, input[111].metric);
    sc_uint<1> comp_r_55_112 = COMP < Q > (input[55].metric, input[112].metric);
    sc_uint<1> comp_r_55_113 = COMP < Q > (input[55].metric, input[113].metric);
    sc_uint<1> comp_r_55_114 = COMP < Q > (input[55].metric, input[114].metric);
    sc_uint<1> comp_r_55_115 = COMP < Q > (input[55].metric, input[115].metric);
    sc_uint<1> comp_r_55_116 = COMP < Q > (input[55].metric, input[116].metric);
    sc_uint<1> comp_r_55_117 = COMP < Q > (input[55].metric, input[117].metric);
    sc_uint<1> comp_r_55_118 = COMP < Q > (input[55].metric, input[118].metric);
    sc_uint<1> comp_r_55_119 = COMP < Q > (input[55].metric, input[119].metric);
    sc_uint<1> comp_r_55_120 = COMP < Q > (input[55].metric, input[120].metric);
    sc_uint<1> comp_r_55_121 = COMP < Q > (input[55].metric, input[121].metric);
    sc_uint<1> comp_r_55_122 = COMP < Q > (input[55].metric, input[122].metric);
    sc_uint<1> comp_r_55_123 = COMP < Q > (input[55].metric, input[123].metric);
    sc_uint<1> comp_r_55_124 = COMP < Q > (input[55].metric, input[124].metric);
    sc_uint<1> comp_r_55_125 = COMP < Q > (input[55].metric, input[125].metric);
    sc_uint<1> comp_r_55_126 = COMP < Q > (input[55].metric, input[126].metric);

    sc_uint<1> comp_r_57_58 = COMP < Q > (input[57].metric, input[58].metric);
    sc_uint<1> comp_r_57_59 = COMP < Q > (input[57].metric, input[59].metric);
    sc_uint<1> comp_r_57_60 = COMP < Q > (input[57].metric, input[60].metric);
    sc_uint<1> comp_r_57_61 = COMP < Q > (input[57].metric, input[61].metric);
    sc_uint<1> comp_r_57_62 = COMP < Q > (input[57].metric, input[62].metric);
    sc_uint<1> comp_r_57_63 = COMP < Q > (input[57].metric, input[63].metric);
    sc_uint<1> comp_r_57_64 = COMP < Q > (input[57].metric, input[64].metric);
    sc_uint<1> comp_r_57_65 = COMP < Q > (input[57].metric, input[65].metric);
    sc_uint<1> comp_r_57_66 = COMP < Q > (input[57].metric, input[66].metric);
    sc_uint<1> comp_r_57_67 = COMP < Q > (input[57].metric, input[67].metric);
    sc_uint<1> comp_r_57_68 = COMP < Q > (input[57].metric, input[68].metric);
    sc_uint<1> comp_r_57_69 = COMP < Q > (input[57].metric, input[69].metric);
    sc_uint<1> comp_r_57_70 = COMP < Q > (input[57].metric, input[70].metric);
    sc_uint<1> comp_r_57_71 = COMP < Q > (input[57].metric, input[71].metric);
    sc_uint<1> comp_r_57_72 = COMP < Q > (input[57].metric, input[72].metric);
    sc_uint<1> comp_r_57_73 = COMP < Q > (input[57].metric, input[73].metric);
    sc_uint<1> comp_r_57_74 = COMP < Q > (input[57].metric, input[74].metric);
    sc_uint<1> comp_r_57_75 = COMP < Q > (input[57].metric, input[75].metric);
    sc_uint<1> comp_r_57_76 = COMP < Q > (input[57].metric, input[76].metric);
    sc_uint<1> comp_r_57_77 = COMP < Q > (input[57].metric, input[77].metric);
    sc_uint<1> comp_r_57_78 = COMP < Q > (input[57].metric, input[78].metric);
    sc_uint<1> comp_r_57_79 = COMP < Q > (input[57].metric, input[79].metric);
    sc_uint<1> comp_r_57_80 = COMP < Q > (input[57].metric, input[80].metric);
    sc_uint<1> comp_r_57_81 = COMP < Q > (input[57].metric, input[81].metric);
    sc_uint<1> comp_r_57_82 = COMP < Q > (input[57].metric, input[82].metric);
    sc_uint<1> comp_r_57_83 = COMP < Q > (input[57].metric, input[83].metric);
    sc_uint<1> comp_r_57_84 = COMP < Q > (input[57].metric, input[84].metric);
    sc_uint<1> comp_r_57_85 = COMP < Q > (input[57].metric, input[85].metric);
    sc_uint<1> comp_r_57_86 = COMP < Q > (input[57].metric, input[86].metric);
    sc_uint<1> comp_r_57_87 = COMP < Q > (input[57].metric, input[87].metric);
    sc_uint<1> comp_r_57_88 = COMP < Q > (input[57].metric, input[88].metric);
    sc_uint<1> comp_r_57_89 = COMP < Q > (input[57].metric, input[89].metric);
    sc_uint<1> comp_r_57_90 = COMP < Q > (input[57].metric, input[90].metric);
    sc_uint<1> comp_r_57_91 = COMP < Q > (input[57].metric, input[91].metric);
    sc_uint<1> comp_r_57_92 = COMP < Q > (input[57].metric, input[92].metric);
    sc_uint<1> comp_r_57_93 = COMP < Q > (input[57].metric, input[93].metric);
    sc_uint<1> comp_r_57_94 = COMP < Q > (input[57].metric, input[94].metric);
    sc_uint<1> comp_r_57_95 = COMP < Q > (input[57].metric, input[95].metric);
    sc_uint<1> comp_r_57_96 = COMP < Q > (input[57].metric, input[96].metric);
    sc_uint<1> comp_r_57_97 = COMP < Q > (input[57].metric, input[97].metric);
    sc_uint<1> comp_r_57_98 = COMP < Q > (input[57].metric, input[98].metric);
    sc_uint<1> comp_r_57_99 = COMP < Q > (input[57].metric, input[99].metric);
    sc_uint<1> comp_r_57_100 = COMP < Q > (input[57].metric, input[100].metric);
    sc_uint<1> comp_r_57_101 = COMP < Q > (input[57].metric, input[101].metric);
    sc_uint<1> comp_r_57_102 = COMP < Q > (input[57].metric, input[102].metric);
    sc_uint<1> comp_r_57_103 = COMP < Q > (input[57].metric, input[103].metric);
    sc_uint<1> comp_r_57_104 = COMP < Q > (input[57].metric, input[104].metric);
    sc_uint<1> comp_r_57_105 = COMP < Q > (input[57].metric, input[105].metric);
    sc_uint<1> comp_r_57_106 = COMP < Q > (input[57].metric, input[106].metric);
    sc_uint<1> comp_r_57_107 = COMP < Q > (input[57].metric, input[107].metric);
    sc_uint<1> comp_r_57_108 = COMP < Q > (input[57].metric, input[108].metric);
    sc_uint<1> comp_r_57_109 = COMP < Q > (input[57].metric, input[109].metric);
    sc_uint<1> comp_r_57_110 = COMP < Q > (input[57].metric, input[110].metric);
    sc_uint<1> comp_r_57_111 = COMP < Q > (input[57].metric, input[111].metric);
    sc_uint<1> comp_r_57_112 = COMP < Q > (input[57].metric, input[112].metric);
    sc_uint<1> comp_r_57_113 = COMP < Q > (input[57].metric, input[113].metric);
    sc_uint<1> comp_r_57_114 = COMP < Q > (input[57].metric, input[114].metric);
    sc_uint<1> comp_r_57_115 = COMP < Q > (input[57].metric, input[115].metric);
    sc_uint<1> comp_r_57_116 = COMP < Q > (input[57].metric, input[116].metric);
    sc_uint<1> comp_r_57_117 = COMP < Q > (input[57].metric, input[117].metric);
    sc_uint<1> comp_r_57_118 = COMP < Q > (input[57].metric, input[118].metric);
    sc_uint<1> comp_r_57_119 = COMP < Q > (input[57].metric, input[119].metric);
    sc_uint<1> comp_r_57_120 = COMP < Q > (input[57].metric, input[120].metric);
    sc_uint<1> comp_r_57_121 = COMP < Q > (input[57].metric, input[121].metric);
    sc_uint<1> comp_r_57_122 = COMP < Q > (input[57].metric, input[122].metric);
    sc_uint<1> comp_r_57_123 = COMP < Q > (input[57].metric, input[123].metric);
    sc_uint<1> comp_r_57_124 = COMP < Q > (input[57].metric, input[124].metric);
    sc_uint<1> comp_r_57_125 = COMP < Q > (input[57].metric, input[125].metric);
    sc_uint<1> comp_r_57_126 = COMP < Q > (input[57].metric, input[126].metric);

    sc_uint<1> comp_r_59_60 = COMP < Q > (input[59].metric, input[60].metric);
    sc_uint<1> comp_r_59_61 = COMP < Q > (input[59].metric, input[61].metric);
    sc_uint<1> comp_r_59_62 = COMP < Q > (input[59].metric, input[62].metric);
    sc_uint<1> comp_r_59_63 = COMP < Q > (input[59].metric, input[63].metric);
    sc_uint<1> comp_r_59_64 = COMP < Q > (input[59].metric, input[64].metric);
    sc_uint<1> comp_r_59_65 = COMP < Q > (input[59].metric, input[65].metric);
    sc_uint<1> comp_r_59_66 = COMP < Q > (input[59].metric, input[66].metric);
    sc_uint<1> comp_r_59_67 = COMP < Q > (input[59].metric, input[67].metric);
    sc_uint<1> comp_r_59_68 = COMP < Q > (input[59].metric, input[68].metric);
    sc_uint<1> comp_r_59_69 = COMP < Q > (input[59].metric, input[69].metric);
    sc_uint<1> comp_r_59_70 = COMP < Q > (input[59].metric, input[70].metric);
    sc_uint<1> comp_r_59_71 = COMP < Q > (input[59].metric, input[71].metric);
    sc_uint<1> comp_r_59_72 = COMP < Q > (input[59].metric, input[72].metric);
    sc_uint<1> comp_r_59_73 = COMP < Q > (input[59].metric, input[73].metric);
    sc_uint<1> comp_r_59_74 = COMP < Q > (input[59].metric, input[74].metric);
    sc_uint<1> comp_r_59_75 = COMP < Q > (input[59].metric, input[75].metric);
    sc_uint<1> comp_r_59_76 = COMP < Q > (input[59].metric, input[76].metric);
    sc_uint<1> comp_r_59_77 = COMP < Q > (input[59].metric, input[77].metric);
    sc_uint<1> comp_r_59_78 = COMP < Q > (input[59].metric, input[78].metric);
    sc_uint<1> comp_r_59_79 = COMP < Q > (input[59].metric, input[79].metric);
    sc_uint<1> comp_r_59_80 = COMP < Q > (input[59].metric, input[80].metric);
    sc_uint<1> comp_r_59_81 = COMP < Q > (input[59].metric, input[81].metric);
    sc_uint<1> comp_r_59_82 = COMP < Q > (input[59].metric, input[82].metric);
    sc_uint<1> comp_r_59_83 = COMP < Q > (input[59].metric, input[83].metric);
    sc_uint<1> comp_r_59_84 = COMP < Q > (input[59].metric, input[84].metric);
    sc_uint<1> comp_r_59_85 = COMP < Q > (input[59].metric, input[85].metric);
    sc_uint<1> comp_r_59_86 = COMP < Q > (input[59].metric, input[86].metric);
    sc_uint<1> comp_r_59_87 = COMP < Q > (input[59].metric, input[87].metric);
    sc_uint<1> comp_r_59_88 = COMP < Q > (input[59].metric, input[88].metric);
    sc_uint<1> comp_r_59_89 = COMP < Q > (input[59].metric, input[89].metric);
    sc_uint<1> comp_r_59_90 = COMP < Q > (input[59].metric, input[90].metric);
    sc_uint<1> comp_r_59_91 = COMP < Q > (input[59].metric, input[91].metric);
    sc_uint<1> comp_r_59_92 = COMP < Q > (input[59].metric, input[92].metric);
    sc_uint<1> comp_r_59_93 = COMP < Q > (input[59].metric, input[93].metric);
    sc_uint<1> comp_r_59_94 = COMP < Q > (input[59].metric, input[94].metric);
    sc_uint<1> comp_r_59_95 = COMP < Q > (input[59].metric, input[95].metric);
    sc_uint<1> comp_r_59_96 = COMP < Q > (input[59].metric, input[96].metric);
    sc_uint<1> comp_r_59_97 = COMP < Q > (input[59].metric, input[97].metric);
    sc_uint<1> comp_r_59_98 = COMP < Q > (input[59].metric, input[98].metric);
    sc_uint<1> comp_r_59_99 = COMP < Q > (input[59].metric, input[99].metric);
    sc_uint<1> comp_r_59_100 = COMP < Q > (input[59].metric, input[100].metric);
    sc_uint<1> comp_r_59_101 = COMP < Q > (input[59].metric, input[101].metric);
    sc_uint<1> comp_r_59_102 = COMP < Q > (input[59].metric, input[102].metric);
    sc_uint<1> comp_r_59_103 = COMP < Q > (input[59].metric, input[103].metric);
    sc_uint<1> comp_r_59_104 = COMP < Q > (input[59].metric, input[104].metric);
    sc_uint<1> comp_r_59_105 = COMP < Q > (input[59].metric, input[105].metric);
    sc_uint<1> comp_r_59_106 = COMP < Q > (input[59].metric, input[106].metric);
    sc_uint<1> comp_r_59_107 = COMP < Q > (input[59].metric, input[107].metric);
    sc_uint<1> comp_r_59_108 = COMP < Q > (input[59].metric, input[108].metric);
    sc_uint<1> comp_r_59_109 = COMP < Q > (input[59].metric, input[109].metric);
    sc_uint<1> comp_r_59_110 = COMP < Q > (input[59].metric, input[110].metric);
    sc_uint<1> comp_r_59_111 = COMP < Q > (input[59].metric, input[111].metric);
    sc_uint<1> comp_r_59_112 = COMP < Q > (input[59].metric, input[112].metric);
    sc_uint<1> comp_r_59_113 = COMP < Q > (input[59].metric, input[113].metric);
    sc_uint<1> comp_r_59_114 = COMP < Q > (input[59].metric, input[114].metric);
    sc_uint<1> comp_r_59_115 = COMP < Q > (input[59].metric, input[115].metric);
    sc_uint<1> comp_r_59_116 = COMP < Q > (input[59].metric, input[116].metric);
    sc_uint<1> comp_r_59_117 = COMP < Q > (input[59].metric, input[117].metric);
    sc_uint<1> comp_r_59_118 = COMP < Q > (input[59].metric, input[118].metric);
    sc_uint<1> comp_r_59_119 = COMP < Q > (input[59].metric, input[119].metric);
    sc_uint<1> comp_r_59_120 = COMP < Q > (input[59].metric, input[120].metric);
    sc_uint<1> comp_r_59_121 = COMP < Q > (input[59].metric, input[121].metric);
    sc_uint<1> comp_r_59_122 = COMP < Q > (input[59].metric, input[122].metric);
    sc_uint<1> comp_r_59_123 = COMP < Q > (input[59].metric, input[123].metric);
    sc_uint<1> comp_r_59_124 = COMP < Q > (input[59].metric, input[124].metric);
    sc_uint<1> comp_r_59_125 = COMP < Q > (input[59].metric, input[125].metric);
    sc_uint<1> comp_r_59_126 = COMP < Q > (input[59].metric, input[126].metric);

    sc_uint<1> comp_r_61_62 = COMP < Q > (input[61].metric, input[62].metric);
    sc_uint<1> comp_r_61_63 = COMP < Q > (input[61].metric, input[63].metric);
    sc_uint<1> comp_r_61_64 = COMP < Q > (input[61].metric, input[64].metric);
    sc_uint<1> comp_r_61_65 = COMP < Q > (input[61].metric, input[65].metric);
    sc_uint<1> comp_r_61_66 = COMP < Q > (input[61].metric, input[66].metric);
    sc_uint<1> comp_r_61_67 = COMP < Q > (input[61].metric, input[67].metric);
    sc_uint<1> comp_r_61_68 = COMP < Q > (input[61].metric, input[68].metric);
    sc_uint<1> comp_r_61_69 = COMP < Q > (input[61].metric, input[69].metric);
    sc_uint<1> comp_r_61_70 = COMP < Q > (input[61].metric, input[70].metric);
    sc_uint<1> comp_r_61_71 = COMP < Q > (input[61].metric, input[71].metric);
    sc_uint<1> comp_r_61_72 = COMP < Q > (input[61].metric, input[72].metric);
    sc_uint<1> comp_r_61_73 = COMP < Q > (input[61].metric, input[73].metric);
    sc_uint<1> comp_r_61_74 = COMP < Q > (input[61].metric, input[74].metric);
    sc_uint<1> comp_r_61_75 = COMP < Q > (input[61].metric, input[75].metric);
    sc_uint<1> comp_r_61_76 = COMP < Q > (input[61].metric, input[76].metric);
    sc_uint<1> comp_r_61_77 = COMP < Q > (input[61].metric, input[77].metric);
    sc_uint<1> comp_r_61_78 = COMP < Q > (input[61].metric, input[78].metric);
    sc_uint<1> comp_r_61_79 = COMP < Q > (input[61].metric, input[79].metric);
    sc_uint<1> comp_r_61_80 = COMP < Q > (input[61].metric, input[80].metric);
    sc_uint<1> comp_r_61_81 = COMP < Q > (input[61].metric, input[81].metric);
    sc_uint<1> comp_r_61_82 = COMP < Q > (input[61].metric, input[82].metric);
    sc_uint<1> comp_r_61_83 = COMP < Q > (input[61].metric, input[83].metric);
    sc_uint<1> comp_r_61_84 = COMP < Q > (input[61].metric, input[84].metric);
    sc_uint<1> comp_r_61_85 = COMP < Q > (input[61].metric, input[85].metric);
    sc_uint<1> comp_r_61_86 = COMP < Q > (input[61].metric, input[86].metric);
    sc_uint<1> comp_r_61_87 = COMP < Q > (input[61].metric, input[87].metric);
    sc_uint<1> comp_r_61_88 = COMP < Q > (input[61].metric, input[88].metric);
    sc_uint<1> comp_r_61_89 = COMP < Q > (input[61].metric, input[89].metric);
    sc_uint<1> comp_r_61_90 = COMP < Q > (input[61].metric, input[90].metric);
    sc_uint<1> comp_r_61_91 = COMP < Q > (input[61].metric, input[91].metric);
    sc_uint<1> comp_r_61_92 = COMP < Q > (input[61].metric, input[92].metric);
    sc_uint<1> comp_r_61_93 = COMP < Q > (input[61].metric, input[93].metric);
    sc_uint<1> comp_r_61_94 = COMP < Q > (input[61].metric, input[94].metric);
    sc_uint<1> comp_r_61_95 = COMP < Q > (input[61].metric, input[95].metric);
    sc_uint<1> comp_r_61_96 = COMP < Q > (input[61].metric, input[96].metric);
    sc_uint<1> comp_r_61_97 = COMP < Q > (input[61].metric, input[97].metric);
    sc_uint<1> comp_r_61_98 = COMP < Q > (input[61].metric, input[98].metric);
    sc_uint<1> comp_r_61_99 = COMP < Q > (input[61].metric, input[99].metric);
    sc_uint<1> comp_r_61_100 = COMP < Q > (input[61].metric, input[100].metric);
    sc_uint<1> comp_r_61_101 = COMP < Q > (input[61].metric, input[101].metric);
    sc_uint<1> comp_r_61_102 = COMP < Q > (input[61].metric, input[102].metric);
    sc_uint<1> comp_r_61_103 = COMP < Q > (input[61].metric, input[103].metric);
    sc_uint<1> comp_r_61_104 = COMP < Q > (input[61].metric, input[104].metric);
    sc_uint<1> comp_r_61_105 = COMP < Q > (input[61].metric, input[105].metric);
    sc_uint<1> comp_r_61_106 = COMP < Q > (input[61].metric, input[106].metric);
    sc_uint<1> comp_r_61_107 = COMP < Q > (input[61].metric, input[107].metric);
    sc_uint<1> comp_r_61_108 = COMP < Q > (input[61].metric, input[108].metric);
    sc_uint<1> comp_r_61_109 = COMP < Q > (input[61].metric, input[109].metric);
    sc_uint<1> comp_r_61_110 = COMP < Q > (input[61].metric, input[110].metric);
    sc_uint<1> comp_r_61_111 = COMP < Q > (input[61].metric, input[111].metric);
    sc_uint<1> comp_r_61_112 = COMP < Q > (input[61].metric, input[112].metric);
    sc_uint<1> comp_r_61_113 = COMP < Q > (input[61].metric, input[113].metric);
    sc_uint<1> comp_r_61_114 = COMP < Q > (input[61].metric, input[114].metric);
    sc_uint<1> comp_r_61_115 = COMP < Q > (input[61].metric, input[115].metric);
    sc_uint<1> comp_r_61_116 = COMP < Q > (input[61].metric, input[116].metric);
    sc_uint<1> comp_r_61_117 = COMP < Q > (input[61].metric, input[117].metric);
    sc_uint<1> comp_r_61_118 = COMP < Q > (input[61].metric, input[118].metric);
    sc_uint<1> comp_r_61_119 = COMP < Q > (input[61].metric, input[119].metric);
    sc_uint<1> comp_r_61_120 = COMP < Q > (input[61].metric, input[120].metric);
    sc_uint<1> comp_r_61_121 = COMP < Q > (input[61].metric, input[121].metric);
    sc_uint<1> comp_r_61_122 = COMP < Q > (input[61].metric, input[122].metric);
    sc_uint<1> comp_r_61_123 = COMP < Q > (input[61].metric, input[123].metric);
    sc_uint<1> comp_r_61_124 = COMP < Q > (input[61].metric, input[124].metric);
    sc_uint<1> comp_r_61_125 = COMP < Q > (input[61].metric, input[125].metric);
    sc_uint<1> comp_r_61_126 = COMP < Q > (input[61].metric, input[126].metric);

    sc_uint<1> comp_r_63_64 = COMP < Q > (input[63].metric, input[64].metric);
    sc_uint<1> comp_r_63_65 = COMP < Q > (input[63].metric, input[65].metric);
    sc_uint<1> comp_r_63_66 = COMP < Q > (input[63].metric, input[66].metric);
    sc_uint<1> comp_r_63_67 = COMP < Q > (input[63].metric, input[67].metric);
    sc_uint<1> comp_r_63_68 = COMP < Q > (input[63].metric, input[68].metric);
    sc_uint<1> comp_r_63_69 = COMP < Q > (input[63].metric, input[69].metric);
    sc_uint<1> comp_r_63_70 = COMP < Q > (input[63].metric, input[70].metric);
    sc_uint<1> comp_r_63_71 = COMP < Q > (input[63].metric, input[71].metric);
    sc_uint<1> comp_r_63_72 = COMP < Q > (input[63].metric, input[72].metric);
    sc_uint<1> comp_r_63_73 = COMP < Q > (input[63].metric, input[73].metric);
    sc_uint<1> comp_r_63_74 = COMP < Q > (input[63].metric, input[74].metric);
    sc_uint<1> comp_r_63_75 = COMP < Q > (input[63].metric, input[75].metric);
    sc_uint<1> comp_r_63_76 = COMP < Q > (input[63].metric, input[76].metric);
    sc_uint<1> comp_r_63_77 = COMP < Q > (input[63].metric, input[77].metric);
    sc_uint<1> comp_r_63_78 = COMP < Q > (input[63].metric, input[78].metric);
    sc_uint<1> comp_r_63_79 = COMP < Q > (input[63].metric, input[79].metric);
    sc_uint<1> comp_r_63_80 = COMP < Q > (input[63].metric, input[80].metric);
    sc_uint<1> comp_r_63_81 = COMP < Q > (input[63].metric, input[81].metric);
    sc_uint<1> comp_r_63_82 = COMP < Q > (input[63].metric, input[82].metric);
    sc_uint<1> comp_r_63_83 = COMP < Q > (input[63].metric, input[83].metric);
    sc_uint<1> comp_r_63_84 = COMP < Q > (input[63].metric, input[84].metric);
    sc_uint<1> comp_r_63_85 = COMP < Q > (input[63].metric, input[85].metric);
    sc_uint<1> comp_r_63_86 = COMP < Q > (input[63].metric, input[86].metric);
    sc_uint<1> comp_r_63_87 = COMP < Q > (input[63].metric, input[87].metric);
    sc_uint<1> comp_r_63_88 = COMP < Q > (input[63].metric, input[88].metric);
    sc_uint<1> comp_r_63_89 = COMP < Q > (input[63].metric, input[89].metric);
    sc_uint<1> comp_r_63_90 = COMP < Q > (input[63].metric, input[90].metric);
    sc_uint<1> comp_r_63_91 = COMP < Q > (input[63].metric, input[91].metric);
    sc_uint<1> comp_r_63_92 = COMP < Q > (input[63].metric, input[92].metric);
    sc_uint<1> comp_r_63_93 = COMP < Q > (input[63].metric, input[93].metric);
    sc_uint<1> comp_r_63_94 = COMP < Q > (input[63].metric, input[94].metric);
    sc_uint<1> comp_r_63_95 = COMP < Q > (input[63].metric, input[95].metric);
    sc_uint<1> comp_r_63_96 = COMP < Q > (input[63].metric, input[96].metric);
    sc_uint<1> comp_r_63_97 = COMP < Q > (input[63].metric, input[97].metric);
    sc_uint<1> comp_r_63_98 = COMP < Q > (input[63].metric, input[98].metric);
    sc_uint<1> comp_r_63_99 = COMP < Q > (input[63].metric, input[99].metric);
    sc_uint<1> comp_r_63_100 = COMP < Q > (input[63].metric, input[100].metric);
    sc_uint<1> comp_r_63_101 = COMP < Q > (input[63].metric, input[101].metric);
    sc_uint<1> comp_r_63_102 = COMP < Q > (input[63].metric, input[102].metric);
    sc_uint<1> comp_r_63_103 = COMP < Q > (input[63].metric, input[103].metric);
    sc_uint<1> comp_r_63_104 = COMP < Q > (input[63].metric, input[104].metric);
    sc_uint<1> comp_r_63_105 = COMP < Q > (input[63].metric, input[105].metric);
    sc_uint<1> comp_r_63_106 = COMP < Q > (input[63].metric, input[106].metric);
    sc_uint<1> comp_r_63_107 = COMP < Q > (input[63].metric, input[107].metric);
    sc_uint<1> comp_r_63_108 = COMP < Q > (input[63].metric, input[108].metric);
    sc_uint<1> comp_r_63_109 = COMP < Q > (input[63].metric, input[109].metric);
    sc_uint<1> comp_r_63_110 = COMP < Q > (input[63].metric, input[110].metric);
    sc_uint<1> comp_r_63_111 = COMP < Q > (input[63].metric, input[111].metric);
    sc_uint<1> comp_r_63_112 = COMP < Q > (input[63].metric, input[112].metric);
    sc_uint<1> comp_r_63_113 = COMP < Q > (input[63].metric, input[113].metric);
    sc_uint<1> comp_r_63_114 = COMP < Q > (input[63].metric, input[114].metric);
    sc_uint<1> comp_r_63_115 = COMP < Q > (input[63].metric, input[115].metric);
    sc_uint<1> comp_r_63_116 = COMP < Q > (input[63].metric, input[116].metric);
    sc_uint<1> comp_r_63_117 = COMP < Q > (input[63].metric, input[117].metric);
    sc_uint<1> comp_r_63_118 = COMP < Q > (input[63].metric, input[118].metric);
    sc_uint<1> comp_r_63_119 = COMP < Q > (input[63].metric, input[119].metric);
    sc_uint<1> comp_r_63_120 = COMP < Q > (input[63].metric, input[120].metric);
    sc_uint<1> comp_r_63_121 = COMP < Q > (input[63].metric, input[121].metric);
    sc_uint<1> comp_r_63_122 = COMP < Q > (input[63].metric, input[122].metric);
    sc_uint<1> comp_r_63_123 = COMP < Q > (input[63].metric, input[123].metric);
    sc_uint<1> comp_r_63_124 = COMP < Q > (input[63].metric, input[124].metric);
    sc_uint<1> comp_r_63_125 = COMP < Q > (input[63].metric, input[125].metric);
    sc_uint<1> comp_r_63_126 = COMP < Q > (input[63].metric, input[126].metric);

    sc_uint<1> comp_r_65_66 = COMP < Q > (input[65].metric, input[66].metric);
    sc_uint<1> comp_r_65_67 = COMP < Q > (input[65].metric, input[67].metric);
    sc_uint<1> comp_r_65_68 = COMP < Q > (input[65].metric, input[68].metric);
    sc_uint<1> comp_r_65_69 = COMP < Q > (input[65].metric, input[69].metric);
    sc_uint<1> comp_r_65_70 = COMP < Q > (input[65].metric, input[70].metric);
    sc_uint<1> comp_r_65_71 = COMP < Q > (input[65].metric, input[71].metric);
    sc_uint<1> comp_r_65_72 = COMP < Q > (input[65].metric, input[72].metric);
    sc_uint<1> comp_r_65_73 = COMP < Q > (input[65].metric, input[73].metric);
    sc_uint<1> comp_r_65_74 = COMP < Q > (input[65].metric, input[74].metric);
    sc_uint<1> comp_r_65_75 = COMP < Q > (input[65].metric, input[75].metric);
    sc_uint<1> comp_r_65_76 = COMP < Q > (input[65].metric, input[76].metric);
    sc_uint<1> comp_r_65_77 = COMP < Q > (input[65].metric, input[77].metric);
    sc_uint<1> comp_r_65_78 = COMP < Q > (input[65].metric, input[78].metric);
    sc_uint<1> comp_r_65_79 = COMP < Q > (input[65].metric, input[79].metric);
    sc_uint<1> comp_r_65_80 = COMP < Q > (input[65].metric, input[80].metric);
    sc_uint<1> comp_r_65_81 = COMP < Q > (input[65].metric, input[81].metric);
    sc_uint<1> comp_r_65_82 = COMP < Q > (input[65].metric, input[82].metric);
    sc_uint<1> comp_r_65_83 = COMP < Q > (input[65].metric, input[83].metric);
    sc_uint<1> comp_r_65_84 = COMP < Q > (input[65].metric, input[84].metric);
    sc_uint<1> comp_r_65_85 = COMP < Q > (input[65].metric, input[85].metric);
    sc_uint<1> comp_r_65_86 = COMP < Q > (input[65].metric, input[86].metric);
    sc_uint<1> comp_r_65_87 = COMP < Q > (input[65].metric, input[87].metric);
    sc_uint<1> comp_r_65_88 = COMP < Q > (input[65].metric, input[88].metric);
    sc_uint<1> comp_r_65_89 = COMP < Q > (input[65].metric, input[89].metric);
    sc_uint<1> comp_r_65_90 = COMP < Q > (input[65].metric, input[90].metric);
    sc_uint<1> comp_r_65_91 = COMP < Q > (input[65].metric, input[91].metric);
    sc_uint<1> comp_r_65_92 = COMP < Q > (input[65].metric, input[92].metric);
    sc_uint<1> comp_r_65_93 = COMP < Q > (input[65].metric, input[93].metric);
    sc_uint<1> comp_r_65_94 = COMP < Q > (input[65].metric, input[94].metric);
    sc_uint<1> comp_r_65_95 = COMP < Q > (input[65].metric, input[95].metric);
    sc_uint<1> comp_r_65_96 = COMP < Q > (input[65].metric, input[96].metric);
    sc_uint<1> comp_r_65_97 = COMP < Q > (input[65].metric, input[97].metric);
    sc_uint<1> comp_r_65_98 = COMP < Q > (input[65].metric, input[98].metric);
    sc_uint<1> comp_r_65_99 = COMP < Q > (input[65].metric, input[99].metric);
    sc_uint<1> comp_r_65_100 = COMP < Q > (input[65].metric, input[100].metric);
    sc_uint<1> comp_r_65_101 = COMP < Q > (input[65].metric, input[101].metric);
    sc_uint<1> comp_r_65_102 = COMP < Q > (input[65].metric, input[102].metric);
    sc_uint<1> comp_r_65_103 = COMP < Q > (input[65].metric, input[103].metric);
    sc_uint<1> comp_r_65_104 = COMP < Q > (input[65].metric, input[104].metric);
    sc_uint<1> comp_r_65_105 = COMP < Q > (input[65].metric, input[105].metric);
    sc_uint<1> comp_r_65_106 = COMP < Q > (input[65].metric, input[106].metric);
    sc_uint<1> comp_r_65_107 = COMP < Q > (input[65].metric, input[107].metric);
    sc_uint<1> comp_r_65_108 = COMP < Q > (input[65].metric, input[108].metric);
    sc_uint<1> comp_r_65_109 = COMP < Q > (input[65].metric, input[109].metric);
    sc_uint<1> comp_r_65_110 = COMP < Q > (input[65].metric, input[110].metric);
    sc_uint<1> comp_r_65_111 = COMP < Q > (input[65].metric, input[111].metric);
    sc_uint<1> comp_r_65_112 = COMP < Q > (input[65].metric, input[112].metric);
    sc_uint<1> comp_r_65_113 = COMP < Q > (input[65].metric, input[113].metric);
    sc_uint<1> comp_r_65_114 = COMP < Q > (input[65].metric, input[114].metric);
    sc_uint<1> comp_r_65_115 = COMP < Q > (input[65].metric, input[115].metric);
    sc_uint<1> comp_r_65_116 = COMP < Q > (input[65].metric, input[116].metric);
    sc_uint<1> comp_r_65_117 = COMP < Q > (input[65].metric, input[117].metric);
    sc_uint<1> comp_r_65_118 = COMP < Q > (input[65].metric, input[118].metric);
    sc_uint<1> comp_r_65_119 = COMP < Q > (input[65].metric, input[119].metric);
    sc_uint<1> comp_r_65_120 = COMP < Q > (input[65].metric, input[120].metric);
    sc_uint<1> comp_r_65_121 = COMP < Q > (input[65].metric, input[121].metric);
    sc_uint<1> comp_r_65_122 = COMP < Q > (input[65].metric, input[122].metric);
    sc_uint<1> comp_r_65_123 = COMP < Q > (input[65].metric, input[123].metric);
    sc_uint<1> comp_r_65_124 = COMP < Q > (input[65].metric, input[124].metric);
    sc_uint<1> comp_r_65_125 = COMP < Q > (input[65].metric, input[125].metric);
    sc_uint<1> comp_r_65_126 = COMP < Q > (input[65].metric, input[126].metric);

    sc_uint<1> comp_r_67_68 = COMP < Q > (input[67].metric, input[68].metric);
    sc_uint<1> comp_r_67_69 = COMP < Q > (input[67].metric, input[69].metric);
    sc_uint<1> comp_r_67_70 = COMP < Q > (input[67].metric, input[70].metric);
    sc_uint<1> comp_r_67_71 = COMP < Q > (input[67].metric, input[71].metric);
    sc_uint<1> comp_r_67_72 = COMP < Q > (input[67].metric, input[72].metric);
    sc_uint<1> comp_r_67_73 = COMP < Q > (input[67].metric, input[73].metric);
    sc_uint<1> comp_r_67_74 = COMP < Q > (input[67].metric, input[74].metric);
    sc_uint<1> comp_r_67_75 = COMP < Q > (input[67].metric, input[75].metric);
    sc_uint<1> comp_r_67_76 = COMP < Q > (input[67].metric, input[76].metric);
    sc_uint<1> comp_r_67_77 = COMP < Q > (input[67].metric, input[77].metric);
    sc_uint<1> comp_r_67_78 = COMP < Q > (input[67].metric, input[78].metric);
    sc_uint<1> comp_r_67_79 = COMP < Q > (input[67].metric, input[79].metric);
    sc_uint<1> comp_r_67_80 = COMP < Q > (input[67].metric, input[80].metric);
    sc_uint<1> comp_r_67_81 = COMP < Q > (input[67].metric, input[81].metric);
    sc_uint<1> comp_r_67_82 = COMP < Q > (input[67].metric, input[82].metric);
    sc_uint<1> comp_r_67_83 = COMP < Q > (input[67].metric, input[83].metric);
    sc_uint<1> comp_r_67_84 = COMP < Q > (input[67].metric, input[84].metric);
    sc_uint<1> comp_r_67_85 = COMP < Q > (input[67].metric, input[85].metric);
    sc_uint<1> comp_r_67_86 = COMP < Q > (input[67].metric, input[86].metric);
    sc_uint<1> comp_r_67_87 = COMP < Q > (input[67].metric, input[87].metric);
    sc_uint<1> comp_r_67_88 = COMP < Q > (input[67].metric, input[88].metric);
    sc_uint<1> comp_r_67_89 = COMP < Q > (input[67].metric, input[89].metric);
    sc_uint<1> comp_r_67_90 = COMP < Q > (input[67].metric, input[90].metric);
    sc_uint<1> comp_r_67_91 = COMP < Q > (input[67].metric, input[91].metric);
    sc_uint<1> comp_r_67_92 = COMP < Q > (input[67].metric, input[92].metric);
    sc_uint<1> comp_r_67_93 = COMP < Q > (input[67].metric, input[93].metric);
    sc_uint<1> comp_r_67_94 = COMP < Q > (input[67].metric, input[94].metric);
    sc_uint<1> comp_r_67_95 = COMP < Q > (input[67].metric, input[95].metric);
    sc_uint<1> comp_r_67_96 = COMP < Q > (input[67].metric, input[96].metric);
    sc_uint<1> comp_r_67_97 = COMP < Q > (input[67].metric, input[97].metric);
    sc_uint<1> comp_r_67_98 = COMP < Q > (input[67].metric, input[98].metric);
    sc_uint<1> comp_r_67_99 = COMP < Q > (input[67].metric, input[99].metric);
    sc_uint<1> comp_r_67_100 = COMP < Q > (input[67].metric, input[100].metric);
    sc_uint<1> comp_r_67_101 = COMP < Q > (input[67].metric, input[101].metric);
    sc_uint<1> comp_r_67_102 = COMP < Q > (input[67].metric, input[102].metric);
    sc_uint<1> comp_r_67_103 = COMP < Q > (input[67].metric, input[103].metric);
    sc_uint<1> comp_r_67_104 = COMP < Q > (input[67].metric, input[104].metric);
    sc_uint<1> comp_r_67_105 = COMP < Q > (input[67].metric, input[105].metric);
    sc_uint<1> comp_r_67_106 = COMP < Q > (input[67].metric, input[106].metric);
    sc_uint<1> comp_r_67_107 = COMP < Q > (input[67].metric, input[107].metric);
    sc_uint<1> comp_r_67_108 = COMP < Q > (input[67].metric, input[108].metric);
    sc_uint<1> comp_r_67_109 = COMP < Q > (input[67].metric, input[109].metric);
    sc_uint<1> comp_r_67_110 = COMP < Q > (input[67].metric, input[110].metric);
    sc_uint<1> comp_r_67_111 = COMP < Q > (input[67].metric, input[111].metric);
    sc_uint<1> comp_r_67_112 = COMP < Q > (input[67].metric, input[112].metric);
    sc_uint<1> comp_r_67_113 = COMP < Q > (input[67].metric, input[113].metric);
    sc_uint<1> comp_r_67_114 = COMP < Q > (input[67].metric, input[114].metric);
    sc_uint<1> comp_r_67_115 = COMP < Q > (input[67].metric, input[115].metric);
    sc_uint<1> comp_r_67_116 = COMP < Q > (input[67].metric, input[116].metric);
    sc_uint<1> comp_r_67_117 = COMP < Q > (input[67].metric, input[117].metric);
    sc_uint<1> comp_r_67_118 = COMP < Q > (input[67].metric, input[118].metric);
    sc_uint<1> comp_r_67_119 = COMP < Q > (input[67].metric, input[119].metric);
    sc_uint<1> comp_r_67_120 = COMP < Q > (input[67].metric, input[120].metric);
    sc_uint<1> comp_r_67_121 = COMP < Q > (input[67].metric, input[121].metric);
    sc_uint<1> comp_r_67_122 = COMP < Q > (input[67].metric, input[122].metric);
    sc_uint<1> comp_r_67_123 = COMP < Q > (input[67].metric, input[123].metric);
    sc_uint<1> comp_r_67_124 = COMP < Q > (input[67].metric, input[124].metric);
    sc_uint<1> comp_r_67_125 = COMP < Q > (input[67].metric, input[125].metric);
    sc_uint<1> comp_r_67_126 = COMP < Q > (input[67].metric, input[126].metric);

    sc_uint<1> comp_r_69_70 = COMP < Q > (input[69].metric, input[70].metric);
    sc_uint<1> comp_r_69_71 = COMP < Q > (input[69].metric, input[71].metric);
    sc_uint<1> comp_r_69_72 = COMP < Q > (input[69].metric, input[72].metric);
    sc_uint<1> comp_r_69_73 = COMP < Q > (input[69].metric, input[73].metric);
    sc_uint<1> comp_r_69_74 = COMP < Q > (input[69].metric, input[74].metric);
    sc_uint<1> comp_r_69_75 = COMP < Q > (input[69].metric, input[75].metric);
    sc_uint<1> comp_r_69_76 = COMP < Q > (input[69].metric, input[76].metric);
    sc_uint<1> comp_r_69_77 = COMP < Q > (input[69].metric, input[77].metric);
    sc_uint<1> comp_r_69_78 = COMP < Q > (input[69].metric, input[78].metric);
    sc_uint<1> comp_r_69_79 = COMP < Q > (input[69].metric, input[79].metric);
    sc_uint<1> comp_r_69_80 = COMP < Q > (input[69].metric, input[80].metric);
    sc_uint<1> comp_r_69_81 = COMP < Q > (input[69].metric, input[81].metric);
    sc_uint<1> comp_r_69_82 = COMP < Q > (input[69].metric, input[82].metric);
    sc_uint<1> comp_r_69_83 = COMP < Q > (input[69].metric, input[83].metric);
    sc_uint<1> comp_r_69_84 = COMP < Q > (input[69].metric, input[84].metric);
    sc_uint<1> comp_r_69_85 = COMP < Q > (input[69].metric, input[85].metric);
    sc_uint<1> comp_r_69_86 = COMP < Q > (input[69].metric, input[86].metric);
    sc_uint<1> comp_r_69_87 = COMP < Q > (input[69].metric, input[87].metric);
    sc_uint<1> comp_r_69_88 = COMP < Q > (input[69].metric, input[88].metric);
    sc_uint<1> comp_r_69_89 = COMP < Q > (input[69].metric, input[89].metric);
    sc_uint<1> comp_r_69_90 = COMP < Q > (input[69].metric, input[90].metric);
    sc_uint<1> comp_r_69_91 = COMP < Q > (input[69].metric, input[91].metric);
    sc_uint<1> comp_r_69_92 = COMP < Q > (input[69].metric, input[92].metric);
    sc_uint<1> comp_r_69_93 = COMP < Q > (input[69].metric, input[93].metric);
    sc_uint<1> comp_r_69_94 = COMP < Q > (input[69].metric, input[94].metric);
    sc_uint<1> comp_r_69_95 = COMP < Q > (input[69].metric, input[95].metric);
    sc_uint<1> comp_r_69_96 = COMP < Q > (input[69].metric, input[96].metric);
    sc_uint<1> comp_r_69_97 = COMP < Q > (input[69].metric, input[97].metric);
    sc_uint<1> comp_r_69_98 = COMP < Q > (input[69].metric, input[98].metric);
    sc_uint<1> comp_r_69_99 = COMP < Q > (input[69].metric, input[99].metric);
    sc_uint<1> comp_r_69_100 = COMP < Q > (input[69].metric, input[100].metric);
    sc_uint<1> comp_r_69_101 = COMP < Q > (input[69].metric, input[101].metric);
    sc_uint<1> comp_r_69_102 = COMP < Q > (input[69].metric, input[102].metric);
    sc_uint<1> comp_r_69_103 = COMP < Q > (input[69].metric, input[103].metric);
    sc_uint<1> comp_r_69_104 = COMP < Q > (input[69].metric, input[104].metric);
    sc_uint<1> comp_r_69_105 = COMP < Q > (input[69].metric, input[105].metric);
    sc_uint<1> comp_r_69_106 = COMP < Q > (input[69].metric, input[106].metric);
    sc_uint<1> comp_r_69_107 = COMP < Q > (input[69].metric, input[107].metric);
    sc_uint<1> comp_r_69_108 = COMP < Q > (input[69].metric, input[108].metric);
    sc_uint<1> comp_r_69_109 = COMP < Q > (input[69].metric, input[109].metric);
    sc_uint<1> comp_r_69_110 = COMP < Q > (input[69].metric, input[110].metric);
    sc_uint<1> comp_r_69_111 = COMP < Q > (input[69].metric, input[111].metric);
    sc_uint<1> comp_r_69_112 = COMP < Q > (input[69].metric, input[112].metric);
    sc_uint<1> comp_r_69_113 = COMP < Q > (input[69].metric, input[113].metric);
    sc_uint<1> comp_r_69_114 = COMP < Q > (input[69].metric, input[114].metric);
    sc_uint<1> comp_r_69_115 = COMP < Q > (input[69].metric, input[115].metric);
    sc_uint<1> comp_r_69_116 = COMP < Q > (input[69].metric, input[116].metric);
    sc_uint<1> comp_r_69_117 = COMP < Q > (input[69].metric, input[117].metric);
    sc_uint<1> comp_r_69_118 = COMP < Q > (input[69].metric, input[118].metric);
    sc_uint<1> comp_r_69_119 = COMP < Q > (input[69].metric, input[119].metric);
    sc_uint<1> comp_r_69_120 = COMP < Q > (input[69].metric, input[120].metric);
    sc_uint<1> comp_r_69_121 = COMP < Q > (input[69].metric, input[121].metric);
    sc_uint<1> comp_r_69_122 = COMP < Q > (input[69].metric, input[122].metric);
    sc_uint<1> comp_r_69_123 = COMP < Q > (input[69].metric, input[123].metric);
    sc_uint<1> comp_r_69_124 = COMP < Q > (input[69].metric, input[124].metric);
    sc_uint<1> comp_r_69_125 = COMP < Q > (input[69].metric, input[125].metric);
    sc_uint<1> comp_r_69_126 = COMP < Q > (input[69].metric, input[126].metric);

    sc_uint<1> comp_r_71_72 = COMP < Q > (input[71].metric, input[72].metric);
    sc_uint<1> comp_r_71_73 = COMP < Q > (input[71].metric, input[73].metric);
    sc_uint<1> comp_r_71_74 = COMP < Q > (input[71].metric, input[74].metric);
    sc_uint<1> comp_r_71_75 = COMP < Q > (input[71].metric, input[75].metric);
    sc_uint<1> comp_r_71_76 = COMP < Q > (input[71].metric, input[76].metric);
    sc_uint<1> comp_r_71_77 = COMP < Q > (input[71].metric, input[77].metric);
    sc_uint<1> comp_r_71_78 = COMP < Q > (input[71].metric, input[78].metric);
    sc_uint<1> comp_r_71_79 = COMP < Q > (input[71].metric, input[79].metric);
    sc_uint<1> comp_r_71_80 = COMP < Q > (input[71].metric, input[80].metric);
    sc_uint<1> comp_r_71_81 = COMP < Q > (input[71].metric, input[81].metric);
    sc_uint<1> comp_r_71_82 = COMP < Q > (input[71].metric, input[82].metric);
    sc_uint<1> comp_r_71_83 = COMP < Q > (input[71].metric, input[83].metric);
    sc_uint<1> comp_r_71_84 = COMP < Q > (input[71].metric, input[84].metric);
    sc_uint<1> comp_r_71_85 = COMP < Q > (input[71].metric, input[85].metric);
    sc_uint<1> comp_r_71_86 = COMP < Q > (input[71].metric, input[86].metric);
    sc_uint<1> comp_r_71_87 = COMP < Q > (input[71].metric, input[87].metric);
    sc_uint<1> comp_r_71_88 = COMP < Q > (input[71].metric, input[88].metric);
    sc_uint<1> comp_r_71_89 = COMP < Q > (input[71].metric, input[89].metric);
    sc_uint<1> comp_r_71_90 = COMP < Q > (input[71].metric, input[90].metric);
    sc_uint<1> comp_r_71_91 = COMP < Q > (input[71].metric, input[91].metric);
    sc_uint<1> comp_r_71_92 = COMP < Q > (input[71].metric, input[92].metric);
    sc_uint<1> comp_r_71_93 = COMP < Q > (input[71].metric, input[93].metric);
    sc_uint<1> comp_r_71_94 = COMP < Q > (input[71].metric, input[94].metric);
    sc_uint<1> comp_r_71_95 = COMP < Q > (input[71].metric, input[95].metric);
    sc_uint<1> comp_r_71_96 = COMP < Q > (input[71].metric, input[96].metric);
    sc_uint<1> comp_r_71_97 = COMP < Q > (input[71].metric, input[97].metric);
    sc_uint<1> comp_r_71_98 = COMP < Q > (input[71].metric, input[98].metric);
    sc_uint<1> comp_r_71_99 = COMP < Q > (input[71].metric, input[99].metric);
    sc_uint<1> comp_r_71_100 = COMP < Q > (input[71].metric, input[100].metric);
    sc_uint<1> comp_r_71_101 = COMP < Q > (input[71].metric, input[101].metric);
    sc_uint<1> comp_r_71_102 = COMP < Q > (input[71].metric, input[102].metric);
    sc_uint<1> comp_r_71_103 = COMP < Q > (input[71].metric, input[103].metric);
    sc_uint<1> comp_r_71_104 = COMP < Q > (input[71].metric, input[104].metric);
    sc_uint<1> comp_r_71_105 = COMP < Q > (input[71].metric, input[105].metric);
    sc_uint<1> comp_r_71_106 = COMP < Q > (input[71].metric, input[106].metric);
    sc_uint<1> comp_r_71_107 = COMP < Q > (input[71].metric, input[107].metric);
    sc_uint<1> comp_r_71_108 = COMP < Q > (input[71].metric, input[108].metric);
    sc_uint<1> comp_r_71_109 = COMP < Q > (input[71].metric, input[109].metric);
    sc_uint<1> comp_r_71_110 = COMP < Q > (input[71].metric, input[110].metric);
    sc_uint<1> comp_r_71_111 = COMP < Q > (input[71].metric, input[111].metric);
    sc_uint<1> comp_r_71_112 = COMP < Q > (input[71].metric, input[112].metric);
    sc_uint<1> comp_r_71_113 = COMP < Q > (input[71].metric, input[113].metric);
    sc_uint<1> comp_r_71_114 = COMP < Q > (input[71].metric, input[114].metric);
    sc_uint<1> comp_r_71_115 = COMP < Q > (input[71].metric, input[115].metric);
    sc_uint<1> comp_r_71_116 = COMP < Q > (input[71].metric, input[116].metric);
    sc_uint<1> comp_r_71_117 = COMP < Q > (input[71].metric, input[117].metric);
    sc_uint<1> comp_r_71_118 = COMP < Q > (input[71].metric, input[118].metric);
    sc_uint<1> comp_r_71_119 = COMP < Q > (input[71].metric, input[119].metric);
    sc_uint<1> comp_r_71_120 = COMP < Q > (input[71].metric, input[120].metric);
    sc_uint<1> comp_r_71_121 = COMP < Q > (input[71].metric, input[121].metric);
    sc_uint<1> comp_r_71_122 = COMP < Q > (input[71].metric, input[122].metric);
    sc_uint<1> comp_r_71_123 = COMP < Q > (input[71].metric, input[123].metric);
    sc_uint<1> comp_r_71_124 = COMP < Q > (input[71].metric, input[124].metric);
    sc_uint<1> comp_r_71_125 = COMP < Q > (input[71].metric, input[125].metric);
    sc_uint<1> comp_r_71_126 = COMP < Q > (input[71].metric, input[126].metric);

    sc_uint<1> comp_r_73_74 = COMP < Q > (input[73].metric, input[74].metric);
    sc_uint<1> comp_r_73_75 = COMP < Q > (input[73].metric, input[75].metric);
    sc_uint<1> comp_r_73_76 = COMP < Q > (input[73].metric, input[76].metric);
    sc_uint<1> comp_r_73_77 = COMP < Q > (input[73].metric, input[77].metric);
    sc_uint<1> comp_r_73_78 = COMP < Q > (input[73].metric, input[78].metric);
    sc_uint<1> comp_r_73_79 = COMP < Q > (input[73].metric, input[79].metric);
    sc_uint<1> comp_r_73_80 = COMP < Q > (input[73].metric, input[80].metric);
    sc_uint<1> comp_r_73_81 = COMP < Q > (input[73].metric, input[81].metric);
    sc_uint<1> comp_r_73_82 = COMP < Q > (input[73].metric, input[82].metric);
    sc_uint<1> comp_r_73_83 = COMP < Q > (input[73].metric, input[83].metric);
    sc_uint<1> comp_r_73_84 = COMP < Q > (input[73].metric, input[84].metric);
    sc_uint<1> comp_r_73_85 = COMP < Q > (input[73].metric, input[85].metric);
    sc_uint<1> comp_r_73_86 = COMP < Q > (input[73].metric, input[86].metric);
    sc_uint<1> comp_r_73_87 = COMP < Q > (input[73].metric, input[87].metric);
    sc_uint<1> comp_r_73_88 = COMP < Q > (input[73].metric, input[88].metric);
    sc_uint<1> comp_r_73_89 = COMP < Q > (input[73].metric, input[89].metric);
    sc_uint<1> comp_r_73_90 = COMP < Q > (input[73].metric, input[90].metric);
    sc_uint<1> comp_r_73_91 = COMP < Q > (input[73].metric, input[91].metric);
    sc_uint<1> comp_r_73_92 = COMP < Q > (input[73].metric, input[92].metric);
    sc_uint<1> comp_r_73_93 = COMP < Q > (input[73].metric, input[93].metric);
    sc_uint<1> comp_r_73_94 = COMP < Q > (input[73].metric, input[94].metric);
    sc_uint<1> comp_r_73_95 = COMP < Q > (input[73].metric, input[95].metric);
    sc_uint<1> comp_r_73_96 = COMP < Q > (input[73].metric, input[96].metric);
    sc_uint<1> comp_r_73_97 = COMP < Q > (input[73].metric, input[97].metric);
    sc_uint<1> comp_r_73_98 = COMP < Q > (input[73].metric, input[98].metric);
    sc_uint<1> comp_r_73_99 = COMP < Q > (input[73].metric, input[99].metric);
    sc_uint<1> comp_r_73_100 = COMP < Q > (input[73].metric, input[100].metric);
    sc_uint<1> comp_r_73_101 = COMP < Q > (input[73].metric, input[101].metric);
    sc_uint<1> comp_r_73_102 = COMP < Q > (input[73].metric, input[102].metric);
    sc_uint<1> comp_r_73_103 = COMP < Q > (input[73].metric, input[103].metric);
    sc_uint<1> comp_r_73_104 = COMP < Q > (input[73].metric, input[104].metric);
    sc_uint<1> comp_r_73_105 = COMP < Q > (input[73].metric, input[105].metric);
    sc_uint<1> comp_r_73_106 = COMP < Q > (input[73].metric, input[106].metric);
    sc_uint<1> comp_r_73_107 = COMP < Q > (input[73].metric, input[107].metric);
    sc_uint<1> comp_r_73_108 = COMP < Q > (input[73].metric, input[108].metric);
    sc_uint<1> comp_r_73_109 = COMP < Q > (input[73].metric, input[109].metric);
    sc_uint<1> comp_r_73_110 = COMP < Q > (input[73].metric, input[110].metric);
    sc_uint<1> comp_r_73_111 = COMP < Q > (input[73].metric, input[111].metric);
    sc_uint<1> comp_r_73_112 = COMP < Q > (input[73].metric, input[112].metric);
    sc_uint<1> comp_r_73_113 = COMP < Q > (input[73].metric, input[113].metric);
    sc_uint<1> comp_r_73_114 = COMP < Q > (input[73].metric, input[114].metric);
    sc_uint<1> comp_r_73_115 = COMP < Q > (input[73].metric, input[115].metric);
    sc_uint<1> comp_r_73_116 = COMP < Q > (input[73].metric, input[116].metric);
    sc_uint<1> comp_r_73_117 = COMP < Q > (input[73].metric, input[117].metric);
    sc_uint<1> comp_r_73_118 = COMP < Q > (input[73].metric, input[118].metric);
    sc_uint<1> comp_r_73_119 = COMP < Q > (input[73].metric, input[119].metric);
    sc_uint<1> comp_r_73_120 = COMP < Q > (input[73].metric, input[120].metric);
    sc_uint<1> comp_r_73_121 = COMP < Q > (input[73].metric, input[121].metric);
    sc_uint<1> comp_r_73_122 = COMP < Q > (input[73].metric, input[122].metric);
    sc_uint<1> comp_r_73_123 = COMP < Q > (input[73].metric, input[123].metric);
    sc_uint<1> comp_r_73_124 = COMP < Q > (input[73].metric, input[124].metric);
    sc_uint<1> comp_r_73_125 = COMP < Q > (input[73].metric, input[125].metric);
    sc_uint<1> comp_r_73_126 = COMP < Q > (input[73].metric, input[126].metric);

    sc_uint<1> comp_r_75_76 = COMP < Q > (input[75].metric, input[76].metric);
    sc_uint<1> comp_r_75_77 = COMP < Q > (input[75].metric, input[77].metric);
    sc_uint<1> comp_r_75_78 = COMP < Q > (input[75].metric, input[78].metric);
    sc_uint<1> comp_r_75_79 = COMP < Q > (input[75].metric, input[79].metric);
    sc_uint<1> comp_r_75_80 = COMP < Q > (input[75].metric, input[80].metric);
    sc_uint<1> comp_r_75_81 = COMP < Q > (input[75].metric, input[81].metric);
    sc_uint<1> comp_r_75_82 = COMP < Q > (input[75].metric, input[82].metric);
    sc_uint<1> comp_r_75_83 = COMP < Q > (input[75].metric, input[83].metric);
    sc_uint<1> comp_r_75_84 = COMP < Q > (input[75].metric, input[84].metric);
    sc_uint<1> comp_r_75_85 = COMP < Q > (input[75].metric, input[85].metric);
    sc_uint<1> comp_r_75_86 = COMP < Q > (input[75].metric, input[86].metric);
    sc_uint<1> comp_r_75_87 = COMP < Q > (input[75].metric, input[87].metric);
    sc_uint<1> comp_r_75_88 = COMP < Q > (input[75].metric, input[88].metric);
    sc_uint<1> comp_r_75_89 = COMP < Q > (input[75].metric, input[89].metric);
    sc_uint<1> comp_r_75_90 = COMP < Q > (input[75].metric, input[90].metric);
    sc_uint<1> comp_r_75_91 = COMP < Q > (input[75].metric, input[91].metric);
    sc_uint<1> comp_r_75_92 = COMP < Q > (input[75].metric, input[92].metric);
    sc_uint<1> comp_r_75_93 = COMP < Q > (input[75].metric, input[93].metric);
    sc_uint<1> comp_r_75_94 = COMP < Q > (input[75].metric, input[94].metric);
    sc_uint<1> comp_r_75_95 = COMP < Q > (input[75].metric, input[95].metric);
    sc_uint<1> comp_r_75_96 = COMP < Q > (input[75].metric, input[96].metric);
    sc_uint<1> comp_r_75_97 = COMP < Q > (input[75].metric, input[97].metric);
    sc_uint<1> comp_r_75_98 = COMP < Q > (input[75].metric, input[98].metric);
    sc_uint<1> comp_r_75_99 = COMP < Q > (input[75].metric, input[99].metric);
    sc_uint<1> comp_r_75_100 = COMP < Q > (input[75].metric, input[100].metric);
    sc_uint<1> comp_r_75_101 = COMP < Q > (input[75].metric, input[101].metric);
    sc_uint<1> comp_r_75_102 = COMP < Q > (input[75].metric, input[102].metric);
    sc_uint<1> comp_r_75_103 = COMP < Q > (input[75].metric, input[103].metric);
    sc_uint<1> comp_r_75_104 = COMP < Q > (input[75].metric, input[104].metric);
    sc_uint<1> comp_r_75_105 = COMP < Q > (input[75].metric, input[105].metric);
    sc_uint<1> comp_r_75_106 = COMP < Q > (input[75].metric, input[106].metric);
    sc_uint<1> comp_r_75_107 = COMP < Q > (input[75].metric, input[107].metric);
    sc_uint<1> comp_r_75_108 = COMP < Q > (input[75].metric, input[108].metric);
    sc_uint<1> comp_r_75_109 = COMP < Q > (input[75].metric, input[109].metric);
    sc_uint<1> comp_r_75_110 = COMP < Q > (input[75].metric, input[110].metric);
    sc_uint<1> comp_r_75_111 = COMP < Q > (input[75].metric, input[111].metric);
    sc_uint<1> comp_r_75_112 = COMP < Q > (input[75].metric, input[112].metric);
    sc_uint<1> comp_r_75_113 = COMP < Q > (input[75].metric, input[113].metric);
    sc_uint<1> comp_r_75_114 = COMP < Q > (input[75].metric, input[114].metric);
    sc_uint<1> comp_r_75_115 = COMP < Q > (input[75].metric, input[115].metric);
    sc_uint<1> comp_r_75_116 = COMP < Q > (input[75].metric, input[116].metric);
    sc_uint<1> comp_r_75_117 = COMP < Q > (input[75].metric, input[117].metric);
    sc_uint<1> comp_r_75_118 = COMP < Q > (input[75].metric, input[118].metric);
    sc_uint<1> comp_r_75_119 = COMP < Q > (input[75].metric, input[119].metric);
    sc_uint<1> comp_r_75_120 = COMP < Q > (input[75].metric, input[120].metric);
    sc_uint<1> comp_r_75_121 = COMP < Q > (input[75].metric, input[121].metric);
    sc_uint<1> comp_r_75_122 = COMP < Q > (input[75].metric, input[122].metric);
    sc_uint<1> comp_r_75_123 = COMP < Q > (input[75].metric, input[123].metric);
    sc_uint<1> comp_r_75_124 = COMP < Q > (input[75].metric, input[124].metric);
    sc_uint<1> comp_r_75_125 = COMP < Q > (input[75].metric, input[125].metric);
    sc_uint<1> comp_r_75_126 = COMP < Q > (input[75].metric, input[126].metric);

    sc_uint<1> comp_r_77_78 = COMP < Q > (input[77].metric, input[78].metric);
    sc_uint<1> comp_r_77_79 = COMP < Q > (input[77].metric, input[79].metric);
    sc_uint<1> comp_r_77_80 = COMP < Q > (input[77].metric, input[80].metric);
    sc_uint<1> comp_r_77_81 = COMP < Q > (input[77].metric, input[81].metric);
    sc_uint<1> comp_r_77_82 = COMP < Q > (input[77].metric, input[82].metric);
    sc_uint<1> comp_r_77_83 = COMP < Q > (input[77].metric, input[83].metric);
    sc_uint<1> comp_r_77_84 = COMP < Q > (input[77].metric, input[84].metric);
    sc_uint<1> comp_r_77_85 = COMP < Q > (input[77].metric, input[85].metric);
    sc_uint<1> comp_r_77_86 = COMP < Q > (input[77].metric, input[86].metric);
    sc_uint<1> comp_r_77_87 = COMP < Q > (input[77].metric, input[87].metric);
    sc_uint<1> comp_r_77_88 = COMP < Q > (input[77].metric, input[88].metric);
    sc_uint<1> comp_r_77_89 = COMP < Q > (input[77].metric, input[89].metric);
    sc_uint<1> comp_r_77_90 = COMP < Q > (input[77].metric, input[90].metric);
    sc_uint<1> comp_r_77_91 = COMP < Q > (input[77].metric, input[91].metric);
    sc_uint<1> comp_r_77_92 = COMP < Q > (input[77].metric, input[92].metric);
    sc_uint<1> comp_r_77_93 = COMP < Q > (input[77].metric, input[93].metric);
    sc_uint<1> comp_r_77_94 = COMP < Q > (input[77].metric, input[94].metric);
    sc_uint<1> comp_r_77_95 = COMP < Q > (input[77].metric, input[95].metric);
    sc_uint<1> comp_r_77_96 = COMP < Q > (input[77].metric, input[96].metric);
    sc_uint<1> comp_r_77_97 = COMP < Q > (input[77].metric, input[97].metric);
    sc_uint<1> comp_r_77_98 = COMP < Q > (input[77].metric, input[98].metric);
    sc_uint<1> comp_r_77_99 = COMP < Q > (input[77].metric, input[99].metric);
    sc_uint<1> comp_r_77_100 = COMP < Q > (input[77].metric, input[100].metric);
    sc_uint<1> comp_r_77_101 = COMP < Q > (input[77].metric, input[101].metric);
    sc_uint<1> comp_r_77_102 = COMP < Q > (input[77].metric, input[102].metric);
    sc_uint<1> comp_r_77_103 = COMP < Q > (input[77].metric, input[103].metric);
    sc_uint<1> comp_r_77_104 = COMP < Q > (input[77].metric, input[104].metric);
    sc_uint<1> comp_r_77_105 = COMP < Q > (input[77].metric, input[105].metric);
    sc_uint<1> comp_r_77_106 = COMP < Q > (input[77].metric, input[106].metric);
    sc_uint<1> comp_r_77_107 = COMP < Q > (input[77].metric, input[107].metric);
    sc_uint<1> comp_r_77_108 = COMP < Q > (input[77].metric, input[108].metric);
    sc_uint<1> comp_r_77_109 = COMP < Q > (input[77].metric, input[109].metric);
    sc_uint<1> comp_r_77_110 = COMP < Q > (input[77].metric, input[110].metric);
    sc_uint<1> comp_r_77_111 = COMP < Q > (input[77].metric, input[111].metric);
    sc_uint<1> comp_r_77_112 = COMP < Q > (input[77].metric, input[112].metric);
    sc_uint<1> comp_r_77_113 = COMP < Q > (input[77].metric, input[113].metric);
    sc_uint<1> comp_r_77_114 = COMP < Q > (input[77].metric, input[114].metric);
    sc_uint<1> comp_r_77_115 = COMP < Q > (input[77].metric, input[115].metric);
    sc_uint<1> comp_r_77_116 = COMP < Q > (input[77].metric, input[116].metric);
    sc_uint<1> comp_r_77_117 = COMP < Q > (input[77].metric, input[117].metric);
    sc_uint<1> comp_r_77_118 = COMP < Q > (input[77].metric, input[118].metric);
    sc_uint<1> comp_r_77_119 = COMP < Q > (input[77].metric, input[119].metric);
    sc_uint<1> comp_r_77_120 = COMP < Q > (input[77].metric, input[120].metric);
    sc_uint<1> comp_r_77_121 = COMP < Q > (input[77].metric, input[121].metric);
    sc_uint<1> comp_r_77_122 = COMP < Q > (input[77].metric, input[122].metric);
    sc_uint<1> comp_r_77_123 = COMP < Q > (input[77].metric, input[123].metric);
    sc_uint<1> comp_r_77_124 = COMP < Q > (input[77].metric, input[124].metric);
    sc_uint<1> comp_r_77_125 = COMP < Q > (input[77].metric, input[125].metric);
    sc_uint<1> comp_r_77_126 = COMP < Q > (input[77].metric, input[126].metric);

    sc_uint<1> comp_r_79_80 = COMP < Q > (input[79].metric, input[80].metric);
    sc_uint<1> comp_r_79_81 = COMP < Q > (input[79].metric, input[81].metric);
    sc_uint<1> comp_r_79_82 = COMP < Q > (input[79].metric, input[82].metric);
    sc_uint<1> comp_r_79_83 = COMP < Q > (input[79].metric, input[83].metric);
    sc_uint<1> comp_r_79_84 = COMP < Q > (input[79].metric, input[84].metric);
    sc_uint<1> comp_r_79_85 = COMP < Q > (input[79].metric, input[85].metric);
    sc_uint<1> comp_r_79_86 = COMP < Q > (input[79].metric, input[86].metric);
    sc_uint<1> comp_r_79_87 = COMP < Q > (input[79].metric, input[87].metric);
    sc_uint<1> comp_r_79_88 = COMP < Q > (input[79].metric, input[88].metric);
    sc_uint<1> comp_r_79_89 = COMP < Q > (input[79].metric, input[89].metric);
    sc_uint<1> comp_r_79_90 = COMP < Q > (input[79].metric, input[90].metric);
    sc_uint<1> comp_r_79_91 = COMP < Q > (input[79].metric, input[91].metric);
    sc_uint<1> comp_r_79_92 = COMP < Q > (input[79].metric, input[92].metric);
    sc_uint<1> comp_r_79_93 = COMP < Q > (input[79].metric, input[93].metric);
    sc_uint<1> comp_r_79_94 = COMP < Q > (input[79].metric, input[94].metric);
    sc_uint<1> comp_r_79_95 = COMP < Q > (input[79].metric, input[95].metric);
    sc_uint<1> comp_r_79_96 = COMP < Q > (input[79].metric, input[96].metric);
    sc_uint<1> comp_r_79_97 = COMP < Q > (input[79].metric, input[97].metric);
    sc_uint<1> comp_r_79_98 = COMP < Q > (input[79].metric, input[98].metric);
    sc_uint<1> comp_r_79_99 = COMP < Q > (input[79].metric, input[99].metric);
    sc_uint<1> comp_r_79_100 = COMP < Q > (input[79].metric, input[100].metric);
    sc_uint<1> comp_r_79_101 = COMP < Q > (input[79].metric, input[101].metric);
    sc_uint<1> comp_r_79_102 = COMP < Q > (input[79].metric, input[102].metric);
    sc_uint<1> comp_r_79_103 = COMP < Q > (input[79].metric, input[103].metric);
    sc_uint<1> comp_r_79_104 = COMP < Q > (input[79].metric, input[104].metric);
    sc_uint<1> comp_r_79_105 = COMP < Q > (input[79].metric, input[105].metric);
    sc_uint<1> comp_r_79_106 = COMP < Q > (input[79].metric, input[106].metric);
    sc_uint<1> comp_r_79_107 = COMP < Q > (input[79].metric, input[107].metric);
    sc_uint<1> comp_r_79_108 = COMP < Q > (input[79].metric, input[108].metric);
    sc_uint<1> comp_r_79_109 = COMP < Q > (input[79].metric, input[109].metric);
    sc_uint<1> comp_r_79_110 = COMP < Q > (input[79].metric, input[110].metric);
    sc_uint<1> comp_r_79_111 = COMP < Q > (input[79].metric, input[111].metric);
    sc_uint<1> comp_r_79_112 = COMP < Q > (input[79].metric, input[112].metric);
    sc_uint<1> comp_r_79_113 = COMP < Q > (input[79].metric, input[113].metric);
    sc_uint<1> comp_r_79_114 = COMP < Q > (input[79].metric, input[114].metric);
    sc_uint<1> comp_r_79_115 = COMP < Q > (input[79].metric, input[115].metric);
    sc_uint<1> comp_r_79_116 = COMP < Q > (input[79].metric, input[116].metric);
    sc_uint<1> comp_r_79_117 = COMP < Q > (input[79].metric, input[117].metric);
    sc_uint<1> comp_r_79_118 = COMP < Q > (input[79].metric, input[118].metric);
    sc_uint<1> comp_r_79_119 = COMP < Q > (input[79].metric, input[119].metric);
    sc_uint<1> comp_r_79_120 = COMP < Q > (input[79].metric, input[120].metric);
    sc_uint<1> comp_r_79_121 = COMP < Q > (input[79].metric, input[121].metric);
    sc_uint<1> comp_r_79_122 = COMP < Q > (input[79].metric, input[122].metric);
    sc_uint<1> comp_r_79_123 = COMP < Q > (input[79].metric, input[123].metric);
    sc_uint<1> comp_r_79_124 = COMP < Q > (input[79].metric, input[124].metric);
    sc_uint<1> comp_r_79_125 = COMP < Q > (input[79].metric, input[125].metric);
    sc_uint<1> comp_r_79_126 = COMP < Q > (input[79].metric, input[126].metric);

    sc_uint<1> comp_r_81_82 = COMP < Q > (input[81].metric, input[82].metric);
    sc_uint<1> comp_r_81_83 = COMP < Q > (input[81].metric, input[83].metric);
    sc_uint<1> comp_r_81_84 = COMP < Q > (input[81].metric, input[84].metric);
    sc_uint<1> comp_r_81_85 = COMP < Q > (input[81].metric, input[85].metric);
    sc_uint<1> comp_r_81_86 = COMP < Q > (input[81].metric, input[86].metric);
    sc_uint<1> comp_r_81_87 = COMP < Q > (input[81].metric, input[87].metric);
    sc_uint<1> comp_r_81_88 = COMP < Q > (input[81].metric, input[88].metric);
    sc_uint<1> comp_r_81_89 = COMP < Q > (input[81].metric, input[89].metric);
    sc_uint<1> comp_r_81_90 = COMP < Q > (input[81].metric, input[90].metric);
    sc_uint<1> comp_r_81_91 = COMP < Q > (input[81].metric, input[91].metric);
    sc_uint<1> comp_r_81_92 = COMP < Q > (input[81].metric, input[92].metric);
    sc_uint<1> comp_r_81_93 = COMP < Q > (input[81].metric, input[93].metric);
    sc_uint<1> comp_r_81_94 = COMP < Q > (input[81].metric, input[94].metric);
    sc_uint<1> comp_r_81_95 = COMP < Q > (input[81].metric, input[95].metric);
    sc_uint<1> comp_r_81_96 = COMP < Q > (input[81].metric, input[96].metric);
    sc_uint<1> comp_r_81_97 = COMP < Q > (input[81].metric, input[97].metric);
    sc_uint<1> comp_r_81_98 = COMP < Q > (input[81].metric, input[98].metric);
    sc_uint<1> comp_r_81_99 = COMP < Q > (input[81].metric, input[99].metric);
    sc_uint<1> comp_r_81_100 = COMP < Q > (input[81].metric, input[100].metric);
    sc_uint<1> comp_r_81_101 = COMP < Q > (input[81].metric, input[101].metric);
    sc_uint<1> comp_r_81_102 = COMP < Q > (input[81].metric, input[102].metric);
    sc_uint<1> comp_r_81_103 = COMP < Q > (input[81].metric, input[103].metric);
    sc_uint<1> comp_r_81_104 = COMP < Q > (input[81].metric, input[104].metric);
    sc_uint<1> comp_r_81_105 = COMP < Q > (input[81].metric, input[105].metric);
    sc_uint<1> comp_r_81_106 = COMP < Q > (input[81].metric, input[106].metric);
    sc_uint<1> comp_r_81_107 = COMP < Q > (input[81].metric, input[107].metric);
    sc_uint<1> comp_r_81_108 = COMP < Q > (input[81].metric, input[108].metric);
    sc_uint<1> comp_r_81_109 = COMP < Q > (input[81].metric, input[109].metric);
    sc_uint<1> comp_r_81_110 = COMP < Q > (input[81].metric, input[110].metric);
    sc_uint<1> comp_r_81_111 = COMP < Q > (input[81].metric, input[111].metric);
    sc_uint<1> comp_r_81_112 = COMP < Q > (input[81].metric, input[112].metric);
    sc_uint<1> comp_r_81_113 = COMP < Q > (input[81].metric, input[113].metric);
    sc_uint<1> comp_r_81_114 = COMP < Q > (input[81].metric, input[114].metric);
    sc_uint<1> comp_r_81_115 = COMP < Q > (input[81].metric, input[115].metric);
    sc_uint<1> comp_r_81_116 = COMP < Q > (input[81].metric, input[116].metric);
    sc_uint<1> comp_r_81_117 = COMP < Q > (input[81].metric, input[117].metric);
    sc_uint<1> comp_r_81_118 = COMP < Q > (input[81].metric, input[118].metric);
    sc_uint<1> comp_r_81_119 = COMP < Q > (input[81].metric, input[119].metric);
    sc_uint<1> comp_r_81_120 = COMP < Q > (input[81].metric, input[120].metric);
    sc_uint<1> comp_r_81_121 = COMP < Q > (input[81].metric, input[121].metric);
    sc_uint<1> comp_r_81_122 = COMP < Q > (input[81].metric, input[122].metric);
    sc_uint<1> comp_r_81_123 = COMP < Q > (input[81].metric, input[123].metric);
    sc_uint<1> comp_r_81_124 = COMP < Q > (input[81].metric, input[124].metric);
    sc_uint<1> comp_r_81_125 = COMP < Q > (input[81].metric, input[125].metric);
    sc_uint<1> comp_r_81_126 = COMP < Q > (input[81].metric, input[126].metric);

    sc_uint<1> comp_r_83_84 = COMP < Q > (input[83].metric, input[84].metric);
    sc_uint<1> comp_r_83_85 = COMP < Q > (input[83].metric, input[85].metric);
    sc_uint<1> comp_r_83_86 = COMP < Q > (input[83].metric, input[86].metric);
    sc_uint<1> comp_r_83_87 = COMP < Q > (input[83].metric, input[87].metric);
    sc_uint<1> comp_r_83_88 = COMP < Q > (input[83].metric, input[88].metric);
    sc_uint<1> comp_r_83_89 = COMP < Q > (input[83].metric, input[89].metric);
    sc_uint<1> comp_r_83_90 = COMP < Q > (input[83].metric, input[90].metric);
    sc_uint<1> comp_r_83_91 = COMP < Q > (input[83].metric, input[91].metric);
    sc_uint<1> comp_r_83_92 = COMP < Q > (input[83].metric, input[92].metric);
    sc_uint<1> comp_r_83_93 = COMP < Q > (input[83].metric, input[93].metric);
    sc_uint<1> comp_r_83_94 = COMP < Q > (input[83].metric, input[94].metric);
    sc_uint<1> comp_r_83_95 = COMP < Q > (input[83].metric, input[95].metric);
    sc_uint<1> comp_r_83_96 = COMP < Q > (input[83].metric, input[96].metric);
    sc_uint<1> comp_r_83_97 = COMP < Q > (input[83].metric, input[97].metric);
    sc_uint<1> comp_r_83_98 = COMP < Q > (input[83].metric, input[98].metric);
    sc_uint<1> comp_r_83_99 = COMP < Q > (input[83].metric, input[99].metric);
    sc_uint<1> comp_r_83_100 = COMP < Q > (input[83].metric, input[100].metric);
    sc_uint<1> comp_r_83_101 = COMP < Q > (input[83].metric, input[101].metric);
    sc_uint<1> comp_r_83_102 = COMP < Q > (input[83].metric, input[102].metric);
    sc_uint<1> comp_r_83_103 = COMP < Q > (input[83].metric, input[103].metric);
    sc_uint<1> comp_r_83_104 = COMP < Q > (input[83].metric, input[104].metric);
    sc_uint<1> comp_r_83_105 = COMP < Q > (input[83].metric, input[105].metric);
    sc_uint<1> comp_r_83_106 = COMP < Q > (input[83].metric, input[106].metric);
    sc_uint<1> comp_r_83_107 = COMP < Q > (input[83].metric, input[107].metric);
    sc_uint<1> comp_r_83_108 = COMP < Q > (input[83].metric, input[108].metric);
    sc_uint<1> comp_r_83_109 = COMP < Q > (input[83].metric, input[109].metric);
    sc_uint<1> comp_r_83_110 = COMP < Q > (input[83].metric, input[110].metric);
    sc_uint<1> comp_r_83_111 = COMP < Q > (input[83].metric, input[111].metric);
    sc_uint<1> comp_r_83_112 = COMP < Q > (input[83].metric, input[112].metric);
    sc_uint<1> comp_r_83_113 = COMP < Q > (input[83].metric, input[113].metric);
    sc_uint<1> comp_r_83_114 = COMP < Q > (input[83].metric, input[114].metric);
    sc_uint<1> comp_r_83_115 = COMP < Q > (input[83].metric, input[115].metric);
    sc_uint<1> comp_r_83_116 = COMP < Q > (input[83].metric, input[116].metric);
    sc_uint<1> comp_r_83_117 = COMP < Q > (input[83].metric, input[117].metric);
    sc_uint<1> comp_r_83_118 = COMP < Q > (input[83].metric, input[118].metric);
    sc_uint<1> comp_r_83_119 = COMP < Q > (input[83].metric, input[119].metric);
    sc_uint<1> comp_r_83_120 = COMP < Q > (input[83].metric, input[120].metric);
    sc_uint<1> comp_r_83_121 = COMP < Q > (input[83].metric, input[121].metric);
    sc_uint<1> comp_r_83_122 = COMP < Q > (input[83].metric, input[122].metric);
    sc_uint<1> comp_r_83_123 = COMP < Q > (input[83].metric, input[123].metric);
    sc_uint<1> comp_r_83_124 = COMP < Q > (input[83].metric, input[124].metric);
    sc_uint<1> comp_r_83_125 = COMP < Q > (input[83].metric, input[125].metric);
    sc_uint<1> comp_r_83_126 = COMP < Q > (input[83].metric, input[126].metric);

    sc_uint<1> comp_r_85_86 = COMP < Q > (input[85].metric, input[86].metric);
    sc_uint<1> comp_r_85_87 = COMP < Q > (input[85].metric, input[87].metric);
    sc_uint<1> comp_r_85_88 = COMP < Q > (input[85].metric, input[88].metric);
    sc_uint<1> comp_r_85_89 = COMP < Q > (input[85].metric, input[89].metric);
    sc_uint<1> comp_r_85_90 = COMP < Q > (input[85].metric, input[90].metric);
    sc_uint<1> comp_r_85_91 = COMP < Q > (input[85].metric, input[91].metric);
    sc_uint<1> comp_r_85_92 = COMP < Q > (input[85].metric, input[92].metric);
    sc_uint<1> comp_r_85_93 = COMP < Q > (input[85].metric, input[93].metric);
    sc_uint<1> comp_r_85_94 = COMP < Q > (input[85].metric, input[94].metric);
    sc_uint<1> comp_r_85_95 = COMP < Q > (input[85].metric, input[95].metric);
    sc_uint<1> comp_r_85_96 = COMP < Q > (input[85].metric, input[96].metric);
    sc_uint<1> comp_r_85_97 = COMP < Q > (input[85].metric, input[97].metric);
    sc_uint<1> comp_r_85_98 = COMP < Q > (input[85].metric, input[98].metric);
    sc_uint<1> comp_r_85_99 = COMP < Q > (input[85].metric, input[99].metric);
    sc_uint<1> comp_r_85_100 = COMP < Q > (input[85].metric, input[100].metric);
    sc_uint<1> comp_r_85_101 = COMP < Q > (input[85].metric, input[101].metric);
    sc_uint<1> comp_r_85_102 = COMP < Q > (input[85].metric, input[102].metric);
    sc_uint<1> comp_r_85_103 = COMP < Q > (input[85].metric, input[103].metric);
    sc_uint<1> comp_r_85_104 = COMP < Q > (input[85].metric, input[104].metric);
    sc_uint<1> comp_r_85_105 = COMP < Q > (input[85].metric, input[105].metric);
    sc_uint<1> comp_r_85_106 = COMP < Q > (input[85].metric, input[106].metric);
    sc_uint<1> comp_r_85_107 = COMP < Q > (input[85].metric, input[107].metric);
    sc_uint<1> comp_r_85_108 = COMP < Q > (input[85].metric, input[108].metric);
    sc_uint<1> comp_r_85_109 = COMP < Q > (input[85].metric, input[109].metric);
    sc_uint<1> comp_r_85_110 = COMP < Q > (input[85].metric, input[110].metric);
    sc_uint<1> comp_r_85_111 = COMP < Q > (input[85].metric, input[111].metric);
    sc_uint<1> comp_r_85_112 = COMP < Q > (input[85].metric, input[112].metric);
    sc_uint<1> comp_r_85_113 = COMP < Q > (input[85].metric, input[113].metric);
    sc_uint<1> comp_r_85_114 = COMP < Q > (input[85].metric, input[114].metric);
    sc_uint<1> comp_r_85_115 = COMP < Q > (input[85].metric, input[115].metric);
    sc_uint<1> comp_r_85_116 = COMP < Q > (input[85].metric, input[116].metric);
    sc_uint<1> comp_r_85_117 = COMP < Q > (input[85].metric, input[117].metric);
    sc_uint<1> comp_r_85_118 = COMP < Q > (input[85].metric, input[118].metric);
    sc_uint<1> comp_r_85_119 = COMP < Q > (input[85].metric, input[119].metric);
    sc_uint<1> comp_r_85_120 = COMP < Q > (input[85].metric, input[120].metric);
    sc_uint<1> comp_r_85_121 = COMP < Q > (input[85].metric, input[121].metric);
    sc_uint<1> comp_r_85_122 = COMP < Q > (input[85].metric, input[122].metric);
    sc_uint<1> comp_r_85_123 = COMP < Q > (input[85].metric, input[123].metric);
    sc_uint<1> comp_r_85_124 = COMP < Q > (input[85].metric, input[124].metric);
    sc_uint<1> comp_r_85_125 = COMP < Q > (input[85].metric, input[125].metric);
    sc_uint<1> comp_r_85_126 = COMP < Q > (input[85].metric, input[126].metric);

    sc_uint<1> comp_r_87_88 = COMP < Q > (input[87].metric, input[88].metric);
    sc_uint<1> comp_r_87_89 = COMP < Q > (input[87].metric, input[89].metric);
    sc_uint<1> comp_r_87_90 = COMP < Q > (input[87].metric, input[90].metric);
    sc_uint<1> comp_r_87_91 = COMP < Q > (input[87].metric, input[91].metric);
    sc_uint<1> comp_r_87_92 = COMP < Q > (input[87].metric, input[92].metric);
    sc_uint<1> comp_r_87_93 = COMP < Q > (input[87].metric, input[93].metric);
    sc_uint<1> comp_r_87_94 = COMP < Q > (input[87].metric, input[94].metric);
    sc_uint<1> comp_r_87_95 = COMP < Q > (input[87].metric, input[95].metric);
    sc_uint<1> comp_r_87_96 = COMP < Q > (input[87].metric, input[96].metric);
    sc_uint<1> comp_r_87_97 = COMP < Q > (input[87].metric, input[97].metric);
    sc_uint<1> comp_r_87_98 = COMP < Q > (input[87].metric, input[98].metric);
    sc_uint<1> comp_r_87_99 = COMP < Q > (input[87].metric, input[99].metric);
    sc_uint<1> comp_r_87_100 = COMP < Q > (input[87].metric, input[100].metric);
    sc_uint<1> comp_r_87_101 = COMP < Q > (input[87].metric, input[101].metric);
    sc_uint<1> comp_r_87_102 = COMP < Q > (input[87].metric, input[102].metric);
    sc_uint<1> comp_r_87_103 = COMP < Q > (input[87].metric, input[103].metric);
    sc_uint<1> comp_r_87_104 = COMP < Q > (input[87].metric, input[104].metric);
    sc_uint<1> comp_r_87_105 = COMP < Q > (input[87].metric, input[105].metric);
    sc_uint<1> comp_r_87_106 = COMP < Q > (input[87].metric, input[106].metric);
    sc_uint<1> comp_r_87_107 = COMP < Q > (input[87].metric, input[107].metric);
    sc_uint<1> comp_r_87_108 = COMP < Q > (input[87].metric, input[108].metric);
    sc_uint<1> comp_r_87_109 = COMP < Q > (input[87].metric, input[109].metric);
    sc_uint<1> comp_r_87_110 = COMP < Q > (input[87].metric, input[110].metric);
    sc_uint<1> comp_r_87_111 = COMP < Q > (input[87].metric, input[111].metric);
    sc_uint<1> comp_r_87_112 = COMP < Q > (input[87].metric, input[112].metric);
    sc_uint<1> comp_r_87_113 = COMP < Q > (input[87].metric, input[113].metric);
    sc_uint<1> comp_r_87_114 = COMP < Q > (input[87].metric, input[114].metric);
    sc_uint<1> comp_r_87_115 = COMP < Q > (input[87].metric, input[115].metric);
    sc_uint<1> comp_r_87_116 = COMP < Q > (input[87].metric, input[116].metric);
    sc_uint<1> comp_r_87_117 = COMP < Q > (input[87].metric, input[117].metric);
    sc_uint<1> comp_r_87_118 = COMP < Q > (input[87].metric, input[118].metric);
    sc_uint<1> comp_r_87_119 = COMP < Q > (input[87].metric, input[119].metric);
    sc_uint<1> comp_r_87_120 = COMP < Q > (input[87].metric, input[120].metric);
    sc_uint<1> comp_r_87_121 = COMP < Q > (input[87].metric, input[121].metric);
    sc_uint<1> comp_r_87_122 = COMP < Q > (input[87].metric, input[122].metric);
    sc_uint<1> comp_r_87_123 = COMP < Q > (input[87].metric, input[123].metric);
    sc_uint<1> comp_r_87_124 = COMP < Q > (input[87].metric, input[124].metric);
    sc_uint<1> comp_r_87_125 = COMP < Q > (input[87].metric, input[125].metric);
    sc_uint<1> comp_r_87_126 = COMP < Q > (input[87].metric, input[126].metric);

    sc_uint<1> comp_r_89_90 = COMP < Q > (input[89].metric, input[90].metric);
    sc_uint<1> comp_r_89_91 = COMP < Q > (input[89].metric, input[91].metric);
    sc_uint<1> comp_r_89_92 = COMP < Q > (input[89].metric, input[92].metric);
    sc_uint<1> comp_r_89_93 = COMP < Q > (input[89].metric, input[93].metric);
    sc_uint<1> comp_r_89_94 = COMP < Q > (input[89].metric, input[94].metric);
    sc_uint<1> comp_r_89_95 = COMP < Q > (input[89].metric, input[95].metric);
    sc_uint<1> comp_r_89_96 = COMP < Q > (input[89].metric, input[96].metric);
    sc_uint<1> comp_r_89_97 = COMP < Q > (input[89].metric, input[97].metric);
    sc_uint<1> comp_r_89_98 = COMP < Q > (input[89].metric, input[98].metric);
    sc_uint<1> comp_r_89_99 = COMP < Q > (input[89].metric, input[99].metric);
    sc_uint<1> comp_r_89_100 = COMP < Q > (input[89].metric, input[100].metric);
    sc_uint<1> comp_r_89_101 = COMP < Q > (input[89].metric, input[101].metric);
    sc_uint<1> comp_r_89_102 = COMP < Q > (input[89].metric, input[102].metric);
    sc_uint<1> comp_r_89_103 = COMP < Q > (input[89].metric, input[103].metric);
    sc_uint<1> comp_r_89_104 = COMP < Q > (input[89].metric, input[104].metric);
    sc_uint<1> comp_r_89_105 = COMP < Q > (input[89].metric, input[105].metric);
    sc_uint<1> comp_r_89_106 = COMP < Q > (input[89].metric, input[106].metric);
    sc_uint<1> comp_r_89_107 = COMP < Q > (input[89].metric, input[107].metric);
    sc_uint<1> comp_r_89_108 = COMP < Q > (input[89].metric, input[108].metric);
    sc_uint<1> comp_r_89_109 = COMP < Q > (input[89].metric, input[109].metric);
    sc_uint<1> comp_r_89_110 = COMP < Q > (input[89].metric, input[110].metric);
    sc_uint<1> comp_r_89_111 = COMP < Q > (input[89].metric, input[111].metric);
    sc_uint<1> comp_r_89_112 = COMP < Q > (input[89].metric, input[112].metric);
    sc_uint<1> comp_r_89_113 = COMP < Q > (input[89].metric, input[113].metric);
    sc_uint<1> comp_r_89_114 = COMP < Q > (input[89].metric, input[114].metric);
    sc_uint<1> comp_r_89_115 = COMP < Q > (input[89].metric, input[115].metric);
    sc_uint<1> comp_r_89_116 = COMP < Q > (input[89].metric, input[116].metric);
    sc_uint<1> comp_r_89_117 = COMP < Q > (input[89].metric, input[117].metric);
    sc_uint<1> comp_r_89_118 = COMP < Q > (input[89].metric, input[118].metric);
    sc_uint<1> comp_r_89_119 = COMP < Q > (input[89].metric, input[119].metric);
    sc_uint<1> comp_r_89_120 = COMP < Q > (input[89].metric, input[120].metric);
    sc_uint<1> comp_r_89_121 = COMP < Q > (input[89].metric, input[121].metric);
    sc_uint<1> comp_r_89_122 = COMP < Q > (input[89].metric, input[122].metric);
    sc_uint<1> comp_r_89_123 = COMP < Q > (input[89].metric, input[123].metric);
    sc_uint<1> comp_r_89_124 = COMP < Q > (input[89].metric, input[124].metric);
    sc_uint<1> comp_r_89_125 = COMP < Q > (input[89].metric, input[125].metric);
    sc_uint<1> comp_r_89_126 = COMP < Q > (input[89].metric, input[126].metric);

    sc_uint<1> comp_r_91_92 = COMP < Q > (input[91].metric, input[92].metric);
    sc_uint<1> comp_r_91_93 = COMP < Q > (input[91].metric, input[93].metric);
    sc_uint<1> comp_r_91_94 = COMP < Q > (input[91].metric, input[94].metric);
    sc_uint<1> comp_r_91_95 = COMP < Q > (input[91].metric, input[95].metric);
    sc_uint<1> comp_r_91_96 = COMP < Q > (input[91].metric, input[96].metric);
    sc_uint<1> comp_r_91_97 = COMP < Q > (input[91].metric, input[97].metric);
    sc_uint<1> comp_r_91_98 = COMP < Q > (input[91].metric, input[98].metric);
    sc_uint<1> comp_r_91_99 = COMP < Q > (input[91].metric, input[99].metric);
    sc_uint<1> comp_r_91_100 = COMP < Q > (input[91].metric, input[100].metric);
    sc_uint<1> comp_r_91_101 = COMP < Q > (input[91].metric, input[101].metric);
    sc_uint<1> comp_r_91_102 = COMP < Q > (input[91].metric, input[102].metric);
    sc_uint<1> comp_r_91_103 = COMP < Q > (input[91].metric, input[103].metric);
    sc_uint<1> comp_r_91_104 = COMP < Q > (input[91].metric, input[104].metric);
    sc_uint<1> comp_r_91_105 = COMP < Q > (input[91].metric, input[105].metric);
    sc_uint<1> comp_r_91_106 = COMP < Q > (input[91].metric, input[106].metric);
    sc_uint<1> comp_r_91_107 = COMP < Q > (input[91].metric, input[107].metric);
    sc_uint<1> comp_r_91_108 = COMP < Q > (input[91].metric, input[108].metric);
    sc_uint<1> comp_r_91_109 = COMP < Q > (input[91].metric, input[109].metric);
    sc_uint<1> comp_r_91_110 = COMP < Q > (input[91].metric, input[110].metric);
    sc_uint<1> comp_r_91_111 = COMP < Q > (input[91].metric, input[111].metric);
    sc_uint<1> comp_r_91_112 = COMP < Q > (input[91].metric, input[112].metric);
    sc_uint<1> comp_r_91_113 = COMP < Q > (input[91].metric, input[113].metric);
    sc_uint<1> comp_r_91_114 = COMP < Q > (input[91].metric, input[114].metric);
    sc_uint<1> comp_r_91_115 = COMP < Q > (input[91].metric, input[115].metric);
    sc_uint<1> comp_r_91_116 = COMP < Q > (input[91].metric, input[116].metric);
    sc_uint<1> comp_r_91_117 = COMP < Q > (input[91].metric, input[117].metric);
    sc_uint<1> comp_r_91_118 = COMP < Q > (input[91].metric, input[118].metric);
    sc_uint<1> comp_r_91_119 = COMP < Q > (input[91].metric, input[119].metric);
    sc_uint<1> comp_r_91_120 = COMP < Q > (input[91].metric, input[120].metric);
    sc_uint<1> comp_r_91_121 = COMP < Q > (input[91].metric, input[121].metric);
    sc_uint<1> comp_r_91_122 = COMP < Q > (input[91].metric, input[122].metric);
    sc_uint<1> comp_r_91_123 = COMP < Q > (input[91].metric, input[123].metric);
    sc_uint<1> comp_r_91_124 = COMP < Q > (input[91].metric, input[124].metric);
    sc_uint<1> comp_r_91_125 = COMP < Q > (input[91].metric, input[125].metric);
    sc_uint<1> comp_r_91_126 = COMP < Q > (input[91].metric, input[126].metric);

    sc_uint<1> comp_r_93_94 = COMP < Q > (input[93].metric, input[94].metric);
    sc_uint<1> comp_r_93_95 = COMP < Q > (input[93].metric, input[95].metric);
    sc_uint<1> comp_r_93_96 = COMP < Q > (input[93].metric, input[96].metric);
    sc_uint<1> comp_r_93_97 = COMP < Q > (input[93].metric, input[97].metric);
    sc_uint<1> comp_r_93_98 = COMP < Q > (input[93].metric, input[98].metric);
    sc_uint<1> comp_r_93_99 = COMP < Q > (input[93].metric, input[99].metric);
    sc_uint<1> comp_r_93_100 = COMP < Q > (input[93].metric, input[100].metric);
    sc_uint<1> comp_r_93_101 = COMP < Q > (input[93].metric, input[101].metric);
    sc_uint<1> comp_r_93_102 = COMP < Q > (input[93].metric, input[102].metric);
    sc_uint<1> comp_r_93_103 = COMP < Q > (input[93].metric, input[103].metric);
    sc_uint<1> comp_r_93_104 = COMP < Q > (input[93].metric, input[104].metric);
    sc_uint<1> comp_r_93_105 = COMP < Q > (input[93].metric, input[105].metric);
    sc_uint<1> comp_r_93_106 = COMP < Q > (input[93].metric, input[106].metric);
    sc_uint<1> comp_r_93_107 = COMP < Q > (input[93].metric, input[107].metric);
    sc_uint<1> comp_r_93_108 = COMP < Q > (input[93].metric, input[108].metric);
    sc_uint<1> comp_r_93_109 = COMP < Q > (input[93].metric, input[109].metric);
    sc_uint<1> comp_r_93_110 = COMP < Q > (input[93].metric, input[110].metric);
    sc_uint<1> comp_r_93_111 = COMP < Q > (input[93].metric, input[111].metric);
    sc_uint<1> comp_r_93_112 = COMP < Q > (input[93].metric, input[112].metric);
    sc_uint<1> comp_r_93_113 = COMP < Q > (input[93].metric, input[113].metric);
    sc_uint<1> comp_r_93_114 = COMP < Q > (input[93].metric, input[114].metric);
    sc_uint<1> comp_r_93_115 = COMP < Q > (input[93].metric, input[115].metric);
    sc_uint<1> comp_r_93_116 = COMP < Q > (input[93].metric, input[116].metric);
    sc_uint<1> comp_r_93_117 = COMP < Q > (input[93].metric, input[117].metric);
    sc_uint<1> comp_r_93_118 = COMP < Q > (input[93].metric, input[118].metric);
    sc_uint<1> comp_r_93_119 = COMP < Q > (input[93].metric, input[119].metric);
    sc_uint<1> comp_r_93_120 = COMP < Q > (input[93].metric, input[120].metric);
    sc_uint<1> comp_r_93_121 = COMP < Q > (input[93].metric, input[121].metric);
    sc_uint<1> comp_r_93_122 = COMP < Q > (input[93].metric, input[122].metric);
    sc_uint<1> comp_r_93_123 = COMP < Q > (input[93].metric, input[123].metric);
    sc_uint<1> comp_r_93_124 = COMP < Q > (input[93].metric, input[124].metric);
    sc_uint<1> comp_r_93_125 = COMP < Q > (input[93].metric, input[125].metric);
    sc_uint<1> comp_r_93_126 = COMP < Q > (input[93].metric, input[126].metric);

    sc_uint<1> comp_r_95_96 = COMP < Q > (input[95].metric, input[96].metric);
    sc_uint<1> comp_r_95_97 = COMP < Q > (input[95].metric, input[97].metric);
    sc_uint<1> comp_r_95_98 = COMP < Q > (input[95].metric, input[98].metric);
    sc_uint<1> comp_r_95_99 = COMP < Q > (input[95].metric, input[99].metric);
    sc_uint<1> comp_r_95_100 = COMP < Q > (input[95].metric, input[100].metric);
    sc_uint<1> comp_r_95_101 = COMP < Q > (input[95].metric, input[101].metric);
    sc_uint<1> comp_r_95_102 = COMP < Q > (input[95].metric, input[102].metric);
    sc_uint<1> comp_r_95_103 = COMP < Q > (input[95].metric, input[103].metric);
    sc_uint<1> comp_r_95_104 = COMP < Q > (input[95].metric, input[104].metric);
    sc_uint<1> comp_r_95_105 = COMP < Q > (input[95].metric, input[105].metric);
    sc_uint<1> comp_r_95_106 = COMP < Q > (input[95].metric, input[106].metric);
    sc_uint<1> comp_r_95_107 = COMP < Q > (input[95].metric, input[107].metric);
    sc_uint<1> comp_r_95_108 = COMP < Q > (input[95].metric, input[108].metric);
    sc_uint<1> comp_r_95_109 = COMP < Q > (input[95].metric, input[109].metric);
    sc_uint<1> comp_r_95_110 = COMP < Q > (input[95].metric, input[110].metric);
    sc_uint<1> comp_r_95_111 = COMP < Q > (input[95].metric, input[111].metric);
    sc_uint<1> comp_r_95_112 = COMP < Q > (input[95].metric, input[112].metric);
    sc_uint<1> comp_r_95_113 = COMP < Q > (input[95].metric, input[113].metric);
    sc_uint<1> comp_r_95_114 = COMP < Q > (input[95].metric, input[114].metric);
    sc_uint<1> comp_r_95_115 = COMP < Q > (input[95].metric, input[115].metric);
    sc_uint<1> comp_r_95_116 = COMP < Q > (input[95].metric, input[116].metric);
    sc_uint<1> comp_r_95_117 = COMP < Q > (input[95].metric, input[117].metric);
    sc_uint<1> comp_r_95_118 = COMP < Q > (input[95].metric, input[118].metric);
    sc_uint<1> comp_r_95_119 = COMP < Q > (input[95].metric, input[119].metric);
    sc_uint<1> comp_r_95_120 = COMP < Q > (input[95].metric, input[120].metric);
    sc_uint<1> comp_r_95_121 = COMP < Q > (input[95].metric, input[121].metric);
    sc_uint<1> comp_r_95_122 = COMP < Q > (input[95].metric, input[122].metric);
    sc_uint<1> comp_r_95_123 = COMP < Q > (input[95].metric, input[123].metric);
    sc_uint<1> comp_r_95_124 = COMP < Q > (input[95].metric, input[124].metric);
    sc_uint<1> comp_r_95_125 = COMP < Q > (input[95].metric, input[125].metric);
    sc_uint<1> comp_r_95_126 = COMP < Q > (input[95].metric, input[126].metric);

    sc_uint<1> comp_r_97_98 = COMP < Q > (input[97].metric, input[98].metric);
    sc_uint<1> comp_r_97_99 = COMP < Q > (input[97].metric, input[99].metric);
    sc_uint<1> comp_r_97_100 = COMP < Q > (input[97].metric, input[100].metric);
    sc_uint<1> comp_r_97_101 = COMP < Q > (input[97].metric, input[101].metric);
    sc_uint<1> comp_r_97_102 = COMP < Q > (input[97].metric, input[102].metric);
    sc_uint<1> comp_r_97_103 = COMP < Q > (input[97].metric, input[103].metric);
    sc_uint<1> comp_r_97_104 = COMP < Q > (input[97].metric, input[104].metric);
    sc_uint<1> comp_r_97_105 = COMP < Q > (input[97].metric, input[105].metric);
    sc_uint<1> comp_r_97_106 = COMP < Q > (input[97].metric, input[106].metric);
    sc_uint<1> comp_r_97_107 = COMP < Q > (input[97].metric, input[107].metric);
    sc_uint<1> comp_r_97_108 = COMP < Q > (input[97].metric, input[108].metric);
    sc_uint<1> comp_r_97_109 = COMP < Q > (input[97].metric, input[109].metric);
    sc_uint<1> comp_r_97_110 = COMP < Q > (input[97].metric, input[110].metric);
    sc_uint<1> comp_r_97_111 = COMP < Q > (input[97].metric, input[111].metric);
    sc_uint<1> comp_r_97_112 = COMP < Q > (input[97].metric, input[112].metric);
    sc_uint<1> comp_r_97_113 = COMP < Q > (input[97].metric, input[113].metric);
    sc_uint<1> comp_r_97_114 = COMP < Q > (input[97].metric, input[114].metric);
    sc_uint<1> comp_r_97_115 = COMP < Q > (input[97].metric, input[115].metric);
    sc_uint<1> comp_r_97_116 = COMP < Q > (input[97].metric, input[116].metric);
    sc_uint<1> comp_r_97_117 = COMP < Q > (input[97].metric, input[117].metric);
    sc_uint<1> comp_r_97_118 = COMP < Q > (input[97].metric, input[118].metric);
    sc_uint<1> comp_r_97_119 = COMP < Q > (input[97].metric, input[119].metric);
    sc_uint<1> comp_r_97_120 = COMP < Q > (input[97].metric, input[120].metric);
    sc_uint<1> comp_r_97_121 = COMP < Q > (input[97].metric, input[121].metric);
    sc_uint<1> comp_r_97_122 = COMP < Q > (input[97].metric, input[122].metric);
    sc_uint<1> comp_r_97_123 = COMP < Q > (input[97].metric, input[123].metric);
    sc_uint<1> comp_r_97_124 = COMP < Q > (input[97].metric, input[124].metric);
    sc_uint<1> comp_r_97_125 = COMP < Q > (input[97].metric, input[125].metric);
    sc_uint<1> comp_r_97_126 = COMP < Q > (input[97].metric, input[126].metric);

    sc_uint<1> comp_r_99_100 = COMP < Q > (input[99].metric, input[100].metric);
    sc_uint<1> comp_r_99_101 = COMP < Q > (input[99].metric, input[101].metric);
    sc_uint<1> comp_r_99_102 = COMP < Q > (input[99].metric, input[102].metric);
    sc_uint<1> comp_r_99_103 = COMP < Q > (input[99].metric, input[103].metric);
    sc_uint<1> comp_r_99_104 = COMP < Q > (input[99].metric, input[104].metric);
    sc_uint<1> comp_r_99_105 = COMP < Q > (input[99].metric, input[105].metric);
    sc_uint<1> comp_r_99_106 = COMP < Q > (input[99].metric, input[106].metric);
    sc_uint<1> comp_r_99_107 = COMP < Q > (input[99].metric, input[107].metric);
    sc_uint<1> comp_r_99_108 = COMP < Q > (input[99].metric, input[108].metric);
    sc_uint<1> comp_r_99_109 = COMP < Q > (input[99].metric, input[109].metric);
    sc_uint<1> comp_r_99_110 = COMP < Q > (input[99].metric, input[110].metric);
    sc_uint<1> comp_r_99_111 = COMP < Q > (input[99].metric, input[111].metric);
    sc_uint<1> comp_r_99_112 = COMP < Q > (input[99].metric, input[112].metric);
    sc_uint<1> comp_r_99_113 = COMP < Q > (input[99].metric, input[113].metric);
    sc_uint<1> comp_r_99_114 = COMP < Q > (input[99].metric, input[114].metric);
    sc_uint<1> comp_r_99_115 = COMP < Q > (input[99].metric, input[115].metric);
    sc_uint<1> comp_r_99_116 = COMP < Q > (input[99].metric, input[116].metric);
    sc_uint<1> comp_r_99_117 = COMP < Q > (input[99].metric, input[117].metric);
    sc_uint<1> comp_r_99_118 = COMP < Q > (input[99].metric, input[118].metric);
    sc_uint<1> comp_r_99_119 = COMP < Q > (input[99].metric, input[119].metric);
    sc_uint<1> comp_r_99_120 = COMP < Q > (input[99].metric, input[120].metric);
    sc_uint<1> comp_r_99_121 = COMP < Q > (input[99].metric, input[121].metric);
    sc_uint<1> comp_r_99_122 = COMP < Q > (input[99].metric, input[122].metric);
    sc_uint<1> comp_r_99_123 = COMP < Q > (input[99].metric, input[123].metric);
    sc_uint<1> comp_r_99_124 = COMP < Q > (input[99].metric, input[124].metric);
    sc_uint<1> comp_r_99_125 = COMP < Q > (input[99].metric, input[125].metric);
    sc_uint<1> comp_r_99_126 = COMP < Q > (input[99].metric, input[126].metric);

    sc_uint<1> comp_r_101_102 = COMP < Q > (input[101].metric, input[102].metric);
    sc_uint<1> comp_r_101_103 = COMP < Q > (input[101].metric, input[103].metric);
    sc_uint<1> comp_r_101_104 = COMP < Q > (input[101].metric, input[104].metric);
    sc_uint<1> comp_r_101_105 = COMP < Q > (input[101].metric, input[105].metric);
    sc_uint<1> comp_r_101_106 = COMP < Q > (input[101].metric, input[106].metric);
    sc_uint<1> comp_r_101_107 = COMP < Q > (input[101].metric, input[107].metric);
    sc_uint<1> comp_r_101_108 = COMP < Q > (input[101].metric, input[108].metric);
    sc_uint<1> comp_r_101_109 = COMP < Q > (input[101].metric, input[109].metric);
    sc_uint<1> comp_r_101_110 = COMP < Q > (input[101].metric, input[110].metric);
    sc_uint<1> comp_r_101_111 = COMP < Q > (input[101].metric, input[111].metric);
    sc_uint<1> comp_r_101_112 = COMP < Q > (input[101].metric, input[112].metric);
    sc_uint<1> comp_r_101_113 = COMP < Q > (input[101].metric, input[113].metric);
    sc_uint<1> comp_r_101_114 = COMP < Q > (input[101].metric, input[114].metric);
    sc_uint<1> comp_r_101_115 = COMP < Q > (input[101].metric, input[115].metric);
    sc_uint<1> comp_r_101_116 = COMP < Q > (input[101].metric, input[116].metric);
    sc_uint<1> comp_r_101_117 = COMP < Q > (input[101].metric, input[117].metric);
    sc_uint<1> comp_r_101_118 = COMP < Q > (input[101].metric, input[118].metric);
    sc_uint<1> comp_r_101_119 = COMP < Q > (input[101].metric, input[119].metric);
    sc_uint<1> comp_r_101_120 = COMP < Q > (input[101].metric, input[120].metric);
    sc_uint<1> comp_r_101_121 = COMP < Q > (input[101].metric, input[121].metric);
    sc_uint<1> comp_r_101_122 = COMP < Q > (input[101].metric, input[122].metric);
    sc_uint<1> comp_r_101_123 = COMP < Q > (input[101].metric, input[123].metric);
    sc_uint<1> comp_r_101_124 = COMP < Q > (input[101].metric, input[124].metric);
    sc_uint<1> comp_r_101_125 = COMP < Q > (input[101].metric, input[125].metric);
    sc_uint<1> comp_r_101_126 = COMP < Q > (input[101].metric, input[126].metric);

    sc_uint<1> comp_r_103_104 = COMP < Q > (input[103].metric, input[104].metric);
    sc_uint<1> comp_r_103_105 = COMP < Q > (input[103].metric, input[105].metric);
    sc_uint<1> comp_r_103_106 = COMP < Q > (input[103].metric, input[106].metric);
    sc_uint<1> comp_r_103_107 = COMP < Q > (input[103].metric, input[107].metric);
    sc_uint<1> comp_r_103_108 = COMP < Q > (input[103].metric, input[108].metric);
    sc_uint<1> comp_r_103_109 = COMP < Q > (input[103].metric, input[109].metric);
    sc_uint<1> comp_r_103_110 = COMP < Q > (input[103].metric, input[110].metric);
    sc_uint<1> comp_r_103_111 = COMP < Q > (input[103].metric, input[111].metric);
    sc_uint<1> comp_r_103_112 = COMP < Q > (input[103].metric, input[112].metric);
    sc_uint<1> comp_r_103_113 = COMP < Q > (input[103].metric, input[113].metric);
    sc_uint<1> comp_r_103_114 = COMP < Q > (input[103].metric, input[114].metric);
    sc_uint<1> comp_r_103_115 = COMP < Q > (input[103].metric, input[115].metric);
    sc_uint<1> comp_r_103_116 = COMP < Q > (input[103].metric, input[116].metric);
    sc_uint<1> comp_r_103_117 = COMP < Q > (input[103].metric, input[117].metric);
    sc_uint<1> comp_r_103_118 = COMP < Q > (input[103].metric, input[118].metric);
    sc_uint<1> comp_r_103_119 = COMP < Q > (input[103].metric, input[119].metric);
    sc_uint<1> comp_r_103_120 = COMP < Q > (input[103].metric, input[120].metric);
    sc_uint<1> comp_r_103_121 = COMP < Q > (input[103].metric, input[121].metric);
    sc_uint<1> comp_r_103_122 = COMP < Q > (input[103].metric, input[122].metric);
    sc_uint<1> comp_r_103_123 = COMP < Q > (input[103].metric, input[123].metric);
    sc_uint<1> comp_r_103_124 = COMP < Q > (input[103].metric, input[124].metric);
    sc_uint<1> comp_r_103_125 = COMP < Q > (input[103].metric, input[125].metric);
    sc_uint<1> comp_r_103_126 = COMP < Q > (input[103].metric, input[126].metric);

    sc_uint<1> comp_r_105_106 = COMP < Q > (input[105].metric, input[106].metric);
    sc_uint<1> comp_r_105_107 = COMP < Q > (input[105].metric, input[107].metric);
    sc_uint<1> comp_r_105_108 = COMP < Q > (input[105].metric, input[108].metric);
    sc_uint<1> comp_r_105_109 = COMP < Q > (input[105].metric, input[109].metric);
    sc_uint<1> comp_r_105_110 = COMP < Q > (input[105].metric, input[110].metric);
    sc_uint<1> comp_r_105_111 = COMP < Q > (input[105].metric, input[111].metric);
    sc_uint<1> comp_r_105_112 = COMP < Q > (input[105].metric, input[112].metric);
    sc_uint<1> comp_r_105_113 = COMP < Q > (input[105].metric, input[113].metric);
    sc_uint<1> comp_r_105_114 = COMP < Q > (input[105].metric, input[114].metric);
    sc_uint<1> comp_r_105_115 = COMP < Q > (input[105].metric, input[115].metric);
    sc_uint<1> comp_r_105_116 = COMP < Q > (input[105].metric, input[116].metric);
    sc_uint<1> comp_r_105_117 = COMP < Q > (input[105].metric, input[117].metric);
    sc_uint<1> comp_r_105_118 = COMP < Q > (input[105].metric, input[118].metric);
    sc_uint<1> comp_r_105_119 = COMP < Q > (input[105].metric, input[119].metric);
    sc_uint<1> comp_r_105_120 = COMP < Q > (input[105].metric, input[120].metric);
    sc_uint<1> comp_r_105_121 = COMP < Q > (input[105].metric, input[121].metric);
    sc_uint<1> comp_r_105_122 = COMP < Q > (input[105].metric, input[122].metric);
    sc_uint<1> comp_r_105_123 = COMP < Q > (input[105].metric, input[123].metric);
    sc_uint<1> comp_r_105_124 = COMP < Q > (input[105].metric, input[124].metric);
    sc_uint<1> comp_r_105_125 = COMP < Q > (input[105].metric, input[125].metric);
    sc_uint<1> comp_r_105_126 = COMP < Q > (input[105].metric, input[126].metric);

    sc_uint<1> comp_r_107_108 = COMP < Q > (input[107].metric, input[108].metric);
    sc_uint<1> comp_r_107_109 = COMP < Q > (input[107].metric, input[109].metric);
    sc_uint<1> comp_r_107_110 = COMP < Q > (input[107].metric, input[110].metric);
    sc_uint<1> comp_r_107_111 = COMP < Q > (input[107].metric, input[111].metric);
    sc_uint<1> comp_r_107_112 = COMP < Q > (input[107].metric, input[112].metric);
    sc_uint<1> comp_r_107_113 = COMP < Q > (input[107].metric, input[113].metric);
    sc_uint<1> comp_r_107_114 = COMP < Q > (input[107].metric, input[114].metric);
    sc_uint<1> comp_r_107_115 = COMP < Q > (input[107].metric, input[115].metric);
    sc_uint<1> comp_r_107_116 = COMP < Q > (input[107].metric, input[116].metric);
    sc_uint<1> comp_r_107_117 = COMP < Q > (input[107].metric, input[117].metric);
    sc_uint<1> comp_r_107_118 = COMP < Q > (input[107].metric, input[118].metric);
    sc_uint<1> comp_r_107_119 = COMP < Q > (input[107].metric, input[119].metric);
    sc_uint<1> comp_r_107_120 = COMP < Q > (input[107].metric, input[120].metric);
    sc_uint<1> comp_r_107_121 = COMP < Q > (input[107].metric, input[121].metric);
    sc_uint<1> comp_r_107_122 = COMP < Q > (input[107].metric, input[122].metric);
    sc_uint<1> comp_r_107_123 = COMP < Q > (input[107].metric, input[123].metric);
    sc_uint<1> comp_r_107_124 = COMP < Q > (input[107].metric, input[124].metric);
    sc_uint<1> comp_r_107_125 = COMP < Q > (input[107].metric, input[125].metric);
    sc_uint<1> comp_r_107_126 = COMP < Q > (input[107].metric, input[126].metric);

    sc_uint<1> comp_r_109_110 = COMP < Q > (input[109].metric, input[110].metric);
    sc_uint<1> comp_r_109_111 = COMP < Q > (input[109].metric, input[111].metric);
    sc_uint<1> comp_r_109_112 = COMP < Q > (input[109].metric, input[112].metric);
    sc_uint<1> comp_r_109_113 = COMP < Q > (input[109].metric, input[113].metric);
    sc_uint<1> comp_r_109_114 = COMP < Q > (input[109].metric, input[114].metric);
    sc_uint<1> comp_r_109_115 = COMP < Q > (input[109].metric, input[115].metric);
    sc_uint<1> comp_r_109_116 = COMP < Q > (input[109].metric, input[116].metric);
    sc_uint<1> comp_r_109_117 = COMP < Q > (input[109].metric, input[117].metric);
    sc_uint<1> comp_r_109_118 = COMP < Q > (input[109].metric, input[118].metric);
    sc_uint<1> comp_r_109_119 = COMP < Q > (input[109].metric, input[119].metric);
    sc_uint<1> comp_r_109_120 = COMP < Q > (input[109].metric, input[120].metric);
    sc_uint<1> comp_r_109_121 = COMP < Q > (input[109].metric, input[121].metric);
    sc_uint<1> comp_r_109_122 = COMP < Q > (input[109].metric, input[122].metric);
    sc_uint<1> comp_r_109_123 = COMP < Q > (input[109].metric, input[123].metric);
    sc_uint<1> comp_r_109_124 = COMP < Q > (input[109].metric, input[124].metric);
    sc_uint<1> comp_r_109_125 = COMP < Q > (input[109].metric, input[125].metric);
    sc_uint<1> comp_r_109_126 = COMP < Q > (input[109].metric, input[126].metric);

    sc_uint<1> comp_r_111_112 = COMP < Q > (input[111].metric, input[112].metric);
    sc_uint<1> comp_r_111_113 = COMP < Q > (input[111].metric, input[113].metric);
    sc_uint<1> comp_r_111_114 = COMP < Q > (input[111].metric, input[114].metric);
    sc_uint<1> comp_r_111_115 = COMP < Q > (input[111].metric, input[115].metric);
    sc_uint<1> comp_r_111_116 = COMP < Q > (input[111].metric, input[116].metric);
    sc_uint<1> comp_r_111_117 = COMP < Q > (input[111].metric, input[117].metric);
    sc_uint<1> comp_r_111_118 = COMP < Q > (input[111].metric, input[118].metric);
    sc_uint<1> comp_r_111_119 = COMP < Q > (input[111].metric, input[119].metric);
    sc_uint<1> comp_r_111_120 = COMP < Q > (input[111].metric, input[120].metric);
    sc_uint<1> comp_r_111_121 = COMP < Q > (input[111].metric, input[121].metric);
    sc_uint<1> comp_r_111_122 = COMP < Q > (input[111].metric, input[122].metric);
    sc_uint<1> comp_r_111_123 = COMP < Q > (input[111].metric, input[123].metric);
    sc_uint<1> comp_r_111_124 = COMP < Q > (input[111].metric, input[124].metric);
    sc_uint<1> comp_r_111_125 = COMP < Q > (input[111].metric, input[125].metric);
    sc_uint<1> comp_r_111_126 = COMP < Q > (input[111].metric, input[126].metric);

    sc_uint<1> comp_r_113_114 = COMP < Q > (input[113].metric, input[114].metric);
    sc_uint<1> comp_r_113_115 = COMP < Q > (input[113].metric, input[115].metric);
    sc_uint<1> comp_r_113_116 = COMP < Q > (input[113].metric, input[116].metric);
    sc_uint<1> comp_r_113_117 = COMP < Q > (input[113].metric, input[117].metric);
    sc_uint<1> comp_r_113_118 = COMP < Q > (input[113].metric, input[118].metric);
    sc_uint<1> comp_r_113_119 = COMP < Q > (input[113].metric, input[119].metric);
    sc_uint<1> comp_r_113_120 = COMP < Q > (input[113].metric, input[120].metric);
    sc_uint<1> comp_r_113_121 = COMP < Q > (input[113].metric, input[121].metric);
    sc_uint<1> comp_r_113_122 = COMP < Q > (input[113].metric, input[122].metric);
    sc_uint<1> comp_r_113_123 = COMP < Q > (input[113].metric, input[123].metric);
    sc_uint<1> comp_r_113_124 = COMP < Q > (input[113].metric, input[124].metric);
    sc_uint<1> comp_r_113_125 = COMP < Q > (input[113].metric, input[125].metric);
    sc_uint<1> comp_r_113_126 = COMP < Q > (input[113].metric, input[126].metric);

    sc_uint<1> comp_r_115_116 = COMP < Q > (input[115].metric, input[116].metric);
    sc_uint<1> comp_r_115_117 = COMP < Q > (input[115].metric, input[117].metric);
    sc_uint<1> comp_r_115_118 = COMP < Q > (input[115].metric, input[118].metric);
    sc_uint<1> comp_r_115_119 = COMP < Q > (input[115].metric, input[119].metric);
    sc_uint<1> comp_r_115_120 = COMP < Q > (input[115].metric, input[120].metric);
    sc_uint<1> comp_r_115_121 = COMP < Q > (input[115].metric, input[121].metric);
    sc_uint<1> comp_r_115_122 = COMP < Q > (input[115].metric, input[122].metric);
    sc_uint<1> comp_r_115_123 = COMP < Q > (input[115].metric, input[123].metric);
    sc_uint<1> comp_r_115_124 = COMP < Q > (input[115].metric, input[124].metric);
    sc_uint<1> comp_r_115_125 = COMP < Q > (input[115].metric, input[125].metric);
    sc_uint<1> comp_r_115_126 = COMP < Q > (input[115].metric, input[126].metric);

    sc_uint<1> comp_r_117_118 = COMP < Q > (input[117].metric, input[118].metric);
    sc_uint<1> comp_r_117_119 = COMP < Q > (input[117].metric, input[119].metric);
    sc_uint<1> comp_r_117_120 = COMP < Q > (input[117].metric, input[120].metric);
    sc_uint<1> comp_r_117_121 = COMP < Q > (input[117].metric, input[121].metric);
    sc_uint<1> comp_r_117_122 = COMP < Q > (input[117].metric, input[122].metric);
    sc_uint<1> comp_r_117_123 = COMP < Q > (input[117].metric, input[123].metric);
    sc_uint<1> comp_r_117_124 = COMP < Q > (input[117].metric, input[124].metric);
    sc_uint<1> comp_r_117_125 = COMP < Q > (input[117].metric, input[125].metric);
    sc_uint<1> comp_r_117_126 = COMP < Q > (input[117].metric, input[126].metric);

    sc_uint<1> comp_r_119_120 = COMP < Q > (input[119].metric, input[120].metric);
    sc_uint<1> comp_r_119_121 = COMP < Q > (input[119].metric, input[121].metric);
    sc_uint<1> comp_r_119_122 = COMP < Q > (input[119].metric, input[122].metric);
    sc_uint<1> comp_r_119_123 = COMP < Q > (input[119].metric, input[123].metric);
    sc_uint<1> comp_r_119_124 = COMP < Q > (input[119].metric, input[124].metric);
    sc_uint<1> comp_r_119_125 = COMP < Q > (input[119].metric, input[125].metric);
    sc_uint<1> comp_r_119_126 = COMP < Q > (input[119].metric, input[126].metric);

    sc_uint<1> comp_r_121_122 = COMP < Q > (input[121].metric, input[122].metric);
    sc_uint<1> comp_r_121_123 = COMP < Q > (input[121].metric, input[123].metric);
    sc_uint<1> comp_r_121_124 = COMP < Q > (input[121].metric, input[124].metric);
    sc_uint<1> comp_r_121_125 = COMP < Q > (input[121].metric, input[125].metric);
    sc_uint<1> comp_r_121_126 = COMP < Q > (input[121].metric, input[126].metric);

    sc_uint<1> comp_r_123_124 = COMP < Q > (input[123].metric, input[124].metric);
    sc_uint<1> comp_r_123_125 = COMP < Q > (input[123].metric, input[125].metric);
    sc_uint<1> comp_r_123_126 = COMP < Q > (input[123].metric, input[126].metric);

    sc_uint<1> comp_r_125_126 = COMP < Q > (input[125].metric, input[126].metric);

// COMPUTE positions

    position[0] = 0;
    position[1] = 1 + comp_r_1_2 + comp_r_1_3 + comp_r_1_4 + comp_r_1_5 + comp_r_1_6 + comp_r_1_7 + comp_r_1_8 + comp_r_1_9 + comp_r_1_10 + comp_r_1_11 + comp_r_1_12 + comp_r_1_13 + comp_r_1_14 + comp_r_1_15 + comp_r_1_16 + comp_r_1_17 + comp_r_1_18 + comp_r_1_19 + comp_r_1_20 + comp_r_1_21 + comp_r_1_22 + comp_r_1_23 + comp_r_1_24 + comp_r_1_25 + comp_r_1_26 + comp_r_1_27 + comp_r_1_28 + comp_r_1_29 + comp_r_1_30 + comp_r_1_31 + comp_r_1_32 + comp_r_1_33 + comp_r_1_34 + comp_r_1_35 + comp_r_1_36 + comp_r_1_37 + comp_r_1_38 + comp_r_1_39 + comp_r_1_40 + comp_r_1_41 + comp_r_1_42 + comp_r_1_43 + comp_r_1_44 + comp_r_1_45 + comp_r_1_46 + comp_r_1_47 + comp_r_1_48 + comp_r_1_49 + comp_r_1_50 + comp_r_1_51 + comp_r_1_52 + comp_r_1_53 + comp_r_1_54 + comp_r_1_55 + comp_r_1_56 + comp_r_1_57 + comp_r_1_58 + comp_r_1_59 + comp_r_1_60 + comp_r_1_61 + comp_r_1_62 + comp_r_1_63 + comp_r_1_64 + comp_r_1_65 + comp_r_1_66 + comp_r_1_67 + comp_r_1_68 + comp_r_1_69 + comp_r_1_70 + comp_r_1_71 + comp_r_1_72 + comp_r_1_73 + comp_r_1_74 + comp_r_1_75 + comp_r_1_76 + comp_r_1_77 + comp_r_1_78 + comp_r_1_79 + comp_r_1_80 + comp_r_1_81 + comp_r_1_82 + comp_r_1_83 + comp_r_1_84 + comp_r_1_85 + comp_r_1_86 + comp_r_1_87 + comp_r_1_88 + comp_r_1_89 + comp_r_1_90 + comp_r_1_91 + comp_r_1_92 + comp_r_1_93 + comp_r_1_94 + comp_r_1_95 + comp_r_1_96 + comp_r_1_97 + comp_r_1_98 + comp_r_1_99 + comp_r_1_100 + comp_r_1_101 + comp_r_1_102 + comp_r_1_103 + comp_r_1_104 + comp_r_1_105 + comp_r_1_106 + comp_r_1_107 + comp_r_1_108 + comp_r_1_109 + comp_r_1_110 + comp_r_1_111 + comp_r_1_112 + comp_r_1_113 + comp_r_1_114 + comp_r_1_115 + comp_r_1_116 + comp_r_1_117 + comp_r_1_118 + comp_r_1_119 + comp_r_1_120 + comp_r_1_121 + comp_r_1_122 + comp_r_1_123 + comp_r_1_124 + comp_r_1_125 + comp_r_1_126;
    position[2] = 1 + (sc_uint<1>) ~(comp_r_1_2);
    position[3] = 2 + (sc_uint<1>) ~(comp_r_1_3) + comp_r_3_4 + comp_r_3_5 + comp_r_3_6 + comp_r_3_7 + comp_r_3_8 + comp_r_3_9 + comp_r_3_10 + comp_r_3_11 + comp_r_3_12 + comp_r_3_13 + comp_r_3_14 + comp_r_3_15 + comp_r_3_16 + comp_r_3_17 + comp_r_3_18 + comp_r_3_19 + comp_r_3_20 + comp_r_3_21 + comp_r_3_22 + comp_r_3_23 + comp_r_3_24 + comp_r_3_25 + comp_r_3_26 + comp_r_3_27 + comp_r_3_28 + comp_r_3_29 + comp_r_3_30 + comp_r_3_31 + comp_r_3_32 + comp_r_3_33 + comp_r_3_34 + comp_r_3_35 + comp_r_3_36 + comp_r_3_37 + comp_r_3_38 + comp_r_3_39 + comp_r_3_40 + comp_r_3_41 + comp_r_3_42 + comp_r_3_43 + comp_r_3_44 + comp_r_3_45 + comp_r_3_46 + comp_r_3_47 + comp_r_3_48 + comp_r_3_49 + comp_r_3_50 + comp_r_3_51 + comp_r_3_52 + comp_r_3_53 + comp_r_3_54 + comp_r_3_55 + comp_r_3_56 + comp_r_3_57 + comp_r_3_58 + comp_r_3_59 + comp_r_3_60 + comp_r_3_61 + comp_r_3_62 + comp_r_3_63 + comp_r_3_64 + comp_r_3_65 + comp_r_3_66 + comp_r_3_67 + comp_r_3_68 + comp_r_3_69 + comp_r_3_70 + comp_r_3_71 + comp_r_3_72 + comp_r_3_73 + comp_r_3_74 + comp_r_3_75 + comp_r_3_76 + comp_r_3_77 + comp_r_3_78 + comp_r_3_79 + comp_r_3_80 + comp_r_3_81 + comp_r_3_82 + comp_r_3_83 + comp_r_3_84 + comp_r_3_85 + comp_r_3_86 + comp_r_3_87 + comp_r_3_88 + comp_r_3_89 + comp_r_3_90 + comp_r_3_91 + comp_r_3_92 + comp_r_3_93 + comp_r_3_94 + comp_r_3_95 + comp_r_3_96 + comp_r_3_97 + comp_r_3_98 + comp_r_3_99 + comp_r_3_100 + comp_r_3_101 + comp_r_3_102 + comp_r_3_103 + comp_r_3_104 + comp_r_3_105 + comp_r_3_106 + comp_r_3_107 + comp_r_3_108 + comp_r_3_109 + comp_r_3_110 + comp_r_3_111 + comp_r_3_112 + comp_r_3_113 + comp_r_3_114 + comp_r_3_115 + comp_r_3_116 + comp_r_3_117 + comp_r_3_118 + comp_r_3_119 + comp_r_3_120 + comp_r_3_121 + comp_r_3_122 + comp_r_3_123 + comp_r_3_124 + comp_r_3_125 + comp_r_3_126;
    position[4] = 2 + (sc_uint<1>) ~(comp_r_1_4) + (sc_uint<1>) ~(comp_r_3_4);
    position[5] = 3 + (sc_uint<1>) ~(comp_r_1_5) + (sc_uint<1>) ~(comp_r_3_5) + comp_r_5_6 + comp_r_5_7 + comp_r_5_8 + comp_r_5_9 + comp_r_5_10 + comp_r_5_11 + comp_r_5_12 + comp_r_5_13 + comp_r_5_14 + comp_r_5_15 + comp_r_5_16 + comp_r_5_17 + comp_r_5_18 + comp_r_5_19 + comp_r_5_20 + comp_r_5_21 + comp_r_5_22 + comp_r_5_23 + comp_r_5_24 + comp_r_5_25 + comp_r_5_26 + comp_r_5_27 + comp_r_5_28 + comp_r_5_29 + comp_r_5_30 + comp_r_5_31 + comp_r_5_32 + comp_r_5_33 + comp_r_5_34 + comp_r_5_35 + comp_r_5_36 + comp_r_5_37 + comp_r_5_38 + comp_r_5_39 + comp_r_5_40 + comp_r_5_41 + comp_r_5_42 + comp_r_5_43 + comp_r_5_44 + comp_r_5_45 + comp_r_5_46 + comp_r_5_47 + comp_r_5_48 + comp_r_5_49 + comp_r_5_50 + comp_r_5_51 + comp_r_5_52 + comp_r_5_53 + comp_r_5_54 + comp_r_5_55 + comp_r_5_56 + comp_r_5_57 + comp_r_5_58 + comp_r_5_59 + comp_r_5_60 + comp_r_5_61 + comp_r_5_62 + comp_r_5_63 + comp_r_5_64 + comp_r_5_65 + comp_r_5_66 + comp_r_5_67 + comp_r_5_68 + comp_r_5_69 + comp_r_5_70 + comp_r_5_71 + comp_r_5_72 + comp_r_5_73 + comp_r_5_74 + comp_r_5_75 + comp_r_5_76 + comp_r_5_77 + comp_r_5_78 + comp_r_5_79 + comp_r_5_80 + comp_r_5_81 + comp_r_5_82 + comp_r_5_83 + comp_r_5_84 + comp_r_5_85 + comp_r_5_86 + comp_r_5_87 + comp_r_5_88 + comp_r_5_89 + comp_r_5_90 + comp_r_5_91 + comp_r_5_92 + comp_r_5_93 + comp_r_5_94 + comp_r_5_95 + comp_r_5_96 + comp_r_5_97 + comp_r_5_98 + comp_r_5_99 + comp_r_5_100 + comp_r_5_101 + comp_r_5_102 + comp_r_5_103 + comp_r_5_104 + comp_r_5_105 + comp_r_5_106 + comp_r_5_107 + comp_r_5_108 + comp_r_5_109 + comp_r_5_110 + comp_r_5_111 + comp_r_5_112 + comp_r_5_113 + comp_r_5_114 + comp_r_5_115 + comp_r_5_116 + comp_r_5_117 + comp_r_5_118 + comp_r_5_119 + comp_r_5_120 + comp_r_5_121 + comp_r_5_122 + comp_r_5_123 + comp_r_5_124 + comp_r_5_125 + comp_r_5_126;
    position[6] = 3 + (sc_uint<1>) ~(comp_r_1_6) + (sc_uint<1>) ~(comp_r_3_6) + (sc_uint<1>) ~(comp_r_5_6);
    position[7] = 4 + (sc_uint<1>) ~(comp_r_1_7) + (sc_uint<1>) ~(comp_r_3_7) + (sc_uint<1>) ~(comp_r_5_7) + comp_r_7_8 + comp_r_7_9 + comp_r_7_10 + comp_r_7_11 + comp_r_7_12 + comp_r_7_13 + comp_r_7_14 + comp_r_7_15 + comp_r_7_16 + comp_r_7_17 + comp_r_7_18 + comp_r_7_19 + comp_r_7_20 + comp_r_7_21 + comp_r_7_22 + comp_r_7_23 + comp_r_7_24 + comp_r_7_25 + comp_r_7_26 + comp_r_7_27 + comp_r_7_28 + comp_r_7_29 + comp_r_7_30 + comp_r_7_31 + comp_r_7_32 + comp_r_7_33 + comp_r_7_34 + comp_r_7_35 + comp_r_7_36 + comp_r_7_37 + comp_r_7_38 + comp_r_7_39 + comp_r_7_40 + comp_r_7_41 + comp_r_7_42 + comp_r_7_43 + comp_r_7_44 + comp_r_7_45 + comp_r_7_46 + comp_r_7_47 + comp_r_7_48 + comp_r_7_49 + comp_r_7_50 + comp_r_7_51 + comp_r_7_52 + comp_r_7_53 + comp_r_7_54 + comp_r_7_55 + comp_r_7_56 + comp_r_7_57 + comp_r_7_58 + comp_r_7_59 + comp_r_7_60 + comp_r_7_61 + comp_r_7_62 + comp_r_7_63 + comp_r_7_64 + comp_r_7_65 + comp_r_7_66 + comp_r_7_67 + comp_r_7_68 + comp_r_7_69 + comp_r_7_70 + comp_r_7_71 + comp_r_7_72 + comp_r_7_73 + comp_r_7_74 + comp_r_7_75 + comp_r_7_76 + comp_r_7_77 + comp_r_7_78 + comp_r_7_79 + comp_r_7_80 + comp_r_7_81 + comp_r_7_82 + comp_r_7_83 + comp_r_7_84 + comp_r_7_85 + comp_r_7_86 + comp_r_7_87 + comp_r_7_88 + comp_r_7_89 + comp_r_7_90 + comp_r_7_91 + comp_r_7_92 + comp_r_7_93 + comp_r_7_94 + comp_r_7_95 + comp_r_7_96 + comp_r_7_97 + comp_r_7_98 + comp_r_7_99 + comp_r_7_100 + comp_r_7_101 + comp_r_7_102 + comp_r_7_103 + comp_r_7_104 + comp_r_7_105 + comp_r_7_106 + comp_r_7_107 + comp_r_7_108 + comp_r_7_109 + comp_r_7_110 + comp_r_7_111 + comp_r_7_112 + comp_r_7_113 + comp_r_7_114 + comp_r_7_115 + comp_r_7_116 + comp_r_7_117 + comp_r_7_118 + comp_r_7_119 + comp_r_7_120 + comp_r_7_121 + comp_r_7_122 + comp_r_7_123 + comp_r_7_124 + comp_r_7_125 + comp_r_7_126;
    position[8] = 4 + (sc_uint<1>) ~(comp_r_1_8) + (sc_uint<1>) ~(comp_r_3_8) + (sc_uint<1>) ~(comp_r_5_8) + (sc_uint<1>) ~(comp_r_7_8);
    position[9] = 5 + (sc_uint<1>) ~(comp_r_1_9) + (sc_uint<1>) ~(comp_r_3_9) + (sc_uint<1>) ~(comp_r_5_9) + (sc_uint<1>) ~(comp_r_7_9) + comp_r_9_10 + comp_r_9_11 + comp_r_9_12 + comp_r_9_13 + comp_r_9_14 + comp_r_9_15 + comp_r_9_16 + comp_r_9_17 + comp_r_9_18 + comp_r_9_19 + comp_r_9_20 + comp_r_9_21 + comp_r_9_22 + comp_r_9_23 + comp_r_9_24 + comp_r_9_25 + comp_r_9_26 + comp_r_9_27 + comp_r_9_28 + comp_r_9_29 + comp_r_9_30 + comp_r_9_31 + comp_r_9_32 + comp_r_9_33 + comp_r_9_34 + comp_r_9_35 + comp_r_9_36 + comp_r_9_37 + comp_r_9_38 + comp_r_9_39 + comp_r_9_40 + comp_r_9_41 + comp_r_9_42 + comp_r_9_43 + comp_r_9_44 + comp_r_9_45 + comp_r_9_46 + comp_r_9_47 + comp_r_9_48 + comp_r_9_49 + comp_r_9_50 + comp_r_9_51 + comp_r_9_52 + comp_r_9_53 + comp_r_9_54 + comp_r_9_55 + comp_r_9_56 + comp_r_9_57 + comp_r_9_58 + comp_r_9_59 + comp_r_9_60 + comp_r_9_61 + comp_r_9_62 + comp_r_9_63 + comp_r_9_64 + comp_r_9_65 + comp_r_9_66 + comp_r_9_67 + comp_r_9_68 + comp_r_9_69 + comp_r_9_70 + comp_r_9_71 + comp_r_9_72 + comp_r_9_73 + comp_r_9_74 + comp_r_9_75 + comp_r_9_76 + comp_r_9_77 + comp_r_9_78 + comp_r_9_79 + comp_r_9_80 + comp_r_9_81 + comp_r_9_82 + comp_r_9_83 + comp_r_9_84 + comp_r_9_85 + comp_r_9_86 + comp_r_9_87 + comp_r_9_88 + comp_r_9_89 + comp_r_9_90 + comp_r_9_91 + comp_r_9_92 + comp_r_9_93 + comp_r_9_94 + comp_r_9_95 + comp_r_9_96 + comp_r_9_97 + comp_r_9_98 + comp_r_9_99 + comp_r_9_100 + comp_r_9_101 + comp_r_9_102 + comp_r_9_103 + comp_r_9_104 + comp_r_9_105 + comp_r_9_106 + comp_r_9_107 + comp_r_9_108 + comp_r_9_109 + comp_r_9_110 + comp_r_9_111 + comp_r_9_112 + comp_r_9_113 + comp_r_9_114 + comp_r_9_115 + comp_r_9_116 + comp_r_9_117 + comp_r_9_118 + comp_r_9_119 + comp_r_9_120 + comp_r_9_121 + comp_r_9_122 + comp_r_9_123 + comp_r_9_124 + comp_r_9_125 + comp_r_9_126;
    position[10] = 5 + (sc_uint<1>) ~(comp_r_1_10) + (sc_uint<1>) ~(comp_r_3_10) + (sc_uint<1>) ~(comp_r_5_10) + (sc_uint<1>) ~(comp_r_7_10) + (sc_uint<1>) ~(comp_r_9_10);
    position[11] = 6 + (sc_uint<1>) ~(comp_r_1_11) + (sc_uint<1>) ~(comp_r_3_11) + (sc_uint<1>) ~(comp_r_5_11) + (sc_uint<1>) ~(comp_r_7_11) + (sc_uint<1>) ~(comp_r_9_11) + comp_r_11_12 + comp_r_11_13 + comp_r_11_14 + comp_r_11_15 + comp_r_11_16 + comp_r_11_17 + comp_r_11_18 + comp_r_11_19 + comp_r_11_20 + comp_r_11_21 + comp_r_11_22 + comp_r_11_23 + comp_r_11_24 + comp_r_11_25 + comp_r_11_26 + comp_r_11_27 + comp_r_11_28 + comp_r_11_29 + comp_r_11_30 + comp_r_11_31 + comp_r_11_32 + comp_r_11_33 + comp_r_11_34 + comp_r_11_35 + comp_r_11_36 + comp_r_11_37 + comp_r_11_38 + comp_r_11_39 + comp_r_11_40 + comp_r_11_41 + comp_r_11_42 + comp_r_11_43 + comp_r_11_44 + comp_r_11_45 + comp_r_11_46 + comp_r_11_47 + comp_r_11_48 + comp_r_11_49 + comp_r_11_50 + comp_r_11_51 + comp_r_11_52 + comp_r_11_53 + comp_r_11_54 + comp_r_11_55 + comp_r_11_56 + comp_r_11_57 + comp_r_11_58 + comp_r_11_59 + comp_r_11_60 + comp_r_11_61 + comp_r_11_62 + comp_r_11_63 + comp_r_11_64 + comp_r_11_65 + comp_r_11_66 + comp_r_11_67 + comp_r_11_68 + comp_r_11_69 + comp_r_11_70 + comp_r_11_71 + comp_r_11_72 + comp_r_11_73 + comp_r_11_74 + comp_r_11_75 + comp_r_11_76 + comp_r_11_77 + comp_r_11_78 + comp_r_11_79 + comp_r_11_80 + comp_r_11_81 + comp_r_11_82 + comp_r_11_83 + comp_r_11_84 + comp_r_11_85 + comp_r_11_86 + comp_r_11_87 + comp_r_11_88 + comp_r_11_89 + comp_r_11_90 + comp_r_11_91 + comp_r_11_92 + comp_r_11_93 + comp_r_11_94 + comp_r_11_95 + comp_r_11_96 + comp_r_11_97 + comp_r_11_98 + comp_r_11_99 + comp_r_11_100 + comp_r_11_101 + comp_r_11_102 + comp_r_11_103 + comp_r_11_104 + comp_r_11_105 + comp_r_11_106 + comp_r_11_107 + comp_r_11_108 + comp_r_11_109 + comp_r_11_110 + comp_r_11_111 + comp_r_11_112 + comp_r_11_113 + comp_r_11_114 + comp_r_11_115 + comp_r_11_116 + comp_r_11_117 + comp_r_11_118 + comp_r_11_119 + comp_r_11_120 + comp_r_11_121 + comp_r_11_122 + comp_r_11_123 + comp_r_11_124 + comp_r_11_125 + comp_r_11_126;
    position[12] = 6 + (sc_uint<1>) ~(comp_r_1_12) + (sc_uint<1>) ~(comp_r_3_12) + (sc_uint<1>) ~(comp_r_5_12) + (sc_uint<1>) ~(comp_r_7_12) + (sc_uint<1>) ~(comp_r_9_12) + (sc_uint<1>) ~(comp_r_11_12);
    position[13] = 7 + (sc_uint<1>) ~(comp_r_1_13) + (sc_uint<1>) ~(comp_r_3_13) + (sc_uint<1>) ~(comp_r_5_13) + (sc_uint<1>) ~(comp_r_7_13) + (sc_uint<1>) ~(comp_r_9_13) + (sc_uint<1>) ~(comp_r_11_13) + comp_r_13_14 + comp_r_13_15 + comp_r_13_16 + comp_r_13_17 + comp_r_13_18 + comp_r_13_19 + comp_r_13_20 + comp_r_13_21 + comp_r_13_22 + comp_r_13_23 + comp_r_13_24 + comp_r_13_25 + comp_r_13_26 + comp_r_13_27 + comp_r_13_28 + comp_r_13_29 + comp_r_13_30 + comp_r_13_31 + comp_r_13_32 + comp_r_13_33 + comp_r_13_34 + comp_r_13_35 + comp_r_13_36 + comp_r_13_37 + comp_r_13_38 + comp_r_13_39 + comp_r_13_40 + comp_r_13_41 + comp_r_13_42 + comp_r_13_43 + comp_r_13_44 + comp_r_13_45 + comp_r_13_46 + comp_r_13_47 + comp_r_13_48 + comp_r_13_49 + comp_r_13_50 + comp_r_13_51 + comp_r_13_52 + comp_r_13_53 + comp_r_13_54 + comp_r_13_55 + comp_r_13_56 + comp_r_13_57 + comp_r_13_58 + comp_r_13_59 + comp_r_13_60 + comp_r_13_61 + comp_r_13_62 + comp_r_13_63 + comp_r_13_64 + comp_r_13_65 + comp_r_13_66 + comp_r_13_67 + comp_r_13_68 + comp_r_13_69 + comp_r_13_70 + comp_r_13_71 + comp_r_13_72 + comp_r_13_73 + comp_r_13_74 + comp_r_13_75 + comp_r_13_76 + comp_r_13_77 + comp_r_13_78 + comp_r_13_79 + comp_r_13_80 + comp_r_13_81 + comp_r_13_82 + comp_r_13_83 + comp_r_13_84 + comp_r_13_85 + comp_r_13_86 + comp_r_13_87 + comp_r_13_88 + comp_r_13_89 + comp_r_13_90 + comp_r_13_91 + comp_r_13_92 + comp_r_13_93 + comp_r_13_94 + comp_r_13_95 + comp_r_13_96 + comp_r_13_97 + comp_r_13_98 + comp_r_13_99 + comp_r_13_100 + comp_r_13_101 + comp_r_13_102 + comp_r_13_103 + comp_r_13_104 + comp_r_13_105 + comp_r_13_106 + comp_r_13_107 + comp_r_13_108 + comp_r_13_109 + comp_r_13_110 + comp_r_13_111 + comp_r_13_112 + comp_r_13_113 + comp_r_13_114 + comp_r_13_115 + comp_r_13_116 + comp_r_13_117 + comp_r_13_118 + comp_r_13_119 + comp_r_13_120 + comp_r_13_121 + comp_r_13_122 + comp_r_13_123 + comp_r_13_124 + comp_r_13_125 + comp_r_13_126;
    position[14] = 7 + (sc_uint<1>) ~(comp_r_1_14) + (sc_uint<1>) ~(comp_r_3_14) + (sc_uint<1>) ~(comp_r_5_14) + (sc_uint<1>) ~(comp_r_7_14) + (sc_uint<1>) ~(comp_r_9_14) + (sc_uint<1>) ~(comp_r_11_14) + (sc_uint<1>) ~(comp_r_13_14);
    position[15] = 8 + (sc_uint<1>) ~(comp_r_1_15) + (sc_uint<1>) ~(comp_r_3_15) + (sc_uint<1>) ~(comp_r_5_15) + (sc_uint<1>) ~(comp_r_7_15) + (sc_uint<1>) ~(comp_r_9_15) + (sc_uint<1>) ~(comp_r_11_15) + (sc_uint<1>) ~(comp_r_13_15) + comp_r_15_16 + comp_r_15_17 + comp_r_15_18 + comp_r_15_19 + comp_r_15_20 + comp_r_15_21 + comp_r_15_22 + comp_r_15_23 + comp_r_15_24 + comp_r_15_25 + comp_r_15_26 + comp_r_15_27 + comp_r_15_28 + comp_r_15_29 + comp_r_15_30 + comp_r_15_31 + comp_r_15_32 + comp_r_15_33 + comp_r_15_34 + comp_r_15_35 + comp_r_15_36 + comp_r_15_37 + comp_r_15_38 + comp_r_15_39 + comp_r_15_40 + comp_r_15_41 + comp_r_15_42 + comp_r_15_43 + comp_r_15_44 + comp_r_15_45 + comp_r_15_46 + comp_r_15_47 + comp_r_15_48 + comp_r_15_49 + comp_r_15_50 + comp_r_15_51 + comp_r_15_52 + comp_r_15_53 + comp_r_15_54 + comp_r_15_55 + comp_r_15_56 + comp_r_15_57 + comp_r_15_58 + comp_r_15_59 + comp_r_15_60 + comp_r_15_61 + comp_r_15_62 + comp_r_15_63 + comp_r_15_64 + comp_r_15_65 + comp_r_15_66 + comp_r_15_67 + comp_r_15_68 + comp_r_15_69 + comp_r_15_70 + comp_r_15_71 + comp_r_15_72 + comp_r_15_73 + comp_r_15_74 + comp_r_15_75 + comp_r_15_76 + comp_r_15_77 + comp_r_15_78 + comp_r_15_79 + comp_r_15_80 + comp_r_15_81 + comp_r_15_82 + comp_r_15_83 + comp_r_15_84 + comp_r_15_85 + comp_r_15_86 + comp_r_15_87 + comp_r_15_88 + comp_r_15_89 + comp_r_15_90 + comp_r_15_91 + comp_r_15_92 + comp_r_15_93 + comp_r_15_94 + comp_r_15_95 + comp_r_15_96 + comp_r_15_97 + comp_r_15_98 + comp_r_15_99 + comp_r_15_100 + comp_r_15_101 + comp_r_15_102 + comp_r_15_103 + comp_r_15_104 + comp_r_15_105 + comp_r_15_106 + comp_r_15_107 + comp_r_15_108 + comp_r_15_109 + comp_r_15_110 + comp_r_15_111 + comp_r_15_112 + comp_r_15_113 + comp_r_15_114 + comp_r_15_115 + comp_r_15_116 + comp_r_15_117 + comp_r_15_118 + comp_r_15_119 + comp_r_15_120 + comp_r_15_121 + comp_r_15_122 + comp_r_15_123 + comp_r_15_124 + comp_r_15_125 + comp_r_15_126;
    position[16] = 8 + (sc_uint<1>) ~(comp_r_1_16) + (sc_uint<1>) ~(comp_r_3_16) + (sc_uint<1>) ~(comp_r_5_16) + (sc_uint<1>) ~(comp_r_7_16) + (sc_uint<1>) ~(comp_r_9_16) + (sc_uint<1>) ~(comp_r_11_16) + (sc_uint<1>) ~(comp_r_13_16) + (sc_uint<1>) ~(comp_r_15_16);
    position[17] = 9 + (sc_uint<1>) ~(comp_r_1_17) + (sc_uint<1>) ~(comp_r_3_17) + (sc_uint<1>) ~(comp_r_5_17) + (sc_uint<1>) ~(comp_r_7_17) + (sc_uint<1>) ~(comp_r_9_17) + (sc_uint<1>) ~(comp_r_11_17) + (sc_uint<1>) ~(comp_r_13_17) + (sc_uint<1>) ~(comp_r_15_17) + comp_r_17_18 + comp_r_17_19 + comp_r_17_20 + comp_r_17_21 + comp_r_17_22 + comp_r_17_23 + comp_r_17_24 + comp_r_17_25 + comp_r_17_26 + comp_r_17_27 + comp_r_17_28 + comp_r_17_29 + comp_r_17_30 + comp_r_17_31 + comp_r_17_32 + comp_r_17_33 + comp_r_17_34 + comp_r_17_35 + comp_r_17_36 + comp_r_17_37 + comp_r_17_38 + comp_r_17_39 + comp_r_17_40 + comp_r_17_41 + comp_r_17_42 + comp_r_17_43 + comp_r_17_44 + comp_r_17_45 + comp_r_17_46 + comp_r_17_47 + comp_r_17_48 + comp_r_17_49 + comp_r_17_50 + comp_r_17_51 + comp_r_17_52 + comp_r_17_53 + comp_r_17_54 + comp_r_17_55 + comp_r_17_56 + comp_r_17_57 + comp_r_17_58 + comp_r_17_59 + comp_r_17_60 + comp_r_17_61 + comp_r_17_62 + comp_r_17_63 + comp_r_17_64 + comp_r_17_65 + comp_r_17_66 + comp_r_17_67 + comp_r_17_68 + comp_r_17_69 + comp_r_17_70 + comp_r_17_71 + comp_r_17_72 + comp_r_17_73 + comp_r_17_74 + comp_r_17_75 + comp_r_17_76 + comp_r_17_77 + comp_r_17_78 + comp_r_17_79 + comp_r_17_80 + comp_r_17_81 + comp_r_17_82 + comp_r_17_83 + comp_r_17_84 + comp_r_17_85 + comp_r_17_86 + comp_r_17_87 + comp_r_17_88 + comp_r_17_89 + comp_r_17_90 + comp_r_17_91 + comp_r_17_92 + comp_r_17_93 + comp_r_17_94 + comp_r_17_95 + comp_r_17_96 + comp_r_17_97 + comp_r_17_98 + comp_r_17_99 + comp_r_17_100 + comp_r_17_101 + comp_r_17_102 + comp_r_17_103 + comp_r_17_104 + comp_r_17_105 + comp_r_17_106 + comp_r_17_107 + comp_r_17_108 + comp_r_17_109 + comp_r_17_110 + comp_r_17_111 + comp_r_17_112 + comp_r_17_113 + comp_r_17_114 + comp_r_17_115 + comp_r_17_116 + comp_r_17_117 + comp_r_17_118 + comp_r_17_119 + comp_r_17_120 + comp_r_17_121 + comp_r_17_122 + comp_r_17_123 + comp_r_17_124 + comp_r_17_125 + comp_r_17_126;
    position[18] = 9 + (sc_uint<1>) ~(comp_r_1_18) + (sc_uint<1>) ~(comp_r_3_18) + (sc_uint<1>) ~(comp_r_5_18) + (sc_uint<1>) ~(comp_r_7_18) + (sc_uint<1>) ~(comp_r_9_18) + (sc_uint<1>) ~(comp_r_11_18) + (sc_uint<1>) ~(comp_r_13_18) + (sc_uint<1>) ~(comp_r_15_18) + (sc_uint<1>) ~(comp_r_17_18);
    position[19] = 10 + (sc_uint<1>) ~(comp_r_1_19) + (sc_uint<1>) ~(comp_r_3_19) + (sc_uint<1>) ~(comp_r_5_19) + (sc_uint<1>) ~(comp_r_7_19) + (sc_uint<1>) ~(comp_r_9_19) + (sc_uint<1>) ~(comp_r_11_19) + (sc_uint<1>) ~(comp_r_13_19) + (sc_uint<1>) ~(comp_r_15_19) + (sc_uint<1>) ~(comp_r_17_19) + comp_r_19_20 + comp_r_19_21 + comp_r_19_22 + comp_r_19_23 + comp_r_19_24 + comp_r_19_25 + comp_r_19_26 + comp_r_19_27 + comp_r_19_28 + comp_r_19_29 + comp_r_19_30 + comp_r_19_31 + comp_r_19_32 + comp_r_19_33 + comp_r_19_34 + comp_r_19_35 + comp_r_19_36 + comp_r_19_37 + comp_r_19_38 + comp_r_19_39 + comp_r_19_40 + comp_r_19_41 + comp_r_19_42 + comp_r_19_43 + comp_r_19_44 + comp_r_19_45 + comp_r_19_46 + comp_r_19_47 + comp_r_19_48 + comp_r_19_49 + comp_r_19_50 + comp_r_19_51 + comp_r_19_52 + comp_r_19_53 + comp_r_19_54 + comp_r_19_55 + comp_r_19_56 + comp_r_19_57 + comp_r_19_58 + comp_r_19_59 + comp_r_19_60 + comp_r_19_61 + comp_r_19_62 + comp_r_19_63 + comp_r_19_64 + comp_r_19_65 + comp_r_19_66 + comp_r_19_67 + comp_r_19_68 + comp_r_19_69 + comp_r_19_70 + comp_r_19_71 + comp_r_19_72 + comp_r_19_73 + comp_r_19_74 + comp_r_19_75 + comp_r_19_76 + comp_r_19_77 + comp_r_19_78 + comp_r_19_79 + comp_r_19_80 + comp_r_19_81 + comp_r_19_82 + comp_r_19_83 + comp_r_19_84 + comp_r_19_85 + comp_r_19_86 + comp_r_19_87 + comp_r_19_88 + comp_r_19_89 + comp_r_19_90 + comp_r_19_91 + comp_r_19_92 + comp_r_19_93 + comp_r_19_94 + comp_r_19_95 + comp_r_19_96 + comp_r_19_97 + comp_r_19_98 + comp_r_19_99 + comp_r_19_100 + comp_r_19_101 + comp_r_19_102 + comp_r_19_103 + comp_r_19_104 + comp_r_19_105 + comp_r_19_106 + comp_r_19_107 + comp_r_19_108 + comp_r_19_109 + comp_r_19_110 + comp_r_19_111 + comp_r_19_112 + comp_r_19_113 + comp_r_19_114 + comp_r_19_115 + comp_r_19_116 + comp_r_19_117 + comp_r_19_118 + comp_r_19_119 + comp_r_19_120 + comp_r_19_121 + comp_r_19_122 + comp_r_19_123 + comp_r_19_124 + comp_r_19_125 + comp_r_19_126;
    position[20] = 10 + (sc_uint<1>) ~(comp_r_1_20) + (sc_uint<1>) ~(comp_r_3_20) + (sc_uint<1>) ~(comp_r_5_20) + (sc_uint<1>) ~(comp_r_7_20) + (sc_uint<1>) ~(comp_r_9_20) + (sc_uint<1>) ~(comp_r_11_20) + (sc_uint<1>) ~(comp_r_13_20) + (sc_uint<1>) ~(comp_r_15_20) + (sc_uint<1>) ~(comp_r_17_20) + (sc_uint<1>) ~(comp_r_19_20);
    position[21] = 11 + (sc_uint<1>) ~(comp_r_1_21) + (sc_uint<1>) ~(comp_r_3_21) + (sc_uint<1>) ~(comp_r_5_21) + (sc_uint<1>) ~(comp_r_7_21) + (sc_uint<1>) ~(comp_r_9_21) + (sc_uint<1>) ~(comp_r_11_21) + (sc_uint<1>) ~(comp_r_13_21) + (sc_uint<1>) ~(comp_r_15_21) + (sc_uint<1>) ~(comp_r_17_21) + (sc_uint<1>) ~(comp_r_19_21) + comp_r_21_22 + comp_r_21_23 + comp_r_21_24 + comp_r_21_25 + comp_r_21_26 + comp_r_21_27 + comp_r_21_28 + comp_r_21_29 + comp_r_21_30 + comp_r_21_31 + comp_r_21_32 + comp_r_21_33 + comp_r_21_34 + comp_r_21_35 + comp_r_21_36 + comp_r_21_37 + comp_r_21_38 + comp_r_21_39 + comp_r_21_40 + comp_r_21_41 + comp_r_21_42 + comp_r_21_43 + comp_r_21_44 + comp_r_21_45 + comp_r_21_46 + comp_r_21_47 + comp_r_21_48 + comp_r_21_49 + comp_r_21_50 + comp_r_21_51 + comp_r_21_52 + comp_r_21_53 + comp_r_21_54 + comp_r_21_55 + comp_r_21_56 + comp_r_21_57 + comp_r_21_58 + comp_r_21_59 + comp_r_21_60 + comp_r_21_61 + comp_r_21_62 + comp_r_21_63 + comp_r_21_64 + comp_r_21_65 + comp_r_21_66 + comp_r_21_67 + comp_r_21_68 + comp_r_21_69 + comp_r_21_70 + comp_r_21_71 + comp_r_21_72 + comp_r_21_73 + comp_r_21_74 + comp_r_21_75 + comp_r_21_76 + comp_r_21_77 + comp_r_21_78 + comp_r_21_79 + comp_r_21_80 + comp_r_21_81 + comp_r_21_82 + comp_r_21_83 + comp_r_21_84 + comp_r_21_85 + comp_r_21_86 + comp_r_21_87 + comp_r_21_88 + comp_r_21_89 + comp_r_21_90 + comp_r_21_91 + comp_r_21_92 + comp_r_21_93 + comp_r_21_94 + comp_r_21_95 + comp_r_21_96 + comp_r_21_97 + comp_r_21_98 + comp_r_21_99 + comp_r_21_100 + comp_r_21_101 + comp_r_21_102 + comp_r_21_103 + comp_r_21_104 + comp_r_21_105 + comp_r_21_106 + comp_r_21_107 + comp_r_21_108 + comp_r_21_109 + comp_r_21_110 + comp_r_21_111 + comp_r_21_112 + comp_r_21_113 + comp_r_21_114 + comp_r_21_115 + comp_r_21_116 + comp_r_21_117 + comp_r_21_118 + comp_r_21_119 + comp_r_21_120 + comp_r_21_121 + comp_r_21_122 + comp_r_21_123 + comp_r_21_124 + comp_r_21_125 + comp_r_21_126;
    position[22] = 11 + (sc_uint<1>) ~(comp_r_1_22) + (sc_uint<1>) ~(comp_r_3_22) + (sc_uint<1>) ~(comp_r_5_22) + (sc_uint<1>) ~(comp_r_7_22) + (sc_uint<1>) ~(comp_r_9_22) + (sc_uint<1>) ~(comp_r_11_22) + (sc_uint<1>) ~(comp_r_13_22) + (sc_uint<1>) ~(comp_r_15_22) + (sc_uint<1>) ~(comp_r_17_22) + (sc_uint<1>) ~(comp_r_19_22) + (sc_uint<1>) ~(comp_r_21_22);
    position[23] = 12 + (sc_uint<1>) ~(comp_r_1_23) + (sc_uint<1>) ~(comp_r_3_23) + (sc_uint<1>) ~(comp_r_5_23) + (sc_uint<1>) ~(comp_r_7_23) + (sc_uint<1>) ~(comp_r_9_23) + (sc_uint<1>) ~(comp_r_11_23) + (sc_uint<1>) ~(comp_r_13_23) + (sc_uint<1>) ~(comp_r_15_23) + (sc_uint<1>) ~(comp_r_17_23) + (sc_uint<1>) ~(comp_r_19_23) + (sc_uint<1>) ~(comp_r_21_23) + comp_r_23_24 + comp_r_23_25 + comp_r_23_26 + comp_r_23_27 + comp_r_23_28 + comp_r_23_29 + comp_r_23_30 + comp_r_23_31 + comp_r_23_32 + comp_r_23_33 + comp_r_23_34 + comp_r_23_35 + comp_r_23_36 + comp_r_23_37 + comp_r_23_38 + comp_r_23_39 + comp_r_23_40 + comp_r_23_41 + comp_r_23_42 + comp_r_23_43 + comp_r_23_44 + comp_r_23_45 + comp_r_23_46 + comp_r_23_47 + comp_r_23_48 + comp_r_23_49 + comp_r_23_50 + comp_r_23_51 + comp_r_23_52 + comp_r_23_53 + comp_r_23_54 + comp_r_23_55 + comp_r_23_56 + comp_r_23_57 + comp_r_23_58 + comp_r_23_59 + comp_r_23_60 + comp_r_23_61 + comp_r_23_62 + comp_r_23_63 + comp_r_23_64 + comp_r_23_65 + comp_r_23_66 + comp_r_23_67 + comp_r_23_68 + comp_r_23_69 + comp_r_23_70 + comp_r_23_71 + comp_r_23_72 + comp_r_23_73 + comp_r_23_74 + comp_r_23_75 + comp_r_23_76 + comp_r_23_77 + comp_r_23_78 + comp_r_23_79 + comp_r_23_80 + comp_r_23_81 + comp_r_23_82 + comp_r_23_83 + comp_r_23_84 + comp_r_23_85 + comp_r_23_86 + comp_r_23_87 + comp_r_23_88 + comp_r_23_89 + comp_r_23_90 + comp_r_23_91 + comp_r_23_92 + comp_r_23_93 + comp_r_23_94 + comp_r_23_95 + comp_r_23_96 + comp_r_23_97 + comp_r_23_98 + comp_r_23_99 + comp_r_23_100 + comp_r_23_101 + comp_r_23_102 + comp_r_23_103 + comp_r_23_104 + comp_r_23_105 + comp_r_23_106 + comp_r_23_107 + comp_r_23_108 + comp_r_23_109 + comp_r_23_110 + comp_r_23_111 + comp_r_23_112 + comp_r_23_113 + comp_r_23_114 + comp_r_23_115 + comp_r_23_116 + comp_r_23_117 + comp_r_23_118 + comp_r_23_119 + comp_r_23_120 + comp_r_23_121 + comp_r_23_122 + comp_r_23_123 + comp_r_23_124 + comp_r_23_125 + comp_r_23_126;
    position[24] = 12 + (sc_uint<1>) ~(comp_r_1_24) + (sc_uint<1>) ~(comp_r_3_24) + (sc_uint<1>) ~(comp_r_5_24) + (sc_uint<1>) ~(comp_r_7_24) + (sc_uint<1>) ~(comp_r_9_24) + (sc_uint<1>) ~(comp_r_11_24) + (sc_uint<1>) ~(comp_r_13_24) + (sc_uint<1>) ~(comp_r_15_24) + (sc_uint<1>) ~(comp_r_17_24) + (sc_uint<1>) ~(comp_r_19_24) + (sc_uint<1>) ~(comp_r_21_24) + (sc_uint<1>) ~(comp_r_23_24);
    position[25] = 13 + (sc_uint<1>) ~(comp_r_1_25) + (sc_uint<1>) ~(comp_r_3_25) + (sc_uint<1>) ~(comp_r_5_25) + (sc_uint<1>) ~(comp_r_7_25) + (sc_uint<1>) ~(comp_r_9_25) + (sc_uint<1>) ~(comp_r_11_25) + (sc_uint<1>) ~(comp_r_13_25) + (sc_uint<1>) ~(comp_r_15_25) + (sc_uint<1>) ~(comp_r_17_25) + (sc_uint<1>) ~(comp_r_19_25) + (sc_uint<1>) ~(comp_r_21_25) + (sc_uint<1>) ~(comp_r_23_25) + comp_r_25_26 + comp_r_25_27 + comp_r_25_28 + comp_r_25_29 + comp_r_25_30 + comp_r_25_31 + comp_r_25_32 + comp_r_25_33 + comp_r_25_34 + comp_r_25_35 + comp_r_25_36 + comp_r_25_37 + comp_r_25_38 + comp_r_25_39 + comp_r_25_40 + comp_r_25_41 + comp_r_25_42 + comp_r_25_43 + comp_r_25_44 + comp_r_25_45 + comp_r_25_46 + comp_r_25_47 + comp_r_25_48 + comp_r_25_49 + comp_r_25_50 + comp_r_25_51 + comp_r_25_52 + comp_r_25_53 + comp_r_25_54 + comp_r_25_55 + comp_r_25_56 + comp_r_25_57 + comp_r_25_58 + comp_r_25_59 + comp_r_25_60 + comp_r_25_61 + comp_r_25_62 + comp_r_25_63 + comp_r_25_64 + comp_r_25_65 + comp_r_25_66 + comp_r_25_67 + comp_r_25_68 + comp_r_25_69 + comp_r_25_70 + comp_r_25_71 + comp_r_25_72 + comp_r_25_73 + comp_r_25_74 + comp_r_25_75 + comp_r_25_76 + comp_r_25_77 + comp_r_25_78 + comp_r_25_79 + comp_r_25_80 + comp_r_25_81 + comp_r_25_82 + comp_r_25_83 + comp_r_25_84 + comp_r_25_85 + comp_r_25_86 + comp_r_25_87 + comp_r_25_88 + comp_r_25_89 + comp_r_25_90 + comp_r_25_91 + comp_r_25_92 + comp_r_25_93 + comp_r_25_94 + comp_r_25_95 + comp_r_25_96 + comp_r_25_97 + comp_r_25_98 + comp_r_25_99 + comp_r_25_100 + comp_r_25_101 + comp_r_25_102 + comp_r_25_103 + comp_r_25_104 + comp_r_25_105 + comp_r_25_106 + comp_r_25_107 + comp_r_25_108 + comp_r_25_109 + comp_r_25_110 + comp_r_25_111 + comp_r_25_112 + comp_r_25_113 + comp_r_25_114 + comp_r_25_115 + comp_r_25_116 + comp_r_25_117 + comp_r_25_118 + comp_r_25_119 + comp_r_25_120 + comp_r_25_121 + comp_r_25_122 + comp_r_25_123 + comp_r_25_124 + comp_r_25_125 + comp_r_25_126;
    position[26] = 13 + (sc_uint<1>) ~(comp_r_1_26) + (sc_uint<1>) ~(comp_r_3_26) + (sc_uint<1>) ~(comp_r_5_26) + (sc_uint<1>) ~(comp_r_7_26) + (sc_uint<1>) ~(comp_r_9_26) + (sc_uint<1>) ~(comp_r_11_26) + (sc_uint<1>) ~(comp_r_13_26) + (sc_uint<1>) ~(comp_r_15_26) + (sc_uint<1>) ~(comp_r_17_26) + (sc_uint<1>) ~(comp_r_19_26) + (sc_uint<1>) ~(comp_r_21_26) + (sc_uint<1>) ~(comp_r_23_26) + (sc_uint<1>) ~(comp_r_25_26);
    position[27] = 14 + (sc_uint<1>) ~(comp_r_1_27) + (sc_uint<1>) ~(comp_r_3_27) + (sc_uint<1>) ~(comp_r_5_27) + (sc_uint<1>) ~(comp_r_7_27) + (sc_uint<1>) ~(comp_r_9_27) + (sc_uint<1>) ~(comp_r_11_27) + (sc_uint<1>) ~(comp_r_13_27) + (sc_uint<1>) ~(comp_r_15_27) + (sc_uint<1>) ~(comp_r_17_27) + (sc_uint<1>) ~(comp_r_19_27) + (sc_uint<1>) ~(comp_r_21_27) + (sc_uint<1>) ~(comp_r_23_27) + (sc_uint<1>) ~(comp_r_25_27) + comp_r_27_28 + comp_r_27_29 + comp_r_27_30 + comp_r_27_31 + comp_r_27_32 + comp_r_27_33 + comp_r_27_34 + comp_r_27_35 + comp_r_27_36 + comp_r_27_37 + comp_r_27_38 + comp_r_27_39 + comp_r_27_40 + comp_r_27_41 + comp_r_27_42 + comp_r_27_43 + comp_r_27_44 + comp_r_27_45 + comp_r_27_46 + comp_r_27_47 + comp_r_27_48 + comp_r_27_49 + comp_r_27_50 + comp_r_27_51 + comp_r_27_52 + comp_r_27_53 + comp_r_27_54 + comp_r_27_55 + comp_r_27_56 + comp_r_27_57 + comp_r_27_58 + comp_r_27_59 + comp_r_27_60 + comp_r_27_61 + comp_r_27_62 + comp_r_27_63 + comp_r_27_64 + comp_r_27_65 + comp_r_27_66 + comp_r_27_67 + comp_r_27_68 + comp_r_27_69 + comp_r_27_70 + comp_r_27_71 + comp_r_27_72 + comp_r_27_73 + comp_r_27_74 + comp_r_27_75 + comp_r_27_76 + comp_r_27_77 + comp_r_27_78 + comp_r_27_79 + comp_r_27_80 + comp_r_27_81 + comp_r_27_82 + comp_r_27_83 + comp_r_27_84 + comp_r_27_85 + comp_r_27_86 + comp_r_27_87 + comp_r_27_88 + comp_r_27_89 + comp_r_27_90 + comp_r_27_91 + comp_r_27_92 + comp_r_27_93 + comp_r_27_94 + comp_r_27_95 + comp_r_27_96 + comp_r_27_97 + comp_r_27_98 + comp_r_27_99 + comp_r_27_100 + comp_r_27_101 + comp_r_27_102 + comp_r_27_103 + comp_r_27_104 + comp_r_27_105 + comp_r_27_106 + comp_r_27_107 + comp_r_27_108 + comp_r_27_109 + comp_r_27_110 + comp_r_27_111 + comp_r_27_112 + comp_r_27_113 + comp_r_27_114 + comp_r_27_115 + comp_r_27_116 + comp_r_27_117 + comp_r_27_118 + comp_r_27_119 + comp_r_27_120 + comp_r_27_121 + comp_r_27_122 + comp_r_27_123 + comp_r_27_124 + comp_r_27_125 + comp_r_27_126;
    position[28] = 14 + (sc_uint<1>) ~(comp_r_1_28) + (sc_uint<1>) ~(comp_r_3_28) + (sc_uint<1>) ~(comp_r_5_28) + (sc_uint<1>) ~(comp_r_7_28) + (sc_uint<1>) ~(comp_r_9_28) + (sc_uint<1>) ~(comp_r_11_28) + (sc_uint<1>) ~(comp_r_13_28) + (sc_uint<1>) ~(comp_r_15_28) + (sc_uint<1>) ~(comp_r_17_28) + (sc_uint<1>) ~(comp_r_19_28) + (sc_uint<1>) ~(comp_r_21_28) + (sc_uint<1>) ~(comp_r_23_28) + (sc_uint<1>) ~(comp_r_25_28) + (sc_uint<1>) ~(comp_r_27_28);
    position[29] = 15 + (sc_uint<1>) ~(comp_r_1_29) + (sc_uint<1>) ~(comp_r_3_29) + (sc_uint<1>) ~(comp_r_5_29) + (sc_uint<1>) ~(comp_r_7_29) + (sc_uint<1>) ~(comp_r_9_29) + (sc_uint<1>) ~(comp_r_11_29) + (sc_uint<1>) ~(comp_r_13_29) + (sc_uint<1>) ~(comp_r_15_29) + (sc_uint<1>) ~(comp_r_17_29) + (sc_uint<1>) ~(comp_r_19_29) + (sc_uint<1>) ~(comp_r_21_29) + (sc_uint<1>) ~(comp_r_23_29) + (sc_uint<1>) ~(comp_r_25_29) + (sc_uint<1>) ~(comp_r_27_29) + comp_r_29_30 + comp_r_29_31 + comp_r_29_32 + comp_r_29_33 + comp_r_29_34 + comp_r_29_35 + comp_r_29_36 + comp_r_29_37 + comp_r_29_38 + comp_r_29_39 + comp_r_29_40 + comp_r_29_41 + comp_r_29_42 + comp_r_29_43 + comp_r_29_44 + comp_r_29_45 + comp_r_29_46 + comp_r_29_47 + comp_r_29_48 + comp_r_29_49 + comp_r_29_50 + comp_r_29_51 + comp_r_29_52 + comp_r_29_53 + comp_r_29_54 + comp_r_29_55 + comp_r_29_56 + comp_r_29_57 + comp_r_29_58 + comp_r_29_59 + comp_r_29_60 + comp_r_29_61 + comp_r_29_62 + comp_r_29_63 + comp_r_29_64 + comp_r_29_65 + comp_r_29_66 + comp_r_29_67 + comp_r_29_68 + comp_r_29_69 + comp_r_29_70 + comp_r_29_71 + comp_r_29_72 + comp_r_29_73 + comp_r_29_74 + comp_r_29_75 + comp_r_29_76 + comp_r_29_77 + comp_r_29_78 + comp_r_29_79 + comp_r_29_80 + comp_r_29_81 + comp_r_29_82 + comp_r_29_83 + comp_r_29_84 + comp_r_29_85 + comp_r_29_86 + comp_r_29_87 + comp_r_29_88 + comp_r_29_89 + comp_r_29_90 + comp_r_29_91 + comp_r_29_92 + comp_r_29_93 + comp_r_29_94 + comp_r_29_95 + comp_r_29_96 + comp_r_29_97 + comp_r_29_98 + comp_r_29_99 + comp_r_29_100 + comp_r_29_101 + comp_r_29_102 + comp_r_29_103 + comp_r_29_104 + comp_r_29_105 + comp_r_29_106 + comp_r_29_107 + comp_r_29_108 + comp_r_29_109 + comp_r_29_110 + comp_r_29_111 + comp_r_29_112 + comp_r_29_113 + comp_r_29_114 + comp_r_29_115 + comp_r_29_116 + comp_r_29_117 + comp_r_29_118 + comp_r_29_119 + comp_r_29_120 + comp_r_29_121 + comp_r_29_122 + comp_r_29_123 + comp_r_29_124 + comp_r_29_125 + comp_r_29_126;
    position[30] = 15 + (sc_uint<1>) ~(comp_r_1_30) + (sc_uint<1>) ~(comp_r_3_30) + (sc_uint<1>) ~(comp_r_5_30) + (sc_uint<1>) ~(comp_r_7_30) + (sc_uint<1>) ~(comp_r_9_30) + (sc_uint<1>) ~(comp_r_11_30) + (sc_uint<1>) ~(comp_r_13_30) + (sc_uint<1>) ~(comp_r_15_30) + (sc_uint<1>) ~(comp_r_17_30) + (sc_uint<1>) ~(comp_r_19_30) + (sc_uint<1>) ~(comp_r_21_30) + (sc_uint<1>) ~(comp_r_23_30) + (sc_uint<1>) ~(comp_r_25_30) + (sc_uint<1>) ~(comp_r_27_30) + (sc_uint<1>) ~(comp_r_29_30);
    position[31] = 16 + (sc_uint<1>) ~(comp_r_1_31) + (sc_uint<1>) ~(comp_r_3_31) + (sc_uint<1>) ~(comp_r_5_31) + (sc_uint<1>) ~(comp_r_7_31) + (sc_uint<1>) ~(comp_r_9_31) + (sc_uint<1>) ~(comp_r_11_31) + (sc_uint<1>) ~(comp_r_13_31) + (sc_uint<1>) ~(comp_r_15_31) + (sc_uint<1>) ~(comp_r_17_31) + (sc_uint<1>) ~(comp_r_19_31) + (sc_uint<1>) ~(comp_r_21_31) + (sc_uint<1>) ~(comp_r_23_31) + (sc_uint<1>) ~(comp_r_25_31) + (sc_uint<1>) ~(comp_r_27_31) + (sc_uint<1>) ~(comp_r_29_31) + comp_r_31_32 + comp_r_31_33 + comp_r_31_34 + comp_r_31_35 + comp_r_31_36 + comp_r_31_37 + comp_r_31_38 + comp_r_31_39 + comp_r_31_40 + comp_r_31_41 + comp_r_31_42 + comp_r_31_43 + comp_r_31_44 + comp_r_31_45 + comp_r_31_46 + comp_r_31_47 + comp_r_31_48 + comp_r_31_49 + comp_r_31_50 + comp_r_31_51 + comp_r_31_52 + comp_r_31_53 + comp_r_31_54 + comp_r_31_55 + comp_r_31_56 + comp_r_31_57 + comp_r_31_58 + comp_r_31_59 + comp_r_31_60 + comp_r_31_61 + comp_r_31_62 + comp_r_31_63 + comp_r_31_64 + comp_r_31_65 + comp_r_31_66 + comp_r_31_67 + comp_r_31_68 + comp_r_31_69 + comp_r_31_70 + comp_r_31_71 + comp_r_31_72 + comp_r_31_73 + comp_r_31_74 + comp_r_31_75 + comp_r_31_76 + comp_r_31_77 + comp_r_31_78 + comp_r_31_79 + comp_r_31_80 + comp_r_31_81 + comp_r_31_82 + comp_r_31_83 + comp_r_31_84 + comp_r_31_85 + comp_r_31_86 + comp_r_31_87 + comp_r_31_88 + comp_r_31_89 + comp_r_31_90 + comp_r_31_91 + comp_r_31_92 + comp_r_31_93 + comp_r_31_94 + comp_r_31_95 + comp_r_31_96 + comp_r_31_97 + comp_r_31_98 + comp_r_31_99 + comp_r_31_100 + comp_r_31_101 + comp_r_31_102 + comp_r_31_103 + comp_r_31_104 + comp_r_31_105 + comp_r_31_106 + comp_r_31_107 + comp_r_31_108 + comp_r_31_109 + comp_r_31_110 + comp_r_31_111 + comp_r_31_112 + comp_r_31_113 + comp_r_31_114 + comp_r_31_115 + comp_r_31_116 + comp_r_31_117 + comp_r_31_118 + comp_r_31_119 + comp_r_31_120 + comp_r_31_121 + comp_r_31_122 + comp_r_31_123 + comp_r_31_124 + comp_r_31_125 + comp_r_31_126;
    position[32] = 16 + (sc_uint<1>) ~(comp_r_1_32) + (sc_uint<1>) ~(comp_r_3_32) + (sc_uint<1>) ~(comp_r_5_32) + (sc_uint<1>) ~(comp_r_7_32) + (sc_uint<1>) ~(comp_r_9_32) + (sc_uint<1>) ~(comp_r_11_32) + (sc_uint<1>) ~(comp_r_13_32) + (sc_uint<1>) ~(comp_r_15_32) + (sc_uint<1>) ~(comp_r_17_32) + (sc_uint<1>) ~(comp_r_19_32) + (sc_uint<1>) ~(comp_r_21_32) + (sc_uint<1>) ~(comp_r_23_32) + (sc_uint<1>) ~(comp_r_25_32) + (sc_uint<1>) ~(comp_r_27_32) + (sc_uint<1>) ~(comp_r_29_32) + (sc_uint<1>) ~(comp_r_31_32);
    position[33] = 17 + (sc_uint<1>) ~(comp_r_1_33) + (sc_uint<1>) ~(comp_r_3_33) + (sc_uint<1>) ~(comp_r_5_33) + (sc_uint<1>) ~(comp_r_7_33) + (sc_uint<1>) ~(comp_r_9_33) + (sc_uint<1>) ~(comp_r_11_33) + (sc_uint<1>) ~(comp_r_13_33) + (sc_uint<1>) ~(comp_r_15_33) + (sc_uint<1>) ~(comp_r_17_33) + (sc_uint<1>) ~(comp_r_19_33) + (sc_uint<1>) ~(comp_r_21_33) + (sc_uint<1>) ~(comp_r_23_33) + (sc_uint<1>) ~(comp_r_25_33) + (sc_uint<1>) ~(comp_r_27_33) + (sc_uint<1>) ~(comp_r_29_33) + (sc_uint<1>) ~(comp_r_31_33) + comp_r_33_34 + comp_r_33_35 + comp_r_33_36 + comp_r_33_37 + comp_r_33_38 + comp_r_33_39 + comp_r_33_40 + comp_r_33_41 + comp_r_33_42 + comp_r_33_43 + comp_r_33_44 + comp_r_33_45 + comp_r_33_46 + comp_r_33_47 + comp_r_33_48 + comp_r_33_49 + comp_r_33_50 + comp_r_33_51 + comp_r_33_52 + comp_r_33_53 + comp_r_33_54 + comp_r_33_55 + comp_r_33_56 + comp_r_33_57 + comp_r_33_58 + comp_r_33_59 + comp_r_33_60 + comp_r_33_61 + comp_r_33_62 + comp_r_33_63 + comp_r_33_64 + comp_r_33_65 + comp_r_33_66 + comp_r_33_67 + comp_r_33_68 + comp_r_33_69 + comp_r_33_70 + comp_r_33_71 + comp_r_33_72 + comp_r_33_73 + comp_r_33_74 + comp_r_33_75 + comp_r_33_76 + comp_r_33_77 + comp_r_33_78 + comp_r_33_79 + comp_r_33_80 + comp_r_33_81 + comp_r_33_82 + comp_r_33_83 + comp_r_33_84 + comp_r_33_85 + comp_r_33_86 + comp_r_33_87 + comp_r_33_88 + comp_r_33_89 + comp_r_33_90 + comp_r_33_91 + comp_r_33_92 + comp_r_33_93 + comp_r_33_94 + comp_r_33_95 + comp_r_33_96 + comp_r_33_97 + comp_r_33_98 + comp_r_33_99 + comp_r_33_100 + comp_r_33_101 + comp_r_33_102 + comp_r_33_103 + comp_r_33_104 + comp_r_33_105 + comp_r_33_106 + comp_r_33_107 + comp_r_33_108 + comp_r_33_109 + comp_r_33_110 + comp_r_33_111 + comp_r_33_112 + comp_r_33_113 + comp_r_33_114 + comp_r_33_115 + comp_r_33_116 + comp_r_33_117 + comp_r_33_118 + comp_r_33_119 + comp_r_33_120 + comp_r_33_121 + comp_r_33_122 + comp_r_33_123 + comp_r_33_124 + comp_r_33_125 + comp_r_33_126;
    position[34] = 17 + (sc_uint<1>) ~(comp_r_1_34) + (sc_uint<1>) ~(comp_r_3_34) + (sc_uint<1>) ~(comp_r_5_34) + (sc_uint<1>) ~(comp_r_7_34) + (sc_uint<1>) ~(comp_r_9_34) + (sc_uint<1>) ~(comp_r_11_34) + (sc_uint<1>) ~(comp_r_13_34) + (sc_uint<1>) ~(comp_r_15_34) + (sc_uint<1>) ~(comp_r_17_34) + (sc_uint<1>) ~(comp_r_19_34) + (sc_uint<1>) ~(comp_r_21_34) + (sc_uint<1>) ~(comp_r_23_34) + (sc_uint<1>) ~(comp_r_25_34) + (sc_uint<1>) ~(comp_r_27_34) + (sc_uint<1>) ~(comp_r_29_34) + (sc_uint<1>) ~(comp_r_31_34) + (sc_uint<1>) ~(comp_r_33_34);
    position[35] = 18 + (sc_uint<1>) ~(comp_r_1_35) + (sc_uint<1>) ~(comp_r_3_35) + (sc_uint<1>) ~(comp_r_5_35) + (sc_uint<1>) ~(comp_r_7_35) + (sc_uint<1>) ~(comp_r_9_35) + (sc_uint<1>) ~(comp_r_11_35) + (sc_uint<1>) ~(comp_r_13_35) + (sc_uint<1>) ~(comp_r_15_35) + (sc_uint<1>) ~(comp_r_17_35) + (sc_uint<1>) ~(comp_r_19_35) + (sc_uint<1>) ~(comp_r_21_35) + (sc_uint<1>) ~(comp_r_23_35) + (sc_uint<1>) ~(comp_r_25_35) + (sc_uint<1>) ~(comp_r_27_35) + (sc_uint<1>) ~(comp_r_29_35) + (sc_uint<1>) ~(comp_r_31_35) + (sc_uint<1>) ~(comp_r_33_35) + comp_r_35_36 + comp_r_35_37 + comp_r_35_38 + comp_r_35_39 + comp_r_35_40 + comp_r_35_41 + comp_r_35_42 + comp_r_35_43 + comp_r_35_44 + comp_r_35_45 + comp_r_35_46 + comp_r_35_47 + comp_r_35_48 + comp_r_35_49 + comp_r_35_50 + comp_r_35_51 + comp_r_35_52 + comp_r_35_53 + comp_r_35_54 + comp_r_35_55 + comp_r_35_56 + comp_r_35_57 + comp_r_35_58 + comp_r_35_59 + comp_r_35_60 + comp_r_35_61 + comp_r_35_62 + comp_r_35_63 + comp_r_35_64 + comp_r_35_65 + comp_r_35_66 + comp_r_35_67 + comp_r_35_68 + comp_r_35_69 + comp_r_35_70 + comp_r_35_71 + comp_r_35_72 + comp_r_35_73 + comp_r_35_74 + comp_r_35_75 + comp_r_35_76 + comp_r_35_77 + comp_r_35_78 + comp_r_35_79 + comp_r_35_80 + comp_r_35_81 + comp_r_35_82 + comp_r_35_83 + comp_r_35_84 + comp_r_35_85 + comp_r_35_86 + comp_r_35_87 + comp_r_35_88 + comp_r_35_89 + comp_r_35_90 + comp_r_35_91 + comp_r_35_92 + comp_r_35_93 + comp_r_35_94 + comp_r_35_95 + comp_r_35_96 + comp_r_35_97 + comp_r_35_98 + comp_r_35_99 + comp_r_35_100 + comp_r_35_101 + comp_r_35_102 + comp_r_35_103 + comp_r_35_104 + comp_r_35_105 + comp_r_35_106 + comp_r_35_107 + comp_r_35_108 + comp_r_35_109 + comp_r_35_110 + comp_r_35_111 + comp_r_35_112 + comp_r_35_113 + comp_r_35_114 + comp_r_35_115 + comp_r_35_116 + comp_r_35_117 + comp_r_35_118 + comp_r_35_119 + comp_r_35_120 + comp_r_35_121 + comp_r_35_122 + comp_r_35_123 + comp_r_35_124 + comp_r_35_125 + comp_r_35_126;
    position[36] = 18 + (sc_uint<1>) ~(comp_r_1_36) + (sc_uint<1>) ~(comp_r_3_36) + (sc_uint<1>) ~(comp_r_5_36) + (sc_uint<1>) ~(comp_r_7_36) + (sc_uint<1>) ~(comp_r_9_36) + (sc_uint<1>) ~(comp_r_11_36) + (sc_uint<1>) ~(comp_r_13_36) + (sc_uint<1>) ~(comp_r_15_36) + (sc_uint<1>) ~(comp_r_17_36) + (sc_uint<1>) ~(comp_r_19_36) + (sc_uint<1>) ~(comp_r_21_36) + (sc_uint<1>) ~(comp_r_23_36) + (sc_uint<1>) ~(comp_r_25_36) + (sc_uint<1>) ~(comp_r_27_36) + (sc_uint<1>) ~(comp_r_29_36) + (sc_uint<1>) ~(comp_r_31_36) + (sc_uint<1>) ~(comp_r_33_36) + (sc_uint<1>) ~(comp_r_35_36);
    position[37] = 19 + (sc_uint<1>) ~(comp_r_1_37) + (sc_uint<1>) ~(comp_r_3_37) + (sc_uint<1>) ~(comp_r_5_37) + (sc_uint<1>) ~(comp_r_7_37) + (sc_uint<1>) ~(comp_r_9_37) + (sc_uint<1>) ~(comp_r_11_37) + (sc_uint<1>) ~(comp_r_13_37) + (sc_uint<1>) ~(comp_r_15_37) + (sc_uint<1>) ~(comp_r_17_37) + (sc_uint<1>) ~(comp_r_19_37) + (sc_uint<1>) ~(comp_r_21_37) + (sc_uint<1>) ~(comp_r_23_37) + (sc_uint<1>) ~(comp_r_25_37) + (sc_uint<1>) ~(comp_r_27_37) + (sc_uint<1>) ~(comp_r_29_37) + (sc_uint<1>) ~(comp_r_31_37) + (sc_uint<1>) ~(comp_r_33_37) + (sc_uint<1>) ~(comp_r_35_37) + comp_r_37_38 + comp_r_37_39 + comp_r_37_40 + comp_r_37_41 + comp_r_37_42 + comp_r_37_43 + comp_r_37_44 + comp_r_37_45 + comp_r_37_46 + comp_r_37_47 + comp_r_37_48 + comp_r_37_49 + comp_r_37_50 + comp_r_37_51 + comp_r_37_52 + comp_r_37_53 + comp_r_37_54 + comp_r_37_55 + comp_r_37_56 + comp_r_37_57 + comp_r_37_58 + comp_r_37_59 + comp_r_37_60 + comp_r_37_61 + comp_r_37_62 + comp_r_37_63 + comp_r_37_64 + comp_r_37_65 + comp_r_37_66 + comp_r_37_67 + comp_r_37_68 + comp_r_37_69 + comp_r_37_70 + comp_r_37_71 + comp_r_37_72 + comp_r_37_73 + comp_r_37_74 + comp_r_37_75 + comp_r_37_76 + comp_r_37_77 + comp_r_37_78 + comp_r_37_79 + comp_r_37_80 + comp_r_37_81 + comp_r_37_82 + comp_r_37_83 + comp_r_37_84 + comp_r_37_85 + comp_r_37_86 + comp_r_37_87 + comp_r_37_88 + comp_r_37_89 + comp_r_37_90 + comp_r_37_91 + comp_r_37_92 + comp_r_37_93 + comp_r_37_94 + comp_r_37_95 + comp_r_37_96 + comp_r_37_97 + comp_r_37_98 + comp_r_37_99 + comp_r_37_100 + comp_r_37_101 + comp_r_37_102 + comp_r_37_103 + comp_r_37_104 + comp_r_37_105 + comp_r_37_106 + comp_r_37_107 + comp_r_37_108 + comp_r_37_109 + comp_r_37_110 + comp_r_37_111 + comp_r_37_112 + comp_r_37_113 + comp_r_37_114 + comp_r_37_115 + comp_r_37_116 + comp_r_37_117 + comp_r_37_118 + comp_r_37_119 + comp_r_37_120 + comp_r_37_121 + comp_r_37_122 + comp_r_37_123 + comp_r_37_124 + comp_r_37_125 + comp_r_37_126;
    position[38] = 19 + (sc_uint<1>) ~(comp_r_1_38) + (sc_uint<1>) ~(comp_r_3_38) + (sc_uint<1>) ~(comp_r_5_38) + (sc_uint<1>) ~(comp_r_7_38) + (sc_uint<1>) ~(comp_r_9_38) + (sc_uint<1>) ~(comp_r_11_38) + (sc_uint<1>) ~(comp_r_13_38) + (sc_uint<1>) ~(comp_r_15_38) + (sc_uint<1>) ~(comp_r_17_38) + (sc_uint<1>) ~(comp_r_19_38) + (sc_uint<1>) ~(comp_r_21_38) + (sc_uint<1>) ~(comp_r_23_38) + (sc_uint<1>) ~(comp_r_25_38) + (sc_uint<1>) ~(comp_r_27_38) + (sc_uint<1>) ~(comp_r_29_38) + (sc_uint<1>) ~(comp_r_31_38) + (sc_uint<1>) ~(comp_r_33_38) + (sc_uint<1>) ~(comp_r_35_38) + (sc_uint<1>) ~(comp_r_37_38);
    position[39] = 20 + (sc_uint<1>) ~(comp_r_1_39) + (sc_uint<1>) ~(comp_r_3_39) + (sc_uint<1>) ~(comp_r_5_39) + (sc_uint<1>) ~(comp_r_7_39) + (sc_uint<1>) ~(comp_r_9_39) + (sc_uint<1>) ~(comp_r_11_39) + (sc_uint<1>) ~(comp_r_13_39) + (sc_uint<1>) ~(comp_r_15_39) + (sc_uint<1>) ~(comp_r_17_39) + (sc_uint<1>) ~(comp_r_19_39) + (sc_uint<1>) ~(comp_r_21_39) + (sc_uint<1>) ~(comp_r_23_39) + (sc_uint<1>) ~(comp_r_25_39) + (sc_uint<1>) ~(comp_r_27_39) + (sc_uint<1>) ~(comp_r_29_39) + (sc_uint<1>) ~(comp_r_31_39) + (sc_uint<1>) ~(comp_r_33_39) + (sc_uint<1>) ~(comp_r_35_39) + (sc_uint<1>) ~(comp_r_37_39) + comp_r_39_40 + comp_r_39_41 + comp_r_39_42 + comp_r_39_43 + comp_r_39_44 + comp_r_39_45 + comp_r_39_46 + comp_r_39_47 + comp_r_39_48 + comp_r_39_49 + comp_r_39_50 + comp_r_39_51 + comp_r_39_52 + comp_r_39_53 + comp_r_39_54 + comp_r_39_55 + comp_r_39_56 + comp_r_39_57 + comp_r_39_58 + comp_r_39_59 + comp_r_39_60 + comp_r_39_61 + comp_r_39_62 + comp_r_39_63 + comp_r_39_64 + comp_r_39_65 + comp_r_39_66 + comp_r_39_67 + comp_r_39_68 + comp_r_39_69 + comp_r_39_70 + comp_r_39_71 + comp_r_39_72 + comp_r_39_73 + comp_r_39_74 + comp_r_39_75 + comp_r_39_76 + comp_r_39_77 + comp_r_39_78 + comp_r_39_79 + comp_r_39_80 + comp_r_39_81 + comp_r_39_82 + comp_r_39_83 + comp_r_39_84 + comp_r_39_85 + comp_r_39_86 + comp_r_39_87 + comp_r_39_88 + comp_r_39_89 + comp_r_39_90 + comp_r_39_91 + comp_r_39_92 + comp_r_39_93 + comp_r_39_94 + comp_r_39_95 + comp_r_39_96 + comp_r_39_97 + comp_r_39_98 + comp_r_39_99 + comp_r_39_100 + comp_r_39_101 + comp_r_39_102 + comp_r_39_103 + comp_r_39_104 + comp_r_39_105 + comp_r_39_106 + comp_r_39_107 + comp_r_39_108 + comp_r_39_109 + comp_r_39_110 + comp_r_39_111 + comp_r_39_112 + comp_r_39_113 + comp_r_39_114 + comp_r_39_115 + comp_r_39_116 + comp_r_39_117 + comp_r_39_118 + comp_r_39_119 + comp_r_39_120 + comp_r_39_121 + comp_r_39_122 + comp_r_39_123 + comp_r_39_124 + comp_r_39_125 + comp_r_39_126;
    position[40] = 20 + (sc_uint<1>) ~(comp_r_1_40) + (sc_uint<1>) ~(comp_r_3_40) + (sc_uint<1>) ~(comp_r_5_40) + (sc_uint<1>) ~(comp_r_7_40) + (sc_uint<1>) ~(comp_r_9_40) + (sc_uint<1>) ~(comp_r_11_40) + (sc_uint<1>) ~(comp_r_13_40) + (sc_uint<1>) ~(comp_r_15_40) + (sc_uint<1>) ~(comp_r_17_40) + (sc_uint<1>) ~(comp_r_19_40) + (sc_uint<1>) ~(comp_r_21_40) + (sc_uint<1>) ~(comp_r_23_40) + (sc_uint<1>) ~(comp_r_25_40) + (sc_uint<1>) ~(comp_r_27_40) + (sc_uint<1>) ~(comp_r_29_40) + (sc_uint<1>) ~(comp_r_31_40) + (sc_uint<1>) ~(comp_r_33_40) + (sc_uint<1>) ~(comp_r_35_40) + (sc_uint<1>) ~(comp_r_37_40) + (sc_uint<1>) ~(comp_r_39_40);
    position[41] = 21 + (sc_uint<1>) ~(comp_r_1_41) + (sc_uint<1>) ~(comp_r_3_41) + (sc_uint<1>) ~(comp_r_5_41) + (sc_uint<1>) ~(comp_r_7_41) + (sc_uint<1>) ~(comp_r_9_41) + (sc_uint<1>) ~(comp_r_11_41) + (sc_uint<1>) ~(comp_r_13_41) + (sc_uint<1>) ~(comp_r_15_41) + (sc_uint<1>) ~(comp_r_17_41) + (sc_uint<1>) ~(comp_r_19_41) + (sc_uint<1>) ~(comp_r_21_41) + (sc_uint<1>) ~(comp_r_23_41) + (sc_uint<1>) ~(comp_r_25_41) + (sc_uint<1>) ~(comp_r_27_41) + (sc_uint<1>) ~(comp_r_29_41) + (sc_uint<1>) ~(comp_r_31_41) + (sc_uint<1>) ~(comp_r_33_41) + (sc_uint<1>) ~(comp_r_35_41) + (sc_uint<1>) ~(comp_r_37_41) + (sc_uint<1>) ~(comp_r_39_41) + comp_r_41_42 + comp_r_41_43 + comp_r_41_44 + comp_r_41_45 + comp_r_41_46 + comp_r_41_47 + comp_r_41_48 + comp_r_41_49 + comp_r_41_50 + comp_r_41_51 + comp_r_41_52 + comp_r_41_53 + comp_r_41_54 + comp_r_41_55 + comp_r_41_56 + comp_r_41_57 + comp_r_41_58 + comp_r_41_59 + comp_r_41_60 + comp_r_41_61 + comp_r_41_62 + comp_r_41_63 + comp_r_41_64 + comp_r_41_65 + comp_r_41_66 + comp_r_41_67 + comp_r_41_68 + comp_r_41_69 + comp_r_41_70 + comp_r_41_71 + comp_r_41_72 + comp_r_41_73 + comp_r_41_74 + comp_r_41_75 + comp_r_41_76 + comp_r_41_77 + comp_r_41_78 + comp_r_41_79 + comp_r_41_80 + comp_r_41_81 + comp_r_41_82 + comp_r_41_83 + comp_r_41_84 + comp_r_41_85 + comp_r_41_86 + comp_r_41_87 + comp_r_41_88 + comp_r_41_89 + comp_r_41_90 + comp_r_41_91 + comp_r_41_92 + comp_r_41_93 + comp_r_41_94 + comp_r_41_95 + comp_r_41_96 + comp_r_41_97 + comp_r_41_98 + comp_r_41_99 + comp_r_41_100 + comp_r_41_101 + comp_r_41_102 + comp_r_41_103 + comp_r_41_104 + comp_r_41_105 + comp_r_41_106 + comp_r_41_107 + comp_r_41_108 + comp_r_41_109 + comp_r_41_110 + comp_r_41_111 + comp_r_41_112 + comp_r_41_113 + comp_r_41_114 + comp_r_41_115 + comp_r_41_116 + comp_r_41_117 + comp_r_41_118 + comp_r_41_119 + comp_r_41_120 + comp_r_41_121 + comp_r_41_122 + comp_r_41_123 + comp_r_41_124 + comp_r_41_125 + comp_r_41_126;
    position[42] = 21 + (sc_uint<1>) ~(comp_r_1_42) + (sc_uint<1>) ~(comp_r_3_42) + (sc_uint<1>) ~(comp_r_5_42) + (sc_uint<1>) ~(comp_r_7_42) + (sc_uint<1>) ~(comp_r_9_42) + (sc_uint<1>) ~(comp_r_11_42) + (sc_uint<1>) ~(comp_r_13_42) + (sc_uint<1>) ~(comp_r_15_42) + (sc_uint<1>) ~(comp_r_17_42) + (sc_uint<1>) ~(comp_r_19_42) + (sc_uint<1>) ~(comp_r_21_42) + (sc_uint<1>) ~(comp_r_23_42) + (sc_uint<1>) ~(comp_r_25_42) + (sc_uint<1>) ~(comp_r_27_42) + (sc_uint<1>) ~(comp_r_29_42) + (sc_uint<1>) ~(comp_r_31_42) + (sc_uint<1>) ~(comp_r_33_42) + (sc_uint<1>) ~(comp_r_35_42) + (sc_uint<1>) ~(comp_r_37_42) + (sc_uint<1>) ~(comp_r_39_42) + (sc_uint<1>) ~(comp_r_41_42);
    position[43] = 22 + (sc_uint<1>) ~(comp_r_1_43) + (sc_uint<1>) ~(comp_r_3_43) + (sc_uint<1>) ~(comp_r_5_43) + (sc_uint<1>) ~(comp_r_7_43) + (sc_uint<1>) ~(comp_r_9_43) + (sc_uint<1>) ~(comp_r_11_43) + (sc_uint<1>) ~(comp_r_13_43) + (sc_uint<1>) ~(comp_r_15_43) + (sc_uint<1>) ~(comp_r_17_43) + (sc_uint<1>) ~(comp_r_19_43) + (sc_uint<1>) ~(comp_r_21_43) + (sc_uint<1>) ~(comp_r_23_43) + (sc_uint<1>) ~(comp_r_25_43) + (sc_uint<1>) ~(comp_r_27_43) + (sc_uint<1>) ~(comp_r_29_43) + (sc_uint<1>) ~(comp_r_31_43) + (sc_uint<1>) ~(comp_r_33_43) + (sc_uint<1>) ~(comp_r_35_43) + (sc_uint<1>) ~(comp_r_37_43) + (sc_uint<1>) ~(comp_r_39_43) + (sc_uint<1>) ~(comp_r_41_43) + comp_r_43_44 + comp_r_43_45 + comp_r_43_46 + comp_r_43_47 + comp_r_43_48 + comp_r_43_49 + comp_r_43_50 + comp_r_43_51 + comp_r_43_52 + comp_r_43_53 + comp_r_43_54 + comp_r_43_55 + comp_r_43_56 + comp_r_43_57 + comp_r_43_58 + comp_r_43_59 + comp_r_43_60 + comp_r_43_61 + comp_r_43_62 + comp_r_43_63 + comp_r_43_64 + comp_r_43_65 + comp_r_43_66 + comp_r_43_67 + comp_r_43_68 + comp_r_43_69 + comp_r_43_70 + comp_r_43_71 + comp_r_43_72 + comp_r_43_73 + comp_r_43_74 + comp_r_43_75 + comp_r_43_76 + comp_r_43_77 + comp_r_43_78 + comp_r_43_79 + comp_r_43_80 + comp_r_43_81 + comp_r_43_82 + comp_r_43_83 + comp_r_43_84 + comp_r_43_85 + comp_r_43_86 + comp_r_43_87 + comp_r_43_88 + comp_r_43_89 + comp_r_43_90 + comp_r_43_91 + comp_r_43_92 + comp_r_43_93 + comp_r_43_94 + comp_r_43_95 + comp_r_43_96 + comp_r_43_97 + comp_r_43_98 + comp_r_43_99 + comp_r_43_100 + comp_r_43_101 + comp_r_43_102 + comp_r_43_103 + comp_r_43_104 + comp_r_43_105 + comp_r_43_106 + comp_r_43_107 + comp_r_43_108 + comp_r_43_109 + comp_r_43_110 + comp_r_43_111 + comp_r_43_112 + comp_r_43_113 + comp_r_43_114 + comp_r_43_115 + comp_r_43_116 + comp_r_43_117 + comp_r_43_118 + comp_r_43_119 + comp_r_43_120 + comp_r_43_121 + comp_r_43_122 + comp_r_43_123 + comp_r_43_124 + comp_r_43_125 + comp_r_43_126;
    position[44] = 22 + (sc_uint<1>) ~(comp_r_1_44) + (sc_uint<1>) ~(comp_r_3_44) + (sc_uint<1>) ~(comp_r_5_44) + (sc_uint<1>) ~(comp_r_7_44) + (sc_uint<1>) ~(comp_r_9_44) + (sc_uint<1>) ~(comp_r_11_44) + (sc_uint<1>) ~(comp_r_13_44) + (sc_uint<1>) ~(comp_r_15_44) + (sc_uint<1>) ~(comp_r_17_44) + (sc_uint<1>) ~(comp_r_19_44) + (sc_uint<1>) ~(comp_r_21_44) + (sc_uint<1>) ~(comp_r_23_44) + (sc_uint<1>) ~(comp_r_25_44) + (sc_uint<1>) ~(comp_r_27_44) + (sc_uint<1>) ~(comp_r_29_44) + (sc_uint<1>) ~(comp_r_31_44) + (sc_uint<1>) ~(comp_r_33_44) + (sc_uint<1>) ~(comp_r_35_44) + (sc_uint<1>) ~(comp_r_37_44) + (sc_uint<1>) ~(comp_r_39_44) + (sc_uint<1>) ~(comp_r_41_44) + (sc_uint<1>) ~(comp_r_43_44);
    position[45] = 23 + (sc_uint<1>) ~(comp_r_1_45) + (sc_uint<1>) ~(comp_r_3_45) + (sc_uint<1>) ~(comp_r_5_45) + (sc_uint<1>) ~(comp_r_7_45) + (sc_uint<1>) ~(comp_r_9_45) + (sc_uint<1>) ~(comp_r_11_45) + (sc_uint<1>) ~(comp_r_13_45) + (sc_uint<1>) ~(comp_r_15_45) + (sc_uint<1>) ~(comp_r_17_45) + (sc_uint<1>) ~(comp_r_19_45) + (sc_uint<1>) ~(comp_r_21_45) + (sc_uint<1>) ~(comp_r_23_45) + (sc_uint<1>) ~(comp_r_25_45) + (sc_uint<1>) ~(comp_r_27_45) + (sc_uint<1>) ~(comp_r_29_45) + (sc_uint<1>) ~(comp_r_31_45) + (sc_uint<1>) ~(comp_r_33_45) + (sc_uint<1>) ~(comp_r_35_45) + (sc_uint<1>) ~(comp_r_37_45) + (sc_uint<1>) ~(comp_r_39_45) + (sc_uint<1>) ~(comp_r_41_45) + (sc_uint<1>) ~(comp_r_43_45) + comp_r_45_46 + comp_r_45_47 + comp_r_45_48 + comp_r_45_49 + comp_r_45_50 + comp_r_45_51 + comp_r_45_52 + comp_r_45_53 + comp_r_45_54 + comp_r_45_55 + comp_r_45_56 + comp_r_45_57 + comp_r_45_58 + comp_r_45_59 + comp_r_45_60 + comp_r_45_61 + comp_r_45_62 + comp_r_45_63 + comp_r_45_64 + comp_r_45_65 + comp_r_45_66 + comp_r_45_67 + comp_r_45_68 + comp_r_45_69 + comp_r_45_70 + comp_r_45_71 + comp_r_45_72 + comp_r_45_73 + comp_r_45_74 + comp_r_45_75 + comp_r_45_76 + comp_r_45_77 + comp_r_45_78 + comp_r_45_79 + comp_r_45_80 + comp_r_45_81 + comp_r_45_82 + comp_r_45_83 + comp_r_45_84 + comp_r_45_85 + comp_r_45_86 + comp_r_45_87 + comp_r_45_88 + comp_r_45_89 + comp_r_45_90 + comp_r_45_91 + comp_r_45_92 + comp_r_45_93 + comp_r_45_94 + comp_r_45_95 + comp_r_45_96 + comp_r_45_97 + comp_r_45_98 + comp_r_45_99 + comp_r_45_100 + comp_r_45_101 + comp_r_45_102 + comp_r_45_103 + comp_r_45_104 + comp_r_45_105 + comp_r_45_106 + comp_r_45_107 + comp_r_45_108 + comp_r_45_109 + comp_r_45_110 + comp_r_45_111 + comp_r_45_112 + comp_r_45_113 + comp_r_45_114 + comp_r_45_115 + comp_r_45_116 + comp_r_45_117 + comp_r_45_118 + comp_r_45_119 + comp_r_45_120 + comp_r_45_121 + comp_r_45_122 + comp_r_45_123 + comp_r_45_124 + comp_r_45_125 + comp_r_45_126;
    position[46] = 23 + (sc_uint<1>) ~(comp_r_1_46) + (sc_uint<1>) ~(comp_r_3_46) + (sc_uint<1>) ~(comp_r_5_46) + (sc_uint<1>) ~(comp_r_7_46) + (sc_uint<1>) ~(comp_r_9_46) + (sc_uint<1>) ~(comp_r_11_46) + (sc_uint<1>) ~(comp_r_13_46) + (sc_uint<1>) ~(comp_r_15_46) + (sc_uint<1>) ~(comp_r_17_46) + (sc_uint<1>) ~(comp_r_19_46) + (sc_uint<1>) ~(comp_r_21_46) + (sc_uint<1>) ~(comp_r_23_46) + (sc_uint<1>) ~(comp_r_25_46) + (sc_uint<1>) ~(comp_r_27_46) + (sc_uint<1>) ~(comp_r_29_46) + (sc_uint<1>) ~(comp_r_31_46) + (sc_uint<1>) ~(comp_r_33_46) + (sc_uint<1>) ~(comp_r_35_46) + (sc_uint<1>) ~(comp_r_37_46) + (sc_uint<1>) ~(comp_r_39_46) + (sc_uint<1>) ~(comp_r_41_46) + (sc_uint<1>) ~(comp_r_43_46) + (sc_uint<1>) ~(comp_r_45_46);
    position[47] = 24 + (sc_uint<1>) ~(comp_r_1_47) + (sc_uint<1>) ~(comp_r_3_47) + (sc_uint<1>) ~(comp_r_5_47) + (sc_uint<1>) ~(comp_r_7_47) + (sc_uint<1>) ~(comp_r_9_47) + (sc_uint<1>) ~(comp_r_11_47) + (sc_uint<1>) ~(comp_r_13_47) + (sc_uint<1>) ~(comp_r_15_47) + (sc_uint<1>) ~(comp_r_17_47) + (sc_uint<1>) ~(comp_r_19_47) + (sc_uint<1>) ~(comp_r_21_47) + (sc_uint<1>) ~(comp_r_23_47) + (sc_uint<1>) ~(comp_r_25_47) + (sc_uint<1>) ~(comp_r_27_47) + (sc_uint<1>) ~(comp_r_29_47) + (sc_uint<1>) ~(comp_r_31_47) + (sc_uint<1>) ~(comp_r_33_47) + (sc_uint<1>) ~(comp_r_35_47) + (sc_uint<1>) ~(comp_r_37_47) + (sc_uint<1>) ~(comp_r_39_47) + (sc_uint<1>) ~(comp_r_41_47) + (sc_uint<1>) ~(comp_r_43_47) + (sc_uint<1>) ~(comp_r_45_47) + comp_r_47_48 + comp_r_47_49 + comp_r_47_50 + comp_r_47_51 + comp_r_47_52 + comp_r_47_53 + comp_r_47_54 + comp_r_47_55 + comp_r_47_56 + comp_r_47_57 + comp_r_47_58 + comp_r_47_59 + comp_r_47_60 + comp_r_47_61 + comp_r_47_62 + comp_r_47_63 + comp_r_47_64 + comp_r_47_65 + comp_r_47_66 + comp_r_47_67 + comp_r_47_68 + comp_r_47_69 + comp_r_47_70 + comp_r_47_71 + comp_r_47_72 + comp_r_47_73 + comp_r_47_74 + comp_r_47_75 + comp_r_47_76 + comp_r_47_77 + comp_r_47_78 + comp_r_47_79 + comp_r_47_80 + comp_r_47_81 + comp_r_47_82 + comp_r_47_83 + comp_r_47_84 + comp_r_47_85 + comp_r_47_86 + comp_r_47_87 + comp_r_47_88 + comp_r_47_89 + comp_r_47_90 + comp_r_47_91 + comp_r_47_92 + comp_r_47_93 + comp_r_47_94 + comp_r_47_95 + comp_r_47_96 + comp_r_47_97 + comp_r_47_98 + comp_r_47_99 + comp_r_47_100 + comp_r_47_101 + comp_r_47_102 + comp_r_47_103 + comp_r_47_104 + comp_r_47_105 + comp_r_47_106 + comp_r_47_107 + comp_r_47_108 + comp_r_47_109 + comp_r_47_110 + comp_r_47_111 + comp_r_47_112 + comp_r_47_113 + comp_r_47_114 + comp_r_47_115 + comp_r_47_116 + comp_r_47_117 + comp_r_47_118 + comp_r_47_119 + comp_r_47_120 + comp_r_47_121 + comp_r_47_122 + comp_r_47_123 + comp_r_47_124 + comp_r_47_125 + comp_r_47_126;
    position[48] = 24 + (sc_uint<1>) ~(comp_r_1_48) + (sc_uint<1>) ~(comp_r_3_48) + (sc_uint<1>) ~(comp_r_5_48) + (sc_uint<1>) ~(comp_r_7_48) + (sc_uint<1>) ~(comp_r_9_48) + (sc_uint<1>) ~(comp_r_11_48) + (sc_uint<1>) ~(comp_r_13_48) + (sc_uint<1>) ~(comp_r_15_48) + (sc_uint<1>) ~(comp_r_17_48) + (sc_uint<1>) ~(comp_r_19_48) + (sc_uint<1>) ~(comp_r_21_48) + (sc_uint<1>) ~(comp_r_23_48) + (sc_uint<1>) ~(comp_r_25_48) + (sc_uint<1>) ~(comp_r_27_48) + (sc_uint<1>) ~(comp_r_29_48) + (sc_uint<1>) ~(comp_r_31_48) + (sc_uint<1>) ~(comp_r_33_48) + (sc_uint<1>) ~(comp_r_35_48) + (sc_uint<1>) ~(comp_r_37_48) + (sc_uint<1>) ~(comp_r_39_48) + (sc_uint<1>) ~(comp_r_41_48) + (sc_uint<1>) ~(comp_r_43_48) + (sc_uint<1>) ~(comp_r_45_48) + (sc_uint<1>) ~(comp_r_47_48);
    position[49] = 25 + (sc_uint<1>) ~(comp_r_1_49) + (sc_uint<1>) ~(comp_r_3_49) + (sc_uint<1>) ~(comp_r_5_49) + (sc_uint<1>) ~(comp_r_7_49) + (sc_uint<1>) ~(comp_r_9_49) + (sc_uint<1>) ~(comp_r_11_49) + (sc_uint<1>) ~(comp_r_13_49) + (sc_uint<1>) ~(comp_r_15_49) + (sc_uint<1>) ~(comp_r_17_49) + (sc_uint<1>) ~(comp_r_19_49) + (sc_uint<1>) ~(comp_r_21_49) + (sc_uint<1>) ~(comp_r_23_49) + (sc_uint<1>) ~(comp_r_25_49) + (sc_uint<1>) ~(comp_r_27_49) + (sc_uint<1>) ~(comp_r_29_49) + (sc_uint<1>) ~(comp_r_31_49) + (sc_uint<1>) ~(comp_r_33_49) + (sc_uint<1>) ~(comp_r_35_49) + (sc_uint<1>) ~(comp_r_37_49) + (sc_uint<1>) ~(comp_r_39_49) + (sc_uint<1>) ~(comp_r_41_49) + (sc_uint<1>) ~(comp_r_43_49) + (sc_uint<1>) ~(comp_r_45_49) + (sc_uint<1>) ~(comp_r_47_49) + comp_r_49_50 + comp_r_49_51 + comp_r_49_52 + comp_r_49_53 + comp_r_49_54 + comp_r_49_55 + comp_r_49_56 + comp_r_49_57 + comp_r_49_58 + comp_r_49_59 + comp_r_49_60 + comp_r_49_61 + comp_r_49_62 + comp_r_49_63 + comp_r_49_64 + comp_r_49_65 + comp_r_49_66 + comp_r_49_67 + comp_r_49_68 + comp_r_49_69 + comp_r_49_70 + comp_r_49_71 + comp_r_49_72 + comp_r_49_73 + comp_r_49_74 + comp_r_49_75 + comp_r_49_76 + comp_r_49_77 + comp_r_49_78 + comp_r_49_79 + comp_r_49_80 + comp_r_49_81 + comp_r_49_82 + comp_r_49_83 + comp_r_49_84 + comp_r_49_85 + comp_r_49_86 + comp_r_49_87 + comp_r_49_88 + comp_r_49_89 + comp_r_49_90 + comp_r_49_91 + comp_r_49_92 + comp_r_49_93 + comp_r_49_94 + comp_r_49_95 + comp_r_49_96 + comp_r_49_97 + comp_r_49_98 + comp_r_49_99 + comp_r_49_100 + comp_r_49_101 + comp_r_49_102 + comp_r_49_103 + comp_r_49_104 + comp_r_49_105 + comp_r_49_106 + comp_r_49_107 + comp_r_49_108 + comp_r_49_109 + comp_r_49_110 + comp_r_49_111 + comp_r_49_112 + comp_r_49_113 + comp_r_49_114 + comp_r_49_115 + comp_r_49_116 + comp_r_49_117 + comp_r_49_118 + comp_r_49_119 + comp_r_49_120 + comp_r_49_121 + comp_r_49_122 + comp_r_49_123 + comp_r_49_124 + comp_r_49_125 + comp_r_49_126;
    position[50] = 25 + (sc_uint<1>) ~(comp_r_1_50) + (sc_uint<1>) ~(comp_r_3_50) + (sc_uint<1>) ~(comp_r_5_50) + (sc_uint<1>) ~(comp_r_7_50) + (sc_uint<1>) ~(comp_r_9_50) + (sc_uint<1>) ~(comp_r_11_50) + (sc_uint<1>) ~(comp_r_13_50) + (sc_uint<1>) ~(comp_r_15_50) + (sc_uint<1>) ~(comp_r_17_50) + (sc_uint<1>) ~(comp_r_19_50) + (sc_uint<1>) ~(comp_r_21_50) + (sc_uint<1>) ~(comp_r_23_50) + (sc_uint<1>) ~(comp_r_25_50) + (sc_uint<1>) ~(comp_r_27_50) + (sc_uint<1>) ~(comp_r_29_50) + (sc_uint<1>) ~(comp_r_31_50) + (sc_uint<1>) ~(comp_r_33_50) + (sc_uint<1>) ~(comp_r_35_50) + (sc_uint<1>) ~(comp_r_37_50) + (sc_uint<1>) ~(comp_r_39_50) + (sc_uint<1>) ~(comp_r_41_50) + (sc_uint<1>) ~(comp_r_43_50) + (sc_uint<1>) ~(comp_r_45_50) + (sc_uint<1>) ~(comp_r_47_50) + (sc_uint<1>) ~(comp_r_49_50);
    position[51] = 26 + (sc_uint<1>) ~(comp_r_1_51) + (sc_uint<1>) ~(comp_r_3_51) + (sc_uint<1>) ~(comp_r_5_51) + (sc_uint<1>) ~(comp_r_7_51) + (sc_uint<1>) ~(comp_r_9_51) + (sc_uint<1>) ~(comp_r_11_51) + (sc_uint<1>) ~(comp_r_13_51) + (sc_uint<1>) ~(comp_r_15_51) + (sc_uint<1>) ~(comp_r_17_51) + (sc_uint<1>) ~(comp_r_19_51) + (sc_uint<1>) ~(comp_r_21_51) + (sc_uint<1>) ~(comp_r_23_51) + (sc_uint<1>) ~(comp_r_25_51) + (sc_uint<1>) ~(comp_r_27_51) + (sc_uint<1>) ~(comp_r_29_51) + (sc_uint<1>) ~(comp_r_31_51) + (sc_uint<1>) ~(comp_r_33_51) + (sc_uint<1>) ~(comp_r_35_51) + (sc_uint<1>) ~(comp_r_37_51) + (sc_uint<1>) ~(comp_r_39_51) + (sc_uint<1>) ~(comp_r_41_51) + (sc_uint<1>) ~(comp_r_43_51) + (sc_uint<1>) ~(comp_r_45_51) + (sc_uint<1>) ~(comp_r_47_51) + (sc_uint<1>) ~(comp_r_49_51) + comp_r_51_52 + comp_r_51_53 + comp_r_51_54 + comp_r_51_55 + comp_r_51_56 + comp_r_51_57 + comp_r_51_58 + comp_r_51_59 + comp_r_51_60 + comp_r_51_61 + comp_r_51_62 + comp_r_51_63 + comp_r_51_64 + comp_r_51_65 + comp_r_51_66 + comp_r_51_67 + comp_r_51_68 + comp_r_51_69 + comp_r_51_70 + comp_r_51_71 + comp_r_51_72 + comp_r_51_73 + comp_r_51_74 + comp_r_51_75 + comp_r_51_76 + comp_r_51_77 + comp_r_51_78 + comp_r_51_79 + comp_r_51_80 + comp_r_51_81 + comp_r_51_82 + comp_r_51_83 + comp_r_51_84 + comp_r_51_85 + comp_r_51_86 + comp_r_51_87 + comp_r_51_88 + comp_r_51_89 + comp_r_51_90 + comp_r_51_91 + comp_r_51_92 + comp_r_51_93 + comp_r_51_94 + comp_r_51_95 + comp_r_51_96 + comp_r_51_97 + comp_r_51_98 + comp_r_51_99 + comp_r_51_100 + comp_r_51_101 + comp_r_51_102 + comp_r_51_103 + comp_r_51_104 + comp_r_51_105 + comp_r_51_106 + comp_r_51_107 + comp_r_51_108 + comp_r_51_109 + comp_r_51_110 + comp_r_51_111 + comp_r_51_112 + comp_r_51_113 + comp_r_51_114 + comp_r_51_115 + comp_r_51_116 + comp_r_51_117 + comp_r_51_118 + comp_r_51_119 + comp_r_51_120 + comp_r_51_121 + comp_r_51_122 + comp_r_51_123 + comp_r_51_124 + comp_r_51_125 + comp_r_51_126;
    position[52] = 26 + (sc_uint<1>) ~(comp_r_1_52) + (sc_uint<1>) ~(comp_r_3_52) + (sc_uint<1>) ~(comp_r_5_52) + (sc_uint<1>) ~(comp_r_7_52) + (sc_uint<1>) ~(comp_r_9_52) + (sc_uint<1>) ~(comp_r_11_52) + (sc_uint<1>) ~(comp_r_13_52) + (sc_uint<1>) ~(comp_r_15_52) + (sc_uint<1>) ~(comp_r_17_52) + (sc_uint<1>) ~(comp_r_19_52) + (sc_uint<1>) ~(comp_r_21_52) + (sc_uint<1>) ~(comp_r_23_52) + (sc_uint<1>) ~(comp_r_25_52) + (sc_uint<1>) ~(comp_r_27_52) + (sc_uint<1>) ~(comp_r_29_52) + (sc_uint<1>) ~(comp_r_31_52) + (sc_uint<1>) ~(comp_r_33_52) + (sc_uint<1>) ~(comp_r_35_52) + (sc_uint<1>) ~(comp_r_37_52) + (sc_uint<1>) ~(comp_r_39_52) + (sc_uint<1>) ~(comp_r_41_52) + (sc_uint<1>) ~(comp_r_43_52) + (sc_uint<1>) ~(comp_r_45_52) + (sc_uint<1>) ~(comp_r_47_52) + (sc_uint<1>) ~(comp_r_49_52) + (sc_uint<1>) ~(comp_r_51_52);
    position[53] = 27 + (sc_uint<1>) ~(comp_r_1_53) + (sc_uint<1>) ~(comp_r_3_53) + (sc_uint<1>) ~(comp_r_5_53) + (sc_uint<1>) ~(comp_r_7_53) + (sc_uint<1>) ~(comp_r_9_53) + (sc_uint<1>) ~(comp_r_11_53) + (sc_uint<1>) ~(comp_r_13_53) + (sc_uint<1>) ~(comp_r_15_53) + (sc_uint<1>) ~(comp_r_17_53) + (sc_uint<1>) ~(comp_r_19_53) + (sc_uint<1>) ~(comp_r_21_53) + (sc_uint<1>) ~(comp_r_23_53) + (sc_uint<1>) ~(comp_r_25_53) + (sc_uint<1>) ~(comp_r_27_53) + (sc_uint<1>) ~(comp_r_29_53) + (sc_uint<1>) ~(comp_r_31_53) + (sc_uint<1>) ~(comp_r_33_53) + (sc_uint<1>) ~(comp_r_35_53) + (sc_uint<1>) ~(comp_r_37_53) + (sc_uint<1>) ~(comp_r_39_53) + (sc_uint<1>) ~(comp_r_41_53) + (sc_uint<1>) ~(comp_r_43_53) + (sc_uint<1>) ~(comp_r_45_53) + (sc_uint<1>) ~(comp_r_47_53) + (sc_uint<1>) ~(comp_r_49_53) + (sc_uint<1>) ~(comp_r_51_53) + comp_r_53_54 + comp_r_53_55 + comp_r_53_56 + comp_r_53_57 + comp_r_53_58 + comp_r_53_59 + comp_r_53_60 + comp_r_53_61 + comp_r_53_62 + comp_r_53_63 + comp_r_53_64 + comp_r_53_65 + comp_r_53_66 + comp_r_53_67 + comp_r_53_68 + comp_r_53_69 + comp_r_53_70 + comp_r_53_71 + comp_r_53_72 + comp_r_53_73 + comp_r_53_74 + comp_r_53_75 + comp_r_53_76 + comp_r_53_77 + comp_r_53_78 + comp_r_53_79 + comp_r_53_80 + comp_r_53_81 + comp_r_53_82 + comp_r_53_83 + comp_r_53_84 + comp_r_53_85 + comp_r_53_86 + comp_r_53_87 + comp_r_53_88 + comp_r_53_89 + comp_r_53_90 + comp_r_53_91 + comp_r_53_92 + comp_r_53_93 + comp_r_53_94 + comp_r_53_95 + comp_r_53_96 + comp_r_53_97 + comp_r_53_98 + comp_r_53_99 + comp_r_53_100 + comp_r_53_101 + comp_r_53_102 + comp_r_53_103 + comp_r_53_104 + comp_r_53_105 + comp_r_53_106 + comp_r_53_107 + comp_r_53_108 + comp_r_53_109 + comp_r_53_110 + comp_r_53_111 + comp_r_53_112 + comp_r_53_113 + comp_r_53_114 + comp_r_53_115 + comp_r_53_116 + comp_r_53_117 + comp_r_53_118 + comp_r_53_119 + comp_r_53_120 + comp_r_53_121 + comp_r_53_122 + comp_r_53_123 + comp_r_53_124 + comp_r_53_125 + comp_r_53_126;
    position[54] = 27 + (sc_uint<1>) ~(comp_r_1_54) + (sc_uint<1>) ~(comp_r_3_54) + (sc_uint<1>) ~(comp_r_5_54) + (sc_uint<1>) ~(comp_r_7_54) + (sc_uint<1>) ~(comp_r_9_54) + (sc_uint<1>) ~(comp_r_11_54) + (sc_uint<1>) ~(comp_r_13_54) + (sc_uint<1>) ~(comp_r_15_54) + (sc_uint<1>) ~(comp_r_17_54) + (sc_uint<1>) ~(comp_r_19_54) + (sc_uint<1>) ~(comp_r_21_54) + (sc_uint<1>) ~(comp_r_23_54) + (sc_uint<1>) ~(comp_r_25_54) + (sc_uint<1>) ~(comp_r_27_54) + (sc_uint<1>) ~(comp_r_29_54) + (sc_uint<1>) ~(comp_r_31_54) + (sc_uint<1>) ~(comp_r_33_54) + (sc_uint<1>) ~(comp_r_35_54) + (sc_uint<1>) ~(comp_r_37_54) + (sc_uint<1>) ~(comp_r_39_54) + (sc_uint<1>) ~(comp_r_41_54) + (sc_uint<1>) ~(comp_r_43_54) + (sc_uint<1>) ~(comp_r_45_54) + (sc_uint<1>) ~(comp_r_47_54) + (sc_uint<1>) ~(comp_r_49_54) + (sc_uint<1>) ~(comp_r_51_54) + (sc_uint<1>) ~(comp_r_53_54);
    position[55] = 28 + (sc_uint<1>) ~(comp_r_1_55) + (sc_uint<1>) ~(comp_r_3_55) + (sc_uint<1>) ~(comp_r_5_55) + (sc_uint<1>) ~(comp_r_7_55) + (sc_uint<1>) ~(comp_r_9_55) + (sc_uint<1>) ~(comp_r_11_55) + (sc_uint<1>) ~(comp_r_13_55) + (sc_uint<1>) ~(comp_r_15_55) + (sc_uint<1>) ~(comp_r_17_55) + (sc_uint<1>) ~(comp_r_19_55) + (sc_uint<1>) ~(comp_r_21_55) + (sc_uint<1>) ~(comp_r_23_55) + (sc_uint<1>) ~(comp_r_25_55) + (sc_uint<1>) ~(comp_r_27_55) + (sc_uint<1>) ~(comp_r_29_55) + (sc_uint<1>) ~(comp_r_31_55) + (sc_uint<1>) ~(comp_r_33_55) + (sc_uint<1>) ~(comp_r_35_55) + (sc_uint<1>) ~(comp_r_37_55) + (sc_uint<1>) ~(comp_r_39_55) + (sc_uint<1>) ~(comp_r_41_55) + (sc_uint<1>) ~(comp_r_43_55) + (sc_uint<1>) ~(comp_r_45_55) + (sc_uint<1>) ~(comp_r_47_55) + (sc_uint<1>) ~(comp_r_49_55) + (sc_uint<1>) ~(comp_r_51_55) + (sc_uint<1>) ~(comp_r_53_55) + comp_r_55_56 + comp_r_55_57 + comp_r_55_58 + comp_r_55_59 + comp_r_55_60 + comp_r_55_61 + comp_r_55_62 + comp_r_55_63 + comp_r_55_64 + comp_r_55_65 + comp_r_55_66 + comp_r_55_67 + comp_r_55_68 + comp_r_55_69 + comp_r_55_70 + comp_r_55_71 + comp_r_55_72 + comp_r_55_73 + comp_r_55_74 + comp_r_55_75 + comp_r_55_76 + comp_r_55_77 + comp_r_55_78 + comp_r_55_79 + comp_r_55_80 + comp_r_55_81 + comp_r_55_82 + comp_r_55_83 + comp_r_55_84 + comp_r_55_85 + comp_r_55_86 + comp_r_55_87 + comp_r_55_88 + comp_r_55_89 + comp_r_55_90 + comp_r_55_91 + comp_r_55_92 + comp_r_55_93 + comp_r_55_94 + comp_r_55_95 + comp_r_55_96 + comp_r_55_97 + comp_r_55_98 + comp_r_55_99 + comp_r_55_100 + comp_r_55_101 + comp_r_55_102 + comp_r_55_103 + comp_r_55_104 + comp_r_55_105 + comp_r_55_106 + comp_r_55_107 + comp_r_55_108 + comp_r_55_109 + comp_r_55_110 + comp_r_55_111 + comp_r_55_112 + comp_r_55_113 + comp_r_55_114 + comp_r_55_115 + comp_r_55_116 + comp_r_55_117 + comp_r_55_118 + comp_r_55_119 + comp_r_55_120 + comp_r_55_121 + comp_r_55_122 + comp_r_55_123 + comp_r_55_124 + comp_r_55_125 + comp_r_55_126;
    position[56] = 28 + (sc_uint<1>) ~(comp_r_1_56) + (sc_uint<1>) ~(comp_r_3_56) + (sc_uint<1>) ~(comp_r_5_56) + (sc_uint<1>) ~(comp_r_7_56) + (sc_uint<1>) ~(comp_r_9_56) + (sc_uint<1>) ~(comp_r_11_56) + (sc_uint<1>) ~(comp_r_13_56) + (sc_uint<1>) ~(comp_r_15_56) + (sc_uint<1>) ~(comp_r_17_56) + (sc_uint<1>) ~(comp_r_19_56) + (sc_uint<1>) ~(comp_r_21_56) + (sc_uint<1>) ~(comp_r_23_56) + (sc_uint<1>) ~(comp_r_25_56) + (sc_uint<1>) ~(comp_r_27_56) + (sc_uint<1>) ~(comp_r_29_56) + (sc_uint<1>) ~(comp_r_31_56) + (sc_uint<1>) ~(comp_r_33_56) + (sc_uint<1>) ~(comp_r_35_56) + (sc_uint<1>) ~(comp_r_37_56) + (sc_uint<1>) ~(comp_r_39_56) + (sc_uint<1>) ~(comp_r_41_56) + (sc_uint<1>) ~(comp_r_43_56) + (sc_uint<1>) ~(comp_r_45_56) + (sc_uint<1>) ~(comp_r_47_56) + (sc_uint<1>) ~(comp_r_49_56) + (sc_uint<1>) ~(comp_r_51_56) + (sc_uint<1>) ~(comp_r_53_56) + (sc_uint<1>) ~(comp_r_55_56);
    position[57] = 29 + (sc_uint<1>) ~(comp_r_1_57) + (sc_uint<1>) ~(comp_r_3_57) + (sc_uint<1>) ~(comp_r_5_57) + (sc_uint<1>) ~(comp_r_7_57) + (sc_uint<1>) ~(comp_r_9_57) + (sc_uint<1>) ~(comp_r_11_57) + (sc_uint<1>) ~(comp_r_13_57) + (sc_uint<1>) ~(comp_r_15_57) + (sc_uint<1>) ~(comp_r_17_57) + (sc_uint<1>) ~(comp_r_19_57) + (sc_uint<1>) ~(comp_r_21_57) + (sc_uint<1>) ~(comp_r_23_57) + (sc_uint<1>) ~(comp_r_25_57) + (sc_uint<1>) ~(comp_r_27_57) + (sc_uint<1>) ~(comp_r_29_57) + (sc_uint<1>) ~(comp_r_31_57) + (sc_uint<1>) ~(comp_r_33_57) + (sc_uint<1>) ~(comp_r_35_57) + (sc_uint<1>) ~(comp_r_37_57) + (sc_uint<1>) ~(comp_r_39_57) + (sc_uint<1>) ~(comp_r_41_57) + (sc_uint<1>) ~(comp_r_43_57) + (sc_uint<1>) ~(comp_r_45_57) + (sc_uint<1>) ~(comp_r_47_57) + (sc_uint<1>) ~(comp_r_49_57) + (sc_uint<1>) ~(comp_r_51_57) + (sc_uint<1>) ~(comp_r_53_57) + (sc_uint<1>) ~(comp_r_55_57) + comp_r_57_58 + comp_r_57_59 + comp_r_57_60 + comp_r_57_61 + comp_r_57_62 + comp_r_57_63 + comp_r_57_64 + comp_r_57_65 + comp_r_57_66 + comp_r_57_67 + comp_r_57_68 + comp_r_57_69 + comp_r_57_70 + comp_r_57_71 + comp_r_57_72 + comp_r_57_73 + comp_r_57_74 + comp_r_57_75 + comp_r_57_76 + comp_r_57_77 + comp_r_57_78 + comp_r_57_79 + comp_r_57_80 + comp_r_57_81 + comp_r_57_82 + comp_r_57_83 + comp_r_57_84 + comp_r_57_85 + comp_r_57_86 + comp_r_57_87 + comp_r_57_88 + comp_r_57_89 + comp_r_57_90 + comp_r_57_91 + comp_r_57_92 + comp_r_57_93 + comp_r_57_94 + comp_r_57_95 + comp_r_57_96 + comp_r_57_97 + comp_r_57_98 + comp_r_57_99 + comp_r_57_100 + comp_r_57_101 + comp_r_57_102 + comp_r_57_103 + comp_r_57_104 + comp_r_57_105 + comp_r_57_106 + comp_r_57_107 + comp_r_57_108 + comp_r_57_109 + comp_r_57_110 + comp_r_57_111 + comp_r_57_112 + comp_r_57_113 + comp_r_57_114 + comp_r_57_115 + comp_r_57_116 + comp_r_57_117 + comp_r_57_118 + comp_r_57_119 + comp_r_57_120 + comp_r_57_121 + comp_r_57_122 + comp_r_57_123 + comp_r_57_124 + comp_r_57_125 + comp_r_57_126;
    position[58] = 29 + (sc_uint<1>) ~(comp_r_1_58) + (sc_uint<1>) ~(comp_r_3_58) + (sc_uint<1>) ~(comp_r_5_58) + (sc_uint<1>) ~(comp_r_7_58) + (sc_uint<1>) ~(comp_r_9_58) + (sc_uint<1>) ~(comp_r_11_58) + (sc_uint<1>) ~(comp_r_13_58) + (sc_uint<1>) ~(comp_r_15_58) + (sc_uint<1>) ~(comp_r_17_58) + (sc_uint<1>) ~(comp_r_19_58) + (sc_uint<1>) ~(comp_r_21_58) + (sc_uint<1>) ~(comp_r_23_58) + (sc_uint<1>) ~(comp_r_25_58) + (sc_uint<1>) ~(comp_r_27_58) + (sc_uint<1>) ~(comp_r_29_58) + (sc_uint<1>) ~(comp_r_31_58) + (sc_uint<1>) ~(comp_r_33_58) + (sc_uint<1>) ~(comp_r_35_58) + (sc_uint<1>) ~(comp_r_37_58) + (sc_uint<1>) ~(comp_r_39_58) + (sc_uint<1>) ~(comp_r_41_58) + (sc_uint<1>) ~(comp_r_43_58) + (sc_uint<1>) ~(comp_r_45_58) + (sc_uint<1>) ~(comp_r_47_58) + (sc_uint<1>) ~(comp_r_49_58) + (sc_uint<1>) ~(comp_r_51_58) + (sc_uint<1>) ~(comp_r_53_58) + (sc_uint<1>) ~(comp_r_55_58) + (sc_uint<1>) ~(comp_r_57_58);
    position[59] = 30 + (sc_uint<1>) ~(comp_r_1_59) + (sc_uint<1>) ~(comp_r_3_59) + (sc_uint<1>) ~(comp_r_5_59) + (sc_uint<1>) ~(comp_r_7_59) + (sc_uint<1>) ~(comp_r_9_59) + (sc_uint<1>) ~(comp_r_11_59) + (sc_uint<1>) ~(comp_r_13_59) + (sc_uint<1>) ~(comp_r_15_59) + (sc_uint<1>) ~(comp_r_17_59) + (sc_uint<1>) ~(comp_r_19_59) + (sc_uint<1>) ~(comp_r_21_59) + (sc_uint<1>) ~(comp_r_23_59) + (sc_uint<1>) ~(comp_r_25_59) + (sc_uint<1>) ~(comp_r_27_59) + (sc_uint<1>) ~(comp_r_29_59) + (sc_uint<1>) ~(comp_r_31_59) + (sc_uint<1>) ~(comp_r_33_59) + (sc_uint<1>) ~(comp_r_35_59) + (sc_uint<1>) ~(comp_r_37_59) + (sc_uint<1>) ~(comp_r_39_59) + (sc_uint<1>) ~(comp_r_41_59) + (sc_uint<1>) ~(comp_r_43_59) + (sc_uint<1>) ~(comp_r_45_59) + (sc_uint<1>) ~(comp_r_47_59) + (sc_uint<1>) ~(comp_r_49_59) + (sc_uint<1>) ~(comp_r_51_59) + (sc_uint<1>) ~(comp_r_53_59) + (sc_uint<1>) ~(comp_r_55_59) + (sc_uint<1>) ~(comp_r_57_59) + comp_r_59_60 + comp_r_59_61 + comp_r_59_62 + comp_r_59_63 + comp_r_59_64 + comp_r_59_65 + comp_r_59_66 + comp_r_59_67 + comp_r_59_68 + comp_r_59_69 + comp_r_59_70 + comp_r_59_71 + comp_r_59_72 + comp_r_59_73 + comp_r_59_74 + comp_r_59_75 + comp_r_59_76 + comp_r_59_77 + comp_r_59_78 + comp_r_59_79 + comp_r_59_80 + comp_r_59_81 + comp_r_59_82 + comp_r_59_83 + comp_r_59_84 + comp_r_59_85 + comp_r_59_86 + comp_r_59_87 + comp_r_59_88 + comp_r_59_89 + comp_r_59_90 + comp_r_59_91 + comp_r_59_92 + comp_r_59_93 + comp_r_59_94 + comp_r_59_95 + comp_r_59_96 + comp_r_59_97 + comp_r_59_98 + comp_r_59_99 + comp_r_59_100 + comp_r_59_101 + comp_r_59_102 + comp_r_59_103 + comp_r_59_104 + comp_r_59_105 + comp_r_59_106 + comp_r_59_107 + comp_r_59_108 + comp_r_59_109 + comp_r_59_110 + comp_r_59_111 + comp_r_59_112 + comp_r_59_113 + comp_r_59_114 + comp_r_59_115 + comp_r_59_116 + comp_r_59_117 + comp_r_59_118 + comp_r_59_119 + comp_r_59_120 + comp_r_59_121 + comp_r_59_122 + comp_r_59_123 + comp_r_59_124 + comp_r_59_125 + comp_r_59_126;
    position[60] = 30 + (sc_uint<1>) ~(comp_r_1_60) + (sc_uint<1>) ~(comp_r_3_60) + (sc_uint<1>) ~(comp_r_5_60) + (sc_uint<1>) ~(comp_r_7_60) + (sc_uint<1>) ~(comp_r_9_60) + (sc_uint<1>) ~(comp_r_11_60) + (sc_uint<1>) ~(comp_r_13_60) + (sc_uint<1>) ~(comp_r_15_60) + (sc_uint<1>) ~(comp_r_17_60) + (sc_uint<1>) ~(comp_r_19_60) + (sc_uint<1>) ~(comp_r_21_60) + (sc_uint<1>) ~(comp_r_23_60) + (sc_uint<1>) ~(comp_r_25_60) + (sc_uint<1>) ~(comp_r_27_60) + (sc_uint<1>) ~(comp_r_29_60) + (sc_uint<1>) ~(comp_r_31_60) + (sc_uint<1>) ~(comp_r_33_60) + (sc_uint<1>) ~(comp_r_35_60) + (sc_uint<1>) ~(comp_r_37_60) + (sc_uint<1>) ~(comp_r_39_60) + (sc_uint<1>) ~(comp_r_41_60) + (sc_uint<1>) ~(comp_r_43_60) + (sc_uint<1>) ~(comp_r_45_60) + (sc_uint<1>) ~(comp_r_47_60) + (sc_uint<1>) ~(comp_r_49_60) + (sc_uint<1>) ~(comp_r_51_60) + (sc_uint<1>) ~(comp_r_53_60) + (sc_uint<1>) ~(comp_r_55_60) + (sc_uint<1>) ~(comp_r_57_60) + (sc_uint<1>) ~(comp_r_59_60);
    position[61] = 31 + (sc_uint<1>) ~(comp_r_1_61) + (sc_uint<1>) ~(comp_r_3_61) + (sc_uint<1>) ~(comp_r_5_61) + (sc_uint<1>) ~(comp_r_7_61) + (sc_uint<1>) ~(comp_r_9_61) + (sc_uint<1>) ~(comp_r_11_61) + (sc_uint<1>) ~(comp_r_13_61) + (sc_uint<1>) ~(comp_r_15_61) + (sc_uint<1>) ~(comp_r_17_61) + (sc_uint<1>) ~(comp_r_19_61) + (sc_uint<1>) ~(comp_r_21_61) + (sc_uint<1>) ~(comp_r_23_61) + (sc_uint<1>) ~(comp_r_25_61) + (sc_uint<1>) ~(comp_r_27_61) + (sc_uint<1>) ~(comp_r_29_61) + (sc_uint<1>) ~(comp_r_31_61) + (sc_uint<1>) ~(comp_r_33_61) + (sc_uint<1>) ~(comp_r_35_61) + (sc_uint<1>) ~(comp_r_37_61) + (sc_uint<1>) ~(comp_r_39_61) + (sc_uint<1>) ~(comp_r_41_61) + (sc_uint<1>) ~(comp_r_43_61) + (sc_uint<1>) ~(comp_r_45_61) + (sc_uint<1>) ~(comp_r_47_61) + (sc_uint<1>) ~(comp_r_49_61) + (sc_uint<1>) ~(comp_r_51_61) + (sc_uint<1>) ~(comp_r_53_61) + (sc_uint<1>) ~(comp_r_55_61) + (sc_uint<1>) ~(comp_r_57_61) + (sc_uint<1>) ~(comp_r_59_61) + comp_r_61_62 + comp_r_61_63 + comp_r_61_64 + comp_r_61_65 + comp_r_61_66 + comp_r_61_67 + comp_r_61_68 + comp_r_61_69 + comp_r_61_70 + comp_r_61_71 + comp_r_61_72 + comp_r_61_73 + comp_r_61_74 + comp_r_61_75 + comp_r_61_76 + comp_r_61_77 + comp_r_61_78 + comp_r_61_79 + comp_r_61_80 + comp_r_61_81 + comp_r_61_82 + comp_r_61_83 + comp_r_61_84 + comp_r_61_85 + comp_r_61_86 + comp_r_61_87 + comp_r_61_88 + comp_r_61_89 + comp_r_61_90 + comp_r_61_91 + comp_r_61_92 + comp_r_61_93 + comp_r_61_94 + comp_r_61_95 + comp_r_61_96 + comp_r_61_97 + comp_r_61_98 + comp_r_61_99 + comp_r_61_100 + comp_r_61_101 + comp_r_61_102 + comp_r_61_103 + comp_r_61_104 + comp_r_61_105 + comp_r_61_106 + comp_r_61_107 + comp_r_61_108 + comp_r_61_109 + comp_r_61_110 + comp_r_61_111 + comp_r_61_112 + comp_r_61_113 + comp_r_61_114 + comp_r_61_115 + comp_r_61_116 + comp_r_61_117 + comp_r_61_118 + comp_r_61_119 + comp_r_61_120 + comp_r_61_121 + comp_r_61_122 + comp_r_61_123 + comp_r_61_124 + comp_r_61_125 + comp_r_61_126;
    position[62] = 31 + (sc_uint<1>) ~(comp_r_1_62) + (sc_uint<1>) ~(comp_r_3_62) + (sc_uint<1>) ~(comp_r_5_62) + (sc_uint<1>) ~(comp_r_7_62) + (sc_uint<1>) ~(comp_r_9_62) + (sc_uint<1>) ~(comp_r_11_62) + (sc_uint<1>) ~(comp_r_13_62) + (sc_uint<1>) ~(comp_r_15_62) + (sc_uint<1>) ~(comp_r_17_62) + (sc_uint<1>) ~(comp_r_19_62) + (sc_uint<1>) ~(comp_r_21_62) + (sc_uint<1>) ~(comp_r_23_62) + (sc_uint<1>) ~(comp_r_25_62) + (sc_uint<1>) ~(comp_r_27_62) + (sc_uint<1>) ~(comp_r_29_62) + (sc_uint<1>) ~(comp_r_31_62) + (sc_uint<1>) ~(comp_r_33_62) + (sc_uint<1>) ~(comp_r_35_62) + (sc_uint<1>) ~(comp_r_37_62) + (sc_uint<1>) ~(comp_r_39_62) + (sc_uint<1>) ~(comp_r_41_62) + (sc_uint<1>) ~(comp_r_43_62) + (sc_uint<1>) ~(comp_r_45_62) + (sc_uint<1>) ~(comp_r_47_62) + (sc_uint<1>) ~(comp_r_49_62) + (sc_uint<1>) ~(comp_r_51_62) + (sc_uint<1>) ~(comp_r_53_62) + (sc_uint<1>) ~(comp_r_55_62) + (sc_uint<1>) ~(comp_r_57_62) + (sc_uint<1>) ~(comp_r_59_62) + (sc_uint<1>) ~(comp_r_61_62);
    position[63] = 32 + (sc_uint<1>) ~(comp_r_1_63) + (sc_uint<1>) ~(comp_r_3_63) + (sc_uint<1>) ~(comp_r_5_63) + (sc_uint<1>) ~(comp_r_7_63) + (sc_uint<1>) ~(comp_r_9_63) + (sc_uint<1>) ~(comp_r_11_63) + (sc_uint<1>) ~(comp_r_13_63) + (sc_uint<1>) ~(comp_r_15_63) + (sc_uint<1>) ~(comp_r_17_63) + (sc_uint<1>) ~(comp_r_19_63) + (sc_uint<1>) ~(comp_r_21_63) + (sc_uint<1>) ~(comp_r_23_63) + (sc_uint<1>) ~(comp_r_25_63) + (sc_uint<1>) ~(comp_r_27_63) + (sc_uint<1>) ~(comp_r_29_63) + (sc_uint<1>) ~(comp_r_31_63) + (sc_uint<1>) ~(comp_r_33_63) + (sc_uint<1>) ~(comp_r_35_63) + (sc_uint<1>) ~(comp_r_37_63) + (sc_uint<1>) ~(comp_r_39_63) + (sc_uint<1>) ~(comp_r_41_63) + (sc_uint<1>) ~(comp_r_43_63) + (sc_uint<1>) ~(comp_r_45_63) + (sc_uint<1>) ~(comp_r_47_63) + (sc_uint<1>) ~(comp_r_49_63) + (sc_uint<1>) ~(comp_r_51_63) + (sc_uint<1>) ~(comp_r_53_63) + (sc_uint<1>) ~(comp_r_55_63) + (sc_uint<1>) ~(comp_r_57_63) + (sc_uint<1>) ~(comp_r_59_63) + (sc_uint<1>) ~(comp_r_61_63) + comp_r_63_64 + comp_r_63_65 + comp_r_63_66 + comp_r_63_67 + comp_r_63_68 + comp_r_63_69 + comp_r_63_70 + comp_r_63_71 + comp_r_63_72 + comp_r_63_73 + comp_r_63_74 + comp_r_63_75 + comp_r_63_76 + comp_r_63_77 + comp_r_63_78 + comp_r_63_79 + comp_r_63_80 + comp_r_63_81 + comp_r_63_82 + comp_r_63_83 + comp_r_63_84 + comp_r_63_85 + comp_r_63_86 + comp_r_63_87 + comp_r_63_88 + comp_r_63_89 + comp_r_63_90 + comp_r_63_91 + comp_r_63_92 + comp_r_63_93 + comp_r_63_94 + comp_r_63_95 + comp_r_63_96 + comp_r_63_97 + comp_r_63_98 + comp_r_63_99 + comp_r_63_100 + comp_r_63_101 + comp_r_63_102 + comp_r_63_103 + comp_r_63_104 + comp_r_63_105 + comp_r_63_106 + comp_r_63_107 + comp_r_63_108 + comp_r_63_109 + comp_r_63_110 + comp_r_63_111 + comp_r_63_112 + comp_r_63_113 + comp_r_63_114 + comp_r_63_115 + comp_r_63_116 + comp_r_63_117 + comp_r_63_118 + comp_r_63_119 + comp_r_63_120 + comp_r_63_121 + comp_r_63_122 + comp_r_63_123 + comp_r_63_124 + comp_r_63_125 + comp_r_63_126;
    position[64] = 32 + (sc_uint<1>) ~(comp_r_1_64) + (sc_uint<1>) ~(comp_r_3_64) + (sc_uint<1>) ~(comp_r_5_64) + (sc_uint<1>) ~(comp_r_7_64) + (sc_uint<1>) ~(comp_r_9_64) + (sc_uint<1>) ~(comp_r_11_64) + (sc_uint<1>) ~(comp_r_13_64) + (sc_uint<1>) ~(comp_r_15_64) + (sc_uint<1>) ~(comp_r_17_64) + (sc_uint<1>) ~(comp_r_19_64) + (sc_uint<1>) ~(comp_r_21_64) + (sc_uint<1>) ~(comp_r_23_64) + (sc_uint<1>) ~(comp_r_25_64) + (sc_uint<1>) ~(comp_r_27_64) + (sc_uint<1>) ~(comp_r_29_64) + (sc_uint<1>) ~(comp_r_31_64) + (sc_uint<1>) ~(comp_r_33_64) + (sc_uint<1>) ~(comp_r_35_64) + (sc_uint<1>) ~(comp_r_37_64) + (sc_uint<1>) ~(comp_r_39_64) + (sc_uint<1>) ~(comp_r_41_64) + (sc_uint<1>) ~(comp_r_43_64) + (sc_uint<1>) ~(comp_r_45_64) + (sc_uint<1>) ~(comp_r_47_64) + (sc_uint<1>) ~(comp_r_49_64) + (sc_uint<1>) ~(comp_r_51_64) + (sc_uint<1>) ~(comp_r_53_64) + (sc_uint<1>) ~(comp_r_55_64) + (sc_uint<1>) ~(comp_r_57_64) + (sc_uint<1>) ~(comp_r_59_64) + (sc_uint<1>) ~(comp_r_61_64) + (sc_uint<1>) ~(comp_r_63_64);
    position[65] = 33 + (sc_uint<1>) ~(comp_r_1_65) + (sc_uint<1>) ~(comp_r_3_65) + (sc_uint<1>) ~(comp_r_5_65) + (sc_uint<1>) ~(comp_r_7_65) + (sc_uint<1>) ~(comp_r_9_65) + (sc_uint<1>) ~(comp_r_11_65) + (sc_uint<1>) ~(comp_r_13_65) + (sc_uint<1>) ~(comp_r_15_65) + (sc_uint<1>) ~(comp_r_17_65) + (sc_uint<1>) ~(comp_r_19_65) + (sc_uint<1>) ~(comp_r_21_65) + (sc_uint<1>) ~(comp_r_23_65) + (sc_uint<1>) ~(comp_r_25_65) + (sc_uint<1>) ~(comp_r_27_65) + (sc_uint<1>) ~(comp_r_29_65) + (sc_uint<1>) ~(comp_r_31_65) + (sc_uint<1>) ~(comp_r_33_65) + (sc_uint<1>) ~(comp_r_35_65) + (sc_uint<1>) ~(comp_r_37_65) + (sc_uint<1>) ~(comp_r_39_65) + (sc_uint<1>) ~(comp_r_41_65) + (sc_uint<1>) ~(comp_r_43_65) + (sc_uint<1>) ~(comp_r_45_65) + (sc_uint<1>) ~(comp_r_47_65) + (sc_uint<1>) ~(comp_r_49_65) + (sc_uint<1>) ~(comp_r_51_65) + (sc_uint<1>) ~(comp_r_53_65) + (sc_uint<1>) ~(comp_r_55_65) + (sc_uint<1>) ~(comp_r_57_65) + (sc_uint<1>) ~(comp_r_59_65) + (sc_uint<1>) ~(comp_r_61_65) + (sc_uint<1>) ~(comp_r_63_65) + comp_r_65_66 + comp_r_65_67 + comp_r_65_68 + comp_r_65_69 + comp_r_65_70 + comp_r_65_71 + comp_r_65_72 + comp_r_65_73 + comp_r_65_74 + comp_r_65_75 + comp_r_65_76 + comp_r_65_77 + comp_r_65_78 + comp_r_65_79 + comp_r_65_80 + comp_r_65_81 + comp_r_65_82 + comp_r_65_83 + comp_r_65_84 + comp_r_65_85 + comp_r_65_86 + comp_r_65_87 + comp_r_65_88 + comp_r_65_89 + comp_r_65_90 + comp_r_65_91 + comp_r_65_92 + comp_r_65_93 + comp_r_65_94 + comp_r_65_95 + comp_r_65_96 + comp_r_65_97 + comp_r_65_98 + comp_r_65_99 + comp_r_65_100 + comp_r_65_101 + comp_r_65_102 + comp_r_65_103 + comp_r_65_104 + comp_r_65_105 + comp_r_65_106 + comp_r_65_107 + comp_r_65_108 + comp_r_65_109 + comp_r_65_110 + comp_r_65_111 + comp_r_65_112 + comp_r_65_113 + comp_r_65_114 + comp_r_65_115 + comp_r_65_116 + comp_r_65_117 + comp_r_65_118 + comp_r_65_119 + comp_r_65_120 + comp_r_65_121 + comp_r_65_122 + comp_r_65_123 + comp_r_65_124 + comp_r_65_125 + comp_r_65_126;
    position[66] = 33 + (sc_uint<1>) ~(comp_r_1_66) + (sc_uint<1>) ~(comp_r_3_66) + (sc_uint<1>) ~(comp_r_5_66) + (sc_uint<1>) ~(comp_r_7_66) + (sc_uint<1>) ~(comp_r_9_66) + (sc_uint<1>) ~(comp_r_11_66) + (sc_uint<1>) ~(comp_r_13_66) + (sc_uint<1>) ~(comp_r_15_66) + (sc_uint<1>) ~(comp_r_17_66) + (sc_uint<1>) ~(comp_r_19_66) + (sc_uint<1>) ~(comp_r_21_66) + (sc_uint<1>) ~(comp_r_23_66) + (sc_uint<1>) ~(comp_r_25_66) + (sc_uint<1>) ~(comp_r_27_66) + (sc_uint<1>) ~(comp_r_29_66) + (sc_uint<1>) ~(comp_r_31_66) + (sc_uint<1>) ~(comp_r_33_66) + (sc_uint<1>) ~(comp_r_35_66) + (sc_uint<1>) ~(comp_r_37_66) + (sc_uint<1>) ~(comp_r_39_66) + (sc_uint<1>) ~(comp_r_41_66) + (sc_uint<1>) ~(comp_r_43_66) + (sc_uint<1>) ~(comp_r_45_66) + (sc_uint<1>) ~(comp_r_47_66) + (sc_uint<1>) ~(comp_r_49_66) + (sc_uint<1>) ~(comp_r_51_66) + (sc_uint<1>) ~(comp_r_53_66) + (sc_uint<1>) ~(comp_r_55_66) + (sc_uint<1>) ~(comp_r_57_66) + (sc_uint<1>) ~(comp_r_59_66) + (sc_uint<1>) ~(comp_r_61_66) + (sc_uint<1>) ~(comp_r_63_66) + (sc_uint<1>) ~(comp_r_65_66);
    position[67] = 34 + (sc_uint<1>) ~(comp_r_1_67) + (sc_uint<1>) ~(comp_r_3_67) + (sc_uint<1>) ~(comp_r_5_67) + (sc_uint<1>) ~(comp_r_7_67) + (sc_uint<1>) ~(comp_r_9_67) + (sc_uint<1>) ~(comp_r_11_67) + (sc_uint<1>) ~(comp_r_13_67) + (sc_uint<1>) ~(comp_r_15_67) + (sc_uint<1>) ~(comp_r_17_67) + (sc_uint<1>) ~(comp_r_19_67) + (sc_uint<1>) ~(comp_r_21_67) + (sc_uint<1>) ~(comp_r_23_67) + (sc_uint<1>) ~(comp_r_25_67) + (sc_uint<1>) ~(comp_r_27_67) + (sc_uint<1>) ~(comp_r_29_67) + (sc_uint<1>) ~(comp_r_31_67) + (sc_uint<1>) ~(comp_r_33_67) + (sc_uint<1>) ~(comp_r_35_67) + (sc_uint<1>) ~(comp_r_37_67) + (sc_uint<1>) ~(comp_r_39_67) + (sc_uint<1>) ~(comp_r_41_67) + (sc_uint<1>) ~(comp_r_43_67) + (sc_uint<1>) ~(comp_r_45_67) + (sc_uint<1>) ~(comp_r_47_67) + (sc_uint<1>) ~(comp_r_49_67) + (sc_uint<1>) ~(comp_r_51_67) + (sc_uint<1>) ~(comp_r_53_67) + (sc_uint<1>) ~(comp_r_55_67) + (sc_uint<1>) ~(comp_r_57_67) + (sc_uint<1>) ~(comp_r_59_67) + (sc_uint<1>) ~(comp_r_61_67) + (sc_uint<1>) ~(comp_r_63_67) + (sc_uint<1>) ~(comp_r_65_67) + comp_r_67_68 + comp_r_67_69 + comp_r_67_70 + comp_r_67_71 + comp_r_67_72 + comp_r_67_73 + comp_r_67_74 + comp_r_67_75 + comp_r_67_76 + comp_r_67_77 + comp_r_67_78 + comp_r_67_79 + comp_r_67_80 + comp_r_67_81 + comp_r_67_82 + comp_r_67_83 + comp_r_67_84 + comp_r_67_85 + comp_r_67_86 + comp_r_67_87 + comp_r_67_88 + comp_r_67_89 + comp_r_67_90 + comp_r_67_91 + comp_r_67_92 + comp_r_67_93 + comp_r_67_94 + comp_r_67_95 + comp_r_67_96 + comp_r_67_97 + comp_r_67_98 + comp_r_67_99 + comp_r_67_100 + comp_r_67_101 + comp_r_67_102 + comp_r_67_103 + comp_r_67_104 + comp_r_67_105 + comp_r_67_106 + comp_r_67_107 + comp_r_67_108 + comp_r_67_109 + comp_r_67_110 + comp_r_67_111 + comp_r_67_112 + comp_r_67_113 + comp_r_67_114 + comp_r_67_115 + comp_r_67_116 + comp_r_67_117 + comp_r_67_118 + comp_r_67_119 + comp_r_67_120 + comp_r_67_121 + comp_r_67_122 + comp_r_67_123 + comp_r_67_124 + comp_r_67_125 + comp_r_67_126;
    position[68] = 34 + (sc_uint<1>) ~(comp_r_1_68) + (sc_uint<1>) ~(comp_r_3_68) + (sc_uint<1>) ~(comp_r_5_68) + (sc_uint<1>) ~(comp_r_7_68) + (sc_uint<1>) ~(comp_r_9_68) + (sc_uint<1>) ~(comp_r_11_68) + (sc_uint<1>) ~(comp_r_13_68) + (sc_uint<1>) ~(comp_r_15_68) + (sc_uint<1>) ~(comp_r_17_68) + (sc_uint<1>) ~(comp_r_19_68) + (sc_uint<1>) ~(comp_r_21_68) + (sc_uint<1>) ~(comp_r_23_68) + (sc_uint<1>) ~(comp_r_25_68) + (sc_uint<1>) ~(comp_r_27_68) + (sc_uint<1>) ~(comp_r_29_68) + (sc_uint<1>) ~(comp_r_31_68) + (sc_uint<1>) ~(comp_r_33_68) + (sc_uint<1>) ~(comp_r_35_68) + (sc_uint<1>) ~(comp_r_37_68) + (sc_uint<1>) ~(comp_r_39_68) + (sc_uint<1>) ~(comp_r_41_68) + (sc_uint<1>) ~(comp_r_43_68) + (sc_uint<1>) ~(comp_r_45_68) + (sc_uint<1>) ~(comp_r_47_68) + (sc_uint<1>) ~(comp_r_49_68) + (sc_uint<1>) ~(comp_r_51_68) + (sc_uint<1>) ~(comp_r_53_68) + (sc_uint<1>) ~(comp_r_55_68) + (sc_uint<1>) ~(comp_r_57_68) + (sc_uint<1>) ~(comp_r_59_68) + (sc_uint<1>) ~(comp_r_61_68) + (sc_uint<1>) ~(comp_r_63_68) + (sc_uint<1>) ~(comp_r_65_68) + (sc_uint<1>) ~(comp_r_67_68);
    position[69] = 35 + (sc_uint<1>) ~(comp_r_1_69) + (sc_uint<1>) ~(comp_r_3_69) + (sc_uint<1>) ~(comp_r_5_69) + (sc_uint<1>) ~(comp_r_7_69) + (sc_uint<1>) ~(comp_r_9_69) + (sc_uint<1>) ~(comp_r_11_69) + (sc_uint<1>) ~(comp_r_13_69) + (sc_uint<1>) ~(comp_r_15_69) + (sc_uint<1>) ~(comp_r_17_69) + (sc_uint<1>) ~(comp_r_19_69) + (sc_uint<1>) ~(comp_r_21_69) + (sc_uint<1>) ~(comp_r_23_69) + (sc_uint<1>) ~(comp_r_25_69) + (sc_uint<1>) ~(comp_r_27_69) + (sc_uint<1>) ~(comp_r_29_69) + (sc_uint<1>) ~(comp_r_31_69) + (sc_uint<1>) ~(comp_r_33_69) + (sc_uint<1>) ~(comp_r_35_69) + (sc_uint<1>) ~(comp_r_37_69) + (sc_uint<1>) ~(comp_r_39_69) + (sc_uint<1>) ~(comp_r_41_69) + (sc_uint<1>) ~(comp_r_43_69) + (sc_uint<1>) ~(comp_r_45_69) + (sc_uint<1>) ~(comp_r_47_69) + (sc_uint<1>) ~(comp_r_49_69) + (sc_uint<1>) ~(comp_r_51_69) + (sc_uint<1>) ~(comp_r_53_69) + (sc_uint<1>) ~(comp_r_55_69) + (sc_uint<1>) ~(comp_r_57_69) + (sc_uint<1>) ~(comp_r_59_69) + (sc_uint<1>) ~(comp_r_61_69) + (sc_uint<1>) ~(comp_r_63_69) + (sc_uint<1>) ~(comp_r_65_69) + (sc_uint<1>) ~(comp_r_67_69) + comp_r_69_70 + comp_r_69_71 + comp_r_69_72 + comp_r_69_73 + comp_r_69_74 + comp_r_69_75 + comp_r_69_76 + comp_r_69_77 + comp_r_69_78 + comp_r_69_79 + comp_r_69_80 + comp_r_69_81 + comp_r_69_82 + comp_r_69_83 + comp_r_69_84 + comp_r_69_85 + comp_r_69_86 + comp_r_69_87 + comp_r_69_88 + comp_r_69_89 + comp_r_69_90 + comp_r_69_91 + comp_r_69_92 + comp_r_69_93 + comp_r_69_94 + comp_r_69_95 + comp_r_69_96 + comp_r_69_97 + comp_r_69_98 + comp_r_69_99 + comp_r_69_100 + comp_r_69_101 + comp_r_69_102 + comp_r_69_103 + comp_r_69_104 + comp_r_69_105 + comp_r_69_106 + comp_r_69_107 + comp_r_69_108 + comp_r_69_109 + comp_r_69_110 + comp_r_69_111 + comp_r_69_112 + comp_r_69_113 + comp_r_69_114 + comp_r_69_115 + comp_r_69_116 + comp_r_69_117 + comp_r_69_118 + comp_r_69_119 + comp_r_69_120 + comp_r_69_121 + comp_r_69_122 + comp_r_69_123 + comp_r_69_124 + comp_r_69_125 + comp_r_69_126;
    position[70] = 35 + (sc_uint<1>) ~(comp_r_1_70) + (sc_uint<1>) ~(comp_r_3_70) + (sc_uint<1>) ~(comp_r_5_70) + (sc_uint<1>) ~(comp_r_7_70) + (sc_uint<1>) ~(comp_r_9_70) + (sc_uint<1>) ~(comp_r_11_70) + (sc_uint<1>) ~(comp_r_13_70) + (sc_uint<1>) ~(comp_r_15_70) + (sc_uint<1>) ~(comp_r_17_70) + (sc_uint<1>) ~(comp_r_19_70) + (sc_uint<1>) ~(comp_r_21_70) + (sc_uint<1>) ~(comp_r_23_70) + (sc_uint<1>) ~(comp_r_25_70) + (sc_uint<1>) ~(comp_r_27_70) + (sc_uint<1>) ~(comp_r_29_70) + (sc_uint<1>) ~(comp_r_31_70) + (sc_uint<1>) ~(comp_r_33_70) + (sc_uint<1>) ~(comp_r_35_70) + (sc_uint<1>) ~(comp_r_37_70) + (sc_uint<1>) ~(comp_r_39_70) + (sc_uint<1>) ~(comp_r_41_70) + (sc_uint<1>) ~(comp_r_43_70) + (sc_uint<1>) ~(comp_r_45_70) + (sc_uint<1>) ~(comp_r_47_70) + (sc_uint<1>) ~(comp_r_49_70) + (sc_uint<1>) ~(comp_r_51_70) + (sc_uint<1>) ~(comp_r_53_70) + (sc_uint<1>) ~(comp_r_55_70) + (sc_uint<1>) ~(comp_r_57_70) + (sc_uint<1>) ~(comp_r_59_70) + (sc_uint<1>) ~(comp_r_61_70) + (sc_uint<1>) ~(comp_r_63_70) + (sc_uint<1>) ~(comp_r_65_70) + (sc_uint<1>) ~(comp_r_67_70) + (sc_uint<1>) ~(comp_r_69_70);
    position[71] = 36 + (sc_uint<1>) ~(comp_r_1_71) + (sc_uint<1>) ~(comp_r_3_71) + (sc_uint<1>) ~(comp_r_5_71) + (sc_uint<1>) ~(comp_r_7_71) + (sc_uint<1>) ~(comp_r_9_71) + (sc_uint<1>) ~(comp_r_11_71) + (sc_uint<1>) ~(comp_r_13_71) + (sc_uint<1>) ~(comp_r_15_71) + (sc_uint<1>) ~(comp_r_17_71) + (sc_uint<1>) ~(comp_r_19_71) + (sc_uint<1>) ~(comp_r_21_71) + (sc_uint<1>) ~(comp_r_23_71) + (sc_uint<1>) ~(comp_r_25_71) + (sc_uint<1>) ~(comp_r_27_71) + (sc_uint<1>) ~(comp_r_29_71) + (sc_uint<1>) ~(comp_r_31_71) + (sc_uint<1>) ~(comp_r_33_71) + (sc_uint<1>) ~(comp_r_35_71) + (sc_uint<1>) ~(comp_r_37_71) + (sc_uint<1>) ~(comp_r_39_71) + (sc_uint<1>) ~(comp_r_41_71) + (sc_uint<1>) ~(comp_r_43_71) + (sc_uint<1>) ~(comp_r_45_71) + (sc_uint<1>) ~(comp_r_47_71) + (sc_uint<1>) ~(comp_r_49_71) + (sc_uint<1>) ~(comp_r_51_71) + (sc_uint<1>) ~(comp_r_53_71) + (sc_uint<1>) ~(comp_r_55_71) + (sc_uint<1>) ~(comp_r_57_71) + (sc_uint<1>) ~(comp_r_59_71) + (sc_uint<1>) ~(comp_r_61_71) + (sc_uint<1>) ~(comp_r_63_71) + (sc_uint<1>) ~(comp_r_65_71) + (sc_uint<1>) ~(comp_r_67_71) + (sc_uint<1>) ~(comp_r_69_71) + comp_r_71_72 + comp_r_71_73 + comp_r_71_74 + comp_r_71_75 + comp_r_71_76 + comp_r_71_77 + comp_r_71_78 + comp_r_71_79 + comp_r_71_80 + comp_r_71_81 + comp_r_71_82 + comp_r_71_83 + comp_r_71_84 + comp_r_71_85 + comp_r_71_86 + comp_r_71_87 + comp_r_71_88 + comp_r_71_89 + comp_r_71_90 + comp_r_71_91 + comp_r_71_92 + comp_r_71_93 + comp_r_71_94 + comp_r_71_95 + comp_r_71_96 + comp_r_71_97 + comp_r_71_98 + comp_r_71_99 + comp_r_71_100 + comp_r_71_101 + comp_r_71_102 + comp_r_71_103 + comp_r_71_104 + comp_r_71_105 + comp_r_71_106 + comp_r_71_107 + comp_r_71_108 + comp_r_71_109 + comp_r_71_110 + comp_r_71_111 + comp_r_71_112 + comp_r_71_113 + comp_r_71_114 + comp_r_71_115 + comp_r_71_116 + comp_r_71_117 + comp_r_71_118 + comp_r_71_119 + comp_r_71_120 + comp_r_71_121 + comp_r_71_122 + comp_r_71_123 + comp_r_71_124 + comp_r_71_125 + comp_r_71_126;
    position[72] = 36 + (sc_uint<1>) ~(comp_r_1_72) + (sc_uint<1>) ~(comp_r_3_72) + (sc_uint<1>) ~(comp_r_5_72) + (sc_uint<1>) ~(comp_r_7_72) + (sc_uint<1>) ~(comp_r_9_72) + (sc_uint<1>) ~(comp_r_11_72) + (sc_uint<1>) ~(comp_r_13_72) + (sc_uint<1>) ~(comp_r_15_72) + (sc_uint<1>) ~(comp_r_17_72) + (sc_uint<1>) ~(comp_r_19_72) + (sc_uint<1>) ~(comp_r_21_72) + (sc_uint<1>) ~(comp_r_23_72) + (sc_uint<1>) ~(comp_r_25_72) + (sc_uint<1>) ~(comp_r_27_72) + (sc_uint<1>) ~(comp_r_29_72) + (sc_uint<1>) ~(comp_r_31_72) + (sc_uint<1>) ~(comp_r_33_72) + (sc_uint<1>) ~(comp_r_35_72) + (sc_uint<1>) ~(comp_r_37_72) + (sc_uint<1>) ~(comp_r_39_72) + (sc_uint<1>) ~(comp_r_41_72) + (sc_uint<1>) ~(comp_r_43_72) + (sc_uint<1>) ~(comp_r_45_72) + (sc_uint<1>) ~(comp_r_47_72) + (sc_uint<1>) ~(comp_r_49_72) + (sc_uint<1>) ~(comp_r_51_72) + (sc_uint<1>) ~(comp_r_53_72) + (sc_uint<1>) ~(comp_r_55_72) + (sc_uint<1>) ~(comp_r_57_72) + (sc_uint<1>) ~(comp_r_59_72) + (sc_uint<1>) ~(comp_r_61_72) + (sc_uint<1>) ~(comp_r_63_72) + (sc_uint<1>) ~(comp_r_65_72) + (sc_uint<1>) ~(comp_r_67_72) + (sc_uint<1>) ~(comp_r_69_72) + (sc_uint<1>) ~(comp_r_71_72);
    position[73] = 37 + (sc_uint<1>) ~(comp_r_1_73) + (sc_uint<1>) ~(comp_r_3_73) + (sc_uint<1>) ~(comp_r_5_73) + (sc_uint<1>) ~(comp_r_7_73) + (sc_uint<1>) ~(comp_r_9_73) + (sc_uint<1>) ~(comp_r_11_73) + (sc_uint<1>) ~(comp_r_13_73) + (sc_uint<1>) ~(comp_r_15_73) + (sc_uint<1>) ~(comp_r_17_73) + (sc_uint<1>) ~(comp_r_19_73) + (sc_uint<1>) ~(comp_r_21_73) + (sc_uint<1>) ~(comp_r_23_73) + (sc_uint<1>) ~(comp_r_25_73) + (sc_uint<1>) ~(comp_r_27_73) + (sc_uint<1>) ~(comp_r_29_73) + (sc_uint<1>) ~(comp_r_31_73) + (sc_uint<1>) ~(comp_r_33_73) + (sc_uint<1>) ~(comp_r_35_73) + (sc_uint<1>) ~(comp_r_37_73) + (sc_uint<1>) ~(comp_r_39_73) + (sc_uint<1>) ~(comp_r_41_73) + (sc_uint<1>) ~(comp_r_43_73) + (sc_uint<1>) ~(comp_r_45_73) + (sc_uint<1>) ~(comp_r_47_73) + (sc_uint<1>) ~(comp_r_49_73) + (sc_uint<1>) ~(comp_r_51_73) + (sc_uint<1>) ~(comp_r_53_73) + (sc_uint<1>) ~(comp_r_55_73) + (sc_uint<1>) ~(comp_r_57_73) + (sc_uint<1>) ~(comp_r_59_73) + (sc_uint<1>) ~(comp_r_61_73) + (sc_uint<1>) ~(comp_r_63_73) + (sc_uint<1>) ~(comp_r_65_73) + (sc_uint<1>) ~(comp_r_67_73) + (sc_uint<1>) ~(comp_r_69_73) + (sc_uint<1>) ~(comp_r_71_73) + comp_r_73_74 + comp_r_73_75 + comp_r_73_76 + comp_r_73_77 + comp_r_73_78 + comp_r_73_79 + comp_r_73_80 + comp_r_73_81 + comp_r_73_82 + comp_r_73_83 + comp_r_73_84 + comp_r_73_85 + comp_r_73_86 + comp_r_73_87 + comp_r_73_88 + comp_r_73_89 + comp_r_73_90 + comp_r_73_91 + comp_r_73_92 + comp_r_73_93 + comp_r_73_94 + comp_r_73_95 + comp_r_73_96 + comp_r_73_97 + comp_r_73_98 + comp_r_73_99 + comp_r_73_100 + comp_r_73_101 + comp_r_73_102 + comp_r_73_103 + comp_r_73_104 + comp_r_73_105 + comp_r_73_106 + comp_r_73_107 + comp_r_73_108 + comp_r_73_109 + comp_r_73_110 + comp_r_73_111 + comp_r_73_112 + comp_r_73_113 + comp_r_73_114 + comp_r_73_115 + comp_r_73_116 + comp_r_73_117 + comp_r_73_118 + comp_r_73_119 + comp_r_73_120 + comp_r_73_121 + comp_r_73_122 + comp_r_73_123 + comp_r_73_124 + comp_r_73_125 + comp_r_73_126;
    position[74] = 37 + (sc_uint<1>) ~(comp_r_1_74) + (sc_uint<1>) ~(comp_r_3_74) + (sc_uint<1>) ~(comp_r_5_74) + (sc_uint<1>) ~(comp_r_7_74) + (sc_uint<1>) ~(comp_r_9_74) + (sc_uint<1>) ~(comp_r_11_74) + (sc_uint<1>) ~(comp_r_13_74) + (sc_uint<1>) ~(comp_r_15_74) + (sc_uint<1>) ~(comp_r_17_74) + (sc_uint<1>) ~(comp_r_19_74) + (sc_uint<1>) ~(comp_r_21_74) + (sc_uint<1>) ~(comp_r_23_74) + (sc_uint<1>) ~(comp_r_25_74) + (sc_uint<1>) ~(comp_r_27_74) + (sc_uint<1>) ~(comp_r_29_74) + (sc_uint<1>) ~(comp_r_31_74) + (sc_uint<1>) ~(comp_r_33_74) + (sc_uint<1>) ~(comp_r_35_74) + (sc_uint<1>) ~(comp_r_37_74) + (sc_uint<1>) ~(comp_r_39_74) + (sc_uint<1>) ~(comp_r_41_74) + (sc_uint<1>) ~(comp_r_43_74) + (sc_uint<1>) ~(comp_r_45_74) + (sc_uint<1>) ~(comp_r_47_74) + (sc_uint<1>) ~(comp_r_49_74) + (sc_uint<1>) ~(comp_r_51_74) + (sc_uint<1>) ~(comp_r_53_74) + (sc_uint<1>) ~(comp_r_55_74) + (sc_uint<1>) ~(comp_r_57_74) + (sc_uint<1>) ~(comp_r_59_74) + (sc_uint<1>) ~(comp_r_61_74) + (sc_uint<1>) ~(comp_r_63_74) + (sc_uint<1>) ~(comp_r_65_74) + (sc_uint<1>) ~(comp_r_67_74) + (sc_uint<1>) ~(comp_r_69_74) + (sc_uint<1>) ~(comp_r_71_74) + (sc_uint<1>) ~(comp_r_73_74);
    position[75] = 38 + (sc_uint<1>) ~(comp_r_1_75) + (sc_uint<1>) ~(comp_r_3_75) + (sc_uint<1>) ~(comp_r_5_75) + (sc_uint<1>) ~(comp_r_7_75) + (sc_uint<1>) ~(comp_r_9_75) + (sc_uint<1>) ~(comp_r_11_75) + (sc_uint<1>) ~(comp_r_13_75) + (sc_uint<1>) ~(comp_r_15_75) + (sc_uint<1>) ~(comp_r_17_75) + (sc_uint<1>) ~(comp_r_19_75) + (sc_uint<1>) ~(comp_r_21_75) + (sc_uint<1>) ~(comp_r_23_75) + (sc_uint<1>) ~(comp_r_25_75) + (sc_uint<1>) ~(comp_r_27_75) + (sc_uint<1>) ~(comp_r_29_75) + (sc_uint<1>) ~(comp_r_31_75) + (sc_uint<1>) ~(comp_r_33_75) + (sc_uint<1>) ~(comp_r_35_75) + (sc_uint<1>) ~(comp_r_37_75) + (sc_uint<1>) ~(comp_r_39_75) + (sc_uint<1>) ~(comp_r_41_75) + (sc_uint<1>) ~(comp_r_43_75) + (sc_uint<1>) ~(comp_r_45_75) + (sc_uint<1>) ~(comp_r_47_75) + (sc_uint<1>) ~(comp_r_49_75) + (sc_uint<1>) ~(comp_r_51_75) + (sc_uint<1>) ~(comp_r_53_75) + (sc_uint<1>) ~(comp_r_55_75) + (sc_uint<1>) ~(comp_r_57_75) + (sc_uint<1>) ~(comp_r_59_75) + (sc_uint<1>) ~(comp_r_61_75) + (sc_uint<1>) ~(comp_r_63_75) + (sc_uint<1>) ~(comp_r_65_75) + (sc_uint<1>) ~(comp_r_67_75) + (sc_uint<1>) ~(comp_r_69_75) + (sc_uint<1>) ~(comp_r_71_75) + (sc_uint<1>) ~(comp_r_73_75) + comp_r_75_76 + comp_r_75_77 + comp_r_75_78 + comp_r_75_79 + comp_r_75_80 + comp_r_75_81 + comp_r_75_82 + comp_r_75_83 + comp_r_75_84 + comp_r_75_85 + comp_r_75_86 + comp_r_75_87 + comp_r_75_88 + comp_r_75_89 + comp_r_75_90 + comp_r_75_91 + comp_r_75_92 + comp_r_75_93 + comp_r_75_94 + comp_r_75_95 + comp_r_75_96 + comp_r_75_97 + comp_r_75_98 + comp_r_75_99 + comp_r_75_100 + comp_r_75_101 + comp_r_75_102 + comp_r_75_103 + comp_r_75_104 + comp_r_75_105 + comp_r_75_106 + comp_r_75_107 + comp_r_75_108 + comp_r_75_109 + comp_r_75_110 + comp_r_75_111 + comp_r_75_112 + comp_r_75_113 + comp_r_75_114 + comp_r_75_115 + comp_r_75_116 + comp_r_75_117 + comp_r_75_118 + comp_r_75_119 + comp_r_75_120 + comp_r_75_121 + comp_r_75_122 + comp_r_75_123 + comp_r_75_124 + comp_r_75_125 + comp_r_75_126;
    position[76] = 38 + (sc_uint<1>) ~(comp_r_1_76) + (sc_uint<1>) ~(comp_r_3_76) + (sc_uint<1>) ~(comp_r_5_76) + (sc_uint<1>) ~(comp_r_7_76) + (sc_uint<1>) ~(comp_r_9_76) + (sc_uint<1>) ~(comp_r_11_76) + (sc_uint<1>) ~(comp_r_13_76) + (sc_uint<1>) ~(comp_r_15_76) + (sc_uint<1>) ~(comp_r_17_76) + (sc_uint<1>) ~(comp_r_19_76) + (sc_uint<1>) ~(comp_r_21_76) + (sc_uint<1>) ~(comp_r_23_76) + (sc_uint<1>) ~(comp_r_25_76) + (sc_uint<1>) ~(comp_r_27_76) + (sc_uint<1>) ~(comp_r_29_76) + (sc_uint<1>) ~(comp_r_31_76) + (sc_uint<1>) ~(comp_r_33_76) + (sc_uint<1>) ~(comp_r_35_76) + (sc_uint<1>) ~(comp_r_37_76) + (sc_uint<1>) ~(comp_r_39_76) + (sc_uint<1>) ~(comp_r_41_76) + (sc_uint<1>) ~(comp_r_43_76) + (sc_uint<1>) ~(comp_r_45_76) + (sc_uint<1>) ~(comp_r_47_76) + (sc_uint<1>) ~(comp_r_49_76) + (sc_uint<1>) ~(comp_r_51_76) + (sc_uint<1>) ~(comp_r_53_76) + (sc_uint<1>) ~(comp_r_55_76) + (sc_uint<1>) ~(comp_r_57_76) + (sc_uint<1>) ~(comp_r_59_76) + (sc_uint<1>) ~(comp_r_61_76) + (sc_uint<1>) ~(comp_r_63_76) + (sc_uint<1>) ~(comp_r_65_76) + (sc_uint<1>) ~(comp_r_67_76) + (sc_uint<1>) ~(comp_r_69_76) + (sc_uint<1>) ~(comp_r_71_76) + (sc_uint<1>) ~(comp_r_73_76) + (sc_uint<1>) ~(comp_r_75_76);
    position[77] = 39 + (sc_uint<1>) ~(comp_r_1_77) + (sc_uint<1>) ~(comp_r_3_77) + (sc_uint<1>) ~(comp_r_5_77) + (sc_uint<1>) ~(comp_r_7_77) + (sc_uint<1>) ~(comp_r_9_77) + (sc_uint<1>) ~(comp_r_11_77) + (sc_uint<1>) ~(comp_r_13_77) + (sc_uint<1>) ~(comp_r_15_77) + (sc_uint<1>) ~(comp_r_17_77) + (sc_uint<1>) ~(comp_r_19_77) + (sc_uint<1>) ~(comp_r_21_77) + (sc_uint<1>) ~(comp_r_23_77) + (sc_uint<1>) ~(comp_r_25_77) + (sc_uint<1>) ~(comp_r_27_77) + (sc_uint<1>) ~(comp_r_29_77) + (sc_uint<1>) ~(comp_r_31_77) + (sc_uint<1>) ~(comp_r_33_77) + (sc_uint<1>) ~(comp_r_35_77) + (sc_uint<1>) ~(comp_r_37_77) + (sc_uint<1>) ~(comp_r_39_77) + (sc_uint<1>) ~(comp_r_41_77) + (sc_uint<1>) ~(comp_r_43_77) + (sc_uint<1>) ~(comp_r_45_77) + (sc_uint<1>) ~(comp_r_47_77) + (sc_uint<1>) ~(comp_r_49_77) + (sc_uint<1>) ~(comp_r_51_77) + (sc_uint<1>) ~(comp_r_53_77) + (sc_uint<1>) ~(comp_r_55_77) + (sc_uint<1>) ~(comp_r_57_77) + (sc_uint<1>) ~(comp_r_59_77) + (sc_uint<1>) ~(comp_r_61_77) + (sc_uint<1>) ~(comp_r_63_77) + (sc_uint<1>) ~(comp_r_65_77) + (sc_uint<1>) ~(comp_r_67_77) + (sc_uint<1>) ~(comp_r_69_77) + (sc_uint<1>) ~(comp_r_71_77) + (sc_uint<1>) ~(comp_r_73_77) + (sc_uint<1>) ~(comp_r_75_77) + comp_r_77_78 + comp_r_77_79 + comp_r_77_80 + comp_r_77_81 + comp_r_77_82 + comp_r_77_83 + comp_r_77_84 + comp_r_77_85 + comp_r_77_86 + comp_r_77_87 + comp_r_77_88 + comp_r_77_89 + comp_r_77_90 + comp_r_77_91 + comp_r_77_92 + comp_r_77_93 + comp_r_77_94 + comp_r_77_95 + comp_r_77_96 + comp_r_77_97 + comp_r_77_98 + comp_r_77_99 + comp_r_77_100 + comp_r_77_101 + comp_r_77_102 + comp_r_77_103 + comp_r_77_104 + comp_r_77_105 + comp_r_77_106 + comp_r_77_107 + comp_r_77_108 + comp_r_77_109 + comp_r_77_110 + comp_r_77_111 + comp_r_77_112 + comp_r_77_113 + comp_r_77_114 + comp_r_77_115 + comp_r_77_116 + comp_r_77_117 + comp_r_77_118 + comp_r_77_119 + comp_r_77_120 + comp_r_77_121 + comp_r_77_122 + comp_r_77_123 + comp_r_77_124 + comp_r_77_125 + comp_r_77_126;
    position[78] = 39 + (sc_uint<1>) ~(comp_r_1_78) + (sc_uint<1>) ~(comp_r_3_78) + (sc_uint<1>) ~(comp_r_5_78) + (sc_uint<1>) ~(comp_r_7_78) + (sc_uint<1>) ~(comp_r_9_78) + (sc_uint<1>) ~(comp_r_11_78) + (sc_uint<1>) ~(comp_r_13_78) + (sc_uint<1>) ~(comp_r_15_78) + (sc_uint<1>) ~(comp_r_17_78) + (sc_uint<1>) ~(comp_r_19_78) + (sc_uint<1>) ~(comp_r_21_78) + (sc_uint<1>) ~(comp_r_23_78) + (sc_uint<1>) ~(comp_r_25_78) + (sc_uint<1>) ~(comp_r_27_78) + (sc_uint<1>) ~(comp_r_29_78) + (sc_uint<1>) ~(comp_r_31_78) + (sc_uint<1>) ~(comp_r_33_78) + (sc_uint<1>) ~(comp_r_35_78) + (sc_uint<1>) ~(comp_r_37_78) + (sc_uint<1>) ~(comp_r_39_78) + (sc_uint<1>) ~(comp_r_41_78) + (sc_uint<1>) ~(comp_r_43_78) + (sc_uint<1>) ~(comp_r_45_78) + (sc_uint<1>) ~(comp_r_47_78) + (sc_uint<1>) ~(comp_r_49_78) + (sc_uint<1>) ~(comp_r_51_78) + (sc_uint<1>) ~(comp_r_53_78) + (sc_uint<1>) ~(comp_r_55_78) + (sc_uint<1>) ~(comp_r_57_78) + (sc_uint<1>) ~(comp_r_59_78) + (sc_uint<1>) ~(comp_r_61_78) + (sc_uint<1>) ~(comp_r_63_78) + (sc_uint<1>) ~(comp_r_65_78) + (sc_uint<1>) ~(comp_r_67_78) + (sc_uint<1>) ~(comp_r_69_78) + (sc_uint<1>) ~(comp_r_71_78) + (sc_uint<1>) ~(comp_r_73_78) + (sc_uint<1>) ~(comp_r_75_78) + (sc_uint<1>) ~(comp_r_77_78);
    position[79] = 40 + (sc_uint<1>) ~(comp_r_1_79) + (sc_uint<1>) ~(comp_r_3_79) + (sc_uint<1>) ~(comp_r_5_79) + (sc_uint<1>) ~(comp_r_7_79) + (sc_uint<1>) ~(comp_r_9_79) + (sc_uint<1>) ~(comp_r_11_79) + (sc_uint<1>) ~(comp_r_13_79) + (sc_uint<1>) ~(comp_r_15_79) + (sc_uint<1>) ~(comp_r_17_79) + (sc_uint<1>) ~(comp_r_19_79) + (sc_uint<1>) ~(comp_r_21_79) + (sc_uint<1>) ~(comp_r_23_79) + (sc_uint<1>) ~(comp_r_25_79) + (sc_uint<1>) ~(comp_r_27_79) + (sc_uint<1>) ~(comp_r_29_79) + (sc_uint<1>) ~(comp_r_31_79) + (sc_uint<1>) ~(comp_r_33_79) + (sc_uint<1>) ~(comp_r_35_79) + (sc_uint<1>) ~(comp_r_37_79) + (sc_uint<1>) ~(comp_r_39_79) + (sc_uint<1>) ~(comp_r_41_79) + (sc_uint<1>) ~(comp_r_43_79) + (sc_uint<1>) ~(comp_r_45_79) + (sc_uint<1>) ~(comp_r_47_79) + (sc_uint<1>) ~(comp_r_49_79) + (sc_uint<1>) ~(comp_r_51_79) + (sc_uint<1>) ~(comp_r_53_79) + (sc_uint<1>) ~(comp_r_55_79) + (sc_uint<1>) ~(comp_r_57_79) + (sc_uint<1>) ~(comp_r_59_79) + (sc_uint<1>) ~(comp_r_61_79) + (sc_uint<1>) ~(comp_r_63_79) + (sc_uint<1>) ~(comp_r_65_79) + (sc_uint<1>) ~(comp_r_67_79) + (sc_uint<1>) ~(comp_r_69_79) + (sc_uint<1>) ~(comp_r_71_79) + (sc_uint<1>) ~(comp_r_73_79) + (sc_uint<1>) ~(comp_r_75_79) + (sc_uint<1>) ~(comp_r_77_79) + comp_r_79_80 + comp_r_79_81 + comp_r_79_82 + comp_r_79_83 + comp_r_79_84 + comp_r_79_85 + comp_r_79_86 + comp_r_79_87 + comp_r_79_88 + comp_r_79_89 + comp_r_79_90 + comp_r_79_91 + comp_r_79_92 + comp_r_79_93 + comp_r_79_94 + comp_r_79_95 + comp_r_79_96 + comp_r_79_97 + comp_r_79_98 + comp_r_79_99 + comp_r_79_100 + comp_r_79_101 + comp_r_79_102 + comp_r_79_103 + comp_r_79_104 + comp_r_79_105 + comp_r_79_106 + comp_r_79_107 + comp_r_79_108 + comp_r_79_109 + comp_r_79_110 + comp_r_79_111 + comp_r_79_112 + comp_r_79_113 + comp_r_79_114 + comp_r_79_115 + comp_r_79_116 + comp_r_79_117 + comp_r_79_118 + comp_r_79_119 + comp_r_79_120 + comp_r_79_121 + comp_r_79_122 + comp_r_79_123 + comp_r_79_124 + comp_r_79_125 + comp_r_79_126;
    position[80] = 40 + (sc_uint<1>) ~(comp_r_1_80) + (sc_uint<1>) ~(comp_r_3_80) + (sc_uint<1>) ~(comp_r_5_80) + (sc_uint<1>) ~(comp_r_7_80) + (sc_uint<1>) ~(comp_r_9_80) + (sc_uint<1>) ~(comp_r_11_80) + (sc_uint<1>) ~(comp_r_13_80) + (sc_uint<1>) ~(comp_r_15_80) + (sc_uint<1>) ~(comp_r_17_80) + (sc_uint<1>) ~(comp_r_19_80) + (sc_uint<1>) ~(comp_r_21_80) + (sc_uint<1>) ~(comp_r_23_80) + (sc_uint<1>) ~(comp_r_25_80) + (sc_uint<1>) ~(comp_r_27_80) + (sc_uint<1>) ~(comp_r_29_80) + (sc_uint<1>) ~(comp_r_31_80) + (sc_uint<1>) ~(comp_r_33_80) + (sc_uint<1>) ~(comp_r_35_80) + (sc_uint<1>) ~(comp_r_37_80) + (sc_uint<1>) ~(comp_r_39_80) + (sc_uint<1>) ~(comp_r_41_80) + (sc_uint<1>) ~(comp_r_43_80) + (sc_uint<1>) ~(comp_r_45_80) + (sc_uint<1>) ~(comp_r_47_80) + (sc_uint<1>) ~(comp_r_49_80) + (sc_uint<1>) ~(comp_r_51_80) + (sc_uint<1>) ~(comp_r_53_80) + (sc_uint<1>) ~(comp_r_55_80) + (sc_uint<1>) ~(comp_r_57_80) + (sc_uint<1>) ~(comp_r_59_80) + (sc_uint<1>) ~(comp_r_61_80) + (sc_uint<1>) ~(comp_r_63_80) + (sc_uint<1>) ~(comp_r_65_80) + (sc_uint<1>) ~(comp_r_67_80) + (sc_uint<1>) ~(comp_r_69_80) + (sc_uint<1>) ~(comp_r_71_80) + (sc_uint<1>) ~(comp_r_73_80) + (sc_uint<1>) ~(comp_r_75_80) + (sc_uint<1>) ~(comp_r_77_80) + (sc_uint<1>) ~(comp_r_79_80);
    position[81] = 41 + (sc_uint<1>) ~(comp_r_1_81) + (sc_uint<1>) ~(comp_r_3_81) + (sc_uint<1>) ~(comp_r_5_81) + (sc_uint<1>) ~(comp_r_7_81) + (sc_uint<1>) ~(comp_r_9_81) + (sc_uint<1>) ~(comp_r_11_81) + (sc_uint<1>) ~(comp_r_13_81) + (sc_uint<1>) ~(comp_r_15_81) + (sc_uint<1>) ~(comp_r_17_81) + (sc_uint<1>) ~(comp_r_19_81) + (sc_uint<1>) ~(comp_r_21_81) + (sc_uint<1>) ~(comp_r_23_81) + (sc_uint<1>) ~(comp_r_25_81) + (sc_uint<1>) ~(comp_r_27_81) + (sc_uint<1>) ~(comp_r_29_81) + (sc_uint<1>) ~(comp_r_31_81) + (sc_uint<1>) ~(comp_r_33_81) + (sc_uint<1>) ~(comp_r_35_81) + (sc_uint<1>) ~(comp_r_37_81) + (sc_uint<1>) ~(comp_r_39_81) + (sc_uint<1>) ~(comp_r_41_81) + (sc_uint<1>) ~(comp_r_43_81) + (sc_uint<1>) ~(comp_r_45_81) + (sc_uint<1>) ~(comp_r_47_81) + (sc_uint<1>) ~(comp_r_49_81) + (sc_uint<1>) ~(comp_r_51_81) + (sc_uint<1>) ~(comp_r_53_81) + (sc_uint<1>) ~(comp_r_55_81) + (sc_uint<1>) ~(comp_r_57_81) + (sc_uint<1>) ~(comp_r_59_81) + (sc_uint<1>) ~(comp_r_61_81) + (sc_uint<1>) ~(comp_r_63_81) + (sc_uint<1>) ~(comp_r_65_81) + (sc_uint<1>) ~(comp_r_67_81) + (sc_uint<1>) ~(comp_r_69_81) + (sc_uint<1>) ~(comp_r_71_81) + (sc_uint<1>) ~(comp_r_73_81) + (sc_uint<1>) ~(comp_r_75_81) + (sc_uint<1>) ~(comp_r_77_81) + (sc_uint<1>) ~(comp_r_79_81) + comp_r_81_82 + comp_r_81_83 + comp_r_81_84 + comp_r_81_85 + comp_r_81_86 + comp_r_81_87 + comp_r_81_88 + comp_r_81_89 + comp_r_81_90 + comp_r_81_91 + comp_r_81_92 + comp_r_81_93 + comp_r_81_94 + comp_r_81_95 + comp_r_81_96 + comp_r_81_97 + comp_r_81_98 + comp_r_81_99 + comp_r_81_100 + comp_r_81_101 + comp_r_81_102 + comp_r_81_103 + comp_r_81_104 + comp_r_81_105 + comp_r_81_106 + comp_r_81_107 + comp_r_81_108 + comp_r_81_109 + comp_r_81_110 + comp_r_81_111 + comp_r_81_112 + comp_r_81_113 + comp_r_81_114 + comp_r_81_115 + comp_r_81_116 + comp_r_81_117 + comp_r_81_118 + comp_r_81_119 + comp_r_81_120 + comp_r_81_121 + comp_r_81_122 + comp_r_81_123 + comp_r_81_124 + comp_r_81_125 + comp_r_81_126;
    position[82] = 41 + (sc_uint<1>) ~(comp_r_1_82) + (sc_uint<1>) ~(comp_r_3_82) + (sc_uint<1>) ~(comp_r_5_82) + (sc_uint<1>) ~(comp_r_7_82) + (sc_uint<1>) ~(comp_r_9_82) + (sc_uint<1>) ~(comp_r_11_82) + (sc_uint<1>) ~(comp_r_13_82) + (sc_uint<1>) ~(comp_r_15_82) + (sc_uint<1>) ~(comp_r_17_82) + (sc_uint<1>) ~(comp_r_19_82) + (sc_uint<1>) ~(comp_r_21_82) + (sc_uint<1>) ~(comp_r_23_82) + (sc_uint<1>) ~(comp_r_25_82) + (sc_uint<1>) ~(comp_r_27_82) + (sc_uint<1>) ~(comp_r_29_82) + (sc_uint<1>) ~(comp_r_31_82) + (sc_uint<1>) ~(comp_r_33_82) + (sc_uint<1>) ~(comp_r_35_82) + (sc_uint<1>) ~(comp_r_37_82) + (sc_uint<1>) ~(comp_r_39_82) + (sc_uint<1>) ~(comp_r_41_82) + (sc_uint<1>) ~(comp_r_43_82) + (sc_uint<1>) ~(comp_r_45_82) + (sc_uint<1>) ~(comp_r_47_82) + (sc_uint<1>) ~(comp_r_49_82) + (sc_uint<1>) ~(comp_r_51_82) + (sc_uint<1>) ~(comp_r_53_82) + (sc_uint<1>) ~(comp_r_55_82) + (sc_uint<1>) ~(comp_r_57_82) + (sc_uint<1>) ~(comp_r_59_82) + (sc_uint<1>) ~(comp_r_61_82) + (sc_uint<1>) ~(comp_r_63_82) + (sc_uint<1>) ~(comp_r_65_82) + (sc_uint<1>) ~(comp_r_67_82) + (sc_uint<1>) ~(comp_r_69_82) + (sc_uint<1>) ~(comp_r_71_82) + (sc_uint<1>) ~(comp_r_73_82) + (sc_uint<1>) ~(comp_r_75_82) + (sc_uint<1>) ~(comp_r_77_82) + (sc_uint<1>) ~(comp_r_79_82) + (sc_uint<1>) ~(comp_r_81_82);
    position[83] = 42 + (sc_uint<1>) ~(comp_r_1_83) + (sc_uint<1>) ~(comp_r_3_83) + (sc_uint<1>) ~(comp_r_5_83) + (sc_uint<1>) ~(comp_r_7_83) + (sc_uint<1>) ~(comp_r_9_83) + (sc_uint<1>) ~(comp_r_11_83) + (sc_uint<1>) ~(comp_r_13_83) + (sc_uint<1>) ~(comp_r_15_83) + (sc_uint<1>) ~(comp_r_17_83) + (sc_uint<1>) ~(comp_r_19_83) + (sc_uint<1>) ~(comp_r_21_83) + (sc_uint<1>) ~(comp_r_23_83) + (sc_uint<1>) ~(comp_r_25_83) + (sc_uint<1>) ~(comp_r_27_83) + (sc_uint<1>) ~(comp_r_29_83) + (sc_uint<1>) ~(comp_r_31_83) + (sc_uint<1>) ~(comp_r_33_83) + (sc_uint<1>) ~(comp_r_35_83) + (sc_uint<1>) ~(comp_r_37_83) + (sc_uint<1>) ~(comp_r_39_83) + (sc_uint<1>) ~(comp_r_41_83) + (sc_uint<1>) ~(comp_r_43_83) + (sc_uint<1>) ~(comp_r_45_83) + (sc_uint<1>) ~(comp_r_47_83) + (sc_uint<1>) ~(comp_r_49_83) + (sc_uint<1>) ~(comp_r_51_83) + (sc_uint<1>) ~(comp_r_53_83) + (sc_uint<1>) ~(comp_r_55_83) + (sc_uint<1>) ~(comp_r_57_83) + (sc_uint<1>) ~(comp_r_59_83) + (sc_uint<1>) ~(comp_r_61_83) + (sc_uint<1>) ~(comp_r_63_83) + (sc_uint<1>) ~(comp_r_65_83) + (sc_uint<1>) ~(comp_r_67_83) + (sc_uint<1>) ~(comp_r_69_83) + (sc_uint<1>) ~(comp_r_71_83) + (sc_uint<1>) ~(comp_r_73_83) + (sc_uint<1>) ~(comp_r_75_83) + (sc_uint<1>) ~(comp_r_77_83) + (sc_uint<1>) ~(comp_r_79_83) + (sc_uint<1>) ~(comp_r_81_83) + comp_r_83_84 + comp_r_83_85 + comp_r_83_86 + comp_r_83_87 + comp_r_83_88 + comp_r_83_89 + comp_r_83_90 + comp_r_83_91 + comp_r_83_92 + comp_r_83_93 + comp_r_83_94 + comp_r_83_95 + comp_r_83_96 + comp_r_83_97 + comp_r_83_98 + comp_r_83_99 + comp_r_83_100 + comp_r_83_101 + comp_r_83_102 + comp_r_83_103 + comp_r_83_104 + comp_r_83_105 + comp_r_83_106 + comp_r_83_107 + comp_r_83_108 + comp_r_83_109 + comp_r_83_110 + comp_r_83_111 + comp_r_83_112 + comp_r_83_113 + comp_r_83_114 + comp_r_83_115 + comp_r_83_116 + comp_r_83_117 + comp_r_83_118 + comp_r_83_119 + comp_r_83_120 + comp_r_83_121 + comp_r_83_122 + comp_r_83_123 + comp_r_83_124 + comp_r_83_125 + comp_r_83_126;
    position[84] = 42 + (sc_uint<1>) ~(comp_r_1_84) + (sc_uint<1>) ~(comp_r_3_84) + (sc_uint<1>) ~(comp_r_5_84) + (sc_uint<1>) ~(comp_r_7_84) + (sc_uint<1>) ~(comp_r_9_84) + (sc_uint<1>) ~(comp_r_11_84) + (sc_uint<1>) ~(comp_r_13_84) + (sc_uint<1>) ~(comp_r_15_84) + (sc_uint<1>) ~(comp_r_17_84) + (sc_uint<1>) ~(comp_r_19_84) + (sc_uint<1>) ~(comp_r_21_84) + (sc_uint<1>) ~(comp_r_23_84) + (sc_uint<1>) ~(comp_r_25_84) + (sc_uint<1>) ~(comp_r_27_84) + (sc_uint<1>) ~(comp_r_29_84) + (sc_uint<1>) ~(comp_r_31_84) + (sc_uint<1>) ~(comp_r_33_84) + (sc_uint<1>) ~(comp_r_35_84) + (sc_uint<1>) ~(comp_r_37_84) + (sc_uint<1>) ~(comp_r_39_84) + (sc_uint<1>) ~(comp_r_41_84) + (sc_uint<1>) ~(comp_r_43_84) + (sc_uint<1>) ~(comp_r_45_84) + (sc_uint<1>) ~(comp_r_47_84) + (sc_uint<1>) ~(comp_r_49_84) + (sc_uint<1>) ~(comp_r_51_84) + (sc_uint<1>) ~(comp_r_53_84) + (sc_uint<1>) ~(comp_r_55_84) + (sc_uint<1>) ~(comp_r_57_84) + (sc_uint<1>) ~(comp_r_59_84) + (sc_uint<1>) ~(comp_r_61_84) + (sc_uint<1>) ~(comp_r_63_84) + (sc_uint<1>) ~(comp_r_65_84) + (sc_uint<1>) ~(comp_r_67_84) + (sc_uint<1>) ~(comp_r_69_84) + (sc_uint<1>) ~(comp_r_71_84) + (sc_uint<1>) ~(comp_r_73_84) + (sc_uint<1>) ~(comp_r_75_84) + (sc_uint<1>) ~(comp_r_77_84) + (sc_uint<1>) ~(comp_r_79_84) + (sc_uint<1>) ~(comp_r_81_84) + (sc_uint<1>) ~(comp_r_83_84);
    position[85] = 43 + (sc_uint<1>) ~(comp_r_1_85) + (sc_uint<1>) ~(comp_r_3_85) + (sc_uint<1>) ~(comp_r_5_85) + (sc_uint<1>) ~(comp_r_7_85) + (sc_uint<1>) ~(comp_r_9_85) + (sc_uint<1>) ~(comp_r_11_85) + (sc_uint<1>) ~(comp_r_13_85) + (sc_uint<1>) ~(comp_r_15_85) + (sc_uint<1>) ~(comp_r_17_85) + (sc_uint<1>) ~(comp_r_19_85) + (sc_uint<1>) ~(comp_r_21_85) + (sc_uint<1>) ~(comp_r_23_85) + (sc_uint<1>) ~(comp_r_25_85) + (sc_uint<1>) ~(comp_r_27_85) + (sc_uint<1>) ~(comp_r_29_85) + (sc_uint<1>) ~(comp_r_31_85) + (sc_uint<1>) ~(comp_r_33_85) + (sc_uint<1>) ~(comp_r_35_85) + (sc_uint<1>) ~(comp_r_37_85) + (sc_uint<1>) ~(comp_r_39_85) + (sc_uint<1>) ~(comp_r_41_85) + (sc_uint<1>) ~(comp_r_43_85) + (sc_uint<1>) ~(comp_r_45_85) + (sc_uint<1>) ~(comp_r_47_85) + (sc_uint<1>) ~(comp_r_49_85) + (sc_uint<1>) ~(comp_r_51_85) + (sc_uint<1>) ~(comp_r_53_85) + (sc_uint<1>) ~(comp_r_55_85) + (sc_uint<1>) ~(comp_r_57_85) + (sc_uint<1>) ~(comp_r_59_85) + (sc_uint<1>) ~(comp_r_61_85) + (sc_uint<1>) ~(comp_r_63_85) + (sc_uint<1>) ~(comp_r_65_85) + (sc_uint<1>) ~(comp_r_67_85) + (sc_uint<1>) ~(comp_r_69_85) + (sc_uint<1>) ~(comp_r_71_85) + (sc_uint<1>) ~(comp_r_73_85) + (sc_uint<1>) ~(comp_r_75_85) + (sc_uint<1>) ~(comp_r_77_85) + (sc_uint<1>) ~(comp_r_79_85) + (sc_uint<1>) ~(comp_r_81_85) + (sc_uint<1>) ~(comp_r_83_85) + comp_r_85_86 + comp_r_85_87 + comp_r_85_88 + comp_r_85_89 + comp_r_85_90 + comp_r_85_91 + comp_r_85_92 + comp_r_85_93 + comp_r_85_94 + comp_r_85_95 + comp_r_85_96 + comp_r_85_97 + comp_r_85_98 + comp_r_85_99 + comp_r_85_100 + comp_r_85_101 + comp_r_85_102 + comp_r_85_103 + comp_r_85_104 + comp_r_85_105 + comp_r_85_106 + comp_r_85_107 + comp_r_85_108 + comp_r_85_109 + comp_r_85_110 + comp_r_85_111 + comp_r_85_112 + comp_r_85_113 + comp_r_85_114 + comp_r_85_115 + comp_r_85_116 + comp_r_85_117 + comp_r_85_118 + comp_r_85_119 + comp_r_85_120 + comp_r_85_121 + comp_r_85_122 + comp_r_85_123 + comp_r_85_124 + comp_r_85_125 + comp_r_85_126;
    position[86] = 43 + (sc_uint<1>) ~(comp_r_1_86) + (sc_uint<1>) ~(comp_r_3_86) + (sc_uint<1>) ~(comp_r_5_86) + (sc_uint<1>) ~(comp_r_7_86) + (sc_uint<1>) ~(comp_r_9_86) + (sc_uint<1>) ~(comp_r_11_86) + (sc_uint<1>) ~(comp_r_13_86) + (sc_uint<1>) ~(comp_r_15_86) + (sc_uint<1>) ~(comp_r_17_86) + (sc_uint<1>) ~(comp_r_19_86) + (sc_uint<1>) ~(comp_r_21_86) + (sc_uint<1>) ~(comp_r_23_86) + (sc_uint<1>) ~(comp_r_25_86) + (sc_uint<1>) ~(comp_r_27_86) + (sc_uint<1>) ~(comp_r_29_86) + (sc_uint<1>) ~(comp_r_31_86) + (sc_uint<1>) ~(comp_r_33_86) + (sc_uint<1>) ~(comp_r_35_86) + (sc_uint<1>) ~(comp_r_37_86) + (sc_uint<1>) ~(comp_r_39_86) + (sc_uint<1>) ~(comp_r_41_86) + (sc_uint<1>) ~(comp_r_43_86) + (sc_uint<1>) ~(comp_r_45_86) + (sc_uint<1>) ~(comp_r_47_86) + (sc_uint<1>) ~(comp_r_49_86) + (sc_uint<1>) ~(comp_r_51_86) + (sc_uint<1>) ~(comp_r_53_86) + (sc_uint<1>) ~(comp_r_55_86) + (sc_uint<1>) ~(comp_r_57_86) + (sc_uint<1>) ~(comp_r_59_86) + (sc_uint<1>) ~(comp_r_61_86) + (sc_uint<1>) ~(comp_r_63_86) + (sc_uint<1>) ~(comp_r_65_86) + (sc_uint<1>) ~(comp_r_67_86) + (sc_uint<1>) ~(comp_r_69_86) + (sc_uint<1>) ~(comp_r_71_86) + (sc_uint<1>) ~(comp_r_73_86) + (sc_uint<1>) ~(comp_r_75_86) + (sc_uint<1>) ~(comp_r_77_86) + (sc_uint<1>) ~(comp_r_79_86) + (sc_uint<1>) ~(comp_r_81_86) + (sc_uint<1>) ~(comp_r_83_86) + (sc_uint<1>) ~(comp_r_85_86);
    position[87] = 44 + (sc_uint<1>) ~(comp_r_1_87) + (sc_uint<1>) ~(comp_r_3_87) + (sc_uint<1>) ~(comp_r_5_87) + (sc_uint<1>) ~(comp_r_7_87) + (sc_uint<1>) ~(comp_r_9_87) + (sc_uint<1>) ~(comp_r_11_87) + (sc_uint<1>) ~(comp_r_13_87) + (sc_uint<1>) ~(comp_r_15_87) + (sc_uint<1>) ~(comp_r_17_87) + (sc_uint<1>) ~(comp_r_19_87) + (sc_uint<1>) ~(comp_r_21_87) + (sc_uint<1>) ~(comp_r_23_87) + (sc_uint<1>) ~(comp_r_25_87) + (sc_uint<1>) ~(comp_r_27_87) + (sc_uint<1>) ~(comp_r_29_87) + (sc_uint<1>) ~(comp_r_31_87) + (sc_uint<1>) ~(comp_r_33_87) + (sc_uint<1>) ~(comp_r_35_87) + (sc_uint<1>) ~(comp_r_37_87) + (sc_uint<1>) ~(comp_r_39_87) + (sc_uint<1>) ~(comp_r_41_87) + (sc_uint<1>) ~(comp_r_43_87) + (sc_uint<1>) ~(comp_r_45_87) + (sc_uint<1>) ~(comp_r_47_87) + (sc_uint<1>) ~(comp_r_49_87) + (sc_uint<1>) ~(comp_r_51_87) + (sc_uint<1>) ~(comp_r_53_87) + (sc_uint<1>) ~(comp_r_55_87) + (sc_uint<1>) ~(comp_r_57_87) + (sc_uint<1>) ~(comp_r_59_87) + (sc_uint<1>) ~(comp_r_61_87) + (sc_uint<1>) ~(comp_r_63_87) + (sc_uint<1>) ~(comp_r_65_87) + (sc_uint<1>) ~(comp_r_67_87) + (sc_uint<1>) ~(comp_r_69_87) + (sc_uint<1>) ~(comp_r_71_87) + (sc_uint<1>) ~(comp_r_73_87) + (sc_uint<1>) ~(comp_r_75_87) + (sc_uint<1>) ~(comp_r_77_87) + (sc_uint<1>) ~(comp_r_79_87) + (sc_uint<1>) ~(comp_r_81_87) + (sc_uint<1>) ~(comp_r_83_87) + (sc_uint<1>) ~(comp_r_85_87) + comp_r_87_88 + comp_r_87_89 + comp_r_87_90 + comp_r_87_91 + comp_r_87_92 + comp_r_87_93 + comp_r_87_94 + comp_r_87_95 + comp_r_87_96 + comp_r_87_97 + comp_r_87_98 + comp_r_87_99 + comp_r_87_100 + comp_r_87_101 + comp_r_87_102 + comp_r_87_103 + comp_r_87_104 + comp_r_87_105 + comp_r_87_106 + comp_r_87_107 + comp_r_87_108 + comp_r_87_109 + comp_r_87_110 + comp_r_87_111 + comp_r_87_112 + comp_r_87_113 + comp_r_87_114 + comp_r_87_115 + comp_r_87_116 + comp_r_87_117 + comp_r_87_118 + comp_r_87_119 + comp_r_87_120 + comp_r_87_121 + comp_r_87_122 + comp_r_87_123 + comp_r_87_124 + comp_r_87_125 + comp_r_87_126;
    position[88] = 44 + (sc_uint<1>) ~(comp_r_1_88) + (sc_uint<1>) ~(comp_r_3_88) + (sc_uint<1>) ~(comp_r_5_88) + (sc_uint<1>) ~(comp_r_7_88) + (sc_uint<1>) ~(comp_r_9_88) + (sc_uint<1>) ~(comp_r_11_88) + (sc_uint<1>) ~(comp_r_13_88) + (sc_uint<1>) ~(comp_r_15_88) + (sc_uint<1>) ~(comp_r_17_88) + (sc_uint<1>) ~(comp_r_19_88) + (sc_uint<1>) ~(comp_r_21_88) + (sc_uint<1>) ~(comp_r_23_88) + (sc_uint<1>) ~(comp_r_25_88) + (sc_uint<1>) ~(comp_r_27_88) + (sc_uint<1>) ~(comp_r_29_88) + (sc_uint<1>) ~(comp_r_31_88) + (sc_uint<1>) ~(comp_r_33_88) + (sc_uint<1>) ~(comp_r_35_88) + (sc_uint<1>) ~(comp_r_37_88) + (sc_uint<1>) ~(comp_r_39_88) + (sc_uint<1>) ~(comp_r_41_88) + (sc_uint<1>) ~(comp_r_43_88) + (sc_uint<1>) ~(comp_r_45_88) + (sc_uint<1>) ~(comp_r_47_88) + (sc_uint<1>) ~(comp_r_49_88) + (sc_uint<1>) ~(comp_r_51_88) + (sc_uint<1>) ~(comp_r_53_88) + (sc_uint<1>) ~(comp_r_55_88) + (sc_uint<1>) ~(comp_r_57_88) + (sc_uint<1>) ~(comp_r_59_88) + (sc_uint<1>) ~(comp_r_61_88) + (sc_uint<1>) ~(comp_r_63_88) + (sc_uint<1>) ~(comp_r_65_88) + (sc_uint<1>) ~(comp_r_67_88) + (sc_uint<1>) ~(comp_r_69_88) + (sc_uint<1>) ~(comp_r_71_88) + (sc_uint<1>) ~(comp_r_73_88) + (sc_uint<1>) ~(comp_r_75_88) + (sc_uint<1>) ~(comp_r_77_88) + (sc_uint<1>) ~(comp_r_79_88) + (sc_uint<1>) ~(comp_r_81_88) + (sc_uint<1>) ~(comp_r_83_88) + (sc_uint<1>) ~(comp_r_85_88) + (sc_uint<1>) ~(comp_r_87_88);
    position[89] = 45 + (sc_uint<1>) ~(comp_r_1_89) + (sc_uint<1>) ~(comp_r_3_89) + (sc_uint<1>) ~(comp_r_5_89) + (sc_uint<1>) ~(comp_r_7_89) + (sc_uint<1>) ~(comp_r_9_89) + (sc_uint<1>) ~(comp_r_11_89) + (sc_uint<1>) ~(comp_r_13_89) + (sc_uint<1>) ~(comp_r_15_89) + (sc_uint<1>) ~(comp_r_17_89) + (sc_uint<1>) ~(comp_r_19_89) + (sc_uint<1>) ~(comp_r_21_89) + (sc_uint<1>) ~(comp_r_23_89) + (sc_uint<1>) ~(comp_r_25_89) + (sc_uint<1>) ~(comp_r_27_89) + (sc_uint<1>) ~(comp_r_29_89) + (sc_uint<1>) ~(comp_r_31_89) + (sc_uint<1>) ~(comp_r_33_89) + (sc_uint<1>) ~(comp_r_35_89) + (sc_uint<1>) ~(comp_r_37_89) + (sc_uint<1>) ~(comp_r_39_89) + (sc_uint<1>) ~(comp_r_41_89) + (sc_uint<1>) ~(comp_r_43_89) + (sc_uint<1>) ~(comp_r_45_89) + (sc_uint<1>) ~(comp_r_47_89) + (sc_uint<1>) ~(comp_r_49_89) + (sc_uint<1>) ~(comp_r_51_89) + (sc_uint<1>) ~(comp_r_53_89) + (sc_uint<1>) ~(comp_r_55_89) + (sc_uint<1>) ~(comp_r_57_89) + (sc_uint<1>) ~(comp_r_59_89) + (sc_uint<1>) ~(comp_r_61_89) + (sc_uint<1>) ~(comp_r_63_89) + (sc_uint<1>) ~(comp_r_65_89) + (sc_uint<1>) ~(comp_r_67_89) + (sc_uint<1>) ~(comp_r_69_89) + (sc_uint<1>) ~(comp_r_71_89) + (sc_uint<1>) ~(comp_r_73_89) + (sc_uint<1>) ~(comp_r_75_89) + (sc_uint<1>) ~(comp_r_77_89) + (sc_uint<1>) ~(comp_r_79_89) + (sc_uint<1>) ~(comp_r_81_89) + (sc_uint<1>) ~(comp_r_83_89) + (sc_uint<1>) ~(comp_r_85_89) + (sc_uint<1>) ~(comp_r_87_89) + comp_r_89_90 + comp_r_89_91 + comp_r_89_92 + comp_r_89_93 + comp_r_89_94 + comp_r_89_95 + comp_r_89_96 + comp_r_89_97 + comp_r_89_98 + comp_r_89_99 + comp_r_89_100 + comp_r_89_101 + comp_r_89_102 + comp_r_89_103 + comp_r_89_104 + comp_r_89_105 + comp_r_89_106 + comp_r_89_107 + comp_r_89_108 + comp_r_89_109 + comp_r_89_110 + comp_r_89_111 + comp_r_89_112 + comp_r_89_113 + comp_r_89_114 + comp_r_89_115 + comp_r_89_116 + comp_r_89_117 + comp_r_89_118 + comp_r_89_119 + comp_r_89_120 + comp_r_89_121 + comp_r_89_122 + comp_r_89_123 + comp_r_89_124 + comp_r_89_125 + comp_r_89_126;
    position[90] = 45 + (sc_uint<1>) ~(comp_r_1_90) + (sc_uint<1>) ~(comp_r_3_90) + (sc_uint<1>) ~(comp_r_5_90) + (sc_uint<1>) ~(comp_r_7_90) + (sc_uint<1>) ~(comp_r_9_90) + (sc_uint<1>) ~(comp_r_11_90) + (sc_uint<1>) ~(comp_r_13_90) + (sc_uint<1>) ~(comp_r_15_90) + (sc_uint<1>) ~(comp_r_17_90) + (sc_uint<1>) ~(comp_r_19_90) + (sc_uint<1>) ~(comp_r_21_90) + (sc_uint<1>) ~(comp_r_23_90) + (sc_uint<1>) ~(comp_r_25_90) + (sc_uint<1>) ~(comp_r_27_90) + (sc_uint<1>) ~(comp_r_29_90) + (sc_uint<1>) ~(comp_r_31_90) + (sc_uint<1>) ~(comp_r_33_90) + (sc_uint<1>) ~(comp_r_35_90) + (sc_uint<1>) ~(comp_r_37_90) + (sc_uint<1>) ~(comp_r_39_90) + (sc_uint<1>) ~(comp_r_41_90) + (sc_uint<1>) ~(comp_r_43_90) + (sc_uint<1>) ~(comp_r_45_90) + (sc_uint<1>) ~(comp_r_47_90) + (sc_uint<1>) ~(comp_r_49_90) + (sc_uint<1>) ~(comp_r_51_90) + (sc_uint<1>) ~(comp_r_53_90) + (sc_uint<1>) ~(comp_r_55_90) + (sc_uint<1>) ~(comp_r_57_90) + (sc_uint<1>) ~(comp_r_59_90) + (sc_uint<1>) ~(comp_r_61_90) + (sc_uint<1>) ~(comp_r_63_90) + (sc_uint<1>) ~(comp_r_65_90) + (sc_uint<1>) ~(comp_r_67_90) + (sc_uint<1>) ~(comp_r_69_90) + (sc_uint<1>) ~(comp_r_71_90) + (sc_uint<1>) ~(comp_r_73_90) + (sc_uint<1>) ~(comp_r_75_90) + (sc_uint<1>) ~(comp_r_77_90) + (sc_uint<1>) ~(comp_r_79_90) + (sc_uint<1>) ~(comp_r_81_90) + (sc_uint<1>) ~(comp_r_83_90) + (sc_uint<1>) ~(comp_r_85_90) + (sc_uint<1>) ~(comp_r_87_90) + (sc_uint<1>) ~(comp_r_89_90);
    position[91] = 46 + (sc_uint<1>) ~(comp_r_1_91) + (sc_uint<1>) ~(comp_r_3_91) + (sc_uint<1>) ~(comp_r_5_91) + (sc_uint<1>) ~(comp_r_7_91) + (sc_uint<1>) ~(comp_r_9_91) + (sc_uint<1>) ~(comp_r_11_91) + (sc_uint<1>) ~(comp_r_13_91) + (sc_uint<1>) ~(comp_r_15_91) + (sc_uint<1>) ~(comp_r_17_91) + (sc_uint<1>) ~(comp_r_19_91) + (sc_uint<1>) ~(comp_r_21_91) + (sc_uint<1>) ~(comp_r_23_91) + (sc_uint<1>) ~(comp_r_25_91) + (sc_uint<1>) ~(comp_r_27_91) + (sc_uint<1>) ~(comp_r_29_91) + (sc_uint<1>) ~(comp_r_31_91) + (sc_uint<1>) ~(comp_r_33_91) + (sc_uint<1>) ~(comp_r_35_91) + (sc_uint<1>) ~(comp_r_37_91) + (sc_uint<1>) ~(comp_r_39_91) + (sc_uint<1>) ~(comp_r_41_91) + (sc_uint<1>) ~(comp_r_43_91) + (sc_uint<1>) ~(comp_r_45_91) + (sc_uint<1>) ~(comp_r_47_91) + (sc_uint<1>) ~(comp_r_49_91) + (sc_uint<1>) ~(comp_r_51_91) + (sc_uint<1>) ~(comp_r_53_91) + (sc_uint<1>) ~(comp_r_55_91) + (sc_uint<1>) ~(comp_r_57_91) + (sc_uint<1>) ~(comp_r_59_91) + (sc_uint<1>) ~(comp_r_61_91) + (sc_uint<1>) ~(comp_r_63_91) + (sc_uint<1>) ~(comp_r_65_91) + (sc_uint<1>) ~(comp_r_67_91) + (sc_uint<1>) ~(comp_r_69_91) + (sc_uint<1>) ~(comp_r_71_91) + (sc_uint<1>) ~(comp_r_73_91) + (sc_uint<1>) ~(comp_r_75_91) + (sc_uint<1>) ~(comp_r_77_91) + (sc_uint<1>) ~(comp_r_79_91) + (sc_uint<1>) ~(comp_r_81_91) + (sc_uint<1>) ~(comp_r_83_91) + (sc_uint<1>) ~(comp_r_85_91) + (sc_uint<1>) ~(comp_r_87_91) + (sc_uint<1>) ~(comp_r_89_91) + comp_r_91_92 + comp_r_91_93 + comp_r_91_94 + comp_r_91_95 + comp_r_91_96 + comp_r_91_97 + comp_r_91_98 + comp_r_91_99 + comp_r_91_100 + comp_r_91_101 + comp_r_91_102 + comp_r_91_103 + comp_r_91_104 + comp_r_91_105 + comp_r_91_106 + comp_r_91_107 + comp_r_91_108 + comp_r_91_109 + comp_r_91_110 + comp_r_91_111 + comp_r_91_112 + comp_r_91_113 + comp_r_91_114 + comp_r_91_115 + comp_r_91_116 + comp_r_91_117 + comp_r_91_118 + comp_r_91_119 + comp_r_91_120 + comp_r_91_121 + comp_r_91_122 + comp_r_91_123 + comp_r_91_124 + comp_r_91_125 + comp_r_91_126;
    position[92] = 46 + (sc_uint<1>) ~(comp_r_1_92) + (sc_uint<1>) ~(comp_r_3_92) + (sc_uint<1>) ~(comp_r_5_92) + (sc_uint<1>) ~(comp_r_7_92) + (sc_uint<1>) ~(comp_r_9_92) + (sc_uint<1>) ~(comp_r_11_92) + (sc_uint<1>) ~(comp_r_13_92) + (sc_uint<1>) ~(comp_r_15_92) + (sc_uint<1>) ~(comp_r_17_92) + (sc_uint<1>) ~(comp_r_19_92) + (sc_uint<1>) ~(comp_r_21_92) + (sc_uint<1>) ~(comp_r_23_92) + (sc_uint<1>) ~(comp_r_25_92) + (sc_uint<1>) ~(comp_r_27_92) + (sc_uint<1>) ~(comp_r_29_92) + (sc_uint<1>) ~(comp_r_31_92) + (sc_uint<1>) ~(comp_r_33_92) + (sc_uint<1>) ~(comp_r_35_92) + (sc_uint<1>) ~(comp_r_37_92) + (sc_uint<1>) ~(comp_r_39_92) + (sc_uint<1>) ~(comp_r_41_92) + (sc_uint<1>) ~(comp_r_43_92) + (sc_uint<1>) ~(comp_r_45_92) + (sc_uint<1>) ~(comp_r_47_92) + (sc_uint<1>) ~(comp_r_49_92) + (sc_uint<1>) ~(comp_r_51_92) + (sc_uint<1>) ~(comp_r_53_92) + (sc_uint<1>) ~(comp_r_55_92) + (sc_uint<1>) ~(comp_r_57_92) + (sc_uint<1>) ~(comp_r_59_92) + (sc_uint<1>) ~(comp_r_61_92) + (sc_uint<1>) ~(comp_r_63_92) + (sc_uint<1>) ~(comp_r_65_92) + (sc_uint<1>) ~(comp_r_67_92) + (sc_uint<1>) ~(comp_r_69_92) + (sc_uint<1>) ~(comp_r_71_92) + (sc_uint<1>) ~(comp_r_73_92) + (sc_uint<1>) ~(comp_r_75_92) + (sc_uint<1>) ~(comp_r_77_92) + (sc_uint<1>) ~(comp_r_79_92) + (sc_uint<1>) ~(comp_r_81_92) + (sc_uint<1>) ~(comp_r_83_92) + (sc_uint<1>) ~(comp_r_85_92) + (sc_uint<1>) ~(comp_r_87_92) + (sc_uint<1>) ~(comp_r_89_92) + (sc_uint<1>) ~(comp_r_91_92);
    position[93] = 47 + (sc_uint<1>) ~(comp_r_1_93) + (sc_uint<1>) ~(comp_r_3_93) + (sc_uint<1>) ~(comp_r_5_93) + (sc_uint<1>) ~(comp_r_7_93) + (sc_uint<1>) ~(comp_r_9_93) + (sc_uint<1>) ~(comp_r_11_93) + (sc_uint<1>) ~(comp_r_13_93) + (sc_uint<1>) ~(comp_r_15_93) + (sc_uint<1>) ~(comp_r_17_93) + (sc_uint<1>) ~(comp_r_19_93) + (sc_uint<1>) ~(comp_r_21_93) + (sc_uint<1>) ~(comp_r_23_93) + (sc_uint<1>) ~(comp_r_25_93) + (sc_uint<1>) ~(comp_r_27_93) + (sc_uint<1>) ~(comp_r_29_93) + (sc_uint<1>) ~(comp_r_31_93) + (sc_uint<1>) ~(comp_r_33_93) + (sc_uint<1>) ~(comp_r_35_93) + (sc_uint<1>) ~(comp_r_37_93) + (sc_uint<1>) ~(comp_r_39_93) + (sc_uint<1>) ~(comp_r_41_93) + (sc_uint<1>) ~(comp_r_43_93) + (sc_uint<1>) ~(comp_r_45_93) + (sc_uint<1>) ~(comp_r_47_93) + (sc_uint<1>) ~(comp_r_49_93) + (sc_uint<1>) ~(comp_r_51_93) + (sc_uint<1>) ~(comp_r_53_93) + (sc_uint<1>) ~(comp_r_55_93) + (sc_uint<1>) ~(comp_r_57_93) + (sc_uint<1>) ~(comp_r_59_93) + (sc_uint<1>) ~(comp_r_61_93) + (sc_uint<1>) ~(comp_r_63_93) + (sc_uint<1>) ~(comp_r_65_93) + (sc_uint<1>) ~(comp_r_67_93) + (sc_uint<1>) ~(comp_r_69_93) + (sc_uint<1>) ~(comp_r_71_93) + (sc_uint<1>) ~(comp_r_73_93) + (sc_uint<1>) ~(comp_r_75_93) + (sc_uint<1>) ~(comp_r_77_93) + (sc_uint<1>) ~(comp_r_79_93) + (sc_uint<1>) ~(comp_r_81_93) + (sc_uint<1>) ~(comp_r_83_93) + (sc_uint<1>) ~(comp_r_85_93) + (sc_uint<1>) ~(comp_r_87_93) + (sc_uint<1>) ~(comp_r_89_93) + (sc_uint<1>) ~(comp_r_91_93) + comp_r_93_94 + comp_r_93_95 + comp_r_93_96 + comp_r_93_97 + comp_r_93_98 + comp_r_93_99 + comp_r_93_100 + comp_r_93_101 + comp_r_93_102 + comp_r_93_103 + comp_r_93_104 + comp_r_93_105 + comp_r_93_106 + comp_r_93_107 + comp_r_93_108 + comp_r_93_109 + comp_r_93_110 + comp_r_93_111 + comp_r_93_112 + comp_r_93_113 + comp_r_93_114 + comp_r_93_115 + comp_r_93_116 + comp_r_93_117 + comp_r_93_118 + comp_r_93_119 + comp_r_93_120 + comp_r_93_121 + comp_r_93_122 + comp_r_93_123 + comp_r_93_124 + comp_r_93_125 + comp_r_93_126;
    position[94] = 47 + (sc_uint<1>) ~(comp_r_1_94) + (sc_uint<1>) ~(comp_r_3_94) + (sc_uint<1>) ~(comp_r_5_94) + (sc_uint<1>) ~(comp_r_7_94) + (sc_uint<1>) ~(comp_r_9_94) + (sc_uint<1>) ~(comp_r_11_94) + (sc_uint<1>) ~(comp_r_13_94) + (sc_uint<1>) ~(comp_r_15_94) + (sc_uint<1>) ~(comp_r_17_94) + (sc_uint<1>) ~(comp_r_19_94) + (sc_uint<1>) ~(comp_r_21_94) + (sc_uint<1>) ~(comp_r_23_94) + (sc_uint<1>) ~(comp_r_25_94) + (sc_uint<1>) ~(comp_r_27_94) + (sc_uint<1>) ~(comp_r_29_94) + (sc_uint<1>) ~(comp_r_31_94) + (sc_uint<1>) ~(comp_r_33_94) + (sc_uint<1>) ~(comp_r_35_94) + (sc_uint<1>) ~(comp_r_37_94) + (sc_uint<1>) ~(comp_r_39_94) + (sc_uint<1>) ~(comp_r_41_94) + (sc_uint<1>) ~(comp_r_43_94) + (sc_uint<1>) ~(comp_r_45_94) + (sc_uint<1>) ~(comp_r_47_94) + (sc_uint<1>) ~(comp_r_49_94) + (sc_uint<1>) ~(comp_r_51_94) + (sc_uint<1>) ~(comp_r_53_94) + (sc_uint<1>) ~(comp_r_55_94) + (sc_uint<1>) ~(comp_r_57_94) + (sc_uint<1>) ~(comp_r_59_94) + (sc_uint<1>) ~(comp_r_61_94) + (sc_uint<1>) ~(comp_r_63_94) + (sc_uint<1>) ~(comp_r_65_94) + (sc_uint<1>) ~(comp_r_67_94) + (sc_uint<1>) ~(comp_r_69_94) + (sc_uint<1>) ~(comp_r_71_94) + (sc_uint<1>) ~(comp_r_73_94) + (sc_uint<1>) ~(comp_r_75_94) + (sc_uint<1>) ~(comp_r_77_94) + (sc_uint<1>) ~(comp_r_79_94) + (sc_uint<1>) ~(comp_r_81_94) + (sc_uint<1>) ~(comp_r_83_94) + (sc_uint<1>) ~(comp_r_85_94) + (sc_uint<1>) ~(comp_r_87_94) + (sc_uint<1>) ~(comp_r_89_94) + (sc_uint<1>) ~(comp_r_91_94) + (sc_uint<1>) ~(comp_r_93_94);
    position[95] = 48 + (sc_uint<1>) ~(comp_r_1_95) + (sc_uint<1>) ~(comp_r_3_95) + (sc_uint<1>) ~(comp_r_5_95) + (sc_uint<1>) ~(comp_r_7_95) + (sc_uint<1>) ~(comp_r_9_95) + (sc_uint<1>) ~(comp_r_11_95) + (sc_uint<1>) ~(comp_r_13_95) + (sc_uint<1>) ~(comp_r_15_95) + (sc_uint<1>) ~(comp_r_17_95) + (sc_uint<1>) ~(comp_r_19_95) + (sc_uint<1>) ~(comp_r_21_95) + (sc_uint<1>) ~(comp_r_23_95) + (sc_uint<1>) ~(comp_r_25_95) + (sc_uint<1>) ~(comp_r_27_95) + (sc_uint<1>) ~(comp_r_29_95) + (sc_uint<1>) ~(comp_r_31_95) + (sc_uint<1>) ~(comp_r_33_95) + (sc_uint<1>) ~(comp_r_35_95) + (sc_uint<1>) ~(comp_r_37_95) + (sc_uint<1>) ~(comp_r_39_95) + (sc_uint<1>) ~(comp_r_41_95) + (sc_uint<1>) ~(comp_r_43_95) + (sc_uint<1>) ~(comp_r_45_95) + (sc_uint<1>) ~(comp_r_47_95) + (sc_uint<1>) ~(comp_r_49_95) + (sc_uint<1>) ~(comp_r_51_95) + (sc_uint<1>) ~(comp_r_53_95) + (sc_uint<1>) ~(comp_r_55_95) + (sc_uint<1>) ~(comp_r_57_95) + (sc_uint<1>) ~(comp_r_59_95) + (sc_uint<1>) ~(comp_r_61_95) + (sc_uint<1>) ~(comp_r_63_95) + (sc_uint<1>) ~(comp_r_65_95) + (sc_uint<1>) ~(comp_r_67_95) + (sc_uint<1>) ~(comp_r_69_95) + (sc_uint<1>) ~(comp_r_71_95) + (sc_uint<1>) ~(comp_r_73_95) + (sc_uint<1>) ~(comp_r_75_95) + (sc_uint<1>) ~(comp_r_77_95) + (sc_uint<1>) ~(comp_r_79_95) + (sc_uint<1>) ~(comp_r_81_95) + (sc_uint<1>) ~(comp_r_83_95) + (sc_uint<1>) ~(comp_r_85_95) + (sc_uint<1>) ~(comp_r_87_95) + (sc_uint<1>) ~(comp_r_89_95) + (sc_uint<1>) ~(comp_r_91_95) + (sc_uint<1>) ~(comp_r_93_95) + comp_r_95_96 + comp_r_95_97 + comp_r_95_98 + comp_r_95_99 + comp_r_95_100 + comp_r_95_101 + comp_r_95_102 + comp_r_95_103 + comp_r_95_104 + comp_r_95_105 + comp_r_95_106 + comp_r_95_107 + comp_r_95_108 + comp_r_95_109 + comp_r_95_110 + comp_r_95_111 + comp_r_95_112 + comp_r_95_113 + comp_r_95_114 + comp_r_95_115 + comp_r_95_116 + comp_r_95_117 + comp_r_95_118 + comp_r_95_119 + comp_r_95_120 + comp_r_95_121 + comp_r_95_122 + comp_r_95_123 + comp_r_95_124 + comp_r_95_125 + comp_r_95_126;
    position[96] = 48 + (sc_uint<1>) ~(comp_r_1_96) + (sc_uint<1>) ~(comp_r_3_96) + (sc_uint<1>) ~(comp_r_5_96) + (sc_uint<1>) ~(comp_r_7_96) + (sc_uint<1>) ~(comp_r_9_96) + (sc_uint<1>) ~(comp_r_11_96) + (sc_uint<1>) ~(comp_r_13_96) + (sc_uint<1>) ~(comp_r_15_96) + (sc_uint<1>) ~(comp_r_17_96) + (sc_uint<1>) ~(comp_r_19_96) + (sc_uint<1>) ~(comp_r_21_96) + (sc_uint<1>) ~(comp_r_23_96) + (sc_uint<1>) ~(comp_r_25_96) + (sc_uint<1>) ~(comp_r_27_96) + (sc_uint<1>) ~(comp_r_29_96) + (sc_uint<1>) ~(comp_r_31_96) + (sc_uint<1>) ~(comp_r_33_96) + (sc_uint<1>) ~(comp_r_35_96) + (sc_uint<1>) ~(comp_r_37_96) + (sc_uint<1>) ~(comp_r_39_96) + (sc_uint<1>) ~(comp_r_41_96) + (sc_uint<1>) ~(comp_r_43_96) + (sc_uint<1>) ~(comp_r_45_96) + (sc_uint<1>) ~(comp_r_47_96) + (sc_uint<1>) ~(comp_r_49_96) + (sc_uint<1>) ~(comp_r_51_96) + (sc_uint<1>) ~(comp_r_53_96) + (sc_uint<1>) ~(comp_r_55_96) + (sc_uint<1>) ~(comp_r_57_96) + (sc_uint<1>) ~(comp_r_59_96) + (sc_uint<1>) ~(comp_r_61_96) + (sc_uint<1>) ~(comp_r_63_96) + (sc_uint<1>) ~(comp_r_65_96) + (sc_uint<1>) ~(comp_r_67_96) + (sc_uint<1>) ~(comp_r_69_96) + (sc_uint<1>) ~(comp_r_71_96) + (sc_uint<1>) ~(comp_r_73_96) + (sc_uint<1>) ~(comp_r_75_96) + (sc_uint<1>) ~(comp_r_77_96) + (sc_uint<1>) ~(comp_r_79_96) + (sc_uint<1>) ~(comp_r_81_96) + (sc_uint<1>) ~(comp_r_83_96) + (sc_uint<1>) ~(comp_r_85_96) + (sc_uint<1>) ~(comp_r_87_96) + (sc_uint<1>) ~(comp_r_89_96) + (sc_uint<1>) ~(comp_r_91_96) + (sc_uint<1>) ~(comp_r_93_96) + (sc_uint<1>) ~(comp_r_95_96);
    position[97] = 49 + (sc_uint<1>) ~(comp_r_1_97) + (sc_uint<1>) ~(comp_r_3_97) + (sc_uint<1>) ~(comp_r_5_97) + (sc_uint<1>) ~(comp_r_7_97) + (sc_uint<1>) ~(comp_r_9_97) + (sc_uint<1>) ~(comp_r_11_97) + (sc_uint<1>) ~(comp_r_13_97) + (sc_uint<1>) ~(comp_r_15_97) + (sc_uint<1>) ~(comp_r_17_97) + (sc_uint<1>) ~(comp_r_19_97) + (sc_uint<1>) ~(comp_r_21_97) + (sc_uint<1>) ~(comp_r_23_97) + (sc_uint<1>) ~(comp_r_25_97) + (sc_uint<1>) ~(comp_r_27_97) + (sc_uint<1>) ~(comp_r_29_97) + (sc_uint<1>) ~(comp_r_31_97) + (sc_uint<1>) ~(comp_r_33_97) + (sc_uint<1>) ~(comp_r_35_97) + (sc_uint<1>) ~(comp_r_37_97) + (sc_uint<1>) ~(comp_r_39_97) + (sc_uint<1>) ~(comp_r_41_97) + (sc_uint<1>) ~(comp_r_43_97) + (sc_uint<1>) ~(comp_r_45_97) + (sc_uint<1>) ~(comp_r_47_97) + (sc_uint<1>) ~(comp_r_49_97) + (sc_uint<1>) ~(comp_r_51_97) + (sc_uint<1>) ~(comp_r_53_97) + (sc_uint<1>) ~(comp_r_55_97) + (sc_uint<1>) ~(comp_r_57_97) + (sc_uint<1>) ~(comp_r_59_97) + (sc_uint<1>) ~(comp_r_61_97) + (sc_uint<1>) ~(comp_r_63_97) + (sc_uint<1>) ~(comp_r_65_97) + (sc_uint<1>) ~(comp_r_67_97) + (sc_uint<1>) ~(comp_r_69_97) + (sc_uint<1>) ~(comp_r_71_97) + (sc_uint<1>) ~(comp_r_73_97) + (sc_uint<1>) ~(comp_r_75_97) + (sc_uint<1>) ~(comp_r_77_97) + (sc_uint<1>) ~(comp_r_79_97) + (sc_uint<1>) ~(comp_r_81_97) + (sc_uint<1>) ~(comp_r_83_97) + (sc_uint<1>) ~(comp_r_85_97) + (sc_uint<1>) ~(comp_r_87_97) + (sc_uint<1>) ~(comp_r_89_97) + (sc_uint<1>) ~(comp_r_91_97) + (sc_uint<1>) ~(comp_r_93_97) + (sc_uint<1>) ~(comp_r_95_97) + comp_r_97_98 + comp_r_97_99 + comp_r_97_100 + comp_r_97_101 + comp_r_97_102 + comp_r_97_103 + comp_r_97_104 + comp_r_97_105 + comp_r_97_106 + comp_r_97_107 + comp_r_97_108 + comp_r_97_109 + comp_r_97_110 + comp_r_97_111 + comp_r_97_112 + comp_r_97_113 + comp_r_97_114 + comp_r_97_115 + comp_r_97_116 + comp_r_97_117 + comp_r_97_118 + comp_r_97_119 + comp_r_97_120 + comp_r_97_121 + comp_r_97_122 + comp_r_97_123 + comp_r_97_124 + comp_r_97_125 + comp_r_97_126;
    position[98] = 49 + (sc_uint<1>) ~(comp_r_1_98) + (sc_uint<1>) ~(comp_r_3_98) + (sc_uint<1>) ~(comp_r_5_98) + (sc_uint<1>) ~(comp_r_7_98) + (sc_uint<1>) ~(comp_r_9_98) + (sc_uint<1>) ~(comp_r_11_98) + (sc_uint<1>) ~(comp_r_13_98) + (sc_uint<1>) ~(comp_r_15_98) + (sc_uint<1>) ~(comp_r_17_98) + (sc_uint<1>) ~(comp_r_19_98) + (sc_uint<1>) ~(comp_r_21_98) + (sc_uint<1>) ~(comp_r_23_98) + (sc_uint<1>) ~(comp_r_25_98) + (sc_uint<1>) ~(comp_r_27_98) + (sc_uint<1>) ~(comp_r_29_98) + (sc_uint<1>) ~(comp_r_31_98) + (sc_uint<1>) ~(comp_r_33_98) + (sc_uint<1>) ~(comp_r_35_98) + (sc_uint<1>) ~(comp_r_37_98) + (sc_uint<1>) ~(comp_r_39_98) + (sc_uint<1>) ~(comp_r_41_98) + (sc_uint<1>) ~(comp_r_43_98) + (sc_uint<1>) ~(comp_r_45_98) + (sc_uint<1>) ~(comp_r_47_98) + (sc_uint<1>) ~(comp_r_49_98) + (sc_uint<1>) ~(comp_r_51_98) + (sc_uint<1>) ~(comp_r_53_98) + (sc_uint<1>) ~(comp_r_55_98) + (sc_uint<1>) ~(comp_r_57_98) + (sc_uint<1>) ~(comp_r_59_98) + (sc_uint<1>) ~(comp_r_61_98) + (sc_uint<1>) ~(comp_r_63_98) + (sc_uint<1>) ~(comp_r_65_98) + (sc_uint<1>) ~(comp_r_67_98) + (sc_uint<1>) ~(comp_r_69_98) + (sc_uint<1>) ~(comp_r_71_98) + (sc_uint<1>) ~(comp_r_73_98) + (sc_uint<1>) ~(comp_r_75_98) + (sc_uint<1>) ~(comp_r_77_98) + (sc_uint<1>) ~(comp_r_79_98) + (sc_uint<1>) ~(comp_r_81_98) + (sc_uint<1>) ~(comp_r_83_98) + (sc_uint<1>) ~(comp_r_85_98) + (sc_uint<1>) ~(comp_r_87_98) + (sc_uint<1>) ~(comp_r_89_98) + (sc_uint<1>) ~(comp_r_91_98) + (sc_uint<1>) ~(comp_r_93_98) + (sc_uint<1>) ~(comp_r_95_98) + (sc_uint<1>) ~(comp_r_97_98);
    position[99] = 50 + (sc_uint<1>) ~(comp_r_1_99) + (sc_uint<1>) ~(comp_r_3_99) + (sc_uint<1>) ~(comp_r_5_99) + (sc_uint<1>) ~(comp_r_7_99) + (sc_uint<1>) ~(comp_r_9_99) + (sc_uint<1>) ~(comp_r_11_99) + (sc_uint<1>) ~(comp_r_13_99) + (sc_uint<1>) ~(comp_r_15_99) + (sc_uint<1>) ~(comp_r_17_99) + (sc_uint<1>) ~(comp_r_19_99) + (sc_uint<1>) ~(comp_r_21_99) + (sc_uint<1>) ~(comp_r_23_99) + (sc_uint<1>) ~(comp_r_25_99) + (sc_uint<1>) ~(comp_r_27_99) + (sc_uint<1>) ~(comp_r_29_99) + (sc_uint<1>) ~(comp_r_31_99) + (sc_uint<1>) ~(comp_r_33_99) + (sc_uint<1>) ~(comp_r_35_99) + (sc_uint<1>) ~(comp_r_37_99) + (sc_uint<1>) ~(comp_r_39_99) + (sc_uint<1>) ~(comp_r_41_99) + (sc_uint<1>) ~(comp_r_43_99) + (sc_uint<1>) ~(comp_r_45_99) + (sc_uint<1>) ~(comp_r_47_99) + (sc_uint<1>) ~(comp_r_49_99) + (sc_uint<1>) ~(comp_r_51_99) + (sc_uint<1>) ~(comp_r_53_99) + (sc_uint<1>) ~(comp_r_55_99) + (sc_uint<1>) ~(comp_r_57_99) + (sc_uint<1>) ~(comp_r_59_99) + (sc_uint<1>) ~(comp_r_61_99) + (sc_uint<1>) ~(comp_r_63_99) + (sc_uint<1>) ~(comp_r_65_99) + (sc_uint<1>) ~(comp_r_67_99) + (sc_uint<1>) ~(comp_r_69_99) + (sc_uint<1>) ~(comp_r_71_99) + (sc_uint<1>) ~(comp_r_73_99) + (sc_uint<1>) ~(comp_r_75_99) + (sc_uint<1>) ~(comp_r_77_99) + (sc_uint<1>) ~(comp_r_79_99) + (sc_uint<1>) ~(comp_r_81_99) + (sc_uint<1>) ~(comp_r_83_99) + (sc_uint<1>) ~(comp_r_85_99) + (sc_uint<1>) ~(comp_r_87_99) + (sc_uint<1>) ~(comp_r_89_99) + (sc_uint<1>) ~(comp_r_91_99) + (sc_uint<1>) ~(comp_r_93_99) + (sc_uint<1>) ~(comp_r_95_99) + (sc_uint<1>) ~(comp_r_97_99) + comp_r_99_100 + comp_r_99_101 + comp_r_99_102 + comp_r_99_103 + comp_r_99_104 + comp_r_99_105 + comp_r_99_106 + comp_r_99_107 + comp_r_99_108 + comp_r_99_109 + comp_r_99_110 + comp_r_99_111 + comp_r_99_112 + comp_r_99_113 + comp_r_99_114 + comp_r_99_115 + comp_r_99_116 + comp_r_99_117 + comp_r_99_118 + comp_r_99_119 + comp_r_99_120 + comp_r_99_121 + comp_r_99_122 + comp_r_99_123 + comp_r_99_124 + comp_r_99_125 + comp_r_99_126;
    position[100] = 50 + (sc_uint<1>) ~(comp_r_1_100) + (sc_uint<1>) ~(comp_r_3_100) + (sc_uint<1>) ~(comp_r_5_100) + (sc_uint<1>) ~(comp_r_7_100) + (sc_uint<1>) ~(comp_r_9_100) + (sc_uint<1>) ~(comp_r_11_100) + (sc_uint<1>) ~(comp_r_13_100) + (sc_uint<1>) ~(comp_r_15_100) + (sc_uint<1>) ~(comp_r_17_100) + (sc_uint<1>) ~(comp_r_19_100) + (sc_uint<1>) ~(comp_r_21_100) + (sc_uint<1>) ~(comp_r_23_100) + (sc_uint<1>) ~(comp_r_25_100) + (sc_uint<1>) ~(comp_r_27_100) + (sc_uint<1>) ~(comp_r_29_100) + (sc_uint<1>) ~(comp_r_31_100) + (sc_uint<1>) ~(comp_r_33_100) + (sc_uint<1>) ~(comp_r_35_100) + (sc_uint<1>) ~(comp_r_37_100) + (sc_uint<1>) ~(comp_r_39_100) + (sc_uint<1>) ~(comp_r_41_100) + (sc_uint<1>) ~(comp_r_43_100) + (sc_uint<1>) ~(comp_r_45_100) + (sc_uint<1>) ~(comp_r_47_100) + (sc_uint<1>) ~(comp_r_49_100) + (sc_uint<1>) ~(comp_r_51_100) + (sc_uint<1>) ~(comp_r_53_100) + (sc_uint<1>) ~(comp_r_55_100) + (sc_uint<1>) ~(comp_r_57_100) + (sc_uint<1>) ~(comp_r_59_100) + (sc_uint<1>) ~(comp_r_61_100) + (sc_uint<1>) ~(comp_r_63_100) + (sc_uint<1>) ~(comp_r_65_100) + (sc_uint<1>) ~(comp_r_67_100) + (sc_uint<1>) ~(comp_r_69_100) + (sc_uint<1>) ~(comp_r_71_100) + (sc_uint<1>) ~(comp_r_73_100) + (sc_uint<1>) ~(comp_r_75_100) + (sc_uint<1>) ~(comp_r_77_100) + (sc_uint<1>) ~(comp_r_79_100) + (sc_uint<1>) ~(comp_r_81_100) + (sc_uint<1>) ~(comp_r_83_100) + (sc_uint<1>) ~(comp_r_85_100) + (sc_uint<1>) ~(comp_r_87_100) + (sc_uint<1>) ~(comp_r_89_100) + (sc_uint<1>) ~(comp_r_91_100) + (sc_uint<1>) ~(comp_r_93_100) + (sc_uint<1>) ~(comp_r_95_100) + (sc_uint<1>) ~(comp_r_97_100) + (sc_uint<1>) ~(comp_r_99_100);
    position[101] = 51 + (sc_uint<1>) ~(comp_r_1_101) + (sc_uint<1>) ~(comp_r_3_101) + (sc_uint<1>) ~(comp_r_5_101) + (sc_uint<1>) ~(comp_r_7_101) + (sc_uint<1>) ~(comp_r_9_101) + (sc_uint<1>) ~(comp_r_11_101) + (sc_uint<1>) ~(comp_r_13_101) + (sc_uint<1>) ~(comp_r_15_101) + (sc_uint<1>) ~(comp_r_17_101) + (sc_uint<1>) ~(comp_r_19_101) + (sc_uint<1>) ~(comp_r_21_101) + (sc_uint<1>) ~(comp_r_23_101) + (sc_uint<1>) ~(comp_r_25_101) + (sc_uint<1>) ~(comp_r_27_101) + (sc_uint<1>) ~(comp_r_29_101) + (sc_uint<1>) ~(comp_r_31_101) + (sc_uint<1>) ~(comp_r_33_101) + (sc_uint<1>) ~(comp_r_35_101) + (sc_uint<1>) ~(comp_r_37_101) + (sc_uint<1>) ~(comp_r_39_101) + (sc_uint<1>) ~(comp_r_41_101) + (sc_uint<1>) ~(comp_r_43_101) + (sc_uint<1>) ~(comp_r_45_101) + (sc_uint<1>) ~(comp_r_47_101) + (sc_uint<1>) ~(comp_r_49_101) + (sc_uint<1>) ~(comp_r_51_101) + (sc_uint<1>) ~(comp_r_53_101) + (sc_uint<1>) ~(comp_r_55_101) + (sc_uint<1>) ~(comp_r_57_101) + (sc_uint<1>) ~(comp_r_59_101) + (sc_uint<1>) ~(comp_r_61_101) + (sc_uint<1>) ~(comp_r_63_101) + (sc_uint<1>) ~(comp_r_65_101) + (sc_uint<1>) ~(comp_r_67_101) + (sc_uint<1>) ~(comp_r_69_101) + (sc_uint<1>) ~(comp_r_71_101) + (sc_uint<1>) ~(comp_r_73_101) + (sc_uint<1>) ~(comp_r_75_101) + (sc_uint<1>) ~(comp_r_77_101) + (sc_uint<1>) ~(comp_r_79_101) + (sc_uint<1>) ~(comp_r_81_101) + (sc_uint<1>) ~(comp_r_83_101) + (sc_uint<1>) ~(comp_r_85_101) + (sc_uint<1>) ~(comp_r_87_101) + (sc_uint<1>) ~(comp_r_89_101) + (sc_uint<1>) ~(comp_r_91_101) + (sc_uint<1>) ~(comp_r_93_101) + (sc_uint<1>) ~(comp_r_95_101) + (sc_uint<1>) ~(comp_r_97_101) + (sc_uint<1>) ~(comp_r_99_101) + comp_r_101_102 + comp_r_101_103 + comp_r_101_104 + comp_r_101_105 + comp_r_101_106 + comp_r_101_107 + comp_r_101_108 + comp_r_101_109 + comp_r_101_110 + comp_r_101_111 + comp_r_101_112 + comp_r_101_113 + comp_r_101_114 + comp_r_101_115 + comp_r_101_116 + comp_r_101_117 + comp_r_101_118 + comp_r_101_119 + comp_r_101_120 + comp_r_101_121 + comp_r_101_122 + comp_r_101_123 + comp_r_101_124 + comp_r_101_125 + comp_r_101_126;
    position[102] = 51 + (sc_uint<1>) ~(comp_r_1_102) + (sc_uint<1>) ~(comp_r_3_102) + (sc_uint<1>) ~(comp_r_5_102) + (sc_uint<1>) ~(comp_r_7_102) + (sc_uint<1>) ~(comp_r_9_102) + (sc_uint<1>) ~(comp_r_11_102) + (sc_uint<1>) ~(comp_r_13_102) + (sc_uint<1>) ~(comp_r_15_102) + (sc_uint<1>) ~(comp_r_17_102) + (sc_uint<1>) ~(comp_r_19_102) + (sc_uint<1>) ~(comp_r_21_102) + (sc_uint<1>) ~(comp_r_23_102) + (sc_uint<1>) ~(comp_r_25_102) + (sc_uint<1>) ~(comp_r_27_102) + (sc_uint<1>) ~(comp_r_29_102) + (sc_uint<1>) ~(comp_r_31_102) + (sc_uint<1>) ~(comp_r_33_102) + (sc_uint<1>) ~(comp_r_35_102) + (sc_uint<1>) ~(comp_r_37_102) + (sc_uint<1>) ~(comp_r_39_102) + (sc_uint<1>) ~(comp_r_41_102) + (sc_uint<1>) ~(comp_r_43_102) + (sc_uint<1>) ~(comp_r_45_102) + (sc_uint<1>) ~(comp_r_47_102) + (sc_uint<1>) ~(comp_r_49_102) + (sc_uint<1>) ~(comp_r_51_102) + (sc_uint<1>) ~(comp_r_53_102) + (sc_uint<1>) ~(comp_r_55_102) + (sc_uint<1>) ~(comp_r_57_102) + (sc_uint<1>) ~(comp_r_59_102) + (sc_uint<1>) ~(comp_r_61_102) + (sc_uint<1>) ~(comp_r_63_102) + (sc_uint<1>) ~(comp_r_65_102) + (sc_uint<1>) ~(comp_r_67_102) + (sc_uint<1>) ~(comp_r_69_102) + (sc_uint<1>) ~(comp_r_71_102) + (sc_uint<1>) ~(comp_r_73_102) + (sc_uint<1>) ~(comp_r_75_102) + (sc_uint<1>) ~(comp_r_77_102) + (sc_uint<1>) ~(comp_r_79_102) + (sc_uint<1>) ~(comp_r_81_102) + (sc_uint<1>) ~(comp_r_83_102) + (sc_uint<1>) ~(comp_r_85_102) + (sc_uint<1>) ~(comp_r_87_102) + (sc_uint<1>) ~(comp_r_89_102) + (sc_uint<1>) ~(comp_r_91_102) + (sc_uint<1>) ~(comp_r_93_102) + (sc_uint<1>) ~(comp_r_95_102) + (sc_uint<1>) ~(comp_r_97_102) + (sc_uint<1>) ~(comp_r_99_102) + (sc_uint<1>) ~(comp_r_101_102);
    position[103] = 52 + (sc_uint<1>) ~(comp_r_1_103) + (sc_uint<1>) ~(comp_r_3_103) + (sc_uint<1>) ~(comp_r_5_103) + (sc_uint<1>) ~(comp_r_7_103) + (sc_uint<1>) ~(comp_r_9_103) + (sc_uint<1>) ~(comp_r_11_103) + (sc_uint<1>) ~(comp_r_13_103) + (sc_uint<1>) ~(comp_r_15_103) + (sc_uint<1>) ~(comp_r_17_103) + (sc_uint<1>) ~(comp_r_19_103) + (sc_uint<1>) ~(comp_r_21_103) + (sc_uint<1>) ~(comp_r_23_103) + (sc_uint<1>) ~(comp_r_25_103) + (sc_uint<1>) ~(comp_r_27_103) + (sc_uint<1>) ~(comp_r_29_103) + (sc_uint<1>) ~(comp_r_31_103) + (sc_uint<1>) ~(comp_r_33_103) + (sc_uint<1>) ~(comp_r_35_103) + (sc_uint<1>) ~(comp_r_37_103) + (sc_uint<1>) ~(comp_r_39_103) + (sc_uint<1>) ~(comp_r_41_103) + (sc_uint<1>) ~(comp_r_43_103) + (sc_uint<1>) ~(comp_r_45_103) + (sc_uint<1>) ~(comp_r_47_103) + (sc_uint<1>) ~(comp_r_49_103) + (sc_uint<1>) ~(comp_r_51_103) + (sc_uint<1>) ~(comp_r_53_103) + (sc_uint<1>) ~(comp_r_55_103) + (sc_uint<1>) ~(comp_r_57_103) + (sc_uint<1>) ~(comp_r_59_103) + (sc_uint<1>) ~(comp_r_61_103) + (sc_uint<1>) ~(comp_r_63_103) + (sc_uint<1>) ~(comp_r_65_103) + (sc_uint<1>) ~(comp_r_67_103) + (sc_uint<1>) ~(comp_r_69_103) + (sc_uint<1>) ~(comp_r_71_103) + (sc_uint<1>) ~(comp_r_73_103) + (sc_uint<1>) ~(comp_r_75_103) + (sc_uint<1>) ~(comp_r_77_103) + (sc_uint<1>) ~(comp_r_79_103) + (sc_uint<1>) ~(comp_r_81_103) + (sc_uint<1>) ~(comp_r_83_103) + (sc_uint<1>) ~(comp_r_85_103) + (sc_uint<1>) ~(comp_r_87_103) + (sc_uint<1>) ~(comp_r_89_103) + (sc_uint<1>) ~(comp_r_91_103) + (sc_uint<1>) ~(comp_r_93_103) + (sc_uint<1>) ~(comp_r_95_103) + (sc_uint<1>) ~(comp_r_97_103) + (sc_uint<1>) ~(comp_r_99_103) + (sc_uint<1>) ~(comp_r_101_103) + comp_r_103_104 + comp_r_103_105 + comp_r_103_106 + comp_r_103_107 + comp_r_103_108 + comp_r_103_109 + comp_r_103_110 + comp_r_103_111 + comp_r_103_112 + comp_r_103_113 + comp_r_103_114 + comp_r_103_115 + comp_r_103_116 + comp_r_103_117 + comp_r_103_118 + comp_r_103_119 + comp_r_103_120 + comp_r_103_121 + comp_r_103_122 + comp_r_103_123 + comp_r_103_124 + comp_r_103_125 + comp_r_103_126;
    position[104] = 52 + (sc_uint<1>) ~(comp_r_1_104) + (sc_uint<1>) ~(comp_r_3_104) + (sc_uint<1>) ~(comp_r_5_104) + (sc_uint<1>) ~(comp_r_7_104) + (sc_uint<1>) ~(comp_r_9_104) + (sc_uint<1>) ~(comp_r_11_104) + (sc_uint<1>) ~(comp_r_13_104) + (sc_uint<1>) ~(comp_r_15_104) + (sc_uint<1>) ~(comp_r_17_104) + (sc_uint<1>) ~(comp_r_19_104) + (sc_uint<1>) ~(comp_r_21_104) + (sc_uint<1>) ~(comp_r_23_104) + (sc_uint<1>) ~(comp_r_25_104) + (sc_uint<1>) ~(comp_r_27_104) + (sc_uint<1>) ~(comp_r_29_104) + (sc_uint<1>) ~(comp_r_31_104) + (sc_uint<1>) ~(comp_r_33_104) + (sc_uint<1>) ~(comp_r_35_104) + (sc_uint<1>) ~(comp_r_37_104) + (sc_uint<1>) ~(comp_r_39_104) + (sc_uint<1>) ~(comp_r_41_104) + (sc_uint<1>) ~(comp_r_43_104) + (sc_uint<1>) ~(comp_r_45_104) + (sc_uint<1>) ~(comp_r_47_104) + (sc_uint<1>) ~(comp_r_49_104) + (sc_uint<1>) ~(comp_r_51_104) + (sc_uint<1>) ~(comp_r_53_104) + (sc_uint<1>) ~(comp_r_55_104) + (sc_uint<1>) ~(comp_r_57_104) + (sc_uint<1>) ~(comp_r_59_104) + (sc_uint<1>) ~(comp_r_61_104) + (sc_uint<1>) ~(comp_r_63_104) + (sc_uint<1>) ~(comp_r_65_104) + (sc_uint<1>) ~(comp_r_67_104) + (sc_uint<1>) ~(comp_r_69_104) + (sc_uint<1>) ~(comp_r_71_104) + (sc_uint<1>) ~(comp_r_73_104) + (sc_uint<1>) ~(comp_r_75_104) + (sc_uint<1>) ~(comp_r_77_104) + (sc_uint<1>) ~(comp_r_79_104) + (sc_uint<1>) ~(comp_r_81_104) + (sc_uint<1>) ~(comp_r_83_104) + (sc_uint<1>) ~(comp_r_85_104) + (sc_uint<1>) ~(comp_r_87_104) + (sc_uint<1>) ~(comp_r_89_104) + (sc_uint<1>) ~(comp_r_91_104) + (sc_uint<1>) ~(comp_r_93_104) + (sc_uint<1>) ~(comp_r_95_104) + (sc_uint<1>) ~(comp_r_97_104) + (sc_uint<1>) ~(comp_r_99_104) + (sc_uint<1>) ~(comp_r_101_104) + (sc_uint<1>) ~(comp_r_103_104);
    position[105] = 53 + (sc_uint<1>) ~(comp_r_1_105) + (sc_uint<1>) ~(comp_r_3_105) + (sc_uint<1>) ~(comp_r_5_105) + (sc_uint<1>) ~(comp_r_7_105) + (sc_uint<1>) ~(comp_r_9_105) + (sc_uint<1>) ~(comp_r_11_105) + (sc_uint<1>) ~(comp_r_13_105) + (sc_uint<1>) ~(comp_r_15_105) + (sc_uint<1>) ~(comp_r_17_105) + (sc_uint<1>) ~(comp_r_19_105) + (sc_uint<1>) ~(comp_r_21_105) + (sc_uint<1>) ~(comp_r_23_105) + (sc_uint<1>) ~(comp_r_25_105) + (sc_uint<1>) ~(comp_r_27_105) + (sc_uint<1>) ~(comp_r_29_105) + (sc_uint<1>) ~(comp_r_31_105) + (sc_uint<1>) ~(comp_r_33_105) + (sc_uint<1>) ~(comp_r_35_105) + (sc_uint<1>) ~(comp_r_37_105) + (sc_uint<1>) ~(comp_r_39_105) + (sc_uint<1>) ~(comp_r_41_105) + (sc_uint<1>) ~(comp_r_43_105) + (sc_uint<1>) ~(comp_r_45_105) + (sc_uint<1>) ~(comp_r_47_105) + (sc_uint<1>) ~(comp_r_49_105) + (sc_uint<1>) ~(comp_r_51_105) + (sc_uint<1>) ~(comp_r_53_105) + (sc_uint<1>) ~(comp_r_55_105) + (sc_uint<1>) ~(comp_r_57_105) + (sc_uint<1>) ~(comp_r_59_105) + (sc_uint<1>) ~(comp_r_61_105) + (sc_uint<1>) ~(comp_r_63_105) + (sc_uint<1>) ~(comp_r_65_105) + (sc_uint<1>) ~(comp_r_67_105) + (sc_uint<1>) ~(comp_r_69_105) + (sc_uint<1>) ~(comp_r_71_105) + (sc_uint<1>) ~(comp_r_73_105) + (sc_uint<1>) ~(comp_r_75_105) + (sc_uint<1>) ~(comp_r_77_105) + (sc_uint<1>) ~(comp_r_79_105) + (sc_uint<1>) ~(comp_r_81_105) + (sc_uint<1>) ~(comp_r_83_105) + (sc_uint<1>) ~(comp_r_85_105) + (sc_uint<1>) ~(comp_r_87_105) + (sc_uint<1>) ~(comp_r_89_105) + (sc_uint<1>) ~(comp_r_91_105) + (sc_uint<1>) ~(comp_r_93_105) + (sc_uint<1>) ~(comp_r_95_105) + (sc_uint<1>) ~(comp_r_97_105) + (sc_uint<1>) ~(comp_r_99_105) + (sc_uint<1>) ~(comp_r_101_105) + (sc_uint<1>) ~(comp_r_103_105) + comp_r_105_106 + comp_r_105_107 + comp_r_105_108 + comp_r_105_109 + comp_r_105_110 + comp_r_105_111 + comp_r_105_112 + comp_r_105_113 + comp_r_105_114 + comp_r_105_115 + comp_r_105_116 + comp_r_105_117 + comp_r_105_118 + comp_r_105_119 + comp_r_105_120 + comp_r_105_121 + comp_r_105_122 + comp_r_105_123 + comp_r_105_124 + comp_r_105_125 + comp_r_105_126;
    position[106] = 53 + (sc_uint<1>) ~(comp_r_1_106) + (sc_uint<1>) ~(comp_r_3_106) + (sc_uint<1>) ~(comp_r_5_106) + (sc_uint<1>) ~(comp_r_7_106) + (sc_uint<1>) ~(comp_r_9_106) + (sc_uint<1>) ~(comp_r_11_106) + (sc_uint<1>) ~(comp_r_13_106) + (sc_uint<1>) ~(comp_r_15_106) + (sc_uint<1>) ~(comp_r_17_106) + (sc_uint<1>) ~(comp_r_19_106) + (sc_uint<1>) ~(comp_r_21_106) + (sc_uint<1>) ~(comp_r_23_106) + (sc_uint<1>) ~(comp_r_25_106) + (sc_uint<1>) ~(comp_r_27_106) + (sc_uint<1>) ~(comp_r_29_106) + (sc_uint<1>) ~(comp_r_31_106) + (sc_uint<1>) ~(comp_r_33_106) + (sc_uint<1>) ~(comp_r_35_106) + (sc_uint<1>) ~(comp_r_37_106) + (sc_uint<1>) ~(comp_r_39_106) + (sc_uint<1>) ~(comp_r_41_106) + (sc_uint<1>) ~(comp_r_43_106) + (sc_uint<1>) ~(comp_r_45_106) + (sc_uint<1>) ~(comp_r_47_106) + (sc_uint<1>) ~(comp_r_49_106) + (sc_uint<1>) ~(comp_r_51_106) + (sc_uint<1>) ~(comp_r_53_106) + (sc_uint<1>) ~(comp_r_55_106) + (sc_uint<1>) ~(comp_r_57_106) + (sc_uint<1>) ~(comp_r_59_106) + (sc_uint<1>) ~(comp_r_61_106) + (sc_uint<1>) ~(comp_r_63_106) + (sc_uint<1>) ~(comp_r_65_106) + (sc_uint<1>) ~(comp_r_67_106) + (sc_uint<1>) ~(comp_r_69_106) + (sc_uint<1>) ~(comp_r_71_106) + (sc_uint<1>) ~(comp_r_73_106) + (sc_uint<1>) ~(comp_r_75_106) + (sc_uint<1>) ~(comp_r_77_106) + (sc_uint<1>) ~(comp_r_79_106) + (sc_uint<1>) ~(comp_r_81_106) + (sc_uint<1>) ~(comp_r_83_106) + (sc_uint<1>) ~(comp_r_85_106) + (sc_uint<1>) ~(comp_r_87_106) + (sc_uint<1>) ~(comp_r_89_106) + (sc_uint<1>) ~(comp_r_91_106) + (sc_uint<1>) ~(comp_r_93_106) + (sc_uint<1>) ~(comp_r_95_106) + (sc_uint<1>) ~(comp_r_97_106) + (sc_uint<1>) ~(comp_r_99_106) + (sc_uint<1>) ~(comp_r_101_106) + (sc_uint<1>) ~(comp_r_103_106) + (sc_uint<1>) ~(comp_r_105_106);
    position[107] = 54 + (sc_uint<1>) ~(comp_r_1_107) + (sc_uint<1>) ~(comp_r_3_107) + (sc_uint<1>) ~(comp_r_5_107) + (sc_uint<1>) ~(comp_r_7_107) + (sc_uint<1>) ~(comp_r_9_107) + (sc_uint<1>) ~(comp_r_11_107) + (sc_uint<1>) ~(comp_r_13_107) + (sc_uint<1>) ~(comp_r_15_107) + (sc_uint<1>) ~(comp_r_17_107) + (sc_uint<1>) ~(comp_r_19_107) + (sc_uint<1>) ~(comp_r_21_107) + (sc_uint<1>) ~(comp_r_23_107) + (sc_uint<1>) ~(comp_r_25_107) + (sc_uint<1>) ~(comp_r_27_107) + (sc_uint<1>) ~(comp_r_29_107) + (sc_uint<1>) ~(comp_r_31_107) + (sc_uint<1>) ~(comp_r_33_107) + (sc_uint<1>) ~(comp_r_35_107) + (sc_uint<1>) ~(comp_r_37_107) + (sc_uint<1>) ~(comp_r_39_107) + (sc_uint<1>) ~(comp_r_41_107) + (sc_uint<1>) ~(comp_r_43_107) + (sc_uint<1>) ~(comp_r_45_107) + (sc_uint<1>) ~(comp_r_47_107) + (sc_uint<1>) ~(comp_r_49_107) + (sc_uint<1>) ~(comp_r_51_107) + (sc_uint<1>) ~(comp_r_53_107) + (sc_uint<1>) ~(comp_r_55_107) + (sc_uint<1>) ~(comp_r_57_107) + (sc_uint<1>) ~(comp_r_59_107) + (sc_uint<1>) ~(comp_r_61_107) + (sc_uint<1>) ~(comp_r_63_107) + (sc_uint<1>) ~(comp_r_65_107) + (sc_uint<1>) ~(comp_r_67_107) + (sc_uint<1>) ~(comp_r_69_107) + (sc_uint<1>) ~(comp_r_71_107) + (sc_uint<1>) ~(comp_r_73_107) + (sc_uint<1>) ~(comp_r_75_107) + (sc_uint<1>) ~(comp_r_77_107) + (sc_uint<1>) ~(comp_r_79_107) + (sc_uint<1>) ~(comp_r_81_107) + (sc_uint<1>) ~(comp_r_83_107) + (sc_uint<1>) ~(comp_r_85_107) + (sc_uint<1>) ~(comp_r_87_107) + (sc_uint<1>) ~(comp_r_89_107) + (sc_uint<1>) ~(comp_r_91_107) + (sc_uint<1>) ~(comp_r_93_107) + (sc_uint<1>) ~(comp_r_95_107) + (sc_uint<1>) ~(comp_r_97_107) + (sc_uint<1>) ~(comp_r_99_107) + (sc_uint<1>) ~(comp_r_101_107) + (sc_uint<1>) ~(comp_r_103_107) + (sc_uint<1>) ~(comp_r_105_107) + comp_r_107_108 + comp_r_107_109 + comp_r_107_110 + comp_r_107_111 + comp_r_107_112 + comp_r_107_113 + comp_r_107_114 + comp_r_107_115 + comp_r_107_116 + comp_r_107_117 + comp_r_107_118 + comp_r_107_119 + comp_r_107_120 + comp_r_107_121 + comp_r_107_122 + comp_r_107_123 + comp_r_107_124 + comp_r_107_125 + comp_r_107_126;
    position[108] = 54 + (sc_uint<1>) ~(comp_r_1_108) + (sc_uint<1>) ~(comp_r_3_108) + (sc_uint<1>) ~(comp_r_5_108) + (sc_uint<1>) ~(comp_r_7_108) + (sc_uint<1>) ~(comp_r_9_108) + (sc_uint<1>) ~(comp_r_11_108) + (sc_uint<1>) ~(comp_r_13_108) + (sc_uint<1>) ~(comp_r_15_108) + (sc_uint<1>) ~(comp_r_17_108) + (sc_uint<1>) ~(comp_r_19_108) + (sc_uint<1>) ~(comp_r_21_108) + (sc_uint<1>) ~(comp_r_23_108) + (sc_uint<1>) ~(comp_r_25_108) + (sc_uint<1>) ~(comp_r_27_108) + (sc_uint<1>) ~(comp_r_29_108) + (sc_uint<1>) ~(comp_r_31_108) + (sc_uint<1>) ~(comp_r_33_108) + (sc_uint<1>) ~(comp_r_35_108) + (sc_uint<1>) ~(comp_r_37_108) + (sc_uint<1>) ~(comp_r_39_108) + (sc_uint<1>) ~(comp_r_41_108) + (sc_uint<1>) ~(comp_r_43_108) + (sc_uint<1>) ~(comp_r_45_108) + (sc_uint<1>) ~(comp_r_47_108) + (sc_uint<1>) ~(comp_r_49_108) + (sc_uint<1>) ~(comp_r_51_108) + (sc_uint<1>) ~(comp_r_53_108) + (sc_uint<1>) ~(comp_r_55_108) + (sc_uint<1>) ~(comp_r_57_108) + (sc_uint<1>) ~(comp_r_59_108) + (sc_uint<1>) ~(comp_r_61_108) + (sc_uint<1>) ~(comp_r_63_108) + (sc_uint<1>) ~(comp_r_65_108) + (sc_uint<1>) ~(comp_r_67_108) + (sc_uint<1>) ~(comp_r_69_108) + (sc_uint<1>) ~(comp_r_71_108) + (sc_uint<1>) ~(comp_r_73_108) + (sc_uint<1>) ~(comp_r_75_108) + (sc_uint<1>) ~(comp_r_77_108) + (sc_uint<1>) ~(comp_r_79_108) + (sc_uint<1>) ~(comp_r_81_108) + (sc_uint<1>) ~(comp_r_83_108) + (sc_uint<1>) ~(comp_r_85_108) + (sc_uint<1>) ~(comp_r_87_108) + (sc_uint<1>) ~(comp_r_89_108) + (sc_uint<1>) ~(comp_r_91_108) + (sc_uint<1>) ~(comp_r_93_108) + (sc_uint<1>) ~(comp_r_95_108) + (sc_uint<1>) ~(comp_r_97_108) + (sc_uint<1>) ~(comp_r_99_108) + (sc_uint<1>) ~(comp_r_101_108) + (sc_uint<1>) ~(comp_r_103_108) + (sc_uint<1>) ~(comp_r_105_108) + (sc_uint<1>) ~(comp_r_107_108);
    position[109] = 55 + (sc_uint<1>) ~(comp_r_1_109) + (sc_uint<1>) ~(comp_r_3_109) + (sc_uint<1>) ~(comp_r_5_109) + (sc_uint<1>) ~(comp_r_7_109) + (sc_uint<1>) ~(comp_r_9_109) + (sc_uint<1>) ~(comp_r_11_109) + (sc_uint<1>) ~(comp_r_13_109) + (sc_uint<1>) ~(comp_r_15_109) + (sc_uint<1>) ~(comp_r_17_109) + (sc_uint<1>) ~(comp_r_19_109) + (sc_uint<1>) ~(comp_r_21_109) + (sc_uint<1>) ~(comp_r_23_109) + (sc_uint<1>) ~(comp_r_25_109) + (sc_uint<1>) ~(comp_r_27_109) + (sc_uint<1>) ~(comp_r_29_109) + (sc_uint<1>) ~(comp_r_31_109) + (sc_uint<1>) ~(comp_r_33_109) + (sc_uint<1>) ~(comp_r_35_109) + (sc_uint<1>) ~(comp_r_37_109) + (sc_uint<1>) ~(comp_r_39_109) + (sc_uint<1>) ~(comp_r_41_109) + (sc_uint<1>) ~(comp_r_43_109) + (sc_uint<1>) ~(comp_r_45_109) + (sc_uint<1>) ~(comp_r_47_109) + (sc_uint<1>) ~(comp_r_49_109) + (sc_uint<1>) ~(comp_r_51_109) + (sc_uint<1>) ~(comp_r_53_109) + (sc_uint<1>) ~(comp_r_55_109) + (sc_uint<1>) ~(comp_r_57_109) + (sc_uint<1>) ~(comp_r_59_109) + (sc_uint<1>) ~(comp_r_61_109) + (sc_uint<1>) ~(comp_r_63_109) + (sc_uint<1>) ~(comp_r_65_109) + (sc_uint<1>) ~(comp_r_67_109) + (sc_uint<1>) ~(comp_r_69_109) + (sc_uint<1>) ~(comp_r_71_109) + (sc_uint<1>) ~(comp_r_73_109) + (sc_uint<1>) ~(comp_r_75_109) + (sc_uint<1>) ~(comp_r_77_109) + (sc_uint<1>) ~(comp_r_79_109) + (sc_uint<1>) ~(comp_r_81_109) + (sc_uint<1>) ~(comp_r_83_109) + (sc_uint<1>) ~(comp_r_85_109) + (sc_uint<1>) ~(comp_r_87_109) + (sc_uint<1>) ~(comp_r_89_109) + (sc_uint<1>) ~(comp_r_91_109) + (sc_uint<1>) ~(comp_r_93_109) + (sc_uint<1>) ~(comp_r_95_109) + (sc_uint<1>) ~(comp_r_97_109) + (sc_uint<1>) ~(comp_r_99_109) + (sc_uint<1>) ~(comp_r_101_109) + (sc_uint<1>) ~(comp_r_103_109) + (sc_uint<1>) ~(comp_r_105_109) + (sc_uint<1>) ~(comp_r_107_109) + comp_r_109_110 + comp_r_109_111 + comp_r_109_112 + comp_r_109_113 + comp_r_109_114 + comp_r_109_115 + comp_r_109_116 + comp_r_109_117 + comp_r_109_118 + comp_r_109_119 + comp_r_109_120 + comp_r_109_121 + comp_r_109_122 + comp_r_109_123 + comp_r_109_124 + comp_r_109_125 + comp_r_109_126;
    position[110] = 55 + (sc_uint<1>) ~(comp_r_1_110) + (sc_uint<1>) ~(comp_r_3_110) + (sc_uint<1>) ~(comp_r_5_110) + (sc_uint<1>) ~(comp_r_7_110) + (sc_uint<1>) ~(comp_r_9_110) + (sc_uint<1>) ~(comp_r_11_110) + (sc_uint<1>) ~(comp_r_13_110) + (sc_uint<1>) ~(comp_r_15_110) + (sc_uint<1>) ~(comp_r_17_110) + (sc_uint<1>) ~(comp_r_19_110) + (sc_uint<1>) ~(comp_r_21_110) + (sc_uint<1>) ~(comp_r_23_110) + (sc_uint<1>) ~(comp_r_25_110) + (sc_uint<1>) ~(comp_r_27_110) + (sc_uint<1>) ~(comp_r_29_110) + (sc_uint<1>) ~(comp_r_31_110) + (sc_uint<1>) ~(comp_r_33_110) + (sc_uint<1>) ~(comp_r_35_110) + (sc_uint<1>) ~(comp_r_37_110) + (sc_uint<1>) ~(comp_r_39_110) + (sc_uint<1>) ~(comp_r_41_110) + (sc_uint<1>) ~(comp_r_43_110) + (sc_uint<1>) ~(comp_r_45_110) + (sc_uint<1>) ~(comp_r_47_110) + (sc_uint<1>) ~(comp_r_49_110) + (sc_uint<1>) ~(comp_r_51_110) + (sc_uint<1>) ~(comp_r_53_110) + (sc_uint<1>) ~(comp_r_55_110) + (sc_uint<1>) ~(comp_r_57_110) + (sc_uint<1>) ~(comp_r_59_110) + (sc_uint<1>) ~(comp_r_61_110) + (sc_uint<1>) ~(comp_r_63_110) + (sc_uint<1>) ~(comp_r_65_110) + (sc_uint<1>) ~(comp_r_67_110) + (sc_uint<1>) ~(comp_r_69_110) + (sc_uint<1>) ~(comp_r_71_110) + (sc_uint<1>) ~(comp_r_73_110) + (sc_uint<1>) ~(comp_r_75_110) + (sc_uint<1>) ~(comp_r_77_110) + (sc_uint<1>) ~(comp_r_79_110) + (sc_uint<1>) ~(comp_r_81_110) + (sc_uint<1>) ~(comp_r_83_110) + (sc_uint<1>) ~(comp_r_85_110) + (sc_uint<1>) ~(comp_r_87_110) + (sc_uint<1>) ~(comp_r_89_110) + (sc_uint<1>) ~(comp_r_91_110) + (sc_uint<1>) ~(comp_r_93_110) + (sc_uint<1>) ~(comp_r_95_110) + (sc_uint<1>) ~(comp_r_97_110) + (sc_uint<1>) ~(comp_r_99_110) + (sc_uint<1>) ~(comp_r_101_110) + (sc_uint<1>) ~(comp_r_103_110) + (sc_uint<1>) ~(comp_r_105_110) + (sc_uint<1>) ~(comp_r_107_110) + (sc_uint<1>) ~(comp_r_109_110);
    position[111] = 56 + (sc_uint<1>) ~(comp_r_1_111) + (sc_uint<1>) ~(comp_r_3_111) + (sc_uint<1>) ~(comp_r_5_111) + (sc_uint<1>) ~(comp_r_7_111) + (sc_uint<1>) ~(comp_r_9_111) + (sc_uint<1>) ~(comp_r_11_111) + (sc_uint<1>) ~(comp_r_13_111) + (sc_uint<1>) ~(comp_r_15_111) + (sc_uint<1>) ~(comp_r_17_111) + (sc_uint<1>) ~(comp_r_19_111) + (sc_uint<1>) ~(comp_r_21_111) + (sc_uint<1>) ~(comp_r_23_111) + (sc_uint<1>) ~(comp_r_25_111) + (sc_uint<1>) ~(comp_r_27_111) + (sc_uint<1>) ~(comp_r_29_111) + (sc_uint<1>) ~(comp_r_31_111) + (sc_uint<1>) ~(comp_r_33_111) + (sc_uint<1>) ~(comp_r_35_111) + (sc_uint<1>) ~(comp_r_37_111) + (sc_uint<1>) ~(comp_r_39_111) + (sc_uint<1>) ~(comp_r_41_111) + (sc_uint<1>) ~(comp_r_43_111) + (sc_uint<1>) ~(comp_r_45_111) + (sc_uint<1>) ~(comp_r_47_111) + (sc_uint<1>) ~(comp_r_49_111) + (sc_uint<1>) ~(comp_r_51_111) + (sc_uint<1>) ~(comp_r_53_111) + (sc_uint<1>) ~(comp_r_55_111) + (sc_uint<1>) ~(comp_r_57_111) + (sc_uint<1>) ~(comp_r_59_111) + (sc_uint<1>) ~(comp_r_61_111) + (sc_uint<1>) ~(comp_r_63_111) + (sc_uint<1>) ~(comp_r_65_111) + (sc_uint<1>) ~(comp_r_67_111) + (sc_uint<1>) ~(comp_r_69_111) + (sc_uint<1>) ~(comp_r_71_111) + (sc_uint<1>) ~(comp_r_73_111) + (sc_uint<1>) ~(comp_r_75_111) + (sc_uint<1>) ~(comp_r_77_111) + (sc_uint<1>) ~(comp_r_79_111) + (sc_uint<1>) ~(comp_r_81_111) + (sc_uint<1>) ~(comp_r_83_111) + (sc_uint<1>) ~(comp_r_85_111) + (sc_uint<1>) ~(comp_r_87_111) + (sc_uint<1>) ~(comp_r_89_111) + (sc_uint<1>) ~(comp_r_91_111) + (sc_uint<1>) ~(comp_r_93_111) + (sc_uint<1>) ~(comp_r_95_111) + (sc_uint<1>) ~(comp_r_97_111) + (sc_uint<1>) ~(comp_r_99_111) + (sc_uint<1>) ~(comp_r_101_111) + (sc_uint<1>) ~(comp_r_103_111) + (sc_uint<1>) ~(comp_r_105_111) + (sc_uint<1>) ~(comp_r_107_111) + (sc_uint<1>) ~(comp_r_109_111) + comp_r_111_112 + comp_r_111_113 + comp_r_111_114 + comp_r_111_115 + comp_r_111_116 + comp_r_111_117 + comp_r_111_118 + comp_r_111_119 + comp_r_111_120 + comp_r_111_121 + comp_r_111_122 + comp_r_111_123 + comp_r_111_124 + comp_r_111_125 + comp_r_111_126;
    position[112] = 56 + (sc_uint<1>) ~(comp_r_1_112) + (sc_uint<1>) ~(comp_r_3_112) + (sc_uint<1>) ~(comp_r_5_112) + (sc_uint<1>) ~(comp_r_7_112) + (sc_uint<1>) ~(comp_r_9_112) + (sc_uint<1>) ~(comp_r_11_112) + (sc_uint<1>) ~(comp_r_13_112) + (sc_uint<1>) ~(comp_r_15_112) + (sc_uint<1>) ~(comp_r_17_112) + (sc_uint<1>) ~(comp_r_19_112) + (sc_uint<1>) ~(comp_r_21_112) + (sc_uint<1>) ~(comp_r_23_112) + (sc_uint<1>) ~(comp_r_25_112) + (sc_uint<1>) ~(comp_r_27_112) + (sc_uint<1>) ~(comp_r_29_112) + (sc_uint<1>) ~(comp_r_31_112) + (sc_uint<1>) ~(comp_r_33_112) + (sc_uint<1>) ~(comp_r_35_112) + (sc_uint<1>) ~(comp_r_37_112) + (sc_uint<1>) ~(comp_r_39_112) + (sc_uint<1>) ~(comp_r_41_112) + (sc_uint<1>) ~(comp_r_43_112) + (sc_uint<1>) ~(comp_r_45_112) + (sc_uint<1>) ~(comp_r_47_112) + (sc_uint<1>) ~(comp_r_49_112) + (sc_uint<1>) ~(comp_r_51_112) + (sc_uint<1>) ~(comp_r_53_112) + (sc_uint<1>) ~(comp_r_55_112) + (sc_uint<1>) ~(comp_r_57_112) + (sc_uint<1>) ~(comp_r_59_112) + (sc_uint<1>) ~(comp_r_61_112) + (sc_uint<1>) ~(comp_r_63_112) + (sc_uint<1>) ~(comp_r_65_112) + (sc_uint<1>) ~(comp_r_67_112) + (sc_uint<1>) ~(comp_r_69_112) + (sc_uint<1>) ~(comp_r_71_112) + (sc_uint<1>) ~(comp_r_73_112) + (sc_uint<1>) ~(comp_r_75_112) + (sc_uint<1>) ~(comp_r_77_112) + (sc_uint<1>) ~(comp_r_79_112) + (sc_uint<1>) ~(comp_r_81_112) + (sc_uint<1>) ~(comp_r_83_112) + (sc_uint<1>) ~(comp_r_85_112) + (sc_uint<1>) ~(comp_r_87_112) + (sc_uint<1>) ~(comp_r_89_112) + (sc_uint<1>) ~(comp_r_91_112) + (sc_uint<1>) ~(comp_r_93_112) + (sc_uint<1>) ~(comp_r_95_112) + (sc_uint<1>) ~(comp_r_97_112) + (sc_uint<1>) ~(comp_r_99_112) + (sc_uint<1>) ~(comp_r_101_112) + (sc_uint<1>) ~(comp_r_103_112) + (sc_uint<1>) ~(comp_r_105_112) + (sc_uint<1>) ~(comp_r_107_112) + (sc_uint<1>) ~(comp_r_109_112) + (sc_uint<1>) ~(comp_r_111_112);
    position[113] = 57 + (sc_uint<1>) ~(comp_r_1_113) + (sc_uint<1>) ~(comp_r_3_113) + (sc_uint<1>) ~(comp_r_5_113) + (sc_uint<1>) ~(comp_r_7_113) + (sc_uint<1>) ~(comp_r_9_113) + (sc_uint<1>) ~(comp_r_11_113) + (sc_uint<1>) ~(comp_r_13_113) + (sc_uint<1>) ~(comp_r_15_113) + (sc_uint<1>) ~(comp_r_17_113) + (sc_uint<1>) ~(comp_r_19_113) + (sc_uint<1>) ~(comp_r_21_113) + (sc_uint<1>) ~(comp_r_23_113) + (sc_uint<1>) ~(comp_r_25_113) + (sc_uint<1>) ~(comp_r_27_113) + (sc_uint<1>) ~(comp_r_29_113) + (sc_uint<1>) ~(comp_r_31_113) + (sc_uint<1>) ~(comp_r_33_113) + (sc_uint<1>) ~(comp_r_35_113) + (sc_uint<1>) ~(comp_r_37_113) + (sc_uint<1>) ~(comp_r_39_113) + (sc_uint<1>) ~(comp_r_41_113) + (sc_uint<1>) ~(comp_r_43_113) + (sc_uint<1>) ~(comp_r_45_113) + (sc_uint<1>) ~(comp_r_47_113) + (sc_uint<1>) ~(comp_r_49_113) + (sc_uint<1>) ~(comp_r_51_113) + (sc_uint<1>) ~(comp_r_53_113) + (sc_uint<1>) ~(comp_r_55_113) + (sc_uint<1>) ~(comp_r_57_113) + (sc_uint<1>) ~(comp_r_59_113) + (sc_uint<1>) ~(comp_r_61_113) + (sc_uint<1>) ~(comp_r_63_113) + (sc_uint<1>) ~(comp_r_65_113) + (sc_uint<1>) ~(comp_r_67_113) + (sc_uint<1>) ~(comp_r_69_113) + (sc_uint<1>) ~(comp_r_71_113) + (sc_uint<1>) ~(comp_r_73_113) + (sc_uint<1>) ~(comp_r_75_113) + (sc_uint<1>) ~(comp_r_77_113) + (sc_uint<1>) ~(comp_r_79_113) + (sc_uint<1>) ~(comp_r_81_113) + (sc_uint<1>) ~(comp_r_83_113) + (sc_uint<1>) ~(comp_r_85_113) + (sc_uint<1>) ~(comp_r_87_113) + (sc_uint<1>) ~(comp_r_89_113) + (sc_uint<1>) ~(comp_r_91_113) + (sc_uint<1>) ~(comp_r_93_113) + (sc_uint<1>) ~(comp_r_95_113) + (sc_uint<1>) ~(comp_r_97_113) + (sc_uint<1>) ~(comp_r_99_113) + (sc_uint<1>) ~(comp_r_101_113) + (sc_uint<1>) ~(comp_r_103_113) + (sc_uint<1>) ~(comp_r_105_113) + (sc_uint<1>) ~(comp_r_107_113) + (sc_uint<1>) ~(comp_r_109_113) + (sc_uint<1>) ~(comp_r_111_113) + comp_r_113_114 + comp_r_113_115 + comp_r_113_116 + comp_r_113_117 + comp_r_113_118 + comp_r_113_119 + comp_r_113_120 + comp_r_113_121 + comp_r_113_122 + comp_r_113_123 + comp_r_113_124 + comp_r_113_125 + comp_r_113_126;
    position[114] = 57 + (sc_uint<1>) ~(comp_r_1_114) + (sc_uint<1>) ~(comp_r_3_114) + (sc_uint<1>) ~(comp_r_5_114) + (sc_uint<1>) ~(comp_r_7_114) + (sc_uint<1>) ~(comp_r_9_114) + (sc_uint<1>) ~(comp_r_11_114) + (sc_uint<1>) ~(comp_r_13_114) + (sc_uint<1>) ~(comp_r_15_114) + (sc_uint<1>) ~(comp_r_17_114) + (sc_uint<1>) ~(comp_r_19_114) + (sc_uint<1>) ~(comp_r_21_114) + (sc_uint<1>) ~(comp_r_23_114) + (sc_uint<1>) ~(comp_r_25_114) + (sc_uint<1>) ~(comp_r_27_114) + (sc_uint<1>) ~(comp_r_29_114) + (sc_uint<1>) ~(comp_r_31_114) + (sc_uint<1>) ~(comp_r_33_114) + (sc_uint<1>) ~(comp_r_35_114) + (sc_uint<1>) ~(comp_r_37_114) + (sc_uint<1>) ~(comp_r_39_114) + (sc_uint<1>) ~(comp_r_41_114) + (sc_uint<1>) ~(comp_r_43_114) + (sc_uint<1>) ~(comp_r_45_114) + (sc_uint<1>) ~(comp_r_47_114) + (sc_uint<1>) ~(comp_r_49_114) + (sc_uint<1>) ~(comp_r_51_114) + (sc_uint<1>) ~(comp_r_53_114) + (sc_uint<1>) ~(comp_r_55_114) + (sc_uint<1>) ~(comp_r_57_114) + (sc_uint<1>) ~(comp_r_59_114) + (sc_uint<1>) ~(comp_r_61_114) + (sc_uint<1>) ~(comp_r_63_114) + (sc_uint<1>) ~(comp_r_65_114) + (sc_uint<1>) ~(comp_r_67_114) + (sc_uint<1>) ~(comp_r_69_114) + (sc_uint<1>) ~(comp_r_71_114) + (sc_uint<1>) ~(comp_r_73_114) + (sc_uint<1>) ~(comp_r_75_114) + (sc_uint<1>) ~(comp_r_77_114) + (sc_uint<1>) ~(comp_r_79_114) + (sc_uint<1>) ~(comp_r_81_114) + (sc_uint<1>) ~(comp_r_83_114) + (sc_uint<1>) ~(comp_r_85_114) + (sc_uint<1>) ~(comp_r_87_114) + (sc_uint<1>) ~(comp_r_89_114) + (sc_uint<1>) ~(comp_r_91_114) + (sc_uint<1>) ~(comp_r_93_114) + (sc_uint<1>) ~(comp_r_95_114) + (sc_uint<1>) ~(comp_r_97_114) + (sc_uint<1>) ~(comp_r_99_114) + (sc_uint<1>) ~(comp_r_101_114) + (sc_uint<1>) ~(comp_r_103_114) + (sc_uint<1>) ~(comp_r_105_114) + (sc_uint<1>) ~(comp_r_107_114) + (sc_uint<1>) ~(comp_r_109_114) + (sc_uint<1>) ~(comp_r_111_114) + (sc_uint<1>) ~(comp_r_113_114);
    position[115] = 58 + (sc_uint<1>) ~(comp_r_1_115) + (sc_uint<1>) ~(comp_r_3_115) + (sc_uint<1>) ~(comp_r_5_115) + (sc_uint<1>) ~(comp_r_7_115) + (sc_uint<1>) ~(comp_r_9_115) + (sc_uint<1>) ~(comp_r_11_115) + (sc_uint<1>) ~(comp_r_13_115) + (sc_uint<1>) ~(comp_r_15_115) + (sc_uint<1>) ~(comp_r_17_115) + (sc_uint<1>) ~(comp_r_19_115) + (sc_uint<1>) ~(comp_r_21_115) + (sc_uint<1>) ~(comp_r_23_115) + (sc_uint<1>) ~(comp_r_25_115) + (sc_uint<1>) ~(comp_r_27_115) + (sc_uint<1>) ~(comp_r_29_115) + (sc_uint<1>) ~(comp_r_31_115) + (sc_uint<1>) ~(comp_r_33_115) + (sc_uint<1>) ~(comp_r_35_115) + (sc_uint<1>) ~(comp_r_37_115) + (sc_uint<1>) ~(comp_r_39_115) + (sc_uint<1>) ~(comp_r_41_115) + (sc_uint<1>) ~(comp_r_43_115) + (sc_uint<1>) ~(comp_r_45_115) + (sc_uint<1>) ~(comp_r_47_115) + (sc_uint<1>) ~(comp_r_49_115) + (sc_uint<1>) ~(comp_r_51_115) + (sc_uint<1>) ~(comp_r_53_115) + (sc_uint<1>) ~(comp_r_55_115) + (sc_uint<1>) ~(comp_r_57_115) + (sc_uint<1>) ~(comp_r_59_115) + (sc_uint<1>) ~(comp_r_61_115) + (sc_uint<1>) ~(comp_r_63_115) + (sc_uint<1>) ~(comp_r_65_115) + (sc_uint<1>) ~(comp_r_67_115) + (sc_uint<1>) ~(comp_r_69_115) + (sc_uint<1>) ~(comp_r_71_115) + (sc_uint<1>) ~(comp_r_73_115) + (sc_uint<1>) ~(comp_r_75_115) + (sc_uint<1>) ~(comp_r_77_115) + (sc_uint<1>) ~(comp_r_79_115) + (sc_uint<1>) ~(comp_r_81_115) + (sc_uint<1>) ~(comp_r_83_115) + (sc_uint<1>) ~(comp_r_85_115) + (sc_uint<1>) ~(comp_r_87_115) + (sc_uint<1>) ~(comp_r_89_115) + (sc_uint<1>) ~(comp_r_91_115) + (sc_uint<1>) ~(comp_r_93_115) + (sc_uint<1>) ~(comp_r_95_115) + (sc_uint<1>) ~(comp_r_97_115) + (sc_uint<1>) ~(comp_r_99_115) + (sc_uint<1>) ~(comp_r_101_115) + (sc_uint<1>) ~(comp_r_103_115) + (sc_uint<1>) ~(comp_r_105_115) + (sc_uint<1>) ~(comp_r_107_115) + (sc_uint<1>) ~(comp_r_109_115) + (sc_uint<1>) ~(comp_r_111_115) + (sc_uint<1>) ~(comp_r_113_115) + comp_r_115_116 + comp_r_115_117 + comp_r_115_118 + comp_r_115_119 + comp_r_115_120 + comp_r_115_121 + comp_r_115_122 + comp_r_115_123 + comp_r_115_124 + comp_r_115_125 + comp_r_115_126;
    position[116] = 58 + (sc_uint<1>) ~(comp_r_1_116) + (sc_uint<1>) ~(comp_r_3_116) + (sc_uint<1>) ~(comp_r_5_116) + (sc_uint<1>) ~(comp_r_7_116) + (sc_uint<1>) ~(comp_r_9_116) + (sc_uint<1>) ~(comp_r_11_116) + (sc_uint<1>) ~(comp_r_13_116) + (sc_uint<1>) ~(comp_r_15_116) + (sc_uint<1>) ~(comp_r_17_116) + (sc_uint<1>) ~(comp_r_19_116) + (sc_uint<1>) ~(comp_r_21_116) + (sc_uint<1>) ~(comp_r_23_116) + (sc_uint<1>) ~(comp_r_25_116) + (sc_uint<1>) ~(comp_r_27_116) + (sc_uint<1>) ~(comp_r_29_116) + (sc_uint<1>) ~(comp_r_31_116) + (sc_uint<1>) ~(comp_r_33_116) + (sc_uint<1>) ~(comp_r_35_116) + (sc_uint<1>) ~(comp_r_37_116) + (sc_uint<1>) ~(comp_r_39_116) + (sc_uint<1>) ~(comp_r_41_116) + (sc_uint<1>) ~(comp_r_43_116) + (sc_uint<1>) ~(comp_r_45_116) + (sc_uint<1>) ~(comp_r_47_116) + (sc_uint<1>) ~(comp_r_49_116) + (sc_uint<1>) ~(comp_r_51_116) + (sc_uint<1>) ~(comp_r_53_116) + (sc_uint<1>) ~(comp_r_55_116) + (sc_uint<1>) ~(comp_r_57_116) + (sc_uint<1>) ~(comp_r_59_116) + (sc_uint<1>) ~(comp_r_61_116) + (sc_uint<1>) ~(comp_r_63_116) + (sc_uint<1>) ~(comp_r_65_116) + (sc_uint<1>) ~(comp_r_67_116) + (sc_uint<1>) ~(comp_r_69_116) + (sc_uint<1>) ~(comp_r_71_116) + (sc_uint<1>) ~(comp_r_73_116) + (sc_uint<1>) ~(comp_r_75_116) + (sc_uint<1>) ~(comp_r_77_116) + (sc_uint<1>) ~(comp_r_79_116) + (sc_uint<1>) ~(comp_r_81_116) + (sc_uint<1>) ~(comp_r_83_116) + (sc_uint<1>) ~(comp_r_85_116) + (sc_uint<1>) ~(comp_r_87_116) + (sc_uint<1>) ~(comp_r_89_116) + (sc_uint<1>) ~(comp_r_91_116) + (sc_uint<1>) ~(comp_r_93_116) + (sc_uint<1>) ~(comp_r_95_116) + (sc_uint<1>) ~(comp_r_97_116) + (sc_uint<1>) ~(comp_r_99_116) + (sc_uint<1>) ~(comp_r_101_116) + (sc_uint<1>) ~(comp_r_103_116) + (sc_uint<1>) ~(comp_r_105_116) + (sc_uint<1>) ~(comp_r_107_116) + (sc_uint<1>) ~(comp_r_109_116) + (sc_uint<1>) ~(comp_r_111_116) + (sc_uint<1>) ~(comp_r_113_116) + (sc_uint<1>) ~(comp_r_115_116);
    position[117] = 59 + (sc_uint<1>) ~(comp_r_1_117) + (sc_uint<1>) ~(comp_r_3_117) + (sc_uint<1>) ~(comp_r_5_117) + (sc_uint<1>) ~(comp_r_7_117) + (sc_uint<1>) ~(comp_r_9_117) + (sc_uint<1>) ~(comp_r_11_117) + (sc_uint<1>) ~(comp_r_13_117) + (sc_uint<1>) ~(comp_r_15_117) + (sc_uint<1>) ~(comp_r_17_117) + (sc_uint<1>) ~(comp_r_19_117) + (sc_uint<1>) ~(comp_r_21_117) + (sc_uint<1>) ~(comp_r_23_117) + (sc_uint<1>) ~(comp_r_25_117) + (sc_uint<1>) ~(comp_r_27_117) + (sc_uint<1>) ~(comp_r_29_117) + (sc_uint<1>) ~(comp_r_31_117) + (sc_uint<1>) ~(comp_r_33_117) + (sc_uint<1>) ~(comp_r_35_117) + (sc_uint<1>) ~(comp_r_37_117) + (sc_uint<1>) ~(comp_r_39_117) + (sc_uint<1>) ~(comp_r_41_117) + (sc_uint<1>) ~(comp_r_43_117) + (sc_uint<1>) ~(comp_r_45_117) + (sc_uint<1>) ~(comp_r_47_117) + (sc_uint<1>) ~(comp_r_49_117) + (sc_uint<1>) ~(comp_r_51_117) + (sc_uint<1>) ~(comp_r_53_117) + (sc_uint<1>) ~(comp_r_55_117) + (sc_uint<1>) ~(comp_r_57_117) + (sc_uint<1>) ~(comp_r_59_117) + (sc_uint<1>) ~(comp_r_61_117) + (sc_uint<1>) ~(comp_r_63_117) + (sc_uint<1>) ~(comp_r_65_117) + (sc_uint<1>) ~(comp_r_67_117) + (sc_uint<1>) ~(comp_r_69_117) + (sc_uint<1>) ~(comp_r_71_117) + (sc_uint<1>) ~(comp_r_73_117) + (sc_uint<1>) ~(comp_r_75_117) + (sc_uint<1>) ~(comp_r_77_117) + (sc_uint<1>) ~(comp_r_79_117) + (sc_uint<1>) ~(comp_r_81_117) + (sc_uint<1>) ~(comp_r_83_117) + (sc_uint<1>) ~(comp_r_85_117) + (sc_uint<1>) ~(comp_r_87_117) + (sc_uint<1>) ~(comp_r_89_117) + (sc_uint<1>) ~(comp_r_91_117) + (sc_uint<1>) ~(comp_r_93_117) + (sc_uint<1>) ~(comp_r_95_117) + (sc_uint<1>) ~(comp_r_97_117) + (sc_uint<1>) ~(comp_r_99_117) + (sc_uint<1>) ~(comp_r_101_117) + (sc_uint<1>) ~(comp_r_103_117) + (sc_uint<1>) ~(comp_r_105_117) + (sc_uint<1>) ~(comp_r_107_117) + (sc_uint<1>) ~(comp_r_109_117) + (sc_uint<1>) ~(comp_r_111_117) + (sc_uint<1>) ~(comp_r_113_117) + (sc_uint<1>) ~(comp_r_115_117) + comp_r_117_118 + comp_r_117_119 + comp_r_117_120 + comp_r_117_121 + comp_r_117_122 + comp_r_117_123 + comp_r_117_124 + comp_r_117_125 + comp_r_117_126;
    position[118] = 59 + (sc_uint<1>) ~(comp_r_1_118) + (sc_uint<1>) ~(comp_r_3_118) + (sc_uint<1>) ~(comp_r_5_118) + (sc_uint<1>) ~(comp_r_7_118) + (sc_uint<1>) ~(comp_r_9_118) + (sc_uint<1>) ~(comp_r_11_118) + (sc_uint<1>) ~(comp_r_13_118) + (sc_uint<1>) ~(comp_r_15_118) + (sc_uint<1>) ~(comp_r_17_118) + (sc_uint<1>) ~(comp_r_19_118) + (sc_uint<1>) ~(comp_r_21_118) + (sc_uint<1>) ~(comp_r_23_118) + (sc_uint<1>) ~(comp_r_25_118) + (sc_uint<1>) ~(comp_r_27_118) + (sc_uint<1>) ~(comp_r_29_118) + (sc_uint<1>) ~(comp_r_31_118) + (sc_uint<1>) ~(comp_r_33_118) + (sc_uint<1>) ~(comp_r_35_118) + (sc_uint<1>) ~(comp_r_37_118) + (sc_uint<1>) ~(comp_r_39_118) + (sc_uint<1>) ~(comp_r_41_118) + (sc_uint<1>) ~(comp_r_43_118) + (sc_uint<1>) ~(comp_r_45_118) + (sc_uint<1>) ~(comp_r_47_118) + (sc_uint<1>) ~(comp_r_49_118) + (sc_uint<1>) ~(comp_r_51_118) + (sc_uint<1>) ~(comp_r_53_118) + (sc_uint<1>) ~(comp_r_55_118) + (sc_uint<1>) ~(comp_r_57_118) + (sc_uint<1>) ~(comp_r_59_118) + (sc_uint<1>) ~(comp_r_61_118) + (sc_uint<1>) ~(comp_r_63_118) + (sc_uint<1>) ~(comp_r_65_118) + (sc_uint<1>) ~(comp_r_67_118) + (sc_uint<1>) ~(comp_r_69_118) + (sc_uint<1>) ~(comp_r_71_118) + (sc_uint<1>) ~(comp_r_73_118) + (sc_uint<1>) ~(comp_r_75_118) + (sc_uint<1>) ~(comp_r_77_118) + (sc_uint<1>) ~(comp_r_79_118) + (sc_uint<1>) ~(comp_r_81_118) + (sc_uint<1>) ~(comp_r_83_118) + (sc_uint<1>) ~(comp_r_85_118) + (sc_uint<1>) ~(comp_r_87_118) + (sc_uint<1>) ~(comp_r_89_118) + (sc_uint<1>) ~(comp_r_91_118) + (sc_uint<1>) ~(comp_r_93_118) + (sc_uint<1>) ~(comp_r_95_118) + (sc_uint<1>) ~(comp_r_97_118) + (sc_uint<1>) ~(comp_r_99_118) + (sc_uint<1>) ~(comp_r_101_118) + (sc_uint<1>) ~(comp_r_103_118) + (sc_uint<1>) ~(comp_r_105_118) + (sc_uint<1>) ~(comp_r_107_118) + (sc_uint<1>) ~(comp_r_109_118) + (sc_uint<1>) ~(comp_r_111_118) + (sc_uint<1>) ~(comp_r_113_118) + (sc_uint<1>) ~(comp_r_115_118) + (sc_uint<1>) ~(comp_r_117_118);
    position[119] = 60 + (sc_uint<1>) ~(comp_r_1_119) + (sc_uint<1>) ~(comp_r_3_119) + (sc_uint<1>) ~(comp_r_5_119) + (sc_uint<1>) ~(comp_r_7_119) + (sc_uint<1>) ~(comp_r_9_119) + (sc_uint<1>) ~(comp_r_11_119) + (sc_uint<1>) ~(comp_r_13_119) + (sc_uint<1>) ~(comp_r_15_119) + (sc_uint<1>) ~(comp_r_17_119) + (sc_uint<1>) ~(comp_r_19_119) + (sc_uint<1>) ~(comp_r_21_119) + (sc_uint<1>) ~(comp_r_23_119) + (sc_uint<1>) ~(comp_r_25_119) + (sc_uint<1>) ~(comp_r_27_119) + (sc_uint<1>) ~(comp_r_29_119) + (sc_uint<1>) ~(comp_r_31_119) + (sc_uint<1>) ~(comp_r_33_119) + (sc_uint<1>) ~(comp_r_35_119) + (sc_uint<1>) ~(comp_r_37_119) + (sc_uint<1>) ~(comp_r_39_119) + (sc_uint<1>) ~(comp_r_41_119) + (sc_uint<1>) ~(comp_r_43_119) + (sc_uint<1>) ~(comp_r_45_119) + (sc_uint<1>) ~(comp_r_47_119) + (sc_uint<1>) ~(comp_r_49_119) + (sc_uint<1>) ~(comp_r_51_119) + (sc_uint<1>) ~(comp_r_53_119) + (sc_uint<1>) ~(comp_r_55_119) + (sc_uint<1>) ~(comp_r_57_119) + (sc_uint<1>) ~(comp_r_59_119) + (sc_uint<1>) ~(comp_r_61_119) + (sc_uint<1>) ~(comp_r_63_119) + (sc_uint<1>) ~(comp_r_65_119) + (sc_uint<1>) ~(comp_r_67_119) + (sc_uint<1>) ~(comp_r_69_119) + (sc_uint<1>) ~(comp_r_71_119) + (sc_uint<1>) ~(comp_r_73_119) + (sc_uint<1>) ~(comp_r_75_119) + (sc_uint<1>) ~(comp_r_77_119) + (sc_uint<1>) ~(comp_r_79_119) + (sc_uint<1>) ~(comp_r_81_119) + (sc_uint<1>) ~(comp_r_83_119) + (sc_uint<1>) ~(comp_r_85_119) + (sc_uint<1>) ~(comp_r_87_119) + (sc_uint<1>) ~(comp_r_89_119) + (sc_uint<1>) ~(comp_r_91_119) + (sc_uint<1>) ~(comp_r_93_119) + (sc_uint<1>) ~(comp_r_95_119) + (sc_uint<1>) ~(comp_r_97_119) + (sc_uint<1>) ~(comp_r_99_119) + (sc_uint<1>) ~(comp_r_101_119) + (sc_uint<1>) ~(comp_r_103_119) + (sc_uint<1>) ~(comp_r_105_119) + (sc_uint<1>) ~(comp_r_107_119) + (sc_uint<1>) ~(comp_r_109_119) + (sc_uint<1>) ~(comp_r_111_119) + (sc_uint<1>) ~(comp_r_113_119) + (sc_uint<1>) ~(comp_r_115_119) + (sc_uint<1>) ~(comp_r_117_119) + comp_r_119_120 + comp_r_119_121 + comp_r_119_122 + comp_r_119_123 + comp_r_119_124 + comp_r_119_125 + comp_r_119_126;
    position[120] = 60 + (sc_uint<1>) ~(comp_r_1_120) + (sc_uint<1>) ~(comp_r_3_120) + (sc_uint<1>) ~(comp_r_5_120) + (sc_uint<1>) ~(comp_r_7_120) + (sc_uint<1>) ~(comp_r_9_120) + (sc_uint<1>) ~(comp_r_11_120) + (sc_uint<1>) ~(comp_r_13_120) + (sc_uint<1>) ~(comp_r_15_120) + (sc_uint<1>) ~(comp_r_17_120) + (sc_uint<1>) ~(comp_r_19_120) + (sc_uint<1>) ~(comp_r_21_120) + (sc_uint<1>) ~(comp_r_23_120) + (sc_uint<1>) ~(comp_r_25_120) + (sc_uint<1>) ~(comp_r_27_120) + (sc_uint<1>) ~(comp_r_29_120) + (sc_uint<1>) ~(comp_r_31_120) + (sc_uint<1>) ~(comp_r_33_120) + (sc_uint<1>) ~(comp_r_35_120) + (sc_uint<1>) ~(comp_r_37_120) + (sc_uint<1>) ~(comp_r_39_120) + (sc_uint<1>) ~(comp_r_41_120) + (sc_uint<1>) ~(comp_r_43_120) + (sc_uint<1>) ~(comp_r_45_120) + (sc_uint<1>) ~(comp_r_47_120) + (sc_uint<1>) ~(comp_r_49_120) + (sc_uint<1>) ~(comp_r_51_120) + (sc_uint<1>) ~(comp_r_53_120) + (sc_uint<1>) ~(comp_r_55_120) + (sc_uint<1>) ~(comp_r_57_120) + (sc_uint<1>) ~(comp_r_59_120) + (sc_uint<1>) ~(comp_r_61_120) + (sc_uint<1>) ~(comp_r_63_120) + (sc_uint<1>) ~(comp_r_65_120) + (sc_uint<1>) ~(comp_r_67_120) + (sc_uint<1>) ~(comp_r_69_120) + (sc_uint<1>) ~(comp_r_71_120) + (sc_uint<1>) ~(comp_r_73_120) + (sc_uint<1>) ~(comp_r_75_120) + (sc_uint<1>) ~(comp_r_77_120) + (sc_uint<1>) ~(comp_r_79_120) + (sc_uint<1>) ~(comp_r_81_120) + (sc_uint<1>) ~(comp_r_83_120) + (sc_uint<1>) ~(comp_r_85_120) + (sc_uint<1>) ~(comp_r_87_120) + (sc_uint<1>) ~(comp_r_89_120) + (sc_uint<1>) ~(comp_r_91_120) + (sc_uint<1>) ~(comp_r_93_120) + (sc_uint<1>) ~(comp_r_95_120) + (sc_uint<1>) ~(comp_r_97_120) + (sc_uint<1>) ~(comp_r_99_120) + (sc_uint<1>) ~(comp_r_101_120) + (sc_uint<1>) ~(comp_r_103_120) + (sc_uint<1>) ~(comp_r_105_120) + (sc_uint<1>) ~(comp_r_107_120) + (sc_uint<1>) ~(comp_r_109_120) + (sc_uint<1>) ~(comp_r_111_120) + (sc_uint<1>) ~(comp_r_113_120) + (sc_uint<1>) ~(comp_r_115_120) + (sc_uint<1>) ~(comp_r_117_120) + (sc_uint<1>) ~(comp_r_119_120);
    position[121] = 61 + (sc_uint<1>) ~(comp_r_1_121) + (sc_uint<1>) ~(comp_r_3_121) + (sc_uint<1>) ~(comp_r_5_121) + (sc_uint<1>) ~(comp_r_7_121) + (sc_uint<1>) ~(comp_r_9_121) + (sc_uint<1>) ~(comp_r_11_121) + (sc_uint<1>) ~(comp_r_13_121) + (sc_uint<1>) ~(comp_r_15_121) + (sc_uint<1>) ~(comp_r_17_121) + (sc_uint<1>) ~(comp_r_19_121) + (sc_uint<1>) ~(comp_r_21_121) + (sc_uint<1>) ~(comp_r_23_121) + (sc_uint<1>) ~(comp_r_25_121) + (sc_uint<1>) ~(comp_r_27_121) + (sc_uint<1>) ~(comp_r_29_121) + (sc_uint<1>) ~(comp_r_31_121) + (sc_uint<1>) ~(comp_r_33_121) + (sc_uint<1>) ~(comp_r_35_121) + (sc_uint<1>) ~(comp_r_37_121) + (sc_uint<1>) ~(comp_r_39_121) + (sc_uint<1>) ~(comp_r_41_121) + (sc_uint<1>) ~(comp_r_43_121) + (sc_uint<1>) ~(comp_r_45_121) + (sc_uint<1>) ~(comp_r_47_121) + (sc_uint<1>) ~(comp_r_49_121) + (sc_uint<1>) ~(comp_r_51_121) + (sc_uint<1>) ~(comp_r_53_121) + (sc_uint<1>) ~(comp_r_55_121) + (sc_uint<1>) ~(comp_r_57_121) + (sc_uint<1>) ~(comp_r_59_121) + (sc_uint<1>) ~(comp_r_61_121) + (sc_uint<1>) ~(comp_r_63_121) + (sc_uint<1>) ~(comp_r_65_121) + (sc_uint<1>) ~(comp_r_67_121) + (sc_uint<1>) ~(comp_r_69_121) + (sc_uint<1>) ~(comp_r_71_121) + (sc_uint<1>) ~(comp_r_73_121) + (sc_uint<1>) ~(comp_r_75_121) + (sc_uint<1>) ~(comp_r_77_121) + (sc_uint<1>) ~(comp_r_79_121) + (sc_uint<1>) ~(comp_r_81_121) + (sc_uint<1>) ~(comp_r_83_121) + (sc_uint<1>) ~(comp_r_85_121) + (sc_uint<1>) ~(comp_r_87_121) + (sc_uint<1>) ~(comp_r_89_121) + (sc_uint<1>) ~(comp_r_91_121) + (sc_uint<1>) ~(comp_r_93_121) + (sc_uint<1>) ~(comp_r_95_121) + (sc_uint<1>) ~(comp_r_97_121) + (sc_uint<1>) ~(comp_r_99_121) + (sc_uint<1>) ~(comp_r_101_121) + (sc_uint<1>) ~(comp_r_103_121) + (sc_uint<1>) ~(comp_r_105_121) + (sc_uint<1>) ~(comp_r_107_121) + (sc_uint<1>) ~(comp_r_109_121) + (sc_uint<1>) ~(comp_r_111_121) + (sc_uint<1>) ~(comp_r_113_121) + (sc_uint<1>) ~(comp_r_115_121) + (sc_uint<1>) ~(comp_r_117_121) + (sc_uint<1>) ~(comp_r_119_121) + comp_r_121_122 + comp_r_121_123 + comp_r_121_124 + comp_r_121_125 + comp_r_121_126;
    position[122] = 61 + (sc_uint<1>) ~(comp_r_1_122) + (sc_uint<1>) ~(comp_r_3_122) + (sc_uint<1>) ~(comp_r_5_122) + (sc_uint<1>) ~(comp_r_7_122) + (sc_uint<1>) ~(comp_r_9_122) + (sc_uint<1>) ~(comp_r_11_122) + (sc_uint<1>) ~(comp_r_13_122) + (sc_uint<1>) ~(comp_r_15_122) + (sc_uint<1>) ~(comp_r_17_122) + (sc_uint<1>) ~(comp_r_19_122) + (sc_uint<1>) ~(comp_r_21_122) + (sc_uint<1>) ~(comp_r_23_122) + (sc_uint<1>) ~(comp_r_25_122) + (sc_uint<1>) ~(comp_r_27_122) + (sc_uint<1>) ~(comp_r_29_122) + (sc_uint<1>) ~(comp_r_31_122) + (sc_uint<1>) ~(comp_r_33_122) + (sc_uint<1>) ~(comp_r_35_122) + (sc_uint<1>) ~(comp_r_37_122) + (sc_uint<1>) ~(comp_r_39_122) + (sc_uint<1>) ~(comp_r_41_122) + (sc_uint<1>) ~(comp_r_43_122) + (sc_uint<1>) ~(comp_r_45_122) + (sc_uint<1>) ~(comp_r_47_122) + (sc_uint<1>) ~(comp_r_49_122) + (sc_uint<1>) ~(comp_r_51_122) + (sc_uint<1>) ~(comp_r_53_122) + (sc_uint<1>) ~(comp_r_55_122) + (sc_uint<1>) ~(comp_r_57_122) + (sc_uint<1>) ~(comp_r_59_122) + (sc_uint<1>) ~(comp_r_61_122) + (sc_uint<1>) ~(comp_r_63_122) + (sc_uint<1>) ~(comp_r_65_122) + (sc_uint<1>) ~(comp_r_67_122) + (sc_uint<1>) ~(comp_r_69_122) + (sc_uint<1>) ~(comp_r_71_122) + (sc_uint<1>) ~(comp_r_73_122) + (sc_uint<1>) ~(comp_r_75_122) + (sc_uint<1>) ~(comp_r_77_122) + (sc_uint<1>) ~(comp_r_79_122) + (sc_uint<1>) ~(comp_r_81_122) + (sc_uint<1>) ~(comp_r_83_122) + (sc_uint<1>) ~(comp_r_85_122) + (sc_uint<1>) ~(comp_r_87_122) + (sc_uint<1>) ~(comp_r_89_122) + (sc_uint<1>) ~(comp_r_91_122) + (sc_uint<1>) ~(comp_r_93_122) + (sc_uint<1>) ~(comp_r_95_122) + (sc_uint<1>) ~(comp_r_97_122) + (sc_uint<1>) ~(comp_r_99_122) + (sc_uint<1>) ~(comp_r_101_122) + (sc_uint<1>) ~(comp_r_103_122) + (sc_uint<1>) ~(comp_r_105_122) + (sc_uint<1>) ~(comp_r_107_122) + (sc_uint<1>) ~(comp_r_109_122) + (sc_uint<1>) ~(comp_r_111_122) + (sc_uint<1>) ~(comp_r_113_122) + (sc_uint<1>) ~(comp_r_115_122) + (sc_uint<1>) ~(comp_r_117_122) + (sc_uint<1>) ~(comp_r_119_122) + (sc_uint<1>) ~(comp_r_121_122);
    position[123] = 62 + (sc_uint<1>) ~(comp_r_1_123) + (sc_uint<1>) ~(comp_r_3_123) + (sc_uint<1>) ~(comp_r_5_123) + (sc_uint<1>) ~(comp_r_7_123) + (sc_uint<1>) ~(comp_r_9_123) + (sc_uint<1>) ~(comp_r_11_123) + (sc_uint<1>) ~(comp_r_13_123) + (sc_uint<1>) ~(comp_r_15_123) + (sc_uint<1>) ~(comp_r_17_123) + (sc_uint<1>) ~(comp_r_19_123) + (sc_uint<1>) ~(comp_r_21_123) + (sc_uint<1>) ~(comp_r_23_123) + (sc_uint<1>) ~(comp_r_25_123) + (sc_uint<1>) ~(comp_r_27_123) + (sc_uint<1>) ~(comp_r_29_123) + (sc_uint<1>) ~(comp_r_31_123) + (sc_uint<1>) ~(comp_r_33_123) + (sc_uint<1>) ~(comp_r_35_123) + (sc_uint<1>) ~(comp_r_37_123) + (sc_uint<1>) ~(comp_r_39_123) + (sc_uint<1>) ~(comp_r_41_123) + (sc_uint<1>) ~(comp_r_43_123) + (sc_uint<1>) ~(comp_r_45_123) + (sc_uint<1>) ~(comp_r_47_123) + (sc_uint<1>) ~(comp_r_49_123) + (sc_uint<1>) ~(comp_r_51_123) + (sc_uint<1>) ~(comp_r_53_123) + (sc_uint<1>) ~(comp_r_55_123) + (sc_uint<1>) ~(comp_r_57_123) + (sc_uint<1>) ~(comp_r_59_123) + (sc_uint<1>) ~(comp_r_61_123) + (sc_uint<1>) ~(comp_r_63_123) + (sc_uint<1>) ~(comp_r_65_123) + (sc_uint<1>) ~(comp_r_67_123) + (sc_uint<1>) ~(comp_r_69_123) + (sc_uint<1>) ~(comp_r_71_123) + (sc_uint<1>) ~(comp_r_73_123) + (sc_uint<1>) ~(comp_r_75_123) + (sc_uint<1>) ~(comp_r_77_123) + (sc_uint<1>) ~(comp_r_79_123) + (sc_uint<1>) ~(comp_r_81_123) + (sc_uint<1>) ~(comp_r_83_123) + (sc_uint<1>) ~(comp_r_85_123) + (sc_uint<1>) ~(comp_r_87_123) + (sc_uint<1>) ~(comp_r_89_123) + (sc_uint<1>) ~(comp_r_91_123) + (sc_uint<1>) ~(comp_r_93_123) + (sc_uint<1>) ~(comp_r_95_123) + (sc_uint<1>) ~(comp_r_97_123) + (sc_uint<1>) ~(comp_r_99_123) + (sc_uint<1>) ~(comp_r_101_123) + (sc_uint<1>) ~(comp_r_103_123) + (sc_uint<1>) ~(comp_r_105_123) + (sc_uint<1>) ~(comp_r_107_123) + (sc_uint<1>) ~(comp_r_109_123) + (sc_uint<1>) ~(comp_r_111_123) + (sc_uint<1>) ~(comp_r_113_123) + (sc_uint<1>) ~(comp_r_115_123) + (sc_uint<1>) ~(comp_r_117_123) + (sc_uint<1>) ~(comp_r_119_123) + (sc_uint<1>) ~(comp_r_121_123) + comp_r_123_124 + comp_r_123_125 + comp_r_123_126;
    position[124] = 62 + (sc_uint<1>) ~(comp_r_1_124) + (sc_uint<1>) ~(comp_r_3_124) + (sc_uint<1>) ~(comp_r_5_124) + (sc_uint<1>) ~(comp_r_7_124) + (sc_uint<1>) ~(comp_r_9_124) + (sc_uint<1>) ~(comp_r_11_124) + (sc_uint<1>) ~(comp_r_13_124) + (sc_uint<1>) ~(comp_r_15_124) + (sc_uint<1>) ~(comp_r_17_124) + (sc_uint<1>) ~(comp_r_19_124) + (sc_uint<1>) ~(comp_r_21_124) + (sc_uint<1>) ~(comp_r_23_124) + (sc_uint<1>) ~(comp_r_25_124) + (sc_uint<1>) ~(comp_r_27_124) + (sc_uint<1>) ~(comp_r_29_124) + (sc_uint<1>) ~(comp_r_31_124) + (sc_uint<1>) ~(comp_r_33_124) + (sc_uint<1>) ~(comp_r_35_124) + (sc_uint<1>) ~(comp_r_37_124) + (sc_uint<1>) ~(comp_r_39_124) + (sc_uint<1>) ~(comp_r_41_124) + (sc_uint<1>) ~(comp_r_43_124) + (sc_uint<1>) ~(comp_r_45_124) + (sc_uint<1>) ~(comp_r_47_124) + (sc_uint<1>) ~(comp_r_49_124) + (sc_uint<1>) ~(comp_r_51_124) + (sc_uint<1>) ~(comp_r_53_124) + (sc_uint<1>) ~(comp_r_55_124) + (sc_uint<1>) ~(comp_r_57_124) + (sc_uint<1>) ~(comp_r_59_124) + (sc_uint<1>) ~(comp_r_61_124) + (sc_uint<1>) ~(comp_r_63_124) + (sc_uint<1>) ~(comp_r_65_124) + (sc_uint<1>) ~(comp_r_67_124) + (sc_uint<1>) ~(comp_r_69_124) + (sc_uint<1>) ~(comp_r_71_124) + (sc_uint<1>) ~(comp_r_73_124) + (sc_uint<1>) ~(comp_r_75_124) + (sc_uint<1>) ~(comp_r_77_124) + (sc_uint<1>) ~(comp_r_79_124) + (sc_uint<1>) ~(comp_r_81_124) + (sc_uint<1>) ~(comp_r_83_124) + (sc_uint<1>) ~(comp_r_85_124) + (sc_uint<1>) ~(comp_r_87_124) + (sc_uint<1>) ~(comp_r_89_124) + (sc_uint<1>) ~(comp_r_91_124) + (sc_uint<1>) ~(comp_r_93_124) + (sc_uint<1>) ~(comp_r_95_124) + (sc_uint<1>) ~(comp_r_97_124) + (sc_uint<1>) ~(comp_r_99_124) + (sc_uint<1>) ~(comp_r_101_124) + (sc_uint<1>) ~(comp_r_103_124) + (sc_uint<1>) ~(comp_r_105_124) + (sc_uint<1>) ~(comp_r_107_124) + (sc_uint<1>) ~(comp_r_109_124) + (sc_uint<1>) ~(comp_r_111_124) + (sc_uint<1>) ~(comp_r_113_124) + (sc_uint<1>) ~(comp_r_115_124) + (sc_uint<1>) ~(comp_r_117_124) + (sc_uint<1>) ~(comp_r_119_124) + (sc_uint<1>) ~(comp_r_121_124) + (sc_uint<1>) ~(comp_r_123_124);
    position[125] = 63 + (sc_uint<1>) ~(comp_r_1_125) + (sc_uint<1>) ~(comp_r_3_125) + (sc_uint<1>) ~(comp_r_5_125) + (sc_uint<1>) ~(comp_r_7_125) + (sc_uint<1>) ~(comp_r_9_125) + (sc_uint<1>) ~(comp_r_11_125) + (sc_uint<1>) ~(comp_r_13_125) + (sc_uint<1>) ~(comp_r_15_125) + (sc_uint<1>) ~(comp_r_17_125) + (sc_uint<1>) ~(comp_r_19_125) + (sc_uint<1>) ~(comp_r_21_125) + (sc_uint<1>) ~(comp_r_23_125) + (sc_uint<1>) ~(comp_r_25_125) + (sc_uint<1>) ~(comp_r_27_125) + (sc_uint<1>) ~(comp_r_29_125) + (sc_uint<1>) ~(comp_r_31_125) + (sc_uint<1>) ~(comp_r_33_125) + (sc_uint<1>) ~(comp_r_35_125) + (sc_uint<1>) ~(comp_r_37_125) + (sc_uint<1>) ~(comp_r_39_125) + (sc_uint<1>) ~(comp_r_41_125) + (sc_uint<1>) ~(comp_r_43_125) + (sc_uint<1>) ~(comp_r_45_125) + (sc_uint<1>) ~(comp_r_47_125) + (sc_uint<1>) ~(comp_r_49_125) + (sc_uint<1>) ~(comp_r_51_125) + (sc_uint<1>) ~(comp_r_53_125) + (sc_uint<1>) ~(comp_r_55_125) + (sc_uint<1>) ~(comp_r_57_125) + (sc_uint<1>) ~(comp_r_59_125) + (sc_uint<1>) ~(comp_r_61_125) + (sc_uint<1>) ~(comp_r_63_125) + (sc_uint<1>) ~(comp_r_65_125) + (sc_uint<1>) ~(comp_r_67_125) + (sc_uint<1>) ~(comp_r_69_125) + (sc_uint<1>) ~(comp_r_71_125) + (sc_uint<1>) ~(comp_r_73_125) + (sc_uint<1>) ~(comp_r_75_125) + (sc_uint<1>) ~(comp_r_77_125) + (sc_uint<1>) ~(comp_r_79_125) + (sc_uint<1>) ~(comp_r_81_125) + (sc_uint<1>) ~(comp_r_83_125) + (sc_uint<1>) ~(comp_r_85_125) + (sc_uint<1>) ~(comp_r_87_125) + (sc_uint<1>) ~(comp_r_89_125) + (sc_uint<1>) ~(comp_r_91_125) + (sc_uint<1>) ~(comp_r_93_125) + (sc_uint<1>) ~(comp_r_95_125) + (sc_uint<1>) ~(comp_r_97_125) + (sc_uint<1>) ~(comp_r_99_125) + (sc_uint<1>) ~(comp_r_101_125) + (sc_uint<1>) ~(comp_r_103_125) + (sc_uint<1>) ~(comp_r_105_125) + (sc_uint<1>) ~(comp_r_107_125) + (sc_uint<1>) ~(comp_r_109_125) + (sc_uint<1>) ~(comp_r_111_125) + (sc_uint<1>) ~(comp_r_113_125) + (sc_uint<1>) ~(comp_r_115_125) + (sc_uint<1>) ~(comp_r_117_125) + (sc_uint<1>) ~(comp_r_119_125) + (sc_uint<1>) ~(comp_r_121_125) + (sc_uint<1>) ~(comp_r_123_125) + comp_r_125_126;
    position[126] = 63 + (sc_uint<1>) ~(comp_r_1_126) + (sc_uint<1>) ~(comp_r_3_126) + (sc_uint<1>) ~(comp_r_5_126) + (sc_uint<1>) ~(comp_r_7_126) + (sc_uint<1>) ~(comp_r_9_126) + (sc_uint<1>) ~(comp_r_11_126) + (sc_uint<1>) ~(comp_r_13_126) + (sc_uint<1>) ~(comp_r_15_126) + (sc_uint<1>) ~(comp_r_17_126) + (sc_uint<1>) ~(comp_r_19_126) + (sc_uint<1>) ~(comp_r_21_126) + (sc_uint<1>) ~(comp_r_23_126) + (sc_uint<1>) ~(comp_r_25_126) + (sc_uint<1>) ~(comp_r_27_126) + (sc_uint<1>) ~(comp_r_29_126) + (sc_uint<1>) ~(comp_r_31_126) + (sc_uint<1>) ~(comp_r_33_126) + (sc_uint<1>) ~(comp_r_35_126) + (sc_uint<1>) ~(comp_r_37_126) + (sc_uint<1>) ~(comp_r_39_126) + (sc_uint<1>) ~(comp_r_41_126) + (sc_uint<1>) ~(comp_r_43_126) + (sc_uint<1>) ~(comp_r_45_126) + (sc_uint<1>) ~(comp_r_47_126) + (sc_uint<1>) ~(comp_r_49_126) + (sc_uint<1>) ~(comp_r_51_126) + (sc_uint<1>) ~(comp_r_53_126) + (sc_uint<1>) ~(comp_r_55_126) + (sc_uint<1>) ~(comp_r_57_126) + (sc_uint<1>) ~(comp_r_59_126) + (sc_uint<1>) ~(comp_r_61_126) + (sc_uint<1>) ~(comp_r_63_126) + (sc_uint<1>) ~(comp_r_65_126) + (sc_uint<1>) ~(comp_r_67_126) + (sc_uint<1>) ~(comp_r_69_126) + (sc_uint<1>) ~(comp_r_71_126) + (sc_uint<1>) ~(comp_r_73_126) + (sc_uint<1>) ~(comp_r_75_126) + (sc_uint<1>) ~(comp_r_77_126) + (sc_uint<1>) ~(comp_r_79_126) + (sc_uint<1>) ~(comp_r_81_126) + (sc_uint<1>) ~(comp_r_83_126) + (sc_uint<1>) ~(comp_r_85_126) + (sc_uint<1>) ~(comp_r_87_126) + (sc_uint<1>) ~(comp_r_89_126) + (sc_uint<1>) ~(comp_r_91_126) + (sc_uint<1>) ~(comp_r_93_126) + (sc_uint<1>) ~(comp_r_95_126) + (sc_uint<1>) ~(comp_r_97_126) + (sc_uint<1>) ~(comp_r_99_126) + (sc_uint<1>) ~(comp_r_101_126) + (sc_uint<1>) ~(comp_r_103_126) + (sc_uint<1>) ~(comp_r_105_126) + (sc_uint<1>) ~(comp_r_107_126) + (sc_uint<1>) ~(comp_r_109_126) + (sc_uint<1>) ~(comp_r_111_126) + (sc_uint<1>) ~(comp_r_113_126) + (sc_uint<1>) ~(comp_r_115_126) + (sc_uint<1>) ~(comp_r_117_126) + (sc_uint<1>) ~(comp_r_119_126) + (sc_uint<1>) ~(comp_r_121_126) + (sc_uint<1>) ~(comp_r_123_126) + (sc_uint<1>) ~(comp_r_125_126);
    position[127] = (sc_uint<7>) 127;

 // Multiplex the inputs in the order
    PS_struct<1,Q,6> temp[127];
#pragma HLS ARRAY_PARTITION variable=temp complete dim=1
    temp[0] = input[0];
    temp[1] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 1);
    temp[2] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 2);
    temp[3] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 3);
    temp[4] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 4);
    temp[5] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 5);
    temp[6] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 6);
    temp[7] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 7);
    temp[8] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 8);
    temp[9] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 9);
    temp[10] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 10);
    temp[11] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 11);
    temp[12] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 12);
    temp[13] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 13);
    temp[14] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 14);
    temp[15] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 15);
    temp[16] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 16);
    temp[17] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 17);
    temp[18] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 18);
    temp[19] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 19);
    temp[20] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 20);
    temp[21] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 21);
    temp[22] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 22);
    temp[23] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 23);
    temp[24] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 24);
    temp[25] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 25);
    temp[26] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 26);
    temp[27] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 27);
    temp[28] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 28);
    temp[29] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 29);
    temp[30] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 30);
    temp[31] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 31);
    temp[32] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 32);
    temp[33] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 33);
    temp[34] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 34);
    temp[35] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 35);
    temp[36] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 36);
    temp[37] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 37);
    temp[38] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 38);
    temp[39] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 39);
    temp[40] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 40);
    temp[41] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 41);
    temp[42] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 42);
    temp[43] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 43);
    temp[44] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 44);
    temp[45] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 45);
    temp[46] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 46);
    temp[47] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 47);
    temp[48] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 48);
    temp[49] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 49);
    temp[50] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 50);
    temp[51] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 51);
    temp[52] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 52);
    temp[53] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 53);
    temp[54] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 54);
    temp[55] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 55);
    temp[56] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 56);
    temp[57] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 57);
    temp[58] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 58);
    temp[59] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 59);
    temp[60] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 60);
    temp[61] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 61);
    temp[62] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 62);
    temp[63] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 63);
    temp[64] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 64);
    temp[65] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 65);
    temp[66] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 66);
    temp[67] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 67);
    temp[68] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 68);
    temp[69] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 69);
    temp[70] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 70);
    temp[71] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 71);
    temp[72] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 72);
    temp[73] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 73);
    temp[74] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 74);
    temp[75] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 75);
    temp[76] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 76);
    temp[77] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 77);
    temp[78] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 78);
    temp[79] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 79);
    temp[80] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 80);
    temp[81] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 81);
    temp[82] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 82);
    temp[83] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 83);
    temp[84] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 84);
    temp[85] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 85);
    temp[86] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 86);
    temp[87] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 87);
    temp[88] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 88);
    temp[89] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 89);
    temp[90] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 90);
    temp[91] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 91);
    temp[92] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 92);
    temp[93] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 93);
    temp[94] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 94);
    temp[95] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 95);
    temp[96] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 96);
    temp[97] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 97);
    temp[98] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 98);
    temp[99] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 99);
    temp[100] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 100);
    temp[101] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 101);
    temp[102] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 102);
    temp[103] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 103);
    temp[104] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 104);
    temp[105] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 105);
    temp[106] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 106);
    temp[107] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 107);
    temp[108] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 108);
    temp[109] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 109);
    temp[110] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 110);
    temp[111] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 111);
    temp[112] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 112);
    temp[113] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 113);
    temp[114] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 114);
    temp[115] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 115);
    temp[116] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 116);
    temp[117] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 117);
    temp[118] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 118);
    temp[119] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 119);
    temp[120] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 120);
    temp[121] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 121);
    temp[122] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 122);
    temp[123] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 123);
    temp[124] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 124);
    temp[125] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 125);
    temp[126] = RO_MUX128 < PS_struct<1,Q,6> > (input, position, (sc_uint<7>) 126);

// output the L smallest input
    if (fb == 0){
    output[0] = temp[63]; 
    output[1] = temp[64]; 
    output[2] = temp[65]; 
    output[3] = temp[66]; 
    output[4] = temp[67]; 
    output[5] = temp[68]; 
    output[6] = temp[69]; 
    output[7] = temp[70]; 
    output[8] = temp[71]; 
    output[9] = temp[72]; 
    output[10] = temp[73]; 
    output[11] = temp[74]; 
    output[12] = temp[75]; 
    output[13] = temp[76]; 
    output[14] = temp[77]; 
    output[15] = temp[78]; 
    output[16] = temp[79]; 
    output[17] = temp[80]; 
    output[18] = temp[81]; 
    output[19] = temp[82]; 
    output[20] = temp[83]; 
    output[21] = temp[84]; 
    output[22] = temp[85]; 
    output[23] = temp[86]; 
    output[24] = temp[87]; 
    output[25] = temp[88]; 
    output[26] = temp[89]; 
    output[27] = temp[90]; 
    output[28] = temp[91]; 
    output[29] = temp[92]; 
    output[30] = temp[93]; 
    output[31] = temp[94]; 
    output[32] = temp[95]; 
    output[33] = temp[96]; 
    output[34] = temp[97]; 
    output[35] = temp[98]; 
    output[36] = temp[99]; 
    output[37] = temp[100]; 
    output[38] = temp[101]; 
    output[39] = temp[102]; 
    output[40] = temp[103]; 
    output[41] = temp[104]; 
    output[42] = temp[105]; 
    output[43] = temp[106]; 
    output[44] = temp[107]; 
    output[45] = temp[108]; 
    output[46] = temp[109]; 
    output[47] = temp[110]; 
    output[48] = temp[111]; 
    output[49] = temp[112]; 
    output[50] = temp[113]; 
    output[51] = temp[114]; 
    output[52] = temp[115]; 
    output[53] = temp[116]; 
    output[54] = temp[117]; 
    output[55] = temp[118]; 
    output[56] = temp[119]; 
    output[57] = temp[120]; 
    output[58] = temp[121]; 
    output[59] = temp[122]; 
    output[60] = temp[123]; 
    output[61] = temp[124]; 
    output[62] = temp[125]; 
    output[63] = temp[126]; 
    }
    else{
    output[0] = temp[0]; 
    output[1] = temp[1]; 
    output[2] = temp[2]; 
    output[3] = temp[3]; 
    output[4] = temp[4]; 
    output[5] = temp[5]; 
    output[6] = temp[6]; 
    output[7] = temp[7]; 
    output[8] = temp[8]; 
    output[9] = temp[9]; 
    output[10] = temp[10]; 
    output[11] = temp[11]; 
    output[12] = temp[12]; 
    output[13] = temp[13]; 
    output[14] = temp[14]; 
    output[15] = temp[15]; 
    output[16] = temp[16]; 
    output[17] = temp[17]; 
    output[18] = temp[18]; 
    output[19] = temp[19]; 
    output[20] = temp[20]; 
    output[21] = temp[21]; 
    output[22] = temp[22]; 
    output[23] = temp[23]; 
    output[24] = temp[24]; 
    output[25] = temp[25]; 
    output[26] = temp[26]; 
    output[27] = temp[27]; 
    output[28] = temp[28]; 
    output[29] = temp[29]; 
    output[30] = temp[30]; 
    output[31] = temp[31]; 
    output[32] = temp[32]; 
    output[33] = temp[33]; 
    output[34] = temp[34]; 
    output[35] = temp[35]; 
    output[36] = temp[36]; 
    output[37] = temp[37]; 
    output[38] = temp[38]; 
    output[39] = temp[39]; 
    output[40] = temp[40]; 
    output[41] = temp[41]; 
    output[42] = temp[42]; 
    output[43] = temp[43]; 
    output[44] = temp[44]; 
    output[45] = temp[45]; 
    output[46] = temp[46]; 
    output[47] = temp[47]; 
    output[48] = temp[48]; 
    output[49] = temp[49]; 
    output[50] = temp[50]; 
    output[51] = temp[51]; 
    output[52] = temp[52]; 
    output[53] = temp[53]; 
    output[54] = temp[54]; 
    output[55] = temp[55]; 
    output[56] = temp[56]; 
    output[57] = temp[57]; 
    output[58] = temp[58]; 
    output[59] = temp[59]; 
    output[60] = temp[60]; 
    output[61] = temp[61]; 
    output[62] = temp[62]; 
    output[63] = temp[63]; 
    }
}

template < int L, int Q, int LOG2L>
void RANKORDER_SORT (PS_struct<1,Q,LOG2L> input[2*L], PS_struct<1,Q,LOG2L> output[L], sc_biguint<1> fb)
{
#pragma HLS INLINE
#if L_SIZE == 2
	RANKORDER_SORT_L2 < Q > ( input, output, fb);
#elif L_SIZE == 4
	RANKORDER_SORT_L4 < Q > ( input, output, fb);
#elif L_SIZE == 8
	RANKORDER_SORT_L8 < Q > ( input, output, fb);
#elif L_SIZE == 16
	RANKORDER_SORT_L16 < Q > ( input, output, fb);
#elif L_SIZE == 32
	RANKORDER_SORT_L32 < Q > ( input, output, fb);
#elif L_SIZE == 64
	RANKORDER_SORT_L64 < Q > ( input, output, fb);
#endif

}

//*************************************************************************//
//** 								FUNCTION							 **//
//*************************************************************************//


template <int L, int Q, int LOG2L, int D>
void List_R_P1 ( LLR_struct<1,Q,LOG2L> input[L], PS_struct<1,Q,LOG2L> output[L], sc_biguint<1> fb, sc_biguint<LOG2L*D> LIST_STACK[L] )
{
#if PARALLEL_TP == 0
	#pragma HLS INLINE off
#else
	#pragma HLS INLINE
#endif
#pragma HLS ARRAY_PARTITION variable=input complete dim=1
#pragma HLS ARRAY_PARTITION variable=output complete dim=1
#pragma HLS ARRAY_PARTITION variable=LIST_STACK complete dim=1

	PS_struct<1,Q,LOG2L> metric_out[2*L];

	// Update the metric and create L new path
	UPDATE_METRIC<L, Q, LOG2L> (input , metric_out, fb);

 /*///////////////   DEBUG   /////////////////////
cout << "%%% OUT METRIC UPDATE %%% "<< endl;
	for(sc_uint<LOG2L+2> j = 0; j < 2*L; j++)
	{
		cout << " #L "<< j << " |" << endl;
		cout << "	bit=  "; SHOW_BITS<PAR >( (sc_biguint<1>) metric_out[j].bit); cout << endl;
		cout << "	met=  "; SHOW_LLRS<PAR, LLR_BITS >( (sc_bigint< Q>) metric_out[j].metric); cout << endl;
		cout << "	pat=  " << metric_out[j].path << endl;
	}
///////////////   DEBUG   /////////////////////*/

	// Sort the 2L path metric and output the L smallest
#if SORTER == 0
	BUBBLE_SORT<L, Q, LOG2L> (metric_out, output);
#elif SORTER == 1
	CUSTOM_SORT<L, Q, LOG2L>(metric_out, output, fb);
#elif SORTER == 2
	RANKORDER_SORT<L, Q, LOG2L>(metric_out, output, fb);
#endif

	/* ///////////////   DEBUG   /////////////////////
	cout << "%%% OUT SORTER %%% "<< endl;
		for(sc_uint<LOG2L+1> j = 0; j < L; j++)
		{
			cout << " #L "<< j << " |" << endl;
			cout << "	bit=  "; SHOW_BITS<PAR >( (sc_biguint<1>) output[j].bit); cout << endl;
			cout << "	met=  "; SHOW_LLRS<PAR, LLR_BITS >( (sc_bigint< Q>) output[j].metric); cout << endl;
			cout << "	pat=  " << output[j].path << endl;
		}
	///////////////   DEBUG   /////////////////////*/

	// List Stack management : make copies of stack according to path
	sc_bigint<LOG2L*D> temp_list_stack[L];
#pragma HLS ARRAY_PARTITION variable=temp_list_stack complete dim=1
	for(sc_uint<LOG2L+1> j = 0; j < L; j++)
	{
	#pragma HLS UNROLL
		temp_list_stack[j] = (sc_bigint<LOG2L*D>)LIST_STACK[j];
	}
	for(sc_uint<LOG2L+1> j = 0; j < L; j++)
	{
	#pragma HLS UNROLL
		sc_uint<LOG2L> at = output[j].path;
		LIST_STACK[j] = (sc_biguint<LOG2L*D>) MUX_LIST <L, LOG2L*D, LOG2L> ( temp_list_stack, at );
		//LIST_STACK[j] = LIST_STACK[at];
	}

}


template <int L, int Q, int LOG2L, int D>
void List_R_P2 ( LLR_struct<2,Q,LOG2L> input[L], PS_struct<2,Q,LOG2L> output[L], sc_biguint<2> fb, sc_biguint<LOG2L*D> LIST_STACK[L] )
{
#if PARALLEL_TP == 0
	#pragma HLS INLINE off
#else
	#pragma HLS INLINE
#endif
#pragma HLS ARRAY_PARTITION variable=input complete dim=1
#pragma HLS ARRAY_PARTITION variable=output complete dim=1

		sc_bigint< (1*Q) > la[L];
		sc_bigint< (1*Q) > lb[L];
		LLR_struct< 1 ,Q,LOG2L> First_R_P1_in[L];
// F
		for(COUNTER j = 0; j < L_SIZE; j++){
		#pragma HLS UNROLL
			la[j] = (input[j].llr).range( (1*Q) -1,0);
			lb[j] = (input[j].llr).range(2* (1*Q) -1, (1*Q) );

			sc_bigint< (1*Q) > la1 = Function_F< 1 ,Q> (la[j],lb[j]);

			First_R_P1_in[j].llr = la1;
			First_R_P1_in[j].metric = input[j].metric;
			First_R_P1_in[j].path = input[j].path;

			//// List Stack management  ////
			LIST_STACK[j] = push_stack <LOG2_L, _DEPTH >( LIST_STACK[j] , (sc_biguint<LOG2_L>) j );
		}

		PS_struct< 1 ,Q,LOG2L> First_R_P1_out[L];
// First P1
		List_R_P1 <L,Q,LOG2L,D> (First_R_P1_in, First_R_P1_out, (sc_biguint<1>)fb[0], LIST_STACK);

		sc_bigint< 1 > sa1[L];
		LLR_struct< 1 ,Q,LOG2L> Second_R_P1_in[L];
// G
		for(COUNTER j = 0; j < L_SIZE; j++){
		#pragma HLS UNROLL
			sa1[j] = (sc_biguint< 1 >) First_R_P1_out[j].bit;

			sc_biguint<LOG2_L> a_l_llr = read_stack< LOG2_L, _DEPTH > (LIST_STACK[j], 2);

			sc_bigint< (1*Q) > i_la = MUX_LIST < L, 1*Q, LOG2L> ( la, a_l_llr);
			sc_bigint< (1*Q) > i_lb = MUX_LIST < L, 1*Q, LOG2L> ( lb, a_l_llr);

			sc_bigint< (1*Q) > lb1 = Function_G< 1 ,Q> (i_la,i_lb,(sc_biguint< 1 >)sa1[j]);

			Second_R_P1_in[j].llr = lb1;
			Second_R_P1_in[j].metric = First_R_P1_out[j].metric;
			Second_R_P1_in[j].path = First_R_P1_out[j].path;

			//// List Stack management  ////
			LIST_STACK[j] = write_stack <LOG2_L, _DEPTH > ( LIST_STACK[j] , (sc_biguint<LOG2_L>) j );
		}

		PS_struct< 1 ,Q,LOG2L> Second_R_P1_out[L];
// Second P1
		List_R_P1 <L,Q,LOG2L,D> (Second_R_P1_in, Second_R_P1_out, (sc_biguint<1>)fb[1], LIST_STACK);

// XOR
		for(COUNTER j = 0; j < L_SIZE; j++){
		#pragma HLS UNROLL
			sc_biguint<LOG2_L> a_list = read_stack< LOG2_L, _DEPTH > (LIST_STACK[j], 1);

			sc_biguint< 1 > i_sa = (sc_biguint< 1 >) MUX_LIST <L, 1, LOG2L> ( sa1, a_list);
			sc_biguint< 1 > i_sb = Second_R_P1_out[j].bit;

			sc_biguint< 1 > sa = VECTOR_XOR< 1 >(i_sa,i_sb);
			sc_biguint< 1 > sb = i_sb;

			output[j].bit = (sb , sa);
			output[j].metric = Second_R_P1_out[j].metric;
			output[j].path = Second_R_P1_out[j].path;

			//// List Stack management  ////
			LIST_STACK[j] = pop_stack <LOG2_L, _DEPTH >( LIST_STACK[j] , (sc_biguint<LOG2_L>) 0 );
		}
}


template <int L, int Q, int LOG2L, int D>
void List_R_P4 ( LLR_struct<4,Q,LOG2L> input[L], PS_struct<4,Q,LOG2L> output[L], sc_biguint<4> fb, sc_biguint<LOG2L*D> LIST_STACK[L] )
{
#if PARALLEL_TP == 0
	#pragma HLS INLINE off
#else
	#pragma HLS INLINE
#endif
#pragma HLS ARRAY_PARTITION variable=input complete dim=1
#pragma HLS ARRAY_PARTITION variable=output complete dim=1

		sc_bigint< (2*Q) > la[L];
		sc_bigint< (2*Q) > lb[L];
		LLR_struct< 2 ,Q,LOG2L> First_R_P1_in[L];
// F
		for(COUNTER j = 0; j < L_SIZE; j++){
		#pragma HLS UNROLL
			la[j] = (input[j].llr).range( (2*Q) -1,0);
			lb[j] = (input[j].llr).range(2* (2*Q) -1, (2*Q) );

			sc_bigint< (2*Q) > la1 = Function_F< 2 ,Q> (la[j],lb[j]);

			First_R_P1_in[j].llr = la1;
			First_R_P1_in[j].metric = input[j].metric;
			First_R_P1_in[j].path = input[j].path;

			//// List Stack management  ////
			LIST_STACK[j] = push_stack <LOG2_L, _DEPTH >( LIST_STACK[j] , (sc_biguint<LOG2_L>) j );
		}

		PS_struct< 2 ,Q,LOG2L> First_R_P1_out[L];
// First P1
		List_R_P2 <L,Q,LOG2L,D> (First_R_P1_in, First_R_P1_out, (sc_biguint< 2 >)fb.range(1, 0), LIST_STACK );

		sc_bigint< 2 > sa1[L];
		LLR_struct< 2 ,Q,LOG2L> Second_R_P1_in[L];
// G
		for(COUNTER j = 0; j < L_SIZE; j++){
		#pragma HLS UNROLL
			sa1[j] = (sc_biguint< 2 >) First_R_P1_out[j].bit;

			sc_biguint<LOG2_L> a_l_llr = read_stack< LOG2_L, _DEPTH > (LIST_STACK[j], 2);

			sc_bigint< (2*Q) > i_la = MUX_LIST < L, 2*Q, LOG2L> ( la, a_l_llr);
			sc_bigint< (2*Q) > i_lb = MUX_LIST < L, 2*Q, LOG2L> ( lb, a_l_llr);

			sc_bigint< (2*Q) > lb1 = Function_G< 2 ,Q> (i_la,i_lb,(sc_biguint< 2 >)sa1[j]);

			Second_R_P1_in[j].llr = lb1;
			Second_R_P1_in[j].metric = First_R_P1_out[j].metric;
			Second_R_P1_in[j].path = First_R_P1_out[j].path;

			//// List Stack management  ////
			LIST_STACK[j] = write_stack <LOG2_L, _DEPTH > ( LIST_STACK[j] , (sc_biguint<LOG2_L>) j );
		}

		PS_struct< 2 ,Q,LOG2L> Second_R_P1_out[L];
// Second P1
		List_R_P2 <L,Q,LOG2L,D> (Second_R_P1_in, Second_R_P1_out, (sc_biguint< 2 >)fb.range(3, 2), LIST_STACK );

// XOR
		for(COUNTER j = 0; j < L_SIZE; j++){
		#pragma HLS UNROLL
			sc_biguint<LOG2_L> a_list = read_stack< LOG2_L, _DEPTH > (LIST_STACK[j], 1);

			sc_biguint< 2 > i_sa = (sc_biguint< 2 >) MUX_LIST <L, 2, LOG2L> ( sa1, a_list);
			sc_biguint< 2 > i_sb = Second_R_P1_out[j].bit;

			sc_biguint< 2 > sa = VECTOR_XOR< 2 >(i_sa,i_sb);
			sc_biguint< 2 > sb = i_sb;

			output[j].bit = (sb , sa);
			output[j].metric = Second_R_P1_out[j].metric;
			output[j].path = Second_R_P1_out[j].path;

			//// List Stack management  ////
			LIST_STACK[j] = pop_stack <LOG2_L, _DEPTH >( LIST_STACK[j] , (sc_biguint<LOG2_L>) 0 );
		}

}


template <int L, int Q, int LOG2L, int D>
void List_R_P8 ( LLR_struct<8,Q,LOG2L> input[L], PS_struct<8,Q,LOG2L> output[L], sc_biguint<8> fb, sc_biguint<LOG2L*D> LIST_STACK[L] )
{
#if PARALLEL_TP == 0
	#pragma HLS INLINE off
#else
	#pragma HLS INLINE
#endif
#pragma HLS ARRAY_PARTITION variable=input complete dim=1
#pragma HLS ARRAY_PARTITION variable=output complete dim=1

		sc_bigint< (4*Q) > la[L];
		sc_bigint< (4*Q) > lb[L];
		LLR_struct< 4 ,Q,LOG2L> First_R_P1_in[L];
// F
		for(COUNTER j = 0; j < L_SIZE; j++){
		#pragma HLS UNROLL
			la[j] = (input[j].llr).range( (4*Q) -1,0);
			lb[j] = (input[j].llr).range(2* (4*Q) -1, (4*Q) );

			sc_bigint< (4*Q) > la1 = Function_F< 4 ,Q> (la[j],lb[j]);

			First_R_P1_in[j].llr = la1;
			First_R_P1_in[j].metric = input[j].metric;
			First_R_P1_in[j].path = input[j].path;

			//// List Stack management  ////
			LIST_STACK[j] = push_stack <LOG2_L, _DEPTH >( LIST_STACK[j] , (sc_biguint<LOG2_L>) j );
		}

		PS_struct< 4 ,Q,LOG2L> First_R_P1_out[L];
// First P1
		List_R_P4 <L,Q,LOG2L,D> (First_R_P1_in, First_R_P1_out, (sc_biguint< 4 >)fb.range(3, 0), LIST_STACK);

		sc_bigint< 4 > sa1[L];
		LLR_struct< 4 ,Q,LOG2L> Second_R_P1_in[L];
// G
		for(COUNTER j = 0; j < L_SIZE; j++){
		#pragma HLS UNROLL
			sa1[j] = First_R_P1_out[j].bit;

			sc_biguint<LOG2_L> a_l_llr = read_stack< LOG2_L, _DEPTH > (LIST_STACK[j], 2);

			sc_bigint< (4*Q) > i_la = MUX_LIST < L, 4*Q, LOG2L> ( la, a_l_llr);
			sc_bigint< (4*Q) > i_lb = MUX_LIST < L, 4*Q, LOG2L> ( lb, a_l_llr);

			sc_bigint< (4*Q) > lb1 = Function_G< 4 ,Q> (i_la,i_lb,sa1[j]);

			Second_R_P1_in[j].llr = lb1;
			Second_R_P1_in[j].metric = First_R_P1_out[j].metric;
			Second_R_P1_in[j].path = First_R_P1_out[j].path;

			//// List Stack management  ////
			LIST_STACK[j] = write_stack <LOG2_L, _DEPTH > ( LIST_STACK[j] , (sc_biguint<LOG2_L>) j );
		}

		PS_struct< 4 ,Q,LOG2L> Second_R_P1_out[L];
// Second P1
		List_R_P4 <L,Q,LOG2L,D> (Second_R_P1_in, Second_R_P1_out, (sc_biguint< 4 >)fb.range(7, 4), LIST_STACK);

// XOR
		for(COUNTER j = 0; j < L_SIZE; j++){
		#pragma HLS UNROLL
			sc_biguint<LOG2_L> a_list = read_stack< LOG2_L, _DEPTH > (LIST_STACK[j], 1);
			sc_biguint< 4 > i_sa = (sc_biguint< 4 >) MUX_LIST <L, 4, LOG2L> ( sa1, a_list);
			sc_biguint< 4 > i_sb = Second_R_P1_out[j].bit;

			sc_biguint< 4 > sa = VECTOR_XOR< 4 >(i_sa,i_sb);
			sc_biguint< 4 > sb = i_sb;

			output[j].bit = (sb , sa);
			output[j].metric = Second_R_P1_out[j].metric;
			output[j].path = Second_R_P1_out[j].path;

			//// List Stack management  ////
			LIST_STACK[j] = pop_stack <LOG2_L, _DEPTH >( LIST_STACK[j] , (sc_biguint<LOG2_L>) 0 );
		}
}


template <int L, int Q, int LOG2L, int D>
void List_R_P16 ( LLR_struct<16,Q,LOG2L> input[L], PS_struct<16,Q,LOG2L> output[L], sc_biguint<16> fb, sc_biguint<LOG2L*D> LIST_STACK[L] )
{
#if PARALLEL_TP == 0
	#pragma HLS INLINE off
#else
	#pragma HLS INLINE
#endif
#pragma HLS ARRAY_PARTITION variable=input complete dim=1
#pragma HLS ARRAY_PARTITION variable=output complete dim=1

		sc_bigint< (8*Q) > la[L];
		sc_bigint< (8*Q) > lb[L];
		LLR_struct< 8 ,Q,LOG2L> First_R_P1_in[L];
// F
		for(COUNTER j = 0; j < L_SIZE; j++){
		#pragma HLS UNROLL
			la[j] = (input[j].llr).range( (8*Q) -1,0);
			lb[j] = (input[j].llr).range(2* (8*Q) -1, (8*Q) );

			sc_bigint< (8*Q) > la1 = Function_F< 8 ,Q> (la[j],lb[j]);

			First_R_P1_in[j].llr = la1;
			First_R_P1_in[j].metric = input[j].metric;
			First_R_P1_in[j].path = input[j].path;

			//// List Stack management  ////
			LIST_STACK[j] = push_stack <LOG2_L, _DEPTH >( LIST_STACK[j] , (sc_biguint<LOG2_L>) j );
		}

		PS_struct< 8 ,Q,LOG2L> First_R_P1_out[L];
// First P1
		List_R_P8 <L,Q,LOG2L,D> (First_R_P1_in, First_R_P1_out, (sc_biguint< 8 >)fb.range(7, 0), LIST_STACK);

		sc_bigint< 8 > sa1[L];
		LLR_struct< 8 ,Q,LOG2L> Second_R_P1_in[L];
// G
		for(COUNTER j = 0; j < L_SIZE; j++){
		#pragma HLS UNROLL
			sa1[j] = First_R_P1_out[j].bit;

			sc_biguint<LOG2_L> a_l_llr = read_stack< LOG2_L, _DEPTH > (LIST_STACK[j], 2);

			sc_bigint< (8*Q) > i_la = MUX_LIST < L, 8*Q, LOG2L> ( la, a_l_llr);
			sc_bigint< (8*Q) > i_lb = MUX_LIST < L, 8*Q, LOG2L> ( lb, a_l_llr);

			sc_bigint< (8*Q) > lb1 = Function_G< 8 ,Q> (i_la,i_lb,sa1[j]);

			Second_R_P1_in[j].llr = lb1;
			Second_R_P1_in[j].metric = First_R_P1_out[j].metric;
			Second_R_P1_in[j].path = First_R_P1_out[j].path;

			//// List Stack management  ////
			LIST_STACK[j] = write_stack <LOG2_L, _DEPTH > ( LIST_STACK[j] , (sc_biguint<LOG2_L>) j );
		}

		PS_struct< 8 ,Q,LOG2L> Second_R_P1_out[L];
// Second P1
		List_R_P8 <L,Q,LOG2L,D> (Second_R_P1_in, Second_R_P1_out, (sc_biguint< 8 >)fb.range(15, 8), LIST_STACK);

// XOR
		for(COUNTER j = 0; j < L_SIZE; j++){
		#pragma HLS UNROLL
			sc_biguint<LOG2_L> a_list = read_stack< LOG2_L, _DEPTH > (LIST_STACK[j], 1);

			sc_biguint< 8 > i_sa = (sc_biguint< 8 >) MUX_LIST <L, 8, LOG2L> ( sa1, a_list);
			sc_biguint< 8 > i_sb = Second_R_P1_out[j].bit;

			sc_biguint< 8 > sa = VECTOR_XOR< 8 >(i_sa,i_sb);
			sc_biguint< 8 > sb = i_sb;

			output[j].bit = (sb , sa);
			output[j].metric = Second_R_P1_out[j].metric;
			output[j].path = Second_R_P1_out[j].path;

			//// List Stack management  ////
			LIST_STACK[j] = pop_stack <LOG2_L, _DEPTH >( LIST_STACK[j] , (sc_biguint<LOG2_L>) 0 );
		}
}


template <int L, int Q, int LOG2L, int D>
void List_R_P32 ( LLR_struct<32,Q,LOG2L> input[L], PS_struct<32,Q,LOG2L> output[L], sc_biguint<32> fb, sc_biguint<LOG2L*D> LIST_STACK[L] )
{
#if PARALLEL_TP == 0
	#pragma HLS INLINE off
#else
	#pragma HLS INLINE
#endif
#pragma HLS ARRAY_PARTITION variable=input complete dim=1
#pragma HLS ARRAY_PARTITION variable=output complete dim=1

		sc_bigint< (16*Q) > la[L];
		sc_bigint< (16*Q) > lb[L];
		LLR_struct< 16 ,Q,LOG2L> First_R_P1_in[L];
// F
		for(COUNTER j = 0; j < L_SIZE; j++){
		#pragma HLS UNROLL
			la[j] = (input[j].llr).range( (16*Q) -1,0);
			lb[j] = (input[j].llr).range(2* (16*Q) -1, (16*Q) );

			sc_bigint< (16*Q) > la1 = Function_F< 16 ,Q> (la[j],lb[j]);

			First_R_P1_in[j].llr = la1;
			First_R_P1_in[j].metric = input[j].metric;
			First_R_P1_in[j].path = input[j].path;

			//// List Stack management  ////
			LIST_STACK[j] = push_stack <LOG2_L, _DEPTH >( LIST_STACK[j] , (sc_biguint<LOG2_L>) j );
		}

		PS_struct< 16 ,Q,LOG2L> First_R_P1_out[L];
// First P1
		List_R_P16 <L,Q,LOG2L,D> (First_R_P1_in, First_R_P1_out, (sc_biguint< 16 >)fb.range(15, 0), LIST_STACK);

		sc_bigint< 16 > sa1[L];
		LLR_struct< 16 ,Q,LOG2L> Second_R_P1_in[L];
// G
		for(COUNTER j = 0; j < L_SIZE; j++){
		#pragma HLS UNROLL
			sa1[j] = First_R_P1_out[j].bit;

			sc_biguint<LOG2_L> a_l_llr = read_stack< LOG2_L, _DEPTH > (LIST_STACK[j], 2);
			sc_bigint< (16*Q) > i_la = MUX_LIST < L, 16*Q, LOG2L> ( la, a_l_llr);
			sc_bigint< (16*Q) > i_lb = MUX_LIST < L, 16*Q, LOG2L> ( lb, a_l_llr);

			sc_bigint< (16*Q) > lb1 = Function_G< 16 ,Q> (i_la,i_lb,sa1[j]);

			Second_R_P1_in[j].llr = lb1;
			Second_R_P1_in[j].metric = First_R_P1_out[j].metric;
			Second_R_P1_in[j].path = First_R_P1_out[j].path;

			//// List Stack management  ////
			LIST_STACK[j] = write_stack <LOG2_L, _DEPTH > ( LIST_STACK[j] , (sc_biguint<LOG2_L>) j );
		}

		PS_struct< 16 ,Q,LOG2L> Second_R_P1_out[L];
// Second P1
		List_R_P16 <L,Q,LOG2L,D> (Second_R_P1_in, Second_R_P1_out, (sc_biguint< 16 >)fb.range(31, 16), LIST_STACK);

// XOR
		for(COUNTER j = 0; j < L_SIZE; j++){
		#pragma HLS UNROLL
			sc_biguint<LOG2_L> a_list = read_stack< LOG2_L, _DEPTH > (LIST_STACK[j], 1);

			sc_biguint< 16 > i_sa = (sc_biguint< 16 >) MUX_LIST <L, 16, LOG2L> ( sa1, a_list);
			sc_biguint< 16 > i_sb = Second_R_P1_out[j].bit;

			sc_biguint< 16 > sa = VECTOR_XOR< 16 >(i_sa,i_sb);
			sc_biguint< 16 > sb = i_sb;

			output[j].bit = (sb , sa);
			output[j].metric = Second_R_P1_out[j].metric;
			output[j].path = Second_R_P1_out[j].path;

			//// List Stack management  ////
			LIST_STACK[j] = pop_stack <LOG2_L, _DEPTH >( LIST_STACK[j] , (sc_biguint<LOG2_L>) 0 );
		}
}


template <int L, int Q, int LOG2L, int D>
void List_R_P64 ( LLR_struct<64,Q,LOG2L> input[L], PS_struct<64,Q,LOG2L> output[L], sc_biguint<64> fb, sc_biguint<LOG2L*D> LIST_STACK[L] )
{
#if PARALLEL_TP == 0
	#pragma HLS INLINE off
#else
	#pragma HLS INLINE
#endif
#pragma HLS ARRAY_PARTITION variable=input complete dim=1
#pragma HLS ARRAY_PARTITION variable=output complete dim=1

		sc_bigint< (32*Q) > la[L];
		sc_bigint< (32*Q) > lb[L];
		LLR_struct< 32 ,Q,LOG2L> First_R_P1_in[L];
// F
		for(COUNTER j = 0; j < L_SIZE; j++){
		#pragma HLS UNROLL
			la[j] = (input[j].llr).range( (32*Q) -1,0);
			lb[j] = (input[j].llr).range(2* (32*Q) -1, (32*Q) );

			sc_bigint< (32*Q) > la1 = Function_F< 32 ,Q> (la[j],lb[j]);

			First_R_P1_in[j].llr = la1;
			First_R_P1_in[j].metric = input[j].metric;
			First_R_P1_in[j].path = input[j].path;

			//// List Stack management  ////
			LIST_STACK[j] = push_stack <LOG2_L, _DEPTH >( LIST_STACK[j] , (sc_biguint<LOG2_L>) j );
		}

		PS_struct< 32 ,Q,LOG2L> First_R_P1_out[L];
// First P1
		List_R_P32 <L,Q,LOG2L,D> (First_R_P1_in, First_R_P1_out, (sc_biguint< 32 >)fb.range(31, 0), LIST_STACK);

		sc_bigint< 32 > sa1[L];
		LLR_struct< 32 ,Q,LOG2L> Second_R_P1_in[L];
// G
		for(COUNTER j = 0; j < L_SIZE; j++){
		#pragma HLS UNROLL
			sa1[j] = First_R_P1_out[j].bit;

			sc_biguint<LOG2_L> a_l_llr = read_stack< LOG2_L, _DEPTH > (LIST_STACK[j], 2);
			sc_bigint< (32*Q) > i_la = MUX_LIST < L, 32*Q, LOG2L> ( la, a_l_llr);
			sc_bigint< (32*Q) > i_lb = MUX_LIST < L, 32*Q, LOG2L> ( lb, a_l_llr);

			sc_bigint< (32*Q) > lb1 = Function_G< 32 ,Q> (i_la,i_lb,sa1[j]);

			Second_R_P1_in[j].llr = lb1;
			Second_R_P1_in[j].metric = First_R_P1_out[j].metric;
			Second_R_P1_in[j].path = First_R_P1_out[j].path;

			//// List Stack management  ////
			LIST_STACK[j] = write_stack <LOG2_L, _DEPTH > ( LIST_STACK[j] , (sc_biguint<LOG2_L>) j );
		}

		PS_struct< 32 ,Q,LOG2L> Second_R_P1_out[L];
// Second P1
		List_R_P32 <L,Q,LOG2L,D> (Second_R_P1_in, Second_R_P1_out, (sc_biguint< 32 >)fb.range(63, 32), LIST_STACK);

// XOR
		for(COUNTER j = 0; j < L_SIZE; j++){
		#pragma HLS UNROLL
			sc_biguint<LOG2_L> a_list = read_stack< LOG2_L, _DEPTH > (LIST_STACK[j], 1);
			sc_biguint< 32 > i_sa = (sc_biguint< 32 >) MUX_LIST <L, 32, LOG2L> ( sa1, a_list);
			sc_biguint< 32 > i_sb = Second_R_P1_out[j].bit;

			sc_biguint< 32 > sa = VECTOR_XOR< 32 >(i_sa,i_sb);
			sc_biguint< 32 > sb = i_sb;

			output[j].bit = (sb , sa);
			output[j].metric = Second_R_P1_out[j].metric;
			output[j].path = Second_R_P1_out[j].path;

			//// List Stack management  ////
			LIST_STACK[j] = pop_stack <LOG2_L, _DEPTH >( LIST_STACK[j] , (sc_biguint<LOG2_L>) 0 );
		}
}


template <int L, int Q, int LOG2L, int D>
void List_R_P128 ( LLR_struct<128,Q,LOG2L> input[L], PS_struct<128,Q,LOG2L> output[L], sc_biguint<128> fb, sc_biguint<LOG2L*D> LIST_STACK[L] )
{
#if PARALLEL_TP == 0
	#pragma HLS INLINE off
#else
	#pragma HLS INLINE
#endif
#pragma HLS ARRAY_PARTITION variable=input complete dim=1
#pragma HLS ARRAY_PARTITION variable=output complete dim=1

		sc_bigint< (64*Q) > la[L];
		sc_bigint< (64*Q) > lb[L];
		LLR_struct< 64 ,Q,LOG2L> First_R_P1_in[L];
// F
		for(COUNTER j = 0; j < L_SIZE; j++){
		#pragma HLS UNROLL
			la[j] = (input[j].llr).range( (64*Q) -1,0);
			lb[j] = (input[j].llr).range(2* (64*Q) -1, (64*Q) );

			sc_bigint< (64*Q) > la1 = Function_F< 64 ,Q> (la[j],lb[j]);

			First_R_P1_in[j].llr = la1;
			First_R_P1_in[j].metric = input[j].metric;
			First_R_P1_in[j].path = input[j].path;

			//// List Stack management  ////
			LIST_STACK[j] = push_stack <LOG2_L, _DEPTH >( LIST_STACK[j] , (sc_biguint<LOG2_L>) j );
		}

		PS_struct< 64 ,Q,LOG2L> First_R_P1_out[L];
// First P1
		List_R_P64 <L,Q,LOG2L,D> (First_R_P1_in, First_R_P1_out, (sc_biguint< 64 >)fb.range(63, 0), LIST_STACK);

		sc_bigint< 64 > sa1[L];
		LLR_struct< 64 ,Q,LOG2L> Second_R_P1_in[L];
// G
		for(COUNTER j = 0; j < L_SIZE; j++){
		#pragma HLS UNROLL
			sa1[j] = First_R_P1_out[j].bit;

			sc_biguint<LOG2_L> a_l_llr = read_stack< LOG2_L, _DEPTH > (LIST_STACK[j], 2);
			sc_bigint< (64*Q) > i_la = MUX_LIST < L, 64*Q, LOG2L> ( la, a_l_llr);
			sc_bigint< (64*Q) > i_lb = MUX_LIST < L, 64*Q, LOG2L> ( lb, a_l_llr);

			sc_bigint< (64*Q) > lb1 = Function_G< 64 ,Q> (i_la,i_lb,sa1[j]);

			Second_R_P1_in[j].llr = lb1;
			Second_R_P1_in[j].metric = First_R_P1_out[j].metric;
			Second_R_P1_in[j].path = First_R_P1_out[j].path;

			//// List Stack management  ////
			LIST_STACK[j] = write_stack <LOG2_L, _DEPTH > ( LIST_STACK[j] , (sc_biguint<LOG2_L>) j );
		}

		PS_struct< 64 ,Q,LOG2L> Second_R_P1_out[L];
// Second P1
		List_R_P64 <L,Q,LOG2L,D> (Second_R_P1_in, Second_R_P1_out, (sc_biguint< 64 >)fb.range(127, 64), LIST_STACK);

// XOR
		for(COUNTER j = 0; j < L_SIZE; j++){
		#pragma HLS UNROLL
			sc_biguint<LOG2_L> a_list = read_stack< LOG2_L, _DEPTH > (LIST_STACK[j], 1);
			sc_biguint< 64 > i_sa = (sc_biguint< 64 >) MUX_LIST <L, 64, LOG2L> ( sa1, a_list);
			sc_biguint< 64 > i_sb = Second_R_P1_out[j].bit;

			sc_biguint< 64 > sa = VECTOR_XOR< 64 >(i_sa,i_sb);
			sc_biguint< 64 > sb = i_sb;

			output[j].bit = (sb , sa);
			output[j].metric = Second_R_P1_out[j].metric;
			output[j].path = Second_R_P1_out[j].path;

			//// List Stack management  ////
			LIST_STACK[j] = pop_stack <LOG2_L, _DEPTH >( LIST_STACK[j] , (sc_biguint<LOG2_L>) 0 );
		}
}


template <int L, int Q, int LOG2L, int D>
void List_R_P256 ( LLR_struct<256,Q,LOG2L> input[L], PS_struct<256,Q,LOG2L> output[L], sc_biguint<256> fb, sc_biguint<LOG2L*D> LIST_STACK[L] )
{
#if PARALLEL_TP == 0
	#pragma HLS INLINE off
#else
	#pragma HLS INLINE
#endif
#pragma HLS ARRAY_PARTITION variable=input complete dim=1
#pragma HLS ARRAY_PARTITION variable=output complete dim=1

		sc_bigint< (128*Q) > la[L];
		sc_bigint< (128*Q) > lb[L];
		LLR_struct< 128 ,Q,LOG2L> First_R_P1_in[L];
// F
		for(COUNTER j = 0; j < L_SIZE; j++){
		#pragma HLS UNROLL
			la[j] = (input[j].llr).range( (128*Q) -1,0);
			lb[j] = (input[j].llr).range(2* (128*Q) -1, (128*Q) );

			sc_bigint< (128*Q) > la1 = Function_F< 128 ,Q> (la[j],lb[j]);

			First_R_P1_in[j].llr = la1;
			First_R_P1_in[j].metric = input[j].metric;
			First_R_P1_in[j].path = input[j].path;

			//// List Stack management  ////
			LIST_STACK[j] = push_stack <LOG2_L, _DEPTH >( LIST_STACK[j] , (sc_biguint<LOG2_L>) j );
		}

		PS_struct< 128 ,Q,LOG2L> First_R_P1_out[L];
// First P1
		List_R_P128 <L,Q,LOG2L,D> (First_R_P1_in, First_R_P1_out, (sc_biguint< 128 >)fb.range(127, 0), LIST_STACK);

		sc_bigint< 128 > sa1[L];
		LLR_struct< 128 ,Q,LOG2L> Second_R_P1_in[L];
// G
		for(COUNTER j = 0; j < L_SIZE; j++){
		#pragma HLS UNROLL
			sa1[j] = First_R_P1_out[j].bit;

			sc_biguint<LOG2_L> a_l_llr = read_stack< LOG2_L, _DEPTH > (LIST_STACK[j], 2);
			sc_bigint< (128*Q) > i_la = MUX_LIST < L, 128*Q, LOG2L> ( la, a_l_llr);
			sc_bigint< (128*Q) > i_lb = MUX_LIST < L, 128*Q, LOG2L> ( lb, a_l_llr);

			sc_bigint< (128*Q) > lb1 = Function_G< 128 ,Q> (i_la,i_lb,sa1[j]);

			Second_R_P1_in[j].llr = lb1;
			Second_R_P1_in[j].metric = First_R_P1_out[j].metric;
			Second_R_P1_in[j].path = First_R_P1_out[j].path;

			//// List Stack management  ////
			LIST_STACK[j] = write_stack <LOG2_L, _DEPTH > ( LIST_STACK[j] , (sc_biguint<LOG2_L>) j );
		}

		PS_struct< 128 ,Q,LOG2L> Second_R_P1_out[L];
// Second P1
		List_R_P128 <L,Q,LOG2L,D> (Second_R_P1_in, Second_R_P1_out, (sc_biguint< 128 >)fb.range(255, 128), LIST_STACK);

// XOR
		for(COUNTER j = 0; j < L_SIZE; j++){
		#pragma HLS UNROLL
			sc_biguint<LOG2_L> a_list = read_stack< LOG2_L, _DEPTH > (LIST_STACK[j], 1);
			sc_biguint< 128 > i_sa = (sc_biguint< 128 >) MUX_LIST <L, 128, LOG2L> ( sa1, a_list);
			sc_biguint< 128 > i_sb = Second_R_P1_out[j].bit;

			sc_biguint< 128 > sa = VECTOR_XOR< 128 >(i_sa,i_sb);
			sc_biguint< 128 > sb = i_sb;

			output[j].bit = (sb , sa);
			output[j].metric = Second_R_P1_out[j].metric;
			output[j].path = Second_R_P1_out[j].path;

			//// List Stack management  ////
			LIST_STACK[j] = pop_stack <LOG2_L, _DEPTH >( LIST_STACK[j] , (sc_biguint<LOG2_L>) 0 );
		}
}

#endif
