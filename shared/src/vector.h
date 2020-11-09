#ifndef VECTOR
#define VECTOR

#include "scalar.h"

template <int P, int Q>
inline void SHOW( sc_bigint<P*Q> a )
{
	for(int i=0; i<P; i+=1)
	{
		sc_bigint<Q> oa = a.range(Q * (i+1) - 1, Q * i);
		printf("%3d ", oa.to_int());
	}
}

template <int P>
inline void SHOW_BITS( sc_biguint<P> a )
{
	for(int i=0; i<P; i+=1)
	{
		sc_biguint<1> oa = (sc_biguint<1>) a[i];
		printf("%3d ", oa.to_int());
	}
}

template <int P, int Q>
inline void SHOW_LLRS( sc_bigint<P*Q> a )
{
#if defined CA2
	return SHOW<P,Q> (a);
#elif defined SIGMAG
	return SHOW_SM<P,Q> ( (sc_biguint<P*Q>)a );
#endif
}

//
////////////////////////////////////////////////////////////////////////////////
//

template <int P, int Q>
inline sc_bigint<P*Q> VECTOR_MIN(sc_bigint<P*Q> a, sc_bigint<P*Q> b)
{
#pragma HLS INLINE
	sc_bigint<P*Q> c;
	for(int i=0; i<P; i+=1)
	{
	#pragma HLS UNROLL
		sc_bigint<Q>  oa = a.range(Q * (i + 1) - 1, Q * i);
		sc_bigint<Q>  ob = b.range(Q * (i + 1) - 1, Q * i);
		c.range(Q * (i+1) - 1, Q * i) = qmin<Q>( oa, ob );
	}
	return c;
}

//
////////////////////////////////////////////////////////////////////////////////
//

template <int P, int Q>
inline sc_biguint<P> VECTOR_IS_MIN(sc_bigint<P*Q> a, sc_bigint<P*Q> b) // comparator : if a < b return 1
{
#pragma HLS INLINE
	sc_biguint<P> c;
	for(int i=0; i<P; i+=1)
	{
	#pragma HLS UNROLL
		sc_bigint<Q>  oa = a.range(Q * (i + 1) - 1, Q * i);
		sc_bigint<Q>  ob = b.range(Q * (i + 1) - 1, Q * i);
		c[i] = ( oa < ob ) ? 1 : 0;
	}
	return c;
}

//
////////////////////////////////////////////////////////////////////////////////
//

template <int P, int Q>
inline sc_bigint<P*Q> VECTOR_SAT(sc_bigint<P*(Q+1)> a)
{
#pragma HLS INLINE
	sc_bigint<P*Q> c;
	for(int i=0; i<P; i+=1)
	{
	#pragma HLS UNROLL
		sc_bigint<Q+1>  oa = a.range((Q+1) * (i + 1) - 1, (Q+1) * i);
		c.range(Q * (i+1) - 1, Q * i) = qsat<Q>( oa );
	}
	return c;
}

//
////////////////////////////////////////////////////////////////////////////////
//

template <int P>
inline sc_biguint<P> VECTOR_INIT_1()
{
#pragma HLS INLINE
	sc_biguint<P> c;
	for(int i=0; i<P; i+=1)
	{
 	#pragma HLS UNROLL
		c.range(i, i) = (sc_biguint<1>)1;
	}
	return c;
}

template <int P, int Q>
inline sc_biguint<P> VECTOR_SIGN(sc_bigint<P*Q> a)
{
#pragma HLS INLINE
	sc_biguint<P> c;
	for(int i=0; i<P; i+=1)
	{
 	#pragma HLS UNROLL
		sc_bigint<Q> llr = a.range(Q * (i + 1) - 1, Q * i);
		c.range(i, i) = qsign<Q>( llr );
	}
	return c;
}

template <int P, int Q>
inline sc_bigint<P*Q> VECTOR_SIGN(sc_bigint<P*Q> a, sc_biguint<P> b)
{
#pragma HLS INLINE
	sc_bigint<P*Q> c;
	for(int i=0; i<P; i+=1)
	{
 	#pragma HLS UNROLL
		sc_bigint<Q> oa = a.range(Q * (i + 1) - 1, Q * i);
		sc_biguint<1> ob = (sc_biguint<1>)b[i];
		c.range(Q * (i+1) - 1, Q * i) = qsign<Q>( oa, ob );
	}
	return c;
}

//
////////////////////////////////////////////////////////////////////////////////
//

template <int P, int Q>
inline sc_bigint<P*Q> VECTOR_ABS(sc_bigint<P*Q> a)
{
#pragma HLS INLINE
	sc_bigint<P*Q> c;
	for(int i=0; i<P; i+=1)
	{
	#pragma HLS UNROLL
		sc_bigint<Q> oa = a.range(Q * (i + 1) - 1, Q * i);
		c.range(Q * (i + 1) - 1, Q * i) = qabs<Q>( oa );
	}
	return c;
}

//
////////////////////////////////////////////////////////////////////////////////
//

template <int P, int Q>
inline sc_bigint<P*Q> VECTOR_ADD(sc_bigint<P*Q> a, sc_bigint<P*Q> b)
{
#pragma HLS INLINE
	sc_bigint<P*Q> c;
	for(int i=0; i<P; i+=1)
	{
	#pragma HLS UNROLL
		sc_bigint<Q> oa = a.range(Q * (i + 1) - 1, Q * i);
		sc_bigint<Q> ob = b.range(Q * (i + 1) - 1, Q * i);
		c.range(Q * (i + 1) - 1, Q * i) = qadd<Q>( oa , ob );
	}
	return c;
}

//
////////////////////////////////////////////////////////////////////////////////
//

template <int P, int Q>
inline sc_bigint<P*(Q+1)> VECTOR_ADD_NOSAT(sc_bigint<P*Q> a, sc_bigint<P*Q> b)
{
#pragma HLS INLINE
	sc_bigint<P*(Q+1)> c;
	for(int i=0; i<P; i+=1)
	{
	#pragma HLS UNROLL
		sc_bigint<Q> oa = a.range(Q * (i + 1) - 1, Q * i);
		sc_bigint<Q> ob = b.range(Q * (i + 1) - 1, Q * i);
		c.range((Q+1) * (i + 1) - 1, (Q+1) * i) = qadd_nosat<Q>( oa , ob );
	}
	return c;
}

//
////////////////////////////////////////////////////////////////////////////////
//

template <int P, int Q>
inline sc_bigint<P*Q> VECTOR_SUB(sc_bigint<P*Q> a, sc_bigint<P*Q> b)
{
#pragma HLS INLINE
	sc_bigint<P*Q> c;
	for(int i=0; i<P; i+=1)
	{
	#pragma HLS UNROLL
		sc_bigint<Q> oa = a.range(Q * (i + 1) - 1, Q * i);
		sc_bigint<Q> ob = b.range(Q * (i + 1) - 1, Q * i);
		c.range(Q * (i + 1) - 1, Q * i) = qsub<Q>( oa , ob );
	}
	return c;
}

//
////////////////////////////////////////////////////////////////////////////////
//

template <int P, int Q>
inline sc_bigint<P*(Q+1)> VECTOR_SUB_NOSAT(sc_bigint<P*Q> a, sc_bigint<P*Q> b)
{
#pragma HLS INLINE
	sc_bigint<P*(Q+1)> c;
	for(int i=0; i<P; i+=1)
	{
	#pragma HLS UNROLL
		sc_bigint<Q> oa = a.range(Q * (i + 1) - 1, Q * i);
		sc_bigint<Q> ob = b.range(Q * (i + 1) - 1, Q * i);
		c.range((Q+1) * (i + 1) - 1, (Q+1) * i) = qsub_nosat<Q>( oa , ob );
	}
	return c;
}

//
////////////////////////////////////////////////////////////////////////////////
//

template <int P, int Q>
inline sc_bigint<P*Q> VECTOR_MUX(sc_bigint<P*Q> a, sc_bigint<P*Q> b, sc_biguint<P> s) //return a when s = 1, b otherwise
{
#pragma HLS INLINE
	sc_bigint<P*Q> c;
	for(int i=0; i<P; i+=1)
	{
	#pragma HLS UNROLL
		sc_bigint<Q> oa = a.range(Q * (i + 1) - 1, Q * i);
		sc_bigint<Q> ob = b.range(Q * (i + 1) - 1, Q * i);
		sc_biguint<1> os = (sc_biguint<1>) s[i];
		c.range(Q * (i + 1) - 1, Q * i) = (sc_uint<1>) os ? oa : ob;
	}
	return c;
}


//*************************************************************************//
//////////////////////////Signe et Magnetude/////////////////////////////////
//*************************************************************************//

template <int P, int Q>
inline void SHOW_SM( sc_biguint<P*Q> a )
{
	for(int i=0; i<P; i+=1)
	{
		sc_bigint<Q> oa = a.range(Q * (i+1) - 1, Q * i);
		sc_bigint<Q> b = (sc_bigint<Q>) qconv_format<Q>(oa);
		printf("%3d ", b.to_int());
	}
}

//
////////////////////////////////////////////////////////////////////////////////
//

template <int P, int Q>
inline sc_biguint<P*Q> VECTOR_MIN_SM(sc_biguint<P*Q> a, sc_biguint<P*Q> b)
{
#pragma HLS INLINE
	sc_biguint<P*Q> c;
	for(int i=0; i<P; i+=1)
	{
	#pragma HLS UNROLL
		sc_biguint<Q>  oa = a.range(Q * (i + 1) - 1, Q * i);
		sc_biguint<Q>  ob = b.range(Q * (i + 1) - 1, Q * i);
		c.range(Q * (i+1) - 1, Q * i) = qmin_sm<Q>( oa, ob );
	}
	return c;
}

//
////////////////////////////////////////////////////////////////////////////////
//

template <int P, int Q>
inline sc_biguint<P*Q> VECTOR_SAT_SM(sc_biguint<P*(Q+1)> a)
{
#pragma HLS INLINE
	sc_biguint<P*Q> c;
	for(int i=0; i<P; i+=1)
	{
	#pragma HLS UNROLL
		sc_biguint<Q+1>  oa = a.range((Q+1) * (i + 1) - 1, (Q+1) * i);
		c.range(Q * (i+1) - 1, Q * i) = qsat_sm<Q>( oa);
	}
	return c;
}

//
////////////////////////////////////////////////////////////////////////////////
//

template <int P, int Q>
inline sc_biguint<P*Q> VECTOR_TRONC_SM(sc_biguint<P*(Q+1)> a)
{
#pragma HLS INLINE
	sc_biguint<P*Q> c;
	for(int i=0; i<P; i+=1)
	{
	#pragma HLS UNROLL
		sc_biguint<Q+1>  oa = a.range((Q+1) * (i + 1) - 1, (Q+1) * i);
		c.range(Q * (i+1) - 1, Q * i) = oa.range(Q-2 , 0);
	}
	return c;
}

//
////////////////////////////////////////////////////////////////////////////////
//

template <int P, int Q>
inline sc_biguint<P> VECTOR_SIGN_SM(sc_biguint<P*Q> a)
{
#pragma HLS INLINE
	sc_biguint<P> c;
	for(int i=0; i<P; i+=1)
	{
	#pragma HLS UNROLL
		sc_biguint<Q>  oa = a.range(Q * (i + 1) - 1, Q * i);
		c.range(i, i) = qsign_sm<Q>( oa);
	}
	return c;
}

//
////////////////////////////////////////////////////////////////////////////////
//

template <int P, int Q>
inline sc_biguint<P*(Q-1)> VECTOR_ABS_SM(sc_biguint<P*Q> a)
{
#pragma HLS INLINE
	sc_biguint<P*(Q-1)> c;
	for(int i=0; i<P; i+=1)
	{
	#pragma HLS UNROLL
		sc_biguint<Q>  oa = a.range(Q * (i + 1) - 1, Q * i);
		c.range((Q-1) * (i+1) - 1, (Q-1) * i) = qabs_sm<Q>( oa );
	}
	return c;
}

//
////////////////////////////////////////////////////////////////////////////////
//

template <int P, int Q>
inline sc_biguint<P*Q> VECTOR_ADD_CARRY_SM(sc_biguint<P*Q> a, sc_biguint<P*Q> b, sc_biguint<P> carry)
{
#pragma HLS INLINE
	sc_biguint<P*Q> c;
	for(int i=0; i<P; i+=1)
	{
	#pragma HLS UNROLL
		sc_biguint<Q>  oa = a.range(Q * (i + 1) - 1, Q * i);
		sc_biguint<Q>  ob = b.range(Q * (i + 1) - 1, Q * i);
		sc_biguint<1> ocarry = (sc_biguint<1>)carry[i];
		c.range(Q * (i+1) - 1, Q * i) = qadd_carry_sm<Q>( oa, ob, ocarry );
	}
	return c;
}

//
////////////////////////////////////////////////////////////////////////////////
//

template <int P, int Q>
inline sc_biguint<P*(Q+1)> VECTOR_FULL_ADDER_SM(sc_biguint<P*Q> a, sc_biguint<P*Q> b)
{
#pragma HLS INLINE
	sc_biguint<P*(Q+1)> c;
	for(int i=0; i<P; i+=1)
	{
	#pragma HLS UNROLL
		sc_biguint<Q>  oa = a.range(Q * (i + 1) - 1, Q * i);
		sc_biguint<Q>  ob = b.range(Q * (i + 1) - 1, Q * i);
		c.range( (Q+1) * (i+1) - 1, (Q+1) * i) = qfull_adder_sm<Q>( oa, ob);
	}
	return c;
}

//
////////////////////////////////////////////////////////////////////////////////
//

template <int P, int Q>
inline sc_biguint<P*Q> VECTOR_FULL_ADDER_SAT_SM(sc_biguint<P*Q> a, sc_biguint<P*Q> b)
{
#pragma HLS INLINE
	sc_biguint<P*(Q+1)> c;
	for(int i=0; i<P; i+=1)
	{
	#pragma HLS UNROLL
		sc_biguint<Q>  oa = a.range(Q * (i + 1) - 1, Q * i);
		sc_biguint<Q>  ob = b.range(Q * (i + 1) - 1, Q * i);
		c.range( Q * (i+1) - 1, Q * i) = qfull_adder_sat_sm<Q>( oa, ob);
	}
	return c;
}

//
////////////////////////////////////////////////////////////////////////////////
//

template <int P, int Q>
inline sc_biguint<P*(Q+1)> VECTOR_FULL_ADD_SUB_SM(sc_biguint<P*Q> a, sc_biguint<P*Q> b, sc_biguint<P> s)  // s = 1 : b-a, s=0 : b+a
{
#pragma HLS INLINE
	sc_biguint<P*(Q+1)> c;
	for(int i=0; i<P; i+=1)
	{
	#pragma HLS UNROLL
		sc_biguint<Q>  oa = a.range(Q * (i + 1) - 1, Q * i);
		sc_biguint<Q>  ob = b.range(Q * (i + 1) - 1, Q * i);
		sc_biguint<1>  os = (sc_biguint<1>)s[i];
		c.range( (Q+1) * (i+1) - 1, (Q+1) * i) = qfull_add_sub_sm<Q>( oa, ob, os);
	}
	return c;
}

//
////////////////////////////////////////////////////////////////////////////////
//

template <int P, int Q>
inline sc_biguint<P*Q> VECTOR_NOT_SM(sc_biguint<P*Q> a)
{
#pragma HLS INLINE
	sc_biguint<P*Q> c;
	for(int i=0; i<P; i+=1)
	{
	#pragma HLS UNROLL
		sc_biguint<Q>  oa = a.range(Q * (i + 1) - 1, Q * i);
		c.range(Q * (i+1) - 1, Q * i) = qnot_sm<Q>( oa );
	}
	return c;
}

//
////////////////////////////////////////////////////////////////////////////////
//

template <int P, int Q>
inline sc_biguint<P*Q> VECTOR_MUX_SM(sc_biguint<P*Q> a, sc_biguint<P*Q> b, sc_biguint<P> s) //return a when s = 1, b otherwise
{
#pragma HLS INLINE
	sc_biguint<P*Q> c;
	for(int i=0; i<P; i+=1)
	{
	#pragma HLS UNROLL
		sc_biguint<Q> oa = a.range(Q * (i + 1) - 1, Q * i);
		sc_biguint<Q> ob = b.range(Q * (i + 1) - 1, Q * i);
		sc_biguint<1> os = (sc_biguint<1>) s[i];
		c.range(Q * (i + 1) - 1, Q * i) = (sc_uint<1>)os ? oa : ob;
	}
	return c;
}

//
////////////////////////////////////////////////////////////////////////////////
//

template <int P, int Q>
inline sc_biguint<P*Q> VECTOR_CONCAT_SM(sc_biguint<P> sign, sc_biguint<P*(Q-1)> mag)
{
#pragma HLS INLINE
	sc_biguint<P*Q> c;
	for(int i=0; i<P; i+=1)
	{
	#pragma HLS UNROLL
		sc_biguint<1> os = (sc_biguint<1>) sign[i];
		sc_biguint<(Q-1)> om = mag.range((Q-1) * (i + 1) - 1, (Q-1) * i);
		c.range(Q * (i+1) - 1, Q * i) = ( os, om );
	}
	return c;
}

//
////////////////////////////////////////////////////////////////////////////////
//

template <int P, int Q>
inline sc_biguint<P*Q> VECTOR_INV_SM(sc_biguint<P*Q> a)  // A = (not A) + 1
{
#pragma HLS INLINE
	sc_biguint<P*Q> c;
	for(int i=0; i<P; i+=1)
	{
	#pragma HLS UNROLL
		sc_biguint<Q>  oa = a.range(Q * (i + 1) - 1, Q * i);
		sc_biguint<Q>  n_oa = (~oa);
		c.range(Q * (i+1) - 1, Q * i) = n_oa + (sc_biguint<1>)(1);
	}
	return c;
}

//
////////////////////////////////////////////////////////////////////////////////
//

template <int P, int Q>
inline sc_biguint<P> VECTOR_IS_MIN_SM(sc_biguint<P*Q> a, sc_biguint<P*Q> b) // comparator : if a < b return 1
{
#pragma HLS INLINE
	sc_biguint<P> c;
	for(int i=0; i<P; i+=1)
	{
	#pragma HLS UNROLL
		sc_biguint<Q>  oa = a.range(Q * (i + 1) - 1, Q * i);
		sc_biguint<Q>  ob = b.range(Q * (i + 1) - 1, Q * i);
		c[i] = ( oa < ob ) ? 1 : 0;
	}
	return c;
}

//
////////////////////////////////////////////////////////////////////////////////
//

template <int P, int Q>
inline sc_biguint<P*(Q+1)> VECTOR_EXT_SM(sc_biguint<P*Q> a)
{
#pragma HLS INLINE
	sc_biguint<P*(Q+1)> c;
	for(int i=0; i<P; i+=1)
	{
	#pragma HLS UNROLL
		sc_biguint<Q>  oa = a.range(Q * (i + 1) - 1, Q * i);
		c.range((Q+1) * (i+1) - 1, (Q+1) * i) = (sc_biguint<Q+1>) oa;
	}
	return c;
}

//*************************************************************************//
/////////////////////////////// BITS ////////////////////////////////////////
//*************************************************************************//

template <int P>
inline sc_biguint<P> VECTOR_XOR(sc_biguint<P> a, sc_biguint<P> b)
{
#pragma HLS INLINE
	sc_biguint<P> c;
	for(int i=0; i<P; i+=1)
	{
	#pragma HLS UNROLL
		sc_biguint<1> oa = (sc_biguint<1>)a[i];
		sc_biguint<1> ob = (sc_biguint<1>)b[i];
		c.range(i,i) = qxor( oa , ob );
	}
	return c;
}

template <int P>
inline sc_uint<P> VECTOR_OR(sc_uint<P> a, sc_uint<P> b)
{
#pragma HLS INLINE
	sc_uint<P> c;
	for(int i=0; i<P; i+=1)
	{
	#pragma HLS UNROLL
		sc_uint<1> oa = (sc_uint<1>)a[i];
		sc_uint<1> ob = (sc_uint<1>)b[i];
		c.range(i,i) = oa | ob;
	}
	return c;
}

template <int P>
inline sc_biguint<P> VECTOR_AND(sc_biguint<P> a, sc_biguint<P> b)
{
#pragma HLS INLINE
	sc_biguint<P> c;
	for(int i=0; i<P; i+=1)
	{
	#pragma HLS UNROLL
		sc_biguint<1> oa = (sc_biguint<1>)a[i];
		sc_biguint<1> ob = (sc_biguint<1>)b[i];
		c.range(i,i) = qand( oa , ob );
	}
	return c;
}

template <int P>
inline sc_uint<P> VECTOR_AND(sc_uint<P> a, sc_uint<P> b)
{
#pragma HLS INLINE
	sc_uint<P> c;
	for(int i=0; i<P; i+=1)
	{
	#pragma HLS UNROLL
		sc_uint<1> oa = (sc_uint<1>)a[i];
		sc_uint<1> ob = (sc_uint<1>)b[i];
		c.range(i,i) = oa & ob;
	}
	return c;
}


//*************************************************************************//
///////////////////////////////    LIST MUX     /////////////////////////////
//*************************************************************************//

template <int L, int Q, int LOG2L>
inline sc_bigint<Q> LIST_MUX2 ( sc_bigint<Q> tab[L], sc_uint<LOG2L> sel)
{
#pragma HLS INLINE
#pragma HLS ARRAY_PARTITION variable=tab complete dim=1
	if( sel == 0 ){
		return tab[0];}
	else{
		return tab[1];}
}

template <int L, int Q, int LOG2L>
inline sc_bigint<Q> LIST_MUX4 ( sc_bigint<Q> tab[L], sc_uint<LOG2L> sel)
{
#pragma HLS INLINE
#pragma HLS ARRAY_PARTITION variable=tab complete dim=1
	if( sel == 0 ){
		return tab[0];}
	else if( sel == 1 ){
		return tab[1];}
	else if( sel == 2 ){
		return tab[2];}
	else{
		return tab[3];}
}

template <int L, int Q, int LOG2L>
inline sc_bigint<Q> LIST_MUX8 ( sc_bigint<Q> tab[L], sc_uint<LOG2L> sel)
{
#pragma HLS INLINE
#pragma HLS ARRAY_PARTITION variable=tab complete dim=1
	if( sel == 0 ){
		return tab[0];}
	else if( sel == 1 ){
		return tab[1];}
	else if( sel == 2 ){
		return tab[2];}
	else if( sel == 3 ){
		return tab[3];}
	else if( sel == 4 ){
		return tab[4];}
	else if( sel == 5 ){
		return tab[5];}
	else if( sel == 6 ){
		return tab[6];}
	else{
		return tab[7];
	}
}

template <int L, int Q, int LOG2L>
inline sc_bigint<Q> LIST_MUX16 ( sc_bigint<Q> tab[L], sc_uint<LOG2L> sel)
{
#pragma HLS INLINE
#pragma HLS ARRAY_PARTITION variable=tab complete dim=1
	if( sel == 0 ){
		return tab[0];}
	else if( sel == 1 ){
		return tab[1];}
	else if( sel == 2 ){
		return tab[2];}
	else if( sel == 3 ){
		return tab[3];}
	else if( sel == 4 ){
		return tab[4];}
	else if( sel == 5 ){
		return tab[5];}
	else if( sel == 6 ){
		return tab[6];}
	else if( sel == 7 ){
		return tab[7];}
	else if( sel == 8 ){
		return tab[8];}
	else if( sel == 9 ){
		return tab[9];}
	else if( sel == 10 ){
		return tab[10];}
	else if( sel == 11 ){
		return tab[11];}
	else if( sel == 12 ){
		return tab[12];}
	else if( sel == 13 ){
		return tab[13];}
	else if( sel == 14 ){
		return tab[14];}
	else{
		return tab[15];
	}
}

template <int L, int Q, int LOG2L>
inline sc_bigint<Q> LIST_MUX32 ( sc_bigint<Q> tab[L], sc_uint<LOG2L> sel)
{
#pragma HLS INLINE
#pragma HLS ARRAY_PARTITION variable=tab complete dim=1
	if( sel == 0 ){
		return tab[0];}
	else if( sel == 1 ){
		return tab[1];}
	else if( sel == 2 ){
		return tab[2];}
	else if( sel == 3 ){
		return tab[3];}
	else if( sel == 4 ){
		return tab[4];}
	else if( sel == 5 ){
		return tab[5];}
	else if( sel == 6 ){
		return tab[6];}
	else if( sel == 7 ){
		return tab[7];}
	else if( sel == 8 ){
		return tab[8];}
	else if( sel == 9 ){
		return tab[9];}
	else if( sel == 10 ){
		return tab[10];}
	else if( sel == 11 ){
		return tab[11];}
	else if( sel == 12 ){
		return tab[12];}
	else if( sel == 13 ){
		return tab[13];}
	else if( sel == 14 ){
		return tab[14];}
	else if( sel == 15 ){
		return tab[15];}
	else if( sel == 16 ){
		return tab[16];}
	else if( sel == 17 ){
		return tab[17];}
	else if( sel == 18 ){
		return tab[18];}
	else if( sel == 19 ){
		return tab[19];}
	else if( sel == 20 ){
		return tab[20];}
	else if( sel == 21 ){
		return tab[21];}
	else if( sel == 22 ){
		return tab[22];}
	else if( sel == 23 ){
		return tab[23];}
	else if( sel == 24 ){
		return tab[24];}
	else if( sel == 25 ){
		return tab[25];}
	else if( sel == 26 ){
		return tab[26];}
	else if( sel == 27 ){
		return tab[27];}
	else if( sel == 28 ){
		return tab[28];}
	else if( sel == 29 ){
		return tab[29];}
	else if( sel == 30 ){
		return tab[30];}
	else{
		return tab[31];
	}
}

template <int L, int Q, int LOG2L>
inline sc_bigint<Q> LIST_MUX64 ( sc_bigint<Q> tab[L], sc_uint<LOG2L> sel)
{
#pragma HLS INLINE
#pragma HLS ARRAY_PARTITION variable=tab complete dim=1
	if( sel == 0 ){
		return tab[0];}
	else if( sel == 1 ){
		return tab[1];}
	else if( sel == 2 ){
		return tab[2];}
	else if( sel == 3 ){
		return tab[3];}
	else if( sel == 4 ){
		return tab[4];}
	else if( sel == 5 ){
		return tab[5];}
	else if( sel == 6 ){
		return tab[6];}
	else if( sel == 7 ){
		return tab[7];}
	else if( sel == 8 ){
		return tab[8];}
	else if( sel == 9 ){
		return tab[9];}
	else if( sel == 10 ){
		return tab[10];}
	else if( sel == 11 ){
		return tab[11];}
	else if( sel == 12 ){
		return tab[12];}
	else if( sel == 13 ){
		return tab[13];}
	else if( sel == 14 ){
		return tab[14];}
	else if( sel == 15 ){
		return tab[15];}
	else if( sel == 16 ){
		return tab[16];}
	else if( sel == 17 ){
		return tab[17];}
	else if( sel == 18 ){
		return tab[18];}
	else if( sel == 19 ){
		return tab[19];}
	else if( sel == 20 ){
		return tab[20];}
	else if( sel == 21 ){
		return tab[21];}
	else if( sel == 22 ){
		return tab[22];}
	else if( sel == 23 ){
		return tab[23];}
	else if( sel == 24 ){
		return tab[24];}
	else if( sel == 25 ){
		return tab[25];}
	else if( sel == 26 ){
		return tab[26];}
	else if( sel == 27 ){
		return tab[27];}
	else if( sel == 28 ){
		return tab[28];}
	else if( sel == 29 ){
		return tab[29];}
	else if( sel == 30 ){
		return tab[30];}
	else if( sel == 31 ){
		return tab[31];}
	else if( sel == 32 ){
		return tab[32];}
	else if( sel == 33 ){
		return tab[33];}
	else if( sel == 34 ){
		return tab[34];}
	else if( sel == 35 ){
		return tab[35];}
	else if( sel == 36 ){
		return tab[36];}
	else if( sel == 37 ){
		return tab[37];}
	else if( sel == 38 ){
		return tab[38];}
	else if( sel == 39 ){
		return tab[39];}
	else if( sel == 40 ){
		return tab[40];}
	else if( sel == 41 ){
		return tab[41];}
	else if( sel == 42 ){
		return tab[42];}
	else if( sel == 43 ){
		return tab[43];}
	else if( sel == 44 ){
		return tab[44];}
	else if( sel == 45 ){
		return tab[45];}
	else if( sel == 46 ){
		return tab[46];}
	else if( sel == 47 ){
		return tab[47];}
	else if( sel == 48 ){
		return tab[48];}
	else if( sel == 49 ){
		return tab[49];}
	else if( sel == 50 ){
		return tab[50];}
	else if( sel == 51 ){
		return tab[51];}
	else if( sel == 52 ){
		return tab[52];}
	else if( sel == 53 ){
		return tab[53];}
	else if( sel == 54 ){
		return tab[54];}
	else if( sel == 55 ){
		return tab[55];}
	else if( sel == 56 ){
		return tab[56];}
	else if( sel == 57 ){
		return tab[57];}
	else if( sel == 58 ){
		return tab[58];}
	else if( sel == 59 ){
		return tab[59];}
	else if( sel == 60 ){
		return tab[60];}
	else if( sel == 61 ){
		return tab[61];}
	else if( sel == 62 ){
		return tab[62];}
	else{
		return tab[63];
	}
}

//*************************************************************************//
////////////////////////// CLASS MUX ////////////////////////////////////////
//*************************************************************************//

template <class T, int L, int Q, int LOG2L>
inline T LIST_MUX2 ( T tab[L], sc_uint<LOG2L> sel)
{
#pragma HLS INLINE
#pragma HLS ARRAY_PARTITION variable=tab complete dim=1
	if( sel == 0 ){
		return tab[0];}
	else{
		return tab[1];}
}

template <class T, int L, int Q, int LOG2L>
inline T LIST_MUX4 ( T tab[L], sc_uint<LOG2L> sel)
{
#pragma HLS INLINE
#pragma HLS ARRAY_PARTITION variable=tab complete dim=1
	if( sel == 0 ){
		return tab[0];}
	else if( sel == 1 ){
		return tab[1];}
	else if( sel == 2 ){
		return tab[2];}
	else{
		return tab[3];}
}

template <class T, int L, int Q, int LOG2L>
inline T LIST_MUX8 ( T tab[L], sc_uint<LOG2L> sel)
{
#pragma HLS INLINE
#pragma HLS ARRAY_PARTITION variable=tab complete dim=1
	if( sel == 0 ){
		return tab[0];}
	else if( sel == 1 ){
		return tab[1];}
	else if( sel == 2 ){
		return tab[2];}
	else if( sel == 3 ){
		return tab[3];}
	else if( sel == 4 ){
		return tab[4];}
	else if( sel == 5 ){
		return tab[5];}
	else if( sel == 6 ){
		return tab[6];}
	else{
		return tab[7];
	}
}

template <class T, int L, int Q, int LOG2L>
inline T LIST_MUX16 ( T tab[L], sc_uint<LOG2L> sel)
{
#pragma HLS INLINE
#pragma HLS ARRAY_PARTITION variable=tab complete dim=1
	if( sel == 0 ){
		return tab[0];}
	else if( sel == 1 ){
		return tab[1];}
	else if( sel == 2 ){
		return tab[2];}
	else if( sel == 3 ){
		return tab[3];}
	else if( sel == 4 ){
		return tab[4];}
	else if( sel == 5 ){
		return tab[5];}
	else if( sel == 6 ){
		return tab[6];}
	else if( sel == 7 ){
		return tab[7];}
	else if( sel == 8 ){
		return tab[8];}
	else if( sel == 9 ){
		return tab[9];}
	else if( sel == 10 ){
		return tab[10];}
	else if( sel == 11 ){
		return tab[11];}
	else if( sel == 12 ){
		return tab[12];}
	else if( sel == 13 ){
		return tab[13];}
	else if( sel == 14 ){
		return tab[14];}
	else{
		return tab[15];
	}
}

template <class T, int L, int Q, int LOG2L>
inline T LIST_MUX32 ( T tab[L], sc_uint<LOG2L> sel)
{
#pragma HLS INLINE
#pragma HLS ARRAY_PARTITION variable=tab complete dim=1
	if( sel == 0 ){
		return tab[0];}
	else if( sel == 1 ){
		return tab[1];}
	else if( sel == 2 ){
		return tab[2];}
	else if( sel == 3 ){
		return tab[3];}
	else if( sel == 4 ){
		return tab[4];}
	else if( sel == 5 ){
		return tab[5];}
	else if( sel == 6 ){
		return tab[6];}
	else if( sel == 7 ){
		return tab[7];}
	else if( sel == 8 ){
		return tab[8];}
	else if( sel == 9 ){
		return tab[9];}
	else if( sel == 10 ){
		return tab[10];}
	else if( sel == 11 ){
		return tab[11];}
	else if( sel == 12 ){
		return tab[12];}
	else if( sel == 13 ){
		return tab[13];}
	else if( sel == 14 ){
		return tab[14];}
	else if( sel == 15 ){
		return tab[15];}
	else if( sel == 16 ){
		return tab[16];}
	else if( sel == 17 ){
		return tab[17];}
	else if( sel == 18 ){
		return tab[18];}
	else if( sel == 19 ){
		return tab[19];}
	else if( sel == 20 ){
		return tab[20];}
	else if( sel == 21 ){
		return tab[21];}
	else if( sel == 22 ){
		return tab[22];}
	else if( sel == 23 ){
		return tab[23];}
	else if( sel == 24 ){
		return tab[24];}
	else if( sel == 25 ){
		return tab[25];}
	else if( sel == 26 ){
		return tab[26];}
	else if( sel == 27 ){
		return tab[27];}
	else if( sel == 28 ){
		return tab[28];}
	else if( sel == 29 ){
		return tab[29];}
	else if( sel == 30 ){
		return tab[30];}
	else{
		return tab[31];
	}
}

template <class T, int L, int Q, int LOG2L>
inline T LIST_MUX64 ( T tab[L], sc_uint<LOG2L> sel)
{
#pragma HLS INLINE
#pragma HLS ARRAY_PARTITION variable=tab complete dim=1
	if( sel == 0 ){
		return tab[0];}
	else if( sel == 1 ){
		return tab[1];}
	else if( sel == 2 ){
		return tab[2];}
	else if( sel == 3 ){
		return tab[3];}
	else if( sel == 4 ){
		return tab[4];}
	else if( sel == 5 ){
		return tab[5];}
	else if( sel == 6 ){
		return tab[6];}
	else if( sel == 7 ){
		return tab[7];}
	else if( sel == 8 ){
		return tab[8];}
	else if( sel == 9 ){
		return tab[9];}
	else if( sel == 10 ){
		return tab[10];}
	else if( sel == 11 ){
		return tab[11];}
	else if( sel == 12 ){
		return tab[12];}
	else if( sel == 13 ){
		return tab[13];}
	else if( sel == 14 ){
		return tab[14];}
	else if( sel == 15 ){
		return tab[15];}
	else if( sel == 16 ){
		return tab[16];}
	else if( sel == 17 ){
		return tab[17];}
	else if( sel == 18 ){
		return tab[18];}
	else if( sel == 19 ){
		return tab[19];}
	else if( sel == 20 ){
		return tab[20];}
	else if( sel == 21 ){
		return tab[21];}
	else if( sel == 22 ){
		return tab[22];}
	else if( sel == 23 ){
		return tab[23];}
	else if( sel == 24 ){
		return tab[24];}
	else if( sel == 25 ){
		return tab[25];}
	else if( sel == 26 ){
		return tab[26];}
	else if( sel == 27 ){
		return tab[27];}
	else if( sel == 28 ){
		return tab[28];}
	else if( sel == 29 ){
		return tab[29];}
	else if( sel == 30 ){
		return tab[30];}
	else if( sel == 31 ){
		return tab[31];}
	else if( sel == 32 ){
		return tab[32];}
	else if( sel == 33 ){
		return tab[33];}
	else if( sel == 34 ){
		return tab[34];}
	else if( sel == 35 ){
		return tab[35];}
	else if( sel == 36 ){
		return tab[36];}
	else if( sel == 37 ){
		return tab[37];}
	else if( sel == 38 ){
		return tab[38];}
	else if( sel == 39 ){
		return tab[39];}
	else if( sel == 40 ){
		return tab[40];}
	else if( sel == 41 ){
		return tab[41];}
	else if( sel == 42 ){
		return tab[42];}
	else if( sel == 43 ){
		return tab[43];}
	else if( sel == 44 ){
		return tab[44];}
	else if( sel == 45 ){
		return tab[45];}
	else if( sel == 46 ){
		return tab[46];}
	else if( sel == 47 ){
		return tab[47];}
	else if( sel == 48 ){
		return tab[48];}
	else if( sel == 49 ){
		return tab[49];}
	else if( sel == 50 ){
		return tab[50];}
	else if( sel == 51 ){
		return tab[51];}
	else if( sel == 52 ){
		return tab[52];}
	else if( sel == 53 ){
		return tab[53];}
	else if( sel == 54 ){
		return tab[54];}
	else if( sel == 55 ){
		return tab[55];}
	else if( sel == 56 ){
		return tab[56];}
	else if( sel == 57 ){
		return tab[57];}
	else if( sel == 58 ){
		return tab[58];}
	else if( sel == 59 ){
		return tab[59];}
	else if( sel == 60 ){
		return tab[60];}
	else if( sel == 61 ){
		return tab[61];}
	else if( sel == 62 ){
		return tab[62];}
	else{
		return tab[63];
	}
}

//*************************************************************************//
///////////////////////// RANK ORDER SORTER MUX /////////////////////////////
//*************************************************************************//

template <class T>
inline T RO_MUX2 ( T tab[2], sc_uint<1> position[2], sc_uint<1> indice)
{
#pragma HLS INLINE off
#pragma HLS ARRAY_PARTITION variable=tab complete dim=1
#pragma HLS ARRAY_PARTITION variable=position complete dim=1
	if(  position[0] == indice ){
		return tab[0];}
	else{
		return tab[1];}
}

template <class T>
inline T RO_MUX4 ( T tab[4], sc_uint<2> position[4], sc_uint<2> indice)
{
#pragma HLS INLINE off
#pragma HLS ARRAY_PARTITION variable=tab complete dim=1
#pragma HLS ARRAY_PARTITION variable=position complete dim=1
	if( position[0] == indice ){
		return tab[0];}
	else if( position[1] == indice ){
		return tab[1];}
	else if( position[2] == indice ){
		return tab[2];}
	else{
		return tab[3];}
}

template <class T>
inline T RO_MUX8 ( T tab[8], sc_uint<3> position[8], sc_uint<3> indice)
{
#pragma HLS INLINE off
#pragma HLS ARRAY_PARTITION variable=tab complete dim=1
#pragma HLS ARRAY_PARTITION variable=position complete dim=1
	if( position[0] == indice ){
		return tab[0];}
	else if( position[1] == indice ){
		return tab[1];}
	else if( position[2] == indice ){
		return tab[2];}
	else if( position[3] == indice ){
		return tab[3];}
	else if( position[4] == indice ){
		return tab[4];}
	else if( position[5] == indice ){
		return tab[5];}
	else if( position[6] == indice ){
		return tab[6];}
	else{
		return tab[7];
	}
}

template <class T>
inline T RO_MUX16 ( T tab[16], sc_uint<4> position[16], sc_uint<4> indice)
{
#pragma HLS INLINE off
#pragma HLS ARRAY_PARTITION variable=tab complete dim=1
#pragma HLS ARRAY_PARTITION variable=position complete dim=1
	if( position[0] == indice ){
		return tab[0];}
	else if( position[1] == indice ){
		return tab[1];}
	else if( position[2] == indice ){
		return tab[2];}
	else if( position[3] == indice ){
		return tab[3];}
	else if( position[4] == indice ){
		return tab[4];}
	else if( position[5] == indice ){
		return tab[5];}
	else if( position[6] == indice ){
		return tab[6];}
	else if( position[7] == indice ){
		return tab[7];}
	else if( position[8] == indice ){
		return tab[8];}
	else if( position[9] == indice ){
		return tab[9];}
	else if( position[10] == indice ){
		return tab[10];}
	else if( position[11] == indice ){
		return tab[11];}
	else if( position[12] == indice ){
		return tab[12];}
	else if( position[13] == indice ){
		return tab[13];}
	else if( position[14] == indice ){
		return tab[14];}
	else{
		return tab[15];
	}
}

template <class T>
inline T RO_MUX32 ( T tab[32], sc_uint<5> position[32], sc_uint<5> indice)
{
#pragma HLS INLINE off
#pragma HLS ARRAY_PARTITION variable=tab complete dim=1
#pragma HLS ARRAY_PARTITION variable=position complete dim=1
     if( position[0] == indice ){
         return tab[0];}
     if( position[1] == indice ){
         return tab[1];}
     if( position[2] == indice ){
         return tab[2];}
     if( position[3] == indice ){
         return tab[3];}
     if( position[4] == indice ){
         return tab[4];}
     if( position[5] == indice ){
         return tab[5];}
     if( position[6] == indice ){
         return tab[6];}
     if( position[7] == indice ){
         return tab[7];}
     if( position[8] == indice ){
         return tab[8];}
     if( position[9] == indice ){
         return tab[9];}
     if( position[10] == indice ){
         return tab[10];}
     if( position[11] == indice ){
         return tab[11];}
     if( position[12] == indice ){
         return tab[12];}
     if( position[13] == indice ){
         return tab[13];}
     if( position[14] == indice ){
         return tab[14];}
     if( position[15] == indice ){
         return tab[15];}
     if( position[16] == indice ){
         return tab[16];}
     if( position[17] == indice ){
         return tab[17];}
     if( position[18] == indice ){
         return tab[18];}
     if( position[19] == indice ){
         return tab[19];}
     if( position[20] == indice ){
         return tab[20];}
     if( position[21] == indice ){
         return tab[21];}
     if( position[22] == indice ){
         return tab[22];}
     if( position[23] == indice ){
         return tab[23];}
     if( position[24] == indice ){
         return tab[24];}
     if( position[25] == indice ){
         return tab[25];}
     if( position[26] == indice ){
         return tab[26];}
     if( position[27] == indice ){
         return tab[27];}
     if( position[28] == indice ){
         return tab[28];}
     if( position[29] == indice ){
         return tab[29];}
     if( position[30] == indice ){
         return tab[30];}
     else{
         return tab[31];}
}

template <class T>
inline T RO_MUX64 ( T tab[64], sc_uint<6> position[64], sc_uint<6> indice)
{
#pragma HLS INLINE off
#pragma HLS ARRAY_PARTITION variable=tab complete dim=1
#pragma HLS ARRAY_PARTITION variable=position complete dim=1
     if( position[0] == indice ){
         return tab[0];}
     if( position[1] == indice ){
         return tab[1];}
     if( position[2] == indice ){
         return tab[2];}
     if( position[3] == indice ){
         return tab[3];}
     if( position[4] == indice ){
         return tab[4];}
     if( position[5] == indice ){
         return tab[5];}
     if( position[6] == indice ){
         return tab[6];}
     if( position[7] == indice ){
         return tab[7];}
     if( position[8] == indice ){
         return tab[8];}
     if( position[9] == indice ){
         return tab[9];}
     if( position[10] == indice ){
         return tab[10];}
     if( position[11] == indice ){
         return tab[11];}
     if( position[12] == indice ){
         return tab[12];}
     if( position[13] == indice ){
         return tab[13];}
     if( position[14] == indice ){
         return tab[14];}
     if( position[15] == indice ){
         return tab[15];}
     if( position[16] == indice ){
         return tab[16];}
     if( position[17] == indice ){
         return tab[17];}
     if( position[18] == indice ){
         return tab[18];}
     if( position[19] == indice ){
         return tab[19];}
     if( position[20] == indice ){
         return tab[20];}
     if( position[21] == indice ){
         return tab[21];}
     if( position[22] == indice ){
         return tab[22];}
     if( position[23] == indice ){
         return tab[23];}
     if( position[24] == indice ){
         return tab[24];}
     if( position[25] == indice ){
         return tab[25];}
     if( position[26] == indice ){
         return tab[26];}
     if( position[27] == indice ){
         return tab[27];}
     if( position[28] == indice ){
         return tab[28];}
     if( position[29] == indice ){
         return tab[29];}
     if( position[30] == indice ){
         return tab[30];}
     if( position[31] == indice ){
         return tab[31];}
     if( position[32] == indice ){
         return tab[32];}
     if( position[33] == indice ){
         return tab[33];}
     if( position[34] == indice ){
         return tab[34];}
     if( position[35] == indice ){
         return tab[35];}
     if( position[36] == indice ){
         return tab[36];}
     if( position[37] == indice ){
         return tab[37];}
     if( position[38] == indice ){
         return tab[38];}
     if( position[39] == indice ){
         return tab[39];}
     if( position[40] == indice ){
         return tab[40];}
     if( position[41] == indice ){
         return tab[41];}
     if( position[42] == indice ){
         return tab[42];}
     if( position[43] == indice ){
         return tab[43];}
     if( position[44] == indice ){
         return tab[44];}
     if( position[45] == indice ){
         return tab[45];}
     if( position[46] == indice ){
         return tab[46];}
     if( position[47] == indice ){
         return tab[47];}
     if( position[48] == indice ){
         return tab[48];}
     if( position[49] == indice ){
         return tab[49];}
     if( position[50] == indice ){
         return tab[50];}
     if( position[51] == indice ){
         return tab[51];}
     if( position[52] == indice ){
         return tab[52];}
     if( position[53] == indice ){
         return tab[53];}
     if( position[54] == indice ){
         return tab[54];}
     if( position[55] == indice ){
         return tab[55];}
     if( position[56] == indice ){
         return tab[56];}
     if( position[57] == indice ){
         return tab[57];}
     if( position[58] == indice ){
         return tab[58];}
     if( position[59] == indice ){
         return tab[59];}
     if( position[60] == indice ){
         return tab[60];}
     if( position[61] == indice ){
         return tab[61];}
     if( position[62] == indice ){
         return tab[62];}
     else{
         return tab[63];}
}

template <class T>
inline T RO_MUX128 ( T tab[128], sc_uint<7> position[128], sc_uint<7> indice)
{
#pragma HLS INLINE off
#pragma HLS ARRAY_PARTITION variable=tab complete dim=1
#pragma HLS ARRAY_PARTITION variable=position complete dim=1
     if( position[0] == indice ){
         return tab[0];}
     if( position[1] == indice ){
         return tab[1];}
     if( position[2] == indice ){
         return tab[2];}
     if( position[3] == indice ){
         return tab[3];}
     if( position[4] == indice ){
         return tab[4];}
     if( position[5] == indice ){
         return tab[5];}
     if( position[6] == indice ){
         return tab[6];}
     if( position[7] == indice ){
         return tab[7];}
     if( position[8] == indice ){
         return tab[8];}
     if( position[9] == indice ){
         return tab[9];}
     if( position[10] == indice ){
         return tab[10];}
     if( position[11] == indice ){
         return tab[11];}
     if( position[12] == indice ){
         return tab[12];}
     if( position[13] == indice ){
         return tab[13];}
     if( position[14] == indice ){
         return tab[14];}
     if( position[15] == indice ){
         return tab[15];}
     if( position[16] == indice ){
         return tab[16];}
     if( position[17] == indice ){
         return tab[17];}
     if( position[18] == indice ){
         return tab[18];}
     if( position[19] == indice ){
         return tab[19];}
     if( position[20] == indice ){
         return tab[20];}
     if( position[21] == indice ){
         return tab[21];}
     if( position[22] == indice ){
         return tab[22];}
     if( position[23] == indice ){
         return tab[23];}
     if( position[24] == indice ){
         return tab[24];}
     if( position[25] == indice ){
         return tab[25];}
     if( position[26] == indice ){
         return tab[26];}
     if( position[27] == indice ){
         return tab[27];}
     if( position[28] == indice ){
         return tab[28];}
     if( position[29] == indice ){
         return tab[29];}
     if( position[30] == indice ){
         return tab[30];}
     if( position[31] == indice ){
         return tab[31];}
     if( position[32] == indice ){
         return tab[32];}
     if( position[33] == indice ){
         return tab[33];}
     if( position[34] == indice ){
         return tab[34];}
     if( position[35] == indice ){
         return tab[35];}
     if( position[36] == indice ){
         return tab[36];}
     if( position[37] == indice ){
         return tab[37];}
     if( position[38] == indice ){
         return tab[38];}
     if( position[39] == indice ){
         return tab[39];}
     if( position[40] == indice ){
         return tab[40];}
     if( position[41] == indice ){
         return tab[41];}
     if( position[42] == indice ){
         return tab[42];}
     if( position[43] == indice ){
         return tab[43];}
     if( position[44] == indice ){
         return tab[44];}
     if( position[45] == indice ){
         return tab[45];}
     if( position[46] == indice ){
         return tab[46];}
     if( position[47] == indice ){
         return tab[47];}
     if( position[48] == indice ){
         return tab[48];}
     if( position[49] == indice ){
         return tab[49];}
     if( position[50] == indice ){
         return tab[50];}
     if( position[51] == indice ){
         return tab[51];}
     if( position[52] == indice ){
         return tab[52];}
     if( position[53] == indice ){
         return tab[53];}
     if( position[54] == indice ){
         return tab[54];}
     if( position[55] == indice ){
         return tab[55];}
     if( position[56] == indice ){
         return tab[56];}
     if( position[57] == indice ){
         return tab[57];}
     if( position[58] == indice ){
         return tab[58];}
     if( position[59] == indice ){
         return tab[59];}
     if( position[60] == indice ){
         return tab[60];}
     if( position[61] == indice ){
         return tab[61];}
     if( position[62] == indice ){
         return tab[62];}
     if( position[63] == indice ){
         return tab[63];}
     if( position[64] == indice ){
         return tab[64];}
     if( position[65] == indice ){
         return tab[65];}
     if( position[66] == indice ){
         return tab[66];}
     if( position[67] == indice ){
         return tab[67];}
     if( position[68] == indice ){
         return tab[68];}
     if( position[69] == indice ){
         return tab[69];}
     if( position[70] == indice ){
         return tab[70];}
     if( position[71] == indice ){
         return tab[71];}
     if( position[72] == indice ){
         return tab[72];}
     if( position[73] == indice ){
         return tab[73];}
     if( position[74] == indice ){
         return tab[74];}
     if( position[75] == indice ){
         return tab[75];}
     if( position[76] == indice ){
         return tab[76];}
     if( position[77] == indice ){
         return tab[77];}
     if( position[78] == indice ){
         return tab[78];}
     if( position[79] == indice ){
         return tab[79];}
     if( position[80] == indice ){
         return tab[80];}
     if( position[81] == indice ){
         return tab[81];}
     if( position[82] == indice ){
         return tab[82];}
     if( position[83] == indice ){
         return tab[83];}
     if( position[84] == indice ){
         return tab[84];}
     if( position[85] == indice ){
         return tab[85];}
     if( position[86] == indice ){
         return tab[86];}
     if( position[87] == indice ){
         return tab[87];}
     if( position[88] == indice ){
         return tab[88];}
     if( position[89] == indice ){
         return tab[89];}
     if( position[90] == indice ){
         return tab[90];}
     if( position[91] == indice ){
         return tab[91];}
     if( position[92] == indice ){
         return tab[92];}
     if( position[93] == indice ){
         return tab[93];}
     if( position[94] == indice ){
         return tab[94];}
     if( position[95] == indice ){
         return tab[95];}
     if( position[96] == indice ){
         return tab[96];}
     if( position[97] == indice ){
         return tab[97];}
     if( position[98] == indice ){
         return tab[98];}
     if( position[99] == indice ){
         return tab[99];}
     if( position[100] == indice ){
         return tab[100];}
     if( position[101] == indice ){
         return tab[101];}
     if( position[102] == indice ){
         return tab[102];}
     if( position[103] == indice ){
         return tab[103];}
     if( position[104] == indice ){
         return tab[104];}
     if( position[105] == indice ){
         return tab[105];}
     if( position[106] == indice ){
         return tab[106];}
     if( position[107] == indice ){
         return tab[107];}
     if( position[108] == indice ){
         return tab[108];}
     if( position[109] == indice ){
         return tab[109];}
     if( position[110] == indice ){
         return tab[110];}
     if( position[111] == indice ){
         return tab[111];}
     if( position[112] == indice ){
         return tab[112];}
     if( position[113] == indice ){
         return tab[113];}
     if( position[114] == indice ){
         return tab[114];}
     if( position[115] == indice ){
         return tab[115];}
     if( position[116] == indice ){
         return tab[116];}
     if( position[117] == indice ){
         return tab[117];}
     if( position[118] == indice ){
         return tab[118];}
     if( position[119] == indice ){
         return tab[119];}
     if( position[120] == indice ){
         return tab[120];}
     if( position[121] == indice ){
         return tab[121];}
     if( position[122] == indice ){
         return tab[122];}
     if( position[123] == indice ){
         return tab[123];}
     if( position[124] == indice ){
         return tab[124];}
     if( position[125] == indice ){
         return tab[125];}
     if( position[126] == indice ){
         return tab[126];}
     else{
         return tab[127];}
}


//*************************************************************************//
///////////////////////// RANK ORDER SORTER MUX 2 ///////////////////////////
//*************************************************************************//

/*
template <class T>
inline T RO_MUX2 ( T tab[2], sc_uint<1> position_0, sc_uint<1> position_1, sc_uint<1> indice)
{
#pragma HLS INLINE
	if(  position_0 == indice ){
		return tab[0];}
	else {
		return tab[1];}
}

template <class T>
inline T RO_MUX4 ( T tab[4], sc_uint<2> position_0, sc_uint<2> position_1, sc_uint<2> position_2, sc_uint<2> position_3, sc_uint<2> indice)
{
#pragma HLS INLINE
	if( position_0 == indice ){
		return tab[0];}
	else if( position_1 == indice ){
		return tab[1];}
	else if( position_2 == indice ){
		return tab[2];}
	else {
		return tab[3];}
}

template <class T>
inline T RO_MUX8 ( T tab[8], sc_uint<3> position_0, sc_uint<3> position_1, sc_uint<3> position_2, sc_uint<3> position_3, sc_uint<3> position_4, sc_uint<3> position_5, sc_uint<3> position_6, sc_uint<3> position_7, sc_uint<3> indice)
{
#pragma HLS INLINE
	if( position_0 == indice ){
		return tab[0];}
	else if( position_1 == indice ){
		return tab[1];}
	else if( position_2 == indice ){
		return tab[2];}
	else if( position_3 == indice ){
		return tab[3];}
	else if( position_4 == indice ){
		return tab[4];}
	else if( position_5 == indice ){
		return tab[5];}
	else if( position_6 == indice ){
		return tab[6];}
	else {
		return tab[7];
	}
}

template <class T>
inline T RO_MUX16 ( T tab[16], sc_uint<4> position_0, sc_uint<4> position_1, sc_uint<4> position_2, sc_uint<4> position_3, sc_uint<4> position_4, sc_uint<4> position_5, sc_uint<4> position_6, sc_uint<4> position_7, sc_uint<4> position_8, sc_uint<4> position_9, sc_uint<4> position_10, sc_uint<4> position_11, sc_uint<4> position_12, sc_uint<4> position_13, sc_uint<4> position_14, sc_uint<4> position_15, sc_uint<4> indice)
{
#pragma HLS INLINE off
	if( position_0 == indice ){
		return tab[0];}
	else if( position_1 == indice ){
		return tab[1];}
	else if( position_2 == indice ){
		return tab[2];}
	else if( position_3 == indice ){
		return tab[3];}
	else if( position_4 == indice ){
		return tab[4];}
	else if( position_5 == indice ){
		return tab[5];}
	else if( position_6 == indice ){
		return tab[6];}
	else if( position_7 == indice ){
		return tab[7];}
	else if( position_8 == indice ){
		return tab[8];}
	else if( position_9 == indice ){
		return tab[9];}
	else if( position_10 == indice ){
		return tab[10];}
	else if( position_11 == indice ){
		return tab[11];}
	else if( position_12 == indice ){
		return tab[12];}
	else if( position_13 == indice ){
		return tab[13];}
	else if( position_14 == indice ){
		return tab[14];}
	else {
		return tab[15];
	}
}

template <class T>
inline T RO_MUX32 ( T tab[32], sc_uint<5> position_0, sc_uint<5> position_1, sc_uint<5> position_2, sc_uint<5> position_3, sc_uint<5> position_4, sc_uint<5> position_5, sc_uint<5> position_6, sc_uint<5> position_7, sc_uint<5> position_8, sc_uint<5> position_9, sc_uint<5> position_10, sc_uint<5> position_11, sc_uint<5> position_12, sc_uint<5> position_13, sc_uint<5> position_14, sc_uint<5> position_15, sc_uint<5> position_16, sc_uint<5> position_17, sc_uint<5> position_18, sc_uint<5> position_19, sc_uint<5> position_20, sc_uint<5> position_21, sc_uint<5> position_22, sc_uint<5> position_23, sc_uint<5> position_24, sc_uint<5> position_25, sc_uint<5> position_26, sc_uint<5> position_27, sc_uint<5> position_28, sc_uint<5> position_29, sc_uint<5> position_30, sc_uint<5> position_31, sc_uint<5> indice)
{
#pragma HLS INLINE
	if( position_0 == indice ){
		return tab[0];}
	else if( position_1 == indice ){
		return tab[1];}
	else if( position_2 == indice ){
		return tab[2];}
	else if( position_3 == indice ){
		return tab[3];}
	else if( position_4 == indice ){
		return tab[4];}
	else if( position_5 == indice ){
		return tab[5];}
	else if( position_6 == indice ){
		return tab[6];}
	else if( position_7 == indice ){
		return tab[7];}
	else if( position_8 == indice ){
		return tab[8];}
	else if( position_9 == indice ){
		return tab[9];}
	else if( position_10 == indice ){
		return tab[10];}
	else if( position_11 == indice ){
		return tab[11];}
	else if( position_12 == indice ){
		return tab[12];}
	else if( position_13 == indice ){
		return tab[13];}
	else if( position_14 == indice ){
		return tab[14];}
	else if( position_15 == indice ){
		return tab[15];}
	else if( position_16 == indice ){
		return tab[16];}
	else if( position_17 == indice ){
		return tab[17];}
	else if( position_18 == indice ){
		return tab[18];}
	else if( position_19 == indice ){
		return tab[19];}
	else if( position_20 == indice ){
		return tab[20];}
	else if( position_21 == indice ){
		return tab[21];}
	else if( position_22 == indice ){
		return tab[22];}
	else if( position_23 == indice ){
		return tab[23];}
	else if( position_24 == indice ){
		return tab[24];}
	else if( position_25 == indice ){
		return tab[25];}
	else if( position_26 == indice ){
		return tab[26];}
	else if( position_27 == indice ){
		return tab[27];}
	else if( position_28 == indice ){
		return tab[28];}
	else if( position_29 == indice ){
		return tab[29];}
	else if( position_30 == indice ){
		return tab[30];}
	else{
		return tab[31];
	}
}
*/
#endif
