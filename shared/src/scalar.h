#ifndef SCALAR
#define SCALAR

template <int Q>
inline void SHOW(sc_bigint<Q> a){
	printf("%3d ", (int)a.to_int());
}

template <int Q>
inline sc_bigint<Q> qmin(sc_bigint<Q> a, sc_bigint<Q> b){
#pragma HLS INLINE
	return (a < b) ? a : b;
}

template <int Q>
inline sc_bigint<Q> qsat(sc_bigint<Q+1> a){
#pragma HLS INLINE
	if     ( a > (sc_bigint<Q>)((sc_uint<1>)0, (sc_uint<Q-1>)0xFFFFFFFF) ) return (sc_bigint<Q>)((sc_uint<1>)0, (sc_uint<Q-1>)0xFFFFFFFF);  // LLR_MAXV
	else if( a < (sc_bigint<Q>)((sc_uint<1>)1, (sc_uint<Q-1>)0x01) ) return (sc_bigint<Q>)((sc_uint<1>)1, (sc_uint<Q-1>)0x01);  // LLR_MINV
	else return (sc_bigint<Q>) a;
}

template <int Q>
inline sc_biguint<1> qsign(sc_bigint<Q> value){  // positif : return 0 ; negatif : 1
#pragma HLS INLINE
	sc_biguint<1> sign = (sc_biguint<1>) value[Q - 1];
	return sign;
}

template <int Q>
inline sc_bigint<Q> qsign(sc_bigint<Q> a, sc_biguint<1> b){
#pragma HLS INLINE
	sc_bigint<Q> positif = +a;
	sc_bigint<Q> negatif = -a;
    return ((sc_uint<1>)b)? negatif : positif ;
}


template <int Q>
inline sc_bigint<Q> qabs(sc_bigint<Q> value){
#pragma HLS INLINE
	sc_bigint<Q> positif = +value;
	sc_bigint<Q> negatif = -value;
	sc_bigint<Q> result  = (value [Q-1]  == 1) ? negatif : positif;
	return result;
}

template <int Q>
inline sc_bigint<Q> qadd(sc_bigint<Q> a, sc_bigint<Q> b){
#pragma HLS INLINE
	sc_bigint<Q+1> r = ((sc_bigint<Q+1>)a) + ((sc_bigint<Q+1>)b);
	sc_bigint<Q> sat = qsat<Q>(r); //if( r > LLR_MAXV ) return LLR_MAXV; else if( r < LLR_MINV ) return LLR_MINV; else return r;
	return sat;
}

template <int Q>
inline sc_bigint<Q+1> qadd_nosat(sc_bigint<Q> a, sc_bigint<Q> b){
#pragma HLS INLINE
	sc_bigint<Q+1> r = ((sc_bigint<Q+1>)a) + ((sc_bigint<Q+1>)b);
	return r;
}

template <int Q>
inline sc_bigint<Q> qsub(sc_bigint<Q> a, sc_bigint<Q> b){
#pragma HLS INLINE
	sc_bigint<Q+1> r = ((sc_bigint<Q+1>)a) - ((sc_bigint<Q+1>)b);
	sc_bigint<Q> sat = qsat<Q>(r); // if ( r > LLR_MAXV ) return LLR_MAXV;else if( r < LLR_MINV ) return LLR_MINV; else return r;
	return r;
}

template <int Q>
inline sc_bigint<Q+1> qsub_nosat(sc_bigint<Q> a, sc_bigint<Q> b){
#pragma HLS INLINE
	sc_bigint<Q+1> r = ((sc_bigint<Q+1>)a) - ((sc_bigint<Q+1>)b);
	return r;
}


//////////////////////////Signe et Magnetude/////////////////////////////////

template <int Q>
inline void SHOW_sm(sc_biguint<Q> a){
	sc_bigint<Q> b = (a[Q-1], (~ (sc_biguint<Q-1>)a.range(Q-2 , 0)) + 1);
	sc_bigint<Q> c = ((sc_uint<1>)a[Q-1]) ? b : a;
	printf("%3d ", (int)c.to_int());
}

template <int Q>
inline sc_biguint<Q> qmin_sm(sc_biguint<Q> a, sc_biguint<Q> b){
#pragma HLS INLINE
	return (a < b) ? a : b;
}

template <int Q>
inline sc_biguint<Q> qsat_sm(sc_biguint<Q+1> a){
#pragma HLS INLINE
	if     ( a > (sc_bigint<Q>)((sc_uint<1>)0, (sc_uint<Q-1>)0xFFFFFF) ) return (sc_bigint<Q>)((sc_uint<1>)0, (sc_uint<Q-1>)0xFFFFFF); // LLR_MAXV
	else return (sc_biguint<Q>) a;
}

template <int Q>
inline sc_biguint<1> qsign_sm(sc_biguint<Q> value){  // positif : return 0 ; negatif : 1
#pragma HLS INLINE
	sc_biguint<1> sign = (sc_biguint<1>) value[Q - 1];
	return sign;
}

template <int Q>
inline sc_biguint<Q-1> qabs_sm(sc_biguint<Q> value){
#pragma HLS INLINE
	sc_biguint<Q-1> result  = (sc_biguint<Q-1>) value.range(Q-2 , 0);
	return result;
}

template <int Q>
inline sc_biguint<Q+1> qadd_nosat_sm(sc_biguint<Q> a, sc_biguint<Q> b){
#pragma HLS INLINE
	sc_biguint<Q+1> r = ((sc_biguint<Q+1>)a) + ((sc_biguint<Q+1>)b);
	return r;
}

template <int Q>
inline sc_biguint<Q> qadd_carry_sm(sc_biguint<Q> a, sc_biguint<Q> b, sc_biguint<1> carry){
#pragma HLS INLINE
	sc_biguint<Q> r = a + b + carry ;
	return r;
}

template <int Q>
inline sc_biguint<Q> qnot_sm(sc_biguint<Q> a){
#pragma HLS INLINE
	return (~a);
}

template <int Q>
inline sc_biguint<Q+1> qfull_adder_sm(sc_biguint<Q> a, sc_biguint<Q> b){
#pragma HLS INLINE

	sc_biguint<1> siga = (sc_biguint<1>) a[Q - 1];
	sc_biguint<1> sigb = (sc_biguint<1>) b[Q - 1];
	sc_biguint<1> Xsig = siga ^ sigb;

	sc_biguint<Q> absla = (sc_biguint<Q>) a.range(Q-2 , 0);
	sc_biguint<Q> invla = ~absla;

	sc_biguint<Q> abslb = (sc_biguint<Q>) b.range(Q-2 , 0);
	sc_biguint<Q> invlb = ~abslb;

	sc_biguint<1> is_min = ( absla < abslb ) ? 1 : 0; // if absla < abslb return 1

	sc_biguint<1> sel_a = Xsig & is_min;
	sc_biguint<1> sel_b = Xsig & ( ~is_min );

	sc_biguint<Q> absA = (sc_uint<1>)sel_a ? invla : absla ;
	sc_biguint<Q> absB = (sc_uint<1>)sel_b ? invlb : abslb ;

	sc_biguint<Q> Somme = absA + absB + Xsig;
	sc_biguint<1> SigSomme = (sc_uint<1>) is_min ? sigb : siga;
	sc_biguint<Q+1> result = (SigSomme , Somme);

	return result;
}

template <int Q>
inline sc_biguint<Q> qfull_adder_sat_sm(sc_biguint<Q> a, sc_biguint<Q> b){
#pragma HLS INLINE

	sc_biguint<1> siga = (sc_biguint<1>) a[Q - 1];
	sc_biguint<1> sigb = (sc_biguint<1>) b[Q - 1];
	sc_biguint<1> Xsig = siga ^ sigb;

	sc_biguint<Q> absla = (sc_biguint<Q>) a.range(Q-2 , 0);
	sc_biguint<Q> invla = ~absla;

	sc_biguint<Q> abslb = (sc_biguint<Q>) b.range(Q-2 , 0);
	sc_biguint<Q> invlb = ~abslb;

	sc_biguint<1> is_min = ( absla < abslb ) ? 1 : 0; // if absla < abslb return 1

	sc_biguint<1> sel_a = Xsig & is_min;
	sc_biguint<1> sel_b = Xsig & ( ~is_min );

	sc_biguint<Q> absA = (sc_uint<1>)sel_a ? invla : absla ;
	sc_biguint<Q> absB = (sc_uint<1>)sel_b ? invlb : abslb ;

	sc_biguint<Q> Somme = absA + absB + Xsig;

	sc_biguint<Q-1> Sat = qsat_sm<Q-1>(Somme);

	sc_biguint<1> SigSomme = (sc_uint<1>) is_min ? sigb : siga;
	sc_biguint<Q> result = (SigSomme , Sat);

	return result;
}

template <int Q>
inline sc_biguint<Q+1> qfull_add_sub_sm(sc_biguint<Q> a, sc_biguint<Q> b, sc_biguint<1> s){  // s = 1 : b-a, s = 0 : b+a
#pragma HLS INLINE

	sc_biguint<1> sla = (sc_biguint<1>) a[Q - 1];
	sc_biguint<1> siga = sla ^ s;

	sc_biguint<1> sigb = (sc_biguint<1>) b[Q - 1];
	sc_biguint<1> Xsig = siga ^ sigb;

	sc_biguint<Q> absla = (sc_biguint<Q>) a.range(Q-2 , 0);
	sc_biguint<Q> invla = ~absla;

	sc_biguint<Q> abslb = (sc_biguint<Q>) b.range(Q-2 , 0);
	sc_biguint<Q> invlb = ~abslb;

	sc_biguint<1> is_min = ( absla < abslb ) ? 1 : 0; // if absla < abslb return 1

	sc_biguint<1> sel_a = Xsig & is_min;
	sc_biguint<1> sel_b = Xsig & ( ~is_min );

	sc_biguint<Q> absA = (sc_uint<1>)sel_a ? invla : absla ;
	sc_biguint<Q> absB = (sc_uint<1>)sel_b ? invlb : abslb ;

	sc_biguint<Q> Somme = absA + absB + Xsig;
	sc_biguint<1> SigSomme = (sc_uint<1>) is_min ? sigb : siga;
	sc_biguint<Q+1> result = (SigSomme , Somme);

	return result;
}

/////////////////////////////////////////////////////////////////////////////

template <int Q>
inline sc_biguint<Q> qconv_format(sc_biguint<Q> a){  // CA2 to SIGMAG / SIGMAG to CA2
#pragma HLS INLINE
	sc_biguint<Q-1> abs = a.range(Q-2, 0);
	sc_biguint< Q > ext = (sc_biguint< Q >) abs;
	sc_biguint< Q > inv = (~ext);
	sc_biguint< Q > add = inv + (sc_biguint<1>) 1;
	sc_uint< 1 > sig = (sc_uint< 1 >) a.range(Q-1, Q-1);
	sc_biguint< Q > res = sig ? add : a;
	return res;
}

/////////////////////////////////////////////////////////////////////////////

inline sc_biguint<1> qxor(sc_biguint<1> a, sc_biguint<1> b){
#pragma HLS INLINE
	return a ^ b;
}

inline sc_biguint<1> qand(sc_biguint<1> a, sc_biguint<1> b){
#pragma HLS INLINE
	return a & b;
}


#endif
