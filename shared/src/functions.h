#ifndef FUNCTIONS
#define FUNCTIONS

#include "systemc.h"
#include "vector.h"

//*************************************************************************//
//** 					STACK MANAGEMENT FUNCTIONS						 **//
//*************************************************************************//

template <int Q, int D>
sc_biguint< Q* D > push_stack ( sc_biguint< Q* D > stack, sc_biguint< Q > value ) { // push a value
#pragma HLS INLINE

	return ( (sc_biguint< Q* (D - 1) >)stack.range( Q* (D - 1) - 1, 0), value);
}

template <int Q, int D>
sc_biguint< Q* D > pop_stack ( sc_biguint< Q* D > stack, sc_biguint< Q > value ) {  // pop with value
#pragma HLS INLINE

	return ( value, (sc_biguint< Q* (D - 1) >)stack.range( (Q* D) - 1, Q));
}

template <int Q, int D>
sc_biguint< Q* D > write_stack ( sc_biguint< Q* D > stack, sc_biguint< Q > value ) {  // write first element
#pragma HLS INLINE

	return ( (sc_biguint< Q* (D - 1) >)stack.range( (Q* D) - 1, Q ), value );
}

template <int Q, int D>
sc_biguint< Q > read_stack ( sc_biguint< Q* D > stack, sc_uint< 8 > adr ) { // read adr th element
#pragma HLS INLINE

	return stack.range( (Q * adr) - 1, Q * (adr -1) );
}


//*************************************************************************//
//** 						POLAR CODE FUNCTIONS						 **//
//*************************************************************************//

//*************************************************************************//
//////////////////////////Complement a 2 ////////////////////////////////////
//*************************************************************************//

template <int P, int Q>  // P : Parallelism , Q : Quantification
sc_bigint<P*Q> F_function_C2(sc_bigint<P*Q> la, sc_bigint<P*Q> lb)
{
#pragma HLS INLINE
	sc_bigint<P*Q> absa = VECTOR_ABS<P,Q> (la);
	sc_bigint<P*Q> absb = VECTOR_ABS<P,Q> (lb);
	sc_bigint<P*Q> min = VECTOR_MIN<P,Q> (absa,absb);
	sc_biguint<P> signa = VECTOR_SIGN<P,Q> (la);
	sc_biguint<P> signb = VECTOR_SIGN<P,Q> (lb);
	sc_biguint<P> sig = VECTOR_XOR<P>(signb,signa);
	sc_bigint<P*Q> result = VECTOR_SIGN<P,Q> (min,sig);

	return result;
}

template <int P, int Q>
sc_bigint<P*Q> G_function_C2(sc_bigint<P*Q> la, sc_bigint<P*Q> lb, sc_biguint<P> sa)
{
#pragma HLS INLINE
	sc_bigint<P*(Q+1)> add =  VECTOR_ADD_NOSAT<P,Q> (lb,la);
	sc_bigint<P*(Q+1)> sub =  VECTOR_SUB_NOSAT<P,Q> (lb,la);

	//sc_bigint<P*(Q+1)> g = VECTOR_MUX<P,(Q+1)> (add, sub, sa); // return add when sa=1, sub otherwise
	sc_bigint<P*(Q+1)> g = VECTOR_MUX<P,(Q+1)> (sub, add, sa);
	sc_bigint<P*Q> result = VECTOR_SAT<P,Q> (g);

	return result;
}

template <int P, int Q>
sc_bigint<P*(Q+1)> G_extended_C2(sc_bigint<P*Q> la, sc_bigint<P*Q> lb, sc_biguint<P> sa) // No saturation
{
#pragma HLS INLINE
	sc_bigint<P*(Q+1)> add =  VECTOR_ADD_NOSAT<P,Q> (lb,la);
	sc_bigint<P*(Q+1)> sub =  VECTOR_SUB_NOSAT<P,Q> (lb,la);

	//sc_bigint<P*(Q+1)> result = VECTOR_MUX<P,(Q+1)> (add, sub, sa); // return add when sa=1, sub otherwise
	sc_bigint<P*(Q+1)> result = VECTOR_MUX<P,(Q+1)> (sub, add, sa);

	return result;
}

template <int P, int Q>
sc_biguint<P> F_simplified_C2(sc_bigint<P*Q> la, sc_bigint<P*Q> lb, sc_biguint<P> fb)
{
#pragma HLS INLINE
	sc_biguint<P> signa = VECTOR_SIGN<P,Q> (la);
	sc_biguint<P> signb = VECTOR_SIGN<P,Q> (lb);
	sc_biguint<P> sig = VECTOR_XOR<P> (signb, signa);

	sc_biguint<P> result = VECTOR_AND<P> (sig, fb); // Partial Sum computation

	return result;
}

template <int P, int Q>
sc_biguint<P> G_simplified_C2(sc_bigint<P*Q> la, sc_bigint<P*Q> lb, sc_biguint<P> sa, sc_biguint<P> fb)
{
#pragma HLS INLINE
	sc_bigint<P*(Q+1)> add =  VECTOR_ADD_NOSAT<P,Q> (lb,la);
	sc_bigint<P*(Q+1)> sub =  VECTOR_SUB_NOSAT<P,Q> (lb,la);

	//sc_bigint<P*(Q+1)> g =  VECTOR_MUX<P,(Q+1)> (add, sub, sa);
	sc_bigint<P*(Q+1)> g =  VECTOR_MUX<P,(Q+1)> (sub, add, sa);

	sc_biguint<P> sig = VECTOR_SIGN<P,(Q+1)> (g);

	sc_biguint<P> result = VECTOR_AND<P> (sig, fb); // Partial Sum computation

	return result;
}

//*************************************************************************//
//////////////////////////Signe et Magnetude/////////////////////////////////
//*************************************************************************//

template <int P, int Q>
sc_bigint<P*Q> F_function_SM(sc_bigint<P*Q> la1, sc_bigint<P*Q> lb1)
{
#pragma HLS INLINE
	sc_biguint<P*Q> la = (sc_biguint<P*Q>) la1; // cast
	sc_biguint<P*Q> lb = (sc_biguint<P*Q>) lb1;

	sc_biguint<P*(Q-1)> absa = VECTOR_ABS_SM<P,Q> (la);
	sc_biguint<P*(Q-1)> absb = VECTOR_ABS_SM<P,Q> (lb);

	sc_biguint<P*(Q-1)> min = VECTOR_MIN_SM<P,Q-1> (absa, absb);

	sc_biguint<P> signa = VECTOR_SIGN_SM<P,Q> (la);
	sc_biguint<P> signb = VECTOR_SIGN_SM<P,Q> (lb);
	sc_biguint<P> sig = VECTOR_XOR<P> (signb,signa);

	sc_biguint<P*Q> result = VECTOR_CONCAT_SM<P,Q> (sig, min);

	sc_bigint<P*Q> result1 = (sc_bigint<P*Q>) result; //cast

	return result1;
}

template <int P, int Q>
sc_bigint<P*Q> G_function_SM(sc_bigint<P*Q> la1, sc_bigint<P*Q> lb1, sc_biguint<P> sa)
{
#pragma HLS INLINE
	sc_biguint<P*Q> la = (sc_biguint<P*Q>) la1; // cast
	sc_biguint<P*Q> lb = (sc_biguint<P*Q>) lb1;

	/*
	sc_biguint<P> sla = VECTOR_SIGN_SM<P,Q> ( la );
	//sc_biguint<P> n_sa = (sc_biguint<P>) VECTOR_NOT_SM<P,1> ( (sc_biguint<P>) sa );
	sc_biguint<P> n_sa =  sa;
	sc_biguint<P> sigla = VECTOR_XOR<P>( sla, n_sa );
	sc_biguint<P> siglb = VECTOR_SIGN_SM<P,Q> (lb);
	sc_biguint<P> Xsig = VECTOR_XOR<P>(sigla, siglb);

	sc_biguint<P*(Q-1)> absla = VECTOR_ABS_SM<P,Q> (la);
	sc_biguint<P*(Q)> absla_e = VECTOR_EXT_SM<P,Q-1> (absla);
	sc_biguint<P*(Q)> invla = VECTOR_NOT_SM<P,Q> (absla_e); // (not absla)

	sc_biguint<P*(Q-1)> abslb = VECTOR_ABS_SM<P,Q> (lb);
	sc_biguint<P*(Q)> abslb_e = VECTOR_EXT_SM<P,Q-1> (abslb);
	sc_biguint<P*(Q)> invlb = VECTOR_NOT_SM<P,Q> (abslb_e); // (not abslb)

	sc_biguint<P> is_min = VECTOR_IS_MIN_SM<P,Q-1> ( absla, abslb ); // if absla < abslb return 1

	sc_biguint<P> sel_a = VECTOR_AND<P> (Xsig, is_min);
	sc_biguint<P> sel_b = VECTOR_AND<P> (Xsig, (sc_biguint<P>) VECTOR_NOT_SM<P,1> ( (sc_biguint<P>) is_min ));

	sc_biguint<P*(Q)> absA = VECTOR_MUX_SM<P,Q> (invla, absla_e, sel_a) ;
	sc_biguint<P*(Q)> absB = VECTOR_MUX_SM<P,Q> (invlb, abslb_e, sel_b) ;

	sc_biguint<P*Q> Somme = VECTOR_ADD_CARRY_SM<P,Q> ( absA, absB, Xsig );

	sc_biguint<P*(Q-1)> Sat = VECTOR_SAT_SM<P,Q-1> (Somme);

	sc_biguint<P> SigSomme = (sc_biguint<P>) VECTOR_MUX_SM<P,1> ((sc_biguint<P>) siglb, (sc_biguint<P>) sigla, is_min);

	sc_biguint<P*Q> result = VECTOR_CONCAT_SM<P,Q> (SigSomme , Sat);
	*/
	sc_biguint<P*(Q+1)> somme = VECTOR_FULL_ADD_SUB_SM< P, Q >(la, lb, sa);
	sc_biguint<P*Q> abs = VECTOR_ABS_SM< P, Q+1 >(somme);
	sc_biguint<P> sign = VECTOR_SIGN_SM< P, Q+1 >(somme);
	sc_biguint<P*(Q-1)> sat = VECTOR_SAT_SM<P,Q-1> (abs);

	sc_biguint<P*Q> result = VECTOR_CONCAT_SM<P,Q> (sign , sat);
	sc_bigint<P*Q> result1 = (sc_bigint<P*Q>) result; //cast

	return result1;
}

template <int P, int Q>
sc_bigint<P*(Q+1)> G_extended_SM(sc_bigint<P*Q> la1, sc_bigint<P*Q> lb1, sc_biguint<P> sa)
{
#pragma HLS INLINE
	sc_biguint<P*Q> la = (sc_biguint<P*Q>) la1; // cast
	sc_biguint<P*Q> lb = (sc_biguint<P*Q>) lb1;

	/*
	sc_biguint<P> sla = VECTOR_SIGN_SM<P,Q> ( la );
	//sc_biguint<P> n_sa = (sc_biguint<P>) VECTOR_NOT_SM<P,1> ( (sc_biguint<P>) sa );
	sc_biguint<P> n_sa =  sa;
	sc_biguint<P> sigla = VECTOR_XOR<P>( sla, n_sa );
	sc_biguint<P> siglb = VECTOR_SIGN_SM<P,Q> (lb);
	sc_biguint<P> Xsig = VECTOR_XOR<P>(sigla, siglb);

	sc_biguint<P*(Q-1)> absla = VECTOR_ABS_SM<P,Q> (la);
	sc_biguint<P*(Q)> absla_e = VECTOR_EXT_SM<P,Q-1> (absla);
	sc_biguint<P*(Q)> invla = VECTOR_NOT_SM<P,Q> (absla_e); // (not absla)

	sc_biguint<P*(Q-1)> abslb = VECTOR_ABS_SM<P,Q> (lb);
	sc_biguint<P*(Q)> abslb_e = VECTOR_EXT_SM<P,Q-1> (abslb);
	sc_biguint<P*(Q)> invlb = VECTOR_NOT_SM<P,Q> (abslb_e); // (not abslb)

	sc_biguint<P> is_min = VECTOR_IS_MIN_SM<P,Q-1> ( absla, abslb ); // if absla < abslb return 1

	sc_biguint<P> sel_a = VECTOR_AND<P> (Xsig, is_min);
	sc_biguint<P> sel_b = VECTOR_AND<P> (Xsig, (sc_biguint<P>) VECTOR_NOT_SM<P,1> ( (sc_biguint<P>) is_min ));

	sc_biguint<P*(Q)> absA = VECTOR_MUX_SM<P,Q> (invla, absla_e, sel_a) ;
	sc_biguint<P*(Q)> absB = VECTOR_MUX_SM<P,Q> (invlb, abslb_e, sel_b) ;

	sc_biguint<P*Q> Somme = VECTOR_ADD_CARRY_SM<P,Q> ( absA, absB, Xsig );

	sc_biguint<P> SigSomme = (sc_biguint<P>) VECTOR_MUX_SM<P,1> ((sc_biguint<P>) siglb, (sc_biguint<P>) sigla, is_min);

	sc_biguint<P*(Q+1)> result = VECTOR_CONCAT_SM<P,(Q+1)> (SigSomme , Somme);
	*/
	sc_biguint<P*(Q+1)> result = VECTOR_FULL_ADD_SUB_SM< P, Q >(la, lb, sa);

	sc_bigint<P*(Q+1)> result1 = (sc_bigint<P*(Q+1)>) result; //cast

	return result1;
}

template <int P, int Q>
sc_biguint<P> F_simplified_SM(sc_bigint<P*Q> la1, sc_bigint<P*Q> lb1, sc_biguint<P> fb)
{
#pragma HLS INLINE
	sc_biguint<P*Q> la = (sc_biguint<P*Q>) la1; // cast
	sc_biguint<P*Q> lb = (sc_biguint<P*Q>) lb1;

	sc_biguint<P> signa = VECTOR_SIGN<P,Q> (la);
	sc_biguint<P> signb = VECTOR_SIGN<P,Q> (lb);
	sc_biguint<P> sig = VECTOR_XOR<P> (signb, signa);

	sc_biguint<P> result = VECTOR_AND<P> (sig, fb); // Partial Sum computation

	return result;
}

template <int P, int Q>
sc_biguint<P> G_simplified_SM(sc_bigint<P*Q> la1, sc_bigint<P*Q> lb1, sc_biguint<P> sa, sc_biguint<P> fb)
{
#pragma HLS INLINE
	sc_biguint<P*Q> la = (sc_biguint<P*Q>) la1; // cast
	sc_biguint<P*Q> lb = (sc_biguint<P*Q>) lb1;

	sc_biguint<P> sla = VECTOR_SIGN_SM<P,Q> ( la );
	//sc_biguint<P> n_sa = (sc_biguint<P>) VECTOR_NOT_SM<P,1> ( (sc_biguint<P>) sa );
	sc_biguint<P> n_sa =  sa;
	sc_biguint<P> sigla = VECTOR_XOR<P>( sla, n_sa );
	sc_biguint<P> siglb = VECTOR_SIGN_SM<P,Q> (lb);

	sc_biguint<P> siga = VECTOR_AND<P>( sigla, fb );
	sc_biguint<P> sigb = VECTOR_AND<P>( siglb, fb );

	sc_biguint<P*(Q-1)> absla = VECTOR_ABS_SM<P,Q> (la);
	sc_biguint<P*(Q-1)> abslb = VECTOR_ABS_SM<P,Q> (lb);

	sc_biguint<P> is_min = VECTOR_IS_MIN_SM<P,Q> ( absla, abslb ); // if absla < abslb return 1

	sc_biguint<P> SigSomme = (sc_biguint<P>) VECTOR_MUX_SM<P,1> ((sc_biguint<P>) sigb, (sc_biguint<P>) siga, is_min);

	return SigSomme;
}

//*************************************************************************//
//** 						POLAR CODE FUNCTIONS						 **//
//*************************************************************************//

template <int P, int Q>  // P : Parallelism , Q : Quantification
sc_bigint<P*Q> Function_F(sc_bigint<P*Q> la, sc_bigint<P*Q> lb){
	#pragma HLS INLINE

#if defined CA2
	return F_function_C2<P,Q> (la,lb);
#elif defined SIGMAG
	return F_function_SM<P,Q> (la,lb);
#endif
}

template <int P, int Q>
sc_bigint<P*Q> Function_G(sc_bigint<P*Q> la, sc_bigint<P*Q> lb, sc_biguint<P> sa){
	#pragma HLS INLINE

#if defined CA2
	return G_function_C2<P,Q> (la,lb,sa);
#elif defined SIGMAG
	return G_function_SM<P,Q> (la,lb,sa);
#endif
}

template <int P, int Q>
sc_bigint<P*(Q+1)> Function_G_ext(sc_bigint<P*Q> la, sc_bigint<P*Q> lb, sc_biguint<P> sa){ // No saturation
	#pragma HLS INLINE

#if defined CA2
	return G_extended_C2<P,Q> (la,lb,sa);
#elif defined SIGMAG
	return G_extended_SM<P,Q> (la,lb,sa);
#endif
}

template <int P, int Q>
sc_biguint<P> Function_F_simp(sc_bigint<P*Q> la, sc_bigint<P*Q> lb, sc_biguint<P> fb){
	#pragma HLS INLINE

#if defined CA2
	return F_simplified_C2<P,Q> (la,lb,fb);
#elif defined SIGMAG
	return F_simplified_SM<P,Q> (la,lb,fb);
#endif
}

template <int P, int Q>
sc_biguint<P> Function_G_simp(sc_bigint<P*Q> la, sc_bigint<P*Q> lb, sc_biguint<P> sa, sc_biguint<P> fb){
	#pragma HLS INLINE

#if defined CA2
	return G_simplified_C2<P,Q> (la,lb,sa,fb);
#elif defined SIGMAG
	return G_simplified_SM<P,Q> (la,lb,sa,fb);
#endif
}

template <int P>
sc_biguint<P> Function_H(sc_biguint<P> sa, sc_biguint<P> sb){
	#pragma HLS INLINE

	return VECTOR_XOR<P> (sa,sb);
}


//*************************************************************************//
//** 					SPECIALIZED POLAR DECODERS						 **//
//*************************************************************************//

template <int Q>
sc_biguint<1> Spec_P1 ( sc_bigint<Q> llr, sc_biguint<1> fb )
{
#pragma HLS INLINE

	sc_biguint<1> sign = VECTOR_SIGN<1,Q>(llr);

	sc_biguint<1> sa = VECTOR_AND<1>(sign,fb);

	return sa;
}

template <int Q>
sc_biguint<2> Spec_P2 ( sc_bigint<2*Q> llr, sc_biguint<2> fb )
{
#pragma HLS INLINE

	sc_bigint<Q> la = llr.range(Q-1,0);
	sc_bigint<Q> lb = llr.range(2*Q-1,Q);

	sc_biguint<1> fba = (sc_biguint<1>) fb[0];
	sc_biguint<1> sa1 = Function_F_simp<1,Q> (la,lb,fba);

	sc_biguint<1> fbb = (sc_biguint<1>) fb[1];
	sc_biguint<1> sb1 = Function_G_simp<1,Q> (la,lb,sa1,fbb);

	sc_biguint<1> sa = VECTOR_XOR<1>(sa1,sb1);
	sc_biguint<1> sb = sb1;

	return (sb , sa);
}

template <int Q>
sc_biguint<4> Spec_P4 ( sc_bigint<4*Q> llr, sc_biguint<4> fb )
{
#pragma HLS INLINE

	sc_bigint<2*Q> la = llr.range((2*Q)-1,0);
	sc_bigint<2*Q> lb = llr.range(2*(2*Q)-1,(2*Q));

	// F
	sc_bigint<2*Q> la1 = Function_F<2,Q> (la,lb);

	// First P2
	sc_biguint<2> sa1 = Spec_P2<Q> ( la1, fb.range(1,0) );

	// G
	sc_bigint<2*Q> lb1 = Function_G<2,Q> (la,lb,sa1);

	// Second P2
	sc_biguint<2> sb1 = Spec_P2<Q> ( lb1 , fb.range(3,2) );

	// XOR
	sc_biguint<2> sa = VECTOR_XOR<2>(sa1,sb1);
	sc_biguint<2> sb = sb1;

	return ( sb, sa );
}

template <int Q>
sc_biguint<4> Spec_P4_ext ( sc_bigint<4*Q> llr, sc_biguint<4> fb )
{
#pragma HLS INLINE

	sc_bigint<2*Q> la = llr.range((2*Q)-1,0);
	sc_bigint<2*Q> lb = llr.range(2*(2*Q)-1,(2*Q));

	// F
	sc_bigint<2*Q> la1 = Function_F<2,Q> (la,lb);

	// First P2
	sc_biguint<2> sa1 = Spec_P2<Q> ( la1, fb.range(1,0) );

	// G
	sc_bigint<2*(Q+1)> lb1 = Function_G_ext<2,Q> (la,lb,sa1);

	// Second P2
	sc_biguint<2> sb1 = Spec_P2<Q+1> ( lb1 , fb.range(3,2) );

	// XOR
	sc_biguint<2> sa = VECTOR_XOR<2>(sa1,sb1);
	sc_biguint<2> sb = sb1;

	return ( sb, sa );
}

template <int Q>
sc_biguint<8> Spec_P8 ( sc_bigint<8*Q> llr, sc_biguint<8> fb )
{
#pragma HLS INLINE

	sc_bigint<4*Q> la = llr.range((4*Q)-1,0);
	sc_bigint<4*Q> lb = llr.range(2*(4*Q)-1,(4*Q));

	// F
	sc_bigint<4*Q> la1 = Function_F<4,Q> (la,lb);

	// First P4
	sc_biguint<4> sa1 = Spec_P4<Q> ( la1, fb.range(3,0) );

	// G
	sc_bigint<4*Q> lb1 = Function_G<4,Q> (la,lb,sa1);

	// Second P2
	sc_biguint<4> sb1 = Spec_P4<Q> ( lb1 , fb.range(7,4) );

	// XOR
	sc_biguint<4> sa = VECTOR_XOR<4>(sa1,sb1);
	sc_biguint<4> sb = sb1;

	return ( sb, sa );
}

template <int Q>
sc_biguint<8> Spec_P8_ext ( sc_bigint<8*Q> llr, sc_biguint<8> fb )
{
#pragma HLS INLINE

	sc_bigint<4*Q> la = llr.range((4*Q)-1,0);
	sc_bigint<4*Q> lb = llr.range(2*(4*Q)-1,(4*Q));

	// F
	sc_bigint<4*Q> la1 = Function_F<4,Q> (la,lb);

	// First P2
	sc_biguint<4> sa1 = Spec_P4_ext<Q> ( la1, fb.range(3,0) );

	// G
	sc_bigint<4*(Q+1)> lb1 = Function_G_ext<4,Q> (la,lb,sa1);

	// Second P2
	sc_biguint<4> sb1 = Spec_P4_ext<Q+1> ( lb1 , fb.range(7,4) );

	// XOR
	sc_biguint<4> sa = VECTOR_XOR<4>(sa1,sb1);
	sc_biguint<4> sb = sb1;

	return ( sb, sa );
}

template <int Q>
sc_biguint<16> Spec_P16 ( sc_bigint<16*Q> llr, sc_biguint<16> fb )
{
#pragma HLS INLINE

	sc_bigint<8*Q> la = llr.range((8*Q)-1,0);
	sc_bigint<8*Q> lb = llr.range(2*(8*Q)-1,(8*Q));

	// F
	sc_bigint<8*Q> la1 = Function_F<8,Q> (la,lb);

	// First P4
	sc_biguint<8> sa1 = Spec_P8<Q> ( la1, fb.range(7,0) );

	// G
	sc_bigint<8*Q> lb1 = Function_G<8,Q> (la,lb,sa1);

	// Second P2
	sc_biguint<8> sb1 = Spec_P8<Q> ( lb1 , fb.range(15,8) );

	// XOR
	sc_biguint<8> sa = VECTOR_XOR<8>(sa1,sb1);
	sc_biguint<8> sb = sb1;

	return ( sb, sa );
}

template <int Q>
sc_biguint<16> Spec_P16_ext ( sc_bigint<16*Q> llr, sc_biguint<16> fb )
{
#pragma HLS INLINE

	sc_bigint<8*Q> la = llr.range((8*Q)-1,0);
	sc_bigint<8*Q> lb = llr.range(2*(8*Q)-1,(8*Q));

	// F
	sc_bigint<8*Q> la1 = Function_F<8,Q> (la,lb);

	// First P2
	sc_biguint<8> sa1 = Spec_P8_ext<Q> ( la1, fb.range(7,0) );

	// G
	sc_bigint<8*(Q+1)> lb1 = Function_G_ext<8,Q> (la,lb,sa1);

	// Second P2
	sc_biguint<8> sb1 = Spec_P8_ext<Q+1> ( lb1 , fb.range(15,8) );

	// XOR
	sc_biguint<8> sa = VECTOR_XOR<8>(sa1,sb1);
	sc_biguint<8> sb = sb1;

	return ( sb, sa );
}

template <int Q>
sc_biguint< 32 > Spec_P32 ( sc_bigint< 32 *Q> llr, sc_biguint< 32 > fb )
{
#pragma HLS INLINE

	sc_bigint< 16 *Q> la = llr.range(( 16 *Q)-1,0);
	sc_bigint< 16 *Q> lb = llr.range(2*( 16 *Q)-1,( 16 *Q));

	// F
	sc_bigint< 16 *Q> la1 = Function_F< 16 ,Q> (la,lb);

	// First P4
	sc_biguint< 16 > sa1 = Spec_P16<Q> ( la1, fb.range( 15, 0 ) );

	// G
	sc_bigint< 16 *Q> lb1 = Function_G< 16 ,Q> (la,lb,sa1);

	// Second P2
	sc_biguint< 16 > sb1 = Spec_P16<Q> ( lb1 , fb.range( 31, 16 ) );

	// XOR
	sc_biguint< 16 > sa = VECTOR_XOR< 16 >(sa1,sb1);
	sc_biguint< 16 > sb = sb1;

	return ( sb, sa );
}

template <int Q>
sc_biguint< 32 > Spec_P32_ext ( sc_bigint< 32 *Q> llr, sc_biguint< 32 > fb )
{
#pragma HLS INLINE

	sc_bigint< 16 *Q> la = llr.range(( 16 *Q)-1,0);
	sc_bigint< 16 *Q> lb = llr.range(2*( 16 *Q)-1,( 16 *Q));

	// F
	sc_bigint< 16 *Q> la1 = Function_F< 16 ,Q> (la,lb);

	// First P2
	sc_biguint< 16 > sa1 = Spec_P16_ext<Q> ( la1, fb.range( 15, 0 ) );

	// G
	sc_bigint< 16 *(Q+1)> lb1 = Function_G_ext< 16 ,Q> (la,lb,sa1);

	// Second P2
	sc_biguint< 16 > sb1 = Spec_P16_ext<Q+1> ( lb1 , fb.range( 31, 16 ) );

	// XOR
	sc_biguint< 16 > sa = VECTOR_XOR< 16 >(sa1,sb1);
	sc_biguint< 16 > sb = sb1;

	return ( sb, sa );
}

template <int Q>
sc_biguint< 64 > Spec_P64 ( sc_bigint< 64 *Q> llr, sc_biguint< 64 > fb )
{
#pragma HLS INLINE

	sc_bigint< 32 *Q> la = llr.range(( 32 *Q)-1,0);
	sc_bigint< 32 *Q> lb = llr.range(2*( 32 *Q)-1,( 32 *Q));

	// F
	sc_bigint< 32 *Q> la1 = Function_F< 32 ,Q> (la,lb);

	// First P4
	sc_biguint< 32 > sa1 = Spec_P32<Q> ( la1, fb.range( 31, 0 ) );

	// G
	sc_bigint< 32 *Q> lb1 = Function_G< 32 ,Q> (la,lb,sa1);

	// Second P2
	sc_biguint< 32 > sb1 = Spec_P32<Q> ( lb1 , fb.range( 63, 32 ) );

	// XOR
	sc_biguint< 32 > sa = VECTOR_XOR< 32 >(sa1,sb1);
	sc_biguint< 32 > sb = sb1;

	return ( sb, sa );
}

template <int Q>
sc_biguint< 64 > Spec_P64_ext ( sc_bigint< 64 *Q> llr, sc_biguint< 64 > fb )
{
#pragma HLS INLINE

	sc_bigint< 32 *Q> la = llr.range(( 32 *Q)-1,0);
	sc_bigint< 32 *Q> lb = llr.range(2*( 32 *Q)-1,( 32 *Q));

	// F
	sc_bigint< 32 *Q> la1 = Function_F< 32 ,Q> (la,lb);

	// First P2
	sc_biguint< 32 > sa1 = Spec_P32_ext<Q> ( la1, fb.range( 31, 0 ) );

	// G
	sc_bigint< 32 *(Q+1)> lb1 = Function_G_ext< 32 ,Q> (la,lb,sa1);

	// Second P2
	sc_biguint< 32 > sb1 = Spec_P32_ext<Q+1> ( lb1 , fb.range( 63, 32 ) );

	// XOR
	sc_biguint< 32 > sa = VECTOR_XOR< 32 >(sa1,sb1);
	sc_biguint< 32 > sb = sb1;

	return ( sb, sa );
}

template <int Q>
sc_biguint< 128 > Spec_P128 ( sc_bigint< 128 *Q> llr, sc_biguint< 128 > fb )
{
#pragma HLS INLINE

	sc_bigint< 64 *Q> la = llr.range(( 64 *Q)-1,0);
	sc_bigint< 64 *Q> lb = llr.range(2*( 64 *Q)-1,( 64 *Q));

	// F
	sc_bigint< 64 *Q> la1 = Function_F< 64 ,Q> (la,lb);

	// First P4
	sc_biguint< 64 > sa1 = Spec_P64<Q> ( la1, fb.range( 63, 0 ) );

	// G
	sc_bigint< 64 *Q> lb1 = Function_G< 64 ,Q> (la,lb,sa1);

	// Second P2
	sc_biguint< 64 > sb1 = Spec_P64<Q> ( lb1 , fb.range( 127, 64 ) );

	// XOR
	sc_biguint< 64 > sa = VECTOR_XOR< 64 >(sa1,sb1);
	sc_biguint< 64 > sb = sb1;

	return ( sb, sa );
}

template <int Q>
sc_biguint< 128 > Spec_P128_ext ( sc_bigint< 128 *Q> llr, sc_biguint< 128 > fb )
{
#pragma HLS INLINE

	sc_bigint< 64 *Q> la = llr.range(( 64 *Q)-1,0);
	sc_bigint< 64 *Q> lb = llr.range(2*( 64 *Q)-1,( 64 *Q));

	// F
	sc_bigint< 64 *Q> la1 = Function_F< 64 ,Q> (la,lb);

	// First P2
	sc_biguint< 64 > sa1 = Spec_P64_ext<Q> ( la1, fb.range( 63, 0 ) );

	// G
	sc_bigint< 64 *(Q+1)> lb1 = Function_G_ext< 64 ,Q> (la,lb,sa1);

	// Second P2
	sc_biguint< 64 > sb1 = Spec_P64_ext<Q+1> ( lb1 , fb.range( 127, 64 ) );

	// XOR
	sc_biguint< 64 > sa = VECTOR_XOR< 64 >(sa1,sb1);
	sc_biguint< 64 > sb = sb1;

	return ( sb, sa );
}

template <int Q>
sc_biguint< 256 > Spec_P256 ( sc_bigint< 256 *Q> llr, sc_biguint< 256 > fb )
{
#pragma HLS INLINE

	sc_bigint< 128 *Q> la = llr.range(( 128 *Q)-1,0);
	sc_bigint< 128 *Q> lb = llr.range(2*( 128 *Q)-1,( 128 *Q));

	// F
	sc_bigint< 128 *Q> la1 = Function_F< 128 ,Q> (la,lb);

	// First P4
	sc_biguint< 128 > sa1 = Spec_P128<Q> ( la1, fb.range( 127, 0 ) );

	// G
	sc_bigint< 128 *Q> lb1 = Function_G< 128 ,Q> (la,lb,sa1);

	// Second P2
	sc_biguint< 128 > sb1 = Spec_P128<Q> ( lb1 , fb.range( 255, 128 ) );

	// XOR
	sc_biguint< 128 > sa = VECTOR_XOR< 128 >(sa1,sb1);
	sc_biguint< 128 > sb = sb1;

	return ( sb, sa );
}

template <int Q>
sc_biguint< 256 > Spec_P256_ext ( sc_bigint< 256 *Q> llr, sc_biguint< 256 > fb )
{
#pragma HLS INLINE

	sc_bigint< 128 *Q> la = llr.range(( 128 *Q)-1,0);
	sc_bigint< 128 *Q> lb = llr.range(2*( 128 *Q)-1,( 128 *Q));

	// F
	sc_bigint< 128 *Q> la1 = Function_F< 128 ,Q> (la,lb);

	// First P2
	sc_biguint< 128 > sa1 = Spec_P128_ext<Q> ( la1, fb.range( 127, 0 ) );

	// G
	sc_bigint< 128 *(Q+1)> lb1 = Function_G_ext< 128 ,Q> (la,lb,sa1);

	// Second P2
	sc_biguint< 128 > sb1 = Spec_P128_ext<Q+1> ( lb1 , fb.range( 255, 128 ) );

	// XOR
	sc_biguint< 128 > sa = VECTOR_XOR< 128 >(sa1,sb1);
	sc_biguint< 128 > sb = sb1;

	return ( sb, sa );
}

//*************************************************************************//
//** 			GLOBAL SPECIALIZED POLAR DECODERS						 **//
//*************************************************************************//
template <int Q>
sc_biguint< 1 > Spec_PolarDec_1 ( sc_bigint< Q> llr, sc_biguint< 1 > fb )
{
#pragma HLS INLINE

	return Spec_P1< Q > (llr, fb);
}


template <int Q>
sc_biguint< 2 > Spec_PolarDec_2 ( sc_bigint< 2 *Q> llr, sc_biguint< 2 > fb )
{
#pragma HLS INLINE

	return Spec_P2< Q > (llr, fb);
}

template <int Q>
sc_biguint< 4 > Spec_PolarDec_4 ( sc_bigint< 4 *Q> llr, sc_biguint< 4 > fb )
{
#pragma HLS INLINE

#if EXTENDED == 1
	return Spec_P4_ext< Q > (llr, fb);
#else
	return Spec_P4< Q > (llr, fb);
#endif
}

template <int Q>
sc_biguint< 8 > Spec_PolarDec_8 ( sc_bigint< 8 *Q> llr, sc_biguint< 8 > fb )
{
#pragma HLS INLINE

#if EXTENDED == 1
	return Spec_P8_ext< Q > (llr, fb);
#else
	return Spec_P8< Q > (llr, fb);
#endif
}

template <int Q>
sc_biguint< 16 > Spec_PolarDec_16 ( sc_bigint< 16 *Q> llr, sc_biguint< 16 > fb )
{
#pragma HLS INLINE

#if EXTENDED == 1
	return Spec_P16_ext< Q > (llr, fb);
#else
	return Spec_P16< Q > (llr, fb);
#endif
}

template <int Q>
sc_biguint< 32 > Spec_PolarDec_32 ( sc_bigint< 32 *Q> llr, sc_biguint< 32 > fb )
{
#pragma HLS INLINE

#if EXTENDED == 1
	return Spec_P32_ext< Q > (llr, fb);
#else
	return Spec_P32< Q > (llr, fb);
#endif
}

template <int Q>
sc_biguint< 64 > Spec_PolarDec_64 ( sc_bigint< 64 *Q> llr, sc_biguint< 64 > fb )
{
#pragma HLS INLINE

#if EXTENDED == 1
	return Spec_P64_ext< Q > (llr, fb);
#else
	return Spec_P64< Q > (llr, fb);
#endif
}

template <int Q>
sc_biguint< 128 > Spec_PolarDec_128 ( sc_bigint< 128 *Q> llr, sc_biguint< 128 > fb )
{
#pragma HLS INLINE

#if EXTENDED == 1
	return Spec_P128_ext< Q > (llr, fb);
#else
	return Spec_P128< Q > (llr, fb);
#endif
}

template <int Q>
sc_biguint< 256 > Spec_PolarDec_256 ( sc_bigint< 256 *Q> llr, sc_biguint< 256 > fb )
{
#pragma HLS INLINE

#if EXTENDED == 1
	return Spec_P256_ext< Q > (llr, fb);
#else
	return Spec_P256< Q > (llr, fb);
#endif
}

//*************************************************************************//
//** 						REP NODE DECODERS							 **//
//*************************************************************************//

/////////////// CA2 ////////////////

template <int Q>
sc_biguint< 2 > REP_2_CA2 (sc_bigint< 2 * Q> llr)
{
#pragma HLS INLINE
	sc_bigint< Q 	 > la  = llr.range( Q - 1, 0 );
	sc_bigint< Q 	 > lb  = llr.range( 2 * Q - 1, Q );
	sc_bigint< Q + 1 > sum = VECTOR_ADD_NOSAT < 1, Q > ( la, lb );

	sc_biguint< 1 > sig = qsign < Q + 1 > (sum);

	sc_biguint< 2 > res = (sig, sig);
	return res;
}

template <int Q>
sc_biguint< 4 > REP_4_CA2 (sc_bigint< 4 * Q> llr)
{
#pragma HLS INLINE
	sc_bigint< 2 * Q 	   > la  = llr.range( (2 * Q) - 1, 0 );
	sc_bigint< 2 * Q 	   > lb  = llr.range( (4 * Q) - 1, (2 * Q) );
	sc_bigint< 2 * (Q + 1) > sum = VECTOR_ADD_NOSAT < 2, Q > ( la, lb );

	sc_biguint< 2 > sig = REP_2_CA2 < Q + 1 > ( sum );

	sc_biguint< 4 > res = (sig, sig);

	return res;
}

template <int Q>
sc_biguint< 8 > REP_8_CA2 (sc_bigint< 8 * Q> llr)
{
#pragma HLS INLINE
	sc_bigint< 4 * Q 	   > la  = llr.range( (4 * Q) - 1, 0 );
	sc_bigint< 4 * Q 	   > lb  = llr.range( (8 * Q) - 1, (4 * Q) );
	sc_bigint< 4 * (Q + 1) > sum = VECTOR_ADD_NOSAT < 4, Q > ( la, lb );

	sc_biguint< 4 > sig = REP_4_CA2 < Q + 1 > ( sum );

	sc_biguint< 8 > res = (sig, sig);

	return res;
}

template <int Q>
sc_biguint< 16 > REP_16_CA2 (sc_bigint< 16 * Q> llr)
{
#pragma HLS INLINE
	sc_bigint< 8 *  Q      > la  = llr.range( (8 * Q) - 1, 0 );
	sc_bigint< 8 *  Q      > lb  = llr.range( (16 * Q) - 1, (8 * Q) );
	sc_bigint< 8 * (Q + 1) > sum = VECTOR_ADD_NOSAT < 8, Q > ( la, lb );

	sc_biguint< 8 > sig = REP_8_CA2 < Q + 1 > ( sum );

	sc_biguint< 16 > res = (sig, sig);

	return res;
}

template <int Q>
sc_biguint< 32 > REP_32_CA2 (sc_bigint< 32 * Q> llr)
{
#pragma HLS INLINE
	sc_bigint< 16 *  Q      > la  = llr.range( (16 * Q) - 1, 0 );
	sc_bigint< 16 *  Q      > lb  = llr.range( (32 * Q) - 1, (16 * Q) );
	sc_bigint< 16 * (Q + 1) > sum = VECTOR_ADD_NOSAT < 16, Q > ( la, lb );

	sc_biguint< 16 > sig = REP_16_CA2 < Q + 1 > ( sum );

	sc_biguint< 32 > res = (sig, sig);

	return res;
}

template <int Q>
sc_biguint< 64 > REP_64_CA2 (sc_bigint< 64 * Q> llr)
{
#pragma HLS INLINE
	sc_bigint< 32 *  Q      > la  = llr.range( (32 * Q) - 1, 0 );
	sc_bigint< 32 *  Q      > lb  = llr.range( (64 * Q) - 1, (32 * Q) );
	sc_bigint< 32 * (Q + 1) > sum = VECTOR_ADD_NOSAT < 32, Q > ( la, lb );

	sc_biguint< 32 > sig = REP_32_CA2 < Q + 1 > ( sum );

	sc_biguint< 64 > res = (sig, sig);

	return res;
}

template <int Q>
sc_biguint< 128 > REP_128_CA2 (sc_bigint< 128 * Q> llr)
{
#pragma HLS INLINE
	sc_bigint< 64 *  Q      > la  = llr.range( (64 * Q) - 1, 0 );
	sc_bigint< 64 *  Q      > lb  = llr.range( (128 * Q) - 1, (64 * Q) );
	sc_bigint< 64 * (Q + 1) > sum = VECTOR_ADD_NOSAT < 64, Q > ( la, lb );

	sc_biguint< 64 > sig = REP_64_CA2 < Q + 1 > ( sum );

	sc_biguint< 128 > res = (sig, sig);

	return res;
}

template <int Q>
sc_biguint< 256 > REP_256_CA2 (sc_bigint< 256 * Q> llr)
{
#pragma HLS INLINE
	sc_bigint< 128 *  Q      > la  = llr.range( (128 * Q) - 1, 0 );
	sc_bigint< 128 *  Q      > lb  = llr.range( (256 * Q) - 1, (128 * Q) );
	sc_bigint< 128 * (Q + 1) > sum = VECTOR_ADD_NOSAT < 128, Q > ( la, lb );

	sc_biguint< 128 > sig = REP_128_CA2 < Q + 1 > ( sum );

	sc_biguint< 256 > res = (sig, sig);

	return res;
}

/////////////// SIGMAG ////////////////

template <int Q>
sc_biguint< 2 > REP_2_SM (sc_biguint< 2 * Q> llr)
{
#pragma HLS INLINE
	sc_biguint< Q 	 > la  = llr.range( Q - 1, 0 );
	sc_biguint< Q 	 > lb  = llr.range( 2 * Q - 1, Q );

	sc_biguint<1> siga = (sc_biguint<1>) la[Q - 1];
	sc_biguint<1> sigb = (sc_biguint<1>) lb[Q - 1];

	sc_biguint<Q-1> absla = la.range(Q-2 , 0);
	sc_biguint<Q-1> abslb = lb.range(Q-2 , 0);

	sc_biguint<1> is_min = ( absla < abslb ) ? 1 : 0; // if absla < abslb return 1

	sc_biguint<1> sig = (sc_uint<1>)is_min ? sigb : siga;

	sc_biguint< 2 > res = (sig, sig);
	return res;
}

template <int Q>
sc_biguint< 4 > REP_4_SM (sc_biguint< 4 * Q> llr)
{
#pragma HLS INLINE
	sc_biguint< 2 * Q 	   > la  = llr.range( (2 * Q) - 1, 0 );
	sc_biguint< 2 * Q 	   > lb  = llr.range( (4 * Q) - 1, (2 * Q) );

	sc_biguint< 2 * (Q + 1) > sum = VECTOR_FULL_ADDER_SM < 2, Q > ( la, lb );

	sc_biguint< 2 > sig = REP_2_SM < Q + 1 > ( sum );

	sc_biguint< 4 > res = (sig, sig);

	return res;
}

template <int Q>
sc_biguint< 8 > REP_8_SM (sc_biguint< 8 * Q> llr)
{
#pragma HLS INLINE
	sc_biguint< 4 * Q 	   > la  = llr.range( (4 * Q) - 1, 0 );
	sc_biguint< 4 * Q 	   > lb  = llr.range( (8 * Q) - 1, (4 * Q) );

	sc_biguint< 4 * (Q + 1) > sum = VECTOR_FULL_ADDER_SM < 4, Q > ( la, lb );

	sc_biguint< 4 > sig = REP_4_SM < Q + 1 > ( sum );

	sc_biguint< 8 > res = (sig, sig);

	return res;
}

template <int Q>
sc_biguint< 16 > REP_16_SM (sc_biguint< 16 * Q> llr)
{
#pragma HLS INLINE
	sc_biguint< 8 * Q 	   > la  = llr.range( (8 * Q) - 1, 0 );
	sc_biguint< 8 * Q 	   > lb  = llr.range( (16 * Q) - 1, (8 * Q) );

	sc_biguint< 8 * (Q + 1) > sum = VECTOR_FULL_ADDER_SM < 8, Q > ( la, lb );

	sc_biguint< 8 > sig = REP_8_SM < Q + 1 > ( sum );

	sc_biguint< 16 > res = (sig, sig);

	return res;
}

template <int Q>
sc_biguint< 32 > REP_32_SM (sc_biguint< 32 * Q> llr)
{
#pragma HLS INLINE
	sc_biguint< 16 * Q 	   > la  = llr.range( (16 * Q) - 1, 0 );
	sc_biguint< 16 * Q 	   > lb  = llr.range( (32 * Q) - 1, (16 * Q) );

	sc_biguint< 16 * (Q + 1) > sum = VECTOR_FULL_ADDER_SM < 16, Q > ( la, lb );

	sc_biguint< 16 > sig = REP_16_SM < Q + 1 > ( sum );

	sc_biguint< 32 > res = (sig, sig);

	return res;
}

template <int Q>
sc_biguint< 64 > REP_64_SM (sc_biguint< 64 * Q> llr)
{
#pragma HLS INLINE
	sc_biguint< 32 * Q 	   > la  = llr.range( (32 * Q) - 1, 0 );
	sc_biguint< 32 * Q 	   > lb  = llr.range( (64 * Q) - 1, (32 * Q) );

	sc_biguint< 32 * (Q + 1) > sum = VECTOR_FULL_ADDER_SM < 32, Q > ( la, lb );

	sc_biguint< 32 > sig = REP_32_SM < Q + 1 > ( sum );

	sc_biguint< 64 > res = (sig, sig);

	return res;
}

template <int Q>
sc_biguint< 128 > REP_128_SM (sc_biguint< 128 * Q> llr)
{
#pragma HLS INLINE
	sc_biguint< 64 * Q 	   > la  = llr.range( (64 * Q) - 1, 0 );
	sc_biguint< 64 * Q 	   > lb  = llr.range( (128 * Q) - 1, (64 * Q) );

	sc_biguint< 64 * (Q + 1) > sum = VECTOR_FULL_ADDER_SM < 64, Q > ( la, lb );

	sc_biguint< 64 > sig = REP_64_SM < Q + 1 > ( sum );

	sc_biguint< 128 > res = (sig, sig);

	return res;
}

template <int Q>
sc_biguint< 256 > REP_256_SM (sc_biguint< 256 * Q> llr)
{
#pragma HLS INLINE
	sc_biguint< 128 * Q 	   > la  = llr.range( (128 * Q) - 1, 0 );
	sc_biguint< 128 * Q 	   > lb  = llr.range( (256 * Q) - 1, (128 * Q) );

	sc_biguint< 128 * (Q + 1) > sum = VECTOR_FULL_ADDER_SM < 128, Q > ( la, lb );

	sc_biguint< 128 > sig = REP_128_SM < Q + 1 > ( sum );

	sc_biguint< 256 > res = (sig, sig);

	return res;
}

/////////////// GLOBAL REP Node ////////////////

template <int Q>
sc_biguint< 2 > REP_2_Node (sc_bigint< 2 * Q> llr)
{
#pragma HLS INLINE

#if defined CA2
	return REP_2_CA2< Q >(llr);
#elif defined SIGMAG
	return REP_2_SM< Q >( (sc_biguint< 2 * Q>) llr);
#endif
}

template <int Q>
sc_biguint< 4 > REP_4_Node (sc_bigint< 4 * Q> llr)
{
#pragma HLS INLINE

#if defined CA2
	return REP_4_CA2< Q >(llr);
#elif defined SIGMAG
	return REP_4_SM< Q >( (sc_biguint< 4 * Q>) llr);
#endif
}

template <int Q>
sc_biguint< 8 > REP_8_Node (sc_bigint< 8 * Q> llr)
{
#pragma HLS INLINE

#if defined CA2
	return REP_8_CA2< Q >(llr);
#elif defined SIGMAG
	return REP_8_SM< Q >( (sc_biguint< 8 * Q>) llr);
#endif
}

template <int Q>
sc_biguint< 16 > REP_16_Node (sc_bigint< 16 * Q> llr)
{
#pragma HLS INLINE

#if defined CA2
	return REP_16_CA2< Q >(llr);
#elif defined SIGMAG
	return REP_16_SM< Q >( (sc_biguint< 16 * Q>) llr);
#endif
}

template <int Q>
sc_biguint< 32 > REP_32_Node (sc_bigint< 32 * Q> llr)
{
#pragma HLS INLINE

#if defined CA2
	return REP_32_CA2< Q >(llr);
#elif defined SIGMAG
	return REP_32_SM< Q >( (sc_biguint< 32 * Q>) llr);
#endif
}

template <int Q>
sc_biguint< 64 > REP_64_Node (sc_bigint< 64 * Q> llr)
{
#pragma HLS INLINE

#if defined CA2
	return REP_64_CA2< Q >(llr);
#elif defined SIGMAG
	return REP_64_SM< Q >( (sc_biguint< 64 * Q>) llr);
#endif
}

template <int Q>
sc_biguint< 128 > REP_128_Node (sc_bigint< 128 * Q> llr)
{
#pragma HLS INLINE

#if defined CA2
	return REP_128_CA2< Q >(llr);
#elif defined SIGMAG
	return REP_128_SM< Q >( (sc_biguint< 128 * Q>) llr);
#endif
}

template <int Q>
sc_biguint< 256 > REP_256_Node (sc_bigint< 256 * Q> llr)
{
#pragma HLS INLINE

#if defined CA2
	return REP_256_CA2< Q >(llr);
#elif defined SIGMAG
	return REP_256_SM< Q >( (sc_biguint< 256 * Q>) llr);
#endif
}
//*************************************************************************//
//** 					REP / REP2 NODE DECODERS						 **//
//*************************************************************************//

/////////////// CA2 ////////////////

template <int Q>
sc_biguint< 2 > REP_REP2_2_CA2 (sc_bigint< 2 * Q> llr, bool sel)  // Sel 0 : REP, Sel 1 : REP2
{
#pragma HLS INLINE
	sc_bigint< Q 	 > la  = llr.range( Q - 1, 0 );
	sc_bigint< Q 	 > lb  = llr.range( 2 * Q - 1, Q );

	sc_biguint< 1 > sla = qsign < Q > (la);
	sc_biguint< 1 > slb = qsign < Q > (lb);
	sc_biguint< 2 > rep2 = (slb, sla);

	sc_bigint< Q + 1 > sum = VECTOR_ADD_NOSAT < 1, Q > ( la, lb );
	sc_biguint< 1 > sign_sum = qsign < Q + 1 > (sum);
	sc_biguint< 2 > rep = (sign_sum, sign_sum);

	if(sel == 0)
		return rep;
	else
		return rep2;
}

template <int Q>
sc_biguint< 4 > REP_REP2_4_CA2 (sc_bigint< 4 * Q> llr, bool sel)
{
#pragma HLS INLINE
	sc_bigint< 2 * Q 	   > la  = llr.range( (2 * Q) - 1, 0 );
	sc_bigint< 2 * Q 	   > lb  = llr.range( (4 * Q) - 1, (2 * Q) );
	sc_bigint< 2 * (Q + 1) > sum = VECTOR_ADD_NOSAT < 2, Q > ( la, lb );

	sc_biguint< 2 > sig = REP_REP2_2_CA2 < Q + 1 > ( sum, sel );

	sc_biguint< 4 > res = (sig, sig);
	return res;
}

template <int Q>
sc_biguint< 8 > REP_REP2_8_CA2 (sc_bigint< 8 * Q> llr, bool sel)
{
#pragma HLS INLINE
	sc_bigint< 4 * Q 	   > la  = llr.range( (4 * Q) - 1, 0 );
	sc_bigint< 4 * Q 	   > lb  = llr.range( (8 * Q) - 1, (4 * Q) );
	sc_bigint< 4 * (Q + 1) > sum = VECTOR_ADD_NOSAT < 4, Q > ( la, lb );

	sc_biguint< 4 > sig = REP_REP2_4_CA2 < Q + 1 > ( sum, sel );

	sc_biguint< 8 > res = (sig, sig);
	return res;
}

template <int Q>
sc_biguint< 16 > REP_REP2_16_CA2 (sc_bigint< 16 * Q> llr, bool sel)
{
#pragma HLS INLINE
	sc_bigint< 8 *  Q      > la  = llr.range( (8 * Q) - 1, 0 );
	sc_bigint< 8 *  Q      > lb  = llr.range( (16 * Q) - 1, (8 * Q) );
	sc_bigint< 8 * (Q + 1) > sum = VECTOR_ADD_NOSAT < 8, Q > ( la, lb );

	sc_biguint< 8 > sig = REP_REP2_8_CA2 < Q + 1 > ( sum, sel );

	sc_biguint< 16 > res = (sig, sig);
	return res;
}

template <int Q>
sc_biguint< 32 > REP_REP2_32_CA2 (sc_bigint< 32 * Q> llr, bool sel)
{
#pragma HLS INLINE
	sc_bigint< 16 *  Q      > la  = llr.range( (16 * Q) - 1, 0 );
	sc_bigint< 16 *  Q      > lb  = llr.range( (32 * Q) - 1, (16 * Q) );
	sc_bigint< 16 * (Q + 1) > sum = VECTOR_ADD_NOSAT < 16, Q > ( la, lb );

	sc_biguint< 16 > sig = REP_REP2_16_CA2 < Q + 1 > ( sum, sel );

	sc_biguint< 32 > res = (sig, sig);
	return res;
}

template <int Q>
sc_biguint< 64 > REP_REP2_64_CA2 (sc_bigint< 64 * Q> llr, bool sel)
{
#pragma HLS INLINE
	sc_bigint< 32 *  Q      > la  = llr.range( (32 * Q) - 1, 0 );
	sc_bigint< 32 *  Q      > lb  = llr.range( (64 * Q) - 1, (32 * Q) );
	sc_bigint< 32 * (Q + 1) > sum = VECTOR_ADD_NOSAT < 32, Q > ( la, lb );

	sc_biguint< 32 > sig = REP_REP2_32_CA2 < Q + 1 > ( sum, sel );

	sc_biguint< 64 > res = (sig, sig);
	return res;
}

template <int Q>
sc_biguint< 128 > REP_REP2_128_CA2 (sc_bigint< 128 * Q> llr, bool sel)
{
#pragma HLS INLINE
	sc_bigint< 64 *  Q      > la  = llr.range( (64 * Q) - 1, 0 );
	sc_bigint< 64 *  Q      > lb  = llr.range( (128 * Q) - 1, (64 * Q) );
	sc_bigint< 64 * (Q + 1) > sum = VECTOR_ADD_NOSAT < 64, Q > ( la, lb );

	sc_biguint< 64 > sig = REP_REP2_64_CA2 < Q + 1 > ( sum, sel );

	sc_biguint< 128 > res = (sig, sig);
	return res;
}

template <int Q>
sc_biguint< 256 > REP_REP2_256_CA2 (sc_bigint< 256 * Q> llr, bool sel)
{
#pragma HLS INLINE
	sc_bigint< 128 *  Q      > la  = llr.range( (128 * Q) - 1, 0 );
	sc_bigint< 128 *  Q      > lb  = llr.range( (256 * Q) - 1, (128 * Q) );
	sc_bigint< 128 * (Q + 1) > sum = VECTOR_ADD_NOSAT < 128, Q > ( la, lb );

	sc_biguint< 128 > sig = REP_REP2_128_CA2 < Q + 1 > ( sum, sel );

	sc_biguint< 256 > res = (sig, sig);
	return res;
}

/////////////// SIGMAG ////////////////

template <int Q>
sc_biguint< 2 > REP_REP2_2_SM (sc_biguint< 2 * Q> llr, bool sel)  // Sel 0 : REP, Sel 1 : REP2
{
#pragma HLS INLINE
	sc_biguint< Q 	 > la  = llr.range( Q - 1, 0 );
	sc_biguint< Q 	 > lb  = llr.range( 2 * Q - 1, Q );

	sc_biguint<1> siga = (sc_biguint<1>) la[Q - 1];
	sc_biguint<1> sigb = (sc_biguint<1>) lb[Q - 1];

	sc_biguint<Q-1> absla = la.range(Q-2 , 0);
	sc_biguint<Q-1> abslb = lb.range(Q-2 , 0);

	sc_biguint<1> is_min = ( absla < abslb ) ? 1 : 0; // if absla < abslb return 1

	sc_biguint<1> sig = (sc_uint<1>)is_min ? sigb : siga;

	sc_biguint< 2 > res1 = (sig, sig);
	sc_biguint< 2 > res2 = (sigb, siga);

	if(sel == 0)
		return res1;
	else
		return res2;
}

template <int Q>
sc_biguint< 4 > REP_REP2_4_SM (sc_biguint< 4 * Q> llr, bool sel)
{
#pragma HLS INLINE
	sc_biguint< 2 * Q 	   > la  = llr.range( (2 * Q) - 1, 0 );
	sc_biguint< 2 * Q 	   > lb  = llr.range( (4 * Q) - 1, (2 * Q) );

	sc_biguint< 2 * (Q + 1) > sum = VECTOR_FULL_ADDER_SM < 2, Q > ( la, lb );

	sc_biguint< 2 > sig = REP_REP2_2_SM < Q + 1 > ( sum, sel );

	sc_biguint< 4 > res = (sig, sig);

	return res;
}

template <int Q>
sc_biguint< 8 > REP_REP2_8_SM (sc_biguint< 8 * Q> llr, bool sel)
{
#pragma HLS INLINE
	sc_biguint< 4 * Q 	   > la  = llr.range( (4 * Q) - 1, 0 );
	sc_biguint< 4 * Q 	   > lb  = llr.range( (8 * Q) - 1, (4 * Q) );

	sc_biguint< 4 * (Q + 1) > sum = VECTOR_FULL_ADDER_SM < 4, Q > ( la, lb );

	sc_biguint< 4 > sig = REP_REP2_4_SM < Q + 1 > ( sum, sel );

	sc_biguint< 8 > res = (sig, sig);

	return res;
}

template <int Q>
sc_biguint< 16 > REP_REP2_16_SM (sc_biguint< 16 * Q> llr, bool sel)
{
#pragma HLS INLINE
	sc_biguint< 8 * Q 	   > la  = llr.range( (8 * Q) - 1, 0 );
	sc_biguint< 8 * Q 	   > lb  = llr.range( (16 * Q) - 1, (8 * Q) );

	sc_biguint< 8 * (Q + 1) > sum = VECTOR_FULL_ADDER_SM < 8, Q > ( la, lb );

	sc_biguint< 8 > sig = REP_REP2_8_SM < Q + 1 > ( sum, sel );

	sc_biguint< 16 > res = (sig, sig);

	return res;
}

template <int Q>
sc_biguint< 32 > REP_REP2_32_SM (sc_biguint< 32 * Q> llr, bool sel)
{
#pragma HLS INLINE
	sc_biguint< 16 * Q 	   > la  = llr.range( (16 * Q) - 1, 0 );
	sc_biguint< 16 * Q 	   > lb  = llr.range( (32 * Q) - 1, (16 * Q) );

	sc_biguint< 16 * (Q + 1) > sum = VECTOR_FULL_ADDER_SM < 16, Q > ( la, lb );

	sc_biguint< 16 > sig = REP_REP2_16_SM < Q + 1 > ( sum, sel );

	sc_biguint< 32 > res = (sig, sig);

	return res;
}

template <int Q>
sc_biguint< 64 > REP_REP2_64_SM (sc_biguint< 64 * Q> llr, bool sel)
{
#pragma HLS INLINE
	sc_biguint< 32 * Q 	   > la  = llr.range( (32 * Q) - 1, 0 );
	sc_biguint< 32 * Q 	   > lb  = llr.range( (64 * Q) - 1, (32 * Q) );

	sc_biguint< 32 * (Q + 1) > sum = VECTOR_FULL_ADDER_SM < 32, Q > ( la, lb );

	sc_biguint< 32 > sig = REP_REP2_32_SM < Q + 1 > ( sum, sel );

	sc_biguint< 64 > res = (sig, sig);

	return res;
}

template <int Q>
sc_biguint< 128 > REP_REP2_128_SM (sc_biguint< 128 * Q> llr, bool sel)
{
#pragma HLS INLINE
	sc_biguint< 64 * Q 	   > la  = llr.range( (64 * Q) - 1, 0 );
	sc_biguint< 64 * Q 	   > lb  = llr.range( (128 * Q) - 1, (64 * Q) );

	sc_biguint< 64 * (Q + 1) > sum = VECTOR_FULL_ADDER_SM < 64, Q > ( la, lb );

	sc_biguint< 64 > sig = REP_REP2_64_SM < Q + 1 > ( sum, sel );

	sc_biguint< 128 > res = (sig, sig);

	return res;
}

template <int Q>
sc_biguint< 256 > REP_REP2_256_SM (sc_biguint< 256 * Q> llr, bool sel)
{
#pragma HLS INLINE
	sc_biguint< 128 * Q 	   > la  = llr.range( (128 * Q) - 1, 0 );
	sc_biguint< 128 * Q 	   > lb  = llr.range( (256 * Q) - 1, (128 * Q) );

	sc_biguint< 128 * (Q + 1) > sum = VECTOR_FULL_ADDER_SM < 128, Q > ( la, lb );

	sc_biguint< 128 > sig = REP_REP2_128_SM < Q + 1 > ( sum, sel );

	sc_biguint< 256 > res = (sig, sig);

	return res;
}


/////////////// GLOBAL REP/REP2 Node ////////////////

template <int Q>
sc_biguint< 2 > REP_REP2_2_Node (sc_bigint< 2 * Q> llr, bool sel)
{
#pragma HLS INLINE

#if defined CA2
	return REP_REP2_2_CA2< Q >(llr, sel);
#elif defined SIGMAG
	return REP_REP2_2_SM< Q >( (sc_biguint< 2 * Q>) llr, sel);
#endif
}

template <int Q>
sc_biguint< 4 > REP_REP2_4_Node (sc_bigint< 4 * Q> llr, bool sel)
{
#pragma HLS INLINE

#if defined CA2
	return REP_REP2_4_CA2< Q >(llr, sel);
#elif defined SIGMAG
	return REP_REP2_4_SM< Q >( (sc_biguint< 4 * Q>) llr, sel);
#endif
}

template <int Q>
sc_biguint< 8 > REP_REP2_8_Node (sc_bigint< 8 * Q> llr, bool sel)
{
#pragma HLS INLINE

#if defined CA2
	return REP_REP2_8_CA2< Q >(llr, sel);
#elif defined SIGMAG
	return REP_REP2_8_SM< Q >( (sc_biguint< 8 * Q>) llr, sel);
#endif
}

template <int Q>
sc_biguint< 16 > REP_REP2_16_Node (sc_bigint< 16 * Q> llr, bool sel)
{
#pragma HLS INLINE

#if defined CA2
	return REP_REP2_16_CA2< Q >(llr, sel);
#elif defined SIGMAG
	return REP_REP2_16_SM< Q >( (sc_biguint< 16 * Q>) llr, sel);
#endif
}

template <int Q>
sc_biguint< 32 > REP_REP2_32_Node (sc_bigint< 32 * Q> llr, bool sel)
{
#pragma HLS INLINE

#if defined CA2
	return REP_REP2_32_CA2< Q >(llr, sel);
#elif defined SIGMAG
	return REP_REP2_32_SM< Q >( (sc_biguint< 32 * Q>) llr, sel);
#endif
}

template <int Q>
sc_biguint< 64 > REP_REP2_64_Node (sc_bigint< 64 * Q> llr, bool sel)
{
#pragma HLS INLINE

#if defined CA2
	return REP_REP2_64_CA2< Q >(llr, sel);
#elif defined SIGMAG
	return REP_REP2_64_SM< Q >( (sc_biguint< 64 * Q>) llr, sel);
#endif
}

template <int Q>
sc_biguint< 128 > REP_REP2_128_Node (sc_bigint< 128 * Q> llr, bool sel)
{
#pragma HLS INLINE

#if defined CA2
	return REP_REP2_128_CA2< Q >(llr, sel);
#elif defined SIGMAG
	return REP_REP2_128_SM< Q >( (sc_biguint< 128 * Q>) llr, sel);
#endif
}

template <int Q>
sc_biguint< 256 > REP_REP2_256_Node (sc_bigint< 256 * Q> llr, bool sel)
{
#pragma HLS INLINE

#if defined CA2
	return REP_REP2_256_CA2< Q >(llr, sel);
#elif defined SIGMAG
	return REP_REP2_256_SM< Q >( (sc_biguint< 256 * Q>) llr, sel);
#endif
}

//*************************************************************************//
//** 						SPC NODE DECODERS							 **//
//*************************************************************************//

inline sc_biguint< 2 > SPC_Parity_2 (sc_biguint< 2 > sign)
{
#pragma HLS INLINE
	sc_biguint< 1 > sa = (sc_biguint< 1 >) sign[0];
	sc_biguint< 1 > sb = (sc_biguint< 1 >) sign[1];
	sc_biguint< 1 > parity = VECTOR_XOR< 1 >(sa, sb);

	//sc_biguint< 1 > n_parity = ~ parity;
	//sc_biguint< 2 > res = (n_parity, n_parity);
	sc_biguint< 2 > res = (parity, parity);

	return res;
}

inline sc_biguint< 4 > SPC_Parity_4 (sc_biguint< 4 > sign)
{
#pragma HLS INLINE
	sc_biguint< 2 > sa = sign.range( 1, 0 );
	sc_biguint< 2 > sb = sign.range( 3, 2 );
	sc_biguint< 2 > s_xor = VECTOR_XOR< 2 >(sa, sb);

	sc_biguint< 2 > parity = SPC_Parity_2 (s_xor);

	sc_biguint< 4 > res = (parity, parity);
	return res;
}

inline sc_biguint< 8 > SPC_Parity_8 (sc_biguint< 8 > sign)
{
#pragma HLS INLINE
	sc_biguint< 4 > sa = sign.range( 3, 0 );
	sc_biguint< 4 > sb = sign.range( 7, 4 );
	sc_biguint< 4 > s_xor = VECTOR_XOR< 4 >(sa, sb);

	sc_biguint< 4 > parity = SPC_Parity_4 (s_xor);

	sc_biguint< 8 > res = (parity, parity);
	return res;
}

inline sc_biguint< 16 > SPC_Parity_16 (sc_biguint< 16 > sign)
{
#pragma HLS INLINE
	sc_biguint< 8 > sa = sign.range( 7, 0 );
	sc_biguint< 8 > sb = sign.range( 15, 8 );
	sc_biguint< 8 > s_xor = VECTOR_XOR< 8 >(sa, sb);

	sc_biguint< 8 > parity = SPC_Parity_8 (s_xor);

	sc_biguint< 16 > res = (parity, parity);
	return res;
}

inline sc_biguint< 32 > SPC_Parity_32 (sc_biguint< 32 > sign)
{
#pragma HLS INLINE
	sc_biguint< 16 > sa = sign.range( 15, 0 );
	sc_biguint< 16 > sb = sign.range( 31, 16 );
	sc_biguint< 16 > s_xor = VECTOR_XOR< 16 >(sa, sb);

	sc_biguint< 16 > parity = SPC_Parity_16 (s_xor);

	sc_biguint< 32 > res = (parity, parity);
	return res;
}

inline sc_biguint< 64 > SPC_Parity_64 (sc_biguint< 64 > sign)
{
#pragma HLS INLINE
	sc_biguint< 32 > sa = sign.range( 31, 0 );
	sc_biguint< 32 > sb = sign.range( 63, 32 );
	sc_biguint< 32 > s_xor = VECTOR_XOR< 32 >(sa, sb);

	sc_biguint< 32 > parity = SPC_Parity_32 (s_xor);

	sc_biguint< 64 > res = (parity, parity);
	return res;
}

inline sc_biguint< 128 > SPC_Parity_128 (sc_biguint< 128 > sign)
{
#pragma HLS INLINE
	sc_biguint< 64 > sa = sign.range( 63, 0 );
	sc_biguint< 64 > sb = sign.range( 127, 64 );
	sc_biguint< 64 > s_xor = VECTOR_XOR< 64 >(sa, sb);

	sc_biguint< 64 > parity = SPC_Parity_64 (s_xor);

	sc_biguint< 128 > res = (parity, parity);
	return res;
}

inline sc_biguint< 256 > SPC_Parity_256 (sc_biguint< 256 > sign)
{
#pragma HLS INLINE
	sc_biguint< 128 > sa = sign.range( 127, 0 );
	sc_biguint< 128 > sb = sign.range( 255, 128 );
	sc_biguint< 128 > s_xor = VECTOR_XOR< 128 >(sa, sb);

	sc_biguint< 128 > parity = SPC_Parity_128 (s_xor);

	sc_biguint< 256 > res = (parity, parity);
	return res;
}

/////////////// CA2 ////////////////

template <int Q>
sc_biguint< 2 > SPC_Min_Mask_2_CA2 (sc_bigint< 2 * Q> Mag)
{
#pragma HLS INLINE
	sc_bigint< Q > ma  = Mag.range( Q - 1, 0 );
	sc_bigint< Q > mb  = Mag.range( 2 * Q - 1, Q );

	sc_biguint< 1 > is_min = VECTOR_IS_MIN< 1, Q >(mb, ma);
	sc_biguint< 1 > n_is_min = (~is_min);

	sc_biguint< 2 > mask = (is_min, n_is_min);
	return mask;
}

template <int Q>
sc_biguint< 4 > SPC_Min_Mask_4_CA2 (sc_bigint< 4 * Q> Mag)
{
#pragma HLS INLINE
	sc_bigint< 2 * Q > ma  = Mag.range( (2 * Q) - 1, 0 );
	sc_bigint< 2 * Q > mb  = Mag.range( (4 * Q) - 1, 2 * Q );

	sc_biguint< 2 > is_min = VECTOR_IS_MIN< 2, Q >( mb, ma );
	sc_biguint< 2 > n_is_min = (~is_min);

	sc_biguint< 4 > mask_A = (is_min, n_is_min );

	sc_bigint< 2 * Q > min = VECTOR_MUX< 2 , Q >( mb, ma, is_min );

	sc_biguint< 2 > i_mask = SPC_Min_Mask_2_CA2< Q >( min );
	sc_biguint< 4 > mask_B = (i_mask, i_mask);

	sc_biguint< 4 > mask = VECTOR_AND< 4 >( mask_A, mask_B );
	return mask;
}

template <int Q>
sc_biguint< 8 > SPC_Min_Mask_8_CA2 (sc_bigint< 8 * Q> Mag)
{
#pragma HLS INLINE
	sc_bigint< 4 * Q > ma  = Mag.range( (4 * Q) - 1, 0 );
	sc_bigint< 4 * Q > mb  = Mag.range( (8 * Q) - 1, 4 * Q );

	sc_biguint< 4 > is_min = VECTOR_IS_MIN< 4, Q >( mb, ma );
	sc_biguint< 4 > n_is_min = (~is_min);

	sc_biguint< 8 > mask_A = (is_min, n_is_min );

	sc_bigint< 4 * Q > min = VECTOR_MUX< 4 , Q >( mb, ma, is_min );

	sc_biguint< 4 > i_mask = SPC_Min_Mask_4_CA2< Q >( min );
	sc_biguint< 8 > mask_B = (i_mask, i_mask);

	sc_biguint< 8 > mask = VECTOR_AND< 8 >( mask_A, mask_B );
	return mask;
}

template <int Q>
sc_biguint< 16 > SPC_Min_Mask_16_CA2 (sc_bigint< 16 * Q> Mag)
{
#pragma HLS INLINE
	sc_bigint< 8 * Q > ma  = Mag.range( (8 * Q) - 1, 0 );
	sc_bigint< 8 * Q > mb  = Mag.range( (16 * Q) - 1, 8 * Q );

	sc_biguint< 8 > is_min = VECTOR_IS_MIN< 8, Q >( mb, ma );
	sc_biguint< 8 > n_is_min = (~is_min);

	sc_biguint< 16 > mask_A = (is_min, n_is_min );

	sc_bigint< 8 * Q > min = VECTOR_MUX< 8 , Q >( mb, ma, is_min );

	sc_biguint< 8 > i_mask = SPC_Min_Mask_8_CA2< Q >( min );
	sc_biguint< 16 > mask_B = (i_mask, i_mask);

	sc_biguint< 16 > mask = VECTOR_AND< 16 >( mask_A, mask_B );
	return mask;
}

template <int Q>
sc_biguint< 32 > SPC_Min_Mask_32_CA2 (sc_bigint< 32 * Q> Mag)
{
#pragma HLS INLINE
	sc_bigint< 16 * Q > ma  = Mag.range( (16 * Q) - 1, 0 );
	sc_bigint< 16 * Q > mb  = Mag.range( (32 * Q) - 1, 16 * Q );

	sc_biguint< 16 > is_min = VECTOR_IS_MIN< 16, Q >( mb, ma );
	sc_biguint< 16 > n_is_min = (~is_min);

	sc_biguint< 32 > mask_A = (is_min, n_is_min );

	sc_bigint< 16 * Q > min = VECTOR_MUX< 16 , Q >( mb, ma, is_min );

	sc_biguint< 16 > i_mask = SPC_Min_Mask_16_CA2< Q >( min );
	sc_biguint< 32 > mask_B = (i_mask, i_mask);

	sc_biguint< 32 > mask = VECTOR_AND< 32 >( mask_A, mask_B );
	return mask;
}

template <int Q>
sc_biguint< 64 > SPC_Min_Mask_64_CA2 (sc_bigint< 64 * Q> Mag)
{
#pragma HLS INLINE
	sc_bigint< 32 * Q > ma  = Mag.range( (32 * Q) - 1, 0 );
	sc_bigint< 32 * Q > mb  = Mag.range( (64 * Q) - 1, 32 * Q );

	sc_biguint< 32 > is_min = VECTOR_IS_MIN< 32, Q >( mb, ma );
	sc_biguint< 32 > n_is_min = (~is_min);

	sc_biguint< 64 > mask_A = (is_min, n_is_min );

	sc_bigint< 32 * Q > min = VECTOR_MUX< 32 , Q >( mb, ma, is_min );

	sc_biguint< 32 > i_mask = SPC_Min_Mask_32_CA2< Q >( min );
	sc_biguint< 64 > mask_B = (i_mask, i_mask);

	sc_biguint< 64 > mask = VECTOR_AND< 64 >( mask_A, mask_B );
	return mask;
}

template <int Q>
sc_biguint< 128 > SPC_Min_Mask_128_CA2 (sc_bigint< 128 * Q> Mag)
{
#pragma HLS INLINE
	sc_bigint< 64 * Q > ma  = Mag.range( (64 * Q) - 1, 0 );
	sc_bigint< 64 * Q > mb  = Mag.range( (128 * Q) - 1, 64 * Q );

	sc_biguint< 64 > is_min = VECTOR_IS_MIN< 64, Q >( mb, ma );
	sc_biguint< 64 > n_is_min = (~is_min);

	sc_biguint< 128 > mask_A = (is_min, n_is_min );

	sc_bigint< 64 * Q > min = VECTOR_MUX< 64 , Q >( mb, ma, is_min );

	sc_biguint< 64 > i_mask = SPC_Min_Mask_64_CA2< Q >( min );
	sc_biguint< 128 > mask_B = (i_mask, i_mask);

	sc_biguint< 128 > mask = VECTOR_AND< 128 >( mask_A, mask_B );
	return mask;
}

template <int Q>
sc_biguint< 256 > SPC_Min_Mask_256_CA2 (sc_bigint< 256 * Q> Mag)
{
#pragma HLS INLINE
	sc_bigint< 128 * Q > ma  = Mag.range( (128 * Q) - 1, 0 );
	sc_bigint< 128 * Q > mb  = Mag.range( (256 * Q) - 1, 128 * Q );

	sc_biguint< 128 > is_min = VECTOR_IS_MIN< 128, Q >( mb, ma );
	sc_biguint< 128 > n_is_min = (~is_min);

	sc_biguint< 256 > mask_A = (is_min, n_is_min );

	sc_bigint< 128 * Q > min = VECTOR_MUX< 128 , Q >( mb, ma, is_min );

	sc_biguint< 128 > i_mask = SPC_Min_Mask_128_CA2< Q >( min );
	sc_biguint< 256 > mask_B = (i_mask, i_mask);

	sc_biguint< 256 > mask = VECTOR_AND< 256 >( mask_A, mask_B );
	return mask;
}
/////////////// SIGMAG ////////////////

template <int Q>
sc_biguint< 2 > SPC_Min_Mask_2_SM (sc_biguint< 2 * Q> Mag)
{
#pragma HLS INLINE
	sc_biguint< Q > ma  = Mag.range( Q - 1, 0 );
	sc_biguint< Q > mb  = Mag.range( 2 * Q - 1, Q );

	sc_biguint< 1 > is_min = VECTOR_IS_MIN_SM< 1, Q >(mb, ma);
	sc_biguint< 1 > n_is_min = (~is_min);

	sc_biguint< 2 > mask = (is_min, n_is_min);
	return mask;
}

template <int Q>
sc_biguint< 4 > SPC_Min_Mask_4_SM (sc_biguint< 4 * Q> Mag)
{
#pragma HLS INLINE
	sc_biguint< 2 * Q > ma  = Mag.range( (2 * Q) - 1, 0 );
	sc_biguint< 2 * Q > mb  = Mag.range( (4 * Q) - 1, 2 * Q );

	sc_biguint< 2 > is_min = VECTOR_IS_MIN_SM< 2, Q >( mb, ma );
	sc_biguint< 2 > n_is_min = (~is_min);

	sc_biguint< 4 > mask_A = (is_min, n_is_min );

	sc_biguint< 2 * Q > min = VECTOR_MUX_SM< 2 , Q >( mb, ma, is_min );

	sc_biguint< 2 > i_mask = SPC_Min_Mask_2_SM< Q >( min );
	sc_biguint< 4 > mask_B = (i_mask, i_mask);

	sc_biguint< 4 > mask = VECTOR_AND< 4 >( mask_A, mask_B );
	return mask;
}

template <int Q>
sc_biguint< 8 > SPC_Min_Mask_8_SM (sc_biguint< 8 * Q> Mag)
{
#pragma HLS INLINE
	sc_biguint< 4 * Q > ma  = Mag.range( (4 * Q) - 1, 0 );
	sc_biguint< 4 * Q > mb  = Mag.range( (8 * Q) - 1, 4 * Q );

	sc_biguint< 4 > is_min = VECTOR_IS_MIN_SM< 4, Q >( mb, ma );
	sc_biguint< 4 > n_is_min = (~is_min);

	sc_biguint< 8 > mask_A = (is_min, n_is_min );

	sc_biguint< 4 * Q > min = VECTOR_MUX_SM< 4 , Q >( mb, ma, is_min );

	sc_biguint< 4 > i_mask = SPC_Min_Mask_4_SM< Q >( min );
	sc_biguint< 8 > mask_B = (i_mask, i_mask);

	sc_biguint< 8 > mask = VECTOR_AND< 8 >( mask_A, mask_B );
	return mask;
}

template <int Q>
sc_biguint< 16 > SPC_Min_Mask_16_SM (sc_biguint< 16 * Q> Mag)
{
#pragma HLS INLINE
	sc_biguint< 8 * Q > ma  = Mag.range( (8 * Q) - 1, 0 );
	sc_biguint< 8 * Q > mb  = Mag.range( (16 * Q) - 1, 8 * Q );

	sc_biguint< 8 > is_min = VECTOR_IS_MIN_SM< 8, Q >( mb, ma );
	sc_biguint< 8 > n_is_min = (~is_min);

	sc_biguint< 16 > mask_A = (is_min, n_is_min );

	sc_biguint< 8 * Q > min = VECTOR_MUX_SM< 8 , Q >( mb, ma, is_min );

	sc_biguint< 8 > i_mask = SPC_Min_Mask_8_SM< Q >( min );
	sc_biguint< 16 > mask_B = (i_mask, i_mask);

	sc_biguint< 16 > mask = VECTOR_AND< 16 >( mask_A, mask_B );
	return mask;
}

template <int Q>
sc_biguint< 32 > SPC_Min_Mask_32_SM (sc_biguint< 32 * Q> Mag)
{
#pragma HLS INLINE
	sc_biguint< 16 * Q > ma  = Mag.range( (16 * Q) - 1, 0 );
	sc_biguint< 16 * Q > mb  = Mag.range( (32 * Q) - 1, 16 * Q );

	sc_biguint< 16 > is_min = VECTOR_IS_MIN_SM< 16, Q >( mb, ma );
	sc_biguint< 16 > n_is_min = (~is_min);

	sc_biguint< 32 > mask_A = (is_min, n_is_min );

	sc_biguint< 16 * Q > min = VECTOR_MUX_SM< 16 , Q >( mb, ma, is_min );

	sc_biguint< 16 > i_mask = SPC_Min_Mask_16_SM< Q >( min );
	sc_biguint< 32 > mask_B = (i_mask, i_mask);

	sc_biguint< 32 > mask = VECTOR_AND< 32 >( mask_A, mask_B );
	return mask;
}

template <int Q>
sc_biguint< 64 > SPC_Min_Mask_64_SM (sc_biguint< 64 * Q> Mag)
{
#pragma HLS INLINE
	sc_biguint< 32 * Q > ma  = Mag.range( (32 * Q) - 1, 0 );
	sc_biguint< 32 * Q > mb  = Mag.range( (64 * Q) - 1, 32 * Q );

	sc_biguint< 32 > is_min = VECTOR_IS_MIN_SM< 32, Q >( mb, ma );
	sc_biguint< 32 > n_is_min = (~is_min);

	sc_biguint< 64 > mask_A = (is_min, n_is_min );

	sc_biguint< 32 * Q > min = VECTOR_MUX_SM< 32 , Q >( mb, ma, is_min );

	sc_biguint< 32 > i_mask = SPC_Min_Mask_32_SM< Q >( min );
	sc_biguint< 64 > mask_B = (i_mask, i_mask);

	sc_biguint< 64 > mask = VECTOR_AND< 64 >( mask_A, mask_B );
	return mask;
}

template <int Q>
sc_biguint< 128 > SPC_Min_Mask_128_SM (sc_biguint< 128 * Q> Mag)
{
#pragma HLS INLINE
	sc_biguint< 64 * Q > ma  = Mag.range( (64 * Q) - 1, 0 );
	sc_biguint< 64 * Q > mb  = Mag.range( (128 * Q) - 1, 64 * Q );

	sc_biguint< 64 > is_min = VECTOR_IS_MIN_SM< 64, Q >( mb, ma );
	sc_biguint< 64 > n_is_min = (~is_min);

	sc_biguint< 128 > mask_A = (is_min, n_is_min );

	sc_biguint< 64 * Q > min = VECTOR_MUX_SM< 64 , Q >( mb, ma, is_min );

	sc_biguint< 64 > i_mask = SPC_Min_Mask_64_SM< Q >( min );
	sc_biguint< 128 > mask_B = (i_mask, i_mask);

	sc_biguint< 128 > mask = VECTOR_AND< 128 >( mask_A, mask_B );
	return mask;
}

template <int Q>
sc_biguint< 256 > SPC_Min_Mask_256_SM (sc_biguint< 256 * Q> Mag)
{
#pragma HLS INLINE
	sc_biguint< 128 * Q > ma  = Mag.range( (128 * Q) - 1, 0 );
	sc_biguint< 128 * Q > mb  = Mag.range( (256 * Q) - 1, 128 * Q );

	sc_biguint< 128 > is_min = VECTOR_IS_MIN_SM< 128, Q >( mb, ma );
	sc_biguint< 128 > n_is_min = (~is_min);

	sc_biguint< 256 > mask_A = (is_min, n_is_min );

	sc_biguint< 128 * Q > min = VECTOR_MUX_SM< 128 , Q >( mb, ma, is_min );

	sc_biguint< 128 > i_mask = SPC_Min_Mask_128_SM< Q >( min );
	sc_biguint< 256 > mask_B = (i_mask, i_mask);

	sc_biguint< 256 > mask = VECTOR_AND< 256 >( mask_A, mask_B );
	return mask;
}

/////////////// GLOBAL SPC Node ////////////////

template <int Q>
sc_biguint< 2 > SPC_Node_2 (sc_bigint< 2 * Q> llr)
{
#pragma HLS INLINE

#if defined CA2
	sc_biguint< 2 > sign = VECTOR_SIGN< 2, Q >(llr);
	sc_bigint < 2 * Q > mag = VECTOR_ABS< 2, Q >(llr);

	sc_biguint< 2 > parity = SPC_Parity_2( sign );
	sc_biguint< 2 > mask = SPC_Min_Mask_2_CA2< Q >( mag );

	sc_biguint< 2 > and_mask = VECTOR_AND< 2 >( parity, mask) ;
	sc_biguint< 2 > bit = VECTOR_XOR< 2 >( sign, and_mask );
	return bit;
#elif defined SIGMAG
	sc_biguint< 2 > sign = VECTOR_SIGN_SM< 2, Q >( (sc_biguint< 2 * Q>) llr);
	sc_biguint< 2 * (Q-1) > mag = VECTOR_ABS_SM< 2, Q >( (sc_biguint< 2 * Q>) llr);

	sc_biguint< 2 > parity = SPC_Parity_2( sign );
	sc_biguint< 2 > mask = SPC_Min_Mask_2_SM< Q-1 >( mag );

	sc_biguint< 2 > and_mask = VECTOR_AND< 2 >( parity, mask) ;
	sc_biguint< 2 > bit = VECTOR_XOR< 2 >( sign, and_mask );
	return bit;
#endif
}

template <int Q>
sc_biguint< 4 > SPC_Node_4 (sc_bigint< 4 * Q> llr)
{
#pragma HLS INLINE

#if defined CA2
	sc_biguint< 4 > sign = VECTOR_SIGN< 4, Q >(llr);
	sc_bigint < 4 * Q > mag = VECTOR_ABS< 4, Q >(llr);

	sc_biguint< 4 > parity = SPC_Parity_4( sign );
	sc_biguint< 4 > mask = SPC_Min_Mask_4_CA2< Q >( mag );

	sc_biguint< 4 > and_mask = VECTOR_AND< 4 >( parity, mask) ;
	sc_biguint< 4 > bit = VECTOR_XOR< 4 >( sign, and_mask );

	return bit;
#elif defined SIGMAG
	sc_biguint< 4 > sign = VECTOR_SIGN_SM< 4, Q >( (sc_biguint< 4 * Q>) llr);
	sc_biguint< 4 * (Q-1) > mag = VECTOR_ABS_SM< 4, Q >( (sc_biguint< 4 * Q>) llr);

	sc_biguint< 4 > parity = SPC_Parity_4( sign );
	sc_biguint< 4 > mask = SPC_Min_Mask_4_SM< Q-1 >( mag );

	sc_biguint< 4 > and_mask = VECTOR_AND< 4 >( parity, mask) ;
	sc_biguint< 4 > bit = VECTOR_XOR< 4 >( sign, and_mask );
	return bit;
#endif
}

template <int Q>
sc_biguint< 8 > SPC_Node_8 (sc_bigint< 8 * Q> llr)
{
#pragma HLS INLINE

#if defined CA2
	sc_biguint< 8 > sign = VECTOR_SIGN< 8, Q >(llr);
	sc_bigint < 8 * Q > mag = VECTOR_ABS< 8, Q >(llr);

	sc_biguint< 8 > parity = SPC_Parity_8( sign );
	sc_biguint< 8 > mask = SPC_Min_Mask_8_CA2< Q >( mag );

	sc_biguint< 8 > and_mask = VECTOR_AND< 8 >( parity, mask) ;
	sc_biguint< 8 > bit = VECTOR_XOR< 8 >( sign, and_mask );
	return bit;
#elif defined SIGMAG
	sc_biguint< 8 > sign = VECTOR_SIGN_SM< 8, Q >( (sc_biguint< 8 * Q>) llr);
	sc_biguint< 8 * (Q-1) > mag = VECTOR_ABS_SM< 8, Q >( (sc_biguint< 8 * Q>) llr);

	sc_biguint< 8 > parity = SPC_Parity_8( sign );
	sc_biguint< 8 > mask = SPC_Min_Mask_8_SM< Q-1 >( mag );

	sc_biguint< 8 > and_mask = VECTOR_AND< 8 >( parity, mask) ;
	sc_biguint< 8 > bit = VECTOR_XOR< 8 >( sign, and_mask );
	return bit;
#endif
}

template <int Q>
sc_biguint< 16 > SPC_Node_16 (sc_bigint< 16 * Q> llr)
{
#pragma HLS INLINE

#if defined CA2
	sc_biguint< 16 > sign = VECTOR_SIGN< 16, Q >(llr);
	sc_bigint < 16 * Q > mag = VECTOR_ABS< 16, Q >(llr);

	sc_biguint< 16 > parity = SPC_Parity_16( sign );
	sc_biguint< 16 > mask = SPC_Min_Mask_16_CA2< Q >( mag );

	sc_biguint< 16 > and_mask = VECTOR_AND< 16 >( parity, mask) ;
	sc_biguint< 16 > bit = VECTOR_XOR< 16 >( sign, and_mask );
	return bit;
#elif defined SIGMAG
	sc_biguint< 16 > sign = VECTOR_SIGN_SM< 16, Q >( (sc_biguint< 16 * Q>) llr);
	sc_biguint< 16 * (Q-1) > mag = VECTOR_ABS_SM< 16, Q >( (sc_biguint< 16 * Q>) llr);

	sc_biguint< 16 > parity = SPC_Parity_16( sign );
	sc_biguint< 16 > mask = SPC_Min_Mask_16_SM< Q-1 >( mag );

	sc_biguint< 16 > and_mask = VECTOR_AND< 16 >( parity, mask) ;
	sc_biguint< 16 > bit = VECTOR_XOR< 16 >( sign, and_mask );
	return bit;
#endif
}

template <int Q>
sc_biguint< 32 > SPC_Node_32 (sc_bigint< 32 * Q> llr)
{
#pragma HLS INLINE

#if defined CA2
	sc_biguint< 32 > sign = VECTOR_SIGN< 32, Q >(llr);
	sc_bigint < 32 * Q > mag = VECTOR_ABS< 32, Q >(llr);

	sc_biguint< 32 > parity = SPC_Parity_32( sign );
	sc_biguint< 32 > mask = SPC_Min_Mask_32_CA2< Q >( mag );

	sc_biguint< 32 > and_mask = VECTOR_AND< 32 >( parity, mask) ;
	sc_biguint< 32 > bit = VECTOR_XOR< 32 >( sign, and_mask );
	return bit;
#elif defined SIGMAG
	sc_biguint< 32 > sign = VECTOR_SIGN_SM< 32, Q >( (sc_biguint< 32 * Q>) llr);
	sc_biguint< 32 * (Q-1) > mag = VECTOR_ABS_SM< 32, Q >( (sc_biguint< 32 * Q>) llr);

	sc_biguint< 32 > parity = SPC_Parity_32( sign );
	sc_biguint< 32 > mask = SPC_Min_Mask_32_SM< Q-1 >( mag );

	sc_biguint< 32 > and_mask = VECTOR_AND< 32 >( parity, mask) ;
	sc_biguint< 32 > bit = VECTOR_XOR< 32 >( sign, and_mask );
	return bit;
#endif
}

template <int Q>
sc_biguint< 64 > SPC_Node_64 (sc_bigint< 64 * Q> llr)
{
#pragma HLS INLINE

#if defined CA2
	sc_biguint< 64 > sign = VECTOR_SIGN< 64, Q >(llr);
	sc_bigint < 64 * Q > mag = VECTOR_ABS< 64, Q >(llr);

	sc_biguint< 64 > parity = SPC_Parity_64( sign );
	sc_biguint< 64 > mask = SPC_Min_Mask_64_CA2< Q >( mag );

	sc_biguint< 64 > and_mask = VECTOR_AND< 64 >( parity, mask) ;
	sc_biguint< 64 > bit = VECTOR_XOR< 64 >( sign, and_mask );
	return bit;
#elif defined SIGMAG
	sc_biguint< 64 > sign = VECTOR_SIGN_SM< 64, Q >( (sc_biguint< 64 * Q>) llr);
	sc_biguint< 64 * (Q-1) > mag = VECTOR_ABS_SM< 64, Q >( (sc_biguint< 64 * Q>) llr);

	sc_biguint< 64 > parity = SPC_Parity_64( sign );
	sc_biguint< 64 > mask = SPC_Min_Mask_64_SM< Q-1 >( mag );

	sc_biguint< 64 > and_mask = VECTOR_AND< 64 >( parity, mask) ;
	sc_biguint< 64 > bit = VECTOR_XOR< 64 >( sign, and_mask );
	return bit;
#endif
}

template <int Q>
sc_biguint< 128 > SPC_Node_128 (sc_bigint< 128 * Q> llr)
{
#pragma HLS INLINE

#if defined CA2
	sc_biguint< 128 > sign = VECTOR_SIGN< 128, Q >(llr);
	sc_bigint < 128 * Q > mag = VECTOR_ABS< 128, Q >(llr);

	sc_biguint< 128 > parity = SPC_Parity_128( sign );
	sc_biguint< 128 > mask = SPC_Min_Mask_128_CA2< Q >( mag );

	sc_biguint< 128 > and_mask = VECTOR_AND< 128 >( parity, mask) ;
	sc_biguint< 128 > bit = VECTOR_XOR< 128 >( sign, and_mask );
	return bit;
#elif defined SIGMAG
	sc_biguint< 128 > sign = VECTOR_SIGN_SM< 128, Q >( (sc_biguint< 128 * Q>) llr);
	sc_biguint< 128 * (Q-1) > mag = VECTOR_ABS_SM< 128, Q >( (sc_biguint< 128 * Q>) llr);

	sc_biguint< 128 > parity = SPC_Parity_128( sign );
	sc_biguint< 128 > mask = SPC_Min_Mask_128_SM< Q-1 >( mag );

	sc_biguint< 128 > and_mask = VECTOR_AND< 128 >( parity, mask) ;
	sc_biguint< 128 > bit = VECTOR_XOR< 128 >( sign, and_mask );
	return bit;
#endif
}

template <int Q>
sc_biguint< 256 > SPC_Node_256 (sc_bigint< 256 * Q> llr)
{
#pragma HLS INLINE

#if defined CA2
	sc_biguint< 256 > sign = VECTOR_SIGN< 256, Q >(llr);
	sc_bigint < 256 * Q > mag = VECTOR_ABS< 256, Q >(llr);

	sc_biguint< 256 > parity = SPC_Parity_256( sign );
	sc_biguint< 256 > mask = SPC_Min_Mask_256_CA2< Q >( mag );

	sc_biguint< 256 > and_mask = VECTOR_AND< 256 >( parity, mask) ;
	sc_biguint< 256 > bit = VECTOR_XOR< 256 >( sign, and_mask );
	return bit;
#elif defined SIGMAG
	sc_biguint< 256 > sign = VECTOR_SIGN_SM< 256, Q >( (sc_biguint< 256 * Q>) llr);
	sc_biguint< 256 * (Q-1) > mag = VECTOR_ABS_SM< 256, Q >( (sc_biguint< 256 * Q>) llr);

	sc_biguint< 256 > parity = SPC_Parity_256( sign );
	sc_biguint< 256 > mask = SPC_Min_Mask_256_SM< Q-1 >( mag );

	sc_biguint< 256 > and_mask = VECTOR_AND< 256 >( parity, mask) ;
	sc_biguint< 256 > bit = VECTOR_XOR< 256 >( sign, and_mask );
	return bit;
#endif
}

//*************************************************************************//
//** 					SPC / SPC2 NODE DECODERS						 **//
//*************************************************************************//

inline sc_biguint< 2 > SPC_SPC2_Parity_2 (sc_biguint< 2 > sign, bool sel) // Sel 0 : SPC, Sel 1 : SPC2
{
#pragma HLS INLINE
	sc_biguint< 1 > sa = (sc_biguint< 1 >) sign[0];
	sc_biguint< 1 > sb = (sc_biguint< 1 >) sign[1];
	sc_biguint< 1 > parity = VECTOR_XOR< 1 >(sa, sb);

	sc_biguint< 2 > res1 = (parity, parity);

	sc_biguint< 2 > res2 = sign;

	if(sel == 0)
		return res1;
	else
		return res2;
}

inline sc_biguint< 4 > SPC_SPC2_Parity_4 (sc_biguint< 4 > sign, bool sel)
{
#pragma HLS INLINE
	sc_biguint< 2 > sa = sign.range( 1, 0 );
	sc_biguint< 2 > sb = sign.range( 3, 2 );
	sc_biguint< 2 > s_xor = VECTOR_XOR< 2 >(sa, sb);

	sc_biguint< 2 > parity = SPC_SPC2_Parity_2 (s_xor, sel);

	sc_biguint< 4 > res = (parity, parity);
	return res;
}

inline sc_biguint< 8 > SPC_SPC2_Parity_8 (sc_biguint< 8 > sign, bool sel)
{
#pragma HLS INLINE
	sc_biguint< 4 > sa = sign.range( 3, 0 );
	sc_biguint< 4 > sb = sign.range( 7, 4 );
	sc_biguint< 4 > s_xor = VECTOR_XOR< 4 >(sa, sb);

	sc_biguint< 4 > parity = SPC_SPC2_Parity_4 (s_xor, sel);

	sc_biguint< 8 > res = (parity, parity);
	return res;
}

inline sc_biguint< 16 > SPC_SPC2_Parity_16 (sc_biguint< 16 > sign, bool sel)
{
#pragma HLS INLINE
	sc_biguint< 8 > sa = sign.range( 7, 0 );
	sc_biguint< 8 > sb = sign.range( 15, 8 );
	sc_biguint< 8 > s_xor = VECTOR_XOR< 8 >(sa, sb);

	sc_biguint< 8 > parity = SPC_SPC2_Parity_8 (s_xor, sel);

	sc_biguint< 16 > res = (parity, parity);
	return res;
}

inline sc_biguint< 32 > SPC_SPC2_Parity_32 (sc_biguint< 32 > sign, bool sel)
{
#pragma HLS INLINE
	sc_biguint< 16 > sa = sign.range( 15, 0 );
	sc_biguint< 16 > sb = sign.range( 31, 16 );
	sc_biguint< 16 > s_xor = VECTOR_XOR< 16 >(sa, sb);

	sc_biguint< 16 > parity = SPC_SPC2_Parity_16 (s_xor, sel);

	sc_biguint< 32 > res = (parity, parity);
	return res;
}

inline sc_biguint< 64 > SPC_SPC2_Parity_64 (sc_biguint< 64 > sign, bool sel)
{
#pragma HLS INLINE
	sc_biguint< 32 > sa = sign.range( 31, 0 );
	sc_biguint< 32 > sb = sign.range( 63, 32 );
	sc_biguint< 32 > s_xor = VECTOR_XOR< 32 >(sa, sb);

	sc_biguint< 32 > parity = SPC_SPC2_Parity_32 (s_xor, sel);

	sc_biguint< 64 > res = (parity, parity);
	return res;
}

inline sc_biguint< 128 > SPC_SPC2_Parity_128 (sc_biguint< 128 > sign, bool sel)
{
#pragma HLS INLINE
	sc_biguint< 64 > sa = sign.range( 63, 0 );
	sc_biguint< 64 > sb = sign.range( 127, 64 );
	sc_biguint< 64 > s_xor = VECTOR_XOR< 64 >(sa, sb);

	sc_biguint< 64 > parity = SPC_SPC2_Parity_64 (s_xor, sel);

	sc_biguint< 128 > res = (parity, parity);
	return res;
}

inline sc_biguint< 256 > SPC_SPC2_Parity_256 (sc_biguint< 256 > sign, bool sel)
{
#pragma HLS INLINE
	sc_biguint< 128 > sa = sign.range( 127, 0 );
	sc_biguint< 128 > sb = sign.range( 255, 128 );
	sc_biguint< 128 > s_xor = VECTOR_XOR< 128 >(sa, sb);

	sc_biguint< 128 > parity = SPC_SPC2_Parity_128 (s_xor, sel);

	sc_biguint< 256 > res = (parity, parity);
	return res;
}

/////////////// CA2 ////////////////

template <int Q>
sc_biguint< 2 > SPC_SPC2_Min_Mask_2_CA2 (sc_bigint< 2 * Q> Mag, bool sel) // Sel 0 : SPC, Sel 1 : SPC2
{
#pragma HLS INLINE
	sc_bigint< Q > ma  = Mag.range( Q - 1, 0 );
	sc_bigint< Q > mb  = Mag.range( 2 * Q - 1, Q );

	sc_biguint< 1 > is_min = VECTOR_IS_MIN< 1, Q >(mb, ma);
	sc_biguint< 1 > n_is_min = (~is_min);

	sc_biguint< 2 > mask1 = (is_min, n_is_min);

	sc_biguint< 2 > mask2 = 0x03; // (11)

	if(sel == 0)
		return mask1;
	else
		return mask2;
}

template <int Q>
sc_biguint< 4 > SPC_SPC2_Min_Mask_4_CA2 (sc_bigint< 4 * Q> Mag, bool sel)
{
#pragma HLS INLINE
	sc_bigint< 2 * Q > ma  = Mag.range( (2 * Q) - 1, 0 );
	sc_bigint< 2 * Q > mb  = Mag.range( (4 * Q) - 1, 2 * Q );

	sc_biguint< 2 > is_min = VECTOR_IS_MIN< 2, Q >( mb, ma );
	sc_biguint< 2 > n_is_min = (~is_min);

	sc_biguint< 4 > mask_A = (is_min, n_is_min );

	sc_bigint< 2 * Q > min = VECTOR_MUX< 2 , Q >( mb, ma, is_min );

	sc_biguint< 2 > i_mask = SPC_SPC2_Min_Mask_2_CA2< Q >( min, sel );
	sc_biguint< 4 > mask_B = (i_mask, i_mask);

	sc_biguint< 4 > mask = VECTOR_AND< 4 >( mask_A, mask_B );
	return mask;
}

template <int Q>
sc_biguint< 8 > SPC_SPC2_Min_Mask_8_CA2 (sc_bigint< 8 * Q> Mag, bool sel)
{
#pragma HLS INLINE
	sc_bigint< 4 * Q > ma  = Mag.range( (4 * Q) - 1, 0 );
	sc_bigint< 4 * Q > mb  = Mag.range( (8 * Q) - 1, 4 * Q );

	sc_biguint< 4 > is_min = VECTOR_IS_MIN< 4, Q >( mb, ma );
	sc_biguint< 4 > n_is_min = (~is_min);

	sc_biguint< 8 > mask_A = (is_min, n_is_min );

	sc_bigint< 4 * Q > min = VECTOR_MUX< 4 , Q >( mb, ma, is_min );

	sc_biguint< 4 > i_mask = SPC_SPC2_Min_Mask_4_CA2< Q >( min, sel );
	sc_biguint< 8 > mask_B = (i_mask, i_mask);

	sc_biguint< 8 > mask = VECTOR_AND< 8 >( mask_A, mask_B );
	return mask;
}

template <int Q>
sc_biguint< 16 > SPC_SPC2_Min_Mask_16_CA2 (sc_bigint< 16 * Q> Mag, bool sel)
{
#pragma HLS INLINE
	sc_bigint< 8 * Q > ma  = Mag.range( (8 * Q) - 1, 0 );
	sc_bigint< 8 * Q > mb  = Mag.range( (16 * Q) - 1, 8 * Q );

	sc_biguint< 8 > is_min = VECTOR_IS_MIN< 8, Q >( mb, ma );
	sc_biguint< 8 > n_is_min = (~is_min);

	sc_biguint< 16 > mask_A = (is_min, n_is_min );

	sc_bigint< 8 * Q > min = VECTOR_MUX< 8 , Q >( mb, ma, is_min );

	sc_biguint< 8 > i_mask = SPC_SPC2_Min_Mask_8_CA2< Q >( min, sel );
	sc_biguint< 16 > mask_B = (i_mask, i_mask);

	sc_biguint< 16 > mask = VECTOR_AND< 16 >( mask_A, mask_B );
	return mask;
}

template <int Q>
sc_biguint< 32 > SPC_SPC2_Min_Mask_32_CA2 (sc_bigint< 32 * Q> Mag, bool sel)
{
#pragma HLS INLINE
	sc_bigint< 16 * Q > ma  = Mag.range( (16 * Q) - 1, 0 );
	sc_bigint< 16 * Q > mb  = Mag.range( (32 * Q) - 1, 16 * Q );

	sc_biguint< 16 > is_min = VECTOR_IS_MIN< 16, Q >( mb, ma );
	sc_biguint< 16 > n_is_min = (~is_min);

	sc_biguint< 32 > mask_A = (is_min, n_is_min );

	sc_bigint< 16 * Q > min = VECTOR_MUX< 16 , Q >( mb, ma, is_min );

	sc_biguint< 16 > i_mask = SPC_SPC2_Min_Mask_16_CA2< Q >( min, sel );
	sc_biguint< 32 > mask_B = (i_mask, i_mask);

	sc_biguint< 32 > mask = VECTOR_AND< 32 >( mask_A, mask_B );
	return mask;
}

template <int Q>
sc_biguint< 64 > SPC_SPC2_Min_Mask_64_CA2 (sc_bigint< 64 * Q> Mag, bool sel)
{
#pragma HLS INLINE
	sc_bigint< 32 * Q > ma  = Mag.range( (32 * Q) - 1, 0 );
	sc_bigint< 32 * Q > mb  = Mag.range( (64 * Q) - 1, 32 * Q );

	sc_biguint< 32 > is_min = VECTOR_IS_MIN< 32, Q >( mb, ma );
	sc_biguint< 32 > n_is_min = (~is_min);

	sc_biguint< 64 > mask_A = (is_min, n_is_min );

	sc_bigint< 32 * Q > min = VECTOR_MUX< 32 , Q >( mb, ma, is_min );

	sc_biguint< 32 > i_mask = SPC_SPC2_Min_Mask_32_CA2< Q >( min, sel );
	sc_biguint< 64 > mask_B = (i_mask, i_mask);

	sc_biguint< 64 > mask = VECTOR_AND< 64 >( mask_A, mask_B );
	return mask;
}

template <int Q>
sc_biguint< 128 > SPC_SPC2_Min_Mask_128_CA2 (sc_bigint< 128 * Q> Mag, bool sel)
{
#pragma HLS INLINE
	sc_bigint< 64 * Q > ma  = Mag.range( (64 * Q) - 1, 0 );
	sc_bigint< 64 * Q > mb  = Mag.range( (128 * Q) - 1, 64 * Q );

	sc_biguint< 64 > is_min = VECTOR_IS_MIN< 64, Q >( mb, ma );
	sc_biguint< 64 > n_is_min = (~is_min);

	sc_biguint< 128 > mask_A = (is_min, n_is_min );

	sc_bigint< 64 * Q > min = VECTOR_MUX< 64 , Q >( mb, ma, is_min );

	sc_biguint< 64 > i_mask = SPC_SPC2_Min_Mask_64_CA2< Q >( min, sel );
	sc_biguint< 128 > mask_B = (i_mask, i_mask);

	sc_biguint< 128 > mask = VECTOR_AND< 128 >( mask_A, mask_B );
	return mask;
}

template <int Q>
sc_biguint< 256 > SPC_SPC2_Min_Mask_256_CA2 (sc_bigint< 256 * Q> Mag, bool sel)
{
#pragma HLS INLINE
	sc_bigint< 128 * Q > ma  = Mag.range( (128 * Q) - 1, 0 );
	sc_bigint< 128 * Q > mb  = Mag.range( (256 * Q) - 1, 128 * Q );

	sc_biguint< 128 > is_min = VECTOR_IS_MIN< 128, Q >( mb, ma );
	sc_biguint< 128 > n_is_min = (~is_min);

	sc_biguint< 256 > mask_A = (is_min, n_is_min );

	sc_bigint< 128 * Q > min = VECTOR_MUX< 128 , Q >( mb, ma, is_min );

	sc_biguint< 128 > i_mask = SPC_SPC2_Min_Mask_128_CA2< Q >( min, sel );
	sc_biguint< 256 > mask_B = (i_mask, i_mask);

	sc_biguint< 256 > mask = VECTOR_AND< 256 >( mask_A, mask_B );
	return mask;
}

/////////////// SIGMAG ////////////////

template <int Q>
sc_biguint< 2 > SPC_SPC2_Min_Mask_2_SM (sc_biguint< 2 * Q> Mag, bool sel) // Sel 0 : SPC, Sel 1 : SPC2
{
#pragma HLS INLINE
	sc_biguint< Q > ma  = Mag.range( Q - 1, 0 );
	sc_biguint< Q > mb  = Mag.range( 2 * Q - 1, Q );

	sc_biguint< 1 > is_min = VECTOR_IS_MIN_SM< 1, Q >(mb, ma);
	sc_biguint< 1 > n_is_min = (~is_min);

	sc_biguint< 2 > mask1 = (is_min, n_is_min);

	sc_biguint< 2 > mask2 = 0x03; // (11)

	if(sel == 0)
		return mask1;
	else
		return mask2;
}

template <int Q>
sc_biguint< 4 > SPC_SPC2_Min_Mask_4_SM (sc_biguint< 4 * Q> Mag, bool sel)
{
#pragma HLS INLINE
	sc_biguint< 2 * Q > ma  = Mag.range( (2 * Q) - 1, 0 );
	sc_biguint< 2 * Q > mb  = Mag.range( (4 * Q) - 1, 2 * Q );

	sc_biguint< 2 > is_min = VECTOR_IS_MIN_SM< 2, Q >( mb, ma );
	sc_biguint< 2 > n_is_min = (~is_min);

	sc_biguint< 4 > mask_A = (is_min, n_is_min );

	sc_biguint< 2 * Q > min = VECTOR_MUX_SM< 2 , Q >( mb, ma, is_min );

	sc_biguint< 2 > i_mask = SPC_SPC2_Min_Mask_2_SM< Q >( min, sel );
	sc_biguint< 4 > mask_B = (i_mask, i_mask);

	sc_biguint< 4 > mask = VECTOR_AND< 4 >( mask_A, mask_B );
	return mask;
}

template <int Q>
sc_biguint< 8 > SPC_SPC2_Min_Mask_8_SM (sc_biguint< 8 * Q> Mag, bool sel)
{
#pragma HLS INLINE
	sc_biguint< 4 * Q > ma  = Mag.range( (4 * Q) - 1, 0 );
	sc_biguint< 4 * Q > mb  = Mag.range( (8 * Q) - 1, 4 * Q );

	sc_biguint< 4 > is_min = VECTOR_IS_MIN_SM< 4, Q >( mb, ma );
	sc_biguint< 4 > n_is_min = (~is_min);

	sc_biguint< 8 > mask_A = (is_min, n_is_min );

	sc_biguint< 4 * Q > min = VECTOR_MUX_SM< 4 , Q >( mb, ma, is_min );

	sc_biguint< 4 > i_mask = SPC_SPC2_Min_Mask_4_SM< Q >( min, sel );
	sc_biguint< 8 > mask_B = (i_mask, i_mask);

	sc_biguint< 8 > mask = VECTOR_AND< 8 >( mask_A, mask_B );
	return mask;
}

template <int Q>
sc_biguint< 16 > SPC_SPC2_Min_Mask_16_SM (sc_biguint< 16 * Q> Mag, bool sel)
{
#pragma HLS INLINE
	sc_biguint< 8 * Q > ma  = Mag.range( (8 * Q) - 1, 0 );
	sc_biguint< 8 * Q > mb  = Mag.range( (16 * Q) - 1, 8 * Q );

	sc_biguint< 8 > is_min = VECTOR_IS_MIN_SM< 8, Q >( mb, ma );
	sc_biguint< 8 > n_is_min = (~is_min);

	sc_biguint< 16 > mask_A = (is_min, n_is_min );

	sc_biguint< 8 * Q > min = VECTOR_MUX_SM< 8 , Q >( mb, ma, is_min );

	sc_biguint< 8 > i_mask = SPC_SPC2_Min_Mask_8_SM< Q >( min, sel );
	sc_biguint< 16 > mask_B = (i_mask, i_mask);

	sc_biguint< 16 > mask = VECTOR_AND< 16 >( mask_A, mask_B );
	return mask;
}

template <int Q>
sc_biguint< 32 > SPC_SPC2_Min_Mask_32_SM (sc_biguint< 32 * Q> Mag, bool sel)
{
#pragma HLS INLINE
	sc_biguint< 16 * Q > ma  = Mag.range( (16 * Q) - 1, 0 );
	sc_biguint< 16 * Q > mb  = Mag.range( (32 * Q) - 1, 16 * Q );

	sc_biguint< 16 > is_min = VECTOR_IS_MIN_SM< 16, Q >( mb, ma );
	sc_biguint< 16 > n_is_min = (~is_min);

	sc_biguint< 32 > mask_A = (is_min, n_is_min );

	sc_biguint< 16 * Q > min = VECTOR_MUX_SM< 16 , Q >( mb, ma, is_min );

	sc_biguint< 16 > i_mask = SPC_SPC2_Min_Mask_16_SM< Q >( min, sel );
	sc_biguint< 32 > mask_B = (i_mask, i_mask);

	sc_biguint< 32 > mask = VECTOR_AND< 32 >( mask_A, mask_B );
	return mask;
}

template <int Q>
sc_biguint< 64 > SPC_SPC2_Min_Mask_64_SM (sc_biguint< 64 * Q> Mag, bool sel)
{
#pragma HLS INLINE
	sc_biguint< 32 * Q > ma  = Mag.range( (32 * Q) - 1, 0 );
	sc_biguint< 32 * Q > mb  = Mag.range( (64 * Q) - 1, 32 * Q );

	sc_biguint< 32 > is_min = VECTOR_IS_MIN_SM< 32, Q >( mb, ma );
	sc_biguint< 32 > n_is_min = (~is_min);

	sc_biguint< 64 > mask_A = (is_min, n_is_min );

	sc_biguint< 32 * Q > min = VECTOR_MUX_SM< 32 , Q >( mb, ma, is_min );

	sc_biguint< 32 > i_mask = SPC_SPC2_Min_Mask_32_SM< Q >( min, sel );
	sc_biguint< 64 > mask_B = (i_mask, i_mask);

	sc_biguint< 64 > mask = VECTOR_AND< 64 >( mask_A, mask_B );
	return mask;
}

template <int Q>
sc_biguint< 128 > SPC_SPC2_Min_Mask_128_SM (sc_biguint< 128 * Q> Mag, bool sel)
{
#pragma HLS INLINE
	sc_biguint< 64 * Q > ma  = Mag.range( (64 * Q) - 1, 0 );
	sc_biguint< 64 * Q > mb  = Mag.range( (128 * Q) - 1, 64 * Q );

	sc_biguint< 64 > is_min = VECTOR_IS_MIN_SM< 64, Q >( mb, ma );
	sc_biguint< 64 > n_is_min = (~is_min);

	sc_biguint< 128 > mask_A = (is_min, n_is_min );

	sc_biguint< 64 * Q > min = VECTOR_MUX_SM< 64 , Q >( mb, ma, is_min );

	sc_biguint< 64 > i_mask = SPC_SPC2_Min_Mask_64_SM< Q >( min, sel );
	sc_biguint< 128 > mask_B = (i_mask, i_mask);

	sc_biguint< 128 > mask = VECTOR_AND< 128 >( mask_A, mask_B );
	return mask;
}

template <int Q>
sc_biguint< 256 > SPC_SPC2_Min_Mask_256_SM (sc_biguint< 256 * Q> Mag, bool sel)
{
#pragma HLS INLINE
	sc_biguint< 128 * Q > ma  = Mag.range( (128 * Q) - 1, 0 );
	sc_biguint< 128 * Q > mb  = Mag.range( (256 * Q) - 1, 128 * Q );

	sc_biguint< 128 > is_min = VECTOR_IS_MIN_SM< 128, Q >( mb, ma );
	sc_biguint< 128 > n_is_min = (~is_min);

	sc_biguint< 256 > mask_A = (is_min, n_is_min );

	sc_biguint< 128 * Q > min = VECTOR_MUX_SM< 128 , Q >( mb, ma, is_min );

	sc_biguint< 128 > i_mask = SPC_SPC2_Min_Mask_128_SM< Q >( min, sel );
	sc_biguint< 256 > mask_B = (i_mask, i_mask);

	sc_biguint< 256 > mask = VECTOR_AND< 256 >( mask_A, mask_B );
	return mask;
}

/////////////// GLOBAL SPC Node ////////////////
template <int Q>
sc_biguint< 2 > SPC_SPC2_Node_2 (sc_bigint< 2 * Q> llr, bool sel)
{
#pragma HLS INLINE

#if defined CA2
	sc_biguint< 2 > sign = VECTOR_SIGN< 2, Q >(llr);
	sc_bigint < 2 * Q > mag = VECTOR_ABS< 2, Q >(llr);

	sc_biguint< 2 > parity = SPC_SPC2_Parity_2( sign, sel );
	sc_biguint< 2 > mask = SPC_SPC2_Min_Mask_2_CA2< Q >( mag, sel );

	sc_biguint< 2 > and_mask = VECTOR_AND< 2 >( parity, mask) ;
	sc_biguint< 2 > bit = VECTOR_XOR< 2 >( sign, and_mask );
	return bit;
#elif defined SIGMAG
	sc_biguint< 2 > sign = VECTOR_SIGN_SM< 2, Q >( (sc_biguint< 2 * Q>) llr);
	sc_biguint< 2 * (Q-1) > mag = VECTOR_ABS_SM< 2, Q >( (sc_biguint< 2 * Q>) llr);

	sc_biguint< 2 > parity = SPC_SPC2_Parity_2( sign, sel );
	sc_biguint< 2 > mask = SPC_SPC2_Min_Mask_2_SM< Q-1 >( mag, sel );

	sc_biguint< 2 > and_mask = VECTOR_AND< 2 >( parity, mask) ;
	sc_biguint< 2 > bit = VECTOR_XOR< 2 >( sign, and_mask );
	return bit;
#endif
}

template <int Q>
sc_biguint< 4 > SPC_SPC2_Node_4 (sc_bigint< 4 * Q> llr, bool sel)
{
#pragma HLS INLINE

#if defined CA2
	sc_biguint< 4 > sign = VECTOR_SIGN< 4, Q >(llr);
	sc_bigint < 4 * Q > mag = VECTOR_ABS< 4, Q >(llr);

	sc_biguint< 4 > parity = SPC_SPC2_Parity_4( sign, sel );
	sc_biguint< 4 > mask = SPC_SPC2_Min_Mask_4_CA2< Q >( mag, sel );

	sc_biguint< 4 > and_mask = VECTOR_AND< 4 >( parity, mask) ;
	sc_biguint< 4 > bit = VECTOR_XOR< 4 >( sign, and_mask );
	return bit;
#elif defined SIGMAG
	sc_biguint< 4 > sign = VECTOR_SIGN_SM< 4, Q >( (sc_biguint< 4 * Q>) llr);
	sc_biguint< 4 * (Q-1) > mag = VECTOR_ABS_SM< 4, Q >( (sc_biguint< 4 * Q>) llr);

	sc_biguint< 4 > parity = SPC_SPC2_Parity_4( sign, sel );
	sc_biguint< 4 > mask = SPC_SPC2_Min_Mask_4_SM< Q-1 >( mag, sel );

	sc_biguint< 4 > and_mask = VECTOR_AND< 4 >( parity, mask) ;
	sc_biguint< 4 > bit = VECTOR_XOR< 4 >( sign, and_mask );
	return bit;
#endif
}

template <int Q>
sc_biguint< 8 > SPC_SPC2_Node_8 (sc_bigint< 8 * Q> llr, bool sel)
{
#pragma HLS INLINE

#if defined CA2
	sc_biguint< 8 > sign = VECTOR_SIGN< 8, Q >(llr);
	sc_bigint < 8 * Q > mag = VECTOR_ABS< 8, Q >(llr);

	sc_biguint< 8 > parity = SPC_SPC2_Parity_8( sign, sel );
	sc_biguint< 8 > mask = SPC_SPC2_Min_Mask_8_CA2< Q >( mag, sel );

	sc_biguint< 8 > and_mask = VECTOR_AND< 8 >( parity, mask) ;
	sc_biguint< 8 > bit = VECTOR_XOR< 8 >( sign, and_mask );
	return bit;
#elif defined SIGMAG
	sc_biguint< 8 > sign = VECTOR_SIGN_SM< 8, Q >( (sc_biguint< 8 * Q>) llr);
	sc_biguint< 8 * (Q-1) > mag = VECTOR_ABS_SM< 8, Q >( (sc_biguint< 8 * Q>) llr);

	sc_biguint< 8 > parity = SPC_SPC2_Parity_8( sign, sel );
	sc_biguint< 8 > mask = SPC_SPC2_Min_Mask_8_SM< Q-1 >( mag, sel );

	sc_biguint< 8 > and_mask = VECTOR_AND< 8 >( parity, mask) ;
	sc_biguint< 8 > bit = VECTOR_XOR< 8 >( sign, and_mask );
	return bit;
#endif
}

template <int Q>
sc_biguint< 16 > SPC_SPC2_Node_16 (sc_bigint< 16 * Q> llr, bool sel)
{
#pragma HLS INLINE

#if defined CA2
	sc_biguint< 16 > sign = VECTOR_SIGN< 16, Q >(llr);
	sc_bigint < 16 * Q > mag = VECTOR_ABS< 16, Q >(llr);

	sc_biguint< 16 > parity = SPC_SPC2_Parity_16( sign, sel );
	sc_biguint< 16 > mask = SPC_SPC2_Min_Mask_16_CA2< Q >( mag, sel );

	sc_biguint< 16 > and_mask = VECTOR_AND< 16 >( parity, mask) ;
	sc_biguint< 16 > bit = VECTOR_XOR< 16 >( sign, and_mask );
	return bit;
#elif defined SIGMAG
	sc_biguint< 16 > sign = VECTOR_SIGN_SM< 16, Q >( (sc_biguint< 16 * Q>) llr);
	sc_biguint< 16 * (Q-1) > mag = VECTOR_ABS_SM< 16, Q >( (sc_biguint< 16 * Q>) llr);

	sc_biguint< 16 > parity = SPC_SPC2_Parity_16( sign, sel );
	sc_biguint< 16 > mask = SPC_SPC2_Min_Mask_16_SM< Q-1 >( mag, sel );

	sc_biguint< 16 > and_mask = VECTOR_AND< 16 >( parity, mask) ;
	sc_biguint< 16 > bit = VECTOR_XOR< 16 >( sign, and_mask );
	return bit;
#endif
}

template <int Q>
sc_biguint< 32 > SPC_SPC2_Node_32 (sc_bigint< 32 * Q> llr, bool sel)
{
#pragma HLS INLINE

#if defined CA2
	sc_biguint< 32 > sign = VECTOR_SIGN< 32, Q >(llr);
	sc_bigint < 32 * Q > mag = VECTOR_ABS< 32, Q >(llr);

	sc_biguint< 32 > parity = SPC_SPC2_Parity_32( sign, sel );
	sc_biguint< 32 > mask = SPC_SPC2_Min_Mask_32_CA2< Q >( mag, sel );

	sc_biguint< 32 > and_mask = VECTOR_AND< 32 >( parity, mask) ;
	sc_biguint< 32 > bit = VECTOR_XOR< 32 >( sign, and_mask );
	return bit;
#elif defined SIGMAG
	sc_biguint< 32 > sign = VECTOR_SIGN_SM< 32, Q >( (sc_biguint< 32 * Q>) llr);
	sc_biguint< 32 * (Q-1) > mag = VECTOR_ABS_SM< 32, Q >( (sc_biguint< 32 * Q>) llr);

	sc_biguint< 32 > parity = SPC_SPC2_Parity_32( sign, sel );
	sc_biguint< 32 > mask = SPC_SPC2_Min_Mask_32_SM< Q-1 >( mag, sel );

	sc_biguint< 32 > and_mask = VECTOR_AND< 32 >( parity, mask) ;
	sc_biguint< 32 > bit = VECTOR_XOR< 32 >( sign, and_mask );
	return bit;
#endif
}

template <int Q>
sc_biguint< 64 > SPC_SPC2_Node_64 (sc_bigint< 64 * Q> llr, bool sel)
{
#pragma HLS INLINE

#if defined CA2
	sc_biguint< 64 > sign = VECTOR_SIGN< 64, Q >(llr);
	sc_bigint < 64 * Q > mag = VECTOR_ABS< 64, Q >(llr);

	sc_biguint< 64 > parity = SPC_SPC2_Parity_64( sign, sel );
	sc_biguint< 64 > mask = SPC_SPC2_Min_Mask_64_CA2< Q >( mag, sel );

	sc_biguint< 64 > and_mask = VECTOR_AND< 64 >( parity, mask) ;
	sc_biguint< 64 > bit = VECTOR_XOR< 64 >( sign, and_mask );
	return bit;
#elif defined SIGMAG
	sc_biguint< 64 > sign = VECTOR_SIGN_SM< 64, Q >( (sc_biguint< 64 * Q>) llr);
	sc_biguint< 64 * (Q-1) > mag = VECTOR_ABS_SM< 64, Q >( (sc_biguint< 64 * Q>) llr);

	sc_biguint< 64 > parity = SPC_SPC2_Parity_64( sign, sel );
	sc_biguint< 64 > mask = SPC_SPC2_Min_Mask_64_SM< Q-1 >( mag, sel );

	sc_biguint< 64 > and_mask = VECTOR_AND< 64 >( parity, mask) ;
	sc_biguint< 64 > bit = VECTOR_XOR< 64 >( sign, and_mask );
	return bit;
#endif
}

template <int Q>
sc_biguint< 128 > SPC_SPC2_Node_128 (sc_bigint< 128 * Q> llr, bool sel)
{
#pragma HLS INLINE

#if defined CA2
	sc_biguint< 128 > sign = VECTOR_SIGN< 128, Q >(llr);
	sc_bigint < 128 * Q > mag = VECTOR_ABS< 128, Q >(llr);

	sc_biguint< 128 > parity = SPC_SPC2_Parity_128( sign, sel );
	sc_biguint< 128 > mask = SPC_SPC2_Min_Mask_128_CA2< Q >( mag, sel );

	sc_biguint< 128 > and_mask = VECTOR_AND< 128 >( parity, mask) ;
	sc_biguint< 128 > bit = VECTOR_XOR< 128 >( sign, and_mask );
	return bit;
#elif defined SIGMAG
	sc_biguint< 128 > sign = VECTOR_SIGN_SM< 128, Q >( (sc_biguint< 128 * Q>) llr);
	sc_biguint< 128 * (Q-1) > mag = VECTOR_ABS_SM< 128, Q >( (sc_biguint< 128 * Q>) llr);

	sc_biguint< 128 > parity = SPC_SPC2_Parity_128( sign, sel );
	sc_biguint< 128 > mask = SPC_SPC2_Min_Mask_128_SM< Q-1 >( mag, sel );

	sc_biguint< 128 > and_mask = VECTOR_AND< 128 >( parity, mask) ;
	sc_biguint< 128 > bit = VECTOR_XOR< 128 >( sign, and_mask );
	return bit;
#endif
}

template <int Q>
sc_biguint< 256 > SPC_SPC2_Node_256 (sc_bigint< 256 * Q> llr, bool sel)
{
#pragma HLS INLINE

#if defined CA2
	sc_biguint< 256 > sign = VECTOR_SIGN< 256, Q >(llr);
	sc_bigint < 256 * Q > mag = VECTOR_ABS< 256, Q >(llr);

	sc_biguint< 256 > parity = SPC_SPC2_Parity_256( sign, sel );
	sc_biguint< 256 > mask = SPC_SPC2_Min_Mask_256_CA2< Q >( mag, sel );

	sc_biguint< 256 > and_mask = VECTOR_AND< 256 >( parity, mask) ;
	sc_biguint< 256 > bit = VECTOR_XOR< 256 >( sign, and_mask );
	return bit;
#elif defined SIGMAG
	sc_biguint< 256 > sign = VECTOR_SIGN_SM< 256, Q >( (sc_biguint< 256 * Q>) llr);
	sc_biguint< 256 * (Q-1) > mag = VECTOR_ABS_SM< 256, Q >( (sc_biguint< 256 * Q>) llr);

	sc_biguint< 256 > parity = SPC_SPC2_Parity_256( sign, sel );
	sc_biguint< 256 > mask = SPC_SPC2_Min_Mask_256_SM< Q-1 >( mag, sel );

	sc_biguint< 256 > and_mask = VECTOR_AND< 256 >( parity, mask) ;
	sc_biguint< 256 > bit = VECTOR_XOR< 256 >( sign, and_mask );
	return bit;
#endif
}

//*************************************************************************//
// 							PRUNING FUNCTIONS REP						   //
//*************************************************************************//

/////////////// CA2 ////////////////

template <int Q>
sc_bigint< Q + 1 > ADD_TREE_2_CA2 (sc_bigint< 2 * Q > llr)
{
#pragma HLS INLINE
	sc_bigint< Q 	 > la  = llr.range( Q - 1, 0 );
	sc_bigint< Q 	 > lb  = llr.range( 2 * Q - 1, Q );
	sc_bigint< Q + 1 > sum = VECTOR_ADD_NOSAT < 1, Q > ( la, lb );

	return sum;
}

template <int Q>
sc_bigint< Q + 2 > ADD_TREE_4_CA2 (sc_bigint< 4 * Q> llr)
{
#pragma HLS INLINE
	sc_bigint< 2 * Q 	   > la  = llr.range( (2 * Q) - 1, 0 );
	sc_bigint< 2 * Q 	   > lb  = llr.range( (4 * Q) - 1, (2 * Q) );
	sc_bigint< 2 * (Q + 1) > sum = VECTOR_ADD_NOSAT < 2, Q > ( la, lb );

	sc_bigint< Q + 2 > res = ADD_TREE_2_CA2 < Q + 1 > ( sum );

	return res;
}

template <int Q>
sc_bigint< Q + 3 > ADD_TREE_8_CA2 (sc_bigint< 8 * Q> llr)
{
#pragma HLS INLINE
	sc_bigint< 4 * Q 	   > la  = llr.range( (4 * Q) - 1, 0 );
	sc_bigint< 4 * Q 	   > lb  = llr.range( (8 * Q) - 1, (4 * Q) );
	sc_bigint< 4 * (Q + 1) > sum = VECTOR_ADD_NOSAT < 4, Q > ( la, lb );

	sc_bigint< Q + 3 > res = ADD_TREE_4_CA2 < Q + 1 > ( sum );

	return res;
}

template <int Q>
sc_bigint< Q + 4 > ADD_TREE_16_CA2 (sc_bigint< 16 * Q> llr)
{
#pragma HLS INLINE
	sc_bigint< 8 * Q 	   > la  = llr.range( (8 * Q) - 1, 0 );
	sc_bigint< 8 * Q 	   > lb  = llr.range( (16 * Q) - 1, (8 * Q) );
	sc_bigint< 8 * (Q + 1) > sum = VECTOR_ADD_NOSAT < 8, Q > ( la, lb );

	sc_bigint< Q + 4 > res = ADD_TREE_8_CA2 < Q + 1 > ( sum );

	return res;
}

template <int Q>
sc_bigint< Q + 5 > ADD_TREE_32_CA2 (sc_bigint< 32 * Q> llr)
{
#pragma HLS INLINE
	sc_bigint< 16 * Q 	   > la  = llr.range( (16 * Q) - 1, 0 );
	sc_bigint< 16 * Q 	   > lb  = llr.range( (32 * Q) - 1, (16 * Q) );
	sc_bigint< 16 * (Q + 1) > sum = VECTOR_ADD_NOSAT < 16, Q > ( la, lb );

	sc_bigint< Q + 5 > res = ADD_TREE_16_CA2 < Q + 1 > ( sum );

	return res;
}

template <int Q>
sc_bigint< Q + 6 > ADD_TREE_64_CA2 (sc_bigint< 64 * Q> llr)
{
#pragma HLS INLINE
	sc_bigint< 32 * Q 	   > la  = llr.range( (32 * Q) - 1, 0 );
	sc_bigint< 32 * Q 	   > lb  = llr.range( (64 * Q) - 1, (32 * Q) );
	sc_bigint< 32 * (Q + 1) > sum = VECTOR_ADD_NOSAT < 32, Q > ( la, lb );

	sc_bigint< Q + 6 > res = ADD_TREE_32_CA2 < Q + 1 > ( sum );

	return res;
}

template <int Q>
sc_bigint< Q + 7 > ADD_TREE_128_CA2 (sc_bigint< 128 * Q> llr)
{
#pragma HLS INLINE
	sc_bigint< 64 * Q 	   > la  = llr.range( (64 * Q) - 1, 0 );
	sc_bigint< 64 * Q 	   > lb  = llr.range( (128 * Q) - 1, (64 * Q) );
	sc_bigint< 64 * (Q + 1) > sum = VECTOR_ADD_NOSAT < 64, Q > ( la, lb );

	sc_bigint< Q + 7 > res = ADD_TREE_64_CA2 < Q + 1 > ( sum );

	return res;
}

template <int Q>
sc_bigint< Q + 8 > ADD_TREE_256_CA2 (sc_bigint< 256 * Q> llr)
{
#pragma HLS INLINE
	sc_bigint< 128 * Q 	   > la  = llr.range( (128 * Q) - 1, 0 );
	sc_bigint< 128 * Q 	   > lb  = llr.range( (256 * Q) - 1, (128 * Q) );
	sc_bigint< 128 * (Q + 1) > sum = VECTOR_ADD_NOSAT < 128, Q > ( la, lb );

	sc_bigint< Q + 8 > res = ADD_TREE_128_CA2 < Q + 1 > ( sum );

	return res;
}

/////////////// SIGMAG ////////////////

template <int Q>
sc_biguint< Q + 1 > ADD_TREE_2_SM (sc_biguint< 2 * Q> llr)
{
#pragma HLS INLINE
	sc_biguint< Q 	 > la  = llr.range( Q - 1, 0 );
	sc_biguint< Q 	 > lb  = llr.range( 2 * Q - 1, Q );
	sc_biguint< Q + 1 > sum  = VECTOR_FULL_ADDER_SM < 1, Q > ( la, lb );

	return sum;
}

template <int Q>
sc_biguint< Q + 2 > ADD_TREE_4_SM (sc_biguint< 4 * Q> llr)
{
#pragma HLS INLINE
	sc_biguint< 2 * Q 	   > la  = llr.range( (2 * Q) - 1, 0 );
	sc_biguint< 2 * Q 	   > lb  = llr.range( (4 * Q) - 1, (2 * Q) );
	sc_biguint< 2 * (Q + 1) > sum = VECTOR_FULL_ADDER_SM < 2, Q > ( la, lb );

	sc_biguint< Q + 2 > res = ADD_TREE_2_SM < Q + 1 > ( sum );

	return res;
}

template <int Q>
sc_biguint< Q + 3 > ADD_TREE_8_SM (sc_biguint< 8 * Q> llr)
{
#pragma HLS INLINE
	sc_biguint< 4 * Q 	   > la  = llr.range( (4 * Q) - 1, 0 );
	sc_biguint< 4 * Q 	   > lb  = llr.range( (8 * Q) - 1, (4 * Q) );
	sc_biguint< 4 * (Q + 1) > sum = VECTOR_FULL_ADDER_SM < 4, Q > ( la, lb );

	sc_biguint< Q + 3 > res = ADD_TREE_4_SM < Q + 1 > ( sum );

	return res;
}

template <int Q>
sc_biguint< Q + 4 > ADD_TREE_16_SM (sc_biguint< 16 * Q> llr)
{
#pragma HLS INLINE
	sc_biguint< 8 * Q 	   > la  = llr.range( (8 * Q) - 1, 0 );
	sc_biguint< 8 * Q 	   > lb  = llr.range( (16 * Q) - 1, (8 * Q) );
	sc_biguint< 8 * (Q + 1) > sum = VECTOR_FULL_ADDER_SM < 8, Q > ( la, lb );

	sc_biguint< Q + 4 > res = ADD_TREE_8_SM < Q + 1 > ( sum );

	return res;
}

template <int Q>
sc_biguint< Q + 5 > ADD_TREE_32_SM (sc_biguint< 32 * Q> llr)
{
#pragma HLS INLINE
	sc_biguint< 16 * Q 	   > la  = llr.range( (16 * Q) - 1, 0 );
	sc_biguint< 16 * Q 	   > lb  = llr.range( (32 * Q) - 1, (16 * Q) );
	sc_biguint< 16 * (Q + 1) > sum = VECTOR_FULL_ADDER_SM < 16, Q > ( la, lb );

	sc_biguint< Q + 5 > res = ADD_TREE_16_SM < Q + 1 > ( sum );

	return res;
}

template <int Q>
sc_biguint< Q + 6 > ADD_TREE_64_SM (sc_biguint< 64 * Q> llr)
{
#pragma HLS INLINE
	sc_biguint< 32 * Q 	   > la  = llr.range( (32 * Q) - 1, 0 );
	sc_biguint< 32 * Q 	   > lb  = llr.range( (64 * Q) - 1, (32 * Q) );
	sc_biguint< 32 * (Q + 1) > sum = VECTOR_FULL_ADDER_SM < 32, Q > ( la, lb );

	sc_biguint< Q + 6 > res = ADD_TREE_32_SM < Q + 1 > ( sum );

	return res;
}

template <int Q>
sc_biguint< Q + 7 > ADD_TREE_128_SM (sc_biguint< 128 * Q> llr)
{
#pragma HLS INLINE
	sc_biguint< 64 * Q 	   > la  = llr.range( (64 * Q) - 1, 0 );
	sc_biguint< 64 * Q 	   > lb  = llr.range( (128 * Q) - 1, (64 * Q) );
	sc_biguint< 64 * (Q + 1) > sum = VECTOR_FULL_ADDER_SM < 64, Q > ( la, lb );

	sc_biguint< Q + 7 > res = ADD_TREE_64_SM < Q + 1 > ( sum );

	return res;
}

template <int Q>
sc_biguint< Q + 8 > ADD_TREE_256_SM (sc_biguint< 256 * Q> llr)
{
#pragma HLS INLINE
	sc_biguint< 128 * Q 	   > la  = llr.range( (128 * Q) - 1, 0 );
	sc_biguint< 128 * Q 	   > lb  = llr.range( (256 * Q) - 1, (128 * Q) );
	sc_biguint< 128 * (Q + 1) > sum = VECTOR_FULL_ADDER_SM < 128, Q > ( la, lb );

	sc_biguint< Q + 8 > res = ADD_TREE_128_SM < Q + 1 > ( sum );

	return res;
}

/////////////// GLOBAL ADDER TREE ////////////////

template <int Q>
sc_bigint< Q + 2 > ADDER_TREE_2 (sc_bigint< 2 * Q > llr, sc_bigint< Q + 2 > old_sum)
{
#pragma HLS INLINE

#if defined CA2
	sc_bigint< Q + 2 > add_tree = (sc_bigint< Q + 2 >) ADD_TREE_2_CA2< Q >(llr);
	sc_bigint< Q + 2 > sum = VECTOR_ADD < 1, Q + 2 >(add_tree, old_sum);
	return sum;
#elif defined SIGMAG
	sc_biguint< Q + 1 > add_tree = ADD_TREE_2_SM< Q >( (sc_biguint< 2 * Q >) llr);
	sc_biguint< Q + 2 > add_tree_ext = ( (sc_biguint<1>)add_tree[Q] , (sc_biguint<Q + 1>)add_tree.range(Q-1, 0) );
	sc_biguint< Q + 2 > sum = VECTOR_FULL_ADDER_SAT_SM < 1, Q + 2 >(add_tree_ext, (sc_biguint< Q + 2 >) old_sum);
	return (sc_bigint< Q + 2 >) sum;
#endif
}

template <int Q>
sc_bigint< Q + 3 > ADDER_TREE_4 (sc_bigint< 4 * Q > llr, sc_bigint< Q + 3 > old_sum)
{
#pragma HLS INLINE

#if defined CA2
	sc_bigint< Q + 3 > add_tree = (sc_bigint< Q + 3 >) ADD_TREE_4_CA2< Q >(llr);
	sc_bigint< Q + 3 > sum = VECTOR_ADD < 1, Q + 3 >(add_tree, old_sum);
	return sum;
#elif defined SIGMAG
	sc_biguint< Q + 2 > add_tree = ADD_TREE_4_SM< Q >( (sc_biguint< 4 * Q >) llr);
	sc_biguint< Q + 3 > add_tree_ext = ( (sc_biguint<1>)add_tree[Q+1] , (sc_biguint<Q + 2>)add_tree.range(Q, 0) );
	sc_biguint< Q + 3 > sum = VECTOR_FULL_ADDER_SAT_SM < 1, Q + 3 >(add_tree_ext, (sc_biguint< Q + 3 >) old_sum);
	return (sc_bigint< Q + 3 >) sum;
#endif
}

template <int Q>
sc_bigint< Q + 4 > ADDER_TREE_8 (sc_bigint< 8 * Q > llr, sc_bigint< Q + 4 > old_sum)
{
#pragma HLS INLINE

#if defined CA2
	sc_bigint< Q + 4 > add_tree = (sc_bigint< Q + 4 >) ADD_TREE_8_CA2< Q >(llr);
	sc_bigint< Q + 4 > sum = VECTOR_ADD < 1, Q + 4 >(add_tree, old_sum);
	return sum;
#elif defined SIGMAG
	sc_biguint< Q + 3 > add_tree = ADD_TREE_8_SM< Q >( (sc_biguint< 8 * Q >) llr);
	sc_biguint< Q + 4 > add_tree_ext = ( (sc_biguint<1>)add_tree[Q+2] , (sc_biguint<Q + 3>)add_tree.range(Q+1, 0) );
	sc_biguint< Q + 4 > sum = VECTOR_FULL_ADDER_SAT_SM < 1, Q + 4 >(add_tree_ext, (sc_biguint< Q + 4 >) old_sum);
	return (sc_bigint< Q + 4 >) sum;
#endif
}

template <int Q>
sc_bigint< Q + 5 > ADDER_TREE_16 (sc_bigint< 16 * Q > llr, sc_bigint< Q + 5 > old_sum)
{
#pragma HLS INLINE

#if defined CA2
	sc_bigint< Q + 5 > add_tree = (sc_bigint< Q + 5 >) ADD_TREE_16_CA2< Q >(llr);
	sc_bigint< Q + 5 > sum = VECTOR_ADD < 1, Q + 5 >(add_tree, old_sum);
	return sum;
#elif defined SIGMAG
	sc_biguint< Q + 4 > add_tree = ADD_TREE_16_SM< Q >( (sc_biguint< 16 * Q >) llr);
	sc_biguint< Q + 5 > add_tree_ext = ( (sc_biguint<1>)add_tree[Q + 3] , (sc_biguint<Q + 4>)add_tree.range(Q + 2, 0) );
	sc_biguint< Q + 5 > sum = VECTOR_FULL_ADDER_SAT_SM < 1, Q + 5 >(add_tree_ext, (sc_biguint< Q + 5 >) old_sum);
	return (sc_bigint< Q + 5 >) sum;
#endif
}

template <int Q>
sc_bigint< Q + 6> ADDER_TREE_32 (sc_bigint< 32 * Q > llr, sc_bigint< Q + 6> old_sum)
{
#pragma HLS INLINE

#if defined CA2
	sc_bigint< Q + 6> add_tree = (sc_bigint< Q + 6>) ADD_TREE_32_CA2< Q >(llr);
	sc_bigint< Q + 6> sum = VECTOR_ADD < 1, Q + 6>(add_tree, old_sum);
	return sum;
#elif defined SIGMAG
	sc_biguint< Q + 5> add_tree = ADD_TREE_32_SM< Q >( (sc_biguint< 32 * Q >) llr);
	sc_biguint< Q + 6> add_tree_ext = ( (sc_biguint<1>)add_tree[Q + 4] , (sc_biguint<Q + 5>)add_tree.range(Q + 3, 0) );
	sc_biguint< Q + 6> sum = VECTOR_FULL_ADDER_SAT_SM < 1, Q + 6>(add_tree_ext, (sc_biguint< Q + 6>) old_sum);
	return (sc_bigint< Q + 6>) sum;
#endif
}

template <int Q>
sc_bigint< Q + 7> ADDER_TREE_64 (sc_bigint< 64 * Q > llr, sc_bigint< Q + 7> old_sum)
{
#pragma HLS INLINE

#if defined CA2
	sc_bigint< Q + 7> add_tree = (sc_bigint< Q + 7>) ADD_TREE_64_CA2< Q >(llr);
	sc_bigint< Q + 7> sum = VECTOR_ADD < 1, Q + 7>(add_tree, old_sum);
	return sum;
#elif defined SIGMAG
	sc_biguint< Q + 6> add_tree = ADD_TREE_64_SM< Q >( (sc_biguint< 64 * Q >) llr);
	sc_biguint< Q + 7> add_tree_ext = ( (sc_biguint<1>)add_tree[Q + 5] , (sc_biguint<Q + 6>)add_tree.range(Q + 4, 0) );
	sc_biguint< Q + 7> sum = VECTOR_FULL_ADDER_SAT_SM < 1, Q + 7>(add_tree_ext, (sc_biguint< Q + 7>) old_sum);
	return (sc_bigint< Q + 7>) sum;
#endif
}

template <int Q>
sc_bigint< Q + 8> ADDER_TREE_128 (sc_bigint< 128 * Q > llr, sc_bigint< Q + 8> old_sum)
{
#pragma HLS INLINE

#if defined CA2
	sc_bigint< Q + 8> add_tree = (sc_bigint< Q + 8>) ADD_TREE_128_CA2< Q >(llr);
	sc_bigint< Q + 8> sum = VECTOR_ADD < 1, Q + 8>(add_tree, old_sum);
	return sum;
#elif defined SIGMAG
	sc_biguint< Q + 7> add_tree = ADD_TREE_128_SM< Q >( (sc_biguint< 128 * Q >) llr);
	sc_biguint< Q + 8> add_tree_ext = ( (sc_biguint<1>)add_tree[Q + 6] , (sc_biguint<Q + 7>)add_tree.range(Q + 5, 0) );
	sc_biguint< Q + 8> sum = VECTOR_FULL_ADDER_SAT_SM < 1, Q + 8>(add_tree_ext, (sc_biguint< Q + 8>) old_sum);
	return (sc_bigint< Q + 8>) sum;
#endif
}

template <int Q>
sc_bigint< Q + 9> ADDER_TREE_256 (sc_bigint< 256 * Q > llr, sc_bigint< Q + 9> old_sum)
{
#pragma HLS INLINE

#if defined CA2
	sc_bigint< Q + 9> add_tree = (sc_bigint< Q + 9>) ADD_TREE_256_CA2< Q >(llr);
	sc_bigint< Q + 9> sum = VECTOR_ADD < 1, Q + 9>(add_tree, old_sum);
	return sum;
#elif defined SIGMAG
	sc_biguint< Q + 8> add_tree = ADD_TREE_256_SM< Q >( (sc_biguint< 256 * Q >) llr);
	sc_biguint< Q + 9> add_tree_ext = ( (sc_biguint<1>)add_tree[Q + 7] , (sc_biguint<Q + 8>)add_tree.range(Q + 6, 0) );
	sc_biguint< Q + 9> sum = VECTOR_FULL_ADDER_SAT_SM < 1, Q + 9>(add_tree_ext, (sc_biguint< Q + 9>) old_sum);
	return (sc_bigint< Q + 9>) sum;
#endif
}

//*************************************************************************//
// 							PRUNING FUNCTIONS SPC						   //
//*************************************************************************//

inline sc_biguint< 1 > xor_TREE_2 (sc_biguint< 2 > sign)
{
#pragma HLS INLINE
	sc_biguint< 1 > sa = (sc_biguint< 1 >) sign[0];
	sc_biguint< 1 > sb = (sc_biguint< 1 >) sign[1];
	sc_biguint< 1 > parity = VECTOR_XOR< 1 >(sa, sb);

	return parity;
}

inline sc_biguint< 1 > xor_TREE_4 (sc_biguint< 4 > sign)
{
#pragma HLS INLINE
	sc_biguint< 2 > sa = sign.range( 1, 0 );
	sc_biguint< 2 > sb = sign.range( 3, 2 );
	sc_biguint< 2 > s_xor = VECTOR_XOR< 2 >(sa, sb);

	sc_biguint< 1 > parity = xor_TREE_2 (s_xor);

	return parity;
}

inline sc_biguint< 1 > xor_TREE_8 (sc_biguint< 8 > sign)
{
#pragma HLS INLINE
	sc_biguint< 4 > sa = sign.range( 3, 0 );
	sc_biguint< 4 > sb = sign.range( 7, 4 );
	sc_biguint< 4 > s_xor = VECTOR_XOR< 4 >(sa, sb);

	sc_biguint< 1 > parity = xor_TREE_4 (s_xor);

	return parity;
}

inline sc_biguint< 1 > xor_TREE_16 (sc_biguint< 16 > sign)
{
#pragma HLS INLINE
	sc_biguint< 8 > sa = sign.range( 7, 0 );
	sc_biguint< 8 > sb = sign.range( 15, 8 );
	sc_biguint< 8 > s_xor = VECTOR_XOR< 8 >(sa, sb);

	sc_biguint< 1 > parity = xor_TREE_8 (s_xor);

	return parity;
}

inline sc_biguint< 1 > xor_TREE_32 (sc_biguint< 32 > sign)
{
#pragma HLS INLINE
	sc_biguint< 16 > sa = sign.range( 15, 0 );
	sc_biguint< 16 > sb = sign.range( 31, 16 );
	sc_biguint< 16 > s_xor = VECTOR_XOR< 16 >(sa, sb);

	sc_biguint< 1 > parity = xor_TREE_16 (s_xor);

	return parity;
}

inline sc_biguint< 1 > xor_TREE_64 (sc_biguint< 64 > sign)
{
#pragma HLS INLINE
	sc_biguint< 32 > sa = sign.range( 31, 0 );
	sc_biguint< 32 > sb = sign.range( 63, 32 );
	sc_biguint< 32 > s_xor = VECTOR_XOR< 32 >(sa, sb);

	sc_biguint< 1 > parity = xor_TREE_32 (s_xor);

	return parity;
}

inline sc_biguint< 1 > xor_TREE_128 (sc_biguint< 128 > sign)
{
#pragma HLS INLINE
	sc_biguint< 64 > sa = sign.range( 63, 0 );
	sc_biguint< 64 > sb = sign.range( 127, 64 );
	sc_biguint< 64 > s_xor = VECTOR_XOR< 64 >(sa, sb);

	sc_biguint< 1 > parity = xor_TREE_64 (s_xor);

	return parity;
}

inline sc_biguint< 1 > xor_TREE_256 (sc_biguint< 256 > sign)
{
#pragma HLS INLINE
	sc_biguint< 128 > sa = sign.range( 127, 0 );
	sc_biguint< 128 > sb = sign.range( 255, 128 );
	sc_biguint< 128 > s_xor = VECTOR_XOR< 128 >(sa, sb);

	sc_biguint< 1 > parity = xor_TREE_128 (s_xor);

	return parity;
}

///////////////  GLOBAL Parity TREE ////////////////

inline sc_biguint< 1 > Parity_TREE_2 (sc_biguint< 2 > sign, sc_biguint< 1 > old_parity)
{
#pragma HLS INLINE

	sc_biguint< 1 > parity = xor_TREE_2 (sign);
	sc_biguint< 1 > s_xor = VECTOR_XOR< 1 >(parity, old_parity);
	return s_xor;
}

inline sc_biguint< 1 > Parity_TREE_4 (sc_biguint< 4 > sign, sc_biguint< 1 > old_parity)
{
#pragma HLS INLINE

	sc_biguint< 1 > parity = xor_TREE_4 (sign);
	sc_biguint< 1 > s_xor = VECTOR_XOR< 1 >(parity, old_parity);
	return s_xor;
}

inline sc_biguint< 1 > Parity_TREE_8 (sc_biguint< 8 > sign, sc_biguint< 1 > old_parity)
{
#pragma HLS INLINE

	sc_biguint< 1 > parity = xor_TREE_8 (sign);
	sc_biguint< 1 > s_xor = VECTOR_XOR< 1 >(parity, old_parity);
	return s_xor;
}

inline sc_biguint< 1 > Parity_TREE_16 (sc_biguint< 16 > sign, sc_biguint< 1 > old_parity)
{
#pragma HLS INLINE

	sc_biguint< 1 > parity = xor_TREE_16 (sign);
	sc_biguint< 1 > s_xor = VECTOR_XOR< 1 >(parity, old_parity);
	return s_xor;
}

inline sc_biguint< 1 > Parity_TREE_32 (sc_biguint< 32 > sign, sc_biguint< 1 > old_parity)
{
#pragma HLS INLINE

	sc_biguint< 1 > parity = xor_TREE_32 (sign);
	sc_biguint< 1 > s_xor = VECTOR_XOR< 1 >(parity, old_parity);
	return s_xor;
}

inline sc_biguint< 1 > Parity_TREE_64 (sc_biguint< 64 > sign, sc_biguint< 1 > old_parity)
{
#pragma HLS INLINE

	sc_biguint< 1 > parity = xor_TREE_64 (sign);
	sc_biguint< 1 > s_xor = VECTOR_XOR< 1 >(parity, old_parity);
	return s_xor;
}

inline sc_biguint< 1 > Parity_TREE_128 (sc_biguint< 128 > sign, sc_biguint< 1 > old_parity)
{
#pragma HLS INLINE

	sc_biguint< 1 > parity = xor_TREE_128 (sign);
	sc_biguint< 1 > s_xor = VECTOR_XOR< 1 >(parity, old_parity);
	return s_xor;
}

inline sc_biguint< 1 > Parity_TREE_256 (sc_biguint< 256 > sign, sc_biguint< 1 > old_parity)
{
#pragma HLS INLINE

	sc_biguint< 1 > parity = xor_TREE_256 (sign);
	sc_biguint< 1 > s_xor = VECTOR_XOR< 1 >(parity, old_parity);
	return s_xor;
}

/////////////// CA2 ////////////////

template <int Q>
sc_bigint< Q + 2 > Min_Mask_2_CA2 (sc_bigint< 2 * Q> Mag)
{
#pragma HLS INLINE
	sc_bigint< Q > ma  = Mag.range( Q - 1, 0 );
	sc_bigint< Q > mb  = Mag.range( 2 * Q - 1, Q );

	sc_biguint< 1 > is_min = VECTOR_IS_MIN< 1, Q >(mb, ma);
	sc_biguint< 1 > n_is_min = (~is_min);

	sc_bigint< Q > min = VECTOR_MUX< 1 , Q >( mb, ma, is_min );

	sc_biguint< 2 > mask = (is_min, n_is_min);

	sc_bigint< Q + 2 > res = ( min , mask);
	return res;
}

template <int Q>
sc_bigint< Q + 4 > Min_Mask_4_CA2 (sc_bigint< 4 * Q> Mag)
{
#pragma HLS INLINE
	sc_bigint< 2 * Q > ma  = Mag.range( (2 * Q) - 1, 0 );
	sc_bigint< 2 * Q > mb  = Mag.range( (4 * Q) - 1, 2 * Q );

	sc_biguint< 2 > is_min = VECTOR_IS_MIN< 2, Q >( mb, ma );
	sc_biguint< 2 > n_is_min = (~is_min);

	sc_biguint< 4 > mask_A = (is_min, n_is_min );

	sc_bigint< 2 * Q > min = VECTOR_MUX< 2 , Q >( mb, ma, is_min );

	sc_bigint< Q + 2 > min_mask = Min_Mask_2_CA2< Q >( min );
	sc_bigint< Q > res_min = min_mask.range(Q + 1 , 2);
	sc_biguint< 2 > i_mask = (sc_biguint< 2 >) min_mask.range(1,0);

	sc_biguint< 4 > mask_B = (i_mask, i_mask);

	sc_biguint< 4 > mask = VECTOR_AND< 4 >( mask_A, mask_B );

	sc_bigint< Q + 4 > res = ( res_min , mask);
	return res;
}


template <int Q>
sc_bigint< Q + 8 > Min_Mask_8_CA2 (sc_bigint< 8 * Q> Mag)
{
#pragma HLS INLINE
	sc_bigint< 4 * Q > ma  = Mag.range( (4 * Q) - 1, 0 );
	sc_bigint< 4 * Q > mb  = Mag.range( (8 * Q) - 1, 4 * Q );

	sc_biguint< 4 > is_min = VECTOR_IS_MIN< 4, Q >( mb, ma );
	sc_biguint< 4 > n_is_min = (~is_min);

	sc_biguint< 8 > mask_A = (is_min, n_is_min );

	sc_bigint< 4 * Q > min = VECTOR_MUX< 4 , Q >( mb, ma, is_min );

	sc_bigint< Q + 4 > min_mask = Min_Mask_4_CA2< Q >( min );
	sc_bigint< Q > res_min = min_mask.range(Q + 3 , 4);
	sc_biguint< 4 > i_mask = (sc_biguint< 4 >) min_mask.range(3,0);

	sc_biguint< 8 > mask_B = (i_mask, i_mask);

	sc_biguint< 8 > mask = VECTOR_AND< 8 >( mask_A, mask_B );

	sc_bigint< Q + 8 > res = ( res_min , mask);
	return res;
}

template <int Q>
sc_bigint< Q + 16 > Min_Mask_16_CA2 (sc_bigint< 16 * Q> Mag)
{
#pragma HLS INLINE
	sc_bigint< 8 * Q > ma  = Mag.range( (8 * Q) - 1, 0 );
	sc_bigint< 8 * Q > mb  = Mag.range( (16 * Q) - 1, 8 * Q );

	sc_biguint< 8 > is_min = VECTOR_IS_MIN< 8, Q >( mb, ma );
	sc_biguint< 8 > n_is_min = (~is_min);

	sc_biguint< 16 > mask_A = (is_min, n_is_min );

	sc_bigint< 8 * Q > min = VECTOR_MUX< 8 , Q >( mb, ma, is_min );

	sc_bigint< Q + 8 > min_mask = Min_Mask_8_CA2< Q >( min );
	sc_bigint< Q > res_min = min_mask.range(Q + 7 , 8);
	sc_biguint< 8 > i_mask = (sc_biguint< 8 >) min_mask.range(7,0);

	sc_biguint< 16 > mask_B = (i_mask, i_mask);

	sc_biguint< 16 > mask = VECTOR_AND< 16 >( mask_A, mask_B );

	sc_bigint< Q + 16 > res = ( res_min , mask);
	return res;
}

template <int Q>
sc_bigint< Q + 32 > Min_Mask_32_CA2 (sc_bigint< 32 * Q> Mag)
{
#pragma HLS INLINE
	sc_bigint< 16 * Q > ma  = Mag.range( (16 * Q) - 1, 0 );
	sc_bigint< 16 * Q > mb  = Mag.range( (32 * Q) - 1, 16 * Q );

	sc_biguint< 16 > is_min = VECTOR_IS_MIN< 16, Q >( mb, ma );
	sc_biguint< 16 > n_is_min = (~is_min);

	sc_biguint< 32 > mask_A = (is_min, n_is_min );

	sc_bigint< 16 * Q > min = VECTOR_MUX< 16 , Q >( mb, ma, is_min );

	sc_bigint< Q + 16 > min_mask = Min_Mask_16_CA2< Q >( min );
	sc_bigint< Q > res_min = min_mask.range(Q + 15 , 16);
	sc_biguint< 16 > i_mask = (sc_biguint< 16 >) min_mask.range(15,0);

	sc_biguint< 32 > mask_B = (i_mask, i_mask);

	sc_biguint< 32 > mask = VECTOR_AND< 32 >( mask_A, mask_B );

	sc_bigint< Q + 32 > res = ( res_min , mask);
	return res;
}

template <int Q>
sc_bigint< Q + 64 > Min_Mask_64_CA2 (sc_bigint< 64 * Q> Mag)
{
#pragma HLS INLINE
	sc_bigint< 32 * Q > ma  = Mag.range( (32 * Q) - 1, 0 );
	sc_bigint< 32 * Q > mb  = Mag.range( (64 * Q) - 1, 32 * Q );

	sc_biguint< 32 > is_min = VECTOR_IS_MIN< 32, Q >( mb, ma );
	sc_biguint< 32 > n_is_min = (~is_min);

	sc_biguint< 64 > mask_A = (is_min, n_is_min );

	sc_bigint< 32 * Q > min = VECTOR_MUX< 32 , Q >( mb, ma, is_min );

	sc_bigint< Q + 32 > min_mask = Min_Mask_32_CA2< Q >( min );
	sc_bigint< Q > res_min = min_mask.range(Q + 31 , 32);
	sc_biguint< 32 > i_mask = (sc_biguint< 32 >) min_mask.range(31,0);

	sc_biguint< 64 > mask_B = (i_mask, i_mask);

	sc_biguint< 64 > mask = VECTOR_AND< 64 >( mask_A, mask_B );

	sc_bigint< Q + 64 > res = ( res_min , mask);
	return res;
}

template <int Q>
sc_bigint< Q + 128 > Min_Mask_128_CA2 (sc_bigint< 128 * Q> Mag)
{
#pragma HLS INLINE
	sc_bigint< 64 * Q > ma  = Mag.range( (64 * Q) - 1, 0 );
	sc_bigint< 64 * Q > mb  = Mag.range( (128 * Q) - 1, 64 * Q );

	sc_biguint< 64 > is_min = VECTOR_IS_MIN< 64, Q >( mb, ma );
	sc_biguint< 64 > n_is_min = (~is_min);

	sc_biguint< 128 > mask_A = (is_min, n_is_min );

	sc_bigint< 64 * Q > min = VECTOR_MUX< 64 , Q >( mb, ma, is_min );

	sc_bigint< Q + 64 > min_mask = Min_Mask_64_CA2< Q >( min );
	sc_bigint< Q > res_min = min_mask.range(Q + 63 , 64);
	sc_biguint< 64 > i_mask = (sc_biguint< 64 >) min_mask.range(63,0);

	sc_biguint< 128 > mask_B = (i_mask, i_mask);

	sc_biguint< 128 > mask = VECTOR_AND< 128 >( mask_A, mask_B );

	sc_bigint< Q + 128 > res = ( res_min , mask);
	return res;
}

template <int Q>
sc_bigint< Q + 256 > Min_Mask_256_CA2 (sc_bigint< 256 * Q> Mag)
{
#pragma HLS INLINE
	sc_bigint< 128 * Q > ma  = Mag.range( (128 * Q) - 1, 0 );
	sc_bigint< 128 * Q > mb  = Mag.range( (256 * Q) - 1, 128 * Q );

	sc_biguint< 128 > is_min = VECTOR_IS_MIN< 128, Q >( mb, ma );
	sc_biguint< 128 > n_is_min = (~is_min);

	sc_biguint< 256 > mask_A = (is_min, n_is_min );

	sc_bigint< 128 * Q > min = VECTOR_MUX< 128 , Q >( mb, ma, is_min );

	sc_bigint< Q + 128 > min_mask = Min_Mask_128_CA2< Q >( min );
	sc_bigint< Q > res_min = min_mask.range(Q + 127 , 128);
	sc_biguint< 128 > i_mask = (sc_biguint< 128 >) min_mask.range(127,0);

	sc_biguint< 256 > mask_B = (i_mask, i_mask);

	sc_biguint< 256 > mask = VECTOR_AND< 256 >( mask_A, mask_B );

	sc_bigint< Q + 256 > res = ( res_min , mask);
	return res;
}

/////////////// SIGMAG ////////////////

template <int Q>
sc_biguint< Q + 2 > Min_Mask_2_SM (sc_biguint< 2 * Q> Mag)
{
#pragma HLS INLINE
	sc_biguint< Q > ma  = Mag.range( Q - 1, 0 );
	sc_biguint< Q > mb  = Mag.range( 2 * Q - 1, Q );

	sc_biguint< 1 > is_min = VECTOR_IS_MIN_SM< 1, Q >(mb, ma);
	sc_biguint< 1 > n_is_min = (~is_min);

	sc_biguint< Q > min = VECTOR_MUX_SM< 1 , Q >( mb, ma, is_min );

	sc_biguint< 2 > mask = (is_min, n_is_min);

	sc_biguint< Q + 2 > res = ( min , mask);
	return res;
}

template <int Q>
sc_biguint< Q + 4 > Min_Mask_4_SM (sc_biguint< 4 * Q> Mag)
{
#pragma HLS INLINE
	sc_biguint< 2 * Q > ma  = Mag.range( (2 * Q) - 1, 0 );
	sc_biguint< 2 * Q > mb  = Mag.range( (4 * Q) - 1, 2 * Q );

	sc_biguint< 2 > is_min = VECTOR_IS_MIN_SM< 2, Q >( mb, ma );
	sc_biguint< 2 > n_is_min = (~is_min);

	sc_biguint< 4 > mask_A = (is_min, n_is_min );

	sc_biguint< 2 * Q > min = VECTOR_MUX_SM< 2 , Q >( mb, ma, is_min );

	sc_biguint< Q + 2 > min_mask = Min_Mask_2_SM< Q >( min );
	sc_biguint< Q > res_min = min_mask.range(Q + 1 , 2);
	sc_biguint< 2 > i_mask = (sc_biguint< 2 >) min_mask.range(1,0);

	sc_biguint< 4 > mask_B = (i_mask, i_mask);

	sc_biguint< 4 > mask = VECTOR_AND< 4 >( mask_A, mask_B );

	sc_biguint< Q + 4 > res = ( res_min , mask);
	return res;
}


template <int Q>
sc_biguint< Q + 8 > Min_Mask_8_SM (sc_biguint< 8 * Q> Mag)
{
#pragma HLS INLINE
	sc_biguint< 4 * Q > ma  = Mag.range( (4 * Q) - 1, 0 );
	sc_biguint< 4 * Q > mb  = Mag.range( (8 * Q) - 1, 4 * Q );

	sc_biguint< 4 > is_min = VECTOR_IS_MIN_SM< 4, Q >( mb, ma );
	sc_biguint< 4 > n_is_min = (~is_min);

	sc_biguint< 8 > mask_A = (is_min, n_is_min );

	sc_biguint< 4 * Q > min = VECTOR_MUX_SM< 4 , Q >( mb, ma, is_min );

	sc_biguint< Q + 4 > min_mask = Min_Mask_4_SM< Q >( min );
	sc_biguint< Q > res_min = min_mask.range(Q + 3 , 4);
	sc_biguint< 4 > i_mask = (sc_biguint< 4 >) min_mask.range(3,0);

	sc_biguint< 8 > mask_B = (i_mask, i_mask);

	sc_biguint< 8 > mask = VECTOR_AND< 8 >( mask_A, mask_B );

	sc_biguint< Q + 8 > res = ( res_min , mask);
	return res;
}

template <int Q>
sc_biguint< Q + 16 > Min_Mask_16_SM (sc_biguint< 16 * Q> Mag)
{
#pragma HLS INLINE
	sc_biguint< 8 * Q > ma  = Mag.range( (8 * Q) - 1, 0 );
	sc_biguint< 8 * Q > mb  = Mag.range( (16 * Q) - 1, 8 * Q );

	sc_biguint< 8 > is_min = VECTOR_IS_MIN_SM< 8, Q >( mb, ma );
	sc_biguint< 8 > n_is_min = (~is_min);

	sc_biguint< 16 > mask_A = (is_min, n_is_min );

	sc_biguint< 8 * Q > min = VECTOR_MUX_SM< 8 , Q >( mb, ma, is_min );

	sc_biguint< Q + 8 > min_mask = Min_Mask_8_SM< Q >( min );
	sc_biguint< Q > res_min = min_mask.range(Q + 7 , 8);
	sc_biguint< 8 > i_mask = (sc_biguint< 8 >) min_mask.range(7,0);

	sc_biguint< 16 > mask_B = (i_mask, i_mask);

	sc_biguint< 16 > mask = VECTOR_AND< 16 >( mask_A, mask_B );

	sc_biguint< Q + 16 > res = ( res_min , mask);
	return res;
}

template <int Q>
sc_biguint< Q + 32 > Min_Mask_32_SM (sc_biguint< 32 * Q> Mag)
{
#pragma HLS INLINE
	sc_biguint< 16 * Q > ma  = Mag.range( (16 * Q) - 1, 0 );
	sc_biguint< 16 * Q > mb  = Mag.range( (32 * Q) - 1, 16 * Q );

	sc_biguint< 16 > is_min = VECTOR_IS_MIN_SM< 16, Q >( mb, ma );
	sc_biguint< 16 > n_is_min = (~is_min);

	sc_biguint< 32 > mask_A = (is_min, n_is_min );

	sc_biguint< 16 * Q > min = VECTOR_MUX_SM< 16 , Q >( mb, ma, is_min );

	sc_biguint< Q + 16 > min_mask = Min_Mask_16_SM< Q >( min );
	sc_biguint< Q > res_min = min_mask.range(Q + 15 , 16);
	sc_biguint< 16 > i_mask = (sc_biguint< 16 >) min_mask.range(15,0);

	sc_biguint< 32 > mask_B = (i_mask, i_mask);

	sc_biguint< 32 > mask = VECTOR_AND< 32 >( mask_A, mask_B );

	sc_biguint< Q + 32 > res = ( res_min , mask);
	return res;
}

template <int Q>
sc_biguint< Q + 64 > Min_Mask_64_SM (sc_biguint< 64 * Q> Mag)
{
#pragma HLS INLINE
	sc_biguint< 32 * Q > ma  = Mag.range( (32 * Q) - 1, 0 );
	sc_biguint< 32 * Q > mb  = Mag.range( (64 * Q) - 1, 32 * Q );

	sc_biguint< 32 > is_min = VECTOR_IS_MIN_SM< 32, Q >( mb, ma );
	sc_biguint< 32 > n_is_min = (~is_min);

	sc_biguint< 64 > mask_A = (is_min, n_is_min );

	sc_biguint< 32 * Q > min = VECTOR_MUX_SM< 32 , Q >( mb, ma, is_min );

	sc_biguint< Q + 32 > min_mask = Min_Mask_32_SM< Q >( min );
	sc_biguint< Q > res_min = min_mask.range(Q + 31 , 32);
	sc_biguint< 32 > i_mask = (sc_biguint< 32 >) min_mask.range(31,0);

	sc_biguint< 64 > mask_B = (i_mask, i_mask);

	sc_biguint< 64 > mask = VECTOR_AND< 64 >( mask_A, mask_B );

	sc_biguint< Q + 64 > res = ( res_min , mask);
	return res;
}

template <int Q>
sc_biguint< Q + 128 > Min_Mask_128_SM (sc_biguint< 128 * Q> Mag)
{
#pragma HLS INLINE
	sc_biguint< 64 * Q > ma  = Mag.range( (64 * Q) - 1, 0 );
	sc_biguint< 64 * Q > mb  = Mag.range( (128 * Q) - 1, 64 * Q );

	sc_biguint< 64 > is_min = VECTOR_IS_MIN_SM< 64, Q >( mb, ma );
	sc_biguint< 64 > n_is_min = (~is_min);

	sc_biguint< 128 > mask_A = (is_min, n_is_min );

	sc_biguint< 64 * Q > min = VECTOR_MUX_SM< 64 , Q >( mb, ma, is_min );

	sc_biguint< Q + 64 > min_mask = Min_Mask_64_SM< Q >( min );
	sc_biguint< Q > res_min = min_mask.range(Q + 63 , 64);
	sc_biguint< 64 > i_mask = (sc_biguint< 64 >) min_mask.range(63,0);

	sc_biguint< 128 > mask_B = (i_mask, i_mask);

	sc_biguint< 128 > mask = VECTOR_AND< 128 >( mask_A, mask_B );

	sc_biguint< Q + 128 > res = ( res_min , mask);
	return res;
}

template <int Q>
sc_biguint< Q + 256 > Min_Mask_256_SM (sc_biguint< 256 * Q> Mag)
{
#pragma HLS INLINE
	sc_biguint< 128 * Q > ma  = Mag.range( (128 * Q) - 1, 0 );
	sc_biguint< 128 * Q > mb  = Mag.range( (256 * Q) - 1, 128 * Q );

	sc_biguint< 128 > is_min = VECTOR_IS_MIN_SM< 128, Q >( mb, ma );
	sc_biguint< 128 > n_is_min = (~is_min);

	sc_biguint< 256 > mask_A = (is_min, n_is_min );

	sc_biguint< 128 * Q > min = VECTOR_MUX_SM< 128 , Q >( mb, ma, is_min );

	sc_biguint< Q + 128 > min_mask = Min_Mask_128_SM< Q >( min );
	sc_biguint< Q > res_min = min_mask.range(Q + 127 , 128);
	sc_biguint< 128 > i_mask = (sc_biguint< 128 >) min_mask.range(127,0);

	sc_biguint< 256 > mask_B = (i_mask, i_mask);

	sc_biguint< 256 > mask = VECTOR_AND< 256 >( mask_A, mask_B );

	sc_biguint< Q + 256 > res = ( res_min , mask);
	return res;
}

/////////////// Global Min Mask ////////////////

template <int Q>
sc_bigint< Q + 2 > Min_Mask_TREE_2 (sc_bigint< 2 * Q> llr)
{
#pragma HLS INLINE

#if defined CA2
	sc_bigint< 2 * Q > abs = VECTOR_ABS< 2, Q >(llr);
	return Min_Mask_2_CA2< Q >( abs );
#elif defined SIGMAG
	sc_biguint< 2 * (Q-1) > abs = VECTOR_ABS_SM< 2, Q >(llr);
	sc_biguint< Q + 2 > res = (sc_biguint< Q + 2 >) Min_Mask_2_SM< Q-1 >( abs );
	return (sc_bigint< Q + 2 >) res;
#endif
}

template <int Q>
sc_bigint< Q + 4 > Min_Mask_TREE_4 (sc_bigint< 4 * Q> llr)
{
#pragma HLS INLINE

#if defined CA2
	sc_bigint< 4 * Q > abs = VECTOR_ABS< 4, Q >(llr);
	return Min_Mask_4_CA2< Q >( abs );
#elif defined SIGMAG
	sc_biguint< 4 * (Q-1) > abs = VECTOR_ABS_SM< 4, Q >(llr);
	sc_biguint< Q + 4 > res = (sc_biguint< Q + 4 >) Min_Mask_4_SM< Q-1 >( abs );
	return (sc_bigint< Q + 4 >) res;
#endif
}

template <int Q>
sc_bigint< Q + 8 > Min_Mask_TREE_8 (sc_bigint< 8 * Q> llr)
{
#pragma HLS INLINE

#if defined CA2
	sc_bigint< 8 * Q > abs = VECTOR_ABS< 8, Q >(llr);
	return Min_Mask_8_CA2< Q >( abs );
#elif defined SIGMAG
	sc_biguint< 8 * (Q-1) > abs = VECTOR_ABS_SM< 8, Q >(llr);
	sc_biguint< Q + 8 > res = (sc_biguint< Q + 8 >) Min_Mask_8_SM< Q-1 >( abs );
	return (sc_bigint< Q + 8 >) res;
#endif
}

template <int Q>
sc_bigint< Q + 16 > Min_Mask_TREE_16 (sc_bigint< 16 * Q> llr)
{
#pragma HLS INLINE

#if defined CA2
	sc_bigint< 16 * Q > abs = VECTOR_ABS< 16, Q >(llr);
	return Min_Mask_16_CA2< Q >( abs );
#elif defined SIGMAG
	sc_biguint< 16 * (Q-1) > abs = VECTOR_ABS_SM< 16, Q >(llr);
	sc_biguint< Q + 16 > res = (sc_biguint< Q + 16 >) Min_Mask_16_SM< Q-1 >( abs );
	return (sc_bigint< Q + 16 >) res;
#endif
}

template <int Q>
sc_bigint< Q + 32 > Min_Mask_TREE_32 (sc_bigint< 32 * Q> llr)
{
#pragma HLS INLINE

#if defined CA2
	sc_bigint< 32 * Q > abs = VECTOR_ABS< 32, Q >(llr);
	return Min_Mask_32_CA2< Q >( abs );
#elif defined SIGMAG
	sc_biguint< 32 * (Q-1) > abs = VECTOR_ABS_SM< 32, Q >(llr);
	sc_biguint< Q + 32 > res = (sc_biguint< Q + 32 >) Min_Mask_32_SM< Q-1 >( abs );
	return (sc_bigint< Q + 32 >) res;
#endif
}

template <int Q>
sc_bigint< Q + 64 > Min_Mask_TREE_64 (sc_bigint< 64 * Q> llr)
{
#pragma HLS INLINE

#if defined CA2
	sc_bigint< 64 * Q > abs = VECTOR_ABS< 64, Q >(llr);
	return Min_Mask_64_CA2< Q >( abs );
#elif defined SIGMAG
	sc_biguint< 64 * (Q-1) > abs = VECTOR_ABS_SM< 64, Q >(llr);
	sc_biguint< Q + 64 > res = (sc_biguint< Q + 64 >) Min_Mask_64_SM< Q-1 >( abs );
	return (sc_bigint< Q + 64 >) res;
#endif
}

template <int Q>
sc_bigint< Q + 128 > Min_Mask_TREE_128 (sc_bigint< 128 * Q> llr)
{
#pragma HLS INLINE

#if defined CA2
	sc_bigint< 128 * Q > abs = VECTOR_ABS< 128, Q >(llr);
	return Min_Mask_128_CA2< Q >( abs );
#elif defined SIGMAG
	sc_biguint< 128 * (Q-1) > abs = VECTOR_ABS_SM< 128, Q >(llr);
	sc_biguint< Q + 128 > res = (sc_biguint< Q + 128 >) Min_Mask_128_SM< Q-1 >( abs );
	return (sc_bigint< Q + 128 >) res;
#endif
}

template <int Q>
sc_bigint< Q + 256 > Min_Mask_TREE_256 (sc_bigint< 256 * Q> llr)
{
#pragma HLS INLINE

#if defined CA2
	sc_bigint< 256 * Q > abs = VECTOR_ABS< 256, Q >(llr);
	return Min_Mask_256_CA2< Q >( abs );
#elif defined SIGMAG
	sc_biguint< 256 * (Q-1) > abs = VECTOR_ABS_SM< 256, Q >(llr);
	sc_biguint< Q + 256 > res = (sc_biguint< Q + 256 >) Min_Mask_256_SM< Q-1 >( abs );
	return (sc_bigint< Q + 256 >) res;
#endif
}

//*************************************************************************//
// 							Recursive Functions 						   //
//*************************************************************************//

template <int P, int Q>
sc_biguint<P> Spec_P_gen ( sc_bigint<P*Q> llr, sc_biguint<P> fb )
{
#pragma HLS INLINE

	if( P == 2)
		return Spec_P2<Q> (llr,fb);
	else
	{
		sc_bigint<(P/2)*Q> la = llr.range( ((P/2) * Q) - 1, 0 );
		sc_bigint<(P/2)*Q> lb = llr.range( 2* ((P/2) * Q) - 1, ((P/2) * Q) );

		// F
		sc_bigint<(P/2)*Q> la1 = Function_F<(P/2),Q> (la,lb);

		// First P
		sc_biguint<(P/2)> sa1 = Spec_P_gen<(P/2),Q> ( la1, fb.range(1,0) );

		// G
		sc_bigint<(P/2)*Q> lb1 = Function_G<(P/2),Q> (la,lb,sa1);

		// Second P
		sc_biguint<(P/2)> sb1 = Spec_P_gen<(P/2),Q> ( lb1 , fb.range(3,2) );

		// XOR
		sc_biguint<(P/2)> sa = VECTOR_XOR<(P/2)>(sa1,sb1);
		sc_biguint<(P/2)> sb = sb1;

		return ( sb, sa );
	}
}

template <int P, int Q>
sc_biguint<P> Spec_P_ext_gen ( sc_bigint<P*Q> llr, sc_biguint<P> fb )
{
#pragma HLS INLINE

	if( P == 2)
		return Spec_P2<Q> (llr,fb);
	else
	{
		sc_bigint<(P/2)*Q> la = llr.range( ((P/2) * Q) - 1, 0 );
		sc_bigint<(P/2)*Q> lb = llr.range( 2* ((P/2) * Q) - 1, ((P/2) * Q) );

		// F
		sc_bigint<(P/2)*Q> la1 = Function_F<(P/2),Q> (la,lb);

		// First P
		sc_biguint<(P/2)> sa1 = Spec_P_ext_gen<(P/2),Q> ( la1, fb.range(1,0) );

		// G
		sc_bigint<(P/2)*(Q+1)> lb1 = Function_G_ext<(P/2),Q> (la,lb,sa1);

		// Second P
		sc_biguint<(P/2)> sb1 = Spec_P_ext_gen<(P/2),(Q+1)> ( lb1 , fb.range(3,2) );

		// XOR
		sc_biguint<(P/2)> sa = VECTOR_XOR<(P/2)>(sa1,sb1);
		sc_biguint<(P/2)> sb = sb1;

		return ( sb, sa );
	}
}

#endif
