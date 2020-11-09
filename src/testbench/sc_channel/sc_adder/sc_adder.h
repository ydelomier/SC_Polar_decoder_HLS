/*
 *  Adder.h
 *  SystemC_SimpleAdder
 *
 *  Created by Le Gal on 07/05/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef _sc_adder_
#define _sc_adder_

#include "systemc.h"

#include <iostream>
using namespace std;
#include <iomanip>

SC_MODULE(sc_adder)
{
  sc_in <bool> clk;
  sc_in <bool> reset;

  //sc_in < sc_uint<5> > snr;
  sc_in < float > sigma;

  sc_fifo_in < float > e1;
  sc_fifo_in < float > e2;
  sc_fifo_out< float > s;

  SC_CTOR(sc_adder)
  {
    SC_CTHREAD(do_gen, clk.pos());
    reset_signal_is(reset, true);
  }

private:
  void do_gen();

};
//
//
////////////////////////////////////////////////////////////////////////////////
//
//
void sc_adder::do_gen()           // process
{
	wait();

/*	const float tab_sigma[ ] =
	{
		1.000000f, //SNR = 0.000000, rate = 0.500000, b = X 000 00
		0.971628f, //SNR = 0.250000, rate = 0.500000, b = X 000 01
		0.944061f, //SNR = 0.500000, rate = 0.500000, b = X 000 10
		0.917276f, //SNR = 0.750000, rate = 0.500000, b = X 000 11
		0.891251f, //SNR = 1.000000, rate = 0.500000, b = X 001 00
		0.865964f, //SNR = 1.250000, rate = 0.500000, b = X 001 01
		0.841395f, //SNR = 1.500000, rate = 0.500000, b = X 001 10
		0.817523f, //SNR = 1.750000, rate = 0.500000, b = X 001 11
		0.794328f, //SNR = 2.000000, rate = 0.500000, b = X 010 00
		0.771792f, //SNR = 2.250000, rate = 0.500000, b = X 010 01
		0.749894f, //SNR = 2.500000, rate = 0.500000, b = X 010 10
		0.728618f, //SNR = 2.750000, rate = 0.500000, b = X 010 11
		0.707946f, //SNR = 3.000000, rate = 0.500000, b = X 011 00
		0.687860f, //SNR = 3.250000, rate = 0.500000, b = X 011 01
		0.668344f, //SNR = 3.500000, rate = 0.500000, b = X 011 10
		0.649382f, //SNR = 3.750000, rate = 0.500000, b = X 011 11
		0.630957f, //SNR = 4.000000, rate = 0.500000, b = X 100 00
		0.613056f, //SNR = 4.250000, rate = 0.500000, b = X 100 01
		0.595662f, //SNR = 4.500000, rate = 0.500000, b = X 100 10
		0.578762f, //SNR = 4.750000, rate = 0.500000, b = X 100 11
		0.562341f, //SNR = 5.000000, rate = 0.500000, b = X 101 00
		0.546387f, //SNR = 5.250000, rate = 0.500000, b = X 101 01
		0.530884f, //SNR = 5.500000, rate = 0.500000, b = X 101 10
		0.515822f, //SNR = 5.750000, rate = 0.500000, b = X 101 11
		0.501187f, //SNR = 6.000000, rate = 0.500000, b = X 110 00
		0.486968f, //SNR = 6.250000, rate = 0.500000, b = X 110 01
		0.473151f, //SNR = 6.500000, rate = 0.500000, b = X 110 10
		0.459727f, //SNR = 6.750000, rate = 0.500000, b = X 110 11
		0.446684f, //SNR = 7.000000, rate = 0.500000, b = X 111 00
		0.434010f, //SNR = 7.250000, rate = 0.500000, b = X 111 01
		0.421697f, //SNR = 7.500000, rate = 0.500000, b = X 111 10
		0.000000f  //SNR = 7.750000, rate = 0.500000, b = X 111 11
	};*/

	//const sc_uint<5> SNR    = snr.read();

	//float const sigma = tab_sigma[ SNR ];


/*#if 0
	#define NB_BITS  (40*1024*1024)
    int* t_bits = new int[NB_BITS];

    for(int i=0; i<NB_BITS; i++)
    	t_bits[i] = i%2;

    printf("SigB = %f ", sigma); printf("\n");

    for(int i=0; i<16; i++) printf("%d ", t_bits[i]); printf("\n");

    float* t_float = new float[NB_BITS];
    for(int i=0; i<NB_BITS; i+=2){
        t_float[i+0] = (t_bits[i+0] == 1 ? 1.0f : -1.0f) + (sigma * e2.read());
        t_float[i+1] = (t_bits[i+1] == 1 ? 1.0f : -1.0f) + (sigma * e2.read());
    }

    for(int i=0; i<16; i++) printf("%f ", t_float[i]); printf("\n");

    int* t_dec = new int[NB_BITS];
    for(int i=0; i<NB_BITS; i++)
    	t_dec[i] = ((int)(t_float[i] * 8.0f)) > 0;

    for(int i=0; i<16; i++) printf("%d ", t_dec[i]); printf("\n");

    int sumE = 0;
    for(int i=0; i<NB_BITS; i++)
    	sumE += (t_dec[i] != t_bits[i]);

    printf("sumE = %d\n", sumE);
    printf("BER = %f\n", ((float)sumE)/((float)NB_BITS));
    exit( 0 );
#endif*/


  int eBits = 0;
  int pBits = 1;
  float sum  = 0.0f;
  float sumA = 0.0f;
	while( true )
	{
		#pragma HLS PIPELINE

		//const sc_uint<5> SNR    = snr.read();
		//float const sigma = tab_sigma[ SNR ];
		const float sig = sigma.read();

		float oBit  = e1.read();					// la valeur du bit module
		float noise = e2.read();
		float n     = noise * sig;  // la valeur du bruit mulitplie par sigma
		float nBit  = oBit + n;

/*
    sumA  += fabs(n);
    sum   += n;
    eBits += (signbit(oBit) != signbit(nBit) );
    printf("%1.3f ", n);
    if( ((pBits++) % 649) == 0 ){
      printf("\n%d %d => %f\n", pBits, eBits, (float)eBits/(float)pBits);
      printf("sumA = %f\n", sumA / pBits);
      printf("sum  = %f\n", sum  / pBits);
      exit( 0 );
    }
*/
		s.write( nBit );							// la valeur du llr post bruitage
	}
}
//
//
////////////////////////////////////////////////////////////////////////////////
//
//
#endif
