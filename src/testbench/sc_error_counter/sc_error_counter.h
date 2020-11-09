/*
 *  Adder.h
 *  SystemC_SimpleAdder
 *
 *  Created by Le Gal on 07/05/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef _sc_error_counter_
#define _sc_error_counter_

#include "systemc.h"

#ifdef __RTL_SIMULATION__
#include "config.h"
#include "polar_parameters.h"
#else
#include "../../module/config.h"
#include "../../module/polar_parameters.h"
#endif


#include <iostream>

SC_MODULE(sc_error_counter)
{
    sc_in <bool> clk;
    sc_in <bool> reset;

    sc_fifo_in< BIT > dec;
    sc_fifo_in< BIT > ref;

    sc_out    < BIT > dec_frame;
    sc_out    < BIT > err_frame;

    sc_out    < sc_uint<64> > err_bits;
    sc_out    < sc_uint<64> > err_frames;

    sc_out    < sc_uint<64> > proc_bits;
    sc_out    < sc_uint<64> > proc_frames;

    SC_CTOR(sc_error_counter)
    {
        SC_CTHREAD(do_gen, clk.pos());
        reset_signal_is(reset, true);
    }

private:

    void do_gen()
    {
    	dec_frame.write  ( 0 );
    	err_frame.write  ( 0 );
    	err_bits.write   ( 0 );
    	err_frames.write ( 0 );
    	proc_bits.write  ( 0 );
    	proc_frames.write( 0 );

    	sc_uint<64> errBE   = 0;
    	sc_uint<64> errFE   = 0;
    	sc_uint<64> pBits   = 0;
    	sc_uint<64> pFrames = 0;

    	wait();

    	int cpti = 0;

    	while( true )
    	{
    		// TODO: ADAPTER LA VALEUR 10 AU CODE !
    		sc_uint<10> err = 0;

    		BIT bit1 [_NBITS];
    		BIT bit2 [_NBITS];

    		for(int i = 0; i < _NBITS;  i += 1)
    		{
    			#pragma HLS PIPELINE
    			BIT bit_1 = dec.read();
    			bit1[i] = bit_1;
    			BIT bit_2 = ref.read();
    			bit2[i] = bit_2;
    			err += (bit_1 != bit_2);
    		}

    		////////////////////////////////////
#ifdef DBG_ERR
    		cpti = cpti + 1;
    		cout << "	[ERR_CNT] "; printf("%5d ", cpti);
    		cout << ", BIT rf = ";
    		for( int i = 0; i < _NBITS; i++)
				printf("%3d ", bit1[i].to_int());
    		cout << endl;
    		cout << "	[ERR_CNT] "; printf("%5d ", cpti);
    		cout << ", BIT dc = ";
    		for( int i = 0; i < _NBITS; i++)
				printf("%3d ", bit2[i].to_int());
    		cout << endl;
#endif
    		////////////////////////////////////

    		//
    		// ON MET A JOUR LES STATISTIQUES
    		//
    		errBE   += err;
    		errFE   += err.or_reduce();
    		pBits   += _NBITS ;
    		pFrames += 1;

    		//
    		// ON RESCALE LES BEs CAR CELA MONTE TRES VITE !
    		//
    		err_bits.write   ( errBE   );
    		err_frames.write ( errFE   );
    		proc_bits.write  ( pBits   );
    		proc_frames.write( pFrames );

    		err_frame.write  ( err.or_reduce() );
    		dec_frame.write  ( 1 );
    		wait( );

    		err_frame.write  ( 0 );
    		dec_frame.write  ( 0 );
    		wait( );
    	}
    }

};

//
//
////////////////////////////////////////////////////////////////////////////////
//
//
#endif
