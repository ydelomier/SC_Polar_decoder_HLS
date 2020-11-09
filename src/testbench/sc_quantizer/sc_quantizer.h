/*
 *  Adder.h
 *  SystemC_SimpleAdder
 *
 *  Created by Le Gal on 07/05/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef _sc_quantizer_
#define _sc_quantizer_

#include "systemc.h"

#ifdef __RTL_SIMULATION__
#include "config.h"
#else
#include "../../module/config.h"
#endif

SC_MODULE(sc_quantizer)          // module (class) declaration
{
    sc_in <bool> clk;
    sc_in <bool> reset;

    sc_in < sc_uint<8> > 	BETA;
    sc_in < sc_int<8> > 	VSATN;
    sc_in < sc_int<8> > 	VSATP;

	sc_fifo_in < float > e;
	sc_fifo_out< LLR   > s;

	SC_CTOR(sc_quantizer)
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
void sc_quantizer::do_gen()           // process
{
	/*
    const float FACTEUR_BETA    = 4.0f;
    const signed short vSAT_NEG = -31;
    const signed short vSAT_POS = +31;
//    const float FACTEUR_BETA    = 64.0f;
//    const signed short vSAT_NEG = -512;
//    const signed short vSAT_POS = +512;


#ifdef _SIMULATION_
    cout << "(II) sc_quantizer :: FACTEUR  = " << FACTEUR_BETA << endl;
    cout << "(II) sc_quantizer :: vSAT_NEG = " << (vSAT_NEG)   << endl;
    cout << "(II) sc_quantizer :: vSAT_POS = " << (vSAT_POS)   << endl;
#endif
*/
    //unsigned int counter = 0;

    sc_uint<8> FACTEUR_BETA   = BETA.read();
    signed short vSAT_NEG 	= (signed short)VSATN.read();
    signed short vSAT_POS 	= (signed short)VSATP.read();
    while( true )
    {
		#pragma HLS PIPELINE

        FACTEUR_BETA 	= BETA.read();
        vSAT_NEG 	= (signed short)VSATN.read();
        vSAT_POS 	= (signed short)VSATP.read();
        const float        fvalue = e.read();
        const signed short ivalue = (signed short)(fvalue * FACTEUR_BETA);
        const signed short mvalue = (ivalue > vSAT_NEG) ? ivalue : vSAT_NEG;
        const signed short rvalue = (mvalue < vSAT_POS) ? mvalue : vSAT_POS;
        s.write( rvalue );
    }
}
//
//
////////////////////////////////////////////////////////////////////////////////
//
//
#endif
