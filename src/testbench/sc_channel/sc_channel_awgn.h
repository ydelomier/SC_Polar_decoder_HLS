/*
 *  Adder.h
 *  SystemC_SimpleAdder
 *
 *  Created by Le Gal on 07/05/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef _sc_top_module_
#define _sc_top_module_

#include "../sc_modulate/sc_bpsk.h"
#include "../sc_channel/sc_xorshift128/sc_xorshift128.h"
#include "../sc_channel/sc_awgn/sc_awgn.h"
#include "../sc_channel/sc_adder/sc_adder.h"
#include "../sc_quantizer/sc_quantizer.h"

SC_MODULE(sc_channel_awgn)
{
    sc_in <bool> clk;
    sc_in <bool> reset;

    sc_in < sc_uint<5> > snr;

    sc_fifo_in < BIT > ivalue;
    sc_fifo_out< LLR > qvalue;

    SC_CTOR(sc_channel_awgn) :
        bpsk   ("bpsk"),
        randc ("rand1"),
        awgn  ("awgn" ),
        adder ("adder"),
        quant ("quant"),
    f2 ("bpsk_2_adder" ,    32),  // TYPE BIT
    f3 ("rand1_2_awgn" ,    32),  // TYPE FLOAT
    f4 ("rand2_2_awgn" ,    32),  // TYPE FLOAT
    f5 ("awgn_2_adder" ,    32),  // TYPE FLOAT
    f6 ("adder_2_quant",    32)
    {
        bpsk.clk  ( clk    );
        bpsk.reset( reset  );
        bpsk.e    ( ivalue );
        bpsk.s    ( f2     );

        randc.clk  ( clk   );
        randc.reset( reset );
        randc.s1   ( f3    );
        randc.s2   ( f4    );

        awgn.clk   ( clk   );
        awgn.reset ( reset );
        awgn.e1    ( f3    );
        awgn.e2    ( f4    );
        awgn.s     ( f5    );

        adder.clk   ( clk   );
        adder.reset ( reset );
        adder.snr   ( snr   );
        adder.e1    ( f2    ); // enc. bits
        adder.e2    ( f5    ); // noise
        adder.s     ( f6    );

        quant.clk   ( clk   );
        quant.reset ( reset );
        quant.e     ( f6    );
        quant.s     ( qvalue);
    }

    sc_bpsk            bpsk;
    sc_xorshift128     randc;
    sc_awgn            awgn;
    sc_adder           adder;
    sc_quantizer       quant;

    sc_fifo< float     > f2;
    sc_fifo< float     > f3;
    sc_fifo< float     > f4;
    sc_fifo< float     > f5;
    sc_fifo< float     > f6;

};

#endif
