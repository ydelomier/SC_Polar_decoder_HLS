/*
 *  sc_awgn.h
 *  sc_awgn
 *
 *  Created by Bertrand LE GAL on 06/11/2016.
 *  Copyright 2016 Bertrand LE GAL. All rights reserved.
 *
 */
#ifndef _sc_awgn_half_
#define _sc_awgn_half_

#include "systemc.h"
#include <hls_half.h>

SC_MODULE(sc_awgn_half)
{
  sc_in <bool> clk;
  sc_in <bool> reset;

  sc_fifo_in < float > e1;
  sc_fifo_in < float > e2;

  sc_fifo_out< float > s;

  SC_CTOR(sc_awgn_half)
  {
    SC_CTHREAD(do_action, clk.pos());
    reset_signal_is(reset,true);
  }

private:
  void do_action();
};

#endif
