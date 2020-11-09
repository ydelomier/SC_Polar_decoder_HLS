#ifndef _sc_top_module_
#define _sc_top_module_

#ifdef __RTL_SIMULATION__
  #include "my_module_rtl_wrapper.h"
  #define sc_fec my_module_rtl_wrapper
#else
  #include "../module/my_module.h"
  #define sc_fec my_module
#endif


#ifdef __RTL_SIMULATION__
	#include "sc_encoder.h"
	#include "sc_bpsk.h"
	#include "sc_xorshift128.h"
	#include "sc_awgn.h"
	#include "sc_adder.h"
	#include "sc_quantizer.h"
	#include "wrapper_in.h"
	#include "wrapper_out.h"
	#include "sc_error_counter.h"
	#include "sc_mbps_compute.h"
	#include "sc_xer_compute.h"
#else
	#include "sc_encoder/sc_encoder.h"
	#include "sc_modulate/sc_bpsk.h"
	#include "sc_channel/sc_xorshift128/sc_xorshift128.h"
	#include "sc_channel/sc_awgn/sc_awgn.h"
	#include "sc_channel/sc_adder/sc_adder.h"
	#include "sc_quantizer/sc_quantizer.h"
	#include "../module/wrapper_in.h"
	#include "../module/wrapper_out.h"
	#include "sc_error_counter/sc_error_counter.h"
	#include "sc_mbps_compute/sc_mbps_compute.h"
	#include "sc_xer_compute/sc_xer_compute.h"
#endif

SC_MODULE(sc_top_module)
{
    sc_in <bool> clk;
    sc_in <bool> reset;

    sc_in <bool> enable;
    sc_in < float > 		SIGMA;
    sc_in < sc_uint<8> > 	BETA;
    sc_in < sc_int<8> > 	VSATN;
    sc_in < sc_int<8> > 	VSATP;
    sc_in < sc_uint<8> > 	NBIT;
    //sc_out< sc_uint<16> > 	CPT;

    sc_in < sc_uint<8> > seed;

    sc_out< sc_uint<64> > err_bits;
    sc_out< sc_uint<64> > err_frames;
    sc_out< sc_uint<64> > proc_bits;
    sc_out< sc_uint<64> > proc_frames;

    sc_out< sc_uint<32> > ber_value;
    sc_out< sc_uint<32> > fer_value;
    sc_out< sc_uint<32> > mbps_value;


   SC_CTOR(sc_top_module) :
	polar  ("polar" ),
	encoder("enc" ),
	bpsk   ("bpsk"),
	randc ("rand1"),
	awgn  ("awgn" ),
	adder ("adder"),
	quant ("quant"),
	wrap_in ("w_in"),
	wrap_out ("w_o"),
	err_counter("err_counter"),
	mbps_compute("mbps_compute"),
	ber_counter("ber_counter"),
	fer_counter("fer_counter"),
	f1 ("enc_2_bpsk"   , 1),  // TYPE BIT
	f2 ("bpsk_2_adder" ,    1),  // TYPE BIT
	f3 ("rand1_2_awgn" ,    1),  // TYPE FLOAT
	f4 ("rand2_2_awgn" ,    1),  // TYPE FLOAT
	f5 ("awgn_2_adder" ,    1),  // TYPE FLOAT
	f6 ("adder_2_quant",    1),  // TYPE FLOAT
	f7 ("quant_2_wrap",  1),  // TYPE LLR
	f8 ("wrap_2_dec",  1),  // TYPE LLR
	f9 ("dec_2_wrap",  1),  // TYPE BITS
	f10("wrap_2_err",    1),  // TYPE BIT
	f11("enc_2_err"    , 32768),  // TYPE BIT
	FB ("fb", 1),
	dec_frame    ("dec_frame"),
	err_frame    ("err_frame"),
	_proc_bits   ("_proc_bits"),
	_proc_frames ("_dec_frame"),
	_err_bits    ("_err_bits"),
	_err_frames  ("_err_frame"),
	__proc_bits  ("__proc_bits"),
	__proc_frames("__dec_frame"),
	__err_bits   ("__err_bits"),
	__err_frames ("__err_frame")
    {

        encoder.clk   ( clk      );
        encoder.reset ( reset    );
        encoder.enable( enable   );
        encoder.s1    ( f1       );
        encoder.s2    ( f11       );
        encoder.FB (FB);

        bpsk.clk  ( clk   );
        bpsk.reset( reset );
        bpsk.e    ( f1    );
        bpsk.s    ( f2    );

        randc.clk  ( clk   );
        randc.reset( reset );
        randc.seed ( seed  );
        randc.s1   ( f3    );
        randc.s2   ( f4    );

        awgn.clk   ( clk   );
        awgn.reset ( reset );
        awgn.e1    ( f3    );
        awgn.e2    ( f4    );
        awgn.s     ( f5    );

        adder.clk   ( clk   );
        adder.reset ( reset );
        adder.sigma ( SIGMA   );
        adder.e1    ( f2    ); // enc. bits
        adder.e2    ( f5    ); // noise
        adder.s     ( f6    );

        quant.clk   ( clk   );
        quant.reset ( reset );
        quant.e     ( f6    );
        quant.s     ( f7    );
        quant.BETA     ( BETA    );
        quant.VSATN     ( VSATN    );
        quant.VSATP     ( VSATP    );

        wrap_in.clk(clk);
        wrap_in.reset(reset);
        wrap_in.e(f7);
        wrap_in.s(f8);

        polar.clk   ( clk      );
        polar.reset ( reset    );
        polar.FB    ( FB       );
        polar.e     ( f8       );
        polar.s     ( f9       );

        wrap_out.clk(clk);
        wrap_out.reset(reset);
        wrap_out.e(f9);
        wrap_out.s(f10);

        err_counter.clk   ( clk   );
        err_counter.reset ( reset );
        err_counter.dec   ( f10   );
        err_counter.ref   ( f11   );
        err_counter.dec_frame  ( dec_frame    );
        err_counter.err_frame  ( err_frame    );
        err_counter.err_bits   ( _err_bits    );
        err_counter.err_frames ( _err_frames  );
        err_counter.proc_bits  ( _proc_bits   );
        err_counter.proc_frames( _proc_frames );

        mbps_compute.clk       ( clk   );
        mbps_compute.reset     ( reset );
        mbps_compute.dec_frame ( dec_frame  );
        mbps_compute.mbps      ( mbps_value );

        ber_counter.clk  ( clk         );
        ber_counter.reset( reset       );
        ber_counter.decs ( __proc_bits );
        ber_counter.errs ( __err_bits  );
        ber_counter.xer  ( ber_value   );

        fer_counter.clk  ( clk           );
        fer_counter.reset( reset         );
        fer_counter.decs ( __proc_frames );
        fer_counter.errs ( __err_frames  );
        fer_counter.xer  ( fer_value     );

        SC_CTHREAD(do_gen, clk.pos());
        reset_signal_is(reset, true);

        SC_METHOD(gen_signals);
        sensitive << seed;
    }

    void gen_signals()
    {
        __proc_bits   = _proc_bits.read().range  (31, 0);
        __proc_frames = _proc_frames.read().range(31, 0);
        __err_bits    = _err_bits.read().range   (31, 0);
        __err_frames  = _err_frames.read().range (31, 0);
    }

    void do_gen()
    {
        proc_bits.write  ( 0 );
        proc_frames.write( 0 );
        err_bits.write   ( 0 );
        err_frames.write ( 0 );
        wait();

        while( true )
        {
            proc_bits.write  (_proc_bits  );
            proc_frames.write(_proc_frames);
            err_bits.write   (_err_bits   );
            err_frames.write (_err_frames );
            wait();
        }
    }

private:

    sc_fec             polar;
    sc_encoder         encoder;
    sc_bpsk            bpsk;
    sc_xorshift128     randc;
    sc_awgn            awgn;
    sc_adder           adder;
    sc_quantizer       quant;
    wrapper_in		   wrap_in;
    wrapper_out		   wrap_out;
    sc_error_counter   err_counter;
    sc_mbps_compute    mbps_compute;
    sc_xer_compute     ber_counter;
    sc_xer_compute     fer_counter;

    sc_fifo< TYPE_LLRS > f8;
    sc_fifo< TYPE_BITS > f9;
    sc_fifo< BIT       > f1;
    sc_fifo< float     > f2;
    sc_fifo< float     > f3;
    sc_fifo< float     > f4;
    sc_fifo< float     > f5;
    sc_fifo< float     > f6;
    sc_fifo< LLR       > f7;
    sc_fifo< BIT       > f10;
    sc_fifo< BIT       > f11;
    sc_fifo< BIT > FB;
    sc_signal< BIT > s1;
    sc_signal< BIT > s2;
    sc_signal< BIT > dec_frame;
    sc_signal< BIT > err_frame;

    sc_signal< sc_uint<64> > _proc_bits;
    sc_signal< sc_uint<64> > _proc_frames;
    sc_signal< sc_uint<64> > _err_bits;
    sc_signal< sc_uint<64> > _err_frames;

    sc_signal< sc_uint<32> > __proc_bits;
    sc_signal< sc_uint<32> > __proc_frames;
    sc_signal< sc_uint<32> > __err_bits;
    sc_signal< sc_uint<32> > __err_frames;

};

#endif
