//#include "./sc_top_module/sc_dual_top_module.h"

#include "sc_top_module.h"
#include <iomanip>

//#define DEBUG_MAIN

#ifdef __RTL_SIMULATION__
#include "sc_terminal.h"
#include "polar_parameters.h"
#else
#include "sc_terminal/sc_terminal.h"
#include "../module/polar_parameters.h"
#endif

#define _VSATN	-31
#define _VSATP	31
#define _BETA	4


int main (int argc, char * argv []){

	cout << "(II) Parametres de simulation : " << endl;
	cout << "     - Frame size: "  << (_NBITS)                       << endl;
	cout << "     - LLR width : "  << (LLR_BITS)                 << endl;
	cout << "     - QUANT  	  : [" << (_BETA) <<", " << (_VSATN) <<", " << (_VSATP)<< "]" << endl;

	sc_clock clock("clock", 10, SC_NS, 0.5);
	sc_signal< bool > reset ("8");

	sc_signal< float > 		sigma ("2");
	sc_signal< sc_uint<8> > beta ("3");
	sc_signal< sc_int<8> > 	vsatn ("4");
	sc_signal< sc_int<8> > 	vsatp ("5");
	sc_signal< sc_uint<8> > nbit ("6");
	sc_signal< sc_uint<16> > cpt ("7");

	sc_signal< bool       > enable ("9");
	//sc_signal< sc_uint<5> > snr;
	sc_signal< sc_uint<8> > seed ("10");

	sc_signal< sc_uint<64> > err_bits ("11");
	sc_signal< sc_uint<64> > err_frames ("12");
	sc_signal< sc_uint<64> > proc_bits ("13");
	sc_signal< sc_uint<64> > proc_frames ("14");
	sc_signal< sc_uint<32> > mbps_value ("15");
	sc_signal< sc_uint<32> > ber_value ("16");
	sc_signal< sc_uint<32> > fer_value ("17");

	sc_signal< BIT > err_frame ("18");

	sc_top_module top ( "top"  );
	top.clk      ( clock     );
	top.reset    ( reset     );

	top.enable   ( enable    );
	top.SIGMA    ( sigma    );
	top.seed     ( seed      );
	top.BETA (beta);
	top.VSATN (vsatn);
	top.VSATP (vsatp);
	top.NBIT (nbit);
	//top.CPT(cpt);

	top.err_bits   ( err_bits    );
	top.err_frames ( err_frames  );
	top.proc_bits  ( proc_bits   );
	top.proc_frames( proc_frames );
	top.ber_value  ( ber_value   );
	top.fer_value  ( fer_value   );
	top.mbps_value ( mbps_value  );


	sc_terminal     term( "term" );
	term.clk        ( clock       );
	term.reset      ( reset       );
	term.err_bits   ( err_bits    );
	term.err_frames ( err_frames  );
	term.proc_bits  ( proc_bits   );
	term.proc_frames( proc_frames );
	term.mbps_value ( mbps_value  );
	term.fer_value  ( fer_value   );

	cout << "(II) Reset the system" << endl;

	beta = _BETA;
	vsatn = _VSATN;
	vsatp = _VSATP;

	//snr    = 2.00f * (4); // x4 a cause du format de codage en vigule fixe
	float snr = 2.5f;
	float _R = 0.5;
#ifdef DEBUG_MAIN
	float sigma1 = 0.0f;
#else
	float sigma1 = 1.0f / sqrtf( 2.f * _R * powf( 10.f , snr/10.f ) );
#endif
	sigma = sigma1;
	cout << "		(!!!) Sigma = " <<std::fixed << std::setprecision(4) << sigma1 << " (!!!)" << endl;

	enable = true;
	reset  = true;
	seed   = 0xF0;
	reset  = false;
	sc_start( 100, SC_NS );

	cout << "(II) Running the system" << endl;
#ifdef DEBUG_MAIN
	sc_start( 500, SC_US );
#else
	#ifdef __RTL_SIMULATION__
		sc_start(1, SC_US);
	#else
		sc_start(3, SC_MS);
#endif
#endif

  return 0;
}
