#ifndef _sc_monitor_
#define _sc_monitor_

#include "systemc.h"
#include "math.h"
#include <iostream>
#include <iomanip>      // std::setw

#ifdef __RTL_SIMULATION__
#include "config.h"
#include "polar_parameters.h"
#else
#include "../../module/config.h"
#include "../../module/polar_parameters.h"
#endif

#define NTH_ITERATION 3

SC_MODULE(sc_monitor)
{
  sc_in <bool> clk;
  sc_in <bool> reset;

  sc_fifo_in < TYPE_BITS > s;

  sc_in <bool> busy;
  sc_in <bool> load;
  sc_in <bool> store;
  sc_in <COUNTER> cnt_i;
  sc_in <int> N_value;
  sc_in < sc_uint<8> > Fct_ID;  // F: 0x88, G: 0x48, H: 0x28, R: 0x18
  	  	  	  	  	  	  	  	// R0: 0x00, R1: 0x0F, REP: 0x02, SPC: 0x04

  SC_CTOR(sc_monitor)
  {
	SC_CTHREAD(read_fifo, clk.pos());
	reset_signal_is(reset, true);

    SC_CTHREAD(do_busy, clk.pos());
    reset_signal_is(reset, true);

    SC_CTHREAD(do_load, clk.pos());
    reset_signal_is(reset, true);

    SC_CTHREAD(do_store, clk.pos());
    reset_signal_is(reset, true);
  }

private:
	void do_busy(){

		int ctn_proc = 0;
		int Lat_level[DEPTH_DIV] = {0};
		int Lat_function[128] = {0};
		int latencies[128][DEPTH_DIV] = {0};
		int nb_occurences[24][DEPTH_DIV] = {0};
		bool one_time = false;

		for(int i = 0; i < DEPTH_DIV; i++)
			Lat_level[i] = 0;
		for(int i = 0; i < 128; i++)
			Lat_function[i] = 0;
		for(int i = 0; i < 128; i++){
			for(int j = 0; j < DEPTH_DIV; j++){
				latencies[i][j] = 0;
				nb_occurences[i][j] = 0;
			}
		}
		wait();
		while(true)
		{
			if(one_time == false)
			{
				while(busy.read() == 0)
					wait();

				int last_F_ID = -1;
				int last_Log2N = -1;

				while(busy.read() == 1){

					int N = N_value.read();
					int Log2N = log2( N );
					Lat_level[ Log2N ] += 1;
					sc_uint<8> F_ID = Fct_ID.read() ;
					Lat_function[ F_ID ] += 1;
					ctn_proc++;
					latencies[F_ID][Log2N] +=1;
					if( (F_ID != last_F_ID) || (last_Log2N != Log2N) )
						nb_occurences[F_ID][Log2N] +=1;
					last_F_ID = F_ID;
					last_Log2N = Log2N;
					wait();
				}

				cout << "[MONITOR] PROCESS latency : " << ctn_proc <<endl<< endl;
				cout << "		*** Latency by Level ***  " << endl;
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << "				level " << DEPTH_DIV-1-i << " : "<< Lat_level[ DEPTH_DIV-1-i] <<endl;
				cout << endl;
				cout << "		*** Latency by Function *** " << endl;
				cout << "				function F     : "<< Lat_function[0x88] <<endl;
				cout << "				function G     : "<< Lat_function[0x48] <<endl;
				cout << "				function H     : "<< Lat_function[0x28] <<endl;
				cout << "				function R     : "<< Lat_function[0x18] <<endl;
				cout << "				......................" << endl;
#if PRUNING_LEVEL == 1
				cout << "				function R_R0  : "<< Lat_function[0x10] <<endl;
	#if ELAG_R1 ==1
				cout << "				function R_R1  : "<< Lat_function[0x1F] <<endl;
	#endif
	#if ELAG_REP ==1
				cout << "				function R_REP : "<< Lat_function[0x12] <<endl;
	#endif
	#if ELAG_REP2 == 1
				cout << "				function R_REP2: "<< Lat_function[0x13] <<endl;
	#endif
	#if ELAG_SPC ==1
				cout << "				function R_SPC : "<< Lat_function[0x14] <<endl;
	#endif
	#if ELAG_SPC2 == 1
				cout << "				function R_SPC2: "<< Lat_function[0x15] <<endl;
	#endif
#elif PRUNING_LEVEL == 2
	#if ELAG_H0 == 0
				cout << "				function F_R0  : "<< Lat_function[0x80] <<endl;
	#else
				cout << "				function H_R0  : "<< Lat_function[0x20] <<endl;
				cout << "				......................" << endl;
	#endif
	#if ELAG_R1 ==1 && ELAG_RARE == 1
				cout << "				function F_R1  : "<< Lat_function[0x8F] <<endl;
	#endif
	#if ELAG_REP ==1
				cout << "				function F_REP : "<< Lat_function[0x82] <<endl;
	#endif
	#if ELAG_REP2 == 1
				cout << "				function F_REP2: "<< Lat_function[0x83] <<endl;
	#endif
	#if ELAG_SPC ==1 && ELAG_RARE == 1
				cout << "				function F_SPC : "<< Lat_function[0x84] <<endl;
	#endif
	#if ELAG_SPC2 == 1 && ELAG_RARE == 1
				cout << "				function F_SPC2: "<< Lat_function[0x85] <<endl;
	#endif
				cout << "				......................" << endl;
	#if ELAG_RARE == 1
				cout << "				function G_R0  : "<< Lat_function[0x40] <<endl;
	#endif
	#if ELAG_R1 ==1
				cout << "				function G_R1  : "<< Lat_function[0x4F] <<endl;
	#endif
	#if ELAG_REP ==1 && ELAG_RARE == 1
				cout << "				function G_REP : "<< Lat_function[0x42] <<endl;
	#endif
	#if ELAG_REP2 == 1 && ELAG_RARE == 1
				cout << "				function G_REP2: "<< Lat_function[0x43] <<endl;
	#endif
	#if ELAG_SPC ==1
				cout << "				function G_SPC : "<< Lat_function[0x44] <<endl;
	#endif
	#if ELAG_SPC2 == 1
				cout << "				function G_SPC2: "<< Lat_function[0x45] <<endl;
	#endif

#endif


				cout << endl;
				cout << "		*** Matrix Level / Function Latency ***  "<< endl;
				cout << "				Level     |";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << std::setw(7) << pow(2,_DEPTH-1-i) ;
				cout << endl << "				__________|";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << "_______";
				cout << endl << "				fct F     |";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << std::setw(7) << latencies[0x88][DEPTH_DIV-1-i] ;
				cout << endl << "				fct G     |";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << std::setw(7) << latencies[0x48][DEPTH_DIV-1-i] ;
				cout << endl << "				fct H     |";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << std::setw(7) << latencies[0x28][DEPTH_DIV-1-i] ;
				cout << endl << "				fct R     |";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << std::setw(7) << latencies[0x18][DEPTH_DIV-1-i] ;
				cout << endl << "				..........|";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << ".......";
#if PRUNING_LEVEL == 1
				cout << endl << "				fct R_R0  |";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << std::setw(7) << latencies[0x10][DEPTH_DIV-1-i] ;
	#if ELAG_R1 ==1
				cout << endl << "				fct R_R1  |";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << std::setw(7) << latencies[0x1F][DEPTH_DIV-1-i] ;
	#endif
	#if ELAG_REP ==1
				cout << endl << "				fct R_REP |";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << std::setw(7) << latencies[0x12][DEPTH_DIV-1-i] ;
	#endif
	#if ELAG_REP2 == 1
				cout << endl << "				fct R_REP2|";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << std::setw(7) << latencies[0x13][DEPTH_DIV-1-i] ;
	#endif
	#if ELAG_SPC ==1
				cout << endl << "				fct R_SPC |";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << std::setw(7) << latencies[0x14][DEPTH_DIV-1-i] ;
	#endif
	#if ELAG_SPC2 == 1
				cout << endl << "				fct R_SPC2|";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << std::setw(7) << latencies[0x15][DEPTH_DIV-1-i] ;
	#endif
#elif PRUNING_LEVEL == 2
	#if ELAG_H0 == 0
				cout << endl << "				fct F_R0  |";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << std::setw(7) << latencies[0x80][DEPTH_DIV-1-i] ;
	#else
				cout << endl << "				fct H_R0  |";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << std::setw(7) << latencies[0x20][DEPTH_DIV-1-i] ;
				cout << endl << "				..........|";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << ".......";
	#endif
	#if ELAG_R1 ==1 && ELAG_RARE == 1
				cout << endl << "				fct F_R1  |";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << std::setw(7) << latencies[0x8F][DEPTH_DIV-1-i] ;
	#endif
	#if ELAG_REP ==1
				cout << endl << "				fct F_REP |";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << std::setw(7) << latencies[0x82][DEPTH_DIV-1-i] ;
	#endif
	#if ELAG_REP2 == 1
				cout << endl << "				fct F_REP2|";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << std::setw(7) << latencies[0x83][DEPTH_DIV-1-i] ;
	#endif
	#if ELAG_SPC ==1 && ELAG_RARE == 1
				cout << endl << "				fct F_SPC |";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << std::setw(7) << latencies[0x84][DEPTH_DIV-1-i] ;
	#endif
	#if ELAG_SPC2 == 1 && ELAG_RARE == 1
				cout << endl << "				fct F_SPC2|";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << std::setw(7) << latencies[0x85][DEPTH_DIV-1-i] ;
	#endif
				cout << endl << "				..........|";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << ".......";
	#if ELAG_RARE == 1
				cout << endl << "				fct G_R0  |";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << std::setw(7) << latencies[0x40][DEPTH_DIV-1-i] ;
	#endif
	#if ELAG_R1 ==1
				cout << endl << "				fct G_R1  |";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << std::setw(7) << latencies[0x4F][DEPTH_DIV-1-i] ;
	#endif
	#if ELAG_REP ==1 && ELAG_RARE == 1
				cout << endl << "				fct G_REP |";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << std::setw(7) << latencies[0x42][DEPTH_DIV-1-i] ;
	#endif
	#if ELAG_REP2 == 1 && ELAG_RARE == 1
				cout << endl << "				fct G_REP2|";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << std::setw(7) << latencies[0x43][DEPTH_DIV-1-i] ;
	#endif
	#if ELAG_SPC ==1
				cout << endl << "				fct G_SPC |";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << std::setw(7) << latencies[0x44][DEPTH_DIV-1-i] ;
	#endif
	#if ELAG_SPC2 == 1
				cout << endl << "				fct G_SPC2|";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << std::setw(7) << latencies[0x45][DEPTH_DIV-1-i] ;
	#endif
#endif
				cout << endl;
				cout << endl;
				cout << "		*** Matrix Level / Function Occurences ***  "<< endl;
				cout << "				Level     |";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << std::setw(7) << pow(2,_DEPTH-1-i) ;
				cout << endl << "				__________|";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << "_______";
				cout << endl << "				fct F     |";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << std::setw(7) << nb_occurences[0x88][DEPTH_DIV-1-i] ;
				cout << endl << "				fct G     |";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << std::setw(7) << nb_occurences[0x48][DEPTH_DIV-1-i] ;
				cout << endl << "				fct H     |";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << std::setw(7) << nb_occurences[0x28][DEPTH_DIV-1-i] ;
				cout << endl << "				fct R     |";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << std::setw(7) << nb_occurences[0x18][DEPTH_DIV-1-i] ;
				cout << endl << "				..........|";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << ".......";
#if PRUNING_LEVEL == 1
				cout << endl << "				fct R_R0  |";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << std::setw(7) << nb_occurences[0x10][DEPTH_DIV-1-i] ;
	#if ELAG_R1 ==1
				cout << endl << "				fct R_R1  |";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << std::setw(7) << nb_occurences[0x1F][DEPTH_DIV-1-i] ;
	#endif
	#if ELAG_REP ==1
				cout << endl << "				fct R_REP |";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << std::setw(7) << nb_occurences[0x12][DEPTH_DIV-1-i] ;
	#endif
	#if ELAG_REP2 == 1
				cout << endl << "				fct R_REP2|";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << std::setw(7) << nb_occurences[0x13][DEPTH_DIV-1-i] ;
	#endif
	#if ELAG_SPC ==1

				cout << endl << "				fct R_SPC |";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << std::setw(7) << nb_occurences[0x14][DEPTH_DIV-1-i] ;
	#endif
	#if ELAG_SPC2 == 1
				cout << endl << "				fct R_SPC2|";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << std::setw(7) << nb_occurences[0x15][DEPTH_DIV-1-i] ;

	#endif
#elif PRUNING_LEVEL == 2
	#if ELAG_H0 == 0
				cout << endl << "				fct F_R0  |";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << std::setw(7) << nb_occurences[0x80][DEPTH_DIV-1-i] ;
	#else
				cout << endl << "				fct H_R0  |";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << std::setw(7) << nb_occurences[0x20][DEPTH_DIV-1-i] ;
				cout << endl << "				..........|";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << ".......";
	#endif
	#if ELAG_R1 ==1 && ELAG_RARE == 1
				cout << endl << "				fct F_R1  |";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << std::setw(7) << nb_occurences[0x8F][DEPTH_DIV-1-i] ;
	#endif
	#if ELAG_REP ==1

				cout << endl << "				fct F_REP |";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << std::setw(7) << nb_occurences[0x82][DEPTH_DIV-1-i] ;
	#endif
	#if ELAG_REP2 == 1
				cout << endl << "				fct F_REP2|";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << std::setw(7) << nb_occurences[0x83][DEPTH_DIV-1-i] ;
	#endif
	#if ELAG_SPC ==1 && ELAG_RARE == 1
				cout << endl << "				fct F_SPC |";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << std::setw(7) << nb_occurences[0x84][DEPTH_DIV-1-i] ;
	#endif
	#if ELAG_SPC2 == 1 && ELAG_RARE == 1
				cout << endl << "				fct F_SPC2|";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << std::setw(7) << nb_occurences[0x85][DEPTH_DIV-1-i] ;
	#endif
				cout << endl << "				..........|";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << ".......";
	#if ELAG_RARE == 1
				cout << endl << "				fct G_R0  |";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << std::setw(7) << nb_occurences[0x40][DEPTH_DIV-1-i] ;
	#endif
	#if ELAG_R1 ==1
				cout << endl << "				fct G_R1  |";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << std::setw(7) << nb_occurences[0x4F][DEPTH_DIV-1-i] ;
	#endif
	#if ELAG_REP ==1 && ELAG_RARE == 1
				cout << endl << "				fct G_REP |";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << std::setw(7) << nb_occurences[0x42][DEPTH_DIV-1-i] ;
	#endif
	#if ELAG_REP2 == 1 && ELAG_RARE == 1
				cout << endl << "				fct G_REP2|";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << std::setw(7) << nb_occurences[0x43][DEPTH_DIV-1-i] ;
	#endif
	#if ELAG_SPC ==1

				cout << endl << "				fct G_SPC |";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << std::setw(7) << nb_occurences[0x44][DEPTH_DIV-1-i] ;
	#endif
	#if ELAG_SPC2 == 1
				cout << endl << "				fct G_SPC2|";
				for(int i = 0; i < DEPTH_DIV; i++)
					cout << std::setw(7) << nb_occurences[0x45][DEPTH_DIV-1-i] ;
	#endif
#endif
				cout << endl;

				ctn_proc = 0;
				for(int i = 0; i < DEPTH_DIV; i++)
					Lat_level[i] = 0;
				for(int i = 0; i < 128; i++)
					Lat_function[i] = 0;
				for(int i = 0; i < 128; i++){
					for(int j = 0; j < DEPTH_DIV; j++){
						latencies[i][j] = 0;
						nb_occurences[i][j] = 0;
					}
				}
			}

			one_time = true;
			wait();
		}

	}

	void do_load(){

		int ctn_load = 0;
		bool one_time = false;
		wait();
		while(true)
		{
			if(one_time == false)
			{
				while(load.read() == 0)
					wait();

				while(load.read() == 1){
					ctn_load++;
					wait();
				}

				cout << "[MONITOR] LOAD latency : " << ctn_load <<endl;
				fflush( stdout );

				ctn_load = 0;
			}

			one_time = true;
			wait();
		}
	}

	void do_store(){

		int ctn_store = 0;
		bool one_time = false;
		wait();
		while(true)
		{
			if(one_time == false)
			{
				while(store.read() == 0)
					wait();

				while(store.read() == 1){
					ctn_store++;
					wait();
				}

				cout << "[MONITOR] STORE latency : " << ctn_store <<endl;
				fflush( stdout );

				ctn_store = 0;
			}

			one_time = true;
			wait();
		}
	}

	void read_fifo(){

		wait();
		while(true)
		{
			TYPE_BITS fifo_in = s.read();
		}
	}

	void read_dbg(){

		wait();
		while(true)
		{
			COUNTER cnt = cnt_i.read();
		}
	}
};

#endif
