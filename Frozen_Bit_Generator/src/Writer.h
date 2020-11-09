//
// Created by ydelomier on 15/11/2018.
//
#ifndef FROZEN_BITS_WRITER_H
#define FROZEN_BITS_WRITER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include "math.h"

#define NMAX 32768  //2097152

using namespace std;

class Writer{
public:

    static void Generate_FB_File(const std::string i_filename, long NBIT, const std::string o_filename, long KBIT, const std::string affect_file, int PAR, bool En, bool i_file_case) {

        ifstream file(i_filename);
        if (file.is_open()) {
            string str;
            long N;
            long Frozen_Bit[NBIT];
            if( i_file_case == 0) {

                long channel_order[NMAX];
                long sub_channel_order[NBIT];
                //
                // ON RECUPERE LA TAILLE N
                //
                getline(file, str);
                istringstream line_1(str);
                line_1 >> N;

                //
                // ON SAUTE LES 2 LIGNES SUIVANTES
                //
                getline(file, str);
                istringstream line_2(str);
                getline(file, str);
                istringstream line_3(str);

                //
                // ON RECUPERE L'ORDRE DES CANAUX PAR ORDRE DE FIABILITE DECROISSANTE
                //
                getline(file, str);
                istringstream channels(str);

                long x = 0;
                int tempvalue;
                while (channels >> tempvalue) {
                    channel_order[x] = tempvalue;
                    x += 1;
                }

                //
                // ON RECUPERE SEULEMENT LES CANAUX INFERIEUR A NBIT
                //
                long k = 0;
                for (long i = 0; i < x; i++) {
                    int value = channel_order[i];
                    if (value < NBIT) {
                        sub_channel_order[k] = value;
                        k += 1;
                    }
                }

                //
                // ON ECRIT LE FICHIER POUR AFFECT
                //
                ofstream affect(affect_file);
                affect << NBIT << endl << "0" << endl << "0" << endl;
                for (long i = 0; i < NBIT; i++) {
                    affect << sub_channel_order[i] << "    ";
                }
                affect.close();

                //
                // ON DEFINIT LES BITS GELE EN FONCTION DE KBIT
                //

                for (long i = 0; i < KBIT; i++) {
                    int indice = sub_channel_order[i];
                    Frozen_Bit[indice] = 1;            // For the KBIT most reliable channels write 1
                }
                for (long i = KBIT; i < NBIT; i++) {
                    int indice = sub_channel_order[i];
                    Frozen_Bit[indice] = 0;           // For the KBIT least reliable channels write 0
                }
            }
            else
            {
                getline(file, str);
                istringstream channels(str);
                long x = 0;
                int tempvalue;
                while (channels >> tempvalue) {
                    Frozen_Bit[x] = tempvalue;
                    x += 1;
                }
            }
            //
            // ON ECRIT LE FICHIER DE SORTIE
            //

            ofstream o_file(o_filename);
            o_file << "#ifndef POLAR_HEADER_H"<< endl;
            o_file << "#define POLAR_HEADER_H"<< endl;
            o_file << endl;
            o_file << "#define _NBITS       " << NBIT << endl;
            o_file << "#define _LOG2N       " << (log2(NBIT)) << endl;
            o_file << "#define _DEPTH       " << (log2(NBIT) + 1) << endl;
            o_file << endl;
            o_file << "#define PAR          " << PAR << endl;
            o_file << "#define LOG2_PAR     " << log2(PAR) << endl;
            o_file << "#define N_DIVIDED    (_NBITS / PAR) "<< endl;
			o_file << "#define DEPTH_DIV    " << (log2(NBIT / PAR) + 1) << endl;
            o_file << endl;
            o_file << "#define COUNTER      sc_uint<_DEPTH>"<< endl;
            o_file << endl;
            if(En == true) {
                o_file << "const sc_bv<PAR> Frozen_Bits[N_DIVIDED] = {" << endl;
                o_file << "   //";
                for (long i = 0; i < NBIT; i++) {
                    o_file << Frozen_Bit[i] << ", ";
                }
                o_file << std::endl << "     \"";

                for (long i = 0; i < (NBIT / PAR); i++) {
                    for (long j = 0; j < PAR; j++) {
                        o_file << Frozen_Bit[((i + 1) * PAR - 1) - j];
                    }
                    o_file << "\", \"";
                }
                o_file.seekp(-3, std::ios_base::cur);
                o_file << std::endl;
                o_file << "};" << std::endl << std::endl;
            }
            else{
                o_file << "const sc_bv<1> Frozen_Bits[_NBITS] = {" << endl;
                o_file << "   // \"";
                for (long i = 0; i < (NBIT / PAR); i++) {
                    for (long j = 0; j < PAR; j++) {
                        o_file << Frozen_Bit[((i + 1) * PAR - 1) - j];
                    }
                    o_file << "\", \"";
                }
                o_file << std::endl << "     ";
                for (long i = 0; i < NBIT; i++) {
                    o_file << Frozen_Bit[i] << ", ";
                }
                o_file.seekp(-2, std::ios_base::cur);
                o_file << std::endl;
                o_file << "};" << std::endl << std::endl;
            }
            o_file << std::endl;
            o_file << "#endif // POLAR_HEADER_H"<< endl;
            o_file.close();
        }
        else{
            cout << "!!! ERROR file does not exist : " << i_filename << " !!!" << endl;
            exit(0);
        }
        file.close();

        cout << "fin" << endl;
    }

};

#endif //FROZEN_BITS_WRITER_H
