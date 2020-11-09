#include <iostream>
#include<string>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "src/Writer.h"

using namespace std;

int main(int argc, char* argv[]) {

    std::cout << "(II) USER GUIDE : command N K P En IFile Input OPath" << std::endl;
    std::cout << "(II)      Parameters : " << std::endl;
    std::cout << "(II)                  - N          : (long) Frame size " << std::endl;
    std::cout << "(II)                  - K          : (long) Number of Information Bits " << std::endl;
    std::cout << "(II)                  - P          : (int) Parallelization factor (Number of parallel PU) " << std::endl;
    std::cout << "(II)                  - En         : (bool) 1 : Write the Frozen Bit table in the (PAR) concatenated form, 0 : otherwise" << std::endl;
    std::cout << "(II)                  - IFile      : (string) File of the reference Channel Fiability Order  (ex : 5G_N1024.txt) " << std::endl;
    std::cout << "(II)                  - Input      : (bool) 1 : if IFile contains frozen bit, 0 : if IFile contains Channel Fiabiliy Order " << std::endl;
    std::cout << "(II)                  - OPath      : (string) Path to the destination file (ex : /c/Polar_decoder_HLS/SC_Polar_decoder_v0/src/module/ ) " << std::endl;

    if (argc < 8)
        return -1;

    long NBIT = atoi(argv [1]);
    long KBIT = atoi(argv [2]);
    int PAR = atoi(argv [3]);
    bool En = (bool) atoi(argv [4]);
    string i_filename = argv [5];
    bool i_file_case = atoi(argv [6]);
    string o_path = argv [7];

    string o_filename = o_path + "polar_parameters.h";

    string affect_file = "../../Frozen_Bit_Tab/FB_N" + to_string(NBIT) + "_K" + to_string(KBIT) + ".txt" ;

    /***************************************************************************/

    Writer::Generate_FB_File(i_filename, NBIT, o_filename, KBIT, affect_file, PAR, En, i_file_case);

    /***************************************************************************/

    std::cout << "  File Generated :" << std::endl;

    std::cout << "  N           : " << NBIT << std::endl;
    std::cout << "  K           : " << KBIT << std::endl;
    std::cout << "  P           : " << PAR << std::endl;
    std::cout << "  En          : " << En << std::endl;
    std::cout << "  iFile       : " << i_filename << std::endl;
    std::cout << "  Input       : " << i_file_case << std::endl;
    std::cout << "  oFile       : " << o_filename << std::endl;
    std::cout << "  AffectFile  : " << affect_file << std::endl;

    exit( 0 );
}

