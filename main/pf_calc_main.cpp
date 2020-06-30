/* 

	Protein packing fraction calc main

	Read in pdb, file location, calculate packing fraction and print to specified file

*/

// include file
#include "proteinPacking.h"

// namespace
using namespace std;

// main function
int main(int argc, char *argv[]){
	// read in data
	string inputFileStr 	= argv[1]; 			// input pdb string
	string outputFileStr 	= argv[2]; 			// output pf file string
	string nmcptsStr 		= argv[3]; 		// number of monte carlo points for particle mass

	// MC variables
	double nmcpts_dbl;
	int nmcpts;

	stringstream nmcptsss(nmcptsStr);
	nmcptsss >> nmcpts_dbl;
	nmcpts = (int)nmcpts_dbl;

	// instantiate object
	cout << "Instantiating object..." << endl;
	proteinPacking pobj(inputFileStr);

	// compute masses
	cout << "Calculating masses..." << endl;
	pobj.calcMasses(nmcpts);

	cout << "Calculating neighbors..." << endl;
	pobj.neighbors();

	// print packing fraction data to file
	cout << "Opening packing file " << outputFileStr << endl;
	pobj.openPackingFile(outputFileStr);

	// print data to file
	cout << "Printing file." << endl;
	pobj.printPackingFraction();

	// end main
	cout << "Ending main." << endl;
	return 0;
}