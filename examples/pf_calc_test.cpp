/* 

protein packing test main

*/

#include "proteinPacking.h"

using namespace std;

int main()
{
	// local variables
	string inputFileStr = "/Users/JackTreado/_pv/pdb/dunbrack/3ccd_H.dat";
	string outputFileStr = "pf.test";

	// instantiate object
	cout << "Instantiating object..." << endl;
	proteinPacking pobj(inputFileStr);

	// print packing fraction data to file
	cout << "Opening packing file " << outputFileStr << endl;
	pobj.openPackingFile(outputFileStr);

	cout << "Printing file." << endl;
	pobj.printPackingFraction();


	cout << "Ending test code." << endl;
	return 0;
}