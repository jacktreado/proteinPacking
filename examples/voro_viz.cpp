/* 

protein packing test main

*/

#include "proteinPacking.h"

using namespace std;

int main()
{
	// local variables
	int i, NR;

	// file strings
	string inputFileStr 	= "/Users/JackTreado/_pv/boxmethod/data/8abp_H.dat";
	string voronoiFileStr 	= "/Users/JackTreado/Jamming/ProteinVoids/pdb/vorodata/data/8abp_H.voro";

	// instantiate object
	cout << "Instantiating object..." << endl;
	proteinPacking pobj(inputFileStr);

	// open voronoi file
	cout << "	** Opening file : " << voronoiFileStr << endl;
	pobj.openVoronoiFile(voronoiFileStr);

	// print neighbors to console
	cout << "Printing neighbor information...." << endl;
	int N = pobj.getN();
	int NN;
	for (int i=0; i<N; i++){
		cout << "neighbors of i = " << i << "  :  ";
		for (int j=0; j<pobj.numResNeighbors(i); j++)
			cout << setw(6) << pobj.resNeighbors(i,j);
		cout << endl;
	}

	// loop over specific residues, print voronoi information
	vector<int> residues;
	// residues.push_back(78);
	// residues.push_back(103);
	for (i=0; i<N; i++){
		if (pobj.resVoro(i) > 0)
			residues.push_back(i);
	}
	NR = residues.size();

	// loop over residues, print voronoi information
	pobj.printToVoroFile(NR);
	pobj.printToVoroFile(pobj.getBound(0,0));
	pobj.printToVoroFile(pobj.getBound(0,1));
	pobj.printToVoroFile(pobj.getBound(1,0));
	pobj.printToVoroFile(pobj.getBound(1,1));
	pobj.printToVoroFile(pobj.getBound(2,0));
	pobj.printToVoroFile(pobj.getBound(2,1));
	for (i=0; i<NR; i++){
		cout << "** printing voronoi cell for " << residues.at(i) << endl;
		pobj.printSingleVoronoiCell(residues.at(i));
	}

	// end
	cout << "Ending voro viz code." << endl;
	return 0;
}