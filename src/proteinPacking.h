#ifndef PROTEINPACKING_H
#define PROTEINPACKING_H

/*
		--- PROTEINPACKING CLASS ---

	Class for protein objects to calculate 
	protein packing fractins using voro++ 
	(radical tessellation)

	saves atomic coordinates, voronoi volumes

	prints to a variety of files
*/


#include "voro++.hh"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include <fstream>
#include <cstdlib>
#include <cmath>

class proteinPacking{
private:
	// number of atoms
	int NA;

	// number of residues
	int N;

	// number of atoms in each residue
	int* residueSize;

	// sequence
	std::string* sequence;

	// atom types
	std::string* atoms;

	// atomic information
	double* x;			// positions
	double* r;			// radii
	double* m;			// masses

	// voronoi information	
	voro::container_poly* voroContainer;
	double* voroVols;
	std::vector<int>* voroNeighbors;

	// file information
	std::ofstream packingFile;
public:
	// constructors
	proteinPacking();
	proteinPacking(std::string& inputFileString, int NMCPTS);

	// destructors
	~proteinPacking();

	// file open
	void openPackingFile(std::string& str){
		packingFile.open(str.c_str());
		if (!packingFile.is_open()){
			std::cout << "	ERROR: packing file not open, file string " << str << " is not functioning, ending." << std::endl; exit(1);
		}
	}

	// setters
	void setPos(int residue, int atom, int d, double val);
	void setPos(int atom, int d, double val);
	void setRad(int residue, int atom, double val);
	void setRad(int atom, double val);
	void setAtom(int atom, std::string& str);
	void setMass(int residue, double val);
	void setResidueSize(int residue, int val);
	void setSequence(int residue, std::string& str);
	void setVoroVol(int residue, double val);

	// getters
	int cumulativeNumberOfAtoms(int residue);
	int atomID(int residue, int atom);
	int resID(int atom);
	double pos(int residue, int atom, int d);
	double pos(int atom, int d);
	double rad(int residue, int atom);
	double rad(int atom);
	double mass(int residue);
	int size(int residue);
	double resVoro(int residue);
	int resNeighbors(int residue, int neighbor);
	int numResNeighbors(int residue);
	std::string getSeq(int residue);
	std::string getAtom(int residue, int atom);
	std::string getAtom(int atom);

	// calculators
	void calcMasses(int NMCPTS);
	void neighbors();

	// printers
	void printPackingFraction();
};

#endif