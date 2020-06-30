/*

	proteinPacking methods

*/

// include
#include "proteinPacking.h"

// namespace
using namespace std;
using namespace voro;

// global constants
#define NCELLS 3
#define NINIT 5
#define NDIM 3
const double PI = 4*atan(1);

/****************************
	
	Initialization

*****************************/

// constructors
proteinPacking::proteinPacking(){
	// set scalars to 0
	NA = 0;
	N = 0;

	// point pointers to null
	residueSize = nullptr;
	sequence = nullptr;
	atoms = nullptr;
	x = nullptr;
	r = nullptr;
	m = nullptr;
	voroContainer = nullptr;
	voroVols = nullptr;
	voroNeighbors = nullptr;
}

proteinPacking::proteinPacking(string& inputFileString){
	// local variables
	int i, j;
	double rtmp,xtmp,ytmp,ztmp,xmin,xmax,ymin,ymax,zmin,zmax,totalmin,totalmax,boxscale;
	int idtmp,atomCounter,resCounter,oldId;
	string restmp,atomtmp;
	char chaintmp;
	voro::voronoicell_neighbor c;
	int residueIndex;
	double boxChangeTol = 1e-4;


	// open input file object
	ifstream inputObject(inputFileString.c_str());
	if (!inputObject.is_open()){
		cout << "	ERROR: cannot open input file object, ending." << endl; exit(1);
	}

	// get file header input information
	inputObject >> NA;
	inputObject >> N;

	// check that numbers are non size
	if (NA == 0){
		cout << "	ERROR: NA read in as 0 from file, ending." << endl; exit(1);
	}
	else if (N == 0){
		cout << "	ERROR: N read in as 0 from file, ending." << endl; exit(1);
	}

	// initialize dynamically allocated memory

	// atomic information
	x 		= new double[NA*NDIM];
	r 		= new double[NA];	
	atoms 	= new string[NA];

	if (!x){
		cout << "	ERROR: cannot allocate memory for x array, ending." << endl; exit(1);
	}
	else if (!r){
		cout << "	ERROR: cannot allocate memory for r array, ending." << endl; exit(1);
	}
	else if (!atoms){
		cout << "	ERROR: cannot allocate memory for atoms array, ending." << endl; exit(1);
	}

	// residue information
	m 				= new double[N];
	residueSize 	= new int[N];
	sequence 		= new string[N];
	voroVols 		= new double[N];
	voroNeighbors 	= new vector<int>[N];

	if (!m){
		cout << "	ERROR: cannot allocate memory for m array, ending." << endl; exit(1);
	}
	else if (!residueSize){
		cout << "	ERROR: cannot allocate memory for residueSize array, ending." << endl; exit(1);
	}
	else if (!sequence){
		cout << "	ERROR: cannot allocate memory for sequence array, ending." << endl; exit(1);
	}
	else if (!voroVols){
		cout << "	ERROR: cannot allocate memory for voroVols array, ending." << endl; exit(1);
	}
	else if (!voroNeighbors){
		cout << "	ERROR: cannot allocate memory for voroNeighbors array, ending." << endl; exit(1);
	}

	// counters
	atomCounter = 0;
	resCounter = 0;
	oldId = 1;
	xmin = 1e3; xmax = -1e3;
	ymin = 1e3; ymax = -1e3;
	zmin = 1e3; zmax = -1e3;

	// loop over file, read in values
	cout << "Reading in values from file " << inputFileString << endl;
	while (!inputObject.eof()){
		// read in data
		inputObject >> idtmp >> restmp >> atomtmp >> chaintmp >> rtmp >> xtmp >> ytmp >> ztmp;
		cout << "atom=" << atomCounter << ";  resname=" << restmp << ";  id=" << idtmp << ";  r=" << rtmp << ";  x=" << xtmp << ";  y=" << ytmp << ";  z=" << ztmp << endl;

		// determine res id
		if (oldId == idtmp-1){
			// store residue information
			setMass(oldId-1,0.0);
			setResidueSize(oldId-1,resCounter);
			setVoroVol(oldId-1,0.0);

			// reset counters
			resCounter = 0;
			oldId = idtmp;			
		}
		if (oldId == idtmp && resCounter == 0)
			setSequence(oldId-1,restmp);

		// save to values
		setPos(atomCounter,0,xtmp);
		setPos(atomCounter,1,ytmp);
		setPos(atomCounter,2,ztmp);
		setRad(atomCounter,rtmp);
		setAtom(atomCounter,atomtmp);
		
		// update atoms on residue counter
		resCounter++;

		// update atom counter
		atomCounter++;
		if (atomCounter == NA){
			cout << "All atoms read in from file, ending loop." << endl;
			setMass(N-1,0.0);
			setResidueSize(N-1,resCounter);
			setVoroVol(N-1,0.0);
			break;
		}

		// determine minimum and maximum values of position
		if (xtmp > xmax)
			xmax = xtmp;
		else if (xtmp < xmin)
			xmin = xtmp;

		if (ytmp > ymax)
			ymax = ytmp;
		else if (ytmp < ymin)
			ymin = ytmp;

		if (ztmp > zmax)
			zmax = ztmp;
		else if (ztmp < zmin)
			zmin = ztmp;
	}
	cout << "Done reading in values from file " << inputFileString << endl;
	cout << "xmin = " << xmin << ", xmax = " << xmax << endl;
	cout << "ymin = " << ymin << ", ymax = " << ymax << endl;
	cout << "zmin = " << zmin << ", zmax = " << zmax << endl;

	// get totalmin
	totalmin = xmin;
	if (xmin < ymin){
		if (xmin < zmin)
			totalmin = xmin;
		else
			totalmin = zmin;
	}
	else{
		if (ymin < zmin)
			totalmin = ymin;
		else
			totalmin = zmin;
	}

	// get totalmax
	totalmax = xmax;
	if (xmax > ymax){
		if (xmax > zmax)
			totalmax = xmax;
		else
			totalmax = zmax;
	}
	else{
		if (ymax > zmax)
			totalmax = ymax;
		else
			totalmax = zmax;
	}

	// scale box sizes so they are a little off of the protein
	if (abs(totalmax) > abs(totalmin))
		boxscale = 0.25*abs(totalmax);
	else
		boxscale = 0.25*abs(totalmin);

	totalmin -= boxscale;
	totalmax += boxscale;
	cout << "totalmin = " << totalmin << ", totalmax = " << totalmax << endl;

	// save box boundaries linearly (xmin, xmax, ymin, ymax, zmin, zmax)
	bounds.resize(2*NDIM);
	bounds.at(0) = totalmin;
	bounds.at(1) = totalmax;
	bounds.at(2) = totalmin;
	bounds.at(3) = totalmax;
	bounds.at(4) = totalmin;
	bounds.at(5) = totalmax;

	// setup container for voronoi
	cout << "Setting up container for voronoi calculation." << endl;

	// instantiate containerp object
	cout << "Getting voronoi volumes for normal box size. " << endl;
	voroContainer = new container_poly(totalmin,totalmax,totalmin,totalmax,totalmin,totalmax,NCELLS,NCELLS,NCELLS,false,false,false,NINIT);

	// instantiate slightly larger box, find surface residues
	totalmin -= 0.25*boxscale;
	totalmax += 0.25*boxscale;
	cout << "Getting voronoi volumes for larger box size. " << endl;
	container_poly largerBox(totalmin,totalmax,totalmin,totalmax,totalmin,totalmax,NCELLS,NCELLS,NCELLS,false,false,false,NINIT);

	// loop over atoms and put into container
	cout << "Putting atoms into both boxes." << endl;
	for (i=0; i<NA; i++){
		voroContainer->put(i,pos(i,0),pos(i,1),pos(i,2),rad(i));
		largerBox.put(i,pos(i,0),pos(i,1),pos(i,2),rad(i));
	}

	// check larger box
	c_loop_all clLarger(largerBox);
	vector<double> lBoxVoro(N,0.0);

	// check larger volumes
	cout << "Looping over volumes in larger box" << endl;
	if (clLarger.start()) do if(largerBox.compute_cell(c,clLarger)) {	
		// get residue index
		residueIndex = resID(clLarger.pid());

		// add to larger box volume
		lBoxVoro.at(residueIndex) += c.volume();
	} while (clLarger.inc());

	// do voro++ loop to get volumes and neighbors
	c_loop_all cl(*voroContainer);

	// get initial volumes
	vector<double> atomicvols(NA,-1.0);
	cout << "Looping over initial volumes in normal box" << endl;
	if (cl.start()) do if(voroContainer->compute_cell(c,cl)) {	
		// get residue index
		residueIndex = resID(cl.pid());

		atomicvols[cl.pid()] = c.volume();

		// increment volume on residueIndex
		setVoroVol(residueIndex,resVoro(residueIndex)+c.volume());
	} while (cl.inc());

	for (i=0; i<NA; i++)
		cout << "atomic voro volume : i = " << i << ", avol = " << atomicvols[i] << endl;

	// compute neighbors for each residue
	neighbors();

	// loop over residues, set volumes to - that share at least one voronoi face with the container
	// OR if the volume changes
	for (i=0; i<N; i++){
		if (lBoxVoro.at(i) - resVoro(i) > boxChangeTol)
			setVoroVol(i,-1.0*resVoro(i));
		else{
			for (j=0; j<numResNeighbors(i); j++){
				if (resNeighbors(i,j) == -1){
					setVoroVol(i,-1.0*resVoro(i));
					break;
				}
			}
		}
	}
}

// destructor
proteinPacking::~proteinPacking(){
	if (residueSize){
		delete [] residueSize;
		residueSize = nullptr;
		cout << "deleting residueSize.\n";
	}
	if (sequence){
		delete [] sequence;
		sequence = nullptr;
		cout << "deleting sequence.\n";
	}
	if (atoms){
		delete [] atoms;
		atoms = nullptr;
		cout << "deleting atoms.\n";
	}
	if (x){
		delete [] x;
		x = nullptr;
		cout << "deleting x.\n";
	}
	if (r){
		delete [] r;
		r = nullptr;
		cout << "deleting r.\n";
	}
	if (m){
		delete [] m;
		m = nullptr;
		cout << "deleting m.\n";
	}
	if (voroContainer){
		delete voroContainer;
		voroContainer = nullptr;
		cout << "deleting voroContainer.\n";
	}
	if (voroVols){
		delete [] voroVols;
		voroVols = nullptr;
		cout << "deleting voroVols.\n";
	}
	if (voroNeighbors){
		delete [] voroNeighbors;
		voroNeighbors = nullptr;
		cout << "deleting voroNeighbors.\n";
	}

	// close file object
	packingFile.close();
}




/****************************
	
	Setters

*****************************/

void proteinPacking::setPos(int residue, int atom, int d, double val){
	// index to access
	int index = NDIM*atomID(residue,atom)+d;

	// bounds checking
	if (index >= NDIM*NA){
		cout << "	ERROR: index in setPos goes beyond NA, residue = " << residue << ", atom = " << atom << ", d = " << d << "." << endl;
		exit(1);
	}

	// set value
	x[index] = val;
}

void proteinPacking::setPos(int atom, int d, double val){
	// index to access
	int index = NDIM*atom+d;

	// bounds checking
	if (index >= NDIM*NA){
		cout << "	ERROR: index in setPos goes beyond NA, atom = " << atom << ", d = " << d << ", index = " << index << "." << endl;
		exit(1);
	}

	// set value
	x[index] = val;
}

void proteinPacking::setRad(int residue, int atom, double val){
	// index to access
	int index = atomID(residue,atom);

	// bounds checking
	if (index >= NA){
		cout << "	ERROR: atom index in setRad goes beyond NA, residue = " << residue << ", atom = " << atom << "." << endl;
		exit(1);
	}

	// set value
	r[atom] = val;
}

void proteinPacking::setRad(int atom, double val){
	// bounds checking
	if (atom >= NA){
		cout << "	ERROR: atom index in setRad goes beyond NA, atom = " << atom << "." << endl;
		exit(1);
	}

	// set value
	r[atom] = val;
}

void proteinPacking::setAtom(int atom, string& str){
	// bounds checking
	if (atom >= NA){
		cout << "	ERROR: atom index in setAtom goes beyond NA, atom = " << atom << "." << endl;
		exit(1);
	}

	// set value
	atoms[atom] = str;
}

void proteinPacking::setMass(int residue, double val){
	// bounds checking
	if (residue >= N){
		cout << "	ERROR: residue index in setMass goes beyond N, residue = " << residue << "." << endl;
		exit(1);
	}

	// set value
	m[residue] = val;
}

void proteinPacking::setResidueSize(int residue, int val){
	// bounds checking
	if (residue >= N){
		cout << "	ERROR: residue index in setResidueSize goes beyond N, residue = " << residue << "." << endl;
		exit(1);
	}

	// set value
	residueSize[residue] = val;
}

void proteinPacking::setSequence(int residue, string& str){
	// bounds checking
	if (residue >= N){
		cout << "	ERROR: residue index in setSequence goes beyond N, residue = " << residue << "." << endl;
		exit(1);
	}

	// set value
	sequence[residue] = str;
}

void proteinPacking::setVoroVol(int residue, double val){
	// bounds checking
	if (residue >= N){
		cout << "	ERROR: residue index in setVoroVol goes beyond N, residue = " << residue << "." << endl;
		exit(1);
	}

	// set value
	voroVols[residue] = val;
}





/****************************
	
	Getters

*****************************/

int proteinPacking::cumulativeNumberOfAtoms(int residue){
	int i;
	int val = 0;

	for (i=0; i<residue; i++)
		val += size(i);

	return val;
}

int proteinPacking::atomID(int residue, int atom){
	return cumulativeNumberOfAtoms(residue) + atom;
}

int proteinPacking::resID(int atom){
	// local variables
	int index = 0;
	int residue = 0;

	// if atom is 0, just return 0 for residue
	if (atom == 0)
		return 0;

	// loop over residues, quit when pass atom
	while (index <= atom){
		index += size(residue);
		residue++;
	}

	return residue-1;
}

double proteinPacking::pos(int residue, int atom, int d){
	// index to access
	int index = NDIM*atomID(residue,atom)+d;

	// bounds checking
	if (index >= NA*NDIM){
		cout << "	ERROR: index in pos goes beyond NA, residue = " << residue << ", atom = " << atom << ", d = " << d << "." << endl;
		exit(1);
	}

	// return position
	return x[index];
}

double proteinPacking::pos(int atom, int d){
	// index to access
	int index = NDIM*atom+d;

	// bounds checking
	if (index >= NA*NDIM){
		cout << "	ERROR: index in pos goes beyond NA, atom = " << atom << ", d = " << d << "." << endl;
		exit(1);
	}

	// return position
	return x[index];
}

double proteinPacking::rad(int residue, int atom){
	// index to access
	int index = atomID(residue,atom);

	// bounds checking
	if (index >= NA){
		cout << "	ERROR: atom index in rad goes beyond NA, residue = " << residue << ", atom = " << atom << "." << endl;
		exit(1);
	}

	// return radius
	return r[index];
}

double proteinPacking::rad(int atom){
	// bounds checking
	if (atom >= NA){
		cout << "	ERROR: atom index in rad goes beyond NA, atom = " << atom << "." << endl;
		exit(1);
	}

	return r[atom];
}

double proteinPacking::mass(int residue){
	// bounds checking
	if (residue >= N){
		cout << "	ERROR: residue index in mass goes beyond N, residue = " << residue << "." << endl;
		exit(1);
	}

	return m[residue];
}

int proteinPacking::size(int residue){
	// bounds checking
	if (residue >= N){
		cout << "	ERROR: residue index in size goes beyond N, residue = " << residue << "." << endl;
		exit(1);
	}

	return residueSize[residue];
}

double proteinPacking::resVoro(int residue){
	// bounds checking
	if (residue >= N){
		cout << "	ERROR: residue index in resVoro goes beyond N, residue = " << residue << "." << endl;
		exit(1);
	}

	return voroVols[residue];
}

int proteinPacking::resNeighbors(int residue, int neighbor){
	// bounds checking
	if (residue >= N){
		cout << "	ERROR: residue index in resNeighbors goes beyond N, residue = " << residue << "." << endl;
		exit(1);
	}

	return voroNeighbors[residue].at(neighbor);
}

int proteinPacking::numResNeighbors(int residue){
	// bounds checking
	if (residue >= N){
		cout << "	ERROR: residue index in resNeighbors goes beyond N, residue = " << residue << "." << endl;
		exit(1);
	}

	return voroNeighbors[residue].size();
}

string proteinPacking::getSeq(int residue){
	// bounds checking
	if (residue >= N){
		cout << "	ERROR: residue index in getSeq goes beyond N, residue = " << residue << "." << endl;
		exit(1);
	}

	return sequence[residue];
}

string proteinPacking::getAtom(int residue, int atom){
	// index to access
	int index = atomID(residue,atom);

	// bounds checking
	if (index >= NA){
		cout << "	ERROR: index in getAtom goes beyond NA, residue = " << residue << ", atom = " << atom << "." << endl;
		exit(1);
	}

	return atoms[index];
}

string proteinPacking::getAtom(int atom){
	// bounds checking
	if (atom >= NA){
		cout << "	ERROR: atom index in getAtom goes beyond NA, atom = " << atom << "." << endl;
		exit(1);
	}

	return atoms[atom];
}


/****************************
	
	Calculators

*****************************/

void proteinPacking::calcMasses(int NMCPTS){
	// local variables
	int i,j,p,hit,tmpOctant;
	int xb,yb,zb;
	double xmin=1e6, xmax=0, ymin=1e6, ymax=0, zmin=1e6, zmax=0;
	double dx,dy,dz,da,daMin,boxVolume;
	double rx,ry,rz;
	double mtmp;
	
	// loop over atoms, get min and max
	for (i=0; i<NA; i++){
		if (pos(i,0) < xmin)
			xmin = pos(i,0);
		else if (pos(i,0) > xmax)
			xmax = pos(i,0);

		if (pos(i,1) < ymin)
			ymin = pos(i,1);
		else if (pos(i,1) > ymax)
			ymax = pos(i,1);

		if (pos(i,2) < zmin)
			zmin = pos(i,2);
		else if (pos(i,2) > zmax)
			zmax = pos(i,2);
	}

	// make box a little bit larger (by 1 Angstroms)
	xmin -= 2.0;
	xmax += 2.0;

	ymin -= 2.0;
	ymax += 2.0;

	zmin -= 2.0;
	zmax += 2.0;

	// get box width
	dx = xmax - xmin;
	dy = ymax - ymin; 
	dz = zmax - zmin;
	boxVolume = dx*dy*dz;

	// label atoms based on octants
	vector<int> octantList[8];
	int xo, yo, zo;
	double skinDepth = 3.0;

	for (i=0; i<NA; i++){

		if (pos(i,0) > (0.5*dx + xmin))
			xo = 1;
		else
			xo = 0;

		if (pos(i,1) > (0.5*dy + ymin))
			yo = 1;
		else
			yo = 0;

		if (pos(i,2) > (0.5*dz + zmin))
			zo = 1;
		else
			zo = 0;

		// pick quadrant
		octantList[4*zo + 2*yo + xo].push_back(i);

		// check to see if near x edge
		xb = 0;
		if (abs(pos(i,0) - (0.5*dx + xmin)) < skinDepth)
			xb = 1;

		// check to see if near y edge
		yb = 0;
		if (abs(pos(i,1) - (0.5*dy + ymin)) < skinDepth)
			yb = 1;

		// check to see if near z edge
		zb = 0;
		if (abs(pos(i,2) - (0.5*dz + zmin)) < skinDepth)
			zb = 1;

		// add to other lists depending on boundary status
		if (xb==1 && yb==1 && zb==1){
			// add to every other octant
			octantList[4*(1-zo) + 2*(1-yo) + 1-xo].push_back(i);
			octantList[4*zo + 2*(1-yo) + 1-xo].push_back(i);
			octantList[4*(1-zo) + 2*yo + 1-xo].push_back(i);
			octantList[4*(1-zo) + 2*(1-yo) + xo].push_back(i);
			octantList[4*(1-zo) + 2*yo + xo].push_back(i);
			octantList[4*zo + 2*(1-yo) + xo].push_back(i);
			octantList[4*zo + 2*yo + 1-xo].push_back(i);
		}
		else if (xb==1 && yb==1 && zb==0){
			// add to 3 other octants, flip x&y
			octantList[4*zo + 2*(1-yo) + 1-xo].push_back(i);
			octantList[4*zo + 2*yo + 1-xo].push_back(i);
			octantList[4*zo + 2*(1-yo) + xo].push_back(i);
		}
		else if (xb==1 && yb==0 && zb==1){
			// add to 3 other octants, flip x&z
			octantList[4*(1-zo) + 2*yo + 1-xo].push_back(i);
			octantList[4*zo + 2*yo + 1-xo].push_back(i);
			octantList[4*(1-zo) + 2*yo + xo].push_back(i);
		}
		else if (xb==0 && yb==1 && zb==1){
			// add to 3 other octants, flip y&z
			octantList[4*(1-zo) + 2*(1-yo) + xo].push_back(i);
			octantList[4*zo + 2*(1-yo) + xo].push_back(i);
			octantList[4*(1-zo) + 2*yo + xo].push_back(i);
		}
		else if (xb==1 && yb==0 && zb==0)
			octantList[4*zo + 2*yo + 1-xo].push_back(i);
		else if (xb==0 && yb==1 && zb==0)
			octantList[4*zo + 2*(1-yo) + xo].push_back(i);
		else if (xb==0 && yb==0 && zb==1)
			octantList[4*(1-zo) + 2*yo + xo].push_back(i);
	}

	// vector of hits
	vector<int> resHits(N,0);

	// run MC integration with octants
	for (p=0; p<NMCPTS; p++){
		if (p % 500 == 0)
			cout << "on mc pt p = " << p << endl;

		// minimum distance reset
		daMin = 1e6;
		hit = -1;

		// generate random point
		rx = dx*drand48() + xmin;
		ry = dy*drand48() + ymin;
		rz = dz*drand48() + zmin;

		// decide on quadrant of point
		if (rx > (0.5*dx + xmin))
			xo = 1;
		else
			xo = 0;

		if (ry > (0.5*dy + ymin))
			yo = 1;
		else
			yo = 0;

		if (rz > (0.5*dz + zmin))
			zo = 1;
		else
			zo = 0;

		// get octant
		tmpOctant = 4*zo + 2*yo + xo;
		
		// find minimum distance, only check octants
		for (j=0; j<octantList[tmpOctant].size(); j++){
			// get atom to check
			i = octantList[tmpOctant][j];

			// distance between atom and pt
			da = sqrt(pow(rx-pos(i,0),2) + pow(ry-pos(i,1),2) + pow(rz-pos(i,2),2));

			// check if smallest distance found
			if (da < daMin){
				daMin = da;
				hit = i;
			}
		}

		// check minimum distance, see if inside residue
		if (hit > 0){
			if (daMin < rad(hit)){
				// minimum distance is inside atom, so incremement hit for associated residue
				resHits.at(resID(hit))++;
			}
		}
	}

	// compute masses
	for (i=0; i<N; i++){
		mtmp = ((double)resHits.at(i)/(double)NMCPTS)*boxVolume;
		setMass(i,mtmp);
	}
}

void proteinPacking::neighbors(){
	// local variables
	int i,j,uqFound;
	voro::voronoicell_neighbor c;
	vector<int> neigh;

	// do voro++ loop to get volumes and neighbors
	c_loop_all cl(*voroContainer);
	int residueIndex = 0;

	// get initial volumes
	if (cl.start()) do if(voroContainer->compute_cell(c,cl)) {	
		// get residue index
		residueIndex = resID(cl.pid());

		// check neighbors of atom in cell c
		c.neighbors(neigh);
	
		// loop over each neighbor of cell c, if new neighbor found add to neighbor list	
		for (i=0; i<neigh.size(); i++){
			// check for new/unique neighbors
			uqFound = 1;

			// if neighbor is part of this residue, skip
			if (resID(neigh.at(i)) == residueIndex)
				continue;

			// loop over neighbors already found
			for (j=0; j<numResNeighbors(residueIndex); j++){
				if (resID(neigh.at(i)) == resNeighbors(residueIndex,j)){
					uqFound = 0;
					break;
				}
			}

			// add residue-mapped neigh[i] to neighbor list
			if (uqFound == 1)
				voroNeighbors[residueIndex].push_back(resID(neigh.at(i)));
		}

	} while (cl.inc());
}






/****************************
	
	Printers

*****************************/

void proteinPacking::printPackingFraction(){
	// local variables
	int i;

	// print header
	packingFile << N << endl;

	// print box volumes
	packingFile << setw(20) << getBound(0,0);
	packingFile << setw(20) << getBound(0,1);
	packingFile << endl;

	packingFile << setw(20) << getBound(1,0);
	packingFile << setw(20) << getBound(1,1);
	packingFile << endl;

	packingFile << setw(20) << getBound(2,0);
	packingFile << setw(20) << getBound(2,1);
	packingFile << endl;

	// print packing fractions for each residue
	for (i=0; i<N; i++){
		packingFile << setw(6) << i+1;
		packingFile << setw(10) << getSeq(i);
		packingFile << setw(25) << mass(i);
		packingFile << setw(25) << resVoro(i);

		// if voronoi < 0, then on surface and has 0 packing fraction
		if (resVoro(i) < 0)
			packingFile << setw(25) << 0.0 << endl;
		else
			packingFile << setw(25) << mass(i)/resVoro(i) << endl;
		
	}
}



void proteinPacking::printSingleVoronoiCell(int i){
	// local variables
	int atot, aloc, atest, w1, w2, w3, n, j, k;
	int natmp, nvtmp;
	double x,y,z;
	voro::voronoicell_neighbor c;
	vector<int> f_vert, neigh;
	vector<double> v;

	// check that voronoi file output object is open
	if (!voronoiFile.is_open()){
		cout << "	** ERROR: voronoi file is not open in printVoronoiCells(), ending." << endl;
		exit(1);
	}

	// output widths
	w1 = 6;
	w2 = 10;
	w3 = 20;

	// print initial residue information
	voronoiFile << setw(w1) << i;
	voronoiFile << setw(w2) << getSeq(i);
	voronoiFile << setw(w2) << size(i);
	voronoiFile << endl;

	// loop over atoms, print voronoi cell information 
	// i.e. vertex positions and which vertex belongs to which face
	c_loop_all cl(*voroContainer);
	if (cl.start()) do if(voroContainer->compute_cell(c,cl)) {
		// get atomic id
		atot = cl.pid();

		// skip if atom id does not coincide with resid i
		if (resID(atot) != i)
			continue;

		// get local label of atom
		aloc = atot - cumulativeNumberOfAtoms(i);

		// check that labelling is done correctly
		atest = atomID(i,aloc);
		if (atest != atot){
			cout << "	** ERROR: atom labelling done incorrectly in printSingleVoronoiCell(), ending. " << endl;
			exit(1);
		}
		else{
			// cout << "-- on residue " << i << ", atomid = " << aloc << endl;
			// cout << "	** resid = " << resID(atot) << endl;
			// cout << "	** cumsum = " << cumulativeNumberOfAtoms(i) << endl;
			// cout << "	** atot = " << atot << endl;
		}

		// print atomic information corresponding to this residue
		voronoiFile << setw(w1) << aloc;
		voronoiFile << setw(w2) << getAtom(atot);
		voronoiFile << setw(w2) << rad(atot);
		voronoiFile << setw(w3) << pos(atot,0);
		voronoiFile << setw(w3) << pos(atot,1);
		voronoiFile << setw(w3) << pos(atot,2);
		voronoiFile << endl;

		// get face vertex information
		cl.pos(x,y,z);
		c.neighbors(neigh);
		c.face_vertices(f_vert);
		c.vertices(x,y,z,v);

		// loop over vertices, print coordinates to file
		nvtmp = v.size()/3;
		voronoiFile << setw(w1) << nvtmp << endl;
		for (j=0; j<nvtmp; j++){
			voronoiFile << setw(w1) << j;
			voronoiFile << setw(w3) << v[3*j];
			voronoiFile << setw(w3) << v[3*j + 1];
			voronoiFile << setw(w3) << v[3*j + 2];
			voronoiFile << endl;
		}

		// Loop over faces, print which vertices map to which face
		// 	NOTE: face info stored as (f, v1, v2, ..., vf), where f is number of vertices
		// 	per face, and v1, ..., vf are the f vertex labels corresponding to this face

		voronoiFile << setw(w1) << neigh.size() << endl;
		j=0;
		for (n=0; n<neigh.size(); n++){
			// output face index and number of vertices / face
			voronoiFile << setw(w1) << f_vert[j];

			// loop over vertices on face, print indices
			for (k=0; k<f_vert[j]; k++)
				voronoiFile << setw(w2) << f_vert[j+k+1];

			// go to next face
			j += f_vert[j]+1;

			// print new line
			voronoiFile << endl;
		}
	} while (cl.inc());
}





















