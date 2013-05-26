
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <new>

#include "CArena.H"
#include "REAL.H"
#include "Utility.H"
#include "IntVect.H"
#include "Box.H"
#include "Geometry.H"
#include "ParmParse.H"
#include "ParallelDescriptor.H"
#include "VisMF.H"
#include "ArrayLim.H"

#include "Amr.H"
#include "AmrLevel.H"

#include "zzz.hpp"

//#include "writePlotFile.H"

//
#if (BL_SPACEDIM != 2)
#	error "Only BL_SPACEDIM=2 is allowed!"
#endif


#define USE_FORT_FUNCTION

#if defined(USE_FORT_FUNCTION)
#	if defined(BL_FORT_USE_UNDERSCORE)
#		define FORT_INIT_PHI init_phi_
#		define FORT_COMPUTE_FLUX compute_flux_
#		define FORT_UPDATE_PHI update_phi_
#	else
#		error "Use BL_FORT_USE_UNDERSCORE or die!"
#	endif

extern "C" {

void FORT_INIT_PHI(Real* data, const int *low, const int *high,
		const int *nghost, const Real *dx,
		const Real *problem_low, const Real *problem_high);

void FORT_COMPUTE_FLUX(Real *phi, const int *nghostPhi,
		Real *fluxx, Real *fluxy, const int *nghostFlux,
		const int *low, const int *high, const Real *dx);

void FORT_UPDATE_PHI(Real *phiOld, Real *phinew, const int *nghostPhi,
		Real *fluxx, Real *fluxy, const int *nghostFlux,
		const int *low, const int *high, const Real *dx, const Real *dt);
}
#endif


static void init_phi(FArrayBox &phi, const int nghost,
		const Box &validBox, const Real *dx, const Geometry &geom)
{
	const int *indexLow = validBox.loVect();
	const int *indexHigh = validBox.hiVect();

	const double *probLow = geom.ProbLo();
	const double *probHigh = geom.ProbHi();

	for(int j=indexLow[1]; j<=indexHigh[1]; j++) {
		double y = probLow[1] + (j+0.5) * dx[1];

		for(int i=indexLow[0]; i<=indexHigh[0]; i++) {
			double x = probLow[0] + (i+0.5) * dx[0];
			double r2 = pow(x-0.25, 2) + pow(y-0.25, 2);
			IntVect idx(i, j);
			phi(idx) = 1.0 + exp(-r2 * 100);
		}
	}

}


//	FORT_COMPUTE_FLUX(data, &nGhostPhi,
//		flux[0][mfi].dataPtr(), flux[1][mfi].dataPtr(),
//		&nGhostFlux, box.loVect(), box.hiVect(), dx);
static void compute_flux(FArrayBox &phi, const int nghostPhi,
		FArrayBox &fluxx, FArrayBox &fluxy, const int nghostFlux,
		const Box &validBox, const Real *dx)
{
	const int *indexLow = validBox.loVect();
	const int *indexHigh = validBox.hiVect();

	// x-flux
	for(int j=indexLow[1]; j<=indexHigh[1]; j++) {
		for(int i=indexLow[0]; i<=indexHigh[0]+1; i++) {
			IntVect idx(i,j);
			IntVect left(i-1,j), right(i,j);
			fluxx(idx) = (phi(right) - phi(left)) / dx[0];
		}
	}

	// y-flux
	for(int j=indexLow[1]; j<=indexHigh[1]+1; j++) {
		for(int i=indexLow[0]; i<=indexHigh[0]; i++) {
			IntVect idx(i,j);
			IntVect up(i,j), down(i,j-1);
			fluxy(idx) = (phi(up) - phi(down)) / dx[1];
		}
	}
}


static Real compute_dt(const Real *dx) {
	double len = dx[0];
	return 0.9 * len*len / (2.0*BL_SPACEDIM);
}

static void advance(MultiFab *oldPhi, MultiFab *newPhi, MultiFab *flux,
		const Real *dx, const Real dt, const Geometry &geom)
{
	// fill ghost cells covered by valid cells
	oldPhi->FillBoundary();
	// fill periodic boundary ghost cells
	geom.FillPeriodicBoundary(*oldPhi);

	int nComp = oldPhi->nComp();
	int nGhostPhi = oldPhi->nGrow();
	int nGhostFlux = flux->nGrow();

	// compute fluxes one grid at a time
	for(MFIter mfi(*oldPhi); mfi.isValid(); ++mfi) {
		const Box &box = mfi.validbox();

		Real *data = (*oldPhi)[mfi].dataPtr();
#ifdef USE_FORT_FUNCTION
		FORT_COMPUTE_FLUX(data, &nGhostPhi,
				flux[0][mfi].dataPtr(), flux[1][mfi].dataPtr(),
				&nGhostFlux, box.loVect(), box.hiVect(), dx);
#else
		compute_flux((*oldPhi)[mfi], nGhostPhi,
				flux[0][mfi], flux[1][mfi], nGhostFlux,
				box, dx);
#endif
	}

	// advance the solution one grid at a time
	for(MFIter mfi(*oldPhi); mfi.isValid(); ++mfi) {
		const Box &box = mfi.validbox();

		Real *oldData = (*oldPhi)[mfi].dataPtr();
		Real *newData = (*newPhi)[mfi].dataPtr();

		FORT_UPDATE_PHI(oldData, newData, &nGhostPhi,
				flux[0][mfi].dataPtr(), flux[1][mfi].dataPtr(), &nGhostFlux,
				box.loVect(), box.hiVect(), dx, &dt);
	}
}

//int main(int argc, char **argv) {
//
//	BoxLib::Initialize(argc, argv);
//
////	Real start_time = ParallelDescriptor::second();
//
//	//
//	ParmParse pp;
//
//	int n_cell = -1;
//	pp.get("n_cell", n_cell);
//
//	int max_grid_size = -1;
//	pp.get("max_grid_size", max_grid_size);
//
//	int plot_int = -1;
//	pp.get("plot_int", plot_int);
//
//	int max_step = -1;
//	pp.get("max_step", max_step);
//
////	Real start_time = -1, stop_time = -1;
////	pp.query("start_time", start_time);
////	pp.query("stop_time", stop_time);
//	Real start_time = 0.0;
//
//	// TODO report
//	std::cout << "n_cell=" << n_cell << ";"
//			<< "max_grid_size=" << max_grid_size << ";"
//			<< "plot_int=" << plot_int << ";"
//			<< "max_step=" << max_step << ";"
//			<< std::endl;
//
//	// integer indexing, both inclusive
//	IntVect domLow(0,0);
//	IntVect domHigh(n_cell-1,n_cell-1);
//	//
//	Box domain(domLow, domHigh);
//	// break the domain into chunks
//	BoxArray bs(domain);
//	bs.maxSize(max_grid_size);
//
//	// physical size
//	RealBox realbox;
//	for(int i=0; i<BL_SPACEDIM; i++) {
//		realbox.setLo(i, -1.0);
//		realbox.setHi(i, 1.0);
//	}
//
//	// use Cartesian coordinates
//	CoordSys::CoordType coord = CoordSys::cartesian;
//	// periodic BC
//	int isPeriodic[BL_SPACEDIM];
//	for(int i=0; i<BL_SPACEDIM; i++) {
//		isPeriodic[i] = 1;
//	}
//
//	// geometry object
//	Geometry geom(domain, &realbox, coord, isPeriodic);
//
//	// mesh spacing
//	Real dx[BL_SPACEDIM];
//	for(int i=0; i<BL_SPACEDIM; i++) {
//		int ni = domain.length(i);
//		dx[i] = (geom.ProbHi(i) - geom.ProbLo(i)) / ni;
//	}
//
//	//
//	const int nGhost = 1;
//	const int nComp = 1;
//	if(nGhost > max_grid_size) {
//		BoxLib::Error("nGhost>max_grid_size, grids are too small.");
//	}
//
////	MultiFab oldPhi(bs, nComp, nGhost);
////	MultiFab newPhi(bs, nComp, nGhost);
////	oldPhi.setVal(0);
////	newPhi.setVal(0);
//	MultiFab *oldPhi = new MultiFab(bs, nComp, nGhost);
//	MultiFab *newPhi = new MultiFab(bs, nComp, nGhost);
//	oldPhi->setVal(0);
//	newPhi->setVal(0);
//
//	// initialize using MultiFab iterator
//	for(MFIter mfi(*newPhi); mfi.isValid(); ++mfi) {
//		const Box &box = mfi.validbox();
//#ifdef USE_FORT_FUNCTION
//		FORT_INIT_PHI(
//				(*newPhi)[mfi].dataPtr(),
//				box.loVect(), box.hiVect(),
//				&nGhost, dx, geom.ProbLo(), geom.ProbHi());
//#else
//		init_phi((*newPhi)[mfi], nGhost, box, dx, geom);
//#endif
//	}
//	std::cout << "Initialized FABs" << std::endl;
//
//	int step = 0;
//	Real dt = compute_dt(dx);
//
//	// write initial output
//	if(plot_int > 0) {
//		std::string pltdir = BoxLib::Concatenate("plt", step, 5);
//		writePlotFile(pltdir, *newPhi, geom);
//	}
//
//	// build flux MultiFabs, staggered on each direction
//	MultiFab *flux = new MultiFab[BL_SPACEDIM];
//	for(int i=0; i<BL_SPACEDIM; i++) {
//		// flux_i has one component, zero ghost cells, nodal in its direction
//		BoxArray edgeGrids(bs);
//		// change to NODE
//		edgeGrids.surroundingNodes(i);
//
//		int nvar = 1, nghost = 0;
//		flux[i].define(edgeGrids, nvar, nghost, Fab_allocate);
//	}
//
//	 // time loop
//	for(int n=1; n<=max_step; n++) {
//		// swap data
//		std::swap(oldPhi, newPhi);
//
//		// new = old + dt * (div flux)
//		advance(oldPhi, newPhi, flux, dx, dt, geom);
//
//		if(ParallelDescriptor::IOProcessor()) {
//			std::cout << "Step " << n << std::endl;
//		}
//
//		// plot file
//		if(plot_int>0 && n%plot_int==0) {
//			std::string pltdir = BoxLib::Concatenate("plt", n, 5);
//			writePlotFile(pltdir, *newPhi, geom);
//		}
//	}
//
////	Real stop_time = ParallelDescriptor::second() - start_time;
////	const int ioProcNumber = ParallelDescriptor::IOProcessorNumber();
////	ParallelDescriptor::ReduceRealMax(stop_time, ioProcNumber);
////	if(ParallelDescriptor::IOProcessor()) {
////		std::cout << "Run time = " << stop_time << std::endl;
////	}
//
//
//
//#	if (0)
//	MultiFab phi(bs, nComp, nGhost);
//	phi.setVal(3.14);
//
//	for(MFIter mfi(phi); mfi.isValid(); ++mfi) {
//		const Box &box = mfi.validbox();
//		FORT_INIT_PHI(phi[mfi].dataPtr(),
//				box.loVect(), box.hiVect(), &nGhost,
//				dx, geom.ProbLo(), geom.ProbHi());
//		std::cout << "Box" << mfi.index() << std::endl;
//	}
//
//	int boxCount = 0;
//	for(MFIter mfi(phi); mfi.isValid(); ++mfi) {
//		const Box& validBox = mfi.validbox();
//		const Box& fabBox = mfi.fabbox();
//
//		std::cout << "Box" << mfi.index() << std::endl;
//		std::cout << "validbox: " << validBox << std::endl;
//		std::cout << "fabbox" << fabBox << std::endl;
//
//		boxCount++;
//	}
//	// initial output
//	if(plot_int > 0) {
//		int n = 0;
//		std::string pltdir = BoxLib::Concatenate("plt", n, 5);
//		writePlotFile(pltdir, phi, geom);
//		std::cout << pltdir << std::endl;
//	}
//#endif
//
//	BoxLib::Finalize();
//
//	return 0;
//}
//

int main(int argc, char **argv) {

	BoxLib::Initialize(argc, argv);

//	Real start_time = ParallelDescriptor::second();

	//
	ParmParse pp;

	int max_step = -1;
	pp.get("max_step", max_step);

	Real start_time = 0.0;
	Real stop_time = -1;
	pp.query("start_time", start_time);
	pp.query("stop_time", stop_time);

	if(ParallelDescriptor::IOProcessor()) {
		std::cout << "max_step=" << max_step << ";"
				<< "stop_time=" << stop_time << ";"
				<< std::endl;
	}

	ParallelDescriptor::Barrier();

	// build AMR object
	Amr *pAmr = new Amr;
	pAmr->init(start_time, stop_time);

	// TODO regrid_on_restart

	// time loop
	while(pAmr->okToContinue() &&
			(pAmr->levelSteps(0)<max_step || max_step<0) &&
			(pAmr->cumTime()<stop_time || stop_time<0))
	{
		pAmr->coarseTimeStep(stop_time);
	}

	delete pAmr;

	// TODO report execute time
	// TODO report CArena memory usage



	BoxLib::Finalize();

	return 0;
}


