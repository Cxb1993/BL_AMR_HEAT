/*
 * zzz.cpp
 *
 *  Created on: 2013/04/08
 *      Author: homu
 */

#include <memory>

#include "LevelBld.H"
#include "winstd.H"
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
#include "BLFort.H"

#include "TagBox.H"
#include "MacBndry.H"
#include "MultiGrid.H"
#include "CGSolver.H"
#include "ABecLaplacian.H"




#include "zzz.hpp"
#include "fortranFunctions.h"
#include "FortranHelper.hpp"


/*
 *
 */
int ZZZAmr::verbose = 1;
Real ZZZAmr::cfl = 0.9;

BCRec ZZZAmr::physBC;

int ZZZAmr::allow_untagging = 0;
ErrorList ZZZAmr::errList;

int ZZZAmr::doReflux = 1;

int ZZZAmr::doImplicitStep = 0;

int ZZZAmr::numState = 1;
int ZZZAmr::numGhost = 1;

int ZZZAmr::StatePhi = -1;

int ZZZAmr::radiusGrow = 1;

/*
 *
 */


/**
 * interior, in-flow, out-flow, symmetry, slip-wall, non-slip-wall
 */
static const int bcmap_scalar[] = {
		INT_DIR, // interior
		EXT_DIR, // in-flow,
		FOEXTRAP, // out-flow
		REFLECT_EVEN, // symmetry
		REFLECT_EVEN, // slip-wall 
		REFLECT_EVEN, // non-slip-wall
};

/**
 * convert BC type for scalars
 */
static void set_scalar_bc(BCRec &mathBC, const BCRec &physBC) {
	const int *lo_bc = physBC.lo();
	for(int i=0; i<BL_SPACEDIM; i++) {
		int low = bcmap_scalar[physBC.lo(i)];
		mathBC.setLo(i, low);

		int high = bcmap_scalar[physBC.hi(i)];
		mathBC.setHi(i, high);
	}
}

static void bndryfunc_fill_phi(
		Real *phi, const int &phiL1, const int &phiL2, const int &phiH1, const int &phiH2,
		const int dlo[], const int dhi[],
		const Real dx[], const Real glo[],
		const Real *time, const int bc[])
{
	BL_FORT_PROC_CALL(FILCC, filcc) (
			phi, phiL1, phiL2, phiH1, phiH2,
			dlo, dhi, dx, glo, bc
	);
}

static void derivefunc_calc_temp(
		Real* tempData, const int &t_l1, const int &t_l2,
		const int &t_h1, const int &t_h2, const int* nvar,
		const Real* stateData, const int &s_l1, const int &s_l2,
		const int &s_h1, const int &s_h2, const int* ncomp,
		const int* lo, const int* hi,
		const int* domain_lo, const int* domain_hi,
		const Real* delta, const Real* xlo,
		const Real* time, const Real* dt,
		const int* bcrec,
		const int* level, const int* grid_no)
{
	// data view
	FortArrayWrap<Real,3> temp(tempData, t_l1,t_h1, t_l2,t_h2, 0,*nvar-1);
	FortArrayWrap<const Real,3> state((stateData), s_l1,s_h1, s_l2,s_h2, 0,*ncomp-1);

	const int IPhi = ZZZAmr::StatePhi;

	for(int j=lo[1]; j<=hi[1]; j++) {
		for(int i=lo[0]; i<=hi[0]; i++) {
			temp(i,j,0) = state(i,j,IPhi) - 1;
		}
	}
}

void ZZZAmr::variableSetUp() {
	if(desc_lst.size() != 0) {
		BoxLib::Error(__FUNCTION__);
	}

	/*
	 * read options, physical BC, etc.
	 */
	readParams();

	/*
	 * setup state variables
	 */
	// set indices for state variables
	StatePhi = 0;

	numState = 1;
	numGhost = 1;

	// add descriptors for state variables
	Interpolater *interp = &cell_cons_interp;
	desc_lst.addDescriptor(State_Type, IndexType::TheCellType(),
			StateDescriptor::Point, numGhost, numState, interp);

	// BC for state variables
	Array<BCRec> bcs(numState);
	Array<std::string> names(numState);
	BCRec bc;

	// scalar quantity PHI
	set_scalar_bc(bc, physBC);
	bcs[StatePhi] = bc; names[StatePhi] = "phi";

	desc_lst.setComponent(State_Type, StatePhi, names[StatePhi], bcs[StatePhi],
			StateDescriptor::BndryFunc(bndryfunc_fill_phi),
//			StateDescriptor::BndryFunc(BL_FORT_PROC_NAME(PHI_FILL, phi_fill)),
			interp);

	/*
	 * add derived variables
	 */
	derive_lst.add("temperature", IndexType::TheCellType(), 1,
			derivefunc_calc_temp, DeriveRec::TheSameBox);
	derive_lst.addComponent("temperature", desc_lst, State_Type, 0, numState);

	ErrorSetUp();
}

void ZZZAmr::variableCleanUp() {
	desc_lst.clear();
}

void ZZZAmr::readParams() {
	static bool done = false;
	if(done) return;

	done = true;

	ParmParse pp("zzz");

	verbose = 1;
	pp.query("verbose", verbose);

	cfl = 0.9;
	pp.get("cfl", cfl);

	doReflux = 1;
	pp.get("doReflux", doReflux);

	doImplicitStep = 0;
	pp.get("doImplicitStep", doImplicitStep);

	// get BC
	Array<int> lowBC(BL_SPACEDIM), highBC(BL_SPACEDIM);
	pp.getarr("lo_bc", lowBC, 0, BL_SPACEDIM);
	pp.getarr("hi_bc", highBC, 0, BL_SPACEDIM);
	for(int i=0; i<BL_SPACEDIM; i++) {
		physBC.setLo(i, lowBC[i]);
		physBC.setHi(i, highBC[i]);
	}

	// check BC consistency
	checkBoundaryConsistency();

	// dump inputs
}

void ZZZAmr::checkBoundaryConsistency() {
	const int *lo_bc = physBC.lo();
	const int *hi_bc = physBC.hi();

	// periodic
	if (Geometry::isAnyPeriodic()) {
		//
		// Do idiot check.  Periodic means interior in those directions.
		//
		for (int dir = 0; dir<BL_SPACEDIM; dir++) {
			if (Geometry::isPeriodic(dir)) {
				if (lo_bc[dir] != Interior) {
					std::cerr << "Castro::read_params:periodic in direction "
							<< dir
							<< " but low BC is not Interior\n";
					BoxLib::Error();
				}
				if (hi_bc[dir] != Interior) {
					std::cerr << "Castro::read_params:periodic in direction "
							<< dir
							<< " but high BC is not Interior\n";
					BoxLib::Error();
				}
			}
		}
	} else {
		//
		// Do idiot check.  If not periodic, should be no interior.
		//
		for (int dir=0; dir<BL_SPACEDIM; dir++) {
			if (lo_bc[dir] == Interior) {
				std::cerr << "Castro::read_params:interior bc in direction "
						<< dir
						<< " but not periodic\n";
				BoxLib::Error();
			}
			if (hi_bc[dir] == Interior) {
				std::cerr << "Castro::read_params:interior bc in direction "
						<< dir
						<< " but not periodic\n";
				BoxLib::Error();
			}
		}
	}

	if ( Geometry::IsRZ() && (lo_bc[0] != Symmetry) ) {
		std::cerr << "ERROR:Castro::read_params: must set r=0 boundary condition to Symmetry for r-z\n";
		BoxLib::Error();
	}

#if (BL_SPACEDIM == 1)
	if ( Geometry::IsSPHERICAL() ) {
		if ( (lo_bc[0] != Symmetry) && (Geometry::ProbLo(0) == 0.0) ) {
			std::cerr << "ERROR:Castro::read_params: must set r=0 boundary condition to Symmetry for spherical\n";
			BoxLib::Error();
		}
	}
#elif (BL_SPACEDIM == 2)
	if ( Geometry::IsSPHERICAL() ) {
		BoxLib::Abort("We don't support spherical coordinate systems in 2D");
	}
#elif (BL_SPACEDIM == 3)
	if ( Geometry::IsRZ() ) {
		BoxLib::Abort("We don't support cylindrical coordinate systems in 3D");
	} else if ( Geometry::IsSPHERICAL() ) {
		BoxLib::Abort("We don't support spherical coordinate systems in 3D");
	}
#endif

}


/**
 * @param level: the refinement level of this array, in C notation (from 0)
 */
static void errorfunc_thres_phi(
		int* tagData, const int &tagl1, const int &tagl2, const int &tagh1, const int &tagh2,
		const int* set, const int* clear,
		Real* phiData, const int &phil1, const int &phil2, const int &phih1, const int &phih2,
		const int* lo, const int * hi, const int* nvar,
		const int* domain_lo, const int* domain_hi,
		const Real* dx, const Real* xlo,
		const Real* prob_lo, const Real* time,
		const int* level)
{
	const int tagVal = *set;
	const int clearVal = *clear;
	const int nvarVal = *nvar;
	const int levelVal = *level;

	// data view
	FortArrayWrap<int,2> tag(tagData, tagl1,tagh1, tagl2,tagh2);
	FortArrayWrap<Real,3> phi(phiData, phil1,phih1, phil2,phih2, 0,nvarVal-1);

	// tag regions of high PHI value
	Real phi_err_thres;
	switch(levelVal) {
	case 0: phi_err_thres = 1.01; break;
	case 1: phi_err_thres = 1.05; break;
	case 2: phi_err_thres = 1.1; break;
	default: phi_err_thres = 1.5; break;
	}

	const int max_phi_err_lev = 4;
	if(levelVal < max_phi_err_lev) {
		for(int j=lo[1]; j<=hi[1]; j++) {
			for(int i=lo[0]; i<=hi[0]; i++) {
				if(phi(i,j,0) >= phi_err_thres) {
					tag(i,j) = tagVal;
				}
			}
		}
	}
}

/**
 * define error estimation quantities
 */
void ZZZAmr::ErrorSetUp() {

	// TODO e.g., error estimation based on Laplacian

	/*
	 * error estimation based on threshold
	 */
//	std::string statePhiName = "phi";
	std::string statePhiName = desc_lst[State_Type].name(StatePhi);
	errList.add(statePhiName, 1, ErrorRec::Special,
			ErrorRec::ErrorFunc(errorfunc_thres_phi));
}



////////////////////////////////////////////////////////////////
// constructors
////////////////////////////////////////////////////////////////

ZZZAmr::ZZZAmr() {
	// nothing
	fluxReg = NULL;
}

ZZZAmr::ZZZAmr(Amr &papa, int level, const Geometry &levelGeom,
		const BoxArray &boxList, Real time)
: AmrLevel(papa, level, levelGeom, boxList, time)
{
	buildMetrics();

	fluxReg = NULL;

	if(level>0 && doReflux) {
		fluxReg = new FluxRegister(grids, crse_ratio, level, numState);
		fluxReg->setVal(0.0);
	}

	// initialize other physics, e.g. gravity, diffusion, reaction, etc.
}

ZZZAmr::~ZZZAmr() {
	if(fluxReg) {
		delete fluxReg;
		fluxReg = NULL;
	}
}

/**
 * build volume, face area, etc.
 */
void ZZZAmr::buildMetrics() {
	// radius, for cylindrical coordinate
	const int numGrids = grids.size();
	radius.resize(numGrids);

	const Real *dx = geom.CellSize();
	const Real *probLo = geom.ProbLo();

	for(int i=0; i<numGrids; i++) {
		const Box &box = grids[i];
		int ilo = box.smallEnd(0) - radiusGrow;
		int ihi = box.bigEnd(0) + radiusGrow;
		int len = ihi - ilo + 1;

		radius[i].resize(len);

		Real *rad = radius[i].dataPtr();

		if(Geometry::IsCartesian()) {
			for(int j=0; j<len; j++) {
				rad[j] = 1.0;
			}
		} else {
			// cylindrical, around y-axis
			RealBox gridloc = RealBox(box, dx, probLo);
			const Real xlo = gridloc.lo(0) - radiusGrow*dx[0];
			for(int j=0; j<len; j++) {
				rad[j] = xlo + (j+0.5) * dx[0];
			}
		}
	}

	volume.clear();
	for(int i=0; i<BL_SPACEDIM; i++) {
		area[i].clear();
	}
	dLogArea[0].clear();

	// volume
	geom.GetVolume(volume, grids, numGhost);

	// faces area
	for(int i=0; i<BL_SPACEDIM; i++) {
		geom.GetFaceArea(area[i], grids, i, numGhost);
	}

	// dLogArea, for 1D/2D only
#if (BL_SPACEDIM <= 2)
	geom.GetDLogA(dLogArea[0], grids, 0, numGhost);
#endif
}

////////////////////////////////////////////////////////////////
// output and plotting
////////////////////////////////////////////////////////////////

std::string ZZZAmr::thePlotFileType() const {
	static const std::string plotFileType("NavierStokes-V1.1");
	return plotFileType;
}


void ZZZAmr::writePlotFile(const std::string &dir, std::ostream &os, VisMF::How how) {
	// list of indices of state variable to plot
	// pair<0> is state_type, pair<1> is component
	std::vector<std::pair<int,int> > plot_var_map;
	for(int k=0; k<desc_lst.size(); k++) {
		for(int comp=0; comp<desc_lst[k].nComp(); comp++) {
			if(parent->isStatePlotVar(desc_lst[k].name(comp)) &&
					desc_lst[k].getType()==IndexType::TheCellType()) {
				plot_var_map.push_back(std::make_pair(k, comp));
			}
		}
	}

	// derived variables
	int num_derive = 0;
	std::list<std::string> derive_names;
	const std::list<DeriveRec> &derivList = derive_lst.dlist();

	for(std::list<DeriveRec>::const_iterator it=derivList.begin();
			it != derivList.end(); ++it)
	{
		if(parent->isDerivePlotVar(it->name())) {
			derive_names.push_back(it->name());
			num_derive++;
		}
	}

	// all variables to plot
	int num_data_items = plot_var_map.size() + num_derive;

	Real curTime = state[State_Type].curTime();

	// Only let 64 CPUs be writing at any one time.
//	VisMF::SetNOutFiles(64);

	// Only the IOProcessor() writes to the header file.
	if (level==0 && ParallelDescriptor::IOProcessor()) {
		// plot file type
		os << thePlotFileType() << '\n';

		if(num_data_items == 0) {
			BoxLib::Error(__FUNCTION__);
		}
		os << num_data_items << '\n';

		// plot variables
		for(int i=0; i<plot_var_map.size(); i++) {
			int type = plot_var_map[i].first;
			int comp = plot_var_map[i].second;
			os << desc_lst[type].name(comp) << '\n';
		}
		for(std::list<std::string>::const_iterator it=derive_names.begin();
				it!=derive_names.end(); ++it) {
			const DeriveRec *rec = derive_lst.get(*it);
			os << rec->variableName(0) << '\n';
		}

		const int finest_lev = parent->finestLevel();

		os << BL_SPACEDIM << '\n';
		os << parent->cumTime() << '\n';
		os << finest_lev << '\n';
		for(int i = 0; i < BL_SPACEDIM; i++) {
			os << geom.ProbLo(i) << ' ';
		}
		os << '\n';
		for (int i = 0; i < BL_SPACEDIM; i++) {
			os << geom.ProbHi(i) << ' ';
		}
		os << '\n';

		for(int i=0; i<finest_lev; i++) {
			os << parent->refRatio(i)[0] << ' ';
		}
		os << '\n';
		for(int i=0; i<=finest_lev; i++) {
			os << parent->Geom(i).Domain() << ' ';
		}
		os << '\n';
		for(int i=0; i<=finest_lev; i++) {
			os << parent->levelSteps(i) << ' ';
		}
		os << '\n';
		for(int i=0; i<=finest_lev; i++) {
			for (int k = 0; k < BL_SPACEDIM; k++) {
				os << parent->Geom(i).CellSize(k) << ' ';
			}
			os << '\n';
		}

		os << (int) Geometry::Coord() << '\n';
		os << "0\n"; // write boundary data
	}

	// Build the directory to hold the MultiFab at this level.
	// The name is relative to the directory containing the Header file.
	static const std::string BaseName = "/Cell";
	std::string Level = BoxLib::Concatenate("Level_", level, 1);

	// Now for the full pathname of that directory.
	std::string FullPath = dir;
	if (!FullPath.empty() && FullPath[FullPath.length()-1] != '/')
		FullPath += '/';
	FullPath += Level;

	// Only the I/O processor makes the directory if it doesn't already exist.
	if (ParallelDescriptor::IOProcessor()) {
		if (!BoxLib::UtilCreateDirectory(FullPath, 0755)) {
			BoxLib::CreateDirectoryFailed(FullPath);
		}
	}

	// Force other processors to wait till directory is built.
	ParallelDescriptor::Barrier();

	if (ParallelDescriptor::IOProcessor()) {

		os << level << ' ' << grids.size() << ' ' << curTime << '\n';
		os << parent->levelSteps(level) << '\n';

		for (int i = 0; i < grids.size(); ++i) {
			RealBox loc = RealBox(grids[i], geom.CellSize(), geom.ProbLo());
			for (int n = 0; n < BL_SPACEDIM; n++) {
				os << loc.lo(n) << ' ' << loc.hi(n) << '\n';
			}
		}

		/*
		 * The full relative pathname of the MultiFabs at this level.
		 * The name is relative to the Header file containing this name.
		 * It's the name that gets written into the Header.
		 */
		if(num_data_items > 0) {
			std::string PathNameInHeader = Level;
			PathNameInHeader += BaseName;
			os << PathNameInHeader << '\n';
		}
	}

	/*
	 * Combine all the MultiFabs (state, derived, etc.) into PlotMF.
	 * NOTE: assuming that each state variable has one component,
	 * but a derived variable is allowed to have multiple components.
	 */
	int cnt = 0;
	const int nonGrow = 0;
	MultiFab plotMF(grids, num_data_items, nonGrow);

	// copy data from state variables, no ghost cells.
	for(int i=0; i<plot_var_map.size(); i++) {
		int type = plot_var_map[i].first;
		int comp = plot_var_map[i].second;

		const MultiFab &this_data = state[type].newData();
		MultiFab::Copy(plotMF, this_data, comp, cnt, 1, nonGrow);

		cnt++;
	}

	// copy data from derived variables
	if(derive_names.size() > 0) {
		for(std::list<std::string>::const_iterator it = derive_names.begin();
				it != derive_names.end(); ++it)
		{
			MultiFab *derive_data = derive(*it, curTime, nonGrow);
			MultiFab::Copy(plotMF, *derive_data, 0, cnt, 1, nonGrow);
			delete derive_data; // this must be done by hand

			cnt++;
		}
	}

	// Use the Full pathname when naming the MultiFab.
	std::string TheFullPath = FullPath;
	TheFullPath += BaseName;

	VisMF::Write(plotMF,TheFullPath);
}

////////////////////////////////////////////////////////////////
// time stepping
////////////////////////////////////////////////////////////////

Real ZZZAmr::initialTimeStep() {
	Real dummy_dt = 0.0;
	Real dt = estimateTimeStep(dummy_dt);

	return dt;
}
Real ZZZAmr::estimateTimeStep(const Real dt_old) {
	const Real *dx = geom.CellSize();
	Real len = dx[0];
	Real dt = cfl * len*len / (2.0*BL_SPACEDIM);

	if(verbose && ParallelDescriptor::IOProcessor()) {
		std::cout << __FUNCTION__
				<< ":level=" << level << ";"
				<< "dt=" << dt << std::endl;
	}

	return dt;
}

void ZZZAmr::computeInitialDt(int finest_level, int sub_cycle,
		Array<int> &n_cycle, const Array<IntVect> &ref_ratio,
		Array<double> &dt_level, Real stop_time)
{
	// grids have been constructed, compute DT for all levels
	if(level > 0) return;

	n_cycle[0] = 1;
	for(int i=1; i<=finest_level; i++) {
		n_cycle[i] = sub_cycle ? parent->MaxRefRatio(i-1) : 1;
	}

	Real dt0 = 1e100;
	int nFactor = 1;
	for(int i=0; i<=finest_level; i++) {
		dt_level[i] = getLevel(i).estimateTimeStep(-1);
		nFactor *= n_cycle[i];
		dt0 = std::min(dt0, nFactor*dt_level[i]);
	}

	// limit DT by the value of stop_time
	const Real eps = 0.0001 * dt0;
	Real cur_time = state[State_Type].curTime();
	if(stop_time >= 0) {
		if(cur_time+dt0 > stop_time-eps) {
			dt0 = stop_time - cur_time;
		}
	}

	nFactor = 1;
	for(int i=0; i<=finest_level; i++) {
		nFactor *= n_cycle[i];
		dt_level[i] = dt0 / nFactor;
	}
}

void ZZZAmr::computeNewDt(int finest_level, int sub_cycle,
		Array<int>& n_cycle, const Array<IntVect>& ref_ratio,
		Array<Real>& dt_min, Array<Real>& dt_level,
		Real stop_time, int post_regrid_flag)
{
	// we are at the end of a coarse grid time cycle.
	// compute the time steps for the next iteration.
	if(level > 0) return;

	n_cycle[0] = 1;
	for(int i=1; i<=finest_level; i++) {
		n_cycle[i] = sub_cycle ? parent->MaxRefRatio(i-1) : 1;
	}

	for(int i=0; i<=finest_level; i++) {
		ZZZAmr &advLevel = getLevel(i);
		dt_min[i] = advLevel.estimateTimeStep(dt_level[i]);
	}

//	if(post_regrid_flag == 1) {
//		// limit dt by pre-regrid values
//		for(int i=0; i<=finest_level; i++) {
//			dt_min[i] = std::min(dt_min[i], dt_level[i]);
//		}
//	}

	// find minimum among all levels
	Real dt0 = 1e100;
	int nFactor = 1;
	for(int i=0; i<=finest_level; i++) {
		nFactor *= n_cycle[i];
		dt0 = std::min(dt0, nFactor*dt_min[i]);
	}

	// limit DT by the value of stop_time
	const Real eps = 0.0001 * dt0;
	Real cur_time = state[State_Type].curTime();
	if(stop_time >= 0) {
		if(cur_time+dt0 > stop_time-eps) {
			dt0 = stop_time - cur_time;
		}
	}

	nFactor = 1;
	for(int i=0; i<=finest_level; i++) {
		nFactor *= n_cycle[i];
		dt_level[i] = dt0 / nFactor;
	}
}

int ZZZAmr::okToContinue() {
	if(level > 0) return 1;

	int ok = 1;
	const double dt_cutoff = 5e-20;
	if(parent->dtLevel(0) < dt_cutoff) {
		ok = 0;
	}

	return ok;
}

////////////////////////////////////////////////////////////////
// advance
////////////////////////////////////////////////////////////////


static void cpp_calc_fluxs(
		Real* phi, const int &phiLo1, const int &phiLo2, const int &phiHi1, const int &phiHi2,
		Real* fluxx, const int &fxLo1, const int &fxLo2, const int &fxHi1, const int &fxHi2,
		Real* fluxy, const int &fyLo1, const int &fyLo2, const int &fyHi1, const int &fyHi2,
		const int lo[], const int hi[], const Real dx[]
) {
//	FORTRAN_ARRAY_2D(state, phi, phiLo1, phiLo2, phiHi1, phiHi2);
//	FORTRAN_ARRAY_2D(fx, fluxx, fxLo1, fxLo2, fxHi1, fxHi2);
	FortArrayWrap<Real, 2> state(phi, phiLo1, phiHi1, phiLo2, phiHi2);

	FortArrayWrap<Real, 2> fx(fluxx, fxLo1, fxHi1, fxLo2, fxHi2);
	FortArrayWrap<Real, 2> fy(fluxy, fyLo1, fyHi1, fyLo2, fyHi2);

	// x-flux
	for(int j=lo[1]; j<=hi[1]; j++) {
		for(int i=lo[0]; i<=hi[0]+1; i++) {
			fx(i,j) = -(state(i,j) - state(i-1,j)) / dx[0];
		}
	}

	// y-flux
	for(int j=lo[1]; j<=hi[1]+1; j++) {
		for(int i=lo[0]; i<=hi[0]; i++) {
			fy(i,j) = -(state(i,j) - state(i,j-1)) / dx[1];
		}
	}
}

static void cpp_update_state(
		Real *oldData, int oldLo1, int oldLo2, int oldHi1, int oldHi2,
		Real *newData, int newLo1, int newLo2, int newHi1, int newHi2,
		Real *fluxx, int fxLo1, int fxLo2, int fxHi1, int fxHi2,
		Real *fluxy, int fyLo1, int fyLo2, int fyHi1, int fyHi2,
		const int lo[], const int hi[],
		const Real dx[], const Real &dt
) {
	// declare Fortran wrappers
	FortArrayWrap<Real,2> oldState(oldData, oldLo1, oldHi1, oldLo2, oldHi2);
	FortArrayWrap<Real,2> newState(newData, newLo1, newHi1, newLo2, newHi2);

	FortArrayWrap<Real,2> fx(fluxx, fxLo1, fxHi1, fxLo2, fxHi2);
	FortArrayWrap<Real,2> fy(fluxy, fyLo1, fyHi1, fyLo2, fyHi2);

	// update with flux divergence
	for(int j=lo[1]; j<=hi[1]; j++) {
		for(int i=lo[0]; i<=hi[0]; i++) {
			double dfx = (fx(i,j) - fx(i+1,j)) / dx[0];
			double dfy = (fy(i,j) - fy(i,j+1)) / dx[1];

			newState(i,j) = oldState(i,j) + dt * (dfx + dfy);
		}
	}
}

Real ZZZAmr::advance_explicitly(Real time, Real dt, int iteration, int ncycle) {
	int finest_level = parent->finestLevel();
	bool is_finest_level = level==finest_level ? true : false;

	if(doReflux && level<finest_level) {
		// reset flux registers to zero
		getFluxReg(level+1).setVal(0.0);
	}

	for(int k=State_Type; k<NUM_STATE_TYPE; k++) {
		state[k].allocOldData();
		state[k].swapTimeLevels(dt);
	}

	const Real prev_time = state[State_Type].prevTime();
	const Real curr_time = state[State_Type].curTime();

	MultiFab &oldState = get_old_data(State_Type);
	MultiFab &newState = get_new_data(State_Type);

	// test of NaN
	if(oldState.contains_nan(0, oldState.nComp(), 0)) {
		for(int i=0; i<oldState.nComp(); i++) {
			if(oldState.contains_nan(i, 1, 0)) {
				std::cout << "Testing component=" << i
						<< " for NaN" << std::endl;
				BoxLib::Abort("oldState has NaNs");
			}
		}
	}

	Real dt_new = dt;

	// get pointers to flux registers, or set them to zero if not there
	FluxRegister *fine = NULL;
	FluxRegister *current = NULL;
	if(doReflux && level<finest_level) {
		fine = &getFluxReg(level+1);
	}
	if(doReflux && level>0) {
		current = &getFluxReg(level);
	}

	AmrLevel &levelData = *this;
	Geometry g = levelData.Geom();

	MultiFab levelVolume;
	g.GetVolume(levelVolume, grids, numGhost);

	MultiFab levelArea[BL_SPACEDIM];
	for(int i=0; i<BL_SPACEDIM; i++) {
		g.GetFaceArea(levelArea[i], grids, i, numGhost);
	}

	// integrate, loop through the grids
	FArrayBox divu;
	FArrayBox grid_volume;
	FArrayBox dloga;
	FArrayBox area[BL_SPACEDIM];
	FArrayBox flux[BL_SPACEDIM];

	const Real *dx = geom.CellSize();
	const int *domain_lo = geom.Domain().loVect();
	const int *domain_hi = geom.Domain().hiVect();
	Real courant = -1e200;

	MultiFab fluxes[BL_SPACEDIM];
	if(doReflux && fine) {
		for(int i=0; i<BL_SPACEDIM; i++) {
			BoxArray edgeBA = newState.boxArray();
			edgeBA.surroundingNodes(i);

			fluxes[i].define(edgeBA, numState, 0, Fab_allocate);
			fluxes[i].setVal(0);
		}
	}

	for(FillPatchIterator fpi(*this, newState, numGhost, time, State_Type, 0, numState);
			fpi.isValid(); ++fpi)
	{
		int mfi_index = fpi.index();

		Box bx(fpi.UngrownBox());
		const int local_grow = 4;
		Box bx_grown(BoxLib::grow(bx, local_grow));

		// create FAB for extended grid values (including boundaries) and fill
		FArrayBox &state = fpi();
		FArrayBox &state_out = newState[fpi];

		grid_volume.resize(bx_grown, 1);
		grid_volume.copy(levelVolume[fpi]);

		for(int i=0; i<BL_SPACEDIM; i++) {
			area[i].resize(BoxLib::surroundingNodes(bx_grown, i));
			area[i].copy(levelArea[i][mfi_index]);
		}

#if (BL_SPACEDIM <= 2)
		dloga.resize(bx_grown);
		dloga.copy(dLogArea[0][mfi_index]);
#endif

		// allocate FABs for fluxes
		for(int i=0; i<BL_SPACEDIM; i++) {
			flux[i].resize(BoxLib::surroundingNodes(bx, i), numState);
		}

		Real cflLoc = -1e200;

		// TODO driver here
		{
			const int *lo = bx.loVect();
			const int *hi = bx.hiVect();

			// calculate flux
			cpp_calc_fluxs(BL_TO_FORTRAN(state),
					BL_TO_FORTRAN(flux[0]),
					BL_TO_FORTRAN(flux[1]),
					lo, hi, dx);

			// update
			cpp_update_state(
					BL_TO_FORTRAN(state), BL_TO_FORTRAN(state_out),
					BL_TO_FORTRAN(flux[0]), BL_TO_FORTRAN(flux[1]),
					lo, hi, dx, dt);
		}

		if(doReflux) {
			// NOTE: fluxes must be scaled by time for correct reflux!!
			for(int i=0; i<BL_SPACEDIM; i++) {
				flux[i].mult(dt);
			}

			if(fine) {
				for(int i=0; i<BL_SPACEDIM; i++) {
					fluxes[i][mfi_index].copy(flux[i]);
				}
			}

			if(current && iteration==ncycle) {
				for(int i=0; i<BL_SPACEDIM; i++) {
//					current->FineAdd(flux[i], i, mfi_index, 0, 0, numState, 1);
					current->FineAdd(flux[i], area[i], i, mfi_index, 0, 0, numState, 1);
				}
			}
		}

		courant = std::max(courant, cflLoc);
	} // end looping FPI

	if(doReflux && fine) {
		for(int i=0; i<BL_SPACEDIM; i++) {
//			fine->CrseInit(fluxes[i], i, 0, 0, numState, -1);
			fine->CrseInit(fluxes[i], levelArea[i], i, 0, 0, numState, -1);
		}
	}

	ParallelDescriptor::ReduceRealMax(courant);
	// check Courant number here
//	dt_new = dt / courant;
	dt_new = dt;

	// test of NaN
	if(newState.contains_nan(0, newState.nComp(), 0)) {
		for(int i=0; i<newState.nComp(); i++) {
			if(newState.contains_nan(i, 1, 0)) {
				std::cout << "Testing component=" << i
						<< " for NaN" << std::endl;
				BoxLib::Abort("newState has NaNs");
			}
		}
	}

	return dt_new;
}

Real ZZZAmr::advance_implicitly(Real time, Real dt, int iteration, int ncycle) {
	// TODO
	BoxLib::Error("Implicit time stepping is not supported!");

	int finest_level = parent->finestLevel();
	int is_finest_level = (level == finest_level);

	if(doReflux && level<finest_level) {
		// reset flux registers to zero
		getFluxReg(level+1).setVal(0.0);
	}

	for(int k=State_Type; k<NUM_STATE_TYPE; k++) {
		state[k].allocOldData();
		state[k].swapTimeLevels(dt);
	}

	const Real prev_time = state[State_Type].prevTime();
	const Real curr_time = state[State_Type].curTime();

	MultiFab &oldState = get_old_data(State_Type);
	MultiFab &newState = get_new_data(State_Type);

	// test of NaN
	if(oldState.contains_nan(0, oldState.nComp(), 0)) {
		for(int i=0; i<oldState.nComp(); i++) {
			if(oldState.contains_nan(i, 1, 0)) {
				std::cout << "Testing component=" << i
						<< " for NaN" << std::endl;
				BoxLib::Abort("oldState has NaNs");
			}
		}
	}

	Real dt_new = dt;

	{
		/*
		 * solver begin
		 */
		MultiFab phi(grids, 1, numGhost, Fab_allocate);
		MultiFab::Copy(phi, oldState, StatePhi, 0, 1, numGhost);

		if(verbose && ParallelDescriptor::IOProcessor()) {
			std::cout << "Solve for PHI at level=" << level << std::endl;
		}

		if(level==0 && !Geometry::isAllPeriodic()) {
			std::cout << "Must be periodic on all directions" << std::endl;
			BoxLib::Error(__FUNCTION__);
		}

		MacBndry bndry(grids, 1, geom);

		IntVect coarsen_ratio = level>0 ? parent->refRatio(level-1)
				: IntVect::TheZeroVector();

		// set Dirichlet BC in ghost cells, use to initialize boundary
		const int src_comp = 0;
		const int dst_comp = 0;
		const int num_comp = 1;
		const int num_grow = 1;

		if(level == 0) {
			bndry.setBndryValues(phi, src_comp, dst_comp, num_comp, physBC);
		} else {
			MultiFab coarsePhi;
			{ // get coarse PHI
				coarsePhi.clear();
				ZZZAmr &coarseLevel = getLevel(level-1);
				const MultiFab &coarseState = coarseLevel.get_new_data(State_Type);

				coarsePhi.define(coarseLevel.boxArray(), 1, 1, Fab_allocate);

				// copy coarser data, but DO NOT trust the ghost cells
				for(MFIter mfi(coarsePhi); mfi.isValid(); ++mfi) {
					coarsePhi[mfi].copy(coarseState[mfi], src_comp, dst_comp, num_comp);
				}
				coarsePhi.FillBoundary();

				const Geometry &coarseGeom = coarseLevel.Geom();
				coarseGeom.FillPeriodicBoundary(coarsePhi, true);
			}

			BoxArray coarsen_boxes = BoxArray(grids).coarsen(coarsen_ratio);

			const int in_rad = 0;
			const int out_rad = 1;
			const int extent_rad = 2;
			BndryRegister coarse_br(coarsen_boxes, in_rad, out_rad, extent_rad, num_comp);
			coarse_br.copyFrom(coarsePhi, coarsePhi.nGrow(), src_comp, dst_comp, num_comp);

			bndry.setBndryValues(coarse_br, src_comp, phi, src_comp, dst_comp, num_comp,
					coarsen_ratio, physBC);
		}

		// setup linear equation
		// alpha * A . phi - beta * div (B . grad phi) = RHS

		// RHS = phi^n
		MultiFab rhs(grids, 1, 0, Fab_allocate);
		MultiFab::Copy(rhs, oldState, StatePhi, 0, 1, 0);

		// we have alpha=1, A=I, beta=dt, B=I
		ABecLaplacian lap(bndry, geom.CellSize());
		const Real alpha = 1.0;
		const Real beta = dt;
		MultiFab acoefs;
		acoefs.define(grids, num_comp, num_grow, Fab_allocate);
		acoefs.setVal(1.0);
		MultiFab bcoefs[BL_SPACEDIM];
		for(int i=0; i<BL_SPACEDIM; i++) {
			BoxArray edgeNodes(grids);
			edgeNodes.surroundingNodes(i);
			bcoefs[i].define(edgeNodes, num_comp, num_grow, Fab_allocate);
			bcoefs[i].setVal(1.0);
		}

		// set coefficients
		lap.setScalars(alpha, beta);
		lap.setCoefficients(acoefs, bcoefs);

		if(verbose && ParallelDescriptor::IOProcessor()) {
			double norm = lap.norm();
			std::cout << "Solver at level=" << level
					<< "; norm=" << norm << std::endl;
		}

		// we solve using CG solver
		const int use_mg_pre = 1;
		CGSolver cg(lap, use_mg_pre);
		const int max_iter = 200;
		const Real tol = 1.0e-9;
		const Real tol_abs = -1;
		cg.setMaxIter(max_iter);

		int stat = cg.solve(phi, rhs, tol, tol_abs);

		if(stat != 0) { // error
			std::cout << "CG solver failed with error=" << stat << std::endl;
			std::cout << "at level=" << level << std::endl;
			BoxLib::Error(__FUNCTION__);
		}

		// assign PHI to state storage
		MultiFab::Copy(newState, phi, src_comp, StatePhi, num_comp, num_grow);
		if(verbose && ParallelDescriptor::IOProcessor()) {
			std::cout << "CG solver returned for level=" << level << std::endl;
		}

		/*
		 * solver end
		 */
	}


	// get pointers to flux registers, or set them to zero if not there
	FluxRegister *fine = NULL;
	FluxRegister *current = NULL;
	if(doReflux && level<finest_level) {
		fine = &getFluxReg(level+1);
	}
	if(doReflux && level>0) {
		current = &getFluxReg();
	}

	AmrLevel &levelData = *this;
	Geometry g = levelData.Geom();

	MultiFab levelVolume;
	g.GetVolume(levelVolume, grids, numGhost);

	MultiFab levelArea[BL_SPACEDIM];
	for(int i=0; i<BL_SPACEDIM; i++) {
		g.GetFaceArea(levelArea[i], grids, i, numGhost);
	}

	// integrate, loop through the grids
	FArrayBox divu, grid_volume;
	FArrayBox dloga, area[BL_SPACEDIM];
	FArrayBox flux[BL_SPACEDIM];

	const Real *dx = geom.CellSize();
	const int *domain_lo = geom.Domain().loVect();
	const int *domain_hi = geom.Domain().hiVect();
	Real courant = -1e200;

	MultiFab fluxes[BL_SPACEDIM];
	if(doReflux && fine!=NULL) {
		for(int i=0; i<BL_SPACEDIM; i++) {
			BoxArray edgeBA = newState.boxArray();
			edgeBA.surroundingNodes(i);

			fluxes[i].define(edgeBA, numState, 0, Fab_allocate);
			fluxes[i].setVal(0);
		}
	}

	for(FillPatchIterator fpi(*this, newState, numGhost, time, State_Type, 0, numState);
			fpi.isValid(); ++fpi)
	{
		int mfi_index = fpi.index();

		Box bx(fpi.UngrownBox());
		const int local_grow = numGhost;
		Box bx_grown(BoxLib::grow(bx, local_grow));

		// create FAB for extended grid values (including boundaries) and fill
		//		FArrayBox &state = fpi();
		//		FArrayBox &state_out = newState[fpi];
		FArrayBox &state = newState[fpi];

		grid_volume.resize(bx_grown, 1);
		grid_volume.copy(levelVolume[fpi]);

		for(int i=0; i<BL_SPACEDIM; i++) {
			area[i].resize(BoxLib::surroundingNodes(bx_grown, i));
			area[i].copy(levelArea[i][mfi_index]);
		}

#if (BL_SPACEDIM <= 2)
		dloga.resize(bx_grown);
		dloga.copy(dLogArea[0][mfi_index]);
#endif

		// allocate FABs for fluxes
		for(int i=0; i<BL_SPACEDIM; i++) {
			flux[i].resize(BoxLib::surroundingNodes(bx_grown, i), numState);
		}

		Real cflLoc = -1e200;

		// TODO driver here
		{
			const int *lo = bx.loVect();
			const int *hi = bx.hiVect();

			// calculate flux
			cpp_calc_fluxs(BL_TO_FORTRAN(state),
					BL_TO_FORTRAN(flux[0]),
					BL_TO_FORTRAN(flux[1]),
					lo, hi, dx);

			//			// update
			//			cpp_update_state(
			//					BL_TO_FORTRAN(state), BL_TO_FORTRAN(state_out),
			//					BL_TO_FORTRAN(flux[0]), BL_TO_FORTRAN(flux[1]),
			//					lo, hi, dx, dt);
		}

		if(doReflux) {
			if(fine != NULL) {
				for(int i=0; i<BL_SPACEDIM; i++) {
					fluxes[i][mfi_index].copy(flux[i]);
				}
			}
			if(current != NULL) {
				for(int i=0; i<BL_SPACEDIM; i++) {
					current->FineAdd(flux[i], i, mfi_index, 0, 0, numState, 1);
				}
			}
		}

		courant = std::max(courant, cflLoc);
	} // end looping FPI

	if(doReflux && fine!=NULL) {
		for(int i=0; i<BL_SPACEDIM; i++) {
			fine->CrseInit(fluxes[i], i, 0, 0, numState, -1);
		}
	}

	ParallelDescriptor::ReduceRealMax(courant);
	// check Courant number here
	//	dt_new = dt / courant;
	dt_new = dt;

	// test of NaN
	if(newState.contains_nan(0, newState.nComp(), 0)) {
		for(int i=0; i<newState.nComp(); i++) {
			if(newState.contains_nan(i, 1, 0)) {
				std::cout << "Testing component=" << i
						<< " for NaN" << std::endl;
				BoxLib::Abort("newState has NaNs");
			}
		}
	}

	return dt_new;
}

/**
 * Do an integration step on this level, and returns maximum safe time step.
 */
Real ZZZAmr::advance(Real time, Real dt, int iteration, int ncycle) {
	Real dt_new = dt;

	const bool isExplicit = !doImplicitStep;
	if(isExplicit) {
		dt_new = advance_explicitly(time, dt, iteration, ncycle);
	} else {
		dt_new = advance_implicitly(time, dt, iteration, ncycle);
	}

	return dt_new;
}


//// the old version, work for level=0
//Real ZZZAmr::advance(Real time, Real dt, int iteration, int ncycle) {
//	MultiFab flux[BL_SPACEDIM];
//
//	for(int i=0; i<BL_SPACEDIM; i++) {
//		BoxArray edgeGrids(grids);
//		edgeGrids.surroundingNodes(i);
//
//		int nvar = numState, nghost = 0;
//		flux[i].define(edgeGrids, nvar, nghost, Fab_allocate);
//		flux[i].setVal(0);
//	}
//
//	for(int k=State_Type; k<NUM_STATE_TYPE; k++) {
//		state[k].allocOldData();
//		state[k].swapTimeLevels(dt);
//	}
//
//	MultiFab &oldData = get_old_data(State_Type);
//	MultiFab &newData = get_new_data(State_Type);
//
//	// fill ghost cells covered by valid cells
//	oldData.FillBoundary();
//	// fill periodic boundary ghost cells
//	geom.FillPeriodicBoundary(oldData);
//
//	const Real *dx = geom.CellSize();
//
////	for(MFIter mfi(oldData); mfi.isValid(); ++mfi) {
////		const Box &validBox = mfi.validbox();
////
////		BL_FORT_PROC_CALL(CALC_FLUX, calc_flux)(
////				BL_TO_FORTRAN(oldData[mfi]),
////				BL_TO_FORTRAN(flux[0][mfi]),
////				BL_TO_FORTRAN(flux[1][mfi]),
////				validBox.loVect(), validBox.hiVect(), dx
////		);
////	}
//
//	FillPatchIterator fpi(*this, oldData, numGhost, time, State_Type, 0, numState);
//	for(; fpi.isValid(); ++fpi) {
//		const Box &validBox = fpi.validbox();
//
//		FArrayBox &oldState = fpi();
//
////		if(ParallelDescriptor::IOProcessor()) {
////			std::cout << "Fluxes for patch " << fpi.index() << std::endl;
////		}
//
////		BL_FORT_PROC_CALL(CALC_FLUX, calc_flux)(
////				BL_TO_FORTRAN(oldState),
////				BL_TO_FORTRAN(flux[0][fpi]),
////				BL_TO_FORTRAN(flux[1][fpi]),
////				validBox.loVect(), validBox.hiVect(), dx
////		);
//		cpp_calc_fluxs(
//				BL_TO_FORTRAN(oldState),
//				BL_TO_FORTRAN(flux[0][fpi]),
//				BL_TO_FORTRAN(flux[1][fpi]),
//				validBox.loVect(), validBox.hiVect(), dx
//		);
//	}
//
//	for(MFIter mfi(oldData); mfi.isValid(); ++mfi) {
//		const Box &validBox = mfi.validbox();
//
////		BL_FORT_PROC_CALL(UPDATE_STATE, update_state)(
////				BL_TO_FORTRAN(oldData[mfi]), BL_TO_FORTRAN(newData[mfi]),
////				BL_TO_FORTRAN(flux[0][mfi]), BL_TO_FORTRAN(flux[1][mfi]),
////				validBox.loVect(), validBox.hiVect(), dx, &dt
////		);
//		cpp_update_state(
//				BL_TO_FORTRAN(oldData[mfi]), BL_TO_FORTRAN(newData[mfi]),
//				BL_TO_FORTRAN(flux[0][mfi]), BL_TO_FORTRAN(flux[1][mfi]),
//				validBox.loVect(), validBox.hiVect(), dx, dt
//		);
//	}
//
//	return dt;
//}



////////////////////////////////////////////////////////////////
// post hooks
////////////////////////////////////////////////////////////////

void ZZZAmr::post_timestep(int iteration) {
	// integration cycle on fine level grids is done
	// do post_timestep stuff here

	const int finestLevel = parent->finestLevel();

	if(doReflux && level<finestLevel) {
		// this level is the coarse one
//		MultiFab &newCoarseState = get_new_data(State_Type);

		reflux();

		// refluxing changes the values of coarse cells underneath fine grids
//		if(level < finestLevel) {
//			averageDown();
//		}
	}

	if(level < finestLevel) {
		averageDown();
	}

	if(level == 0) {
		int nstep = parent->levelSteps(0);
		// sum integrated quantities to check conservation?
	}
}

void ZZZAmr::post_restart() {
	BoxLib::Error(__FUNCTION__);
	// do not support restart
}


void ZZZAmr::post_regrid(int lbase, int new_finest) {
	// do nothing
}

/*
 * Operations after initialization.
 */
void ZZZAmr::post_init(Real stop_time) {
	if(level > 0) return;

	/*
	 * average data down from finer levels,
	 * so that conserved data is consistent between levels
	 */
	int finest_level = parent->finestLevel();
	// NOTE: must be done in a reversed order!
	for(int k=finest_level-1; k>=0; k--) {
		getLevel(k).averageDown();
	}
}



static void avgdown_xmajor(
		FortArrayWrap<Real,3> &crse, FortArrayWrap<Real,2> &cv,
		FortArrayWrap<Real,3> &fine, FortArrayWrap<Real,2> &fv,
		const int &nvar,
		const int lo[], const int hi[],
		const int &lratx, const int &lraty)
{
	for(int n=0; n<nvar; n++) {
		// set coarse grid to zero on overlap
		for(int jc=lo[1]; jc<=hi[1]; jc++) {
			for(int ic=lo[0]; ic<=hi[0]; ic++) {
				crse(ic,jc,n) = 0.0;
			}
		}

		// sum fine data
		for(int joff=0; joff<lraty; joff++) {
			for(int jc=lo[1]; jc<=hi[1]; jc++) {
				int j = jc * lraty + joff;

				for(int ioff=0; ioff<lratx; ioff++) {
					for(int ic=lo[0]; ic<=hi[0]; ic++) {
						int i = ic * lratx + ioff;

						crse(ic,jc,n) += fv(i,j) * fine(i,j,n);
					}
				}
			}
		}

		// divide by volume weight
		for(int jc=lo[1]; jc<=hi[1]; jc++) {
			for(int ic=lo[0]; ic<=hi[0]; ic++) {
				crse(ic,jc,n) /= cv(ic,jc);
			}
		}
	}
}
static void avgdown_ymajor(
		FortArrayWrap<Real,3> &crse, FortArrayWrap<Real,2> &cv,
		FortArrayWrap<Real,3> &fine, FortArrayWrap<Real,2> &fv,
		const int &nvar,
		const int lo[], const int hi[],
		const int &lratx, const int &lraty)
{
	for(int n=0; n<nvar; n++) {
		// set coarse grid to zero on overlap
		for(int ic=lo[0]; ic<=hi[0]; ic++) {
			for(int jc=lo[1]; jc<=hi[1]; jc++) {
				crse(ic,jc,n) = 0.0;
			}
		}

		// sum fine data
		for(int ioff=0; ioff<lratx; ioff++) {
			for(int ic=lo[0]; ic<=hi[0]; ic++) {
				int i = ic * lratx + ioff;

				for(int joff=0; joff<lraty; joff++) {
					for(int jc=lo[1]; jc<=hi[1]; jc++) {
						int j = jc * lraty + joff;

						crse(ic,jc,n) += fv(i,j) * fine(i,j,n);
					}
				}
			}
		}

		// divide by volume weight
		for(int ic=lo[0]; ic<=hi[0]; ic++) {
			for(int jc=lo[1]; jc<=hi[1]; jc++) {
				crse(ic,jc,n) /= cv(ic,jc);
			}
		}
	}
}
/**
 * Volume-weight average the fine grid data onto the coarse
 * grid.  Overlap is given in coarse grid coordinates.
 * INPUTS / OUTPUTS:
! ::  crse      <=  coarse grid data
! ::  clo,chi    => index limits of crse array interior
! ::  ngc        => number of ghost cells in coarse array
! ::  nvar	 => number of components in arrays
! ::  fine       => fine grid data
! ::  flo,fhi    => index limits of fine array interior
! ::  ngf        => number of ghost cells in fine array
! ::  rfine      => (ignore) used in 2-D RZ calc
! ::  lo,hi      => index limits of overlap (crse grid)
! ::  lrat       => refinement ratio
! ::
 *
 * NOTE: Assumes all data cell centered
 */
static void avgdown_volumeWeighted(
		Real *coarseData, // coarse grid data
		const int &clo1, const int &clo2,
		const int &chi1, const int &chi2,
		const int &nvar,
		Real *coarseVolData, // coarse volume
		const int &cvlo1, const int &cvlo2,
		const int &cvhi1, const int &cvhi2,
		Real *fineData, // fine grid data
		const int &flo1, const int &flo2,
		const int &fhi1, const int &fhi2,
		Real *fineVolData, // fine volume
		const int &fvlo1, const int &fvlo2,
		const int &fvhi1, const int &fvhi2,
		const int lo[], const int hi[], // limits for coarse grid
		const int lrat[])
{
	// data view
	FortArrayWrap<Real,3> crse(coarseData, clo1,chi1, clo2,chi2, 0,nvar-1);
	FortArrayWrap<Real,2> cv(coarseVolData, cvlo1,cvhi1, cvlo2,cvhi2);
	FortArrayWrap<Real,3> fine(fineData, flo1,fhi1, flo2,fhi2, 0, nvar-1);
	FortArrayWrap<Real,2> fv(fineVolData, fvlo1,fvhi1, fvlo2,fvhi2);

	const int lratx = lrat[0];
	const int lraty = lrat[1];
	const int lenx = hi[0] - lo[0] + 1;
	const int leny = hi[1] - lo[1] + 1;
	const int maxlen = std::max(lenx, leny);

	if(lenx == maxlen) {
		avgdown_xmajor(crse, cv, fine, fv, nvar,
				lo, hi, lratx, lraty);
	} else {
		avgdown_ymajor(crse, cv, fine, fv, nvar,
				lo, hi, lratx, lraty);
	}
}

void ZZZAmr::averageDown() {
	const int finest_level = parent->finestLevel();
	if(level == finest_level) return;

	averageDown(State_Type);
}
void ZZZAmr::averageDown(int state_index) {
	if(level == parent->finestLevel()) return;

	ZZZAmr &fineLevel = getLevel(level+1);
	MultiFab &coarseState = get_new_data(state_index);
	MultiFab &fineState = fineLevel.get_new_data(state_index);
	MultiFab &fineVolume = fineLevel.volume;

	const int ncomp = fineState.nComp();

	BL_ASSERT(coarseState.boxArray() == volume.boxArray());
	BL_ASSERT(fineState.boxArray() == fineVolume.boxArray());

	// coarsen the finer data
	BoxArray coarsenFineStateBoxes(fineState.boxArray().size());
	for(int i=0; i<fineState.boxArray().size(); i++) {
		coarsenFineStateBoxes.set(i,
				BoxLib::coarsen(fineState.box(i), fine_ratio));
	}

	MultiFab coarsenFineState(coarsenFineStateBoxes, ncomp, 0);
	MultiFab coarsenFineVolume(coarsenFineStateBoxes, 1, 0);
	coarsenFineVolume.copy(volume);

	for(MFIter mfi(fineState); mfi.isValid(); ++mfi) {
		const int i = mfi.index();

		const Box &ovlp = coarsenFineStateBoxes[i];

		FArrayBox &coarsenFab = coarsenFineState[i];
		FArrayBox &coarsenVol = coarsenFineVolume[i];

		FArrayBox &fineFab = fineState[i];
		FArrayBox &fineVol = fineVolume[i];

		// average-down implementation
		avgdown_volumeWeighted(
				BL_TO_FORTRAN(coarsenFab), ncomp,
				BL_TO_FORTRAN(coarsenVol),
				BL_TO_FORTRAN(fineFab),
				BL_TO_FORTRAN(fineVol),
				ovlp.loVect(), ovlp.hiVect(),
				fine_ratio.getVect()
		);
	}

	// assign coarsen values
	coarseState.copy(coarsenFineState);
}

void ZZZAmr::reflux() {
	if(level >= parent->finestLevel()) {
		BoxLib::Error(__FUNCTION__);
	}

	getFluxReg(level+1).Reflux(
			get_new_data(State_Type),
			volume, 1.0,
			0, 0, numState, geom);

	if(verbose && ParallelDescriptor::IOProcessor()) {
		std::cout << __FUNCTION__ << " at level=" << level << std::endl;
	}
}

//void ZZZAmr::syncInterp(
//		MultiFab &coarseSync, int coarse_level,
//		MultiFab &fineSync, int fine_level,
//		IntVect &ratio,
//		int src_comp, int dst_comp, int num_comp, int increment,
//		Real coarseDt, int **bc_orig_qty,
//		SyncInterpType which_interp,
//		int state_comp)
//{
//}




////////////////////////////////////////////////////////////////
// error estimation
////////////////////////////////////////////////////////////////

void ZZZAmr::errorEst(TagBoxArray &tags, int clearval, int tagval,
		Real time, int n_error_buf, int ngrow)
{
	const int *domain_lo = geom.Domain().loVect();
	const int *domain_hi = geom.Domain().hiVect();
	const Real *dx = geom.CellSize();
	const Real *prob_lo = geom.ProbLo();


	for(int j=0; j<errList.size(); j++) {
		MultiFab *mf_ptr = derive(errList[j].name(), time, errList[j].nGrow());
		if(mf_ptr == NULL) {
			std::string msg = "ZZZAmr::errorEst: failed to derive ";
			msg += errList[j].name();
			BoxLib::Error(msg.c_str());
		}
		std::auto_ptr<MultiFab> mf(mf_ptr);

		for(MFIter mfi(*mf); mfi.isValid(); ++mfi) {
			const int idx = mfi.index();
			RealBox gridloc(grids[idx], dx, prob_lo);

			TagBox &tb = tags[idx];

			Array<int> itags = tb.tags();
			int *tptr = itags.dataPtr();
			const int *tlo = tb.loVect();
			const int *thi = tb.hiVect();

			const int *lo = mfi.validbox().loVect();
			const int *hi = mfi.validbox().hiVect();
			const Real *xlo = gridloc.lo();

			FArrayBox &fab = (*mf)[mfi];
			Real *data = fab.dataPtr();
			const int *dlo = fab.loVect();
			const int *dhi = fab.hiVect();
			const int ncomp = fab.nComp();

			errList[j].errFunc()(
					tptr, ARLIM(tlo), ARLIM(thi),
					&tagval, &clearval,
					data, ARLIM(dlo), ARLIM(dhi),
					lo, hi, &ncomp, domain_lo, domain_hi,
					dx, xlo, prob_lo, &time, &level
			);

			// set the tags in the TagBox
			if(allow_untagging == 1) {
				tb.tags_and_untags(itags);
			} else {
				tb.tags(itags);
			}
		}
	}
}


////////////////////////////////////////////////////////////////
// initialization routines
////////////////////////////////////////////////////////////////

static void cpp_init_data(
		const int &level, const Real &time,
		const int lo[], const int hi[],
		const int &nstate,
		Real *stateData,
		const int &sl1, const int &sl2,
		const int &sh1, const int &sh2,
		const Real dx[],
		const Real xlow[], const Real xhi[])
{
	// data wrapper for Fortran array
	FortArrayWrap<Real,3> state(stateData, sl1,sh1, sl2,sh2, 0,nstate-1);

	const int IPhi = ZZZAmr::StatePhi;

    for(int j=lo[1]; j<=hi[1]; j++) {
    	for(int i=lo[0]; i<=hi[0]; i++) {
    		double x = xlow[0] + (i-lo[0]+0.5) * dx[0];
    		double y = xlow[1] + (j-lo[1]+0.5) * dx[1];
    		double r2 = pow(x-0.25,2) + pow(y-0.25,2);

    		state(i,j,IPhi) = 1.0 + exp(-r2*100);
    	}
    }
}

/**
 * Initialize grid data at problem start-up.
 */
void ZZZAmr::initData() {
	const int nState = numState;
	const Real *dx = geom.CellSize();
	MultiFab &newState = get_new_data(State_Type);
	const Real currentTime = state[State_Type].curTime();

	if(verbose && ParallelDescriptor::IOProcessor()) {
		std::cout << "Initializing data at level=" << level << std::endl;
	}

	// initialize data
	for(MFIter mfi(newState); mfi.isValid(); ++mfi) {
		const int fabIndex = mfi.index();

		RealBox gridLoc(grids[fabIndex], geom.CellSize(), geom.ProbLo());

		const Box &validBox = mfi.validbox();
		const int *lo = validBox.loVect();
		const int *hi = validBox.hiVect();

		FArrayBox &fab = newState[mfi];
		fab.setVal(0.0);

		cpp_init_data(level, currentTime, lo, hi, nState,
				BL_TO_FORTRAN(fab),
				dx, gridLoc.lo(), gridLoc.hi());
	}

	if(verbose && ParallelDescriptor::IOProcessor()) {
		std::cout << "Initialized data at level=" << level << std::endl;
	}
}


/**
 * Called for re-grid, to build based on the old level.
 */
void ZZZAmr::init(AmrLevel &old) {
	ZZZAmr *oldLevel = static_cast<ZZZAmr*>(&old);

	// create new grid data by fill-patching from the old
	Real dtNew = parent->dtLevel(level);
	Real curTime = oldLevel->state[State_Type].curTime();
	Real prevTime = oldLevel->state[State_Type].prevTime();
	Real dtOld = curTime - prevTime;
	setTimeLevel(curTime, dtOld, dtNew);

	MultiFab &newState = get_new_data(State_Type);
	for(FillPatchIterator fpi(old, newState, 0, curTime, State_Type, 0, numState);
			fpi.isValid(); ++fpi) {
		const FArrayBox &fabOld = fpi();
		newState[fpi].copy(fabOld);
	}
}
/**
 * Called for re-grid, to initialize the data on a new level
 * that did not exist before regridding.
 */
void ZZZAmr::init() {
	ZZZAmr &baseLevel = getLevel(level-1);

	Real dt = parent->dtLevel(level);
	Real curTime = baseLevel.state[State_Type].curTime();
	Real prevTime = baseLevel.state[State_Type].prevTime();

	const int maxRefineRatio = parent->MaxRefRatio(level-1);
	Real dtOld = (curTime - prevTime) / maxRefineRatio;

	setTimeLevel(curTime, dtOld, dt);

	MultiFab &newState = get_new_data(State_Type);
	FillCoarsePatch(newState, 0, curTime, State_Type, 0, numState);
}





/**
 * factory class for ZZZAmr
 */
class ZZZAmrBluid : public LevelBld {
    //
    // Perform any problem-dependent setup such as physical
    // boundary condition and derived quantities.
    // This is a pure virtual function and hence MUST
    // be implemented by derived classes.
    //
    virtual void variableSetUp () {
    	ZZZAmr::variableSetUp();
    }
    //
    // Perform any problem-dependent cleanup.
    // This is a pure virtual function and hence MUST
    // be implemented by derived classes.
    //
    virtual void variableCleanUp () {
    	ZZZAmr::variableCleanUp();
    }
    //
    // This is a virtual constructor for types derived
    // from AmrLevel.  The derived type is initialized
    // with the default constructor.
    // This is a pure virtual function and hence MUST
    // be implemented by derived classes.
    //
    virtual AmrLevel* operator() () {
    	return new ZZZAmr;
    }
    //
    // This is a virtual constructor for types derived
    // from AmrLevel.  The derived type is initialized
    // with the five specified variables.
    // This is a pure virtual function and hence MUST
    // be implemented by derived classes.
    //
    virtual AmrLevel* operator() (Amr&            papa,
                                  int             lev,
                                  const Geometry& geom_lev,
                                  const BoxArray& ba,
                                  Real            time)
    {
    	return new ZZZAmr(papa, lev, geom_lev, ba, time);
    }

public:
    static ZZZAmrBluid instance;
};

ZZZAmrBluid ZZZAmrBluid::instance;


/**
 *
 */
LevelBld* getLevelBld() {
	return &(ZZZAmrBluid::instance);
}


