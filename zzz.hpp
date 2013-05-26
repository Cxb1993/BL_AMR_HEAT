/*
 * zzz.hpp
 *
 *  Created on: 2013/04/08
 *      Author: homu
 */

#ifndef ZZZ_HPP_
#define ZZZ_HPP_

#include <iostream>


#include "BC_TYPES.H"
#include "AmrLevel.H"
#include "Amr.H"
#include "ErrorList.H"
#include "FluxRegister.H"


enum StateType {
	State_Type = 0,

	NUM_STATE_TYPE
};


class ZZZAmr : public AmrLevel
{
public:
	// TODO
	ZZZAmr();
	ZZZAmr(Amr &papa, int level, const Geometry &levelGeom,
			const BoxArray &boxList, Real time);

	virtual ~ZZZAmr();


    //
    // A string written as the first item in writePlotFile() at
    // level zero. This MUST be defined by each derived class.
    // It is so we can distinguish between different types of
    // plot files.  This is a pure virtual function and hence MUST
    // be implemented by derived classes.
    //
    virtual std::string thePlotFileType () const;
    //
    // Write plot file stuff to specified directory.  This is a
    // pure virtual function and hence MUST be implemented by
    // derived classes.
    //
    virtual void writePlotFile (const std::string& dir,
                                std::ostream&      os,
                                VisMF::How         how = VisMF::OneFilePerCPU);


    //
    // Compute the initial time step.  This is a pure virtual function
    // and hence MUST be implemented by derived classes.
    //
    virtual void computeInitialDt (int                   finest_level,
                                   int                   sub_cycle,
                                   Array<int>&           n_cycle,
                                   const Array<IntVect>& ref_ratio,
                                   Array<Real>&          dt_level,
                                   Real                  stop_time);
    //
    // Compute the next time step.  This is a pure virtual function
    // and hence MUST be implemented by derived classes.
    //
    virtual void computeNewDt (int                   finest_level,
                               int                   sub_cycle,
                               Array<int>&           n_cycle,
                               const Array<IntVect>& ref_ratio,
                               Array<Real>&          dt_min,
                               Array<Real>&          dt_level,
                               Real                  stop_time,
                               int                   post_regrid_flag);


    //
    // Do an integration step on this level.  Returns maximum safe
    // time step.  This is a pure virtual function and hence MUST
    // be implemented by derived classes.
    //
    virtual Real advance (Real time,
                          Real dt,
                          int  iteration,
                          int  ncycle);

    //
    // Contains operations to be done after a timestep.  This is a
    // pure virtual function and hence MUST be implemented by derived
    // classes.
    //
    virtual  void post_timestep (int iteration);
    //
    // Operations to be done after restart.  This is a pure virtual
    // function and hence MUST be implemented by derived classes.
    //
    virtual  void post_restart ();
    //
    // Operations to be done after regridding (like avgDown).
    // This is a pure virtual function and hence MUST be
    // implemented by derived classes.
    //
    virtual  void post_regrid (int lbase,
                               int new_finest);
    //
    // Operations to be done after initialization.
    // This is a pure virtual function and hence MUST be
    // implemented by derived classes.
    //
    virtual  void post_init (Real stop_time);

    //
    // Is it ok to continue the calculation?
    // This is a pure virtual function and hence MUST be
    // implemented by derived classes.
    //
    virtual  int okToContinue ();

    //
    // Init grid data at problem start-up.
    // This is a pure virtual function and hence MUST be
    // implemented by derived classes.
    //
    virtual void initData ();
    //
    // Init data on this level from another AmrLevel (during regrid).
    // This is a pure virtual function and hence MUST be
    // implemented by derived classes.
    //
    virtual void init (AmrLevel &old);
    //
    // Init data on this level after regridding if old AmrLevel
    // did not previously exist. This is a pure virtual function
    // and hence MUST be implemented by derived classes.
    //
    virtual void init ();

    //
    // Error estimation for regridding. This is a pure virtual
    // function and hence MUST be implemented by derived classes.
    //
    virtual void errorEst (TagBoxArray& tb,
                           int          clearval,
                           int          tagval,
                           Real         time,
                           int          n_error_buf = 0,
                           int          ngrow = 0);


	// define data descriptors
	static void variableSetUp();
	static void readParams();
	static void checkBoundaryConsistency();

	// clean up data descriptors at the end of run
	static void variableCleanUp();

	// define tagging functions
	static void ErrorSetUp();

protected:

	inline ZZZAmr& getLevel(int iLevel) {
		AmrLevel *pLevel = &parent->getLevel(iLevel);
		return *static_cast<ZZZAmr*>(pLevel);
	}

    inline FluxRegister& getFluxReg() {
    	BL_ASSERT(fluxReg);
    	return *fluxReg;
    }
    inline FluxRegister& getFluxReg(int lev) {
    	return getLevel(lev).getFluxReg();
    }

    inline MultiFab& getVolume() {
    	return volume;
    }
    inline MultiFab* getArea() {
    	return area;
    }
    inline MultiFab& getArea(int dir) {
    	return area[dir];
    }

    void buildMetrics();

    Real initialTimeStep();
	Real estimateTimeStep(const Real dt_old);

	/**
	 * take the average from finer level,
	 * this level is the coarse one
	 */
	void averageDown();
	void averageDown(int state_index);

	void reflux();

//	// interpolate cell-centered Sync. correction
//	// coarse -> fine
//	enum SyncInterpType {
//		PC_T, CellCons_T, CellConsLin_T, CellConsProt_T,
//	};
//
//	void syncInterp(
//			MultiFab &coarseSync, int coarse_level,
//			MultiFab &fineSync, int fine_level,
//			IntVect &ratio,
//			int src_comp, int dst_comp, int num_comp, int increment,
//			Real coarseDt, int **bc_orig_qty,
//			SyncInterpType which_interp = CellCons_T,
//			int state_comp = -1);

	Real advance_explicitly(Real time, Real dt, int iteration, int ncycle);
	Real advance_implicitly(Real time, Real dt, int iteration, int ncycle);




	FluxRegister *fluxReg;

	MultiFab volume;
	MultiFab area[BL_SPACEDIM];
	MultiFab dLogArea[1];

	Array<Array<Real> > radius; // e.g. for cylindrical coordinates

public:
	/*
	 * static members
	 */
	static int verbose;

	static Real cfl;

	static int doReflux;

	static int doImplicitStep;

	//
	static BCRec physBC;

	static int numGhost;
	static int numState;

	static int radiusGrow;

	static int StatePhi;

	static int allow_untagging;
	static ErrorList errList;
};





#endif /* ZZZ_HPP_ */
