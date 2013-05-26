/*
 * fortranFunctions.h
 *
 *  Created on: 2013/04/09
 *      Author: homu
 */

#ifndef FORTRANFUNCTIONS_H_
#define FORTRANFUNCTIONS_H_


#include <BLFort.H>
#include <REAL.H>


/**
 * default Fortran filler function
 */
BL_FORT_PROC_DECL(FILCC, filcc) (
		BL_FORT_FAB_ARG(state),
		const int dlo[], const int dhi[],
		const Real dx[], const Real glo[],
		const int bc[]
);


/**
 * fills the boundary data for quantity PHI
 */
BL_FORT_PROC_DECL(PHI_FILL, phi_fill) (
		BL_FORT_FAB_ARG(state),
		const int dlo[], const int dhi[],
		const Real dx[], const Real glo[],
		const Real *time, const int bc[]
);


/**
 * the routine to initialize the problem
 */
BL_FORT_PROC_DECL(INIT_DATA, init_data) (
		const int &level, const Real &time,
		const int lo[], const int hi[],
		const int &numState,
		BL_FORT_FAB_ARG(state),
		const Real dx[],
		const Real xlow[], const Real xhi[]
);


/**
 *
 */
BL_FORT_PROC_DECL(CALC_FLUX, calc_flux) (
		BL_FORT_FAB_ARG(old),
		BL_FORT_FAB_ARG(fluxx), BL_FORT_FAB_ARG(fluxy),
		const int lo[], const int hi[],
		const Real dx[]
);

BL_FORT_PROC_DECL(UPDATE_STATE, update_state) (
		BL_FORT_FAB_ARG(oldData), BL_FORT_FAB_ARG(newData),
		BL_FORT_FAB_ARG(fluxx), BL_FORT_FAB_ARG(fluxy),
		const int lo[], const int hi[],
		const Real dx[], const Real *dt
);

#endif /* FORTRANFUNCTIONS_H_ */



