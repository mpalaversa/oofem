/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2013   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */
 
#ifndef pltqd4dkt_h
#define pltqd4dkt_h

#include "qdplate.h"
#include "femcmpnn.h"
#include "error.h"
#include "floatmatrix.h"
#include "floatarray.h"

#include <cstdio>
#include <vector>
#include <memory>

#define _IFT_PltQd4DKT_Name "pltqd4dkt"

namespace oofem {
	class PltQd4DKT : public QdPlate
	{
	public:
		PltQd4DKT(int n, Domain* d);
		virtual ~PltQd4DKT() {}

		void computeBmatrixAt(double xi, double eta, FloatMatrix& answer) override { }
		void computeStrainVectorAt(FloatArray& answer, double xi, double eta, TimeStep* tStep) override { }

		const char* giveClassName() const override { return "PltQd4DKT"; }
		const char* giveInputRecordName() const override { return _IFT_PltQd4DKT_Name; }
		void computeStressVector(FloatArray& answer, const FloatArray& strain, GaussPoint* gp, TimeStep* tStep) override { }
		void computeConstitutiveMatrixAt(FloatMatrix& answer, MatResponseMode rMode, GaussPoint* gp, TimeStep* tStep) override {};

		// giveInternalForcesVector is used only in non-linear analysis. This should be changed when non-linear analysis capabilities are implemented.
		void giveInternalForcesVector(FloatArray& answer, TimeStep* tStep, int useUpdatedGpRecord) override { answer.resize(12); }
	};
}
#endif