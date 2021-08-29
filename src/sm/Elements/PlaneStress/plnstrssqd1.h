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
 
#ifndef plnstrssqd1_h
#define plnstrssqd1_h

#include "sm/Elements/PlaneStress/qdmembrane.h"
#include "femcmpnn.h"
#include "error.h"
#include "floatmatrix.h"
#include "floatarray.h"

#include <cstdio>
#include <vector>
#include <memory>

#define _IFT_PlnStrssQd1_Name "plnstrssqd1"
#define _IFT_PlnStrssQd1_csClass "csclass"
#define _IFT_PlnStrssQd1_outputAtXY "outputatxy"
#define _IFT_PlnStrssQd1_outputType "outputtype"

namespace oofem {
	class FEI2dQuadLin;

	class PlnStrssQd1 : public QdMembrane
	{
		friend class ShellQd41;
	protected:
		static FEI2dQuadLin interpolation;

	public:
		PlnStrssQd1(int n, Domain* d);
		virtual ~PlnStrssQd1() {}

		bool computeGtoLRotationMatrix(FloatMatrix& answer) override;
		void giveDofManDofIDMask(int inode, IntArray&) const override;
		FEInterpolation* giveInterpolation() const override;
		MaterialMode giveMaterialMode() override { return _PlaneStress; }

		void computeBmatrixAt(double xi, double eta, FloatMatrix& answer) override;
		const char* giveClassName() const override { return "PlnStrssQd1"; }
		const char* giveInputRecordName() const override { return _IFT_PlnStrssQd1_Name; }
		void initializeFrom(InputRecord& ir) override;

		// giveInternalForcesVector is used only in non-linear analysis. This should be changed when non-linear analysis capabilities are implemented.
		//void giveInternalForcesVector(FloatArray& answer, TimeStep* tStep, int useUpdatedGpRecord) override { answer.resize(8); }
	};
}
#endif