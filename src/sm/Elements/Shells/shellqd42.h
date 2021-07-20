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
 
#ifndef shellqd42_h
#define shellqd42_h

#include "sm/Elements/PlaneStress/plnstrssqd1rot.h"
#include "sm/Elements/Plates/pltqd4dkt.h"

#include "qdshell.h"
#include "femcmpnn.h"
#include "error.h"
#include "floatmatrix.h"
#include "floatarray.h"

#include <cstdio>
#include <vector>
#include <memory>

#define _IFT_ShellQd42_Name "shellqd42"
#define _IFT_ShellQd42_outputAtXY "outputatxy"
#define _IFT_ShellQd42_outputType "outputtype"
#define _IFT_ShellQd42_outputAtZ "outputatz"
#define _IFT_ShellQd42_outputCategory "outputcategory"

namespace oofem {
	class ShellQd42 : public QdShell
	{
		PlnStrssQd1Rot* membrane;
		PltQd4DKT* plate;

	private:
		IntArray positionVectorMemb, positionVectorPlate;
	
	protected:
		// This should be moved to the QdShell class (with casted plate and membrane objects).
		void computeGaussPoints() override;
	
	public:
		ShellQd42(int n, Domain* d);
		virtual ~ShellQd42() { delete[] membrane; delete[] plate; }

		void computeStiffnessMatrix(FloatMatrix& answer, MatResponseMode rMode, TimeStep* tStep) override;
		void computeStrainVectorAt(FloatArray& answer, double xi, double eta, TimeStep* tStep) override;
		// This method is used for stress calculation at Gauss points only for this element.
		void computeStressVector(FloatArray& answer, const FloatArray& strain, GaussPoint* gp, TimeStep* tStep) override;
		void computeStressVectorAtCentre(FloatArray& answer, TimeStep* tStep, const FloatArray& strain = 0) override;
		void getStressesTopBottom(FloatArray& answer, TimeStep* tStep) override;
		const char* giveClassName() const override { return "ShellQd42"; }
		int giveDefaultIntegrationRule() const override { return plate->giveDefaultIntegrationRule(); }
		void giveInternalForcesVector(FloatArray& answer, TimeStep* tStep, int useUpdatedGpRecord) override;
		FEInterpolation* giveInterpolation() const override { return plate->giveInterpolation(); }
		const char* giveInputRecordName() const override { return _IFT_ShellQd42_Name; }
		MaterialMode giveMaterialMode() override { return _PlaneStress; }
		void initializeFrom(InputRecord& ir) override;
		void setCrossSection(int csIndx) override;
		void updateLocalNumbering(EntityRenumberingFunctor& f) override;

		void computeBodyLoadVectorAt(FloatArray& answer, Load* forLoad, TimeStep* tStep, ValueModeType mode) override;
		void computeBoundarySurfaceLoadVector(FloatArray& answer, BoundaryLoad* load, int boundary, CharType type, ValueModeType mode, TimeStep* tStep, bool global = true) override;
		int computeLoadGToLRotationMtrx(FloatMatrix& answer) override { return membrane->computeLoadGToLRotationMtrx(answer); }
		int computeLoadLSToLRotationMatrix(FloatMatrix& answer, int iSurf, GaussPoint* gp) override { return 0; }
		void computeSurfaceNMatrix(FloatMatrix& answer, int boundaryID, const FloatArray& lcoords) override;

		// This method should never be called.
		void computeBmatrixAt(double xi, double eta, FloatMatrix& answer) override { }
		// This method should never be called.
		void computeConstitutiveMatrixAt(FloatMatrix& answer, MatResponseMode rMode, GaussPoint* gp, TimeStep* tStep) override {};
	};
}
#endif