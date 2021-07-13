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
 
#ifndef qdmembrane_h
#define qdmembrane_h


#include "sm/Elements/qdelement.h"
#include "femcmpnn.h"
#include "error.h"
#include "floatmatrix.h"
#include "floatarray.h"

#include <cstdio>
#include <vector>
#include <memory>

namespace oofem {
	class QdMembrane : public QdElement
	{
	public:
		QdMembrane(int n, Domain* d);

		void computeConstitutiveMatrixAt(FloatMatrix& answer, MatResponseMode rMode, GaussPoint* gp, TimeStep* tStep) override;
		virtual bool computeGtoLRotationMatrix(FloatMatrix& answer) override = 0;
		void computeStrainVectorAt(FloatArray& answer, double xi, double eta, TimeStep* tStep) override;
		void computeStressVector(FloatArray& answer, const FloatArray& strain, GaussPoint* gp, TimeStep* tStep) override;
		virtual double computeVolumeAround(GaussPoint* gp) override;
		bool giveRotationMatrix(FloatMatrix& answer) override;
		virtual void giveSurfaceDofMapping(IntArray& answer, int iSurf) const override;
		void postInitialize() override;
		void printOutputAt(FILE* file, TimeStep* tStep) override;

	protected:
		void computeGaussPoints() override;
	};
}
#endif