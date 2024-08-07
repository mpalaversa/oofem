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
 
#ifndef qdplate_h
#define qdplate_h

#include "sm/Elements/qdelement.h"
#include "femcmpnn.h"
#include "error.h"
#include "floatmatrix.h"
#include "floatarray.h"

#include <cstdio>
#include <vector>
#include <memory>

namespace oofem {
	class QdPlate : public QdElement
	{
	protected:
		double outputAtZ;

		void computeGaussPoints() override;

	public:
		QdPlate(int n, Domain* d);
		virtual ~QdPlate() = default;

		void computeBodyLoadVectorAt(FloatArray& answer, Load* forLoad, TimeStep* tStep, ValueModeType mode) override;
		virtual void computeCurvaturesAt(FloatArray& answer, double xi, double eta, TimeStep* tStep) = 0;
		void computeStressVectorAtCentre(FloatArray& answer, TimeStep* tStep, const FloatArray& strain = 0) override;
		double getOutputLocationInZ() { return outputAtZ; }
		void getStressesTopBottom(FloatArray& answer, TimeStep* tStep) { }
		void giveDofManDofIDMask(int inode, IntArray&) const override;
		void giveSurfaceDofMapping(IntArray& answer, int iSurf) const override;
	};
}
#endif