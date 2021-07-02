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
 
#ifndef qdelement_h
#define qdelement_h

#include "femcmpnn.h"
#include "error.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "sm/Elements/nlstructuralelement.h"

#include <cstdio>
#include <vector>
#include <memory>

#define _IFT_QdElement_Name "qdelement"

namespace oofem {
	class QdElement : public NLStructuralElement
	{
	public:
		QdElement(int n, Domain* d);
		virtual ~QdElement() {}

		void computeBmatrixAt(GaussPoint* gp, FloatMatrix& answer,
			int lowerIndx = 1, int upperIndx = ALL_STRAINS) override { }
		void computeStressVector(FloatArray& answer, const FloatArray& strain, GaussPoint* gp, TimeStep* tStep) override { }
		void computeConstitutiveMatrixAt(FloatMatrix& answer,
			MatResponseMode rMode, GaussPoint* gp,
			TimeStep* tStep) override { }
		const char* giveClassName() const override { return "QdElement"; }
		const char* giveInputRecordName() const override { return _IFT_QdElement_Name; }
	};
}
#endif