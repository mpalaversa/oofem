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
 
#include "sm/Elements/Plates/qdplate.h"

#include "classfactory.h"
#include "gaussintegrationrule.h"
#include "sm/CrossSections/structuralcrosssection.h"

namespace oofem {
    QdPlate::QdPlate(int n, Domain* aDomain) : QdElement(n, aDomain)
    {
        outputAtZ = 0.0;
    }

    void
    QdPlate::computeGaussPoints()
        // Sets up the array containing the four Gauss points of the receiver.
    {
        if (integrationRulesArray.size() == 0) {
            integrationRulesArray.resize(1);
            integrationRulesArray[0] = std::make_unique<GaussIntegrationRule>(1, this, 1, 5);
            this->giveCrossSection()->setupIntegrationPoints(*integrationRulesArray[0], numberOfGaussPoints, this);
        }
    }

    void
    QdPlate::giveDofManDofIDMask(int inode, IntArray& answer) const
    {
        answer = { D_u, D_v, D_w, R_u, R_v, R_w };
    }

    void
    QdPlate::giveSurfaceDofMapping(IntArray& answer, int iSurf) const
    {
        if (iSurf == 1 || iSurf == 2) {
            answer.enumerate(24);
        }
        else {
            OOFEM_ERROR("wrong surface number");
        }
    }
}