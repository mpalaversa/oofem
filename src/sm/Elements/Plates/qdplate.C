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
#include "load.h"
#include "sm/CrossSections/structuralcrosssection.h"

namespace oofem {
    QdPlate::QdPlate(int n, Domain* aDomain) : QdElement(n, aDomain)
    {
        outputAtZ = 0.0;
    }

    void
    QdPlate::computeBodyLoadVectorAt(FloatArray& answer, Load* forLoad, TimeStep* tStep, ValueModeType mode) {
        double density, dV;
        FloatArray acceleration, distributedAcceleration;
        FloatMatrix T;

        if ((forLoad->giveBCGeoType() != BodyLoadBGT) || (forLoad->giveBCValType() != ForceLoadBVT)) {
            OOFEM_ERROR("Unknown load type.");
        }

        forLoad->computeComponentArrayAt(acceleration, tStep, mode);
        // The acceleration vector is given in the global coordinate system and needs to be transformed to the element local c.s.
        if (this->computeLoadGToLRotationMtrx(T))
            acceleration.rotatedWith(T, 'n');

        // A general membrane element can take up two forces (only forces in plane of the element).
        FloatArray localAccelerationVector;
        localAccelerationVector.resize(3);
        int j = 3;
        for (int i = 1; i <= 3; i++) {
            localAccelerationVector.at(i) = acceleration.at(j);
            j++;
        }

        FloatArray NMatrixTemp;
        FloatMatrix NMatrix;
        for (GaussPoint* gp : *this->giveDefaultIntegrationRulePtr()) {
            giveInterpolation()->evalN(NMatrixTemp, gp->giveSubPatchCoordinates(), *giveCellGeometryWrapper());
            NMatrix.beNMatrixOf(NMatrixTemp, 3);
            dV = computeSurfaceVolumeAround(gp, 1) * this->giveCrossSection()->give(CS_Thickness, gp);
            density = this->giveCrossSection()->give('d', gp);
            distributedAcceleration.beTProductOf(NMatrix, localAccelerationVector);
            answer.add(dV * density, distributedAcceleration);
        }
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
    QdPlate::computeStressVectorAtCentre(FloatArray& answer, TimeStep* tStep, const FloatArray& strain) {
        computeStressVector(answer, strain, this->giveIntegrationRulesArray()[0]->getIntegrationPoint(0), tStep);
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