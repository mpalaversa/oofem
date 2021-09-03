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
 
#include "sm/Elements/PlaneStress/qdmembrane.h"

#include "classfactory.h"
#include "gaussintegrationrule.h"
#include "gausspoint.h"
#include "load.h"
#include "sm/Materials/structuralms.h"
#include "sm/CrossSections/structuralcrosssection.h"

namespace oofem {
QdMembrane::QdMembrane(int n, Domain* aDomain) : QdElement(n, aDomain)
{

}

void
QdMembrane::computeBodyLoadVectorAt(FloatArray& answer, Load* forLoad, TimeStep* tStep, ValueModeType mode) {
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
    localAccelerationVector.resize(2);
    localAccelerationVector.at(1) = acceleration.at(1);
    localAccelerationVector.at(2) = acceleration.at(2);
    
    FloatArray NMatrixTemp;
    FloatMatrix NMatrix;
    for (GaussPoint* gp : *this->giveDefaultIntegrationRulePtr()) {
        giveInterpolation()->evalN(NMatrixTemp, gp->giveSubPatchCoordinates(), *giveCellGeometryWrapper());
        NMatrix.beNMatrixOf(NMatrixTemp, 2);
        dV = computeSurfaceVolumeAround(gp, 1) * this->giveCrossSection()->give(CS_Thickness, gp);
        density = this->giveCrossSection()->give('d', gp);
        distributedAcceleration.beTProductOf(NMatrix, localAccelerationVector);
        answer.add(dV * density, distributedAcceleration);
    }
}

void
QdMembrane::computeConstitutiveMatrixAt(FloatMatrix& answer, MatResponseMode rMode, GaussPoint* gp, TimeStep* tStep)
{
    answer = this->giveStructuralCrossSection()->giveStiffnessMatrix_PlaneStress(rMode, gp, tStep);
}

void
QdMembrane::computeGaussPoints()
{
    // Sets up the integration rule array which contains all the Gauss points
    // Default: create one integration rule
    if (integrationRulesArray.size() == 0) {
        integrationRulesArray.resize(1);
        integrationRulesArray[0] = std::make_unique<GaussIntegrationRule>(1, this, 1, 3);
        this->giveCrossSection()->setupIntegrationPoints(*integrationRulesArray[0], this->numberOfGaussPoints, this);
    }
}

void
QdMembrane::computeStrainVectorAt(FloatArray& answer, double xi, double eta, TimeStep* tStep) {
    FloatMatrix b;
    FloatArray u;

    switch (outputType) {
    case OutputType::Standard:
        this->computeVectorOf(VM_Total, tStep, u);
        /*
        if (initialDisplacements) {
            u.subtract(*initialDisplacements);
        }*/
        computeBmatrixAt(xi, eta, b);
        answer.beProductOf(b, u);
        break;
    case OutputType::Principal:
        OOFEM_ERROR("Not yet implemented for a %s element.", giveClassName());
        break;
    case OutputType::VM:
        OOFEM_ERROR("Not yet implemented for a %s element.", giveClassName());
        break;
    case OutputType::All:
        OOFEM_ERROR("Not yet implemented for a %s element.", giveClassName());
        break;
    default:
        OOFEM_ERROR("Something went wrong. An unknown output category requested for element %d.", giveGlobalNumber());
        break;
    }
}

void
QdMembrane::computeStressVector(FloatArray& answer, const FloatArray& strain, GaussPoint* gp, TimeStep* tStep) {
    answer = this->giveStructuralCrossSection()->giveRealStress_PlaneStress(strain, gp, tStep);
}

void
QdMembrane::computeStressVectorAtCentre(FloatArray& answer, TimeStep* tStep, const FloatArray& strain) {
    computeStressVector(answer, strain, this->giveIntegrationRulesArray()[0]->getIntegrationPoint(0), tStep);
}

double
QdMembrane::computeVolumeAround(GaussPoint* gp)
{
    // Computes the volume element dV associated with the given gp.

    double weight = gp->giveWeight();
    const FloatArray& lCoords = gp->giveNaturalCoordinates(); // local/natural coords of the gp (parent domain)
    double detJ = fabs(this->giveInterpolation()->giveTransformationJacobian(lCoords, *this->giveCellGeometryWrapper()));
    double thickness = this->giveCrossSection()->give(CS_Thickness, gp); // the cross section keeps track of the thickness

    return detJ * thickness * weight; // dV
}

bool
QdMembrane::giveRotationMatrix(FloatMatrix& answer)
{
    bool is_GtoL, is_NtoG;
    FloatMatrix GtoL, NtoG;
    IntArray nodes;
    nodes.enumerate(this->giveNumberOfDofManagers());

    is_GtoL = this->computeGtoLRotationMatrix(GtoL);
    is_NtoG = this->computeDofTransformationMatrix(NtoG, nodes, true);

    if (is_GtoL && NtoG.isNotEmpty()) {
        answer.beProductOf(GtoL, NtoG);
    }
    else if (is_GtoL) {
        answer = GtoL;
    }
    else if (is_NtoG) {
        answer = NtoG;
    }
    else {
        answer.clear();
        return false;
    }
    return true;
}

void
QdMembrane::giveSurfaceDofMapping(IntArray& answer, int iSurf) const
{
    if (iSurf == 1 || iSurf == 2) {
        answer.enumerate(8);
    }
    else {
        OOFEM_ERROR("wrong surface number");
    }
}

void
QdMembrane::postInitialize()
{
    // Element must be created before giveNumberOfNodes can be called
    StructuralElement::postInitialize();
    this->numberOfDofMans = this->giveNumberOfNodes();
}

void
QdMembrane::printOutputAt(FILE* file, TimeStep* tStep)
{
    FloatArray v;
    fprintf(file, "element %d (%8d):\n", this->giveLabel(), number);
    int nPointsXY{ 0 };
    switch (outputAtXY)
    {
    case OutputLocationXY::GaussPoints:
        nPointsXY = giveIntegrationRulesArray()[0]->giveNumberOfIntegrationPoints();
        break;
    case OutputLocationXY::Centre:
        nPointsXY = 1;
        break;
    case OutputLocationXY::Corners:
        nPointsXY = giveNumberOfDofManagers();
        break;
    case OutputLocationXY::All:
        nPointsXY = giveIntegrationRulesArray()[0]->giveNumberOfIntegrationPoints() + giveNumberOfDofManagers() + 1;
        break;
    default:
        OOFEM_ERROR("Something went wrong. The following options for output location at XY plane can be selected for a membrane element: 'All', 'Centroid', 'Gauss points' and 'Corners'.");
        break;
    }

    StructuralMaterialStatus* ms;
    GaussPoint* gp;
    for (int i = 0; i < nPointsXY; i++) {
        switch (outputAtXY)
        {
        case OutputLocationXY::GaussPoints:
            fprintf(file, "  GP %d :", i + 1);
            gp = integrationRulesArray[0]->getIntegrationPoint(i);
            ms = static_cast<StructuralMaterialStatus*>(gp->giveMaterialStatus());

            /*
            fprintf(file, "  forces     ");
            for (auto& val : this->giveMidplaneIPValue(i, IST_ShellForceTensor, tStep)) {
                fprintf(file, " %.4e", val);
            }

            fprintf(file, "\n          moments    ");
            for (auto& val : this->giveMidplaneIPValue(i, IST_ShellMomentTensor, tStep)) {
                fprintf(file, " %.4e", val);
            }

            fprintf(file, "\n          strains    ");
            for (auto& val : ms->giveStrainVector()) {
                fprintf(file, " %.4e", val);
            }

            fprintf(file, "\n          curvatures ");
            for (auto& val : this->giveMidplaneIPValue(i, IST_CurvatureTensor, tStep)) {
                fprintf(file, " %.4e", val);
            }*/

            fprintf(file, "\n          strains    ");
            for (auto& val : ms->giveStrainVector())
                fprintf(file, " %.4e", val);

            fprintf(file, "\n          stresses    ");
            for (auto& val : ms->giveStressVector())
                fprintf(file, " %.4e", val);

            fprintf(file, "\n");
            break;
        case OutputLocationXY::Centre:
            fprintf(file, "  Centre");
            gp = integrationRulesArray[0]->getIntegrationPoint(0);
            ms = static_cast<StructuralMaterialStatus*>(gp->giveMaterialStatus());

            /*
            fprintf(file, "  forces     ");
            for (auto& val : this->giveMidplaneIPValue(i, IST_ShellForceTensor, tStep)) {
                fprintf(file, " %.4e", val);
            }

            fprintf(file, "\n          moments    ");
            for (auto& val : this->giveMidplaneIPValue(i, IST_ShellMomentTensor, tStep)) {
                fprintf(file, " %.4e", val);
            }

            fprintf(file, "\n          strains    ");
            for (auto& val : ms->giveStrainVector()) {
                fprintf(file, " %.4e", val);
            }

            fprintf(file, "\n          curvatures ");
            for (auto& val : this->giveMidplaneIPValue(i, IST_CurvatureTensor, tStep)) {
                fprintf(file, " %.4e", val);
            }*/

            fprintf(file, "\n          strains    ");
            for (auto& val : ms->giveStrainVector())
                fprintf(file, " %.4e", val);

            fprintf(file, "\n          stresses    ");
            for (auto& val : ms->giveStressVector())
                fprintf(file, " %.4e", val);

            fprintf(file, "\n");
            break;
        case OutputLocationXY::Corners:
            fprintf(file, "  Node %d :", i + 1);
            break;
        case OutputLocationXY::All:
            if (i < 4)
                fprintf(file, "  GP %d :", i + 1);
            else if (i >= 4 && i < 8)
                fprintf(file, "  Node %d :", i - 3);
            else
                fprintf(file, "  Centroid :");
            break;
        default:
            OOFEM_ERROR("Something went wrong. The following options for output location at XY plane can be selected for a membrane element: 'All', 'Centroid', 'Gauss points' and 'Corners'.");
            break;
        }
    }
}
}