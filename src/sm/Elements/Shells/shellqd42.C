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
 
#include "sm/Elements/Shells/shellqd42.h"

#include "boundaryload.h"
#include "classfactory.h"
#include "crosssection.h"
#include "gaussintegrationrule.h"
#include "gausspoint.h"
#include "load.h"
#include "sm/CrossSections/structuralcrosssection.h"

namespace oofem {
REGISTER_Element(ShellQd42);

ShellQd42 :: ShellQd42(int n, Domain* aDomain) : QdShell(n, aDomain)
{
    numberOfDofMans = 4;
    numberOfGaussPoints = 4;

    membrane = new PlnStrssQd1Rot(n, aDomain);
    plate = new PltQd4DKT(n, aDomain);

    positionVectorMemb.resize(12);
    positionVectorMemb = { 1, 2, 6, 7, 8, 12, 13, 14, 18, 19, 20, 24 };
    positionVectorPlate.resize(12);
    positionVectorPlate = { 3, 4, 5, 9, 10, 11, 15, 16, 17, 21, 22, 23 };
}

void
ShellQd42::computeBodyLoadVectorAt(FloatArray& answer, Load* forLoad, TimeStep* tStep, ValueModeType mode) {
    FloatArray loadFromMembrane, loadFromPlate;
    FloatMatrix rotMat;

    // Calculate distribution of loads for the membrane and plate part. Note that they are all given in the element coord. sys.
    membrane->computeBodyLoadVectorAt(loadFromMembrane, forLoad, tStep, mode);
    plate->computeBodyLoadVectorAt(loadFromPlate, forLoad, tStep, mode);

    // Assemble shell load vector from the membrane and plate parts.
    answer.resize(24);
    int j = 1;
    for (int i = 1; i <= 12; i += 3) {
        answer.at(j) = loadFromMembrane.at(i);
        answer.at(j + 1) = loadFromMembrane.at(i + 1);
        answer.at(j + 2) = loadFromPlate.at(i);
        answer.at(j + 3) = loadFromPlate.at(i + 1);
        answer.at(j + 4) = loadFromPlate.at(i + 2);
        answer.at(j + 5) = loadFromMembrane.at(i + 2);
        j += 6;
    }
}

void
ShellQd42::computeBoundarySurfaceLoadVector(FloatArray& answer, BoundaryLoad* load, int boundary, CharType type, ValueModeType mode, TimeStep* tStep, bool global)
{
    FloatArray loadFromMembrane, loadFromPlate;
    // Calculate loads from the membrane and plate part.
    membrane->computeBoundarySurfaceLoadVector(loadFromMembrane, load, boundary, type, mode, tStep, global);
    plate->computeBoundarySurfaceLoadVector(loadFromPlate, load, boundary, type, mode, tStep, global);

    // Add the membrane and the plate contribution.
    answer.add(loadFromMembrane);
    answer.add(loadFromPlate);
}

void
ShellQd42::computeGaussPoints()
{
    // Sets up the integration rule array which contains all the Gauss points
    // Default: create one integration rule
    
    if (integrationRulesArray.size() == 0) {
        integrationRulesArray.resize(1);
        integrationRulesArray[0] = std::make_unique<GaussIntegrationRule>(1, this, 1, 3);
        this->giveCrossSection()->setupIntegrationPoints(*integrationRulesArray[0], this->numberOfGaussPoints, this);
    }

    membrane->computeGaussPoints();
    plate->computeGaussPoints();
}

void
ShellQd42::computeStiffnessMatrix(FloatMatrix& answer, MatResponseMode rMode, TimeStep* tStep)
{
    FloatMatrix stiffMatPlate;
    plate->computeStiffnessMatrix(stiffMatPlate, rMode, tStep);
    FloatMatrix stiffMatMembrane;
    membrane->computeStiffnessMatrix(stiffMatMembrane, rMode, tStep);
    
    answer.clear();
    answer.resize(24, 24);

    int startCol{ 1 };
    for (int i = 1; i <= 12; i++) {
        for (int j = startCol; j <= 12; j++) {
            answer.at(positionVectorMemb.at(i), positionVectorMemb.at(j)) = stiffMatMembrane.at(i, j);
            answer.at(positionVectorPlate.at(i), positionVectorPlate.at(j)) = stiffMatPlate.at(i, j);
        }
        startCol++;
    }

    answer.symmetrized();
}

void
ShellQd42::computeStrainVectorAt(FloatArray& answer, double xi, double eta, TimeStep* tStep) {
    FloatArray membraneStrains, plateStrains;

    switch (outputCategory) {
    case OutputCategory::Membrane:
        membrane->computeStrainVectorAt(answer, xi, eta, tStep);
        break;
    case OutputCategory::Plate:
        plate->computeStrainVectorAt(answer, xi, eta, tStep);
        break;
    case OutputCategory::Combined:
        membrane->computeStrainVectorAt(membraneStrains, xi, eta, tStep);
        plate->computeStrainVectorAt(plateStrains, xi, eta, tStep);
        membraneStrains.add(plateStrains);
        answer.resize(3);
        answer = membraneStrains;
        break;
    default:
        OOFEM_ERROR("Something went wrong. An unknown output category requested for the element %d.", giveGlobalNumber());
        break;
    }
}

void
ShellQd42::computeStressVector(FloatArray& answer, const FloatArray& strain, GaussPoint* gp, TimeStep* tStep) {
    FloatArray membraneStrains, membraneStresses, plateStrains, plateStresses;
    switch (outputCategory) {
    case OutputCategory::Membrane:
        membrane->computeStrainVectorAt(membraneStrains, gp->giveNaturalCoordinate(1), gp->giveNaturalCoordinate(2), tStep);
        answer.resize(3);
        membrane->computeStressVector(answer, membraneStrains, gp, tStep);
        break;
    case OutputCategory::Plate:
        plate->computeStrainVectorAt(plateStrains, gp->giveNaturalCoordinate(1), gp->giveNaturalCoordinate(2), tStep);
        answer.resize(3);
        plate->computeStressVector(answer, plateStrains, gp, tStep);
        break;
    case OutputCategory::Combined:
        membrane->computeStrainVectorAt(membraneStrains, gp->giveNaturalCoordinate(1), gp->giveNaturalCoordinate(2), tStep);
        plate->computeStrainVectorAt(plateStrains, gp->giveNaturalCoordinate(1), gp->giveNaturalCoordinate(2), tStep);
        answer = giveStructuralCrossSection()->giveRealStress_Shell(membraneStrains, plateStrains, gp, tStep);
        break;
    default:
        OOFEM_ERROR("Something went wrong. An unknown output category requested for the element %d.", giveGlobalNumber());
        break;
    }
}

void
ShellQd42::computeStressVectorAtCentre(FloatArray& answer, TimeStep* tStep, const FloatArray& strain) {
    FloatArray membraneStrains, plateStrains, membraneStresses, plateStresses;

    switch (outputCategory) {
    case OutputCategory::Membrane:
        membrane->computeStressVector(answer, strain, this->giveIntegrationRulesArray()[0]->getIntegrationPoint(0), tStep);
        break;
    case OutputCategory::Plate:
        plate->computeStressVector(answer, strain, this->giveIntegrationRulesArray()[0]->getIntegrationPoint(0), tStep);
        break;
    case OutputCategory::Combined:
        membrane->computeStrainVectorAt(membraneStrains, 0.0, 0.0, tStep);
        plate->computeStrainVectorAt(plateStrains, 0.0, 0.0, tStep);
        answer = this->giveStructuralCrossSection()->giveRealStress_Shell(membraneStrains, plateStrains, this->giveIntegrationRulesArray()[0]->getIntegrationPoint(0), tStep);
        break;
    default:
        OOFEM_ERROR("Something went wrong. An unknown output category requested for element %d.", giveGlobalNumber());
        break;
    }
}

void
ShellQd42::computeSurfaceNMatrix(FloatMatrix& answer, int boundaryID, const FloatArray& lcoords)
{
    FloatArray n_vec;
    this->giveInterpolation()->boundarySurfaceEvalN(n_vec, boundaryID, lcoords, *giveCellGeometryWrapper());
    answer.beNMatrixOf(n_vec, 6);
}

void
ShellQd42::getStressesTopBottom(FloatArray& answer, TimeStep* tStep) {
    // Remove the following 4 lines of code when the method is considered generic.
    outputAtXY = OutputLocationXY::Centre;
    outputCategory = OutputCategory::Combined;
    outputType = OutputType::Standard;
    outputAtZ = this->giveStructuralCrossSection()->give(CS_Thickness, this->giveDefaultIntegrationRulePtr()->getIntegrationPoint(0)) / 2;

    // Initialise output options for the membrane and the plate part.
    membrane->outputAtXY = outputAtXY;
    membrane->outputType = outputType;
    plate->outputAtXY = outputAtXY;
    plate->outputType = outputType;
    plate->outputAtZ = outputAtZ;

    computeStressVectorAtCentre(answer, tStep);
}

void
ShellQd42::giveInternalForcesVector(FloatArray& answer, TimeStep* tStep, int useUpdatedGpRecord) {
    FloatArray internalForcesMemb, internalForcesPlate;

    membrane->giveInternalForcesVector(internalForcesMemb, tStep, 0);
    plate->giveInternalForcesVector(internalForcesPlate, tStep, 0);

    answer.resize(24);
    for (int i = 1; i <= 12; i++) {
        answer.at(positionVectorMemb.at(i)) = internalForcesMemb.at(i);
        answer.at(positionVectorPlate.at(i)) = internalForcesPlate.at(i);
    }
}

void
ShellQd42::initializeFrom(InputRecord& ir)
{
    StructuralElement::initializeFrom(ir);
    plate->initializeFrom(ir);
    membrane->initializeFrom(ir);

    int outputAtXYTemp, outputTypeTemp, outputAtZTemp, outputCategoryTemp;
    outputAtXYTemp = outputTypeTemp = outputAtZTemp = outputCategoryTemp = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, outputAtXYTemp, _IFT_ShellQd42_outputAtXY);
    IR_GIVE_OPTIONAL_FIELD(ir, outputTypeTemp, _IFT_ShellQd42_outputType);
    IR_GIVE_OPTIONAL_FIELD(ir, outputCategoryTemp, _IFT_ShellQd42_outputCategory);
    IR_GIVE_OPTIONAL_FIELD(ir, outputAtZ, _IFT_ShellQd42_outputAtZ);

    switch (outputAtXYTemp) {
    case 1:
        outputAtXY = OutputLocationXY::GaussPoints;
        break;
    case 2:
        outputAtXY = OutputLocationXY::Centre;
        break;
    case 3:
        outputAtXY = OutputLocationXY::Corners;
        break;
    case 4:
        outputAtXY = OutputLocationXY::All;
        break;
    default:
        outputAtXY = OutputLocationXY::GaussPoints;
        break;
    }
    switch (outputCategoryTemp) {
    case 1:
        outputCategory = OutputCategory::Membrane;
        break;
    case 2:
        outputCategory = OutputCategory::Plate;
        break;
    case 3:
        outputCategory = OutputCategory::Combined;
        break;
    case 4:
        outputCategory = OutputCategory::All;
        break;
    default:
        outputCategory = OutputCategory::Membrane;
        break;
    }
    switch (outputTypeTemp) {
    case 1:
        outputType = OutputType::Standard;
        break;
    case 2:
        outputType = OutputType::Principal;
        break;
    case 3:
        outputType = OutputType::VM;
        break;
    case 4:
        outputType = OutputType::All;
        break;
    default:
        outputType = OutputType::Standard;
        break;
    }

    int csClassTemp = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, csClassTemp, _IFT_ShellQd42_csClass);
    switch (csClassTemp) {
    case 1:
        csClass = CSClass::OOFEM;
        break;
    case 2:
        csClass = CSClass::Nastran;
        break;
    default:
        csClass = CSClass::OOFEM;
        break;
    }

    // Initialise output options for the membrane and the plate part.
    membrane->csClass = csClass;
    membrane->outputAtXY = outputAtXY;
    membrane->outputType = outputType;
    plate->csClass = csClass;
    plate->outputAtXY = outputAtXY;
    plate->outputType = outputType;
    plate->outputAtZ = outputAtZ;
}

void
ShellQd42::setCrossSection(int csIndx)
{
    StructuralElement::setCrossSection(csIndx);
    plate->setCrossSection(csIndx);
    membrane->setCrossSection(csIndx);
}

void
ShellQd42::updateLocalNumbering(EntityRenumberingFunctor& f)
{
    // Update numbering of the related DOFs for the membrane and plate part.
    membrane->updateLocalNumbering(f);
    plate->updateLocalNumbering(f);
    // Update numbering of the related DOFs for ShellQd41.
    for (auto& dnum : dofManArray)
        dnum = f(dnum, ERS_DofManager);
}
}