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

#include "classfactory.h"
#include "gaussintegrationrule.h"
#include "gausspoint.h"
#include "load.h"

namespace oofem {
REGISTER_Element(ShellQd42);

ShellQd42 :: ShellQd42(int n, Domain* aDomain) : QdShell(n, aDomain)
{
    numberOfDofMans = 4;
    numberOfGaussPoints = 4;

    membrane = new PlnStrssQd1Rot(n, aDomain);
    plate = new PltQd4DKT(n, aDomain);
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
/*
void
ShellQd42::computeStressVector(FloatArray& answer, const FloatArray& strain, GaussPoint* gp, TimeStep* tStep) {
    FloatArray membraneStrains, membraneStresses, plateStrains, plateStresses;
    // IT SEEMS THAT THIS CANNOT WORK FOR SHELL ELEMENTS. HOW DO YOU CALCULATE COMBINED STRESSES AT ELEMENT'S CENTRE WITH THIS METHOD? 
    switch (outputCategory) {
    case OutputCategory::Membrane:
        membrane->computeStrainVectorAt(membraneStrains, gp->giveNaturalCoordinate(1), gp->giveNaturalCoordinate(2), tStep);
        membrane->computeStressVector(membraneStresses, membraneStrains, gp, tStep);
        break;
    case OutputCategory::Plate:
        plate->computeStrainVectorAt(plateStrains, gp->giveNaturalCoordinate(1), gp->giveNaturalCoordinate(2), tStep);
        plate->computeStressVector(plateStresses, plateStrains, gp, tStep);
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
}*/

void
ShellQd42::computeStiffnessMatrix(FloatMatrix& answer, MatResponseMode rMode, TimeStep* tStep)
{
    StructuralCrossSection* cs = this->giveStructuralCrossSection();

    FloatMatrix stiffMatPlate;
    plate->computeStiffnessMatrix(stiffMatPlate, rMode, tStep);
    FloatMatrix stiffMatMembrane;
    membrane->computeStiffnessMatrix(stiffMatMembrane, rMode, tStep);

    answer.clear();

    // Compute stiffness matrix of the shell element by manipulating relevant terms of its constituent parts' stiffness matrices.
    if (integrationRulesArray.size() == 1) {

        answer.resize(24, 24);

        answer = stiffMatMembrane;
        answer.add(stiffMatPlate);
    }
    else
        OOFEM_ERROR("It is not allowed to define more than 1 integration rule for this element.");
}

void
ShellQd42::initializeFrom(InputRecord& ir)
{
    StructuralElement::initializeFrom(ir);
    plate->initializeFrom(ir);
    membrane->initializeFrom(ir);

    int outputAtXYTemp, outputTypeTemp, outputAtZTemp, outputCategoryTemp;
    outputAtXYTemp = outputTypeTemp = outputAtZTemp = outputCategoryTemp = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, outputAtXYTemp, _IFT_ShellQd41_outputAtXY);
    IR_GIVE_OPTIONAL_FIELD(ir, outputTypeTemp, _IFT_ShellQd41_outputType);
    IR_GIVE_OPTIONAL_FIELD(ir, outputCategoryTemp, _IFT_ShellQd41_outputCategory);
    IR_GIVE_OPTIONAL_FIELD(ir, outputAtZ, _IFT_ShellQd41_outputAtZ);

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

    // Initialise output options for the membrane and the plate part.
    membrane->outputAtXY = outputAtXY;
    membrane->outputType = outputType;
    plate->outputAtXY = outputAtXY;
    plate->outputType = outputType;
    plate->outputAtZ = outputAtZ;
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