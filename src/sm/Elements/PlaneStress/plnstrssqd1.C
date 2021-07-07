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
 
#include "sm/Elements/PlaneStress/plnstrssqd1.h"
#include "classfactory.h"
#include "fei2dquadlin.h"

namespace oofem {
REGISTER_Element(PlnStrssQd1);
FEI2dQuadLin PlnStrssQd1::interpolation(1, 2);

PlnStrssQd1 :: PlnStrssQd1(int n, Domain* aDomain) : QdMembrane(n, aDomain)
{
    OutputLocationXY outputAtXY = OutputLocationXY::GaussPoints;
    OutputType outputType = OutputType::Standard;

    numberOfDofMans = 4;
    numberOfGaussPoints = 4;
}

void
PlnStrssQd1::computeBmatrixAt(double xi, double eta, FloatMatrix& answer) {
    FloatArray naturalCoordinates;
    naturalCoordinates.resize(2);
    naturalCoordinates.at(1) = xi;
    naturalCoordinates.at(2) = eta;

    FloatMatrix dnx;
    this->interpolation.evaldNdx(dnx, naturalCoordinates, *this->giveCellGeometryWrapper());

    answer.resize(3, 8);
    answer.zero();

    for (int i = 1; i <= 4; i++) {
        // Assembles terms in B matrix related to normal strains.
        answer.at(1, 2 * i - 1) = dnx.at(i, 1);
        answer.at(2, 2 * i - 0) = dnx.at(i, 2);
    }

    this->interpolation.evaldNdx(dnx, { 0., 0. }, *this->giveCellGeometryWrapper());

    for (int i = 1; i <= 4; i++) {
        // Assembles terms in B matrix related to shear strains.
        answer.at(3, 2 * i - 1) = dnx.at(i, 2);
        answer.at(3, 2 * i - 0) = dnx.at(i, 1);
    }
}

bool
PlnStrssQd1::computeGtoLRotationMatrix(FloatMatrix& answer)
// Returns the rotation matrix of the receiver of the size [8,12]
// r(local) = T * r(global)
// for one node (r written transposed): {u,v} = T * {u,v,w}
{
    // test if pereviously computed
    if (GtoLRotationMatrix == NULL)
        QdElement::computeGtoLRotationMatrix();

    answer.resize(8, 12);
    answer.zero();

    for (int i = 1; i <= 3; i++) {
        answer.at(1, i) = answer.at(3, i + 3) = answer.at(5, i + 6) = answer.at(7, i + 9) = GtoLRotationMatrix->at(1, i);
        answer.at(2, i) = answer.at(4, i + 3) = answer.at(6, i + 6) = answer.at(8, i + 9) = GtoLRotationMatrix->at(2, i);
    }

    return 1;
}

void
PlnStrssQd1::giveDofManDofIDMask(int inode, IntArray& answer) const
{
    answer = { D_u, D_v, D_w };
}

FEInterpolation*
PlnStrssQd1::giveInterpolation() const {
    return &interpolation;
}

bool
PlnStrssQd1::giveRotationMatrix(FloatMatrix& answer)
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
}