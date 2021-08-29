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
 
#include "sm/Elements/PlaneStress/plnstrssqd1rot.h"
#include "classfactory.h"
#include "fei2dquadquad.h"

namespace oofem {
REGISTER_Element(PlnStrssQd1Rot);
FEI2dQuadQuad PlnStrssQd1Rot::interpolation(1, 2);

PlnStrssQd1Rot :: PlnStrssQd1Rot(int n, Domain* aDomain) : QdMembrane(n, aDomain)
{
    numberOfDofMans = 4;
    numberOfGaussPoints = 4;
}

void
PlnStrssQd1Rot::computeBmatrixAt(double xi, double eta, FloatMatrix& answer)
{
    FloatArray naturalCoordinates;
    naturalCoordinates.resize(2);
    naturalCoordinates.at(1) = xi;
    naturalCoordinates.at(2) = eta;

    FEInterpolation* interp = this->giveInterpolation();
    FloatMatrix dNdx;
    interp->evaldNdx(dNdx, naturalCoordinates, *this->giveCellGeometryWrapper());

    // B-matrix of the 8-node plane stress quadrilateral.
    FloatMatrix BMatrix8;
    BMatrix8.resize(3, dNdx.giveNumberOfRows() * 2);

    for (int i = 1; i <= dNdx.giveNumberOfRows(); i++) {
        BMatrix8.at(1, i * 2 - 1) = dNdx.at(i, 1);
        BMatrix8.at(2, i * 2 - 0) = dNdx.at(i, 2);

        BMatrix8.at(3, 2 * i - 1) = dNdx.at(i, 2);
        BMatrix8.at(3, 2 * i - 0) = dNdx.at(i, 1);
    }

    // Fetch nodes coordinates in element's coord. system.
    std::vector< FloatArray > locCoords = giveNodeCoordinates();

    // Matrix that transforms the displacement vector of the 8-node quad to the displacement vector of the 4-node quad with rotational DOFs.
    FloatMatrix TMatrix;
    TMatrix.resize(16, 12);
    // u and v at the vertices of the 8-node element correspond to those of the 4-node element.
    TMatrix.at(1, 1) = TMatrix.at(2, 2) = TMatrix.at(3, 4) = TMatrix.at(4, 5) = TMatrix.at(5, 7) = TMatrix.at(6, 8) = TMatrix.at(7, 10) = TMatrix.at(8, 11) = 1;
    
    // Set transformation coefficients for the DOFs at the midside nodes.
    int midsideNode{ 5 };
    IntArray vertexNodes;
    for (int i = 9; i <= 16; i+=2) {
        getVertexNodes(vertexNodes, midsideNode);
        // Set values for transformation of u at a midside node.
        TMatrix.at(i, (vertexNodes.at(1) - 1) * 3 + 1) = 0.5;
        TMatrix.at(i, (vertexNodes.at(1) - 1) * 3 + 3) = - (locCoords.at(vertexNodes.at(2) - 1)[1] - locCoords.at(vertexNodes.at(1) - 1)[1])/8;
        TMatrix.at(i, (vertexNodes.at(2) - 1) * 3 + 1) = 0.5;
        TMatrix.at(i, (vertexNodes.at(2) - 1) * 3 + 3) = (locCoords.at(vertexNodes.at(2) - 1)[1] - locCoords.at(vertexNodes.at(1) - 1)[1]) / 8;
        // Set values for transformation of v at a midside node.
        TMatrix.at(i + 1, (vertexNodes.at(1) - 1) * 3 + 2) = 0.5;
        TMatrix.at(i + 1, (vertexNodes.at(1) - 1) * 3 + 3) = (locCoords.at(vertexNodes.at(2) - 1)[0] - locCoords.at(vertexNodes.at(1) - 1)[0]) / 8;
        TMatrix.at(i + 1, (vertexNodes.at(2) - 1) * 3 + 2) = 0.5;
        TMatrix.at(i + 1, (vertexNodes.at(2) - 1) * 3 + 3) = - (locCoords.at(vertexNodes.at(2) - 1)[0] - locCoords.at(vertexNodes.at(1) - 1)[0]) / 8;
        midsideNode++;
    }

    answer.beProductOf(BMatrix8, TMatrix);
}

bool
PlnStrssQd1Rot::computeGtoLRotationMatrix(FloatMatrix& answer)
// Returns the rotation matrix of the receiver of the size [24,24]
// r(local) = T * r(global)
// for one node (r written transposed): {u,v,R_w} = T * {u,v,w,R_u,R_v,R_w}
{
    if (GtoLRotationMatrix == NULL)
        QdElement::computeGtoLRotationMatrix();

    answer.resize(12, 24);
    answer.zero();

    for (int i = 1; i <= 3; i++) {
        answer.at(1, i) = answer.at(4, i + 6) = answer.at(7, i + 12) = answer.at(10, i + 18) = GtoLRotationMatrix->at(1, i);
        answer.at(2, i) = answer.at(5, i + 6) = answer.at(8, i + 12) = answer.at(11, i + 18) = GtoLRotationMatrix->at(2, i);
        answer.at(3, i + 3) = answer.at(6, i + 9) = answer.at(9, i + 15) = answer.at(12, i + 21) = GtoLRotationMatrix->at(3, i);
    }

    return 1;
}

void
PlnStrssQd1Rot::getVertexNodes(IntArray &answer, int midsideNode)
// Returns nodes associated with the element's vertices next to the given midside node of a corresponding 8-node quad.
{
    answer.resize(2);

    int smallerVertex = midsideNode - 4;
    int greaterVertex = smallerVertex + 1;
    if (greaterVertex == 5)
        greaterVertex = 1;

    answer.at(1) = smallerVertex;
    answer.at(2) = greaterVertex;
}

void
PlnStrssQd1Rot::giveDofManDofIDMask(int inode, IntArray& answer) const
{
	answer = { D_u, D_v, D_w, R_u, R_v, R_w };
}

FEInterpolation*
PlnStrssQd1Rot::giveInterpolation() const {
    return &interpolation;
}

void
PlnStrssQd1Rot::initializeFrom(InputRecord& ir)
{
    StructuralElement::initializeFrom(ir);

    int outputAtXYTemp, outputTypeTemp;
    outputAtXYTemp = outputTypeTemp = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, outputAtXYTemp, _IFT_PlnStrssQd1Rot_outputAtXY);
    IR_GIVE_OPTIONAL_FIELD(ir, outputTypeTemp, _IFT_PlnStrssQd1Rot_outputType);

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
    IR_GIVE_OPTIONAL_FIELD(ir, csClassTemp, _IFT_PlnStrssQd1Rot_csClass);
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
}
}