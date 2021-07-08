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
    
    // Group N,x and N,y into respective matrices.
    FloatMatrix NXDiff;
    NXDiff.resize(1, 4);
    FloatMatrix NYDiff;
    NYDiff.resize(1, 4);
    int k{ 9 };
    for (int i = 1; i <= 4; i++) {
        NXDiff.at(1, i) = BMatrix8.at(1, k);
        NYDiff.at(1, i) = BMatrix8.at(2, k + 1);
        k += 2;
    }

    // Fetch nodes coordinates in element's coord. system.
    std::vector< FloatArray > locCoords = giveNodeCoordinates();

    FloatMatrix transformedRows;
    transformedRows.resize(3, 4);
    int j{ 1 };
    k = 1;
    double xPlus{0.0};
    double yPlus{ 0.0 };
    double xMinus{ 0.0 };
    double yMinus{ 0.0 };
    double x{ 0.0 };
    double y{ 0.0 };
    for (int i = 1; i <= 4; i++) {
        j = i - 1;
        if (j < 1)
            j = 4;
        k = i + 1;
        if (k > 4)
            k = 1;
        
        x = locCoords.at(i - 1).at(1);
        y = locCoords.at(i - 1).at(2);
        xPlus = locCoords.at(k - 1).at(1);
        yPlus = locCoords.at(k - 1).at(2);
        xMinus = locCoords.at(j - 1).at(1);
        yMinus = locCoords.at(j - 1).at(2);
        // Get the derivatives that correspond to the rotational DOFs at vertices.
        transformedRows.at(1, i) = (y - yPlus) * NXDiff.at(1, i) + (y - yMinus) * NXDiff.at(1, j);
        transformedRows.at(2, i) = (xPlus - x) * NYDiff.at(1, i) + (xMinus - x) * NYDiff.at(1, j);
        transformedRows.at(3, i) = (xPlus - x) * NXDiff.at(1, i) + (xMinus - x) * NXDiff.at(1, j) + (y - yPlus) * NYDiff.at(1, i) + (y - yMinus) * NYDiff.at(1, j);
    }

    // Form B-matrix of size 3x12 for the element with rotational DOFs.
    answer.resize(3, 12);
    k = 9;
    int l{ 1 };
    for (int i = 1; i <= 3; i++) {
        l = 1;
        k = 9;
        for (int j = 1; j <= 12; j++) {
            if (j % 3 == 0) {
                answer.at(i, j) = transformedRows.at(i, l);
                l++;
            }
            else {
                answer.at(i, j) = BMatrix8.at(i, k);
                k++;
            }
        }
    }
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
PlnStrssQd1Rot::giveDofManDofIDMask(int inode, IntArray& answer) const
{
	answer = { D_u, D_v, R_w }; // Does the plane R_w rotation transform into all three rotations in space?
}

FEInterpolation*
PlnStrssQd1Rot::giveInterpolation() const {
    return &interpolation;
}
}