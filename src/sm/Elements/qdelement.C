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
 
#include <cmath>

#include "sm/Elements/qdelement.h"

#include "classfactory.h"
#include "gausspoint.h"

namespace oofem {
QdElement::QdElement(int n, Domain* aDomain) : NLStructuralElement(n, aDomain)
{
        outputAtXY = OutputLocationXY::GaussPoints;
        outputType = OutputType::Standard;
        csClass = CSClass::OOFEM;

        GtoLRotationMatrix = NULL;
        cellGeometryWrapper = NULL;
}

void
QdElement::computeBmatrixAt(GaussPoint* gp, FloatMatrix& answer, int lowerIndx, int upperIndx) {
    computeBmatrixAt(gp->giveNaturalCoordinate(1), gp->giveNaturalCoordinate(2), answer);
}

const FloatMatrix*
QdElement::computeGtoLRotationMatrix()
// Returns the translation matrix of the receiver of the size [3,3]
// coords(local) = T * coords(global)
//
// OOFEM's local coordinate system (described by vector triplet e1',e2',e3') is defined as follows (for Nastran's coordinate system description, see the header file):
//
// e1'    : [N2-N1]    Ni - means i - th node
// help   : [N3-N1]
// e3'    : e1' x help
// e2'    : e3' x e1'
{
    if (GtoLRotationMatrix == NULL) {
        FloatArray node1 = this->giveNode(1)->giveCoordinates();
        FloatArray node2 = this->giveNode(2)->giveCoordinates();
        FloatArray node3 = this->giveNode(3)->giveCoordinates();
        FloatArray node4 = this->giveNode(4)->giveCoordinates();
        FloatArray e1, e2, e3;

        switch (csClass)
        {
        case CSClass::Nastran:
        {
            // Vectors spanned bewteen nodes 1 and 2, 1 and 3 and 2 and 4 respectively.
            FloatArray v12, v13, v24;
            v12.beDifferenceOf(node2, node1);
            v13.beDifferenceOf(node3, node1);
            v24.beDifferenceOf(node4, node2);

            // Find angles between vectors v12 and v13 (beta) and between vectors v12 and v24 (gamma).
            double cosBeta{ 0.0 }, cosGamma{ 0.0 };
            FloatArray v12Temp;
            v12Temp.beScaled(-1.0, v12);
            cosBeta = (v12.dotProduct(v13)) / (v12.computeNorm() * v13.computeNorm());
            cosGamma = (v12Temp.dotProduct(v24)) / (v12.computeNorm() * v24.computeNorm());

            // Calculate angle between the x-axis and vC3.
            double alpha = (acos(cosBeta) + acos(cosGamma)) / 2;

            // Calculate vector that is coincident with the x-axis.
            FloatMatrix rotMat;
            rotMat.resize(3, 3);
            rotMat.at(1, 1) = cos(-alpha);
            rotMat.at(1, 2) = -sin(-alpha);
            rotMat.at(1, 3) = 0.0;
            rotMat.at(2, 1) = sin(-alpha);
            rotMat.at(2, 2) = cos(-alpha);
            rotMat.at(2, 3) = 0.0;
            rotMat.at(3, 1) = 0.0;
            rotMat.at(3, 2) = 0.0;
            rotMat.at(3, 3) = 1.0;

            FloatArray vX;
            vX.beProductOf(rotMat, v13);

            // Calculate vector coincident with the z-axis.
            FloatArray vZ;
            vX.resizeWithValues(3);
            vZ.beVectorProductOf(vX, v13);

            // Calculate vector coincident with th y-axis.
            FloatArray vY;
            vY.beVectorProductOf(vZ, vX);

            // Calculate unit vectors in the direction of the x, y and z axis.
            vX.normalize();
            e1.copySubVector(vX, 1);
            vY.normalize();
            e2.copySubVector(vY, 1);
            vZ.normalize();
            e3.copySubVector(vZ, 1);

            break;
        }
        default:
        {
            FloatArray help;

            // compute e1' = [N2-N1]  and  help = [N3-N1]
            e1.beDifferenceOf(node2, node1);
            help.beDifferenceOf(node3, node1);

            // let us normalize e1'
            e1.normalize();

            // compute e3' : vector product of e1' x help
            e3.beVectorProductOf(e1, help);
            // let us normalize
            e3.normalize();

            // now from e3' x e1' compute e2'
            e2.beVectorProductOf(e3, e1);

            break;
        }
        }

        GtoLRotationMatrix = new FloatMatrix(3, 3);

        for (int i = 1; i <= 3; i++) {
                GtoLRotationMatrix->at(1, i) = e1.at(i);
                GtoLRotationMatrix->at(2, i) = e2.at(i);
                GtoLRotationMatrix->at(3, i) = e3.at(i);
        }
    }

    return GtoLRotationMatrix;
}

int
QdElement::computeLoadGToLRotationMtrx(FloatMatrix& answer)
// Returns the rotation matrix of the receiver of the size [6,6]
// f(local) = T * f(global)
{
    if (GtoLRotationMatrix == NULL) {
        this->computeGtoLRotationMatrix();
    }

    answer.resize(6, 6);
    answer.zero();

    for (int i = 1; i <= 3; i++) {
        answer.at(1, i) = answer.at(4, i + 3) = GtoLRotationMatrix->at(1, i);
        answer.at(2, i) = answer.at(5, i + 3) = GtoLRotationMatrix->at(2, i);
        answer.at(3, i) = answer.at(6, i + 3) = GtoLRotationMatrix->at(3, i);
    }

    return 1;
}

void
QdElement::computeLocalNodalCoordinates(std::vector< FloatArray >& lxy)
// Returns nodal coordinates in element's (local) coordinate system
{
    if (GtoLRotationMatrix == NULL) {
        this->computeGtoLRotationMatrix();
    }

    lxy.resize(4);
    for (int i = 0; i < 4; i++) {
        const auto& nc = this->giveNode(i + 1)->giveCoordinates();
        lxy[i].beProductOf(*GtoLRotationMatrix, nc);
    }
}

int
QdElement::computeNumberOfDofs()
{
    ///@todo move one hiearchy up and generalize
    IntArray dofIdMask;
    this->giveDofManDofIDMask(-1, dofIdMask); // ok for standard elements
    return this->giveInterpolation()->giveNumberOfNodes() * dofIdMask.giveSize();
}

void
QdElement::computeStrainVector(FloatArray& answer, GaussPoint* gp, TimeStep* tStep) {
    computeStrainVectorAt(answer, gp->giveNaturalCoordinate(1), gp->giveNaturalCoordinate(2), tStep);
}

void
QdElement::computeSurfaceNMatrix(FloatMatrix& answer, int boundaryID, const FloatArray& lcoords)
{
    FloatArray n_vec;
    this->giveInterpolation()->boundarySurfaceEvalN(n_vec, boundaryID, lcoords, *giveCellGeometryWrapper());
    answer.beNMatrixOf(n_vec, 6);
}

double
QdElement::computeSurfaceVolumeAround(GaussPoint* gp, int iSurf)
// Returns the portion of the receiver which is attached to gp.
{
    double detJ, weight;

    weight = gp->giveWeight();
    detJ = fabs(giveInterpolation()->giveTransformationJacobian(gp->giveNaturalCoordinates(), *giveCellGeometryWrapper()));
    return detJ * weight;
}

FEICellGeometry*
QdElement::giveCellGeometryWrapper()
{
    if (cellGeometryWrapper) {
        return cellGeometryWrapper;
    }
    else {
        this->computeLocalNodalCoordinates(localCoords);
        return (cellGeometryWrapper = new FEIVertexListGeometryWrapper(localCoords));
    }
}

std::vector< FloatArray >
QdElement::giveNodeCoordinates()
{
    IntArray newDofMans(this->dofManArray.giveSize());
    for (int i = 1; i <= this->dofManArray.giveSize(); i++)
        newDofMans.at(i) = this->dofManArray.at(i);

    setDofManagers(newDofMans);

    std::vector< FloatArray > c;
    computeLocalNodalCoordinates(c);
    return c;
}

void
QdElement::updateInternalState(TimeStep* tStep)
// Updates receiver at the end of the step.
{
    FloatArray stress, strain;
    // force updating strains & stresses
    switch (outputAtXY) {
    case OutputLocationXY::GaussPoints:
        for (auto& iRule : integrationRulesArray) {
            for (GaussPoint* gp : *iRule) {
                this->computeStrainVector(strain, gp, tStep);
                this->computeStressVector(stress, strain, gp, tStep);
            }
        }
        break;
    case OutputLocationXY::Centre:
        this->computeStrainVectorAt(strain, 0.0, 0.0, tStep);
        this->computeStressVectorAtCentre(stress, tStep, strain);
        break;
    case OutputLocationXY::Corners:
        OOFEM_ERROR("Not implemented for this element");
        break;
    case OutputLocationXY::All:
        OOFEM_ERROR("Not implemented for this element");
        break;
    deafult:
        OOFEM_ERROR("Something went wrong. Output requested at an unknown location in x-y plane of element %d.", giveGlobalNumber());
        break;
    }
}
}