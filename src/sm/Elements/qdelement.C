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

        //std::vector< FloatArray > locCoords = giveNodeCoordinates();

        switch (csClass)
        {
        case CSClass::Nastran:
        {
            /*
                In order to obtain three unit vectors (e1, e2 and e3) according to definition of the "Nastran's coordinate system", we start by first defining
                an element coordinate system in the "standard OOFEM" manner. This coord. sys. will be referred to as the "native OOFEM" coord. sys. and all variables
                related with such a system will be denoted by "O".
            */
            // Vectors spanned bewteen nodes 1 and 2 and 1 and 3 respectively.
            FloatArray v12, v13;
            // v12 is a vector in the x-direction of the element coord. sys. when it is defined as the native OOFEM coord. sys. (e.g. we could call it vXO).
            v12.beDifferenceOf(node2, node1);
            v13.beDifferenceOf(node3, node1);

            // Vectors in direction of y and z axes of the element coord. sys. fixed to the plane of the element and defined as the native OOFEM coord. sys.
            FloatArray vZO, vYO;
            vZO.beVectorProductOf(v12, v13);
            vYO.beVectorProductOf(vZO, v12);

            /*
                Another coordinate system that we will be making use of here is the global coordinate system. Variables related to vectors defined in this system
                will not have any special endings in this code (e.g. x, y, z, etc.).
            */

            // Unit vectors in direction of the axes of the global coord. syst.
            FloatArray x, y, z;
            x.resize(3);
            x.at(1) = 1.0;
            y.resize(3);
            y.at(2) = 1.0;
            z.resize(3);
            z.at(3) = 1.0;

            // Transformation matrix from the global to the coord. sys. fixed to the element's plane.
            FloatMatrix transMatGlobToEl;
            transMatGlobToEl.resize(3, 3);
            transMatGlobToEl.at(1, 1) = v12.dotProduct(x) / v12.computeNorm();
            transMatGlobToEl.at(1, 2) = v12.dotProduct(y) / v12.computeNorm();
            transMatGlobToEl.at(1, 3) = v12.dotProduct(z) / v12.computeNorm();
            transMatGlobToEl.at(2, 1) = vYO.dotProduct(x) / vYO.computeNorm();
            transMatGlobToEl.at(2, 2) = vYO.dotProduct(y) / vYO.computeNorm();
            transMatGlobToEl.at(2, 3) = vYO.dotProduct(z) / vYO.computeNorm();
            transMatGlobToEl.at(3, 1) = vZO.dotProduct(x) / vZO.computeNorm();
            transMatGlobToEl.at(3, 2) = vZO.dotProduct(y) / vZO.computeNorm();
            transMatGlobToEl.at(3, 3) = vZO.dotProduct(z) / vZO.computeNorm();

            /*
                The third set of variable names that we will be using here is related to vectors given in the element's coordinate system when it is defined
                according to the native OOFEM coord. sys. These will end with "El".
            */
            
            // Vectors spanned bewteen nodes 1 and 2, 1 and 3 and 2 and 4 respectively in the coord. sys. fixed to the element's plane.
            FloatArray v12El, v13El, v24, v24El;
            v12El.beProductOf(transMatGlobToEl, v12);
            v13El.beProductOf(transMatGlobToEl, v13);
            v24.beDifferenceOf(node4, node2);
            v24El.beProductOf(transMatGlobToEl, v24);

            // Find angles between vectors v12El and v13El (\beta) and between vectors v12EL and v24El (\gamma).
            double cosBeta{ 0.0 }, cosGamma{ 0.0 };
            FloatArray v12ElTemp;
            v12ElTemp.beScaled(-1.0, v12El);
            cosBeta = (v12El.dotProduct(v13El)) / (v12El.computeNorm() * v13El.computeNorm());
            cosGamma = (v12ElTemp.dotProduct(v24El)) / (v12El.computeNorm() * v24El.computeNorm());

            // Calculate angle between the x-axis of the elment coord. sys. and vC3El.
            double alpha = (acos(cosBeta) + acos(cosGamma)) / 2;

            // Calculate vector that is coincident with the x-axis of the element coord. sys.
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

            // Vector vXEl is coincident with the x-axis of the new ("Nastran's") coord. sys.
            FloatArray vXEl;
            vXEl.beProductOf(rotMat, v13El);

            // Calculate vector coincident with the z-axis of the element ("Nastran's") coord. sys.
            FloatArray vZEl;
            vZEl.beVectorProductOf(vXEl, v13El);

            // Calculate vector coincident with th y-axis of the element ("Nastran's") coord. sys.
            FloatArray vYEl;
            vYEl.beVectorProductOf(vZEl, vXEl);

            // Vectors in direction of the element ("Nastran's") coord. sys. axes given in the global coordinate system.
            FloatArray vX, vY, vZ;
            vX.beTProductOf(transMatGlobToEl, vXEl);
            vY.beTProductOf(transMatGlobToEl, vYEl);
            vZ.beTProductOf(transMatGlobToEl, vZEl);

            // Calculate unit vectors in the direction of the x, y and z axis of the element.
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