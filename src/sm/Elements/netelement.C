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
#include "fei3dlinelin.h"
#include "sm/Elements/netelement.h"
#include "../sm/CrossSections/structuralcrosssection.h"
#include "gaussintegrationrule.h"

#include "classfactory.h"
#include "gausspoint.h"

namespace oofem {
FEI3dLineLin NetElement ::interp;

NetElement::NetElement(int n, Domain* aDomain) : NLStructuralElement(n, aDomain)
{
        GtoLRotationMatrix = NULL;
        cellGeometryWrapper = NULL;
}

void
NetElement::computeBmatrixAt(GaussPoint* gp, FloatMatrix& answer, int lowerIndx, int upperIndx) {
    computeBmatrixAt(gp->giveNaturalCoordinate(1), gp->giveNaturalCoordinate(2), answer);
}

void
NetElement ::computeConstitutiveMatrixAt( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep )
{
    answer = this->giveStructuralCrossSection()->giveStiffnessMatrix_1d( rMode, gp, tStep );
}

void
NetElement ::computeGaussPoints()
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize( 1 );
        integrationRulesArray[0] = std::make_unique<GaussIntegrationRule>( 1, this, 1, 2 );
        this->giveCrossSection()->setupIntegrationPoints( *integrationRulesArray[0], 1, this );
    }
}

const FloatMatrix*
NetElement::computeGtoLRotationMatrix()
// Returns the rotation matrix of the receiver of the size [3,3]
// coords(local) = T * coords(global)
//
// OOFEM's local coordinate system (described by vector triplet e1',e2',e3') is defined as follows:
// e1'    : [N2-N1]    Ni - means i - th node
// v31   : [N3-N1]
// e3'    : e1' x v31
// e2'    : e3' x e1'
{
    if (GtoLRotationMatrix == NULL) {
        FloatArray node1 = this->giveNode(1)->giveCoordinates();
        FloatArray node2 = this->giveNode(2)->giveCoordinates();
        FloatArray node3 = this->giveNode(3)->giveCoordinates();
        FloatArray node4 = this->giveNode(4)->giveCoordinates();
        
        FloatArray e1, e2, e3;
        e1.beDifferenceOf( node2, node1 );
        e1.normalize();
        FloatArray v31;
        v31.beDifferenceOf( node3, node1 );
        e3.beVectorProductOf( e1, v31 );
        e3.normalize();
        e2.beVectorProductOf( e3, e1 );

        GtoLRotationMatrix = new FloatMatrix(3, 3);
        for (int i = 1; i <= 3; i++) {
                GtoLRotationMatrix->at(1, i) = e1.at(i);
                GtoLRotationMatrix->at(2, i) = e2.at(i);
                GtoLRotationMatrix->at(3, i) = e3.at(i);
        }
    }

    return GtoLRotationMatrix;
}

double
NetElement::computeLength()
{
    return this->interp.giveLength( FEIElementGeometryWrapper( this ) );
}

int
NetElement::computeLoadGToLRotationMtrx(FloatMatrix& answer)
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
NetElement::computeLocalNodalCoordinates(std::vector< FloatArray >& lxy)
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

void
NetElement::computeSurfaceNMatrix(FloatMatrix& answer, int boundaryID, const FloatArray& lcoords)
{
    FloatArray n_vec;
    this->giveInterpolation()->boundarySurfaceEvalN(n_vec, boundaryID, lcoords, *giveCellGeometryWrapper());
    answer.beNMatrixOf(n_vec, 6);
}

double
NetElement::computeSurfaceVolumeAround(GaussPoint* gp, int iSurf)
// Returns the portion of the receiver which is attached to gp.
{
    double detJ, weight;

    weight = gp->giveWeight();
    detJ = fabs(giveInterpolation()->giveTransformationJacobian(gp->giveNaturalCoordinates(), *giveCellGeometryWrapper()));
    return detJ * weight;
}

FEICellGeometry*
NetElement::giveCellGeometryWrapper()
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
NetElement::giveNodeCoordinates()
{
    IntArray newDofMans(this->dofManArray.giveSize());
    for (int i = 1; i <= this->dofManArray.giveSize(); i++)
        newDofMans.at(i) = this->dofManArray.at(i);

    setDofManagers(newDofMans);

    std::vector< FloatArray > c;
    computeLocalNodalCoordinates(c);
    return c;
}

FEInterpolation*
NetElement::giveInterpolation() const { return &interp; }

void
NetElement::updateInternalState(TimeStep* tStep)
// Updates receiver at the end of the step.
{
    FloatArray stress, strain;
    // force updating strains & stresses
}
}