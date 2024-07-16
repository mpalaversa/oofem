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

#include "sm/Elements/Bars/truss3d.h"
#include "sm/CrossSections/structuralcrosssection.h"
#include "sm/CrossSections/decoupledcrosssection.h"
#include "fei3dlinelin.h"
#include "node.h"
#include "material.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "mathfem.h"
#include "classfactory.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
REGISTER_Element(Truss3d);

FEI3dLineLin Truss3d :: interp;

Truss3d :: Truss3d(int n, Domain *aDomain) :
    NLStructuralElement(n, aDomain), ZZNodalRecoveryModelInterface(this)
{
    numberOfDofMans = 2;
}


FEInterpolation *Truss3d :: giveInterpolation() const { return & interp; }


Interface *
Truss3d :: giveInterface(InterfaceType interface)
{
    if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return static_cast< ZZNodalRecoveryModelInterface * >(this);
    } else if ( interface == NodalAveragingRecoveryModelInterfaceType ) {
        return static_cast< NodalAveragingRecoveryModelInterface * >(this);
    }

    //OOFEM_LOG_INFO("Interface on Truss3d element not supported");
    return NULL;
}


void
Truss3d :: NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node, InternalStateType type, TimeStep *tStep)
{
    answer.clear();
    OOFEM_WARNING("IP values will not be transferred to nodes. Use ZZNodalRecovery instead (parameter stype 1)");
}


void
Truss3d :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
//
// Returns linear part of geometrical equations of the receiver at gp.
// Returns the linear part of the B matrix
//
{
    FloatMatrix dN;
    this->interp.evaldNdx( dN, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(1, 6);
    answer.at(1, 1) = dN.at(1, 1);
    answer.at(1, 2) = dN.at(1, 2);
    answer.at(1, 3) = dN.at(1, 3);
    answer.at(1, 4) = dN.at(2, 1);
    answer.at(1, 5) = dN.at(2, 2);
    answer.at(1, 6) = dN.at(2, 3);
}


  void
Truss3d :: computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer)
{
    this->computeBmatrixAt(gp, answer);
}

void Truss3d ::computeHydrodynamicLoadVector( FloatArray &answer, FloatArray velocity, TimeStep *tStep )
{
    // Length of the element
    double l     = this->computeLength();
    // The associated decoupled cross-section
    DecoupledCrossSection *cs = this->giveDecoupledCrossSectionOfType( DecoupledMaterial::DecoupledMaterialType::DecoupledFluidMaterial );
    // Characteristic dimension of the cross-section
    double characteristicDim = cs->giveCharacteristicDimension();
    // Density of the fluid the element is immersed in
    double density = cs->giveMagnitudeOfMaterialProperty( 'd' );
    // Dynamic viscosity of the fluid material
    double mu      = cs->giveMaterial()->giveDynamicViscosity();

    FloatArray ex, u, currentNode1Coordinates, currentNode2Coordinates;
    this->computeVectorOf( VM_Total, tStep, u );
    // Form two vectors, one for current position of node 1 and the other for the node 2
    // by storing current displacements into the corresponding vectors.
    currentNode1Coordinates.resize( 3 );
    currentNode1Coordinates.at( 1 ) = u.at( 1 );
    currentNode1Coordinates.at( 2 )  = u.at( 2 );
    currentNode1Coordinates.at( 3 )  = u.at( 3 );
    currentNode2Coordinates.resize( 3 );
    currentNode2Coordinates.at( 1 )  = u.at( 4 );
    currentNode2Coordinates.at( 2 )  = u.at( 5 );
    currentNode2Coordinates.at( 3 )  = u.at( 6 );
    // Add coordinates of the initial position of the nodes
    currentNode1Coordinates.add( this->giveNode( 1 )->giveCoordinates() );
    currentNode2Coordinates.add( this->giveNode( 2 )->giveCoordinates());
    // Calculate unit vector along the element's longitudinal axis
    ex.beDifferenceOf( currentNode2Coordinates, currentNode1Coordinates );
    ex.normalize();

    // Get velocity of element nodes in the current time step
    FloatArray currentNodalVelocity;
    this->computeVectorOf( VM_Velocity, tStep, currentNodalVelocity );
    // Calculate velocity of the fluid relative to the element
    FloatArray relativeVelocity;
    relativeVelocity.resize( 3 );
    relativeVelocity.at( 1 ) = velocity.at( 1 ) - ( currentNodalVelocity.at( 1 ) + currentNodalVelocity.at( 4 ) ) / 2;
    relativeVelocity.at( 2 ) = velocity.at( 2 ) - ( currentNodalVelocity.at( 2 ) + currentNodalVelocity.at( 5 ) ) / 2;
    relativeVelocity.at( 3 ) = velocity.at( 3 ) - ( currentNodalVelocity.at( 3 ) + currentNodalVelocity.at( 6 ) ) / 2;
    // Calculate tangential component of the relative velocity
    double tangentialVelocityMagnitude = ex.dotProduct( relativeVelocity );
    // Form the tangential velocity component vector
    FloatArray tangentialVelocity, normalVelocity;
    tangentialVelocity.beScaled( tangentialVelocityMagnitude, ex );
    // Calculate normal velocity component, i.e. velocity component orthogonal
    // to the element, that contributes to the drag force (see Morison's equation)
    normalVelocity.beDifferenceOf(relativeVelocity, tangentialVelocity);
    double normalVelocityMagnitude = normalVelocity.computeNorm();

    FloatArray dragCoeffs;
    double userDefinedDragCoeff = cs->giveDragCoefficient();
    if ( userDefinedDragCoeff == 0.0 )
        dragCoeffs = computeDragCoefficients( density, mu, characteristicDim, normalVelocityMagnitude );
    else {
        dragCoeffs.resize( 2 );
        dragCoeffs.at( 1 ) = userDefinedDragCoeff;
        dragCoeffs.at( 2 ) = 0;
    }
    
    // Drag force on the element based on Morison's equation
    double normalDragForceMagnitude = 0.5 * density * dragCoeffs.at(1) *characteristicDim *l * pow(normalVelocityMagnitude, 2);
    double tangentialDragForceMagnitude = dragCoeffs.at(2) * l * pow( tangentialVelocityMagnitude, 2 );
    FloatArray dragForceVector, normalDragForceVector, tangentialDragForceVector;
    // Form the drag force vector (this vector is orthogonal to the element)
    normalVelocity.normalize();
    normalDragForceVector.beScaled( normalDragForceMagnitude / 2, normalVelocity );
    tangentialDragForceVector.beScaled( tangentialDragForceMagnitude / 2, ex );
    dragForceVector.add(normalDragForceVector);
    dragForceVector.add( tangentialDragForceVector );
    // The force is equally distributed among the element's nodes
    answer.resize( 6 );
    answer.at( 1 ) = answer.at( 4 ) = dragForceVector.at( 1 );
    answer.at( 2 ) = answer.at( 5 ) = dragForceVector.at( 2 );
    answer.at( 3 ) = answer.at( 6 ) = dragForceVector.at( 3 );

    // Get acceleration of element nodes in the current time step
    FloatArray currentNodalAcceleration;
    this->computeVectorOf( VM_Acceleration, tStep, currentNodalAcceleration );
    if ( currentNodalAcceleration.computeNorm() != 0 ) {
        // Calculate an average fluid acceleration on the element
        FloatArray averageAcceleration;
        averageAcceleration.resize( 3 );
        averageAcceleration.at( 1 ) = ( currentNodalAcceleration.at( 1 ) + currentNodalAcceleration.at( 4 ) ) / 2;
        averageAcceleration.at( 2 ) = ( currentNodalAcceleration.at( 2 ) + currentNodalAcceleration.at( 5 ) ) / 2;
        averageAcceleration.at( 3 ) = ( currentNodalAcceleration.at( 3 ) + currentNodalAcceleration.at( 6 ) ) / 2;
        // Calculate tangential component of the nodal acceleration
        double tangentialAccelerationMagnitude = ex.dotProduct( averageAcceleration );
        // Form the tangential acceleration component vector
        FloatArray tangentialAcceleration, normalAcceleration;
        tangentialAcceleration.beScaled( tangentialAccelerationMagnitude, ex );
        // Calculate normal acceleration component, i.e. acceleration component orthogonal
        // to the element, that contributes to the added-mass force (see Morison's equation)
        normalAcceleration.beDifferenceOf( averageAcceleration, tangentialAcceleration );

        // The added-mass force on the element based on Morison's equation
        double amForce = -density * ( characteristicDim * characteristicDim * 3.14 / 4 ) * l * normalAcceleration.computeNorm();
        FloatArray amForceVector;
        // The added-mass force vector (this vector is orthogonal to the element)
        normalAcceleration.normalize();
        amForceVector.beScaled( amForce / 2, normalAcceleration );
        // The force is equally distributed among the element's nodes
        answer.at( 1 ) = answer.at( 1 ) + amForceVector.at( 1 );
        answer.at( 4 ) = answer.at( 4 ) + amForceVector.at( 1 );
        answer.at( 2 ) = answer.at( 2 ) + amForceVector.at( 2 );
        answer.at( 5 ) = answer.at( 5 ) + amForceVector.at( 2 );
        answer.at( 3 ) = answer.at( 3 ) + amForceVector.at( 3 );
        answer.at( 6 ) = answer.at( 6 ) + amForceVector.at( 3 );
    }
}

FloatArray
Truss3d::computeDragCoefficients( double density, double mu, double characteristicDim, double relativeNormalVelocity ) {
    double reynoldsNo = density * characteristicDim * relativeNormalVelocity / mu;
    double s   = -0.077215665 + log( 8 / reynoldsNo );

    FloatArray dragCoeffs;
    dragCoeffs.resize( 2 );
    if ( reynoldsNo <= 1 )
        dragCoeffs.at(1) = 8 * 3.14 / ( reynoldsNo * s ) * ( 1 - 0.87 * pow(s,-2));
    else if ( reynoldsNo <= 30 )
        dragCoeffs.at( 1 ) = 1.45 + 8.55 * pow( reynoldsNo, -0.9 );
    else
        dragCoeffs.at( 1 ) = 1.1 + 4 * pow( reynoldsNo, -0.5 );

    dragCoeffs.at( 2 ) = 3.14 * mu * ( 0.55 * pow( reynoldsNo, 0.5 ) + 0.084 * pow( reynoldsNo, 2 / 3 ) );

    return dragCoeffs;
}

void
Truss3d :: computeGaussPoints()
// Sets up the array of Gauss Points of the receiver.
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize( 1 );
        integrationRulesArray [ 0 ] = std::make_unique<GaussIntegrationRule>(1, this, 1, 2);
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], 1, this);
    }
}


double
Truss3d :: computeLength()
{
    return this->interp.giveLength( FEIElementGeometryWrapper(this) );
}


void
Truss3d :: computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep)
// Returns the lumped mass matrix of the receiver. This expression is
// valid in both local and global axes.
{
    answer.resize(6, 6);
    answer.zero();
    if ( !this->isActivated(tStep) ) {
        return;
    }

    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    double density = this->giveStructuralCrossSection()->give('d', gp);
    double halfMass = density * this->giveCrossSection()->give(CS_Area, gp) * this->computeLength() * 0.5;
    answer.at(1, 1) = halfMass;
    answer.at(2, 2) = halfMass;
    answer.at(3, 3) = halfMass;
    answer.at(4, 4) = halfMass;
    answer.at(5, 5) = halfMass;
    answer.at(6, 6) = halfMass;
}


void
Truss3d :: computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer)
// Returns the displacement interpolation matrix {N} of the receiver, eva-
// luated at gp.
{
    FloatArray n;
    this->interp.evalN( n, iLocCoord, FEIElementGeometryWrapper(this) );
    answer.beNMatrixOf(n, 3);
}


double
Truss3d :: computeVolumeAround(GaussPoint *gp)
// Returns the length of the receiver. This method is valid only if 1
// Gauss point is used.
{
    double detJ = this->interp.giveTransformationJacobian( gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    double weight  = gp->giveWeight();
    return detJ *weight *this->giveCrossSection()->give(CS_Area, gp);
}


int
Truss3d :: giveLocalCoordinateSystem(FloatMatrix &answer)
//
// returns a unit vectors of local coordinate system at element
// stored rowwise (mainly used by some materials with ortho and anisotrophy)
//
{
    FloatArray lx, ly(3), lz;

    lx.beDifferenceOf( this->giveNode(2)->giveCoordinates(), this->giveNode(1)->giveCoordinates() );
    lx.normalize();

    ly(0) = lx(1);
    ly(1) = -lx(2);
    ly(2) = lx(0);

    // Construct orthogonal vector
    double npn = ly.dotProduct(lx);
    ly.add(-npn, lx);
    ly.normalize();
    lz.beVectorProductOf(ly, lx);

    answer.resize(3, 3);
    for ( int i = 1; i <= 3; i++ ) {
        answer.at(1, i) = lx.at(i);
        answer.at(2, i) = ly.at(i);
        answer.at(3, i) = lz.at(i);
    }

    return 1;
}

void
Truss3d :: initializeFrom(InputRecord &ir)
{
    NLStructuralElement :: initializeFrom(ir);
}


void
Truss3d :: computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    answer = this->giveStructuralCrossSection()->giveRealStress_1d(strain, gp, tStep);
}

void
Truss3d :: computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    answer = this->giveStructuralCrossSection()->giveStiffnessMatrix_1d(rMode, gp, tStep);
}

void
Truss3d :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {D_u, D_v, D_w};
}


void
Truss3d :: giveEdgeDofMapping(IntArray &answer, int iEdge) const
{
    /*
     * provides dof mapping of local edge dofs (only nonzero are taken into account)
     * to global element dofs
     */

    if ( iEdge != 1 ) {
        OOFEM_ERROR("wrong edge number");
    }

    answer.resize(6);
    answer.at(1) = 1;
    answer.at(2) = 2;
    answer.at(3) = 3;
    answer.at(4) = 4;
    answer.at(5) = 5;
    answer.at(6) = 6;
}


double
Truss3d :: computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
    if ( iEdge != 1 ) { // edge between nodes 1 2
        OOFEM_ERROR("wrong edge number");
    }

    double weight = gp->giveWeight();
    return this->interp.giveTransformationJacobian( gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) ) * weight;
}


int
Truss3d :: computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp)
{
    // returns transformation matrix from
    // edge local coordinate system
    // to element local c.s
    // (same as global c.s in this case)
    //
    // i.e. f(element local) = T * f(edge local)
    //
    FloatMatrix lcs;
    this->giveLocalCoordinateSystem(lcs);
    answer.beTranspositionOf(lcs);

    return 1;
}


#ifdef __OOFEG
void Truss3d :: drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep)
{
    GraphicObj *go;
    //  if (!go) { // create new one
    WCRec p [ 2 ]; /* point */
    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getElementColor() );
    EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveCoordinate(1);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveCoordinate(2);
    p [ 0 ].z = ( FPNum ) this->giveNode(1)->giveCoordinate(3);
    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveCoordinate(1);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveCoordinate(2);
    p [ 1 ].z = ( FPNum ) this->giveNode(2)->giveCoordinate(3);
    go = CreateLine3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
}


void Truss3d :: drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType type)
{
    GraphicObj *go;
    double defScale = gc.getDefScale();
    //  if (!go) { // create new one
    WCRec p [ 2 ]; /* point */
    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getDeformedElementColor() );
    EASValsSetLayer(OOFEG_DEFORMED_GEOMETRY_LAYER);
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 0 ].z = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(3, tStep, defScale);

    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 1 ].z = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(3, tStep, defScale);
    go = CreateLine3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
    EMAddGraphicsToModel(ESIModel(), go);
}
#endif
} // end namespace oofem
