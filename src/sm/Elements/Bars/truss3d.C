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

void Truss3d ::computeHydrodynamicLoadVector( FloatArray &answer, FloatArray flowCharacteristics, TimeStep *tStep )
{
    FloatArray et, u, currentNode1Coordinates, currentNode2Coordinates;
    this->computeVectorOf( VM_Total, tStep, u );
    // Form two vectors, one for current position of node 1 and the other for the node 2
    // by storing current displacements into the corresponding vectors.
    currentNode1Coordinates.resize( 3 );
    currentNode1Coordinates.at( 1 ) = u.at( 1 );
    currentNode1Coordinates.at( 2 ) = u.at( 2 );
    currentNode1Coordinates.at( 3 ) = u.at( 3 );
    currentNode2Coordinates.resize( 3 );
    currentNode2Coordinates.at( 1 ) = u.at( 4 );
    currentNode2Coordinates.at( 2 ) = u.at( 5 );
    currentNode2Coordinates.at( 3 ) = u.at( 6 );
    // Add coordinates of the initial position of the nodes
    currentNode1Coordinates.add( this->giveNode( 1 )->giveCoordinates() );
    currentNode2Coordinates.add( this->giveNode( 2 )->giveCoordinates() );
    // Calculate unit vector along the element's longitudinal axis (tangential direction)
    et.beDifferenceOf( currentNode2Coordinates, currentNode1Coordinates );
    et.normalize();

    // Form the fluid velocity vector
    FloatArray velocity;
    velocity.resize( 3 );
    for ( int i = 1; i <= 3; i++ )
        velocity.at( i ) = flowCharacteristics.at( i );
    
    // Check if the element is downstream relative to another element
    if ( this->isDownstream )
        // Reduce the inflow velocity by the given velocity reduction factor (r)
        velocity.beScaled( this->velocityReductionFactor, velocity );

    // Get velocity and acceleration of element nodes in the current time step
    FloatArray currentNodalVelocity, currentNodalAcceleration;
    this->computeVectorOf( VM_Velocity, tStep, currentNodalVelocity );
    this->computeVectorOf( VM_Acceleration, tStep, currentNodalAcceleration );

    // Calculate velocity of the fluid relative to the element
    FloatArray relativeVelocity;
    relativeVelocity.resize( 3 );
    relativeVelocity.at( 1 ) = velocity.at( 1 ) - ( currentNodalVelocity.at( 1 ) + currentNodalVelocity.at( 4 ) ) / 2;
    relativeVelocity.at( 2 ) = velocity.at( 2 ) - ( currentNodalVelocity.at( 2 ) + currentNodalVelocity.at( 5 ) ) / 2;
    relativeVelocity.at( 3 ) = velocity.at( 3 ) - ( currentNodalVelocity.at( 3 ) + currentNodalVelocity.at( 6 ) ) / 2;

    // Calculate tangential component of the relative velocity
    FloatArray tangentialRelativeVelocity;
    tangentialRelativeVelocity.beScaled( relativeVelocity.dotProduct( et ), et );

    // Calculate normal component of the relative velocity
    FloatArray normalRelativeVelocity;
    normalRelativeVelocity.beDifferenceOf( relativeVelocity, tangentialRelativeVelocity );

    // Fetch length of the element, characteristic dimension of the cross-section and density and dynamic viscosity of the fluid
    double l = this->computeLength();
    DecoupledCrossSection *cs = this->giveDecoupledCrossSectionOfType( DecoupledMaterial::DecoupledMaterialType::DecoupledFluidMaterial );
    double characteristicDim = cs->giveCharacteristicDimension();
    double density = cs->giveMagnitudeOfMaterialProperty( 'd' );
    double mu = cs->giveMaterial()->giveDynamicViscosity();
    
    // Caluclate drag coefficients
    FloatArray dragCoeffs;
    // Check if the user has defined a constant value of the drag coefficient
    double userDefinedDragCoeff = cs->giveDragCoefficient();
    if ( userDefinedDragCoeff == 0.0 )
        // If no, calculate the drag coefficients.
        dragCoeffs = computeDragCoefficients( density, mu, characteristicDim, normalRelativeVelocity.computeNorm() );
    else {
        // If yes, use the user defined drag coefficient for the normal viscous force component
        // and disregard the tangential one.
        dragCoeffs.resize( 2 );
        dragCoeffs.at( 1 ) = userDefinedDragCoeff;
        dragCoeffs.at( 2 ) = 0;
    }

    // Caluclate viscous force
    FloatArray normalViscousForce, tangentialViscousForce;
    normalViscousForce.beScaled( 0.5 * density * dragCoeffs.at( 1 ) * characteristicDim * l * normalRelativeVelocity.computeNorm(), normalRelativeVelocity );
    tangentialViscousForce.beScaled( dragCoeffs.at( 2 ) * l, tangentialRelativeVelocity );
    this->viscousForce.zero();
    this->viscousForce.add( normalViscousForce );
    this->viscousForce.add( tangentialViscousForce );

    // The force is equally distributed among the element's nodes
    answer.resize( 6 );
    answer.at( 1 ) = answer.at( 4 ) = viscousForce.at( 1 ) / 2;
    answer.at( 2 ) = answer.at( 5 ) = viscousForce.at( 2 ) / 2;
    answer.at( 3 ) = answer.at( 6 ) = viscousForce.at( 3 ) / 2;

    if ( currentNodalAcceleration.computeNorm() != 0 ) {
        // Form the fluid acceleration vector
        FloatArray acceleration;
        acceleration.resize( 3 );
        for ( int i = 4; i <= 6; i++ )
            acceleration.at( i - 3 ) = flowCharacteristics.at( i );

        // Calculate an average fluid acceleration on the element - this should be changed when the fluid acceleration becomes available as an input quantity
        FloatArray relativeAcceleration;
        relativeAcceleration.resize( 3 );
        relativeAcceleration.at( 1 ) = acceleration.at( 1 ) - ( currentNodalAcceleration.at( 1 ) + currentNodalAcceleration.at( 4 ) ) / 2;
        relativeAcceleration.at( 2 ) = acceleration.at( 2 ) - ( currentNodalAcceleration.at( 2 ) + currentNodalAcceleration.at( 5 ) ) / 2;
        relativeAcceleration.at( 3 ) = acceleration.at( 3 ) - ( currentNodalAcceleration.at( 3 ) + currentNodalAcceleration.at( 6 ) ) / 2;

        // Calculate tangential component of the relative acceleration
        FloatArray tangentialRelativeAcceleration;
        tangentialRelativeAcceleration.beScaled( relativeAcceleration.dotProduct( et ), et );

        // Calculate normal component of the relative acceleration
        FloatArray normalRelativeAcceleration;
        normalRelativeAcceleration.beDifferenceOf( relativeAcceleration, tangentialRelativeAcceleration );

        // Calculate tangential component of the fluid acceleration
        FloatArray tangentialAcceleration;
        tangentialAcceleration.beScaled( acceleration.dotProduct( et ), et );

        // Calculate normal component of the fluid acceleration
        FloatArray normalAcceleration;
        normalAcceleration.beDifferenceOf( acceleration, tangentialAcceleration );
        
        // Added-mass coefficient
        double cm = cs->giveAddedMassCoefficient();

        // Calculate the added-mass force. [CURRENTLY ASSUMES A CIRCULAR CROSS-SECTION.]
        FloatArray addedMassForce;
        addedMassForce.beScaled( density * ( pow( characteristicDim, 2 ) * 3.14 / 4 ) * l, normalAcceleration );
        addedMassForce.add( density * ( pow( characteristicDim, 2 ) * 3.14 / 4 ) * l * cm, normalRelativeAcceleration );

        // The force is equally distributed among the element's nodes
        answer.at( 1 ) = answer.at( 1 ) + addedMassForce.at( 1 ) / 2;
        answer.at( 4 ) = answer.at( 4 ) + addedMassForce.at( 1 ) / 2;
        answer.at( 2 ) = answer.at( 2 ) + addedMassForce.at( 2 ) / 2;
        answer.at( 5 ) = answer.at( 5 ) + addedMassForce.at( 2 ) / 2;
        answer.at( 3 ) = answer.at( 3 ) + addedMassForce.at( 3 ) / 2;
        answer.at( 6 ) = answer.at( 6 ) + addedMassForce.at( 3 ) / 2;
    }
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
