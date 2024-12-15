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

#include "sm/Elements/NetElements/netrec4trla.h"
#include "../sm/CrossSections/decoupledcrosssection.h"
#include "../sm/CrossSections/structuralcrosssection.h"
#include "sm/Materials/structuralms.h"

#include "classfactory.h"
#include "gaussintegrationrule.h"
#include "gausspoint.h"
#include "fei2dquadlin.h"

namespace oofem {
REGISTER_Element( NetRec4TrLa );

FEI2dQuadLin NetRec4TrLa::interpolation( 1, 2 );

NetRec4TrLa::NetRec4TrLa(int n, Domain* aDomain) : NetElement(n, aDomain)
{
    numberOfDofMans = 4;
    numberOfGaussPoints = 1;
    
    mask = -1;
    a0   = -1;
    dx0  = -1;
    dy0  = -1;
    nx   = -1;
    ny   = -1;
    initialDimensions.resize( 2 );
}

FloatArray
NetRec4TrLa::calculateCurrentUnitNormalToElement(TimeStep* tStep) {
    FloatArray node1 = this->giveNode( 1 )->giveCoordinates();
    FloatArray node2 = this->giveNode( 2 )->giveCoordinates();
    FloatArray node3 = this->giveNode( 3 )->giveCoordinates();

    if ( !tStep->isTheFirstStep() ) {
        // Fetch total displacements in the current configuration
        FloatArray u;
        u.resize( 12 );
        this->computeVectorOf( VM_Total, tStep, u );
        // Current coordinates of node 1
        double x      = node1.at( 1 ) + u.at( 1 );
        double y      = node1.at( 2 ) + u.at( 2 );
        double z      = node1.at( 3 ) + u.at( 3 );
        node1.at( 1 ) = x;
        node1.at( 2 ) = y;
        node1.at( 3 ) = z;
        // Current coordinates of node 2
        x             = node2.at( 1 ) + u.at( 4 );
        y             = node2.at( 2 ) + u.at( 5 );
        z             = node2.at( 3 ) + u.at( 6 );
        node2.at( 1 ) = x;
        node2.at( 2 ) = y;
        node2.at( 3 ) = z;
        // Current coordinates of node 3
        x             = node3.at( 1 ) + u.at( 7 );
        y             = node3.at( 2 ) + u.at( 8 );
        z             = node3.at( 3 ) + u.at( 9 );
        node3.at( 1 ) = x;
        node3.at( 2 ) = y;
        node3.at( 3 ) = z;
    }
    
    FloatArray e1, e2, e3;
    e1.beDifferenceOf( node2, node1 );
    e1.normalize();
    FloatArray v31;
    v31.beDifferenceOf( node3, node1 );
    e3.beVectorProductOf( e1, v31 );
    e3.normalize();

    return e3;
}

FloatArray
NetRec4TrLa::calculateElementDimensions( TimeStep *tStep )
{
    // Fetch initial position of the nodes
    FloatArray node1 = this->giveNode( 1 )->giveCoordinates();
    FloatArray node2 = this->giveNode( 2 )->giveCoordinates();
    FloatArray node3 = this->giveNode( 3 )->giveCoordinates();
    FloatArray node4 = this->giveNode( 4 )->giveCoordinates();

    if ( !tStep->isTheFirstStep() ) {
        // Fetch total displacements in the current configuration
        FloatArray u;
        u.resize( 12 );
        this->computeVectorOf( VM_Total, tStep, u );

        // Update nodes' positions
        double x = node1.at( 1 ) + u.at( 1 );
        double y = node1.at( 2 ) + u.at( 2 );
        double z      = node1.at( 3 ) + u.at( 3 );
        node1.at( 1 ) = x;
        node1.at( 2 ) = y;
        node1.at( 3 ) = z;

        x = node2.at( 1 ) + u.at( 4 );
        y = node2.at( 2 ) + u.at( 5 );
        z             = node2.at( 3 ) + u.at( 6 );
        node2.at( 1 ) = x;
        node2.at( 2 ) = y;
        node2.at( 3 ) = z;

        x             = node3.at( 1 ) + u.at( 7 );
        y             = node3.at( 2 ) + u.at( 8 );
        z             = node3.at( 3 ) + u.at( 9 );
        node3.at( 1 ) = x;
        node3.at( 2 ) = y;
        node3.at( 3 ) = z;
        
        x = node4.at( 1 ) + u.at( 10 );
        y = node4.at( 2 ) + u.at( 11 );
        z = node4.at( 3 ) + u.at( 12 );
        node4.at( 1 ) = x;
        node4.at( 2 ) = y;
        node4.at( 3 ) = z;
    }

    // Calculate length of the element's sides in the current configuration
    FloatArray side12, side23, side34, side14;
    side12.beDifferenceOf( node2, node1 );
    side14.beDifferenceOf( node4, node1 );
    side23.beDifferenceOf( node3, node2 );
    side34.beDifferenceOf( node4, node3 );
    FloatArray dimensions;
    dimensions.resize( 4 );
    dimensions.at( 1 ) = side12.computeNorm();
    dimensions.at( 2 ) = side14.computeNorm();
    dimensions.at( 3 ) = side34.computeNorm();
    dimensions.at( 4 ) = side23.computeNorm();
        
    return dimensions;
}

void
NetRec4TrLa::calculateEquivalentLumpedNodalValues( FloatArray &answer, FloatArray vector )
{
    answer.resize( 12 );
    answer.at( 1 ) = answer.at( 4 ) = answer.at( 7 ) = answer.at( 10 ) = vector.at( 1 ) / 4;
    answer.at( 2 ) = answer.at( 5 ) = answer.at( 8 ) = answer.at( 11 ) = vector.at( 2 ) / 4;
    answer.at( 3 ) = answer.at( 6 ) = answer.at( 9 ) = answer.at( 12 ) = vector.at( 3 ) / 4;
}

FloatArray
NetRec4TrLa::calculateInternalDisplacements( TimeStep *tStep )
{
    FloatArray currentDimensions;
    currentDimensions = calculateElementDimensions( tStep );
    
    FloatArray uInternal;
    uInternal.resize( 2 );
    uInternal.at( 1 ) = currentDimensions.at( 1 ) - initialDimensions.at( 1 );
    uInternal.at( 2 ) = currentDimensions.at( 2 ) - initialDimensions.at( 2 );

    return uInternal;
}

FloatArray
NetRec4TrLa::calculateRelativeAcceleration( FloatArray acceleration, TimeStep *tStep )
{
    // Get acceleration of element nodes in the current time step
    FloatArray currentNodalAcceleration;
    this->computeVectorOf( VM_Acceleration, tStep, currentNodalAcceleration );

    // Calculate an average fluid acceleration on the element - this should be changed when the fluid acceleration becomes available as an input quantity
    FloatArray relativeAcceleration;
    relativeAcceleration.resize( 3 );
    relativeAcceleration.at( 1 ) = acceleration.at( 1 ) - ( currentNodalAcceleration.at( 1 ) + currentNodalAcceleration.at( 4 ) + currentNodalAcceleration.at( 7 ) + currentNodalAcceleration.at( 10 ) ) / 4;
    relativeAcceleration.at( 2 ) = acceleration.at( 2 ) - ( currentNodalAcceleration.at( 2 ) + currentNodalAcceleration.at( 5 ) + currentNodalAcceleration.at( 8 ) + currentNodalAcceleration.at( 11 ) ) / 4;
    relativeAcceleration.at( 3 ) = acceleration.at( 3 ) - ( currentNodalAcceleration.at( 3 ) + currentNodalAcceleration.at( 6 ) + currentNodalAcceleration.at( 9 ) + currentNodalAcceleration.at( 12 ) ) / 4;

    return relativeAcceleration;
}

FloatArray
NetRec4TrLa::calculateRelativeVelocity( FloatArray velocity, TimeStep *tStep )
{
    // Get velocity and acceleration of element nodes in the current time step
    FloatArray currentNodalVelocity;
    this->computeVectorOf( VM_Velocity, tStep, currentNodalVelocity );

    // Calculate velocity of the fluid relative to the element
    FloatArray relativeVelocity;
    relativeVelocity.resize( 3 );
    relativeVelocity.at( 1 ) = velocity.at( 1 ) - ( currentNodalVelocity.at( 1 ) + currentNodalVelocity.at( 4 ) + currentNodalVelocity.at( 7 ) + currentNodalVelocity.at( 10 ) ) / 4;
    relativeVelocity.at( 2 ) = velocity.at( 2 ) - ( currentNodalVelocity.at( 2 ) + currentNodalVelocity.at( 5 ) + currentNodalVelocity.at( 8 ) + currentNodalVelocity.at( 11 ) ) / 4;
    relativeVelocity.at( 3 ) = velocity.at( 3 ) - ( currentNodalVelocity.at( 3 ) + currentNodalVelocity.at( 6 ) + currentNodalVelocity.at( 9 ) + currentNodalVelocity.at( 12 ) ) / 4;

    return relativeVelocity;
}

void
NetRec4TrLa::computeGaussPoints()
{
    // Sets up the integration rule array which contains all the Gauss points
    // Default: create one integration rule
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize( 1 );
        integrationRulesArray[0] = std::make_unique<GaussIntegrationRule>( 1, this, 1, 3 );
        this->giveCrossSection()->setupIntegrationPoints( *integrationRulesArray[0], this->numberOfGaussPoints, this );
    }
}

FloatMatrix
NetRec4TrLa::calculateTransformationMatrix( TimeStep *tStep )
// Returns the transformation matrix of the receiver of the size [8,12]
// r(local) = T * r(global)
// for one node (r written transposed): {u,v} = T * {u,v,w}
{
    FloatArray node1 = this->giveNode( 1 )->giveCoordinates();
    FloatArray node2 = this->giveNode( 2 )->giveCoordinates();
    FloatArray node3 = this->giveNode( 3 )->giveCoordinates();

    if ( tStep != 0 ) {
        // Fetch total displacements in the current configuration
        FloatArray u;
        u.resize( 12 );
        if ( !tStep->isTheFirstStep() )
            this->computeVectorOf( VM_Total, tStep, u );
        // Current coordinates of node 1
        double x      = node1.at( 1 ) + u.at( 1 );
        double y      = node1.at( 2 ) + u.at( 2 );
        double z      = node1.at( 3 ) + u.at( 3 );
        node1.at( 1 ) = x;
        node1.at( 2 ) = y;
        node1.at( 3 ) = z;
        // Current coordinates of node 2
        x             = node2.at( 1 ) + u.at( 4 );
        y             = node2.at( 2 ) + u.at( 5 );
        z             = node2.at( 3 ) + u.at( 6 );
        node2.at( 1 ) = x;
        node2.at( 2 ) = y;
        node2.at( 3 ) = z;
        // Current coordinates of node 3
        x             = node3.at( 1 ) + u.at( 7 );
        y             = node3.at( 2 ) + u.at( 8 );
        z             = node3.at( 3 ) + u.at( 9 );
        node3.at( 1 ) = x;
        node3.at( 2 ) = y;
        node3.at( 3 ) = z;
    }

    FloatArray e1, e2, e3;
    e1.beDifferenceOf( node2, node1 );
    e1.normalize();
    FloatArray v31;
    v31.beDifferenceOf( node3, node1 );
    e3.beVectorProductOf( e1, v31 );
    e3.normalize();
    e2.beVectorProductOf( e3, e1 );

    FloatMatrix GtoLRotationMatrix;
    GtoLRotationMatrix.resize( 3, 3 );
    for ( int i = 1; i <= 3; i++ ) {
        GtoLRotationMatrix.at( 1, i ) = e1.at( i );
        GtoLRotationMatrix.at( 2, i ) = e2.at( i );
        GtoLRotationMatrix.at( 3, i ) = e3.at( i );
    }

    FloatMatrix answer;
    answer.resize( 8, 12 );
    answer.zero();

    for ( int i = 1; i <= 3; i++ ) {
        answer.at( 1, i ) = answer.at( 3, i + 3 ) = answer.at( 5, i + 6 ) = answer.at( 7, i + 9 ) = GtoLRotationMatrix.at( 1, i );
        answer.at( 2, i ) = answer.at( 4, i + 3 ) = answer.at( 6, i + 6 ) = answer.at( 8, i + 9 ) = GtoLRotationMatrix.at( 2, i );
    }

    return answer;
}

void
NetRec4TrLa ::computeLumpedMassMatrix( FloatMatrix &answer, TimeStep *tStep )
// Returns the lumped mass matrix of the receiver. This expression is
// valid in both local and global axes.
{
    answer.resize( 12, 12 );
    if ( !this->isActivated( tStep ) ) {
        return;
    }

    if ( initialDimensions.at( 1 ) == 0 && tStep->isTheFirstStep() ) {
        if ( nx == -1 )
            initialDimensions = calculateElementDimensions( tStep );
        else {
            initialDimensions.at( 1 ) = nx * a0;
            initialDimensions.at( 2 ) = ny * a0;
        }
    }
    // Calculate number of twines within the element
    double noOfTwines = giveNumberOfTwines();

    GaussPoint *gp    = integrationRulesArray[0]->getIntegrationPoint( 0 );
    double density    = this->giveStructuralCrossSection()->give( 'd', gp );
    double quarterMass   = 0.25 * noOfTwines * density * this->giveCrossSection()->give( CS_Area, gp ) * a0;
    
    answer.at( 1, 1 ) = answer.at( 2, 2 ) = answer.at( 3, 3 ) = answer.at( 4, 4 ) = answer.at( 5, 5 ) = answer.at( 6, 6 ) = answer.at( 7, 7 ) = answer.at( 8, 8 ) = 
        answer.at( 9, 9 ) = answer.at( 10, 10 ) = answer.at( 11, 11 ) = answer.at( 12, 12 ) = quarterMass;
}

void
NetRec4TrLa::computeStiffnessMatrix(FloatMatrix& answer, MatResponseMode rMode, TimeStep* tStep) {
    if ( initialDimensions.at( 1 ) == 0 && tStep->isTheFirstStep() ) {
        if (nx == -1)
            initialDimensions = calculateElementDimensions( tStep );
        else {
            initialDimensions.at( 1 ) = nx * a0;
            initialDimensions.at( 2 ) = ny * a0;
        }
    }
    if ( initialDimensions.at( 1 ) == 0)
        OOFEM_ERROR( "Element dimensions are not calculated in the first step (NetRec4TrLa)." );

    FloatMatrix D;
    // Area and modulus of elasticity of a twine
    double At, Et;
    // Fetch At and Et from the user's input data
    for ( auto gp : *this->giveDefaultIntegrationRulePtr() ) {
        this->computeConstitutiveMatrixAt( D, rMode, gp, tStep );
        At = this->giveCrossSection()->give( CS_Area, gp );
    }
    Et = D.at( 1, 1 );

    // Get internal displacements
    FloatArray internalDisplacements;
    internalDisplacements.resize( 2 );
    internalDisplacements = calculateInternalDisplacements( tStep );
    
    FloatMatrix ke;
    if ( mask == 1 ) {
        // Calculate element's stiffness matrix
        ke.resize( 8, 8 );
        ke.at( 1, 1 ) = ke.at( 3, 3 ) = ke.at( 5, 5 ) = ke.at( 7, 7 ) = ( ny / 2 ) * ( At * Et / ( initialDimensions.at( 1 ) + internalDisplacements.at( 1 ) ) + a0 * At * Et * ( 1 / a0 - initialDimensions.at( 1 ) / ( a0 * ( initialDimensions.at( 1 ) + internalDisplacements.at( 1 ) ) ) ) / initialDimensions.at( 1 ) );
        ke.at( 1, 3 ) = ke.at( 3, 1 ) = ke.at( 5, 7 ) = ke.at( 7, 5 ) = - ( ny / 2 ) * ( At * Et / ( initialDimensions.at( 1 ) + internalDisplacements.at( 1 ) ) + a0 * At * Et * ( 1 / a0 - initialDimensions.at( 1 ) / ( a0 * ( initialDimensions.at( 1 ) + internalDisplacements.at( 1 ) ) ) ) / initialDimensions.at( 1 ) );
        ke.at( 2, 2 ) = ke.at( 4, 4 ) = ke.at( 6, 6 ) = ke.at( 8, 8 ) = ( nx / 2 ) * ( At * Et / ( initialDimensions.at( 2 ) + internalDisplacements.at( 2 ) ) + a0 * At * Et * ( 1 / a0 - initialDimensions.at( 2 ) / ( a0 * ( initialDimensions.at( 2 ) + internalDisplacements.at( 2 ) ) ) ) / initialDimensions.at( 2 ) );
        ke.at( 2, 8 ) = ke.at( 4, 6 ) = ke.at( 6, 4 ) = ke.at( 8, 4 ) = -( nx / 2 ) * ( At * Et / ( initialDimensions.at( 2 ) + internalDisplacements.at( 2 ) ) + a0 * At * Et * ( 1 / a0 - initialDimensions.at( 2 ) / ( a0 * ( initialDimensions.at( 2 ) + internalDisplacements.at( 2 ) ) ) ) / initialDimensions.at( 2 ) );
    }

    FloatMatrix T, TTke;
    T = calculateTransformationMatrix( tStep );
    TTke.beTProductOf( T, ke );
    answer.beProductOf( TTke, T );
}

void
NetRec4TrLa::giveDofManDofIDMask( int inode, IntArray &answer ) const
{
    answer = { D_u, D_v, D_w };
}

void
NetRec4TrLa::giveInternalForcesVector( FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord )
{
    FloatArray strainsInSidesRed, strainsInSides, stresses;
    double At = 0;
    // Calculate element's current dimensions in the local x- and y- direction
    FloatArray currentDimensions;
    currentDimensions = calculateElementDimensions( tStep );

    // Calculate strains in element's sides
    strainsInSides.resize( 4 );
    strainsInSides.at( 1 ) = ( currentDimensions.at( 1 ) - initialDimensions.at( 1 ) ) / initialDimensions.at( 1 );
    strainsInSides.at( 2 ) = ( currentDimensions.at( 2 ) - initialDimensions.at( 2 ) ) / initialDimensions.at( 2 );
    strainsInSides.at( 3 ) = ( currentDimensions.at( 3 ) - initialDimensions.at( 1 ) ) / initialDimensions.at( 1 );
    strainsInSides.at( 4 ) = ( currentDimensions.at( 4 ) - initialDimensions.at( 2 ) ) / initialDimensions.at( 2 );

    strainsInSidesRed.resize( 2 );
    //Find stresses in side 1-2
    strainsInSidesRed.at( 1 ) = strainsInSides.at( 1 );
    strainsInSidesRed.at( 2 ) = 0.0;
    for ( auto gp : *this->giveDefaultIntegrationRulePtr() ) {
        computeStressVector( stresses, strainsInSidesRed, gp, tStep, false );
        At = this->giveCrossSection()->give( CS_Area, gp );
    }
    // Internal force in twines coincident with side 1-2
    double Fx12 = At * stresses.at( 1 );

    // Find stresses in side 3-4
    strainsInSidesRed.at( 1 ) = strainsInSides.at( 3 );
    strainsInSidesRed.at( 2 ) = 0.0;
    for ( auto gp : *this->giveDefaultIntegrationRulePtr() )
        computeStressVector( stresses, strainsInSidesRed, gp, tStep, false );
    
    // Internal force in twines coincident with side 3-4
    double Fx34 = At * stresses.at( 1 );

    // Find stresses in side 1-4
    strainsInSidesRed.at( 1 ) = 0.0;
    strainsInSidesRed.at( 2 ) = strainsInSides.at( 2 );
    for ( auto gp : *this->giveDefaultIntegrationRulePtr() )
        computeStressVector( stresses, strainsInSidesRed, gp, tStep, false );

    // Internal force in twines coincident with side 1-4
    double Fy14 = At * stresses.at( 2 );

    // Find stresses in side 2-3
    strainsInSidesRed.at( 1 ) = 0.0;
    strainsInSidesRed.at( 2 ) = strainsInSides.at( 4 );
    for ( auto gp : *this->giveDefaultIntegrationRulePtr() )
        computeStressVector( stresses, strainsInSidesRed, gp, tStep, false );

    // Internal force in twines coincident with side 1-4
    double Fy23 = At * stresses.at( 2 );

    // Calculate internal forces at the nodes in the local coord. system
    FloatArray fe;
    fe.resize( 8 );
    fe.at( 1 ) = -Fx12 * ny / 2;
    fe.at( 3 ) = Fx12 * ny / 2;
    fe.at( 5 ) = Fx34 * ny / 2;
    fe.at( 7 ) = -Fx34 * ny / 2;
    fe.at( 2 ) = -Fy14 * nx / 2;
    fe.at( 8 ) = Fy14 * nx / 2;
    fe.at( 4 ) = -Fy23 * nx / 2;
    fe.at( 6 ) = Fy23 * nx / 2;

    FloatMatrix T;
    T = calculateTransformationMatrix( tStep );

    // Transform into the global coord. system
    answer.beTProductOf( T, fe );

    // Calculate strains and stresses to be saved
    FloatArray strains;
    for ( auto gp : *this->giveDefaultIntegrationRulePtr() ) {
        computeStrainVector( strains, gp, tStep );
        computeStressVector( stresses, strains, gp, tStep );
    }
}

FEInterpolation *
NetRec4TrLa::giveInterpolation() const
{
    return &interpolation;
}

double
NetRec4TrLa::giveNumberOfTwines() {
    if ( mask == 1 )
        return ( 2 * initialDimensions.at( 1 ) / a0 * initialDimensions.at( 2 ) / a0 );
    else
        OOFEM_ERROR( "Only mask 1 is implemented for NetRec4TrLa at this moment." );

    return 0;
}

double
NetRec4TrLa::giveTwineLength() {
    if ( mask == 1 )
        return a0;
    else
        OOFEM_ERROR( "Only mask 1 is implemented for NetRec4TrLa at this moment." );
    
    return 0;
}

void
NetRec4TrLa ::computeStrainVector( FloatArray &answer, GaussPoint *gp, TimeStep *tStep )
{
    answer.resize( 2 );
    // Calculate element's current dimensions in the local x- and y- direction
    FloatArray currentDimensions;
    currentDimensions = calculateElementDimensions( tStep );

    // Calculate strains in x- and y-twines
    answer.at( 1 ) = ( currentDimensions.at( 1 ) - initialDimensions.at( 1 ) ) / initialDimensions.at( 1 );
    answer.at( 2 ) = ( currentDimensions.at( 2 ) - initialDimensions.at( 2 ) ) / initialDimensions.at( 2 );
}

void
NetRec4TrLa ::computeStressVector( FloatArray &answer, const FloatArray &strains, GaussPoint *gp, TimeStep *tStep )
{
    answer = this->giveStructuralCrossSection()->giveRealStress_Netting( strains, gp, tStep );
}

void NetRec4TrLa ::computeStressVector( FloatArray &answer, const FloatArray &strains, GaussPoint *gp, TimeStep *tStep, bool saveContext )
{
    answer = this->giveStructuralCrossSection()->giveRealStress_Netting( strains, gp, tStep, saveContext );
}

void
NetRec4TrLa::initializeFrom( InputRecord &ir )
{
    StructuralElement::initializeFrom( ir );

    IR_GIVE_FIELD( ir, mask, _IFT_NetRec4TrLa_Mask );

    IR_GIVE_FIELD( ir, dx0, _IFT_NetRec4TrLa_dx0 );
    IR_GIVE_OPTIONAL_FIELD( ir, dy0, _IFT_NetRec4TrLa_dy0 );
    if ( dy0 = -1 ) {
        a0 = dx0;
        dx0 = -1;
    }

    IR_GIVE_OPTIONAL_FIELD( ir, nx, _IFT_NetRec4TrLa_nx );
    IR_GIVE_OPTIONAL_FIELD( ir, ny, _IFT_NetRec4TrLa_ny );
}
}