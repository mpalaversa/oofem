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

#include "sm/Elements/NetElements/nettr3pr.h"
#include "../sm/CrossSections/decoupledcrosssection.h"
#include "../sm/CrossSections/structuralcrosssection.h"
#include "sm/Materials/structuralms.h"

#include "boundaryload.h"
#include "classfactory.h"
#include "fei2dtrlin.h"
#include "gaussintegrationrule.h"
#include "gausspoint.h"
#include "load.h"

namespace oofem {
REGISTER_Element( NetTr3Pr );
FEI2dTrLin NetTr3Pr ::interp( 1, 2 );
NetTr3Pr::NetTr3Pr(int n, Domain* aDomain) : NetElement(n, aDomain)
{
    numberOfDofMans = 3;
    numberOfGaussPoints = 1;
    
    d  = -1;
    L0 = -1;
	U1 = -1;
	V1 = -1;
	U2 = -1;
	V2 = -1;
	U3 = -1;
	V3 = -1;

    srf = 1.0;
}

FloatArray
NetTr3Pr::calculateCurrentUnitNormalToElement( TimeStep *tStep )
{
    FloatArray U = computeUTwine( tStep );
    FloatArray V = computeVTwine( tStep );

    FloatArray en;
    en.beVectorProductOf( U, V );
    en.normalize();
    return en;
}

FloatArray
NetTr3Pr::calculateRelativeAcceleration( FloatArray acceleration, TimeStep *tStep )
{
    // Get acceleration of element nodes in the current time step
    FloatArray currentNodalAcceleration;
    this->computeVectorOf( VM_Acceleration, tStep, currentNodalAcceleration );

    // Calculate an average fluid acceleration on the element - this should be changed when the fluid acceleration becomes available as an input quantity
    FloatArray relativeAcceleration;
    relativeAcceleration.resize( 3 );
    relativeAcceleration.at( 1 ) = acceleration.at( 1 ) - ( currentNodalAcceleration.at( 1 ) + currentNodalAcceleration.at( 4 ) + currentNodalAcceleration.at( 7 ) ) / 3;
    relativeAcceleration.at( 2 ) = acceleration.at( 2 ) - ( currentNodalAcceleration.at( 2 ) + currentNodalAcceleration.at( 5 ) + currentNodalAcceleration.at( 8 ) ) / 3;
    relativeAcceleration.at( 3 ) = acceleration.at( 3 ) - ( currentNodalAcceleration.at( 3 ) + currentNodalAcceleration.at( 6 ) + currentNodalAcceleration.at( 9 ) ) / 3;

    return relativeAcceleration;
}

FloatArray
NetTr3Pr::calculateRelativeVelocity(FloatArray velocity, TimeStep* tStep) {
    // Get velocity and acceleration of element nodes in the current time step
    FloatArray currentNodalVelocity;
    this->computeVectorOf( VM_Velocity, tStep, currentNodalVelocity );

    // Calculate velocity of the fluid relative to the element
    FloatArray relativeVelocity;
    relativeVelocity.resize( 3 );
    relativeVelocity.at( 1 ) = velocity.at( 1 ) - ( currentNodalVelocity.at( 1 ) + currentNodalVelocity.at( 4 ) + currentNodalVelocity.at( 7 ) ) / 3;
    relativeVelocity.at( 2 ) = velocity.at( 2 ) - ( currentNodalVelocity.at( 2 ) + currentNodalVelocity.at( 5 ) + currentNodalVelocity.at( 8 ) ) / 3;
    relativeVelocity.at( 3 ) = velocity.at( 3 ) - ( currentNodalVelocity.at( 3 ) + currentNodalVelocity.at( 6 ) + currentNodalVelocity.at( 9 ) ) / 3;

    return relativeVelocity;
}

void
NetTr3Pr ::computeGaussPoints()
{
    // Sets up the integration rule array which contains all the Gauss points
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize( 1 );
        integrationRulesArray[0] = std::make_unique<GaussIntegrationRule>( 1, this, 1, 6 );
        this->giveCrossSection()->setupIntegrationPoints( *integrationRulesArray[0], this->numberOfGaussPoints, this );
    }
}

void NetTr3Pr ::computeHydrodynamicLoadVector( FloatArray &answer, FloatArray flowCharacteristics, TimeStep *tStep )
{
    // Form the fluid velocity vector
    FloatArray velocity;
    velocity.resize( 3 );
    for ( int i = 1; i <= 3; i++ )
        velocity.at( i ) = flowCharacteristics.at( i );

    DecoupledCrossSection *cs = this->giveDecoupledCrossSectionOfType( DecoupledMaterial::DecoupledMaterialType::DecoupledFluidMaterial );
    // Check if the element is downstream relative to another element
    if ( this->isDownstream ) {
        double sn = cs->giveSolidityRatio();
        // Reduce the inflow velocity by the velocity reduction coefficient as given in Loland, G. Current forces on and flow thorugh fish farms
        if ( sn > 0 )
            velocity.beScaled( 1 - 0.46 * ( 0.33 * sn + 6.54 * pow( sn, 2 ) - 4.88 * pow( sn, 3 ) ), velocity );
        else
            OOFEM_ERROR( "Element %d is denoted as downstream, but the solidity ratio is not specified.", this->giveNumber() );
    }

    FloatArray relativeVelocity = calculateRelativeVelocity( velocity, tStep );

    // Get nodal coordinates in the Oxyz and the O'UVW coord. system
    // Fetch initial coordinates of element's nodes
    FloatArray node1 = this->giveNode( 1 )->giveCoordinates();
    FloatArray node2 = this->giveNode( 2 )->giveCoordinates();
    FloatArray node3 = this->giveNode( 3 )->giveCoordinates();
    // Fetch total displacements in the current configuration
    FloatArray u;
    u.resize( 9 );
    if ( !tStep->isTheFirstStep() )
        this->computeVectorOf( VM_Total, tStep, u );
    // Current coordinates of node 1
    double x1 = node1.at( 1 ) + u.at( 1 );
    double y1 = node1.at( 2 ) + u.at( 2 );
    double z1 = node1.at( 3 ) + u.at( 3 );
    // Current coordinates of node 2
    double x2 = node2.at( 1 ) + u.at( 4 );
    double y2 = node2.at( 2 ) + u.at( 5 );
    double z2 = node2.at( 3 ) + u.at( 6 );
    // Current coordinates of node 3
    double x3 = node3.at( 1 ) + u.at( 7 );
    double y3 = node3.at( 2 ) + u.at( 8 );
    double z3 = node3.at( 3 ) + u.at( 9 );

    // Side vectors
    FloatArray s12, s13;
    s12.resize( 3 );
    s13.resize( 3 );

    s12.at( 1 ) = x2 - x1;
    s12.at( 2 ) = y2 - y1;
    s12.at( 3 ) = z2 - z1;

    s13.at( 1 ) = x3 - x1;
    s13.at( 2 ) = y3 - y1;
    s13.at( 3 ) = z3 - z1;

    // Twine vectors
    FloatArray U, U12, U13;
    U12.beScaled( ( V3 - V1 ) / d, s12 );
    U13.beScaled( ( V2 - V1 ) / d, s13 );
    U.beDifferenceOf( U12, U13 );

    FloatArray V, V12, V13;
    V12.beScaled( ( U3 - U1 ) / d, s12 );
    V13.beScaled( ( U2 - U1 ) / d, s13 );
    V.beDifferenceOf( V13, V12 );
    // Unit twine vectors
    FloatArray EU, EV;
    EU.beScaled( 1 / U.computeNorm(), U );
    EV.beScaled( 1 / V.computeNorm(), V );

    // Tangential relative velocity component in U- and V-twines
    FloatArray vRtU, vRtV;
    vRtU.beScaled( relativeVelocity.dotProduct( EU ), EU );
    vRtV.beScaled( relativeVelocity.dotProduct( EV ), EV );
    // Normal relative velocity component in U- and V-twines
    FloatArray vRnU, vRnV;
    vRnU.beDifferenceOf( relativeVelocity, vRtU );
    vRnV.beDifferenceOf( relativeVelocity, vRtV );

    // Fetch length of a twine, characteristic dimension of the cross-section and density and dynamic viscosity of the fluid
    double userDefinedDragCoeff = cs->giveDragCoefficient();
    double l                    = giveTwineLength();
    double density              = cs->giveMagnitudeOfMaterialProperty( 'd' );
    double mu                   = cs->giveMaterial()->giveDynamicViscosity();
    double characteristicDim    = cs->giveCharacteristicDimension();

    // Caluclate drag coefficients
    FloatArray dragCoeffsOnU, dragCoeffsOnV;
    if ( userDefinedDragCoeff == 0.0 ) {
        // If no, calculate the drag coefficients.
        dragCoeffsOnU = computeDragCoefficients( density, mu, characteristicDim, vRnU.computeNorm() );
        dragCoeffsOnV = computeDragCoefficients( density, mu, characteristicDim, vRnV.computeNorm() );
    } else {
        // If yes, use the user defined drag coefficient for the normal viscous force component
        // and disrfegard the tangential one.
        dragCoeffsOnU.resize( 2 );
        dragCoeffsOnV.resize( 2 );
        dragCoeffsOnU.at( 1 ) = dragCoeffsOnV.at( 1 ) = userDefinedDragCoeff;
    }

    // Calculate normal and tangential viscous force components on a U- and a V-twine
    FloatArray FRnU, FRtU, FRnV, FRtV;
    FRnU.beScaled( 0.5 * density * dragCoeffsOnU.at( 1 ) * l * characteristicDim * vRnU.computeNorm(), vRnU );
    FRnV.beScaled( 0.5 * density * dragCoeffsOnV.at( 1 ) * l * characteristicDim * vRnV.computeNorm(), vRnV );
    FRtU.beScaled( dragCoeffsOnU.at( 2 ) * l, vRtU );
    FRtV.beScaled( dragCoeffsOnV.at( 2 ) * l, vRtV );

    // The total viscous force on the element
    FRnV.add( FRtV );
    FRnU.add( FRtU );
    FRnU.add( FRnV );
    FloatArray Fv;
    Fv.beScaled( d / 2, FRnU );

    // Distribute Fv to the nodes (calculate nodal force contribution due to the viscous force)
    calculateEquivalentLumpedNodalValues( answer, Fv );
    
    // Form the fluid acceleration vector
    FloatArray acceleration;
    acceleration.resize( 3 );
    for ( int i = 4; i <= 6; i++ )
        acceleration.at( i - 3 ) = flowCharacteristics.at( i );

    FloatArray relativeAcceleration = calculateRelativeAcceleration( acceleration, tStep );

    if ( relativeAcceleration.computeNorm() != 0 ) {
        // Tangential relative and absolute acceleration component in U- and V-twines
        FloatArray aRtU, aRtV, atU, atV;
        aRtU.beScaled( relativeAcceleration.dotProduct( EU ), EU );
        aRtV.beScaled( relativeAcceleration.dotProduct( EV ), EV );
        atU.beScaled( acceleration.dotProduct( EU ), EU );
        atV.beScaled( acceleration.dotProduct( EV ), EV );

        // Normal relative and absolute acceleration component in U- and V-twines
        FloatArray aRnU, aRnV, anU, anV;
        aRnU.beDifferenceOf( relativeAcceleration, aRtU );
        aRnV.beDifferenceOf( relativeAcceleration, aRtV );
        anU.beDifferenceOf( acceleration, atU );
        anV.beDifferenceOf( acceleration, atV );

        // Added-mass coefficient
        double cm = cs->giveAddedMassCoefficient();

        // Added-mass force on a U- and V-twine
        FloatArray FAU, FAV, accelerationComponent;
        accelerationComponent.beScaled( cm, aRnU );
        anU.add( accelerationComponent );
        FAU.beScaled( density * ( pow( characteristicDim, 2 ) * 3.14 / 4 ) * l, anU );
        accelerationComponent.beScaled( cm, aRnV );
        anV.add( accelerationComponent );
        FAV.beScaled( density * ( pow( characteristicDim, 2 ) * 3.14 / 4 ) * l, anV );

        // Total added-mass force on the element
        FloatArray FA;
        FAU.add( FAV );
        FA.beScaled( d / 2, FAU );

        // The force is equally distributed to the element's nodes
        FloatArray distributedAddedMassForce;
        calculateEquivalentLumpedNodalValues( distributedAddedMassForce, FA );
        answer.add( distributedAddedMassForce );
    }
}

void
NetTr3Pr::computeLumpedMassMatrix( FloatMatrix &answer, TimeStep *tStep )
{
    answer.resize( 9, 9 );
    if ( !this->isActivated( tStep ) ) {
        return;
    }

    GaussPoint *gp    = integrationRulesArray[0]->getIntegrationPoint( 0 );
    double density    = this->giveStructuralCrossSection()->give( 'd', gp );
    // Number of twines in each direction is d/2. Thus, the total number of twines within the element is d.
    double thirdMass   = 0.333 * d * density * this->giveCrossSection()->give( CS_Area, gp ) * L0;
    answer.at( 1, 1 ) = answer.at( 2, 2 ) = answer.at( 3, 3 ) = answer.at( 4, 4 ) = answer.at( 5, 5 ) = answer.at( 6, 6 ) = answer.at( 7, 7 ) = answer.at( 8, 8 ) = answer.at( 9, 9 ) = thirdMass;
}

void
NetTr3Pr::computeStiffnessMatrix(FloatMatrix& answer, MatResponseMode rMode, TimeStep* tStep) {
    FloatMatrix D;
    // Area and modulus of elasticity of a twine
    double At, Et;
    // Fetch At and Et from the user's input data
    for ( auto gp : *this->giveDefaultIntegrationRulePtr() ) {
        this->computeConstitutiveMatrixAt( D, rMode, gp, tStep );
        At = this->giveCrossSection()->give( CS_Area, gp );
    }
    Et = D.at( 1, 1 );

    // Fetch initial coordinates of element's nodes
    FloatArray node1 = this->giveNode( 1 )->giveCoordinates();
    FloatArray node2 = this->giveNode( 2 )->giveCoordinates();
    FloatArray node3 = this->giveNode( 3 )->giveCoordinates();
    // Fetch total displacements in the current configuration
    FloatArray u;
    u.resize( 9 );
    if ( !tStep->isTheFirstStep() )
        this->computeVectorOf( VM_Total, tStep, u );
    // Current coordinates of node 1
    double x1 = node1.at( 1 ) + u.at( 1 );
    double y1 = node1.at( 2 ) + u.at( 2 );
    double z1 = node1.at( 3 ) + u.at( 3 );
    // Current coordinates of node 2
    double x2 = node2.at( 1 ) + u.at( 4 );
    double y2 = node2.at( 2 ) + u.at( 5 );
    double z2 = node2.at( 3 ) + u.at( 6 );
    // Current coordinates of node 3
    double x3 = node3.at( 1 ) + u.at( 7 );
    double y3 = node3.at( 2 ) + u.at( 8 );
    double z3 = node3.at( 3 ) + u.at( 9 );

    FloatArray U = computeUTwine( tStep );
    FloatArray V = computeVTwine( tStep );

    double dUxdx1, dUydy1, dUzdz1, dUxdx2, dUydy2, dUzdz2, dUxdx3, dUydy3, dUzdz3, dUxdy1, dUxdy2, dUxdy3, dUxdz1, dUxdz2, dUxdz3, dUydz1, dUydz2, dUydz3, dUydx1, dUydx2, dUydx3, dUzdx1, dUzdx2, dUzdx3, dUzdy1, dUzdy2, dUzdy3;
    dUxdx1 = dUydy1 = dUzdz1 = ( V2 - V3 ) / d;
    dUxdx2 = dUydy2 = dUzdz2 = ( V3 - V1 ) / d;
    dUxdx3 = dUydy3 = dUzdz3 = ( V1 - V2 ) / d;
    dUxdy1 = dUxdy2 = dUxdy3 = dUxdz1 = dUxdz2 = dUxdz3 = dUydz1 = dUydz2 = dUydz3 = dUydx1 = dUydx2 = dUydx3 = dUzdx1 = dUzdx2 = dUzdx3 = dUzdy1 = dUzdy2 = dUzdy3 = 0;

    double dVxdx1, dVydy1, dVzdz1, dVxdx2, dVydy2, dVzdz2, dVxdx3, dVydy3, dVzdz3, dVxdy1, dVxdy2, dVxdy3, dVxdz1, dVxdz2, dVxdz3, dVydx1, dVydx2, dVydx3, dVydz1, dVydz2, dVydz3, dVzdx1, dVzdx2, dVzdx3, dVzdy1, dVzdy2, dVzdy3;
    dVxdx1 = dVydy1 = dVzdz1 = ( U3 - U2 ) / d;
    dVxdx2 = dVydy2 = dVzdz2 = ( U1 - U3 ) / d;
    dVxdx3 = dVydy3 = dVzdz3 = ( U2 - U1 ) / d;
    dVxdy1 = dVxdy2 = dVxdy3 = dVxdz1 = dVxdz2 = dVxdz3 = dVydx1 = dVydx2 = dVydx3 = dVydz1 = dVydz2 = dVydz3 = dVzdx1 = dVzdx2 = dVzdx3 = dVzdy1 = dVzdy2 = dVzdy3 = 0;

    double dUdx1, dUdx2, dUdx3, dUdy1, dUdy2, dUdy3, dUdz1, dUdz2, dUdz3; 
    dUdx1 = ( V2 - V3 ) / pow( d, 2 ) * ( ( x2 - x1 ) * ( V3 - V1 ) - ( x3 - x1 ) * ( V2 - V1 ) );
    dUdx2 = ( V3 - V1 ) / pow( d, 2 ) * ( ( x2 - x1 ) * ( V3 - V1 ) - ( x3 - x1 ) * ( V2 - V1 ) );
    dUdx3 = ( V1 - V2 ) / pow( d, 2 ) * ( ( x2 - x1 ) * ( V3 - V1 ) - ( x3 - x1 ) * ( V2 - V1 ) );
    dUdy1 = ( V2 - V3 ) / pow( d, 2 ) * ( ( y2 - y1 ) * ( V3 - V1 ) - ( y3 - y1 ) * ( V2 - V1 ) );
    dUdy2 = ( V3 - V1 ) / pow( d, 2 ) * ( ( y2 - y1 ) * ( V3 - V1 ) - ( y3 - y1 ) * ( V2 - V1 ) );
    dUdy3 = ( V1 - V2 ) / pow( d, 2 ) * ( ( y2 - y1 ) * ( V3 - V1 ) - ( y3 - y1 ) * ( V2 - V1 ) );
    dUdz1 = ( V2 - V3 ) / pow( d, 2 ) * ( ( z2 - z1 ) * ( V3 - V1 ) - ( z3 - z1 ) * ( V2 - V1 ) );
    dUdz2 = ( V3 - V1 ) / pow( d, 2 ) * ( ( z2 - z1 ) * ( V3 - V1 ) - ( z3 - z1 ) * ( V2 - V1 ) );
    dUdz3 = ( V1 - V2 ) / pow( d, 2 ) * ( ( z2 - z1 ) * ( V3 - V1 ) - ( z3 - z1 ) * ( V2 - V1 ) );

    double dVdx1, dVdx2, dVdx3, dVdy1, dVdy2, dVdy3, dVdz1, dVdz2, dVdz3;
    dVdx1 = ( U2 - U3 ) / pow( d, 2 ) * ( ( x2 - x1 ) * ( U3 - U1 ) - ( x3 - x1 ) * ( U2 - U1 ) );
    dVdx2 = ( U3 - U1 ) / pow( d, 2 ) * ( ( x2 - x1 ) * ( U3 - U1 ) - ( x3 - x1 ) * ( U2 - U1 ) );
    dVdx3 = ( U1 - U2 ) / pow( d, 2 ) * ( ( x2 - x1 ) * ( U3 - U1 ) - ( x3 - x1 ) * ( U2 - U1 ) );
    dVdy1 = ( U2 - U3 ) / pow( d, 2 ) * ( ( y2 - y1 ) * ( U3 - U1 ) - ( y3 - y1 ) * ( U2 - U1 ) );
    dVdy2 = ( U3 - U1 ) / pow( d, 2 ) * ( ( y2 - y1 ) * ( U3 - U1 ) - ( y3 - y1 ) * ( U2 - U1 ) );
    dVdy3 = ( U1 - U2 ) / pow( d, 2 ) * ( ( y2 - y1 ) * ( U3 - U1 ) - ( y3 - y1 ) * ( U2 - U1 ) );
    dVdz1 = ( U2 - U3 ) / pow( d, 2 ) * ( ( z2 - z1 ) * ( U3 - U1 ) - ( z3 - z1 ) * ( U2 - U1 ) );
    dVdz2 = ( U3 - U1 ) / pow( d, 2 ) * ( ( z2 - z1 ) * ( U3 - U1 ) - ( z3 - z1 ) * ( U2 - U1 ) );
    dVdz3 = ( U1 - U2 ) / pow( d, 2 ) * ( ( z2 - z1 ) * ( U3 - U1 ) - ( z3 - z1 ) * ( U2 - U1 ) );

    bool reduceUStiffness = false, reduceVStiffness = false;
    if ( L0 > U.computeNorm() )
        reduceUStiffness = true;
    if ( L0 > V.computeNorm() )
        reduceVStiffness = true;

    // Calculate elements of the stiffness matrix. This stiffness matrix is already in the global coord. system
    answer.resize( 9, 9 );
    //First row in the matrix
        //dFx1dx1
    if ( srf == 1.0 || ( !reduceUStiffness && !reduceVStiffness ) ) {
        answer.at( 1, 1 ) = Et * At * ( V3 - V2 ) / 2 * ( dUxdx1 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx1 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVxdx1 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx1 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness && reduceVStiffness ) {
        answer.at( 1, 1 ) = ( Et / srf ) * At * ( V3 - V2 ) / 2 * ( dUxdx1 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx1 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U2 - U3 ) / 2 * ( dVxdx1 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx1 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceVStiffness ) {
        answer.at( 1, 1 ) = Et * At * ( V3 - V2 ) / 2 * ( dUxdx1 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx1 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U2 - U3 ) / 2 * ( dVxdx1 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx1 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness) {
        answer.at( 1, 1 ) = ( Et / srf ) * At * ( V3 - V2 ) / 2 * ( dUxdx1 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx1 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVxdx1 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx1 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else {
        OOFEM_ERROR( "No stiffness is associated with element (1, 1) of the stiffness matrix of the finite element %d.", this->giveNumber() );
    }
        //dFx1dy1
    //answer.at( 1, 2 ) = Et * At * ( V3 - V2 ) / 2 * ( dUxdy1 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy1 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVxdy1 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy1 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    if ( srf == 1.0 || ( !reduceUStiffness && !reduceVStiffness ) ) {
        answer.at( 1, 2 ) = Et * At * ( V3 - V2 ) / 2 * ( dUxdy1 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy1 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVxdy1 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy1 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness && reduceVStiffness ) {
        answer.at( 1, 2 ) = ( Et / srf ) * At * ( V3 - V2 ) / 2 * ( dUxdy1 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy1 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U2 - U3 ) / 2 * ( dVxdy1 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy1 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceVStiffness ) {
        answer.at( 1, 2 ) = Et * At * ( V3 - V2 ) / 2 * ( dUxdy1 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy1 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U2 - U3 ) / 2 * ( dVxdy1 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy1 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness ) {
        answer.at( 1, 2 ) = ( Et / srf ) * At * ( V3 - V2 ) / 2 * ( dUxdy1 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy1 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVxdy1 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy1 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else {
        OOFEM_ERROR( "No stiffness is associated with element (1, 2) of the stiffness matrix of the finite element %d.", this->giveNumber() );
    }
        //dFx1dz1
    //answer.at( 1, 3 ) = Et * At * ( V3 - V2 ) / 2 * ( dUxdz1 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz1 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVxdz1 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz1 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    if ( srf == 1.0 || ( !reduceUStiffness && !reduceVStiffness ) ) {
        answer.at( 1, 3 ) = Et * At * ( V3 - V2 ) / 2 * ( dUxdz1 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz1 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVxdz1 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz1 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness && reduceVStiffness ) {
        answer.at( 1, 3 ) = ( Et / srf ) * At * ( V3 - V2 ) / 2 * ( dUxdz1 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz1 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U2 - U3 ) / 2 * ( dVxdz1 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz1 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceVStiffness ) {
        answer.at( 1, 3 ) = Et * At * ( V3 - V2 ) / 2 * ( dUxdz1 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz1 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U2 - U3 ) / 2 * ( dVxdz1 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz1 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness ) {
        answer.at( 1, 3 ) = ( Et / srf ) * At * ( V3 - V2 ) / 2 * ( dUxdz1 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz1 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVxdz1 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz1 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else {
        OOFEM_ERROR( "No stiffness is associated with element (1, 3) of the stiffness matrix of the finite element %d.", this->giveNumber() );
    }
        //dFx1dx2
    //answer.at( 1, 4 ) = Et * At * ( V3 - V2 ) / 2 * ( dUxdx2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx2 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVxdx2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx2 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    if ( srf == 1.0 || ( !reduceUStiffness && !reduceVStiffness ) ) {
        answer.at( 1, 4 ) = Et * At * ( V3 - V2 ) / 2 * ( dUxdx2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx2 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVxdx2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx2 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness && reduceVStiffness ) {
        answer.at( 1, 4 ) = ( Et / srf ) * At * ( V3 - V2 ) / 2 * ( dUxdx2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx2 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U2 - U3 ) / 2 * ( dVxdx2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx2 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceVStiffness ) {
        answer.at( 1, 4 ) = Et * At * ( V3 - V2 ) / 2 * ( dUxdx2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx2 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U2 - U3 ) / 2 * ( dVxdx2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx2 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness ) {
        answer.at( 1, 4 ) = ( Et / srf ) * At * ( V3 - V2 ) / 2 * ( dUxdx2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx2 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVxdx2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx2 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else {
        OOFEM_ERROR( "No stiffness is associated with element (1, 4) of the stiffness matrix of the finite element %d.", this->giveNumber() );
    }
        //dFx1dy2
    //answer.at( 1, 5 ) = Et * At * ( V3 - V2 ) / 2 * ( dUxdy2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy2 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVxdy2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy2 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    if ( srf == 1.0 || ( !reduceUStiffness && !reduceVStiffness ) ) {
        answer.at( 1, 5 ) = Et * At * ( V3 - V2 ) / 2 * ( dUxdy2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy2 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVxdy2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy2 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness && reduceVStiffness ) {
        answer.at( 1, 5 ) = ( Et / srf ) * At * ( V3 - V2 ) / 2 * ( dUxdy2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy2 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U2 - U3 ) / 2 * ( dVxdy2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy2 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceVStiffness ) {
        answer.at( 1, 5 ) = Et * At * ( V3 - V2 ) / 2 * ( dUxdy2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy2 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U2 - U3 ) / 2 * ( dVxdy2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy2 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness ) {
        answer.at( 1, 5 ) = ( Et / srf ) * At * ( V3 - V2 ) / 2 * ( dUxdy2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy2 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVxdy2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy2 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else {
        OOFEM_ERROR( "No stiffness is associated with element (1, 5) of the stiffness matrix of the finite element %d.", this->giveNumber() );
    }
        //dFx1dz2
    //answer.at( 1, 6 ) = Et * At * ( V3 - V2 ) / 2 * ( dUxdz2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz2 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVxdz2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz2 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    if ( srf == 1.0 || ( !reduceUStiffness && !reduceVStiffness ) ) {
        answer.at( 1, 6 ) = Et * At * ( V3 - V2 ) / 2 * ( dUxdz2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz2 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVxdz2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz2 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness && reduceVStiffness ) {
        answer.at( 1, 6 ) = ( Et / srf ) * At * ( V3 - V2 ) / 2 * ( dUxdz2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz2 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U2 - U3 ) / 2 * ( dVxdz2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz2 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceVStiffness ) {
        answer.at( 1, 6 ) = Et * At * ( V3 - V2 ) / 2 * ( dUxdz2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz2 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U2 - U3 ) / 2 * ( dVxdz2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz2 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness ) {
        answer.at( 1, 6 ) = ( Et / srf ) * At * ( V3 - V2 ) / 2 * ( dUxdz2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz2 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVxdz2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz2 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else {
        OOFEM_ERROR( "No stiffness is associated with element (1, 6) of the stiffness matrix of the finite element %d.", this->giveNumber() );
    }
        //dFx1dx3
    //answer.at( 1, 7 ) = Et * At * ( V3 - V2 ) / 2 * ( dUxdx3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx3 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVxdx3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx3 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    if ( srf == 1.0 || ( !reduceUStiffness && !reduceVStiffness ) ) {
        answer.at( 1, 7 ) = Et * At * ( V3 - V2 ) / 2 * ( dUxdx3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx3 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVxdx3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx3 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness && reduceVStiffness ) {
        answer.at( 1, 7 ) = ( Et / srf ) * At * ( V3 - V2 ) / 2 * ( dUxdx3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx3 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U2 - U3 ) / 2 * ( dVxdx3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx3 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceVStiffness ) {
        answer.at( 1, 7 ) = Et * At * ( V3 - V2 ) / 2 * ( dUxdx3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx3 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U2 - U3 ) / 2 * ( dVxdx3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx3 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness ) {
        answer.at( 1, 7 ) = ( Et / srf ) * At * ( V3 - V2 ) / 2 * ( dUxdx3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx3 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVxdx3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx3 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else {
        OOFEM_ERROR( "No stiffness is associated with element (1, 7) of the stiffness matrix of the finite element %d.", this->giveNumber() );
    }
        //dFx1dy3
    //answer.at( 1, 8 ) = Et * At * ( V3 - V2 ) / 2 * ( dUxdy3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy3 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVxdy3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy3 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    if ( srf == 1.0 || ( !reduceUStiffness && !reduceVStiffness ) ) {
        answer.at( 1, 8 ) = Et * At * ( V3 - V2 ) / 2 * ( dUxdy3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy3 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVxdy3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy3 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness && reduceVStiffness ) {
        answer.at( 1, 8 ) = ( Et / srf ) * At * ( V3 - V2 ) / 2 * ( dUxdy3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy3 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U2 - U3 ) / 2 * ( dVxdy3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy3 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceVStiffness ) {
        answer.at( 1, 8 ) = Et * At * ( V3 - V2 ) / 2 * ( dUxdy3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy3 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U2 - U3 ) / 2 * ( dVxdy3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy3 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness ) {
        answer.at( 1, 8 ) = ( Et / srf ) * At * ( V3 - V2 ) / 2 * ( dUxdy3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy3 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVxdy3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy3 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else {
        OOFEM_ERROR( "No stiffness is associated with element (1, 8) of the stiffness matrix of the finite element %d.", this->giveNumber() );
    }
        //dFx1dz3
    //answer.at( 1, 9 ) = Et * At * ( V3 - V2 ) / 2 * ( dUxdz3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz3 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVxdz3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz3 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    if ( srf == 1.0 || ( !reduceUStiffness && !reduceVStiffness ) ) {
        answer.at( 1, 9 ) = Et * At * ( V3 - V2 ) / 2 * ( dUxdz3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz3 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVxdz3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz3 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness && reduceVStiffness ) {
        answer.at( 1, 9 ) = ( Et / srf ) * At * ( V3 - V2 ) / 2 * ( dUxdz3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz3 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U2 - U3 ) / 2 * ( dVxdz3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz3 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceVStiffness ) {
        answer.at( 1, 9 ) = Et * At * ( V3 - V2 ) / 2 * ( dUxdz3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz3 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U2 - U3 ) / 2 * ( dVxdz3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz3 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness ) {
        answer.at( 1, 9 ) = ( Et / srf ) * At * ( V3 - V2 ) / 2 * ( dUxdz3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz3 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVxdz3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz3 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else {
        OOFEM_ERROR( "No stiffness is associated with element (1, 9) of the stiffness matrix of the finite element %d.", this->giveNumber() );
    }

    //Second row in the matrix
        //dFy1dy1
    //answer.at( 2, 2 ) = Et * At * ( V3 - V2 ) / 2 * ( dUydy1 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy1 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVydy1 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy1 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    if ( srf == 1.0 || ( !reduceUStiffness && !reduceVStiffness ) ) {
        answer.at( 2, 2 ) = Et * At * ( V3 - V2 ) / 2 * ( dUydy1 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy1 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVydy1 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy1 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness && reduceVStiffness ) {
        answer.at( 2, 2 ) = ( Et / srf ) * At * ( V3 - V2 ) / 2 * ( dUydy1 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy1 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U2 - U3 ) / 2 * ( dVydy1 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy1 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceVStiffness ) {
        answer.at( 2, 2 ) = Et * At * ( V3 - V2 ) / 2 * ( dUydy1 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy1 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U2 - U3 ) / 2 * ( dVydy1 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy1 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness ) {
        answer.at( 2, 2 ) = ( Et / srf ) * At * ( V3 - V2 ) / 2 * ( dUydy1 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy1 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVydy1 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy1 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else {
        OOFEM_ERROR( "No stiffness is associated with element (2, 2) of the stiffness matrix of the finite element %d.", this->giveNumber() );
    }
        //dFy1dz1
    //answer.at( 2, 3 ) = Et * At * ( V3 - V2 ) / 2 * ( dUydz1 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz1 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVydz1 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz1 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    if ( srf == 1.0 || ( !reduceUStiffness && !reduceVStiffness ) ) {
        answer.at( 2, 3 ) = Et * At * ( V3 - V2 ) / 2 * ( dUydz1 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz1 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVydz1 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz1 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness && reduceVStiffness ) {
        answer.at( 2, 3 ) = ( Et / srf ) * At * ( V3 - V2 ) / 2 * ( dUydz1 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz1 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U2 - U3 ) / 2 * ( dVydz1 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz1 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceVStiffness ) {
        answer.at( 2, 3 ) = Et * At * ( V3 - V2 ) / 2 * ( dUydz1 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz1 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U2 - U3 ) / 2 * ( dVydz1 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz1 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness ) {
        answer.at( 2, 3 ) = ( Et / srf ) * At * ( V3 - V2 ) / 2 * ( dUydz1 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz1 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVydz1 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz1 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else {
        OOFEM_ERROR( "No stiffness is associated with element (2, 3) of the stiffness matrix of the finite element %d.", this->giveNumber() );
    }
        //dFy1dx2
    //answer.at( 2, 4 ) = Et * At * ( V3 - V2 ) / 2 * ( dUydx2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx2 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVydx2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx2 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    if ( srf == 1.0 || ( !reduceUStiffness && !reduceVStiffness ) ) {
        answer.at( 2, 4 ) = Et * At * ( V3 - V2 ) / 2 * ( dUydx2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx2 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVydx2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx2 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness && reduceVStiffness ) {
        answer.at( 2, 4 ) = ( Et / srf ) * At * ( V3 - V2 ) / 2 * ( dUydx2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx2 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U2 - U3 ) / 2 * ( dVydx2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx2 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceVStiffness ) {
        answer.at( 2, 4 ) = Et * At * ( V3 - V2 ) / 2 * ( dUydx2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx2 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U2 - U3 ) / 2 * ( dVydx2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx2 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness ) {
        answer.at( 2, 4 ) = ( Et / srf ) * At * ( V3 - V2 ) / 2 * ( dUydx2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx2 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVydx2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx2 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else {
        OOFEM_ERROR( "No stiffness is associated with element (2, 4) of the stiffness matrix of the finite element %d.", this->giveNumber() );
    }
        //dFy1dy2
    //answer.at( 2, 5 ) = Et * At * ( V3 - V2 ) / 2 * ( dUydy2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy2 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVydy2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy2 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    if ( srf == 1.0 || ( !reduceUStiffness && !reduceVStiffness ) ) {
        answer.at( 2, 5 ) = Et * At * ( V3 - V2 ) / 2 * ( dUydy2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy2 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVydy2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy2 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness && reduceVStiffness ) {
        answer.at( 2, 5 ) = ( Et / srf ) * At * ( V3 - V2 ) / 2 * ( dUydy2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy2 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U2 - U3 ) / 2 * ( dVydy2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy2 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceVStiffness ) {
        answer.at( 2, 5 ) = Et * At * ( V3 - V2 ) / 2 * ( dUydy2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy2 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U2 - U3 ) / 2 * ( dVydy2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy2 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness ) {
        answer.at( 2, 5 ) = ( Et / srf ) * At * ( V3 - V2 ) / 2 * ( dUydy2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy2 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVydy2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy2 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else {
        OOFEM_ERROR( "No stiffness is associated with element (2, 5) of the stiffness matrix of the finite element %d.", this->giveNumber() );
    }
        //dFy1dz2
    //answer.at( 2, 6 ) = Et * At * ( V3 - V2 ) / 2 * ( dUydz2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz2 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVydz2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz2 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    if ( srf == 1.0 || ( !reduceUStiffness && !reduceVStiffness ) ) {
        answer.at( 2, 6 ) = Et * At * ( V3 - V2 ) / 2 * ( dUydz2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz2 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVydz2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz2 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness && reduceVStiffness ) {
        answer.at( 2, 6 ) = ( Et / srf ) * At * ( V3 - V2 ) / 2 * ( dUydz2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz2 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U2 - U3 ) / 2 * ( dVydz2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz2 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceVStiffness ) {
        answer.at( 2, 6 ) = Et * At * ( V3 - V2 ) / 2 * ( dUydz2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz2 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U2 - U3 ) / 2 * ( dVydz2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz2 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness ) {
        answer.at( 2, 6 ) = ( Et / srf ) * At * ( V3 - V2 ) / 2 * ( dUydz2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz2 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVydz2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz2 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else {
        OOFEM_ERROR( "No stiffness is associated with element (2, 6) of the stiffness matrix of the finite element %d.", this->giveNumber() );
    }
        //dFy1dx3
    //answer.at( 2, 7 ) = Et * At * ( V3 - V2 ) / 2 * ( dUydx3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx3 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVydx3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx3 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    if ( srf == 1.0 || ( !reduceUStiffness && !reduceVStiffness ) ) {
        answer.at( 2, 7 ) = Et * At * ( V3 - V2 ) / 2 * ( dUydx3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx3 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVydx3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx3 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness && reduceVStiffness ) {
        answer.at( 2, 7 ) = ( Et / srf ) * At * ( V3 - V2 ) / 2 * ( dUydx3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx3 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U2 - U3 ) / 2 * ( dVydx3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx3 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceVStiffness ) {
        answer.at( 2, 7 ) = Et * At * ( V3 - V2 ) / 2 * ( dUydx3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx3 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U2 - U3 ) / 2 * ( dVydx3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx3 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness ) {
        answer.at( 2, 7 ) = ( Et / srf ) * At * ( V3 - V2 ) / 2 * ( dUydx3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx3 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVydx3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx3 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else {
        OOFEM_ERROR( "No stiffness is associated with element (2, 7) of the stiffness matrix of the finite element %d.", this->giveNumber() );
    }
        //dFy1dy3
    //answer.at( 2, 8 ) = Et * At * ( V3 - V2 ) / 2 * ( dUydy3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy3 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVydy3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy3 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    if ( srf == 1.0 || ( !reduceUStiffness && !reduceVStiffness ) ) {
        answer.at( 2, 8 ) = Et * At * ( V3 - V2 ) / 2 * ( dUydy3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy3 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVydy3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy3 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness && reduceVStiffness ) {
        answer.at( 2, 8 ) = ( Et / srf ) * At * ( V3 - V2 ) / 2 * ( dUydy3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy3 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U2 - U3 ) / 2 * ( dVydy3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy3 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceVStiffness ) {
        answer.at( 2, 8 ) = Et * At * ( V3 - V2 ) / 2 * ( dUydy3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy3 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U2 - U3 ) / 2 * ( dVydy3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy3 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness ) {
        answer.at( 2, 8 ) = ( Et / srf ) * At * ( V3 - V2 ) / 2 * ( dUydy3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy3 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVydy3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy3 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else {
        OOFEM_ERROR( "No stiffness is associated with element (2, 8) of the stiffness matrix of the finite element %d.", this->giveNumber() );
    }
        //dFy1dz3
    //answer.at( 2, 9 ) = Et * At * ( V3 - V2 ) / 2 * ( dUydz3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz3 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVydz3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz3 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    if ( srf == 1.0 || ( !reduceUStiffness && !reduceVStiffness ) ) {
        answer.at( 2, 9 ) = Et * At * ( V3 - V2 ) / 2 * ( dUydz3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz3 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVydz3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz3 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness && reduceVStiffness ) {
        answer.at( 2, 9 ) = ( Et / srf ) * At * ( V3 - V2 ) / 2 * ( dUydz3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz3 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U2 - U3 ) / 2 * ( dVydz3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz3 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceVStiffness ) {
        answer.at( 2, 9 ) = Et * At * ( V3 - V2 ) / 2 * ( dUydz3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz3 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U2 - U3 ) / 2 * ( dVydz3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz3 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness ) {
        answer.at( 2, 9 ) = ( Et / srf ) * At * ( V3 - V2 ) / 2 * ( dUydz3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz3 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVydz3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz3 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else {
        OOFEM_ERROR( "No stiffness is associated with element (2, 9) of the stiffness matrix of the finite element %d.", this->giveNumber() );
    }

    //Third row in the matrix
        //dFz1dz1
    //answer.at( 3, 3 ) = Et * At * ( V3 - V2 ) / 2 * ( dUzdz1 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz1 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVzdz1 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz1 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    if ( srf == 1.0 || ( !reduceUStiffness && !reduceVStiffness ) ) {
        answer.at( 3, 3 ) = Et * At * ( V3 - V2 ) / 2 * ( dUzdz1 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz1 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVzdz1 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz1 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness && reduceVStiffness ) {
        answer.at( 3, 3 ) = ( Et / srf ) * At * ( V3 - V2 ) / 2 * ( dUzdz1 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz1 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U2 - U3 ) / 2 * ( dVzdz1 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz1 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceVStiffness ) {
        answer.at( 3, 3 ) = Et * At * ( V3 - V2 ) / 2 * ( dUzdz1 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz1 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U2 - U3 ) / 2 * ( dVzdz1 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz1 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness ) {
        answer.at( 3, 3 ) = ( Et / srf ) * At * ( V3 - V2 ) / 2 * ( dUzdz1 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz1 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVzdz1 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz1 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    } else {
        OOFEM_ERROR( "No stiffness is associated with element (3, 3) of the stiffness matrix of the finite element %d.", this->giveNumber() );
    }
        //dFz1dx2
    //answer.at( 3, 4 ) = Et * At * ( V3 - V2 ) / 2 * ( dUzdx2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx2 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVzdx2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx2 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    if ( srf == 1.0 || ( !reduceUStiffness && !reduceVStiffness ) ) {
        answer.at( 3, 4 ) = Et * At * ( V3 - V2 ) / 2 * ( dUzdx2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx2 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVzdx2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx2 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness && reduceVStiffness ) {
        answer.at( 3, 4 ) = ( Et / srf ) * At * ( V3 - V2 ) / 2 * ( dUzdx2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx2 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U2 - U3 ) / 2 * ( dVzdx2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx2 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceVStiffness ) {
        answer.at( 3, 4 ) = Et * At * ( V3 - V2 ) / 2 * ( dUzdx2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx2 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U2 - U3 ) / 2 * ( dVzdx2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx2 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness ) {
        answer.at( 3, 4 ) = ( Et / srf ) * At * ( V3 - V2 ) / 2 * ( dUzdx2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx2 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVzdx2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx2 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    } else {
        OOFEM_ERROR( "No stiffness is associated with element (3, 4) of the stiffness matrix of the finite element %d.", this->giveNumber() );
    }
        //dFz1dy2
    //answer.at( 3, 5 ) = Et * At * ( V3 - V2 ) / 2 * ( dUzdy2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy2 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVzdy2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy2 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    if ( srf == 1.0 || ( !reduceUStiffness && !reduceVStiffness ) ) {
        answer.at( 3, 5 ) = Et * At * ( V3 - V2 ) / 2 * ( dUzdy2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy2 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVzdy2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy2 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness && reduceVStiffness ) {
        answer.at( 3, 5 ) = ( Et / srf ) * At * ( V3 - V2 ) / 2 * ( dUzdy2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy2 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U2 - U3 ) / 2 * ( dVzdy2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy2 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceVStiffness ) {
        answer.at( 3, 5 ) = Et * At * ( V3 - V2 ) / 2 * ( dUzdy2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy2 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U2 - U3 ) / 2 * ( dVzdy2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy2 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness ) {
        answer.at( 3, 5 ) = ( Et / srf ) * At * ( V3 - V2 ) / 2 * ( dUzdy2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy2 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVzdy2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy2 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    } else {
        OOFEM_ERROR( "No stiffness is associated with element (3, 5) of the stiffness matrix of the finite element %d.", this->giveNumber() );
    }
        //dFz1dz2
    //answer.at( 3, 6 ) = Et * At * ( V3 - V2 ) / 2 * ( dUzdz2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz2 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVzdz2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz2 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    if ( srf == 1.0 || ( !reduceUStiffness && !reduceVStiffness ) ) {
        answer.at( 3, 6 ) = Et * At * ( V3 - V2 ) / 2 * ( dUzdz2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz2 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVzdz2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz2 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness && reduceVStiffness ) {
        answer.at( 3, 6 ) = ( Et / srf ) * At * ( V3 - V2 ) / 2 * ( dUzdz2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz2 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U2 - U3 ) / 2 * ( dVzdz2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz2 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceVStiffness ) {
        answer.at( 3, 6 ) = Et * At * ( V3 - V2 ) / 2 * ( dUzdz2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz2 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U2 - U3 ) / 2 * ( dVzdz2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz2 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness ) {
        answer.at( 3, 6 ) = ( Et / srf ) * At * ( V3 - V2 ) / 2 * ( dUzdz2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz2 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVzdz2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz2 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    } else {
        OOFEM_ERROR( "No stiffness is associated with element (3, 6) of the stiffness matrix of the finite element %d.", this->giveNumber() );
    }
        //dFz1dx3
    //answer.at( 3, 7 ) = Et * At * ( V3 - V2 ) / 2 * ( dUzdx3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx3 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVzdx3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx3 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    if ( srf == 1.0 || ( !reduceUStiffness && !reduceVStiffness ) ) {
        answer.at( 3, 7 ) = Et * At * ( V3 - V2 ) / 2 * ( dUzdx3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx3 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVzdx3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx3 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness && reduceVStiffness ) {
        answer.at( 3, 7 ) = ( Et / srf ) * At * ( V3 - V2 ) / 2 * ( dUzdx3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx3 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U2 - U3 ) / 2 * ( dVzdx3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx3 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceVStiffness ) {
        answer.at( 3, 7 ) = Et * At * ( V3 - V2 ) / 2 * ( dUzdx3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx3 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U2 - U3 ) / 2 * ( dVzdx3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx3 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness ) {
        answer.at( 3, 7 ) = ( Et / srf ) * At * ( V3 - V2 ) / 2 * ( dUzdx3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx3 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVzdx3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx3 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    } else {
        OOFEM_ERROR( "No stiffness is associated with element (3, 7) of the stiffness matrix of the finite element %d.", this->giveNumber() );
    }
        //dFz1dy3
    //answer.at( 3, 8 ) = Et * At * ( V3 - V2 ) / 2 * ( dUzdy3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy3 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVzdy3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy3 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    if ( srf == 1.0 || ( !reduceUStiffness && !reduceVStiffness ) ) {
        answer.at( 3, 8 ) = Et * At * ( V3 - V2 ) / 2 * ( dUzdy3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy3 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVzdy3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy3 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness && reduceVStiffness ) {
        answer.at( 3, 8 ) = ( Et / srf ) * At * ( V3 - V2 ) / 2 * ( dUzdy3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy3 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U2 - U3 ) / 2 * ( dVzdy3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy3 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceVStiffness ) {
        answer.at( 3, 8 ) = Et * At * ( V3 - V2 ) / 2 * ( dUzdy3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy3 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U2 - U3 ) / 2 * ( dVzdy3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy3 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness ) {
        answer.at( 3, 8 ) = ( Et / srf ) * At * ( V3 - V2 ) / 2 * ( dUzdy3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy3 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVzdy3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy3 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    } else {
        OOFEM_ERROR( "No stiffness is associated with element (3, 8) of the stiffness matrix of the finite element %d.", this->giveNumber() );
    }
        //dFz1dz3
    //answer.at( 3, 9 ) = Et * At * ( V3 - V2 ) / 2 * ( dUzdz3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz3 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVzdz3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz3 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    if ( srf == 1.0 || ( !reduceUStiffness && !reduceVStiffness ) ) {
        answer.at( 3, 9 ) = Et * At * ( V3 - V2 ) / 2 * ( dUzdz3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz3 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVzdz3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz3 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness && reduceVStiffness ) {
        answer.at( 3, 9 ) = ( Et / srf ) * At * ( V3 - V2 ) / 2 * ( dUzdz3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz3 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U2 - U3 ) / 2 * ( dVzdz3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz3 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceVStiffness ) {
        answer.at( 3, 9 ) = Et * At * ( V3 - V2 ) / 2 * ( dUzdz3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz3 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U2 - U3 ) / 2 * ( dVzdz3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz3 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness ) {
        answer.at( 3, 9 ) = ( Et / srf ) * At * ( V3 - V2 ) / 2 * ( dUzdz3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz3 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U2 - U3 ) / 2 * ( dVzdz3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz3 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    } else {
        OOFEM_ERROR( "No stiffness is associated with element (3, 9) of the stiffness matrix of the finite element %d.", this->giveNumber() );
    }

    //Fourth row in the matrix
        //dFx2dx2
    //answer.at( 4, 4 ) = Et * At * ( V1 - V3 ) / 2 * ( dUxdx2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx2 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U3 - U1 ) / 2 * ( dVxdx2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx2 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    if ( srf == 1.0 || ( !reduceUStiffness && !reduceVStiffness ) ) {
        answer.at( 4, 4 ) = Et * At * ( V1 - V3 ) / 2 * ( dUxdx2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx2 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U3 - U1 ) / 2 * ( dVxdx2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx2 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness && reduceVStiffness ) {
        answer.at( 4, 4 ) = ( Et / srf ) * At * ( V1 - V3 ) / 2 * ( dUxdx2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx2 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U3 - U1 ) / 2 * ( dVxdx2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx2 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceVStiffness ) {
        answer.at( 4, 4 ) = Et * At * ( V1 - V3 ) / 2 * ( dUxdx2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx2 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U3 - U1 ) / 2 * ( dVxdx2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx2 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness ) {
        answer.at( 4, 4 ) = ( Et / srf ) * At * ( V1 - V3 ) / 2 * ( dUxdx2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx2 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U3 - U1 ) / 2 * ( dVxdx2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx2 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else {
        OOFEM_ERROR( "No stiffness is associated with element (4, 4) of the stiffness matrix of the finite element %d.", this->giveNumber() );
    }
        //dFx2dy2
    //answer.at( 4, 5 ) = Et * At * ( V1 - V3 ) / 2 * ( dUxdy2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy2 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U3 - U1 ) / 2 * ( dVxdy2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy2 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    if ( srf == 1.0 || ( !reduceUStiffness && !reduceVStiffness ) ) {
        answer.at( 4, 5 ) = Et * At * ( V1 - V3 ) / 2 * ( dUxdy2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy2 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U3 - U1 ) / 2 * ( dVxdy2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy2 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness && reduceVStiffness ) {
        answer.at( 4, 5 ) = ( Et / srf ) * At * ( V1 - V3 ) / 2 * ( dUxdy2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy2 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U3 - U1 ) / 2 * ( dVxdy2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy2 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceVStiffness ) {
        answer.at( 4, 5 ) = Et * At * ( V1 - V3 ) / 2 * ( dUxdy2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy2 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U3 - U1 ) / 2 * ( dVxdy2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy2 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness ) {
        answer.at( 4, 5 ) = ( Et / srf ) * At * ( V1 - V3 ) / 2 * ( dUxdy2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy2 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U3 - U1 ) / 2 * ( dVxdy2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy2 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else {
        OOFEM_ERROR( "No stiffness is associated with element (4, 5) of the stiffness matrix of the finite element %d.", this->giveNumber() );
    }
        //dFx2dz2
    //answer.at( 4, 6 ) = Et * At * ( V1 - V3 ) / 2 * ( dUxdz2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz2 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U3 - U1 ) / 2 * ( dVxdz2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz2 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    if ( srf == 1.0 || ( !reduceUStiffness && !reduceVStiffness ) ) {
        answer.at( 4, 6 ) = Et * At * ( V1 - V3 ) / 2 * ( dUxdz2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz2 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U3 - U1 ) / 2 * ( dVxdz2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz2 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness && reduceVStiffness ) {
        answer.at( 4, 6 ) = ( Et / srf ) * At * ( V1 - V3 ) / 2 * ( dUxdz2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz2 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U3 - U1 ) / 2 * ( dVxdz2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz2 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceVStiffness ) {
        answer.at( 4, 6 ) = Et * At * ( V1 - V3 ) / 2 * ( dUxdz2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz2 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U3 - U1 ) / 2 * ( dVxdz2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz2 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness ) {
        answer.at( 4, 6 ) = ( Et / srf ) * At * ( V1 - V3 ) / 2 * ( dUxdz2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz2 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U3 - U1 ) / 2 * ( dVxdz2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz2 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else {
        OOFEM_ERROR( "No stiffness is associated with element (4, 6) of the stiffness matrix of the finite element %d.", this->giveNumber() );
    }
        //dFx2dx3
    //answer.at( 4, 7 ) = Et * At * ( V1 - V3 ) / 2 * ( dUxdx3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx3 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U3 - U1 ) / 2 * ( dVxdx3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx3 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    if ( srf == 1.0 || ( !reduceUStiffness && !reduceVStiffness ) ) {
        answer.at( 4, 7 ) = Et * At * ( V1 - V3 ) / 2 * ( dUxdx3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx3 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U3 - U1 ) / 2 * ( dVxdx3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx3 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness && reduceVStiffness ) {
        answer.at( 4, 7 ) = ( Et / srf ) * At * ( V1 - V3 ) / 2 * ( dUxdx3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx3 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U3 - U1 ) / 2 * ( dVxdx3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx3 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceVStiffness ) {
        answer.at( 4, 7 ) = Et * At * ( V1 - V3 ) / 2 * ( dUxdx3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx3 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U3 - U1 ) / 2 * ( dVxdx3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx3 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness ) {
        answer.at( 4, 7 ) = ( Et / srf ) * At * ( V1 - V3 ) / 2 * ( dUxdx3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx3 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U3 - U1 ) / 2 * ( dVxdx3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx3 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else {
        OOFEM_ERROR( "No stiffness is associated with element (4, 7) of the stiffness matrix of the finite element %d.", this->giveNumber() );
    }
        //dFx2dy3
    //answer.at( 4, 8 ) = Et * At * ( V1 - V3 ) / 2 * ( dUxdy3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy3 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U3 - U1 ) / 2 * ( dVxdy3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy3 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    if ( srf == 1.0 || ( !reduceUStiffness && !reduceVStiffness ) ) {
        answer.at( 4, 8 ) = Et * At * ( V1 - V3 ) / 2 * ( dUxdy3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy3 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U3 - U1 ) / 2 * ( dVxdy3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy3 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness && reduceVStiffness ) {
        answer.at( 4, 8 ) = ( Et / srf ) * At * ( V1 - V3 ) / 2 * ( dUxdy3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy3 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U3 - U1 ) / 2 * ( dVxdy3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy3 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceVStiffness ) {
        answer.at( 4, 8 ) = Et * At * ( V1 - V3 ) / 2 * ( dUxdy3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy3 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U3 - U1 ) / 2 * ( dVxdy3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy3 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness ) {
        answer.at( 4, 8 ) = ( Et / srf ) * At * ( V1 - V3 ) / 2 * ( dUxdy3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy3 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U3 - U1 ) / 2 * ( dVxdy3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy3 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else {
        OOFEM_ERROR( "No stiffness is associated with element (4, 8) of the stiffness matrix of the finite element %d.", this->giveNumber() );
    }
        //dFx2dz3
    //answer.at( 4, 9 ) = Et * At * ( V1 - V3 ) / 2 * ( dUxdz3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz3 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U3 - U1 ) / 2 * ( dVxdz3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz3 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    if ( srf == 1.0 || ( !reduceUStiffness && !reduceVStiffness ) ) {
        answer.at( 4, 9 ) = Et * At * ( V1 - V3 ) / 2 * ( dUxdz3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz3 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U3 - U1 ) / 2 * ( dVxdz3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz3 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness && reduceVStiffness ) {
        answer.at( 4, 9 ) = ( Et / srf ) * At * ( V1 - V3 ) / 2 * ( dUxdz3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz3 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U3 - U1 ) / 2 * ( dVxdz3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz3 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceVStiffness ) {
        answer.at( 4, 9 ) = Et * At * ( V1 - V3 ) / 2 * ( dUxdz3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz3 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U3 - U1 ) / 2 * ( dVxdz3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz3 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness ) {
        answer.at( 4, 9 ) = ( Et / srf ) * At * ( V1 - V3 ) / 2 * ( dUxdz3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz3 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U3 - U1 ) / 2 * ( dVxdz3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz3 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else {
        OOFEM_ERROR( "No stiffness is associated with element (4, 9) of the stiffness matrix of the finite element %d.", this->giveNumber() );
    }

    //Fifth row in the matrix
        //dFy2dy2
    //answer.at( 5, 5 ) = Et * At * ( V1 - V3 ) / 2 * ( dUydy2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy2 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U3 - U1 ) / 2 * ( dVydy2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy2 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    if ( srf == 1.0 || ( !reduceUStiffness && !reduceVStiffness ) ) {
        answer.at( 5, 5 ) = Et * At * ( V1 - V3 ) / 2 * ( dUydy2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy2 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U3 - U1 ) / 2 * ( dVydy2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy2 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness && reduceVStiffness ) {
        answer.at( 5, 5 ) = ( Et / srf ) * At * ( V1 - V3 ) / 2 * ( dUydy2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy2 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U3 - U1 ) / 2 * ( dVydy2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy2 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceVStiffness ) {
        answer.at( 5, 5 ) = Et * At * ( V1 - V3 ) / 2 * ( dUydy2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy2 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U3 - U1 ) / 2 * ( dVydy2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy2 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness ) {
        answer.at( 5, 5 ) = ( Et / srf ) * At * ( V1 - V3 ) / 2 * ( dUydy2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy2 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U3 - U1 ) / 2 * ( dVydy2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy2 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else {
        OOFEM_ERROR( "No stiffness is associated with element (5, 5) of the stiffness matrix of the finite element %d.", this->giveNumber() );
    }
        //dFy2dz2
    //answer.at( 5, 6 ) = Et * At * ( V1 - V3 ) / 2 * ( dUydz2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz2 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U3 - U1 ) / 2 * ( dVydz2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz2 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    if ( srf == 1.0 || ( !reduceUStiffness && !reduceVStiffness ) ) {
        answer.at( 5, 6 ) = Et * At * ( V1 - V3 ) / 2 * ( dUydz2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz2 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U3 - U1 ) / 2 * ( dVydz2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz2 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness && reduceVStiffness ) {
        answer.at( 5, 6 ) = ( Et / srf ) * At * ( V1 - V3 ) / 2 * ( dUydz2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz2 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U3 - U1 ) / 2 * ( dVydz2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz2 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceVStiffness ) {
        answer.at( 5, 6 ) = Et * At * ( V1 - V3 ) / 2 * ( dUydz2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz2 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U3 - U1 ) / 2 * ( dVydz2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz2 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness ) {
        answer.at( 5, 6 ) = ( Et / srf ) * At * ( V1 - V3 ) / 2 * ( dUydz2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz2 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U3 - U1 ) / 2 * ( dVydz2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz2 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else {
        OOFEM_ERROR( "No stiffness is associated with element (5, 6) of the stiffness matrix of the finite element %d.", this->giveNumber() );
    }
        //dFy2dx3
    //answer.at( 5, 7 ) = Et * At * ( V1 - V3 ) / 2 * ( dUydx3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx3 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U3 - U1 ) / 2 * ( dVydx3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx3 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    if ( srf == 1.0 || ( !reduceUStiffness && !reduceVStiffness ) ) {
        answer.at( 5, 7 ) = Et * At * ( V1 - V3 ) / 2 * ( dUydx3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx3 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U3 - U1 ) / 2 * ( dVydx3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx3 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness && reduceVStiffness ) {
        answer.at( 5, 7 ) = ( Et / srf ) * At * ( V1 - V3 ) / 2 * ( dUydx3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx3 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U3 - U1 ) / 2 * ( dVydx3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx3 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceVStiffness ) {
        answer.at( 5, 7 ) = Et * At * ( V1 - V3 ) / 2 * ( dUydx3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx3 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U3 - U1 ) / 2 * ( dVydx3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx3 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness ) {
        answer.at( 5, 7 ) = ( Et / srf ) * At * ( V1 - V3 ) / 2 * ( dUydx3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx3 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U3 - U1 ) / 2 * ( dVydx3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx3 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else {
        OOFEM_ERROR( "No stiffness is associated with element (5, 7) of the stiffness matrix of the finite element %d.", this->giveNumber() );
    }
        //dFy2dy3
    //answer.at( 5, 8 ) = Et * At * ( V1 - V3 ) / 2 * ( dUydy3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy3 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U3 - U1 ) / 2 * ( dVydy3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy3 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    if ( srf == 1.0 || ( !reduceUStiffness && !reduceVStiffness ) ) {
        answer.at( 5, 8 ) = Et * At * ( V1 - V3 ) / 2 * ( dUydy3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy3 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U3 - U1 ) / 2 * ( dVydy3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy3 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness && reduceVStiffness ) {
        answer.at( 5, 8 ) = ( Et / srf ) * At * ( V1 - V3 ) / 2 * ( dUydy3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy3 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U3 - U1 ) / 2 * ( dVydy3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy3 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceVStiffness ) {
        answer.at( 5, 8 ) = Et * At * ( V1 - V3 ) / 2 * ( dUydy3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy3 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U3 - U1 ) / 2 * ( dVydy3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy3 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness ) {
        answer.at( 5, 8 ) = ( Et / srf ) * At * ( V1 - V3 ) / 2 * ( dUydy3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy3 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U3 - U1 ) / 2 * ( dVydy3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy3 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else {
        OOFEM_ERROR( "No stiffness is associated with element (5, 8) of the stiffness matrix of the finite element %d.", this->giveNumber() );
    }
        //dFy2dz3
    //answer.at( 5, 9 ) = Et * At * ( V1 - V3 ) / 2 * ( dUydz3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz3 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U3 - U1 ) / 2 * ( dVydz3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz3 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    if ( srf == 1.0 || ( !reduceUStiffness && !reduceVStiffness ) ) {
        answer.at( 5, 9 ) = Et * At * ( V1 - V3 ) / 2 * ( dUydz3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz3 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U3 - U1 ) / 2 * ( dVydz3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz3 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness && reduceVStiffness ) {
        answer.at( 5, 9 ) = ( Et / srf ) * At * ( V1 - V3 ) / 2 * ( dUydz3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz3 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U3 - U1 ) / 2 * ( dVydz3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz3 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceVStiffness ) {
        answer.at( 5, 9 ) = Et * At * ( V1 - V3 ) / 2 * ( dUydz3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz3 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U3 - U1 ) / 2 * ( dVydz3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz3 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness ) {
        answer.at( 5, 9 ) = ( Et / srf ) * At * ( V1 - V3 ) / 2 * ( dUydz3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz3 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U3 - U1 ) / 2 * ( dVydz3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz3 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else {
        OOFEM_ERROR( "No stiffness is associated with element (5, 9) of the stiffness matrix of the finite element %d.", this->giveNumber() );
    }

    //Sixth row in the matrix
        //dFz2dz2
    //answer.at( 6, 6 ) = Et * At * ( V1 - V3 ) / 2 * ( dUzdz2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz2 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U3 - U1 ) / 2 * ( dVzdz2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz2 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    if ( srf == 1.0 || ( !reduceUStiffness && !reduceVStiffness ) ) {
        answer.at( 6, 6 ) = Et * At * ( V1 - V3 ) / 2 * ( dUzdz2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz2 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U3 - U1 ) / 2 * ( dVzdz2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz2 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness && reduceVStiffness ) {
        answer.at( 6, 6 ) = ( Et / srf ) * At * ( V1 - V3 ) / 2 * ( dUzdz2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz2 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U3 - U1 ) / 2 * ( dVzdz2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz2 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceVStiffness ) {
        answer.at( 6, 6 ) = Et * At * ( V1 - V3 ) / 2 * ( dUzdz2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz2 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U3 - U1 ) / 2 * ( dVzdz2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz2 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness ) {
        answer.at( 6, 6 ) = ( Et / srf ) * At * ( V1 - V3 ) / 2 * ( dUzdz2 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz2 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U3 - U1 ) / 2 * ( dVzdz2 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz2 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    } else {
        OOFEM_ERROR( "No stiffness is associated with element (6, 6) of the stiffness matrix of the finite element %d.", this->giveNumber() );
    }
        //dFz2dx3
    //answer.at( 6, 7 ) = Et * At * ( V1 - V3 ) / 2 * ( dUzdx3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx3 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U3 - U1 ) / 2 * ( dVzdx3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx3 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    if ( srf == 1.0 || ( !reduceUStiffness && !reduceVStiffness ) ) {
        answer.at( 6, 7 ) = Et * At * ( V1 - V3 ) / 2 * ( dUzdx3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx3 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U3 - U1 ) / 2 * ( dVzdx3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx3 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness && reduceVStiffness ) {
        answer.at( 6, 7 ) = ( Et / srf ) * At * ( V1 - V3 ) / 2 * ( dUzdx3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx3 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U3 - U1 ) / 2 * ( dVzdx3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx3 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceVStiffness ) {
        answer.at( 6, 7 ) = Et * At * ( V1 - V3 ) / 2 * ( dUzdx3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx3 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U3 - U1 ) / 2 * ( dVzdx3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx3 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness ) {
        answer.at( 6, 7 ) = ( Et / srf ) * At * ( V1 - V3 ) / 2 * ( dUzdx3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx3 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U3 - U1 ) / 2 * ( dVzdx3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx3 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    } else {
        OOFEM_ERROR( "No stiffness is associated with element (6, 7) of the stiffness matrix of the finite element %d.", this->giveNumber() );
    }
        //dFz2dy3
    //answer.at( 6, 8 ) = Et * At * ( V1 - V3 ) / 2 * ( dUzdy3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy3 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U3 - U1 ) / 2 * ( dVzdy3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy3 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    if ( srf == 1.0 || ( !reduceUStiffness && !reduceVStiffness ) ) {
        answer.at( 6, 8 ) = Et * At * ( V1 - V3 ) / 2 * ( dUzdy3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy3 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U3 - U1 ) / 2 * ( dVzdy3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy3 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness && reduceVStiffness ) {
        answer.at( 6, 8 ) = ( Et / srf ) * At * ( V1 - V3 ) / 2 * ( dUzdy3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy3 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U3 - U1 ) / 2 * ( dVzdy3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy3 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceVStiffness ) {
        answer.at( 6, 8 ) = Et * At * ( V1 - V3 ) / 2 * ( dUzdy3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy3 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U3 - U1 ) / 2 * ( dVzdy3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy3 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness ) {
        answer.at( 6, 8 ) = ( Et / srf ) * At * ( V1 - V3 ) / 2 * ( dUzdy3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy3 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U3 - U1 ) / 2 * ( dVzdy3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy3 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    } else {
        OOFEM_ERROR( "No stiffness is associated with element (6, 8) of the stiffness matrix of the finite element %d.", this->giveNumber() );
    }
        //dFz2dz3
    //answer.at( 6, 9 ) = Et * At * ( V1 - V3 ) / 2 * ( dUzdz3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz3 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U3 - U1 ) / 2 * ( dVzdz3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz3 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    if ( srf == 1.0 || ( !reduceUStiffness && !reduceVStiffness ) ) {
        answer.at( 6, 9 ) = Et * At * ( V1 - V3 ) / 2 * ( dUzdz3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz3 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U3 - U1 ) / 2 * ( dVzdz3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz3 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness && reduceVStiffness ) {
        answer.at( 6, 9 ) = ( Et / srf ) * At * ( V1 - V3 ) / 2 * ( dUzdz3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz3 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U3 - U1 ) / 2 * ( dVzdz3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz3 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceVStiffness ) {
        answer.at( 6, 9 ) = Et * At * ( V1 - V3 ) / 2 * ( dUzdz3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz3 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U3 - U1 ) / 2 * ( dVzdz3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz3 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness ) {
        answer.at( 6, 9 ) = ( Et / srf ) * At * ( V1 - V3 ) / 2 * ( dUzdz3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz3 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U3 - U1 ) / 2 * ( dVzdz3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz3 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    } else {
        OOFEM_ERROR( "No stiffness is associated with element (6, 9) of the stiffness matrix of the finite element %d.", this->giveNumber() );
    }

    //Seventh row in the matrix
        //dFx3dx3
    //answer.at( 7, 7 ) = Et * At * ( V2 - V1 ) / 2 * ( dUxdx3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx3 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U1 - U2 ) / 2 * ( dVxdx3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx3 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    if ( srf == 1.0 || ( !reduceUStiffness && !reduceVStiffness ) ) {
        answer.at( 7, 7 ) = Et * At * ( V2 - V1 ) / 2 * ( dUxdx3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx3 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U1 - U2 ) / 2 * ( dVxdx3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx3 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness && reduceVStiffness ) {
        answer.at( 7, 7 ) = ( Et / srf ) * At * ( V2 - V1 ) / 2 * ( dUxdx3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx3 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U1 - U2 ) / 2 * ( dVxdx3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx3 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceVStiffness ) {
        answer.at( 7, 7 ) = Et * At * ( V2 - V1 ) / 2 * ( dUxdx3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx3 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U1 - U2 ) / 2 * ( dVxdx3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx3 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness ) {
        answer.at( 7, 7 ) = ( Et / srf ) * At * ( V2 - V1 ) / 2 * ( dUxdx3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdx3 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U1 - U2 ) / 2 * ( dVxdx3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdx3 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else {
        OOFEM_ERROR( "No stiffness is associated with element (7, 7) of the stiffness matrix of the finite element %d.", this->giveNumber() );
    }
        //dFx3dy3
    //answer.at( 7, 8 ) = Et * At * ( V2 - V1 ) / 2 * ( dUxdy3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy3 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U1 - U2 ) / 2 * ( dVxdy3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy3 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    if ( srf == 1.0 || ( !reduceUStiffness && !reduceVStiffness ) ) {
        answer.at( 7, 8 ) = Et * At * ( V2 - V1 ) / 2 * ( dUxdy3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy3 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U1 - U2 ) / 2 * ( dVxdy3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy3 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness && reduceVStiffness ) {
        answer.at( 7, 8 ) = ( Et / srf ) * At * ( V2 - V1 ) / 2 * ( dUxdy3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy3 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U1 - U2 ) / 2 * ( dVxdy3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy3 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceVStiffness ) {
        answer.at( 7, 8 ) = Et * At * ( V2 - V1 ) / 2 * ( dUxdy3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy3 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U1 - U2 ) / 2 * ( dVxdy3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy3 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness ) {
        answer.at( 7, 8 ) = ( Et / srf ) * At * ( V2 - V1 ) / 2 * ( dUxdy3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy3 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U1 - U2 ) / 2 * ( dVxdy3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy3 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else {
        OOFEM_ERROR( "No stiffness is associated with element (7, 8) of the stiffness matrix of the finite element %d.", this->giveNumber() );
    }
        //dFx3dz3
    //answer.at( 7, 9 ) = Et * At * ( V2 - V1 ) / 2 * ( dUxdz3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz3 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U1 - U2 ) / 2 * ( dVxdz3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz3 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    if ( srf == 1.0 || ( !reduceUStiffness && !reduceVStiffness ) ) {
        answer.at( 7, 9 ) = Et * At * ( V2 - V1 ) / 2 * ( dUxdz3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz3 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U1 - U2 ) / 2 * ( dVxdz3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz3 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness && reduceVStiffness ) {
        answer.at( 7, 9 ) = ( Et / srf ) * At * ( V2 - V1 ) / 2 * ( dUxdz3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz3 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U1 - U2 ) / 2 * ( dVxdz3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz3 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceVStiffness ) {
        answer.at( 7, 9 ) = Et * At * ( V2 - V1 ) / 2 * ( dUxdz3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz3 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U1 - U2 ) / 2 * ( dVxdz3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz3 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness ) {
        answer.at( 7, 9 ) = ( Et / srf ) * At * ( V2 - V1 ) / 2 * ( dUxdz3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz3 * U.at( 1 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U1 - U2 ) / 2 * ( dVxdz3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz3 * V.at( 1 ) / pow( V.computeNorm(), 2 ) );
    } else {
        OOFEM_ERROR( "No stiffness is associated with element (7, 9) of the stiffness matrix of the finite element %d.", this->giveNumber() );
    }

    //Eigth row in the matrix
        //dFy3dy3
    //answer.at( 8, 8 ) = Et * At * ( V2 - V1 ) / 2 * ( dUydy3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy3 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U1 - U2 ) / 2 * ( dVydy3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy3 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    if ( srf == 1.0 || ( !reduceUStiffness && !reduceVStiffness ) ) {
        answer.at( 8, 8 ) = Et * At * ( V2 - V1 ) / 2 * ( dUydy3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy3 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U1 - U2 ) / 2 * ( dVydy3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy3 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness && reduceVStiffness ) {
        answer.at( 8, 8 ) = ( Et / srf ) * At * ( V2 - V1 ) / 2 * ( dUydy3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy3 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U1 - U2 ) / 2 * ( dVydy3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy3 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceVStiffness ) {
        answer.at( 8, 8 ) = Et * At * ( V2 - V1 ) / 2 * ( dUydy3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy3 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U1 - U2 ) / 2 * ( dVydy3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy3 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness ) {
        answer.at( 8, 8 ) = ( Et / srf ) * At * ( V2 - V1 ) / 2 * ( dUydy3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdy3 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U1 - U2 ) / 2 * ( dVydy3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdy3 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else {
        OOFEM_ERROR( "No stiffness is associated with element (8, 8) of the stiffness matrix of the finite element %d.", this->giveNumber() );
    }
        //dFy3dz3
    //answer.at( 8, 9 ) = Et * At * ( V2 - V1 ) / 2 * ( dUydz3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz3 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U1 - U2 ) / 2 * ( dVydz3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz3 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    if ( srf == 1.0 || ( !reduceUStiffness && !reduceVStiffness ) ) {
        answer.at( 8, 9 ) = Et * At * ( V2 - V1 ) / 2 * ( dUydz3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz3 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U1 - U2 ) / 2 * ( dVydz3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz3 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness && reduceVStiffness ) {
        answer.at( 8, 9 ) = ( Et / srf ) * At * ( V2 - V1 ) / 2 * ( dUydz3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz3 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U1 - U2 ) / 2 * ( dVydz3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz3 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceVStiffness ) {
        answer.at( 8, 9 ) = Et * At * ( V2 - V1 ) / 2 * ( dUydz3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz3 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U1 - U2 ) / 2 * ( dVydz3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz3 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness ) {
        answer.at( 8, 9 ) = ( Et / srf ) * At * ( V2 - V1 ) / 2 * ( dUydz3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz3 * U.at( 2 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U1 - U2 ) / 2 * ( dVydz3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz3 * V.at( 2 ) / pow( V.computeNorm(), 2 ) );
    } else {
        OOFEM_ERROR( "No stiffness is associated with element (8, 9) of the stiffness matrix of the finite element %d.", this->giveNumber() );
    }

    //Ninth row in the matrix
        //dFz3dz3
    //answer.at( 9, 9 ) = Et * At * ( V2 - V1 ) / 2 * ( dUzdz3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz3 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U1 - U2 ) / 2 * ( dVzdz3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz3 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    if ( srf == 1.0 || ( !reduceUStiffness && !reduceVStiffness ) ) {
        answer.at( 9, 9 ) = Et * At * ( V2 - V1 ) / 2 * ( dUzdz3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz3 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U1 - U2 ) / 2 * ( dVzdz3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz3 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness && reduceVStiffness ) {
        answer.at( 9, 9 ) = ( Et / srf ) * At * ( V2 - V1 ) / 2 * ( dUzdz3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz3 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U1 - U2 ) / 2 * ( dVzdz3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz3 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceVStiffness ) {
        answer.at( 9, 9 ) = Et * At * ( V2 - V1 ) / 2 * ( dUzdz3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz3 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + ( Et / srf ) * At * ( U1 - U2 ) / 2 * ( dVzdz3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz3 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    } else if ( reduceUStiffness ) {
        answer.at( 9, 9 ) = ( Et / srf ) * At * ( V2 - V1 ) / 2 * ( dUzdz3 * ( 1 / L0 - 1 / U.computeNorm() ) + dUdz3 * U.at( 3 ) / pow( U.computeNorm(), 2 ) ) + Et * At * ( U1 - U2 ) / 2 * ( dVzdz3 * ( 1 / L0 - 1 / V.computeNorm() ) + dVdz3 * V.at( 3 ) / pow( V.computeNorm(), 2 ) );
    } else {
        OOFEM_ERROR( "No stiffness is associated with element (9, 9) of the stiffness matrix of the finite element %d.", this->giveNumber() );
    }
    
    answer.symmetrized();
    for ( int i = 1; i <= 9; i++ ) {
        for ( int j = 1; j <= 9; j++ )
            answer.at( i, j ) = -1 * answer.at( i, j );
    }
}

FloatArray
NetTr3Pr::computeUTwine( TimeStep *tStep )
{
    // Fetch coordinates of element's nodes in initial configuration
    FloatArray node1 = this->giveNode( 1 )->giveCoordinates();
    FloatArray node2 = this->giveNode( 2 )->giveCoordinates();
    FloatArray node3 = this->giveNode( 3 )->giveCoordinates();
    // Fetch total displacements in the current configuration
    FloatArray u;
    u.resize( 9 );
    if ( !tStep->isTheFirstStep() )
        this->computeVectorOf( VM_Total, tStep, u );
    // Current coordinates of node 1
    double x1 = node1.at( 1 ) + u.at( 1 );
    double y1 = node1.at( 2 ) + u.at( 2 );
    double z1 = node1.at( 3 ) + u.at( 3 );
    // Current coordinates of node 2
    double x2 = node2.at( 1 ) + u.at( 4 );
    double y2 = node2.at( 2 ) + u.at( 5 );
    double z2 = node2.at( 3 ) + u.at( 6 );
    // Current coordinates of node 3
    double x3 = node3.at( 1 ) + u.at( 7 );
    double y3 = node3.at( 2 ) + u.at( 8 );
    double z3 = node3.at( 3 ) + u.at( 9 );

    // Express unit vectors of the twine coord. system in the global coord. syst.
    FloatArray U;
    U.resize( 3 );
    U.at( 1 ) = ( V3 - V1 ) * ( x2 - x1 ) / d - ( V2 - V1 ) * ( x3 - x1 ) / d;
    U.at( 2 ) = ( V3 - V1 ) * ( y2 - y1 ) / d - ( V2 - V1 ) * ( y3 - y1 ) / d;
    U.at( 3 ) = ( V3 - V1 ) * ( z2 - z1 ) / d - ( V2 - V1 ) * ( z3 - z1 ) / d;
    return U;
}

FloatArray
NetTr3Pr::computeVTwine( TimeStep *tStep )
{
    // Fetch initial coordinates of element's nodes
    FloatArray node1 = this->giveNode( 1 )->giveCoordinates();
    FloatArray node2 = this->giveNode( 2 )->giveCoordinates();
    FloatArray node3 = this->giveNode( 3 )->giveCoordinates();
    // Fetch total displacements in the current configuration
    FloatArray u;
    u.resize( 9 );
    if ( !tStep->isTheFirstStep() )
        this->computeVectorOf( VM_Total, tStep, u );
    // Current coordinates of node 1
    double x1 = node1.at( 1 ) + u.at( 1 );
    double y1 = node1.at( 2 ) + u.at( 2 );
    double z1 = node1.at( 3 ) + u.at( 3 );
    // Current coordinates of node 2
    double x2 = node2.at( 1 ) + u.at( 4 );
    double y2 = node2.at( 2 ) + u.at( 5 );
    double z2 = node2.at( 3 ) + u.at( 6 );
    // Current coordinates of node 3
    double x3 = node3.at( 1 ) + u.at( 7 );
    double y3 = node3.at( 2 ) + u.at( 8 );
    double z3 = node3.at( 3 ) + u.at( 9 );

    // Express unit vectors of the twine coord. system in the global coord. syst.
    FloatArray V;
    V.resize( 3 );
    V.at( 1 ) = ( U2 - U1 ) * ( x3 - x1 ) / d - ( U3 - U1 ) * ( x2 - x1 ) / d;
    V.at( 2 ) = ( U2 - U1 ) * ( y3 - y1 ) / d - ( U3 - U1 ) * ( y2 - y1 ) / d;
    V.at( 3 ) = ( U2 - U1 ) * ( z3 - z1 ) / d - ( U3 - U1 ) * ( z2 - z1 ) / d;
    return V;
}

std::unique_ptr<IntegrationRule>
NetTr3Pr::giveBoundarySurfaceIntegrationRule( int order, int boundary )
{
    return this->giveInterpolation()->giveIntegrationRule( order );
}

void
NetTr3Pr::giveDofManDofIDMask( int inode, IntArray &answer ) const
{
    answer = { D_u, D_v, D_w };
}

void
NetTr3Pr::giveInternalForcesVector( FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord )
{
    FloatArray strains, stresses;
    double At = 0;
    for ( auto gp : *this->giveDefaultIntegrationRulePtr() ) {
        // Calculate strains and stresses in U- and V-twines
        computeStrainVector( strains, gp, tStep );
        computeStressVector( stresses, strains, gp, tStep );
        At = this->giveCrossSection()->give( CS_Area, gp );
    }

    // Calculate tension forces in U- and V-twines
    double Tu = At * stresses.at( 1 );
    double Tv = At * stresses.at( 2 );

    // Calculate current length of a twine in U- and in V-direction
    FloatArray U = computeUTwine( tStep );
    FloatArray V = computeVTwine( tStep );

    // Calculate internal nodal forces    
    answer.resize( 9 );
        // Node 1
    answer.at( 1 ) = ( V3 - V2 ) * Tu * U.at( 1 ) / ( 2 * U.computeNorm() ) + ( U2 - U3 ) * Tv * V.at( 1 ) / ( 2 * V.computeNorm() );
    answer.at( 2 ) = ( V3 - V2 ) * Tu * U.at( 2 ) / ( 2 * U.computeNorm() ) + ( U2 - U3 ) * Tv * V.at( 2 ) / ( 2 * V.computeNorm() );
    answer.at( 3 ) = ( V3 - V2 ) * Tu * U.at( 3 ) / ( 2 * U.computeNorm() ) + ( U2 - U3 ) * Tv * V.at( 3 ) / ( 2 * V.computeNorm() );
        // Node 2
    answer.at( 4 ) = ( V1 - V3 ) * Tu * U.at( 1 ) / ( 2 * U.computeNorm() ) + ( U3 - U1 ) * Tv * V.at( 1 ) / ( 2 * V.computeNorm() );     
    answer.at( 5 ) = ( V1 - V3 ) * Tu * U.at( 2 ) / ( 2 * U.computeNorm() ) + ( U3 - U1 ) * Tv * V.at( 2 ) / ( 2 * V.computeNorm() );  
    answer.at( 6 ) = ( V1 - V3 ) * Tu * U.at( 3 ) / ( 2 * U.computeNorm() ) + ( U3 - U1 ) * Tv * V.at( 3 ) / ( 2 * V.computeNorm() );
        // Node 3   
    answer.at( 7 ) = ( V2 - V1 ) * Tu * U.at( 1 ) / ( 2 * U.computeNorm() ) + ( U1 - U2 ) * Tv * V.at( 1 ) / ( 2 * V.computeNorm() );
    answer.at( 8 ) = ( V2 - V1 ) * Tu * U.at( 2 ) / ( 2 * U.computeNorm() ) + ( U1 - U2 ) * Tv * V.at( 2 ) / ( 2 * V.computeNorm() ); 
    answer.at( 9 ) = ( V2 - V1 ) * Tu * U.at( 3 ) / ( 2 * U.computeNorm() ) + ( U1 - U2 ) * Tv * V.at( 3 ) / ( 2 * V.computeNorm() );

    for ( int i = 1; i <= 9; i++ )
        answer.at( i ) = -1 * answer.at( i );
}

FEInterpolation *
NetTr3Pr::giveInterpolation() const { return &interp; }

void
NetTr3Pr ::computeStrainVector( FloatArray &answer, GaussPoint *gp, TimeStep *tStep )
{
    answer.resize( 2 );
    // Calculate current dimension of a U- and a V- twine
    double lengthU, lengthV;
    lengthU = computeUTwine(tStep).computeNorm();
    lengthV = computeVTwine(tStep).computeNorm();

    if ( srf == 1 || ( lengthU >= ( 0.995 * L0 ) && lengthV >= ( 0.995 * L0 ) ) ) {
        // Calculate strains in U- and V-twines
        answer.at( 1 ) = ( lengthU - L0 ) / L0;
        answer.at( 2 ) = ( lengthV - L0 ) / L0;
    } else if ( lengthU < ( 0.995 * L0 ) && lengthV >= ( 0.995 * L0 ) ) {
        // Calculate strains in V-twines. Reduce strains in U-twines since these are in compression
        answer.at( 1 ) = ( lengthU - L0 ) / ( srf * L0 );
        answer.at( 2 ) = ( lengthV - L0 ) / L0;
    } else if ( lengthV < ( 0.995 * L0 ) && lengthU >= ( 0.995 * L0 ) ) {
        // Calculate strains in U-twines. Reduce strains in V-twines since these are in compression
        answer.at( 1 ) = ( lengthU - L0 ) / L0;
        answer.at( 2 ) = ( lengthV - L0 ) / ( srf * L0 );
    } else {
        // Reduce strains in U- and V-twines since both are in compression
        answer.at( 1 ) = ( lengthU - L0 ) / ( srf * L0 );
        answer.at( 2 ) = ( lengthV - L0 ) / ( srf * L0 );
    }
}

void
NetTr3Pr ::computeStressVector( FloatArray &answer, const FloatArray &strains, GaussPoint *gp, TimeStep *tStep )
{
    answer = this->giveStructuralCrossSection()->giveRealStress_Netting( strains, gp, tStep );
}

void
NetTr3Pr::calculateEquivalentLumpedNodalValues( FloatArray &answer, FloatArray vector )
{
    answer.resize( 9 );
    answer.at( 1 ) = answer.at( 4 ) = answer.at( 7 ) = vector.at( 1 ) / 3;
    answer.at( 2 ) = answer.at( 5 ) = answer.at( 8 ) = vector.at( 2 ) / 3;
    answer.at( 3 ) = answer.at( 6 ) = answer.at( 9 ) = vector.at( 3 ) / 3;
}

void
NetTr3Pr::initializeFrom( InputRecord &ir )
{
    StructuralElement::initializeFrom( ir );

    IR_GIVE_FIELD( ir, L0, _IFT_NetTr3Pr_L0 );
	IR_GIVE_FIELD( ir, U1, _IFT_NetTr3Pr_U1 );
	IR_GIVE_FIELD( ir, V1, _IFT_NetTr3Pr_V1 );
	IR_GIVE_FIELD( ir, U2, _IFT_NetTr3Pr_U2 );
	IR_GIVE_FIELD( ir, V2, _IFT_NetTr3Pr_V2 );
	IR_GIVE_FIELD( ir, U3, _IFT_NetTr3Pr_U3 );
	IR_GIVE_FIELD( ir, V3, _IFT_NetTr3Pr_V3 );
    IR_GIVE_OPTIONAL_FIELD( ir, srf, _IFT_NetTr3Pr_srf );

    d = ( U2 - U1 ) * ( V3 - V1 ) - ( U3 - U1 ) * ( V2 - V1 );
}
}