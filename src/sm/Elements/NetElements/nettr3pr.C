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

    // Calculate elements of the stiffness matrix. This stiffness matrix is already in the global coord. system
    answer.resize( 9, 9 );

    //First row in the matrix
        //dFx1dx1
    answer.at( 1, 1 ) = Et * At * ( V3 - V2 ) / 2 * ( ( V2 - V3 ) / d * ( 1 / L0 - 1 / U.computeNorm() ) + U.at( 1 ) / pow( U.computeNorm(), 2 ) * ( V2 - V3 ) / pow( d, 2 ) * ( ( node2.at( 1 ) - node1.at( 1 ) ) * ( V3 - V1 ) - ( node3.at( 1 ) - node1.at( 1 ) ) * ( V2 - V1 ) ) )
        + Et * At * ( U2 - U3 ) / 2 * ( ( U3 - U2 ) / d * ( 1 / L0 - 1 / V.computeNorm() ) + V.at( 1 ) / pow( V.computeNorm(), 2 ) * ( U2 - U3 ) / pow( d, 2 ) * ( ( node2.at( 1 ) - node1.at( 1 ) ) * ( U3 - U1 ) - ( node3.at( 1 ) - node1.at( 1 ) ) * ( U2 - U1 ) ) );
        //dFx1dy1
    answer.at( 1, 2 ) = Et * At * ( V3 - V2 ) / 2 * ( U.at( 1 ) * ( V2 - V3 ) / ( pow( U.computeNorm(), 2 ) * pow( d, 2 ) ) * ( ( node2.at( 2 ) - node1.at( 2 ) ) * ( V3 - V1 ) - ( node3.at( 2 ) - node1.at( 2 ) ) * ( V2 - V1 ) ) )
        + Et * At * ( U2 - U3 ) / 2 * ( V.at( 1 ) * ( U2 - U3 ) / ( pow( V.computeNorm(), 2 ) * pow( d, 2 ) ) * ( ( node2.at( 2 ) - node1.at( 2 ) ) * ( U3 - U1 ) - ( node3.at( 2 ) - node1.at( 2 ) ) * ( U2 - U1 ) ) );
        //dFx1dz1
    answer.at( 1, 3 ) = Et * At * ( V3 - V2 ) / 2 * ( U.at( 1 ) * ( V2 - V3 ) / ( pow( U.computeNorm(), 2 ) * pow( d, 2 ) ) * ( ( z2 - z1 ) * ( V3 - V1 ) - ( z3 - z1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U2 - U3 ) / 2 * ( V.at( 1 ) * ( U2 - U3 ) / ( pow( V.computeNorm(), 2 ) * pow( d, 2 ) ) * ( ( z2 - z1 ) * ( U3 - U1 ) - ( z3 - z1 ) * ( U2 - U1 ) ) );
        //dFx1dx2
    answer.at( 1, 4 ) = Et * At * ( V3 - V2 ) / 2 * ( ( V3 - V1 ) / d * ( 1 / L0 - 1 / U.computeNorm() ) + U.at( 1 ) / pow( U.computeNorm(), 2 ) * ( V3 - V1 ) / pow( U.computeNorm(), 2 ) * ( ( x2 - x1 ) * ( V3 - V1 ) - ( x3 - x1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U2 - U3 ) / 2 * ( ( U1 - U3 ) / d * ( 1 / L0 - 1 / V.computeNorm() ) + V.at( 1 ) / pow( V.computeNorm(), 2 ) * ( U3 - U1 ) / pow( U.computeNorm(), 2 ) * ( ( x2 - x1 ) * ( U3 - U1 ) - ( x3 - x1 ) * ( U2 - U1 ) ) );
        //dFx1dy2
    answer.at( 1, 5 ) = Et * At * ( V3 - V2 ) / 2 * ( U.at( 1 ) * ( V3 - V1 ) / ( pow( U.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( y2 - y1 ) * ( V3 - V1 ) - ( y3 - y1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U2 - U3 ) / 2 * ( V.at( 1 ) * ( U3 - U1 ) / ( pow( V.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( y2 - y1 ) * ( U3 - U1 ) - ( y3 - y1 ) * ( U2 - U1 ) ) );
        //dFx1dz2
    answer.at( 1, 6 ) = Et * At * ( V3 - V2 ) / 2 * ( U.at( 1 ) * ( V3 - V1 ) / ( pow( U.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( z2 - z1 ) * ( V3 - V1 ) - ( z3 - z1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U2 - U3 ) / 2 * ( V.at( 1 ) * ( U3 - U1 ) / ( pow( V.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( z2 - z1 ) * ( U3 - U1 ) - ( z3 - z1 ) * ( U2 - U1 ) ) );
        //dFx1dx3
    answer.at( 1, 7 ) = Et * At * ( V3 - V2 ) / 2 * ( ( V1 - V2 ) / d * ( 1 / L0 - 1 / U.computeNorm() ) + U.at( 1 ) / pow( U.computeNorm(), 2 ) * ( V1 - V2 ) / pow( U.computeNorm(), 2 ) * ( ( x2 - x1 ) * ( V3 - V1 ) - ( x3 - x1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U2 - U3 ) / 2 * ( ( U2 - U1 ) / d * ( 1 / L0 - 1 / V.computeNorm() ) + V.at( 1 ) / pow( V.computeNorm(), 2 ) * ( U1 - U2 ) / pow( U.computeNorm(), 2 ) * ( ( x2 - x1 ) * ( U3 - U1 ) - ( x3 - x1 ) * ( U2 - U1 ) ) );
        //dFx1dy3
    answer.at( 1, 8 ) = Et * At * ( V3 - V2 ) / 2 * ( U.at( 1 ) * ( V1 - V2 ) / ( pow( U.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( y2 - y1 ) * ( V3 - V1 ) - ( y3 - y1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U2 - U3 ) / 2 * ( V.at( 1 ) * ( U1 - U2 ) / ( pow( V.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( y2 - y1 ) * ( U3 - U1 ) - ( y3 - y1 ) * ( U2 - U1 ) ) );
        //dFx1dz3
    answer.at( 1, 9 ) = Et * At * ( V3 - V2 ) / 2 * ( U.at( 1 ) * ( V1 - V2 ) / ( pow( U.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( z2 - z1 ) * ( V3 - V1 ) - ( z3 - z1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U2 - U3 ) / 2 * ( V.at( 1 ) * ( U1 - U2 ) / ( pow( V.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( z2 - z1 ) * ( U3 - U1 ) - ( z3 - z1 ) * ( U2 - U1 ) ) );

    //Second row in the matrix
        //dFy1dx1
    answer.at( 2, 1 ) = Et * At * ( V3 - V2 ) / 2 * ( U.at(2) * ( V2 - V3 ) / ( pow( U.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( x2 - x1 ) * ( V3 - V1 ) - ( x3 - x1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U2 - U3 ) / 2 * ( V.at(2) * ( U2 - U3 ) / ( pow( V.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( x2 - x1 ) * ( U3 - U1 ) - ( x3 - x1 ) * ( U2 - U1 ) ) );
        //dFy1dy1
    answer.at( 2, 2 ) = Et * At * ( V3 - V2 ) / 2 * ( ( V2 - V3 ) / d * ( 1 / L0 - 1 / U.computeNorm() ) + U.at( 2 ) / pow( U.computeNorm(), 2 ) * ( V2 - V3 ) / pow( U.computeNorm(), 2 ) * ( ( y2 - y1 ) * ( V3 - V1 ) - ( y3 - y1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U2 - U3 ) / 2 * ( ( U3 - U2 ) / d * ( 1 / L0 - 1 / V.computeNorm() ) + V.at( 2 ) / pow( V.computeNorm(), 2 ) * ( U2 - U3 ) / pow( U.computeNorm(), 2 ) * ( ( y2 - y1 ) * ( U3 - U1 ) - ( y3 - y1 ) * ( U2 - U1 ) ) );
        //dFy1dz1
    answer.at( 2, 3 ) = Et * At * ( V3 - V2 ) / 2 * ( U.at( 2 ) * ( V2 - V3 ) / ( pow( U.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( z2 - z1 ) * ( V3 - V1 ) - ( z3 - z1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U2 - U3 ) / 2 * ( V.at( 2 ) * ( U2 - U3 ) / ( pow( V.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( z2 - z1 ) * ( U3 - U1 ) - ( z3 - z1 ) * ( U2 - U1 ) ) );
        //dFy1dx2
    answer.at( 2, 4 ) = Et * At * ( V3 - V2 ) / 2 * ( U.at( 2 ) * ( V3 - V1 ) / ( pow( U.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( x2 - x1 ) * ( V3 - V1 ) - ( x3 - x1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U2 - U3 ) / 2 * ( V.at( 2 ) * ( U3 - U1 ) / ( pow( V.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( x2 - x1 ) * ( U3 - U1 ) - ( x3 - x1 ) * ( U2 - U1 ) ) );
        //dFy1dy2
    answer.at( 2, 5 ) = Et * At * ( V3 - V2 ) / 2 * ( ( V3 - V1 ) / d * ( 1 / L0 - 1 / U.computeNorm() ) + U.at( 2 ) / pow( U.computeNorm(), 2 ) * ( V3 - V1 ) / pow( U.computeNorm(), 2 ) * ( ( y2 - y1 ) * ( V3 - V1 ) - ( y3 - y1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U2 - U3 ) / 2 * ( ( U1 - U3 ) / d * ( 1 / L0 - 1 / V.computeNorm() ) + V.at( 2 ) / pow( V.computeNorm(), 2 ) * ( U3 - U1 ) / pow( U.computeNorm(), 2 ) * ( ( y2 - y1 ) * ( U3 - U1 ) - ( y3 - y1 ) * ( U2 - U1 ) ) );
        //dFy1dz2
    answer.at( 2, 6 ) = Et * At * ( V3 - V2 ) / 2 * ( U.at( 2 ) * ( V3 - V1 ) / ( pow( U.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( z2 - z1 ) * ( V3 - V1 ) - ( z3 - z1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U2 - U3 ) / 2 * ( V.at( 2 ) * ( U3 - U1 ) / ( pow( V.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( z2 - z1 ) * ( U3 - U1 ) - ( z3 - z1 ) * ( U2 - U1 ) ) );
        //dFy1dx3
    answer.at( 2, 7 ) = Et * At * ( V3 - V2 ) / 2 * ( U.at( 2 ) * ( V1 - V2 ) / ( pow( U.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( x2 - x1 ) * ( V3 - V1 ) - ( x3 - x1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U2 - U3 ) / 2 * ( V.at( 2 ) * ( U1 - U2 ) / ( pow( V.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( x2 - x1 ) * ( U3 - U1 ) - ( x3 - x1 ) * ( U2 - U1 ) ) );
        //dFy1dy3
    answer.at( 2, 8 ) = Et * At * ( V3 - V2 ) / 2 * ( ( V1 - V2 ) / d * ( 1 / L0 - 1 / U.computeNorm() ) + U.at( 2 ) / pow( U.computeNorm(), 2 ) * ( V1 - V2 ) / pow( U.computeNorm(), 2 ) * ( ( y2 - y1 ) * ( V3 - V1 ) - ( y3 - y1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U2 - U3 ) / 2 * ( ( U2 - U1 ) / d * ( 1 / L0 - 1 / V.computeNorm() ) + V.at( 2 ) / pow( V.computeNorm(), 2 ) * ( U1 - U2 ) / pow( U.computeNorm(), 2 ) * ( ( y2 - y1 ) * ( U3 - U1 ) - ( y3 - y1 ) * ( U2 - U1 ) ) );
        //dFy1dz3
    answer.at( 2, 9 ) = Et * At * ( V3 - V2 ) / 2 * ( U.at( 2 ) * ( V1 - V2 ) / ( pow( U.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( z2 - z1 ) * ( V3 - V1 ) - ( z3 - z1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U2 - U3 ) / 2 * ( V.at( 2 ) * ( U1 - U2 ) / ( pow( V.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( z2 - z1 ) * ( U3 - U1 ) - ( z3 - z1 ) * ( U2 - U1 ) ) );

    //Third row in the matrix
        //dFz1dx1
    answer.at( 3, 1 ) = Et * At * ( V3 - V2 ) / 2 * ( U.at( 3 ) * ( V2 - V3 ) / ( pow( U.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( x2 - x1 ) * ( V3 - V1 ) - ( x3 - x1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U2 - U3 ) / 2 * ( V.at( 3 ) * ( U2 - U3 ) / ( pow( V.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( x2 - x1 ) * ( U3 - U1 ) - ( x3 - x1 ) * ( U2 - U1 ) ) );
        //dFz1dy1
    answer.at( 3, 2 ) = Et * At * ( V3 - V2 ) / 2 * ( U.at( 3 ) * ( V2 - V3 ) / ( pow( U.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( y2 - y1 ) * ( V3 - V1 ) - ( y3 - y1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U2 - U3 ) / 2 * ( V.at( 3 ) * ( U2 - U3 ) / ( pow( V.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( y2 - y1 ) * ( U3 - U1 ) - ( y3 - y1 ) * ( U2 - U1 ) ) );
        //dFz1dz1
    answer.at( 3, 3 ) = Et * At * ( V3 - V2 ) / 2 * ( ( V2 - V3 ) / d * ( 1 / L0 - 1 / U.computeNorm() ) + U.at( 3 ) / pow( U.computeNorm(), 2 ) * ( V2 - V3 ) / pow( U.computeNorm(), 2 ) * ( ( z2 - z1 ) * ( V3 - V1 ) - ( z3 - z1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U2 - U3 ) / 2 * ( ( U3 - U2 ) / d * ( 1 / L0 - 1 / V.computeNorm() ) + V.at( 3 ) / pow( V.computeNorm(), 2 ) * ( U2 - U3 ) / pow( U.computeNorm(), 2 ) * ( ( z2 - z1 ) * ( U3 - U1 ) - ( z3 - z1 ) * ( U2 - U1 ) ) );
        //dFz1dx2
    answer.at( 3, 4 ) = Et * At * ( V3 - V2 ) / 2 * ( U.at( 3 ) * ( V3 - V1 ) / ( pow( U.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( x2 - x1 ) * ( V3 - V1 ) - ( x3 - x1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U2 - U3 ) / 2 * ( V.at( 3 ) * ( U3 - U1 ) / ( pow( V.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( x2 - x1 ) * ( U3 - U1 ) - ( x3 - x1 ) * ( U2 - U1 ) ) );
        //dFz1dy2
    answer.at( 3, 5 ) = Et * At * ( V3 - V2 ) / 2 * ( U.at( 3 ) * ( V3 - V1 ) / ( pow( U.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( y2 - y1 ) * ( V3 - V1 ) - ( y3 - y1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U2 - U3 ) / 2 * ( V.at( 3 ) * ( U3 - U1 ) / ( pow( V.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( y2 - y1 ) * ( U3 - U1 ) - ( y3 - y1 ) * ( U2 - U1 ) ) );
        //dFz1dz2
    answer.at( 3, 6 ) = Et * At * ( V3 - V2 ) / 2 * ( ( V3 - V1 ) / d * ( 1 / L0 - 1 / U.computeNorm() ) + U.at( 3 ) / pow( U.computeNorm(), 2 ) * ( V3 - V1 ) / pow( U.computeNorm(), 2 ) * ( ( z2 - z1 ) * ( V3 - V1 ) - ( z3 - z1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U2 - U3 ) / 2 * ( ( U1 - U3 ) / d * ( 1 / L0 - 1 / V.computeNorm() ) + V.at( 3 ) / pow( V.computeNorm(), 2 ) * ( U3 - U1 ) / pow( U.computeNorm(), 2 ) * ( ( z2 - z1 ) * ( U3 - U1 ) - ( z3 - z1 ) * ( U2 - U1 ) ) );
        //dFz1dx3
    answer.at( 3, 7 ) = Et * At * ( V3 - V2 ) / 2 * ( U.at( 3 ) * ( V1 - V2 ) / ( pow( U.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( x2 - x1 ) * ( V3 - V1 ) - ( x3 - x1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U2 - U3 ) / 2 * ( V.at( 3 ) * ( U1 - U2 ) / ( pow( V.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( x2 - x1 ) * ( U3 - U1 ) - ( x3 - x1 ) * ( U2 - U1 ) ) );
        //dFz1dy3
    answer.at( 3, 8 ) = Et * At * ( V3 - V2 ) / 2 * ( U.at( 3 ) * ( V1 - V2 ) / ( pow( U.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( y2 - y1 ) * ( V3 - V1 ) - ( y3 - y1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U2 - U3 ) / 2 * ( V.at( 3 ) * ( U1 - U2 ) / ( pow( V.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( y2 - y1 ) * ( U3 - U1 ) - ( y3 - y1 ) * ( U2 - U1 ) ) );

        //dFz1dz3
    answer.at( 3, 9 ) = Et * At * ( V3 - V2 ) / 2 * ( ( V1 - V2 ) / d * ( 1 / L0 - 1 / U.computeNorm() ) + U.at( 3 ) / pow( U.computeNorm(), 2 ) * ( V1 - V2 ) / pow( U.computeNorm(), 2 ) * ( ( z2 - z1 ) * ( V3 - V1 ) - ( z3 - z1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U2 - U3 ) / 2 * ( ( U2 - U1 ) / d * ( 1 / L0 - 1 / V.computeNorm() ) + V.at( 3 ) / pow( V.computeNorm(), 2 ) * ( U1 - U2 ) / pow( U.computeNorm(), 2 ) * ( ( z2 - z1 ) * ( U3 - U1 ) - ( z3 - z1 ) * ( U2 - U1 ) ) );

    //Fourth row in the matrix
        //dFx2dx1
    answer.at( 4, 1 ) = Et * At * ( V1 - V3 ) / 2 * ( ( V2 - V3 ) / d * ( 1 / L0 - 1 / U.computeNorm() ) + U.at( 1 ) / pow( U.computeNorm(), 2 ) * ( V2 - V3 ) / pow( U.computeNorm(), 2 ) * ( ( x2 - x1 ) * ( V3 - V1 ) - ( x3 - x1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U3 - U1 ) / 2 * ( ( U3 - U2 ) / d * ( 1 / L0 - 1 / V.computeNorm() ) + V.at( 1 ) / pow( V.computeNorm(), 2 ) * ( U2 - U3 ) / pow( U.computeNorm(), 2 ) * ( ( x2 - x1 ) * ( U3 - U1 ) - ( x3 - x1 ) * ( U2 - U1 ) ) );
        //dFx2dy1
    answer.at( 4, 2 ) = Et * At * ( V1 - V3 ) / 2 * ( U.at( 1 ) * ( V2 - V3 ) / ( pow( U.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( y2 - y1 ) * ( V3 - V1 ) - ( y3 - y1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U3 - U1 ) / 2 * ( V.at( 1 ) * ( U2 - U3 ) / ( pow( V.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( y2 - y1 ) * ( U3 - U1 ) - ( y3 - y1 ) * ( U2 - U1 ) ) );
        //dFx2dz1
    answer.at( 4, 3 ) = Et * At * ( V1 - V3 ) / 2 * ( U.at( 1 ) * ( V2 - V3 ) / ( pow( U.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( z2 - z1 ) * ( V3 - V1 ) - ( z3 - z1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U3 - U1 ) / 2 * ( V.at( 1 ) * ( U2 - U3 ) / ( pow( V.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( z2 - z1 ) * ( U3 - U1 ) - ( z3 - z1 ) * ( U2 - U1 ) ) );
        //dFx2dx2
    answer.at( 4, 4 ) = Et * At * ( V1 - V3 ) / 2 * ( ( V3 - V1 ) / d * ( 1 / L0 - 1 / U.computeNorm() ) + U.at( 1 ) / pow( U.computeNorm(), 2 ) * ( V3 - V1 ) / pow( U.computeNorm(), 2 ) * ( ( x2 - x1 ) * ( V3 - V1 ) - ( x3 - x1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U3 - U1 ) / 2 * ( ( U1 - U3 ) / d * ( 1 / L0 - 1 / V.computeNorm() ) + V.at( 1 ) / pow( V.computeNorm(), 2 ) * ( U3 - U1 ) / pow( U.computeNorm(), 2 ) * ( ( x2 - x1 ) * ( U3 - U1 ) - ( x3 - x1 ) * ( U2 - U1 ) ) );
        //dFx2dy2
    answer.at( 4, 5 ) = Et * At * ( V1 - V3 ) / 2 * ( U.at( 1 ) * ( V3 - V1 ) / ( pow( U.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( y2 - y1 ) * ( V3 - V1 ) - ( y3 - y1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U3 - U1 ) / 2 * ( V.at( 1 ) * ( U3 - U1 ) / ( pow( V.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( y2 - y1 ) * ( U3 - U1 ) - ( y3 - y1 ) * ( U2 - U1 ) ) );
        //dFx2dz2
    answer.at( 4, 6 ) = Et * At * ( V1 - V3 ) / 2 * ( U.at( 1 ) * ( V3 - V1 ) / ( pow( U.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( z2 - z1 ) * ( V3 - V1 ) - ( z3 - z1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U3 - U1 ) / 2 * ( V.at( 1 ) * ( U3 - U1 ) / ( pow( V.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( z2 - z1 ) * ( U3 - U1 ) - ( z3 - z1 ) * ( U2 - U1 ) ) );
        //dFx2dx3
    answer.at( 4, 7 ) = Et * At * ( V1 - V3 ) / 2 * ( ( V1 - V2 ) / d * ( 1 / L0 - 1 / U.computeNorm() ) + U.at( 1 ) / pow( U.computeNorm(), 2 ) * ( V1 - V2 ) / pow( U.computeNorm(), 2 ) * ( ( x2 - x1 ) * ( V3 - V1 ) - ( x3 - x1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U3 - U1 ) / 2 * ( ( U2 - U1 ) / d * ( 1 / L0 - 1 / V.computeNorm() ) + V.at( 1 ) / pow( V.computeNorm(), 2 ) * ( U1 - U2 ) / pow( U.computeNorm(), 2 ) * ( ( x2 - x1 ) * ( U3 - U1 ) - ( x3 - x1 ) * ( U2 - U1 ) ) );
        //dFx2dy3
    answer.at( 4, 8 ) = Et * At * ( V1 - V3 ) / 2 * ( U.at( 1 ) * ( V1 - V2 ) / ( pow( U.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( y2 - y1 ) * ( V3 - V1 ) - ( y3 - y1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U3 - U1 ) / 2 * ( V.at( 1 ) * ( U1 - U2 ) / ( pow( V.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( y2 - y1 ) * ( U3 - U1 ) - ( y3 - y1 ) * ( U2 - U1 ) ) );
        //dFx2dz3
    answer.at( 4, 9 ) = Et * At * ( V1 - V3 ) / 2 * ( U.at( 1 ) * ( V1 - V2 ) / ( pow( U.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( z2 - z1 ) * ( V3 - V1 ) - ( z3 - z1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U3 - U1 ) / 2 * ( V.at( 1 ) * ( U1 - U2 ) / ( pow( V.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( z2 - z1 ) * ( U3 - U1 ) - ( z3 - z1 ) * ( U2 - U1 ) ) );

    //Fifth row in the matrix
        //dFy2dx1
    answer.at( 5, 1 ) = Et * At * ( V1 - V3 ) / 2 * ( U.at( 2 ) * ( V2 - V3 ) / ( pow( U.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( x2 - x1 ) * ( V3 - V1 ) - ( x3 - x1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U3 - U1 ) / 2 * ( V.at( 2 ) * ( U2 - U3 ) / ( pow( V.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( x2 - x1 ) * ( U3 - U1 ) - ( x3 - x1 ) * ( U2 - U1 ) ) );
        //dFy2dy1
    answer.at( 5, 2 ) = Et * At * ( V1 - V3 ) / 2 * ( ( V2 - V3 ) / d * ( 1 / L0 - 1 / U.computeNorm() ) + U.at( 2 ) / pow( U.computeNorm(), 2 ) * ( V2 - V3 ) / pow( U.computeNorm(), 2 ) * ( ( y2 - y1 ) * ( V3 - V1 ) - ( y3 - y1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U3 - U1 ) / 2 * ( ( U3 - U2 ) / d * ( 1 / L0 - 1 / V.computeNorm() ) + V.at( 2 ) / pow( V.computeNorm(), 2 ) * ( U2 - U3 ) / pow( U.computeNorm(), 2 ) * ( ( y2 - y1 ) * ( U3 - U1 ) - ( y3 - y1 ) * ( U2 - U1 ) ) );
        //dFy2dz1
    answer.at( 5, 3 ) = Et * At * ( V1 - V3 ) / 2 * ( U.at( 2 ) * ( V2 - V3 ) / ( pow( U.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( z2 - z1 ) * ( V3 - V1 ) - ( z3 - z1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U3 - U1 ) / 2 * ( V.at( 2 ) * ( U2 - U3 ) / ( pow( V.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( z2 - z1 ) * ( U3 - U1 ) - ( z3 - z1 ) * ( U2 - U1 ) ) );
        //dFy2dx2
    answer.at( 5, 4 ) = Et * At * ( V1 - V3 ) / 2 * ( U.at( 2 ) * ( V3 - V1 ) / ( pow( U.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( x2 - x1 ) * ( V3 - V1 ) - ( x3 - x1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U3 - U1 ) / 2 * ( V.at( 2 ) * ( U3 - U1 ) / ( pow( V.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( x2 - x1 ) * ( U3 - U1 ) - ( x3 - x1 ) * ( U2 - U1 ) ) );
        //dFy2dy2
    answer.at( 5, 5 ) = Et * At * ( V1 - V3 ) / 2 * ( ( V3 - V1 ) / d * ( 1 / L0 - 1 / U.computeNorm() ) + U.at( 2 ) / pow( U.computeNorm(), 2 ) * ( V3 - V1 ) / pow( U.computeNorm(), 2 ) * ( ( y2 - y1 ) * ( V3 - V1 ) - ( y3 - y1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U3 - U1 ) / 2 * ( ( U1 - U3 ) / d * ( 1 / L0 - 1 / V.computeNorm() ) + V.at( 2 ) / pow( V.computeNorm(), 2 ) * ( U3 - U1 ) / pow( U.computeNorm(), 2 ) * ( ( y2 - y1 ) * ( U3 - U1 ) - ( y3 - y1 ) * ( U2 - U1 ) ) );
        //dFy2dz3
    answer.at( 5, 6 ) = Et * At * ( V1 - V3 ) / 2 * ( U.at( 2 ) * ( V3 - V1 ) / ( pow( U.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( z2 - z1 ) * ( V3 - V1 ) - ( z3 - z1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U3 - U1 ) / 2 * ( V.at( 2 ) * ( U3 - U1 ) / ( pow( V.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( z2 - z1 ) * ( U3 - U1 ) - ( z3 - z1 ) * ( U2 - U1 ) ) );
        //dFy2dx3
    answer.at( 5, 7 ) = Et * At * ( V1 - V3 ) / 2 * ( U.at( 2 ) * ( V1 - V2 ) / ( pow( U.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( x2 - x1 ) * ( V3 - V1 ) - ( x3 - x1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U3 - U1 ) / 2 * ( V.at( 2 ) * ( U1 - U2 ) / ( pow( V.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( x2 - x1 ) * ( U3 - U1 ) - ( x3 - x1 ) * ( U2 - U1 ) ) );
        //dFy2dy3
    answer.at( 5, 8 ) = Et * At * ( V1 - V3 ) / 2 * ( ( V1 - V2 ) / d * ( 1 / L0 - 1 / U.computeNorm() ) + U.at( 2 ) / pow( U.computeNorm(), 2 ) * ( V1 - V2 ) / pow( U.computeNorm(), 2 ) * ( ( y2 - y1 ) * ( V3 - V1 ) - ( y3 - y1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U3 - U1 ) / 2 * ( ( U2 - U1 ) / d * ( 1 / L0 - 1 / V.computeNorm() ) + V.at( 2 ) / pow( V.computeNorm(), 2 ) * ( U1 - U2 ) / pow( U.computeNorm(), 2 ) * ( ( y2 - y1 ) * ( U3 - U1 ) - ( y3 - y1 ) * ( U2 - U1 ) ) );
        //dFy2dz3
    answer.at( 5, 9 ) = Et * At * ( V1 - V3 ) / 2 * ( U.at( 2 ) * ( V1 - V2 ) / ( pow( U.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( z2 - z1 ) * ( V3 - V1 ) - ( z3 - z1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U3 - U1 ) / 2 * ( V.at( 2 ) * ( U1 - U2 ) / ( pow( V.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( z2 - z1 ) * ( U3 - U1 ) - ( z3 - z1 ) * ( U2 - U1 ) ) );

    //Sixth row in the matrix
        //dFz2dx1
    answer.at( 6, 1 ) = Et * At * ( V1 - V3 ) / 2 * ( U.at( 3 ) * ( V2 - V3 ) / ( pow( U.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( x2 - x1 ) * ( V3 - V1 ) - ( x3 - x1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U3 - U1 ) / 2 * ( V.at( 3 ) * ( U2 - U3 ) / ( pow( V.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( x2 - x1 ) * ( U3 - U1 ) - ( x3 - x1 ) * ( U2 - U1 ) ) );
        //dFz2dy1
    answer.at( 6, 2 ) = Et * At * ( V1 - V3 ) / 2 * ( U.at( 3 ) * ( V2 - V3 ) / ( pow( U.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( y2 - y1 ) * ( V3 - V1 ) - ( y3 - y1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U3 - U1 ) / 2 * ( V.at( 3 ) * ( U2 - U3 ) / ( pow( V.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( y2 - y1 ) * ( U3 - U1 ) - ( y3 - y1 ) * ( U2 - U1 ) ) );
        //dFz2dz1
    answer.at( 6, 3 ) = Et * At * ( V1 - V3 ) / 2 * ( ( V2 - V3 ) / d * ( 1 / L0 - 1 / U.computeNorm() ) + U.at( 3 ) / pow( U.computeNorm(), 2 ) * ( V2 - V3 ) / pow( U.computeNorm(), 2 ) * ( ( z2 - z1 ) * ( V3 - V1 ) - ( z3 - z1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U3 - U1 ) / 2 * ( ( U3 - U2 ) / d * ( 1 / L0 - 1 / V.computeNorm() ) + V.at( 3 ) / pow( V.computeNorm(), 2 ) * ( U2 - U3 ) / pow( U.computeNorm(), 2 ) * ( ( z2 - z1 ) * ( U3 - U1 ) - ( z3 - z1 ) * ( U2 - U1 ) ) );
        //dFz2dx2
    answer.at( 6, 4 ) = Et * At * ( V1 - V3 ) / 2 * ( U.at( 3 ) * ( V3 - V1 ) / ( pow( U.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( x2 - x1 ) * ( V3 - V1 ) - ( x3 - x1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U3 - U1 ) / 2 * ( V.at( 3 ) * ( U3 - U1 ) / ( pow( V.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( x2 - x1 ) * ( U3 - U1 ) - ( x3 - x1 ) * ( U2 - U1 ) ) );
        //dFz2dy2
    answer.at( 6, 5 ) = Et * At * ( V1 - V3 ) / 2 * ( U.at( 3 ) * ( V3 - V1 ) / ( pow( U.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( y2 - y1 ) * ( V3 - V1 ) - ( y3 - y1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U3 - U1 ) / 2 * ( V.at( 3 ) * ( U3 - U1 ) / ( pow( V.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( y2 - y1 ) * ( U3 - U1 ) - ( y3 - y1 ) * ( U2 - U1 ) ) );
        //dFz2dz2
    answer.at( 6, 6 ) = Et * At * ( V1 - V3 ) / 2 * ( ( V3 - V1 ) / d * ( 1 / L0 - 1 / U.computeNorm() ) + U.at( 3 ) / pow( U.computeNorm(), 2 ) * ( V3 - V1 ) / pow( U.computeNorm(), 2 ) * ( ( z2 - z1 ) * ( V3 - V1 ) - ( z3 - z1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U3 - U1 ) / 2 * ( ( U1 - U3 ) / d * ( 1 / L0 - 1 / V.computeNorm() ) + V.at( 3 ) / pow( V.computeNorm(), 2 ) * ( U3 - U1 ) / pow( U.computeNorm(), 2 ) * ( ( z2 - z1 ) * ( U3 - U1 ) - ( z3 - z1 ) * ( U2 - U1 ) ) );
        //dFz2dx3
    answer.at( 6, 7 ) = Et * At * ( V1 - V3 ) / 2 * ( U.at( 3 ) * ( V1 - V2 ) / ( pow( U.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( x2 - x1 ) * ( V3 - V1 ) - ( x3 - x1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U3 - U1 ) / 2 * ( V.at( 3 ) * ( U1 - U2 ) / ( pow( V.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( x2 - x1 ) * ( U3 - U1 ) - ( x3 - x1 ) * ( U2 - U1 ) ) );
        //dFz2dy3
    answer.at( 6, 8 ) = Et * At * ( V1 - V3 ) / 2 * ( U.at( 3 ) * ( V1 - V2 ) / ( pow( U.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( y2 - y1 ) * ( V3 - V1 ) - ( y3 - y1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U3 - U1 ) / 2 * ( V.at( 3 ) * ( U1 - U2 ) / ( pow( V.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( y2 - y1 ) * ( U3 - U1 ) - ( y3 - y1 ) * ( U2 - U1 ) ) );
        //dFz2dz3
    answer.at( 6, 9 ) = Et * At * ( V1 - V3 ) / 2 * ( ( V1 - V2 ) / d * ( 1 / L0 - 1 / U.computeNorm() ) + U.at( 3 ) / pow( U.computeNorm(), 2 ) * ( V1 - V2 ) / pow( U.computeNorm(), 2 ) * ( ( z2 - z1 ) * ( V3 - V1 ) - ( z3 - z1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U3 - U1 ) / 2 * ( ( U2 - U1 ) / d * ( 1 / L0 - 1 / V.computeNorm() ) + V.at( 3 ) / pow( V.computeNorm(), 2 ) * ( U1 - U2 ) / pow( U.computeNorm(), 2 ) * ( ( z2 - z1 ) * ( U3 - U1 ) - ( z3 - z1 ) * ( U2 - U1 ) ) );

    //Seventh row in the matrix
        //dFx3dx1
    answer.at( 7, 1 ) = Et * At * ( V2 - V1 ) / 2 * ( ( V2 - V3 ) / d * ( 1 / L0 - 1 / U.computeNorm() ) + U.at( 1 ) / pow( U.computeNorm(), 2 ) * ( V2 - V3 ) / pow( U.computeNorm(), 2 ) * ( ( x2 - x1 ) * ( V3 - V1 ) - ( x3 - x1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U1 - U2 ) / 2 * ( ( U3 - U2 ) / d * ( 1 / L0 - 1 / V.computeNorm() ) + V.at( 1 ) / pow( V.computeNorm(), 2 ) * ( U2 - U3 ) / pow( U.computeNorm(), 2 ) * ( ( x2 - x1 ) * ( U3 - U1 ) - ( x3 - x1 ) * ( U2 - U1 ) ) );
        //dFx3dy1
    answer.at( 7, 2 ) = Et * At * ( V2 - V1 ) / 2 * ( U.at( 1 ) * ( V2 - V3 ) / ( pow( U.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( y2 - y1 ) * ( V3 - V1 ) - ( y3 - y1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U1 - U2 ) / 2 * ( V.at( 1 ) * ( U2 - U3 ) / ( pow( V.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( y2 - y1 ) * ( U3 - U1 ) - ( y3 - y1 ) * ( U2 - U1 ) ) );
        //dFx3dz1
    answer.at( 7, 3 ) = Et * At * ( V2 - V1 ) / 2 * ( U.at( 1 ) * ( V2 - V3 ) / ( pow( U.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( z2 - z1 ) * ( V3 - V1 ) - ( z3 - z1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U1 - U2 ) / 2 * ( V.at( 1 ) * ( U2 - U3 ) / ( pow( V.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( z2 - z1 ) * ( U3 - U1 ) - ( z3 - z1 ) * ( U2 - U1 ) ) );
        //dFx3dx2
    answer.at( 7, 4 ) = Et * At * ( V2 - V1 ) / 2 * ( ( V3 - V1 ) / d * ( 1 / L0 - 1 / U.computeNorm() ) + U.at( 1 ) / pow( U.computeNorm(), 2 ) * ( V3 - V1 ) / pow( U.computeNorm(), 2 ) * ( ( x2 - x1 ) * ( V3 - V1 ) - ( x3 - x1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U1 - U2 ) / 2 * ( ( U1 - U3 ) / d * ( 1 / L0 - 1 / V.computeNorm() ) + V.at( 1 ) / pow( V.computeNorm(), 2 ) * ( U3 - U1 ) / pow( U.computeNorm(), 2 ) * ( ( x2 - x1 ) * ( U3 - U1 ) - ( x3 - x1 ) * ( U2 - U1 ) ) );
        //dFx3dy2
    answer.at( 7, 5 ) = Et * At * ( V2 - V1 ) / 2 * ( U.at( 1 ) * ( V3 - V1 ) / ( pow( U.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( y2 - y1 ) * ( V3 - V1 ) - ( y3 - y1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U1 - U2 ) / 2 * ( V.at( 1 ) * ( U3 - U1 ) / ( pow( V.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( y2 - y1 ) * ( U3 - U1 ) - ( y3 - y1 ) * ( U2 - U1 ) ) );
        //dFx3dz2
    answer.at( 7, 6 ) = Et * At * ( V2 - V1 ) / 2 * ( U.at( 1 ) * ( V3 - V1 ) / ( pow( U.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( z2 - z1 ) * ( V3 - V1 ) - ( z3 - z1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U1 - U2 ) / 2 * ( V.at( 1 ) * ( U3 - U1 ) / ( pow( V.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( z2 - z1 ) * ( U3 - U1 ) - ( z3 - z1 ) * ( U2 - U1 ) ) );
        //dFx3dx3
    answer.at( 7, 7 ) = Et * At * ( V2 - V1 ) / 2 * ( ( V1 - V2 ) / d * ( 1 / L0 - 1 / U.computeNorm() ) + U.at( 1 ) / pow( U.computeNorm(), 2 ) * ( V1 - V2 ) / pow( U.computeNorm(), 2 ) * ( ( x2 - x1 ) * ( V3 - V1 ) - ( x3 - x1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U1 - U2 ) / 2 * ( ( U2 - U1 ) / d * ( 1 / L0 - 1 / V.computeNorm() ) + V.at( 1 ) / pow( V.computeNorm(), 2 ) * ( U1 - U2 ) / pow( U.computeNorm(), 2 ) * ( ( x2 - x1 ) * ( U3 - U1 ) - ( x3 - x1 ) * ( U2 - U1 ) ) );
        //dFx3dy3
    answer.at( 7, 8 ) = Et * At * ( V2 - V1 ) / 2 * ( U.at( 1 ) * ( V1 - V2 ) / ( pow( U.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( y2 - y1 ) * ( V3 - V1 ) - ( y3 - y1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U1 - U2 ) / 2 * ( V.at( 1 ) * ( U1 - U2 ) / ( pow( V.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( y2 - y1 ) * ( U3 - U1 ) - ( y3 - y1 ) * ( U2 - U1 ) ) );
        //dFx3dz3
    answer.at( 7, 9 ) = Et * At * ( V2 - V1 ) / 2 * ( U.at( 1 ) * ( V1 - V2 ) / ( pow( U.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( z2 - z1 ) * ( V3 - V1 ) - ( z3 - z1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U1 - U2 ) / 2 * ( V.at( 1 ) * ( U1 - U2 ) / ( pow( V.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( z2 - z1 ) * ( U3 - U1 ) - ( z3 - z1 ) * ( U2 - U1 ) ) );

    //Eigth row in the matrix
        //dFy3dx1
    answer.at( 8, 1 ) = Et * At * ( V2 - V1 ) / 2 * ( U.at( 2 ) * ( V2 - V3 ) / ( pow( U.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( x2 - x1 ) * ( V3 - V1 ) - ( x3 - x1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U1 - U2 ) / 2 * ( V.at( 2 ) * ( U2 - U3 ) / ( pow( V.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( x2 - x1 ) * ( U3 - U1 ) - ( x3 - x1 ) * ( U2 - U1 ) ) );
        //dFy3dy1
    answer.at( 8, 2 ) = Et * At * ( V2 - V1 ) / 2 * ( ( V2 - V3 ) / d * ( 1 / L0 - 1 / U.computeNorm() ) + U.at( 2 ) / pow( U.computeNorm(), 2 ) * ( V2 - V3 ) / pow( U.computeNorm(), 2 ) * ( ( y2 - y1 ) * ( V3 - V1 ) - ( y3 - y1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U1 - U2 ) / 2 * ( ( U3 - U2 ) / d * ( 1 / L0 - 1 / V.computeNorm() ) + V.at( 2 ) / pow( V.computeNorm(), 2 ) * ( U2 - U3 ) / pow( U.computeNorm(), 2 ) * ( ( y2 - y1 ) * ( U3 - U1 ) - ( y3 - y1 ) * ( U2 - U1 ) ) );
        //dFy3dz1
    answer.at( 8, 3 ) = Et * At * ( V2 - V1 ) / 2 * ( U.at( 2 ) * ( V2 - V3 ) / ( pow( U.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( z2 - z1 ) * ( V3 - V1 ) - ( z3 - z1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U1 - U2 ) / 2 * ( V.at( 2 ) * ( U2 - U3 ) / ( pow( V.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( z2 - z1 ) * ( U3 - U1 ) - ( z3 - z1 ) * ( U2 - U1 ) ) );
        //dFy3dx2
    answer.at( 8, 4 ) = Et * At * ( V2 - V1 ) / 2 * ( U.at( 2 ) * ( V3 - V1 ) / ( pow( U.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( x2 - x1 ) * ( V3 - V1 ) - ( x3 - x1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U1 - U2 ) / 2 * ( V.at( 2 ) * ( U3 - U1 ) / ( pow( V.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( x2 - x1 ) * ( U3 - U1 ) - ( x3 - x1 ) * ( U2 - U1 ) ) );
        //dFy3dy2
    answer.at( 8, 5 ) = Et * At * ( V2 - V1 ) / 2 * ( ( V3 - V1 ) / d * ( 1 / L0 - 1 / U.computeNorm() ) + U.at( 2 ) / pow( U.computeNorm(), 2 ) * ( V3 - V1 ) / pow( U.computeNorm(), 2 ) * ( ( y2 - y1 ) * ( V3 - V1 ) - ( y3 - y1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U1 - U2 ) / 2 * ( ( U1 - U3 ) / d * ( 1 / L0 - 1 / V.computeNorm() ) + V.at( 2 ) / pow( V.computeNorm(), 2 ) * ( U3 - U1 ) / pow( U.computeNorm(), 2 ) * ( ( y2 - y1 ) * ( U3 - U1 ) - ( y3 - y1 ) * ( U2 - U1 ) ) );
        //dFy3dz2
    answer.at( 8, 6 ) = Et * At * ( V2 - V1 ) / 2 * ( U.at( 2 ) * ( V3 - V1 ) / ( pow( U.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( z2 - z1 ) * ( V3 - V1 ) - ( z3 - z1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U1 - U2 ) / 2 * ( V.at( 2 ) * ( U3 - U1 ) / ( pow( V.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( z2 - z1 ) * ( U3 - U1 ) - ( z3 - z1 ) * ( U2 - U1 ) ) );
        //dFy3dx3
    answer.at( 8, 7 ) = Et * At * ( V2 - V1 ) / 2 * ( U.at( 2 ) * ( V1 - V2 ) / ( pow( U.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( x2 - x1 ) * ( V3 - V1 ) - ( x3 - x1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U1 - U2 ) / 2 * ( V.at( 2 ) * ( U1 - U2 ) / ( pow( V.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( x2 - x1 ) * ( U3 - U1 ) - ( x3 - x1 ) * ( U2 - U1 ) ) );
        //dFy3dy3
    answer.at( 8, 8 ) = Et * At * ( V2 - V1 ) / 2 * ( ( V1 - V2 ) / d * ( 1 / L0 - 1 / U.computeNorm() ) + U.at( 2 ) / pow( U.computeNorm(), 2 ) * ( V1 - V2 ) / pow( U.computeNorm(), 2 ) * ( ( y2 - y1 ) * ( V3 - V1 ) - ( y3 - y1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U1 - U2 ) / 2 * ( ( U2 - U1 ) / d * ( 1 / L0 - 1 / V.computeNorm() ) + V.at( 2 ) / pow( V.computeNorm(), 2 ) * ( U1 - U2 ) / pow( U.computeNorm(), 2 ) * ( ( y2 - y1 ) * ( U3 - U1 ) - ( y3 - y1 ) * ( U2 - U1 ) ) );
        //dFy3dz3
    answer.at( 8, 9 ) = Et * At * ( V2 - V1 ) / 2 * ( U.at( 2 ) * ( V1 - V2 ) / ( pow( U.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( z2 - z1 ) * ( V3 - V1 ) - ( z3 - z1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U1 - U2 ) / 2 * ( V.at( 2 ) * ( U1 - U2 ) / ( pow( V.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( z2 - z1 ) * ( U3 - U1 ) - ( z3 - z1 ) * ( U2 - U1 ) ) );

    //Ninth row in the matrix
        //dFz3dx1
    answer.at( 9, 1 ) = Et * At * ( V2 - V1 ) / 2 * ( U.at( 3 ) * ( V2 - V3 ) / ( pow( U.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( x2 - x1 ) * ( V3 - V1 ) - ( x3 - x1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U1 - U2 ) / 2 * ( V.at( 3 ) * ( U2 - U3 ) / ( pow( V.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( x2 - x1 ) * ( U3 - U1 ) - ( x3 - x1 ) * ( U2 - U1 ) ) );
        //dFz3dy1
    answer.at( 9, 2 ) = Et * At * ( V2 - V1 ) / 2 * ( U.at( 3 ) * ( V2 - V3 ) / ( pow( U.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( y2 - y1 ) * ( V3 - V1 ) - ( y3 - y1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U1 - U2 ) / 2 * ( V.at( 3 ) * ( U2 - U3 ) / ( pow( V.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( y2 - y1 ) * ( U3 - U1 ) - ( y3 - y1 ) * ( U2 - U1 ) ) );
        //dFz3dz1
    answer.at( 9, 3 ) = Et * At * ( V2 - V1 ) / 2 * ( ( V2 - V3 ) / d * ( 1 / L0 - 1 / U.computeNorm() ) + U.at( 3 ) / pow( U.computeNorm(), 2 ) * ( V2 - V3 ) / pow( U.computeNorm(), 2 ) * ( ( z2 - z1 ) * ( V3 - V1 ) - ( z3 - z1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U1 - U2 ) / 2 * ( ( U3 - U1 ) / d * ( 1 / L0 - 1 / V.computeNorm() ) + V.at( 3 ) / pow( V.computeNorm(), 2 ) * ( U2 - U3 ) / pow( U.computeNorm(), 2 ) * ( ( z2 - z1 ) * ( U3 - U1 ) - ( z3 - z1 ) * ( U2 - U1 ) ) );
        //dFz3dx2
    answer.at( 9, 4 ) = Et * At * ( V2 - V1 ) / 2 * ( U.at( 3 ) * ( V3 - V1 ) / ( pow( U.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( x2 - x1 ) * ( V3 - V1 ) - ( x3 - x1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U1 - U2 ) / 2 * ( V.at( 3 ) * ( U3 - U1 ) / ( pow( V.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( x2 - x1 ) * ( U3 - U1 ) - ( x3 - x1 ) * ( U2 - U1 ) ) );
        //dFz3dy2
    answer.at( 9, 5 ) = Et * At * ( V2 - V1 ) / 2 * ( U.at( 3 ) * ( V3 - V1 ) / ( pow( U.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( y2 - y1 ) * ( V3 - V1 ) - ( y3 - y1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U1 - U2 ) / 2 * ( V.at( 3 ) * ( U3 - U1 ) / ( pow( V.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( y2 - y1 ) * ( U3 - U1 ) - ( y3 - y1 ) * ( U2 - U1 ) ) );
        //dFz3dz2
    answer.at( 9, 6 ) = Et * At * ( V2 - V1 ) / 2 * ( ( V3 - V1 ) / d * ( 1 / L0 - 1 / U.computeNorm() ) + U.at( 3 ) / pow( U.computeNorm(), 2 ) * ( V3 - V1 ) / pow( U.computeNorm(), 2 ) * ( ( z2 - z1 ) * ( V3 - V1 ) - ( z3 - z1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U1 - U2 ) / 2 * ( ( U1 - U3 ) / d * ( 1 / L0 - 1 / V.computeNorm() ) + V.at( 3 ) / pow( V.computeNorm(), 2 ) * ( U3 - U1 ) / pow( U.computeNorm(), 2 ) * ( ( z2 - z1 ) * ( U3 - U1 ) - ( z3 - z1 ) * ( U2 - U1 ) ) );
        //dFz3dx3
    answer.at( 9, 7 ) = Et * At * ( V2 - V1 ) / 2 * ( U.at( 3 ) * ( V1 - V2 ) / ( pow( U.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( x2 - x1 ) * ( V3 - V1 ) - ( x3 - x1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U1 - U2 ) / 2 * ( V.at( 3 ) * ( U1 - U2 ) / ( pow( V.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( x2 - x1 ) * ( U3 - U1 ) - ( x3 - x1 ) * ( U2 - U1 ) ) );
        //dFz3dy3
    answer.at( 9, 8 ) = Et * At * ( V2 - V1 ) / 2 * ( U.at( 3 ) * ( V1 - V2 ) / ( pow( U.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( y2 - y1 ) * ( V3 - V1 ) - ( y3 - y1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U1 - U2 ) / 2 * ( V.at( 3 ) * ( U1 - U2 ) / ( pow( V.computeNorm(), 2 ) * pow( U.computeNorm(), 2 ) ) * ( ( y2 - y1 ) * ( U3 - U1 ) - ( y3 - y1 ) * ( U2 - U1 ) ) );
        //dFz3dz3
    answer.at( 9, 9 ) = Et * At * ( V2 - V1 ) / 2 * ( ( V1 - V2 ) / d * ( 1 / L0 - 1 / U.computeNorm() ) + U.at( 3 ) / pow( U.computeNorm(), 2 ) * ( V1 - V2 ) / pow( U.computeNorm(), 2 ) * ( ( z2 - z1 ) * ( V3 - V1 ) - ( z3 - z1 ) * ( V2 - V1 ) ) )
        + Et * At * ( U1 - U2 ) / 2 * ( ( U2 - U1 ) / d * ( 1 / L0 - 1 / V.computeNorm() ) + V.at( 3 ) / pow( V.computeNorm(), 2 ) * ( U1 - U2 ) / pow( U.computeNorm(), 2 ) * ( ( z2 - z1 ) * ( U3 - U1 ) - ( z3 - z1 ) * ( U2 - U1 ) ) );

    for ( int i = 1; i <= 9; i++ ) {
        for ( int j = 1; j <= 9; j++ )
            answer.at( i, j ) = -1 * answer.at( i, j );
    }
}

FloatArray
NetTr3Pr::computeUTwine( TimeStep *tStep )
{
    // Fetch current coordinates of element's nodes
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

    // Calculate strains in U- and V-twines
    answer.at( 1 ) = ( lengthU - L0 ) / L0;
    answer.at( 2 ) = ( lengthV - L0 ) / L0;
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

    d = ( U2 - U1 ) * ( V3 - V1 ) - ( U3 - U1 ) * ( V2 - V1 );
}
}