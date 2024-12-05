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
NetRec4TrLa::calculateElementDimensions( TimeStep *tStep )
{
    // Fetch initial position of the nodes
    FloatArray node1 = this->giveNode( 1 )->giveCoordinates();
    FloatArray node2 = this->giveNode( 2 )->giveCoordinates();
    FloatArray node4 = this->giveNode( 4 )->giveCoordinates();
    
    if ( !tStep->isTheFirstStep() ) {
        // Fetch total displacements in the current configuration
        FloatArray u;
        this->computeVectorOf( VM_Total, tStep, u );

        // Update nodes' positions
        node1.at( 1 ) = node1.at( 1 ) + u.at( 1 );
        node1.at( 2 ) = node1.at( 2 ) + u.at( 2 );
        node2.at( 1 ) = node2.at( 1 ) + u.at( 3 );
        node2.at( 2 ) = node2.at( 2 ) + u.at( 4 );
        node4.at( 1 ) = node4.at( 1 ) + u.at( 7 );
        node4.at( 2 ) = node4.at( 2 ) + u.at( 8 );
    }

    // Calculate length of the element's sides in the current configuration
    FloatArray xSide, ySide;
    xSide.beDifferenceOf( node2, node1 );
    ySide.beDifferenceOf( node4, node1 );
    FloatArray dimensions;
    dimensions.resize( 2 );
    dimensions.at( 1 ) = xSide.computeNorm();
    dimensions.at( 2 ) = ySide.computeNorm();

    return dimensions;
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

bool
NetRec4TrLa::computeGtoLRotationMatrix( FloatMatrix &answer )
// Returns the rotation matrix of the receiver of the size [8,12]
// r(local) = T * r(global)
// for one node (r written transposed): {u,v} = T * {u,v,w}
{
    if ( GtoLRotationMatrix == NULL)
        NetElement::computeGtoLRotationMatrix();

    answer.resize( 8, 12 );
    answer.zero();

    for ( int i = 1; i <= 3; i++ ) {
        answer.at( 1, i ) = answer.at( 3, i + 3 ) = answer.at( 5, i + 6 ) = answer.at( 7, i + 9 ) = GtoLRotationMatrix->at( 1, i );
        answer.at( 2, i ) = answer.at( 4, i + 3 ) = answer.at( 6, i + 6 ) = answer.at( 8, i + 9 ) = GtoLRotationMatrix->at( 2, i );
    }

    return 1;
}

void
NetRec4TrLa ::computeLumpedMassMatrix( FloatMatrix &answer, TimeStep *tStep )
// Returns the lumped mass matrix of the receiver. This expression is
// valid in both local and global axes.
{
    answer.resize( 8, 8 );
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
    double noOfTwines = 2 * initialDimensions.at( 1 ) / a0 * initialDimensions.at( 2 ) / a0;

    GaussPoint *gp    = integrationRulesArray[0]->getIntegrationPoint( 0 );
    double density    = this->giveStructuralCrossSection()->give( 'd', gp );
    double quarterMass   = 0.25 * noOfTwines * density * this->giveCrossSection()->give( CS_Area, gp ) * a0;
    answer.at( 1, 1 ) = answer.at( 2, 2 ) = answer.at( 3, 3 ) = answer.at( 4, 4 ) = answer.at( 5, 5 ) = answer.at( 6, 6 ) = answer.at( 7, 7 ) = answer.at( 8, 8 ) = quarterMass;
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
    
    if ( mask == 1 ) {
        // Calculate element's stiffness matrix
        answer.resize( 8, 8 );
        answer.at( 1, 1 ) = answer.at( 3, 3 ) = answer.at( 5, 5 ) = answer.at( 7, 7 ) = At * Et / ( initialDimensions.at( 1 ) + internalDisplacements.at( 1 ) ) + a0 * At * Et * ( 1 / a0 - initialDimensions.at( 1 ) / ( a0 * ( initialDimensions.at( 1 ) + internalDisplacements.at( 1 ) ) ) ) / initialDimensions.at( 1 );
        answer.at( 1, 3 ) = answer.at( 3, 1 ) = answer.at( 5, 7 ) = answer.at( 7, 5 ) = -1 * ( At * Et / ( initialDimensions.at( 1 ) + internalDisplacements.at( 1 ) ) + a0 * At * Et * ( 1 / a0 - initialDimensions.at( 1 ) / ( a0 * ( initialDimensions.at( 1 ) + internalDisplacements.at( 1 ) ) ) ) / initialDimensions.at( 1 ) );
        answer.at( 2, 2 ) = answer.at( 4, 4 ) = answer.at( 6, 6 ) = answer.at( 8, 8 ) = At * Et / ( initialDimensions.at( 2 ) + internalDisplacements.at( 2 ) ) + a0 * At * Et * ( 1 / a0 - initialDimensions.at( 2 ) / ( a0 * ( initialDimensions.at( 2 ) + internalDisplacements.at( 2 ) ) ) ) / initialDimensions.at( 2 );
        answer.at( 2, 8 ) = answer.at( 4, 6 ) = answer.at( 6, 4 ) = answer.at( 8, 4 ) = -1 * ( At * Et / ( initialDimensions.at( 2 ) + internalDisplacements.at( 2 ) ) + a0 * At * Et * ( 1 / a0 - initialDimensions.at( 2 ) / ( a0 * ( initialDimensions.at( 2 ) + internalDisplacements.at( 2 ) ) ) ) / initialDimensions.at( 2 ) );
    }
}

void
NetRec4TrLa::giveDofManDofIDMask( int inode, IntArray &answer ) const
{
    answer = { D_u, D_v, D_w };
}

void
NetRec4TrLa::giveInternalForcesVector( FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord )
{
    FloatArray strains, stresses;
    double At = 0;
    for ( auto gp : *this->giveDefaultIntegrationRulePtr() ) {
        // Calculate strains and stresses in x- and y-twines
        computeStrainVector( strains, gp, tStep );
        computeStressVector( stresses, strains, gp, tStep );
        At = this->giveCrossSection()->give( CS_Area, gp );
    }

    // Calculate internal forces in x- and y-direction
    double Fx = At * stresses.at( 1 );
    double Fy = At * stresses.at( 2 );

    // Calculate internal forces at the nodes
    answer.resize( 8 );
    answer.at( 1 ) = answer.at( 7 ) = -Fx;
    answer.at( 3 ) = answer.at( 5 ) = Fx;
    answer.at( 2 ) = answer.at( 4 ) = -Fy;
    answer.at( 6 ) = answer.at( 8 ) = Fy;
}

FEInterpolation *
NetRec4TrLa::giveInterpolation() const
{
    return &interpolation;
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