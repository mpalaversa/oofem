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

#include "nodalload.h"
#include "classfactory.h"
#include "datastream.h"
#include "dynamicinputrecord.h"
#include "contextioerr.h"
#include "CrossSections/decoupledcrosssection.h"

namespace oofem {
REGISTER_BoundaryCondition(NodalLoad);

NodalLoad :: NodalLoad( int n, Domain *d ) : Load( n, d ) {
    sn = 0.0;
    a3 = 0.05;
    b4 = 0.05;
}

void NodalLoad ::computeHydrodynamicKFLoad( FloatArray &load, std::vector<Element*> associatedElements, DofManager *node, TimeStep *tStep )
{
    // Check if solidity ratio is defined. This hydrodynamic model cannot work without it.
    if ( this->sn == 0 )
        OOFEM_ERROR( "The Kristiansen-Faltinsen hydrodynamic model cannot work without a non-zero value for the net solidity ratio being specified." );
    // This currently works only for elements with 2 nodes in 3D.
    std::vector<FloatArray> elementVectors;
    FloatArray nodes, vTemp;
    int nodeNumber = node->giveNumber();
    // The following block of code organizes edges of the panel in elementVectors such that
    // the current node (given by nodeNumber) is always the first point of the edge.
    for ( int i = 0; i < associatedElements.size(); i++ ) {
        nodes.resize( associatedElements.at( i )->giveNumberOfNodes() * 3 );
        for ( int j = 0; j < associatedElements.at( i )->giveNumberOfNodes(); j++ ) {
            Node *elNode = associatedElements.at( i )->giveNode( j + 1 );
            if ( elNode->giveNumber() == nodeNumber ) {
                nodes.at( 1 ) = elNode->giveCoordinate( 1 );
                nodes.at( 2 ) = elNode->giveCoordinate( 2 );
                nodes.at( 3 ) = elNode->giveCoordinate( 3 );
            } else {
                nodes.at( 4 ) = elNode->giveCoordinate( 1 );
                nodes.at( 5 ) = elNode->giveCoordinate( 2 );
                nodes.at( 6 ) = elNode->giveCoordinate( 3 );
            }
        }
        vTemp.resize( 6 );
        vTemp.at( 1 ) = nodes.at( 4 ) - nodes.at( 1 );
        vTemp.at( 2 ) = nodes.at( 5 ) - nodes.at( 2 );
        vTemp.at( 3 ) = nodes.at( 6 ) - nodes.at( 3 );
        DecoupledCrossSection *dcs = associatedElements.at( i )->giveDecoupledCrossSectionOfType( DecoupledMaterial::DecoupledMaterialType::DecoupledFluidMaterial );
        double characteristicDim   = dcs->giveCharacteristicDimension();
        vTemp.at( 4 ) = characteristicDim;
        double density             = dcs->giveMagnitudeOfMaterialProperty( 'd' );
        vTemp.at( 5 )              = density;
        double mu                  = dcs->giveMaterial()->giveDynamicViscosity();
        vTemp.at( 6 )              = mu;
        elementVectors.push_back( vTemp );
    }
    // This block of code finds edges that belong to the same panel. These edge pairs are stored in panelEdges.
    std::vector<IntArray> panelEdges;
    IntArray panelEdgesTemp;
    // Iterate thorugh all potential panel edges stored in elementVectors
    for ( int i = 0; i < elementVectors.size(); i++ ) {
        // For each potential panel edge, iterate thorugh all other potential panel edges
        for ( int j = 0; j < elementVectors.size(); j++ ) {
            // If the potential panel edge j is not i, if the potential edges are not vectors ling on the same line
            // and if j is not already paired with i
            FloatArray basisVector1, basisVector2, basisVector3;
            basisVector1.resize( 3 );
            basisVector2.resize( 3 );
            basisVector1.at( 1 ) = elementVectors.at( i ).at( 1 );
            basisVector1.at( 2 ) = elementVectors.at( i ).at( 2 );
            basisVector1.at( 3 ) = elementVectors.at( i ).at( 3 );
            basisVector2.at( 1 ) = elementVectors.at( j ).at( 1 );
            basisVector2.at( 2 ) = elementVectors.at( j ).at( 2 );
            basisVector2.at( 3 ) = elementVectors.at( j ).at( 3 );
            basisVector3.beVectorProductOf( basisVector1, basisVector2 );
            bool areParallel = basisVector1.isParallelTo( basisVector2 );
            if ( i != j && basisVector3.computeNorm() != 0 && isEdgePairAlreadyFound( panelEdges, j, i ) == false ) {
                FloatMatrix matOfBasisVectors;
                matOfBasisVectors.resize( 3, 3 );
                matOfBasisVectors.at( 1, 1 ) = elementVectors.at( i ).at( 1 );
                matOfBasisVectors.at( 2, 1 ) = elementVectors.at( i ).at( 2 );
                matOfBasisVectors.at( 3, 1 ) = elementVectors.at( i ).at( 3 );
                matOfBasisVectors.at( 1, 2 ) = elementVectors.at( j ).at( 1 );
                matOfBasisVectors.at( 2, 2 ) = elementVectors.at( j ).at( 2 );
                matOfBasisVectors.at( 3, 2 ) = elementVectors.at( j ).at( 3 );
                matOfBasisVectors.at( 1, 3 ) = basisVector3.at( 1 );
                matOfBasisVectors.at( 2, 3 ) = basisVector3.at( 2 );
                matOfBasisVectors.at( 3, 3 ) = basisVector3.at( 3 );
                // Check whether they are edges of a net panel, i.e. there are no other elements between the two.
                FloatArray thirdVector, coefficients;
                thirdVector.resize( 3 );
                panelEdgesTemp.resize( 2 );
                // Iterate thorugh all edges/elements that are neither element i nor element j
                for ( int k = 0; k < elementVectors.size(); k++ ) {
                    if ( i != k && j != k ) {
                        // Element k is the third vector - the first two being element i and element j stored in matOfBasisVectors
                        thirdVector.at( 1 ) = elementVectors.at( k ).at( 1 );
                        thirdVector.at( 2 ) = elementVectors.at( k ).at( 2 );
                        thirdVector.at( 3 ) = elementVectors.at( k ).at( 3 );
                        // Solves for the coefficients/projections of thirdVector onto elementVectors.at(i) and elementVectors.at(j) 
                        matOfBasisVectors.solveForRhs( thirdVector, coefficients );
                        // If the resulting coefficients are not all positive, this means that thirdVector doesn't lie in the first quadrant
                        // of a vector space spanned by elementVectors.at(i) and elementVectors.at(j), i.e. the element represented by
                        // thirdVector is not between elements given by vectors elementVectors.at(i) and elementVectors.at(j).
                        if ( !( coefficients.at( 1 ) > 0 && coefficients.at( 2 ) > 0 && coefficients.at( 3 ) == 0 ) ) {
                            panelEdgesTemp.at( 1 ) = i;
                            panelEdgesTemp.at( 2 ) = j;
                            panelEdges.push_back( panelEdgesTemp );
                            break;
                        }
                    } else if ( k == elementVectors.size() - 1 ) {
                        panelEdgesTemp.at( 1 ) = i;
                        panelEdgesTemp.at( 2 ) = j;
                        panelEdges.push_back( panelEdgesTemp );
                    }
                }
            }
        }
    }

    // Iterate thorugh panelEdges. Number of edge pairs in panelEdges is equal to number of panels associated with the node given by nodeNumber.
    FloatArray v1, v2, nodalVelocity;
    // This works just for Truss3d
    IntArray dofIDMask = { D_u, D_v, D_w };
    node->giveUnknownVector( nodalVelocity, dofIDMask, VM_Velocity, tStep, true );

    FloatArray relativeVelocity;
    relativeVelocity.beDifferenceOf( this->giveComponentArray(), nodalVelocity );
    
    double inflowAngle, characteristicDim, area, density, mu, cn, ct, cdc;
    FloatArray v1xv2, en, et, normalViscousForce, tangentialViscousForce, viscousForce;
    load.resize( 3 );
    for ( int i = 0; i < panelEdges.size(); i++ ) {
        // Form vectors v1 and v2 representing panel edges
        v1.resize( 3 );
        v1.at( 1 ) = elementVectors.at( panelEdges.at( i ).at( 1 ) ).at( 1 );
        v1.at( 2 ) = elementVectors.at( panelEdges.at( i ).at( 1 ) ).at( 2 );
        v1.at( 3 ) = elementVectors.at( panelEdges.at( i ).at( 1 ) ).at( 3 );
        v2.resize( 3 );
        v2.at( 1 ) = elementVectors.at( panelEdges.at( i ).at( 2 ) ).at( 1 );
        v2.at( 2 ) = elementVectors.at( panelEdges.at( i ).at( 2 ) ).at( 2 );
        v2.at( 3 ) = elementVectors.at( panelEdges.at( i ).at( 2 ) ).at( 3 );

        characteristicDim = ( elementVectors.at( panelEdges.at( i ).at( 1 ) ).at( 4 ) + elementVectors.at( panelEdges.at( i ).at( 2 ) ).at( 4 ) ) / 2;
        density = ( elementVectors.at( panelEdges.at( i ).at( 1 ) ).at( 5 ) + elementVectors.at( panelEdges.at( i ).at( 2 ) ).at( 5 ) ) / 2;
        mu           = ( elementVectors.at( panelEdges.at( i ).at( 1 ) ).at( 6 ) + elementVectors.at( panelEdges.at( i ).at( 2 ) ).at( 6 ) ) / 2;

        v1xv2.beVectorProductOf( v1, v2 );
        area = 0.25 * v1xv2.computeNorm();

        // Calculate normal to the current panel
        en.beScaled( 1 / v1xv2.computeNorm(), v1xv2 );
        // Assign the sign to the normal such that it points in the same falf-space as the relative velocity (see the paper for more info)
        en.beScaled( relativeVelocity.dotProduct( en ) / abs( relativeVelocity.dotProduct( en ) ), en );

        cdc         = computeDragCoefficientOfCircularCylinder( density, mu, characteristicDim, relativeVelocity.computeNorm(), this->sn );
        inflowAngle = acos( relativeVelocity.dotProduct( en ) / relativeVelocity.computeNorm() );

        if ( inflowAngle != 0 )
            et.beScaled( 1 / ( relativeVelocity - relativeVelocity.dotProduct( en ) * en ).computeNorm(), relativeVelocity - relativeVelocity.dotProduct( en ) * en );
        else
            et.resize( 3 );

        if ( inflowAngle <= 3.15 / 4 ) {
            cn  = cdc * this->sn * ( 2 - this->sn ) / ( 2 * pow( 1 - this->sn, 2 ) ) * pow( cos( inflowAngle ), 2 );
            ct  = 4 * cn * inflowAngle / ( 8 + cn );
        } else if ( inflowAngle <= 3.15 / 2) {
            double cd, cl;
            cd = cdc * this->sn * ( 2 - this->sn ) / ( 2 * pow( 1 - this->sn, 2 ) ) * ( cos( inflowAngle ) + this->a3 * ( cos( 3 * inflowAngle ) - cos( inflowAngle ) ) );
            double cn45 = cdc * this->sn * ( 2 - this->sn ) / ( 2 * pow( 1 - this->sn, 2 ) ) * pow( cos( 3.14 / 4 ), 2 );
            cl = ( 0.5 * cdc * this->sn * ( 2 - this->sn ) / ( 2 * pow( 1 - this->sn, 2 ) ) - 4 * cn45 * ( 3.14 / 4 ) / ( 8 + cn45 ) ) / sqrt( 2 ) * ( sin( 2 * inflowAngle ) + this->b4 * sin( 4 * inflowAngle ) );

            FloatMatrix A;
            A.resize( 2, 2 );
            A.at( 1, 1 ) = cos( inflowAngle );
            A.at( 1, 2 ) = sin( inflowAngle );
            A.at( 2, 1 ) = sin( inflowAngle );
            A.at( 2, 2 ) = -cos( inflowAngle );
            FloatArray b;
            b.resize( 2 );
            b.at( 1 ) = cd;
            b.at( 2 ) = cl;
            FloatArray x;
            A.solveForRhs( b, x );

            cn = x.at( 1 );
            ct = x.at( 2 );
        } else
            OOFEM_ERROR( "Inflow angle is greater than 90 degrees. The Kristiansen-Faltinsen hydrodynamic model is not applicable to inflow angles greater than 90 degrees." );

        normalViscousForce.beScaled( 0.5 * cn * density * area * relativeVelocity.dotProduct( relativeVelocity ), en );
        tangentialViscousForce.beScaled( 0.5 * ct * density * area * relativeVelocity.dotProduct( relativeVelocity ), et );

        viscousForce.resize( 3 );
        viscousForce.add( normalViscousForce );
        viscousForce.add( tangentialViscousForce );
        load.add( viscousForce );
    }
}

bool NodalLoad::isEdgePairAlreadyFound( std::vector<IntArray> panelEdges, int edge1, int edge2 )
{
    if ( panelEdges.size() > 0 ) {
        for ( int i = 0; i < panelEdges.size(); i++ ) {
            if ( panelEdges.at( i ).at( 1 ) == edge1 && panelEdges.at( i ).at( 2 ) == edge2 )
                return true;
        }
        return false;
    } else
        return false;
}

double NodalLoad::computeDragCoefficientOfCircularCylinder( double density, double mu, double dt, double relativeVelocity, double sn )
{
    double rn = density * dt * relativeVelocity / ( mu * ( 1 - sn ) );

    return -78.46675 + 254.73873 * log10( rn ) - 327.8864 * pow( log10( rn ), 2 ) + 223.64577 * pow( log10( rn ), 4 ) + 20.00769 * pow( log10( rn ), 5 ) - 2.44894 * pow( log10( rn ), 6 ) + 0.12479 * pow( log10( rn ), 7 );
}


void
NodalLoad :: initializeFrom(InputRecord &ir)
{
    Load :: initializeFrom(ir);
    int value = 1;
    IR_GIVE_OPTIONAL_FIELD(ir, value, _IFT_NodalLoad_cstype);
    coordSystemType = ( CoordSystType ) value;

    value = 0;
    IR_GIVE_OPTIONAL_FIELD( ir, value, _IFT_NodalLoad_loadtype );
    lType = (bcType)value;

    IR_GIVE_OPTIONAL_FIELD( ir, sn, _IFT_NodalLoad_sn );
    IR_GIVE_OPTIONAL_FIELD( ir, a3, _IFT_NodalLoad_a3 );
    IR_GIVE_OPTIONAL_FIELD( ir, b4, _IFT_NodalLoad_b4 );
}


void NodalLoad :: giveInputRecord(DynamicInputRecord &input)
{
    Load :: giveInputRecord(input);
    input.setField(this->coordSystemType, _IFT_NodalLoad_cstype);
}


void
NodalLoad :: saveContext(DataStream &stream, ContextMode mode)
{
    Load :: saveContext(stream, mode);

    if ( mode & CM_Definition ) {
        if ( !stream.write(coordSystemType) ) {
          THROW_CIOERR(CIO_IOERR);
        }
    }
}


void
NodalLoad :: restoreContext(DataStream &stream, ContextMode mode)
{
    Load :: restoreContext(stream, mode);

    if ( mode & CM_Definition ) {
        int _val;
        if ( !stream.read(_val) ) {
          THROW_CIOERR(CIO_IOERR);
        }
        coordSystemType = (CoordSystType) _val;
    }
}


} // end namespace oofem
