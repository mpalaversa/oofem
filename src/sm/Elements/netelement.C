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
#include "sm/Elements/netelement.h"
#include "../sm/CrossSections/decoupledcrosssection.h"
#include "../sm/CrossSections/structuralcrosssection.h"
#include "gaussintegrationrule.h"

#include "bodyload.h"
#include "boundaryload.h"
#include "classfactory.h"
#include "gausspoint.h"
#include "load.h"

namespace oofem {
NetElement::NetElement(int n, Domain* aDomain) : NLStructuralElement(n, aDomain)
{
    sr                  = 1;    
    GtoLRotationMatrix = NULL;
    cellGeometryWrapper = NULL;
}

void
NetElement::computeBmatrixAt(GaussPoint* gp, FloatMatrix& answer, int lowerIndx, int upperIndx) {
    computeBmatrixAt(gp->giveNaturalCoordinate(1), gp->giveNaturalCoordinate(2), answer);
}

void
NetElement::computeBoundarySurfaceLoadVector( FloatArray &answer, BoundaryLoad *boundaryLoad, int boundary, CharType type, ValueModeType mode, TimeStep *tStep, bool global )
{
    answer.clear();
    if ( type != ExternalForcesVector ) {
        return;
    }
    FloatArray force;
    if ( boundaryLoad->giveType() == bcType::HydrodynamicMorison || boundaryLoad->giveType() == bcType::HydrodynamicKF ) {
        Load *load = dynamic_cast<BoundaryLoad *>( boundaryLoad );
        load->computeComponentArrayAt( force, tStep, mode );
        this->computeHydrodynamicLoadVector( answer, force, tStep );
    } else
        OOFEM_ERROR( "Only hydrodynamic surface load can be defined for element %d at this point.", giveGlobalNumber() );
}

void
NetElement ::computeConstitutiveMatrixAt( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep )
{
    answer = this->giveStructuralCrossSection()->giveStiffnessMatrix_1d( rMode, gp, tStep );
}

FloatMatrix*
NetElement::computeBasicGtoLRotationMatrix( TimeStep *tStep )
// Returns the rotation matrix of the receiver of the size [3,3]
// coords(local) = T * coords(global)
//
// OOFEM's local coordinate system (described by vector triplet e1',e2',e3') is defined as follows:
// e1'    : [N2-N1]    Ni - means i - th node
// v31   : [N3-N1]
// e3'    : e1' x v31
// e2'    : e3' x e1'
{
    FloatArray node1 = this->giveNode( 1 )->giveCoordinates();
    FloatArray node2 = this->giveNode( 2 )->giveCoordinates();
    FloatArray node3 = this->giveNode( 3 )->giveCoordinates();

    if ( tStep != 0 ) {
        // Fetch total displacements in the current configuration
        FloatArray u;
        u.resize( 9 );
        if ( !tStep->isTheFirstStep() )
            this->computeVectorOf( VM_Total, tStep, u );
        // Current coordinates of node 1
        double x = node1.at( 1 ) + u.at( 1 );
        double y = node1.at( 2 ) + u.at( 2 );
        double z = node1.at( 3 ) + u.at( 3 );
        node1.at( 1 ) = x;
        node1.at( 2 ) = y;
        node1.at( 3 ) = z;
        // Current coordinates of node 2
        x = node2.at( 1 ) + u.at( 4 );
        y = node2.at( 2 ) + u.at( 5 );
        z = node2.at( 3 ) + u.at( 6 );
        node2.at( 1 ) = x;
        node2.at( 2 ) = y;
        node2.at( 3 ) = z;
        // Current coordinates of node 3
        x = node3.at( 1 ) + u.at( 7 );
        y = node3.at( 2 ) + u.at( 8 );
        z = node3.at( 3 ) + u.at( 9 );
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

    GtoLRotationMatrix = new FloatMatrix( 3, 3 );
    for ( int i = 1; i <= 3; i++ ) {
        GtoLRotationMatrix->at( 1, i ) = e1.at( i );
        GtoLRotationMatrix->at( 2, i ) = e2.at( i );
        GtoLRotationMatrix->at( 3, i ) = e3.at( i );
    }

    return GtoLRotationMatrix;
}

void NetElement ::computeBodyLoadVectorAt( FloatArray &answer, Load *forLoad, TimeStep *tStep, ValueModeType mode )
// Computes numerically the load vector of the receiver due to the body
// loads, at tStep.
// load is assumed to be in global cs.
// load vector is then transformed to coordinate system in each node.
// (should be global coordinate system, but there may be defined
//  different coordinate system in each node)
{
    double dens, dV;
    FloatArray force;

    if ( ( forLoad->giveBCGeoType() != BodyLoadBGT ) || ( forLoad->giveBCValType() != ForceLoadBVT ) ) {
        OOFEM_ERROR( "unknown load type" );
    }

    // note: force is assumed to be in global coordinate system.
    forLoad->computeComponentArrayAt( force, tStep, mode );

    FloatArray bodyLoad;
    bodyLoad.resize( 3 );
    BodyLoad *load = dynamic_cast<BodyLoad *>( forLoad );
    if ( force.giveSize() ) {
        for ( GaussPoint *gp : *this->giveDefaultIntegrationRulePtr() ) {
            dV = this->computeVolumeAround( gp );

            if ( load->giveBodyLoadType() == BodyLoad::BodyLoadType::BuoyantForce ) {
                DecoupledCrossSection *cs = this->giveDecoupledCrossSectionOfType( DecoupledMaterial::DecoupledMaterialType::DecoupledFluidMaterial );
                dens                      = cs->giveMagnitudeOfMaterialProperty( 'd' );
            } else
                dens = this->giveCrossSection()->give( 'd', gp );

            bodyLoad.add( dV * dens, force );
        }

        this->calculateEquivalentLumpedNodalValues( answer, bodyLoad );
    } else {
        return;
    }
}

void
NetElement::computeLocalNodalCoordinates(std::vector< FloatArray >& lxy)
// Returns nodal coordinates in element's (local) coordinate system
{
    if (GtoLRotationMatrix == NULL) {
        this->computeBasicGtoLRotationMatrix();
    }

    lxy.resize(4);
    for (int i = 0; i < 4; i++) {
        const auto& nc = this->giveNode(i + 1)->giveCoordinates();
        lxy[i].beProductOf(*GtoLRotationMatrix, nc);
    }
}

double
NetElement::computeVolumeAround( GaussPoint *gp )
{
    double noOfTwines  = this->giveNumberOfTwines();
    double twineLength = this->giveTwineLength();
    double At          = this->giveCrossSection()->give( CS_Area, gp );

    return noOfTwines * At * twineLength;
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

void
NetElement::updateInternalState(TimeStep* tStep)
// Updates receiver at the end of the step.
{
    FloatArray stress, strain;
    // force updating strains & stresses
}
}