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

#include "../sm/Elements/NetElements/consistentnetelement.h"
#include "../sm/CrossSections/structuralcrosssection.h"
#include "../sm/CrossSections/decoupledcrosssection.h"
#include "../sm/Materials/structuralms.h"
#include "bodyload.h"
#include "fei3dlinelin.h"
#include "node.h"
#include "material.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "load.h"
#include "mathfem.h"
#include "classfactory.h"
#include <math.h>


namespace oofem {
REGISTER_Element(ConsistentNetElement);

ConsistentNetElement :: ConsistentNetElement(int n, Domain *aDomain) : Truss3dnl(n, aDomain)
{
    l0 = -1;
    gf = -1;
    viscousForce.resize( 0 );
}

void
ConsistentNetElement ::computeBodyLoadVectorAt( FloatArray &answer, Load *forLoad, TimeStep *tStep, ValueModeType mode )
// Computes numerically the load vector of the receiver due to the body
// loads, at tStep.
// load is assumed to be in global cs.
// load vector is then transformed to coordinate system in each node.
// (should be global coordinate system, but there may be defined
//  different coordinate system in each node)
{
    double dens, dV, detJ, weight;
    FloatArray force, ntf;
    FloatMatrix n, T;

    if ( ( forLoad->giveBCGeoType() != BodyLoadBGT ) || ( forLoad->giveBCValType() != ForceLoadBVT ) ) {
        OOFEM_ERROR( "unknown load type" );
    }

    // note: force is assumed to be in global coordinate system.
    forLoad->computeComponentArrayAt( force, tStep, mode );
    // transform from global to element local c.s
    if ( this->computeLoadGToLRotationMtrx( T ) ) {
        force.rotatedWith( T, 'n' );
    }

    answer.clear();
    BodyLoad *load = dynamic_cast<BodyLoad *>( forLoad );
    if ( force.giveSize() ) {
        for ( GaussPoint *gp : *this->giveDefaultIntegrationRulePtr() ) {
            this->computeNmatrixAt( gp->giveSubPatchCoordinates(), n );
            detJ   = this->interp.giveTransformationJacobian( gp->giveNaturalCoordinates(), FEIElementGeometryWrapper( this ) );
            weight = gp->giveWeight();
            double area = this->giveCrossSection()->give( CS_Area, gp );

            if ( load->giveBodyLoadType() == BodyLoad::BodyLoadType::BuoyantForce ) {
                DecoupledCrossSection *cs = this->giveDecoupledCrossSectionOfType( DecoupledMaterial::DecoupledMaterialType::DecoupledFluidMaterial );
                dens                      = cs->giveMagnitudeOfMaterialProperty( 'd' );
            } else
                dens = this->giveCrossSection()->give( 'd', gp );

            ntf.beTProductOf( n, force );
            answer.add( detJ * weight * area * dens, ntf );
        }
    } else {
        return;
    }
}

double ConsistentNetElement::giveCharacteristicHydrodynamicDimension()
{
    DecoupledCrossSection *cs = this->giveDecoupledCrossSectionOfType( DecoupledMaterial::DecoupledMaterialType::DecoupledFluidMaterial );
    return gf * cs->giveCharacteristicDimension();

}

double ConsistentNetElement::giveCharacteristicWeightDimension() {
    DecoupledCrossSection *cs = this->giveDecoupledCrossSectionOfType( DecoupledMaterial::DecoupledMaterialType::DecoupledFluidMaterial );
    return cs->giveCharacteristicDimension() * sqrt( gf );
}


void
ConsistentNetElement :: initializeFrom(InputRecord &ir)
{
  Truss3dnl :: initializeFrom(ir);
  IR_GIVE_OPTIONAL_FIELD( ir, gf, _IFT_ConsistentNetElement_gf );
  IR_GIVE_OPTIONAL_FIELD( ir, l0, _IFT_ConsistentNetElement_l0 );
  if ( ( l0 > 0 && gf < 0 ) || ( gf > 0 && l0 == -1 ) )
      OOFEM_ERROR( "Both the globalization factor(s) and the undeformed twine lenght must be defined for element %d.", this->giveNumber() );
}

} // end namespace oofem