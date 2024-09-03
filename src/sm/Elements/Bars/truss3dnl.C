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

#include "../sm/Elements/Bars/truss3dnl.h"
#include "../sm/CrossSections/structuralcrosssection.h"
#include "../sm/CrossSections/decoupledcrosssection.h"
#include "../sm/Materials/structuralms.h"
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
#include <math.h>


namespace oofem {
REGISTER_Element(Truss3dnl);

Truss3dnl :: Truss3dnl(int n, Domain *aDomain) : Truss3d(n, aDomain)
{
    
}


void
Truss3dnl :: initializeFrom(InputRecord &ir)
{
  Truss3d :: initializeFrom(ir);
  initialStretch = 1;
  IR_GIVE_OPTIONAL_FIELD(ir, initialStretch, _IFT_Truss3dnl_initialStretch);
  gf.resize( 1 );
  IR_GIVE_OPTIONAL_FIELD( ir, gf, _IFT_Truss3dnl_gf );
  l0 = -1;
  IR_GIVE_OPTIONAL_FIELD( ir, l0, _IFT_Truss3dnl_l0 );
  if ( ( l0 > 0 && gf.computeNorm() == 0 ) || ( gf.computeNorm() != 0 && l0 == -1 ) )
      OOFEM_ERROR( "Both the globalization factor(s) and the undeformed twine lenght must be defined for element %d.", this->giveNumber() );
}

  
void
Truss3dnl :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
  FloatMatrix B, Be;
  FloatArray vStress, vStrain, u;
  
  // This function can be quite costly to do inside the loops when one has many slave dofs.
  this->computeVectorOf(VM_Total, tStep, u);
  // subtract initial displacements, if defined
  if ( initialDisplacements ) {
    u.subtract(* initialDisplacements);
  }
  
  // zero answer will resize accordingly when adding first contribution
  answer.clear();
  
  for ( auto &gp: *this->giveDefaultIntegrationRulePtr() ) {
    StructuralMaterialStatus *matStat = static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() );
    this->computeBmatrixAt(gp, B, tStep, true);
    this->computeBmatrixAt(gp, Be, tStep);
    if ( useUpdatedGpRecord == 1 ) {
      vStress = matStat->giveStressVector();
    } else {
      ///@todo Is this really what we should do for inactive elements?
      if ( !this->isActivated(tStep) ) {
	vStrain.resize( StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() ) );
	vStrain.zero();
      }
      vStrain.beProductOf(Be, u);
      // add influence of initial stress/stretch
      double l2 = initialStretch*initialStretch;
      vStrain.times(l2);
      FloatArray E0(1);
      E0.at(1) = (l2-1.)/2.;
      vStrain.add(E0);
      //
      this->computeStressVector(vStress, vStrain, gp, tStep);
    }
    
    if ( vStress.giveSize() == 0 ) { /// @todo is this check really necessary?
      break;
    }
    
    // Compute nodal internal forces at nodes as f = B^T*Stress dV
    double dV  = this->computeVolumeAround(gp);
    
    if ( vStress.giveSize() == 6 ) {
      // It may happen that e.g. plane strain is computed
      // using the default 3D implementation. If so,
      // the stress needs to be reduced.
      // (Note that no reduction will take place if
      //  the simulation is actually 3D.)
      FloatArray stressTemp;
      StructuralMaterial :: giveReducedSymVectorForm( stressTemp, vStress, gp->giveMaterialMode() );
      answer.plusProduct(B, stressTemp, dV);
    } else   {
      answer.plusProduct(B, vStress, dV);
    }
    
    
    // If inactive: update fields but do not give any contribution to the internal forces
    if ( !this->isActivated(tStep) ) {
      answer.zero();
      return;
    }
  }
}
  
  
  
  
void
Truss3dnl :: computeStiffnessMatrix(FloatMatrix &answer,
				    MatResponseMode rMode, TimeStep *tStep)
{
  StructuralCrossSection *cs = this->giveStructuralCrossSection();
  bool matStiffSymmFlag = cs->isCharacteristicMtrxSymmetric(rMode);
  
  answer.clear();
  
  if ( !this->isActivated(tStep) ) {
    return;
  }
  
  // Compute matrix from material stiffness (total stiffness for small def.) - B^T * dS/dE * B
  if ( integrationRulesArray.size() == 1 ) {
    FloatMatrix B, D, DB, Ksigma;
    for ( auto &gp : *this->giveDefaultIntegrationRulePtr() ) {
      this->computeBmatrixAt(gp, B, tStep, true);
      this->computeConstitutiveMatrixAt(D, rMode, gp, tStep);
      double dV = this->computeVolumeAround(gp);
      DB.beProductOf(D, B);
      if ( matStiffSymmFlag ) {
	answer.plusProductSymmUpper(B, DB, dV);
      } else {
	answer.plusProductUnsym(B, DB, dV);
      }
      this->computeInitialStressStiffness(Ksigma, gp, tStep);
      Ksigma.times(dV);
      answer.add(Ksigma);
      
    }
    
    if ( matStiffSymmFlag ) {
      answer.symmetrized();
    }
  }
}
  
double
Truss3dnl ::computeVolumeAround( GaussPoint *gp )
// Returns the length of the receiver. This method is valid only if 1
// Gauss point is used.
{
    double detJ   = this->interp.giveTransformationJacobian( gp->giveNaturalCoordinates(), FEIElementGeometryWrapper( this ) );
    double weight = gp->giveWeight();
    double area   = 0.0;
    if ( gf.computeNorm() == 0 )
        area = this->giveCrossSection()->give( CS_Area, gp );
    else {
        // The following block of code is used for equivalent numerical twines
        double de = 0.0;
        if ( gf.giveSize() == 1 )
            de = this->giveDecoupledCrossSectionOfType( DecoupledMaterial::DecoupledMaterialType::DecoupledFluidMaterial )->giveCharacteristicDimension() * sqrt( gf.at( 1 ) * this->computeLength() / l0 );
        else
            de = this->giveDecoupledCrossSectionOfType( DecoupledMaterial::DecoupledMaterialType::DecoupledFluidMaterial )->giveCharacteristicDimension() * sqrt( gf.at( 1 ) * gf.at( 2 ) );
        // Calculate area of element's cross-section (assumes circular shape)
        area = pow( de, 2) * 3.14 / 4;
    }
    return detJ * weight * area;
}

void
Truss3dnl :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, TimeStep *tStep, bool lin)
{
  FloatMatrix Bl, Bnl;
  this->computeBlMatrixAt(gp, Bl);
  this->computeBnlMatrixAt(gp, Bnl, tStep, lin);
  answer = Bl;
  answer.add(Bnl);
}



void
Truss3dnl :: computeBlMatrixAt(GaussPoint *gp, FloatMatrix &answer)
//
// Returns linear part of geometrical equations of the receiver at gp.
// Returns the linear part of the B matrix
//
{
  Truss3d::computeBmatrixAt(gp, answer);
}



void
Truss3dnl :: computeBnlMatrixAt(GaussPoint *gp, FloatMatrix &answer, TimeStep *tStep, bool lin)
//
// Returns linear part of geometrical equations of the receiver at gp.
// Returns the linear part of the B matrix
//
{
  FloatArray d;
  this->computeVectorOf(VM_Total, tStep, d);
    
  FloatMatrix Bnl, A(6,6);
  A.at(1,1) = A.at(2,2) = A.at(3,3) = A.at(4,4) = A.at(5,5) = A.at(6,6) =  1.0;
  A.at(1,4) = A.at(2,5) = A.at(3,6) = A.at(4,1) = A.at(5,2) = A.at(6,3) = -1.0;
  double l0 = this->computeLength();
  double factor = 1/l0/l0;
  if(!lin) {
    factor /= 2;
  } 
  Bnl.beProductOf(A,d);
  Bnl.times(factor);
  answer.beTranspositionOf(Bnl);
  
}

void Truss3dnl ::computeHydrodynamicLoadVector( FloatArray &answer, FloatArray flowCharacteristics, TimeStep *tStep )
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

    DecoupledCrossSection *cs = this->giveDecoupledCrossSectionOfType( DecoupledMaterial::DecoupledMaterialType::DecoupledFluidMaterial );
    // Check if the element is downstream relative to another element
    if ( this->isDownstream ) {
        double sn = cs->giveSolidityRatio();
        // Reduce the inflow velocity by the velocity reduction coefficient as given in Loland, G. Current forces on and flow thorugh fish farms
        if ( sn > 0 )
            velocity.beScaled( 1 - 0.46 * ( 0.33 * sn + 6.54 * pow( sn, 2 ) - 4.88 * pow( sn, 3 ) ), velocity );
        else
            OOFEM_ERROR( "Element %d is denoted as downstream, but the solidty ratio is not specified.", this->giveNumber() );
    }

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
    double userDefinedDragCoeff = cs->giveDragCoefficient();
    double l                  = this->computeLength();
    double density            = cs->giveMagnitudeOfMaterialProperty( 'd' );
    double mu                 = cs->giveMaterial()->giveDynamicViscosity();
    double characteristicDim  = 0.0;
    if ( this->gf.computeNorm() == 0 )
        characteristicDim = cs->giveCharacteristicDimension();
    // If the globalization factor(s) are defined, the equivalent numerical twines are used and the characteristic dimension must be calculated
    else {
        if ( gf.giveSize() == 1 )
            characteristicDim = gf.at( 1 ) * cs->giveCharacteristicDimension();
        else
            characteristicDim = 2 * gf.at( 1 ) * gf.at( 2 ) * l0 * cs->giveCharacteristicDimension() / l;
    }

    // Caluclate drag coefficients
    FloatArray dragCoeffs;
    if ( userDefinedDragCoeff == 0.0 ) {
        // If no, calculate the drag coefficients.
        if ( this->gf.computeNorm() == 0 )
            dragCoeffs = computeDragCoefficients( density, mu, characteristicDim, normalRelativeVelocity.computeNorm() );
        else
            dragCoeffs = computeDragCoefficients( density, mu, cs->giveCharacteristicDimension(), normalRelativeVelocity.computeNorm() );
    }    
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

        if ( gf.computeNorm() != 0 ) {
            if ( gf.giveSize() == 1 )
                characteristicDim = cs->giveCharacteristicDimension() * sqrt( gf.at( 1 ) );
            else
                characteristicDim = cs->giveCharacteristicDimension() * sqrt( 2 * gf.at( 1 ) * gf.at( 2 ) * l0 / l );
        }
        
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
Truss3dnl :: computeInitialStressStiffness(FloatMatrix &answer, GaussPoint *gp, TimeStep *tStep)
{
    answer.resize(6,6);
    answer.at(1,1) = answer.at(2,2) = answer.at(3,3) = answer.at(4,4) = answer.at(5,5) = answer.at(6,6) =  1.0;
    answer.at(1,4) = answer.at(2,5) = answer.at(3,6) = answer.at(4,1) = answer.at(5,2) = answer.at(6,3) = -1.0;
    
    FloatArray d, strain;
    FloatMatrix B;
    this->computeVectorOf(VM_Total, tStep, d);
    this->computeBmatrixAt(gp, B, tStep);	  
    strain.beProductOf(B, d);
    // add influence of initial stress/stretch
    double l2 = initialStretch*initialStretch;
    strain.times(l2);
    FloatArray E0(1);
    E0.at(1) = (l2-1.)/2;
    strain.add(E0);
    /////////////////////////////////////////////////////////////////////////////////////////
    auto stress = this->giveStructuralCrossSection()->giveRealStress_1d(strain, gp, tStep);
    double l0 = this->computeLength();	
    double factor = 1/l0/l0;
    // prevent zero initial stress stiffness
    if ( stress.at(1) == 0 ) {
        stress.at(1) = 1;
    }
    answer.times(stress.at(1));
    answer.times(factor);
}
    

} // end namespace oofem

