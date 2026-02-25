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
    speedUp = false;
}


void
Truss3dnl :: initializeFrom(InputRecord &ir)
{
  Truss3d :: initializeFrom(ir);
  initialStretch = 1;
  IR_GIVE_OPTIONAL_FIELD(ir, initialStretch, _IFT_Truss3dnl_initialStretch);

  int sup = 0;
  IR_GIVE_OPTIONAL_FIELD( ir, sup, _IFT_ConsistentNetElement_sup );
  if ( sup == 1 )
      speedUp = true;
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
    double area   = this->giveCrossSection()->give( CS_Area, gp );

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

double Truss3dnl::giveCharacteristicHydrodynamicDimension() {
    DecoupledCrossSection *cs = this->giveDecoupledCrossSectionOfType( DecoupledMaterial::DecoupledMaterialType::DecoupledFluidMaterial );
    return cs->giveCharacteristicDimension();
}

double Truss3dnl::giveCharacteristicWeightDimension()
{
    return giveCharacteristicHydrodynamicDimension();
}

void Truss3dnl ::computeHydrodynamicLoadMorison( FloatArray &answer, FloatArray flowCharacteristics, TimeStep *tStep, bool knotted )
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

    // Account for the local speed-up if needed
    if ( speedUp ) {
        double sn = cs->giveSolidityRatio();
        for ( int i = 1; i <= 3; i++ )
            relativeVelocity.at( i ) = relativeVelocity.at( i ) / ( 1 - sn );
    }

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
    double characteristicDim  = this->giveCharacteristicHydrodynamicDimension();

    // Get drag coefficients
    FloatArray dragCoeffs;
    if ( userDefinedDragCoeff == 0.0 ) {
        // If no, calculate the drag coefficients.
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
    //if ( this->giveNumber() == 145 || this->giveNumber() == 146 || this->giveNumber() == 147 || this->giveNumber() == 148 || this->giveNumber() == 149 || this->giveNumber() == 150 || this->giveNumber() == 151 || this->giveNumber() == 152 ) {
        normalViscousForce.beScaled( 0.5 * density * dragCoeffs.at( 1 ) * characteristicDim * l * normalRelativeVelocity.computeNorm(), normalRelativeVelocity );
   // } else {
        //normalViscousForce.beScaled( density * dragCoeffs.at( 1 ) * characteristicDim * l * normalRelativeVelocity.computeNorm(), normalRelativeVelocity );
   // }
    tangentialViscousForce.beScaled( dragCoeffs.at( 2 ) * l, tangentialRelativeVelocity );
    this->viscousForce.zero();
    this->viscousForce.add( normalViscousForce );
    this->viscousForce.add( tangentialViscousForce );

    // Check if the net is knotted and add contribution of the knot
    FloatArray dragForceOnKnots;
    if ( knotted ) {
        computeDragForceOnKnots( dragForceOnKnots, density, relativeVelocity );
        this->viscousForce.add( dragForceOnKnots );
    }

    // The force is equally distributed among the element's nodes
    answer.resize( 6 );
    answer.at( 1 ) = answer.at( 4 ) = viscousForce.at( 1 ) / 2;
    answer.at( 2 ) = answer.at( 5 ) = viscousForce.at( 2 ) / 2;
    answer.at( 3 ) = answer.at( 6 ) = viscousForce.at( 3 ) / 2;

    // Form the fluid acceleration vector
    FloatArray acceleration;
    acceleration.resize( 3 );
    for ( int i = 4; i <= 6; i++ )
        //acceleration.at( i - 3 ) = 0;
        acceleration.at( i - 3 ) = flowCharacteristics.at( i );

    
    // Calculate acceleration of the fluid relative to the element
    FloatArray relativeAcceleration;
    relativeAcceleration.resize( 3 );
    relativeAcceleration.at( 1 ) = acceleration.at( 1 ) - ( currentNodalAcceleration.at( 1 ) + currentNodalAcceleration.at( 4 ) ) / 2;
    relativeAcceleration.at( 2 ) = acceleration.at( 2 ) - ( currentNodalAcceleration.at( 2 ) + currentNodalAcceleration.at( 5 ) ) / 2;
    relativeAcceleration.at( 3 ) = acceleration.at( 3 ) - ( currentNodalAcceleration.at( 3 ) + currentNodalAcceleration.at( 6 ) ) / 2;
    
    if ( relativeAcceleration.computeNorm() != 0 ) {
        
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

        characteristicDim = this->giveCharacteristicWeightDimension();
        
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

    //if ( this->giveNumber() == 148 || this->giveNumber() == 105 || this->giveNumber() == 36 || this->giveNumber() == 109 || this->giveNumber() == 112 )
        //OOFEM_LOG_RELEVANT( "Element %d. z: %e. Cdn: %e. Cdt: %e. u: %e. w: %e. ax: %e. az: %e \n", this->giveNumber(), (currentNode1Coordinates.at(3) + currentNode2Coordinates.at(3))/2, dragCoeffs.at(1), dragCoeffs.at(2), flowCharacteristics.at(1), flowCharacteristics.at(3), flowCharacteristics.at(4), flowCharacteristics.at(6) );
}

void
Truss3dnl::computeHydrodynamicLoadVector( FloatArray &answer, FloatArray loadInputData, bcType loadType, TimeStep *tStep )
{
    FloatArray currentLoadsMorison, waveLoadsStokes2;

    if ( loadType == bcType::HydrodynamicMorison ) {
        computeHydrodynamicLoadMorison( answer, loadInputData, tStep );
    } else if ( loadType == bcType::HydrodynamicWaveStokes2 ) {
        computeHydrodynamicLoadFromWavesStokes2( answer, loadInputData, tStep );
    } else
        OOFEM_ERROR( "The following hydrodynamic loads are implemented at the moment: current loads according to the Morison's equation, wave loads according to the Stokes 2nd-order wave theory." );
}

void
Truss3dnl ::computeHydrodynamicLoadFromWavesStokes2( FloatArray &answer, FloatArray waveCharacteristics, TimeStep *tStep, bool knotted )
{
    // User-defined wave height (H), wave period (T), wave direction (beta in deg), water depth (h) and the convergence criterion for wave number (kErr)
    double H    = waveCharacteristics.at( 1 );
    double T    = waveCharacteristics.at( 2 );
    double beta = waveCharacteristics.at( 3 ) * 3.14 / 180;
    double h    = waveCharacteristics.at( 4 );
    double kErr = waveCharacteristics.at( 5 );

    // Current time
    double t = tStep->giveTargetTime();

    // Determine k by means of the secant method
    /*
    FloatArray kTemp;
    kTemp.resize( 3 );
    kTemp.at( 1 ) = 4 * pow( 3.14, 2 ) / ( 9.81 * pow( T, 2 ) );
    kTemp.at( 2 ) = 1.1 * kTemp.at( 1 );
    int steps     = 3;
    double omega  = 2 * 3.14 / T;
    while ( ( abs( kTemp.at( 2 ) - kTemp.at( 1 ) ) > kErr ) && ( steps <= 100 ) ) {
        if ( steps > 3 ) {
            kTemp.at( 1 ) = kTemp.at( 2 );
            kTemp.at( 2 ) = kTemp.at( 3 );
        }
        kTemp.at( 3 ) = kTemp.at( 2 ) - ( kTemp.at( 2 ) - kTemp.at( 1 ) ) * ( pow( omega, 2 ) - kTemp.at( 2 ) * 9.81 * tanh( kTemp.at( 2 ) * h ) ) / ( ( pow( omega, 2 ) - kTemp.at( 2 ) * 9.81 * tanh( kTemp.at( 2 ) * h ) ) - ( pow( omega, 2 ) - kTemp.at( 1 ) * 9.81 * tanh( kTemp.at( 1 ) * h ) ) );
        steps++;
    }
    if ( steps > 100 )
        OOFEM_ERROR( "\n Wave number k was not found in 100 steps for element %d.", this->giveNumber() );
    double k = kTemp.at( 3 );
    */
    double k = 4 * pow( 3.14 , 2 ) / ( 9.81 * pow( T, 2 ) );
    double omega    = 2 * 3.14 / T;
    // To find the velocities, we need to distance from the origin (x) and depth (z) at which the current finite element is
    // Fetch initial coordinates of element's nodes in Oxyz
    FloatArray node1 = this->giveNode( 1 )->giveCoordinates();
    FloatArray node2 = this->giveNode( 2 )->giveCoordinates();
    // Fetch total displacements in the current configuration
    FloatArray u;
    u.resize( 12 );
    if ( !tStep->isTheFirstStep() )
        this->computeVectorOf( VM_Total, tStep, u );
    // Calculate the average x- and z-coordinate of the element
    double x = ( node1.at( 1 ) + u.at( 1 ) + node2.at( 1 ) + u.at( 4 ) ) / 2;
    double z = ( node1.at( 3 ) + u.at( 3 ) + node2.at( 3 ) + u.at( 6 ) ) / 2;

    FloatArray flowCharacteristics;
    flowCharacteristics.resize( 6 );
    // Fluid velocity in x and z direction
    //flowCharacteristics.at( 1 ) = H / 2 * ( 9.81 * k / omega ) * cosh( k * ( h + z ) ) / cosh( k * h ) * cos( k * x - cos( beta ) * omega * t ) + 3 / 16 * pow( H, 2 ) * omega * k * cosh( 2 * k * ( h + z ) ) / pow( sinh( k * h ), 4 ) * cos( 2 * ( k * x - cos( beta ) * omega * t ) );
    //flowCharacteristics.at( 3 ) = H / 2 * ( 9.81 * k / omega ) * sinh( k * ( h + z ) ) / cosh( k * h ) * sin( k * x - cos( beta ) * omega * t ) + 3 / 16 * pow( H, 2 ) * omega * k * sinh( 2 * k * ( h + z ) ) / pow( sinh( k * h ), 4 ) * sin( 2 * ( k * x - cos( beta ) * omega * t ) );
    flowCharacteristics.at( 1 ) = H / 2 * omega * cosh( k * ( h + z ) ) / sinh( k * h ) * cos( k * x - cos( beta ) * omega * t );
    flowCharacteristics.at( 3 ) = H / 2 * omega * sinh( k * ( h + z ) ) / sinh( k * h ) * sin( k * x - cos( beta ) * omega * t );

    // Fluid acceleration in x and z direction
    //flowCharacteristics.at( 4 ) = H / 2 * 9.81 * k * cosh( k * ( h + z ) ) / cosh( k * h ) * sin( k * x - cos( beta ) * omega * t ) - pow( H, 2 ) / 4 * 9.81 * pow( k, 2 ) * sin( 2 * ( k * x - cos( beta ) * omega * t ) ) / sinh( 2 * k * h ) + 3 / 8 * pow( H, 2 ) * pow( omega, 2 ) * k * cosh( 2 * k * ( h + z ) ) / pow( sinh( k * h ), 4 ) * sin( 2 * ( k * x - cos( beta ) * omega * t ) );
    //flowCharacteristics.at( 6 ) = -H / 2 * 9.81 * k * sinh( k * ( h + z ) ) / cosh( k * h ) * cos( k * x - cos( beta ) * omega * t ) + pow( H, 2 ) / 4 * 9.81 * pow( k, 2 ) * sinh( 2 * k * ( h + z ) ) / sinh( 2 * k * h ) - 3 / 8 * pow( H, 2 ) * pow( omega, 2 ) * k * sinh( 2 * k * ( h + z ) ) / pow( sinh( k * h ), 4 ) * sin( 2 * ( k * x - cos( beta ) * omega * t ) );
    flowCharacteristics.at( 4 ) = H / 2 * pow( omega, 2 ) * cosh( k * ( h + z ) ) / sinh( k * h ) * sin ( k * x - cos( beta ) * omega * t );
    flowCharacteristics.at( 6 ) = - H / 2 * pow( omega, 2 ) * sinh( k * ( h + z ) ) / sinh( k * h ) * cos( k * x - cos( beta ) * omega * t );


    computeHydrodynamicLoadMorison( answer, flowCharacteristics, tStep, knotted );
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

