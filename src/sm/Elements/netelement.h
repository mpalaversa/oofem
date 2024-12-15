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
 
#ifndef netelement_h
#define netelement_h

#include "error.h"
#include "feinterpol2d.h"
#include "femcmpnn.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "sm/Elements/nlstructuralelement.h"

#include <cstdio>
#include <vector>
#include <memory>

namespace oofem {
	class NetElement : public NLStructuralElement
	{
	protected:
        FEICellGeometry* cellGeometryWrapper;
        FloatMatrix* GtoLRotationMatrix;
        // Local vertex coordinates
        std::vector< FloatArray > localCoords;

        virtual FloatArray calculateCurrentUnitNormalToElement( TimeStep *tStep ) = 0;
        virtual FloatArray calculateRelativeVelocity( FloatArray velocity, TimeStep *tStep ) = 0;
        virtual FloatArray calculateRelativeAcceleration( FloatArray acceleration, TimeStep *tStep ) = 0;
        virtual void computeBmatrixAt(double xi, double eta, FloatMatrix& answer) = 0;
        void computeBmatrixAt(GaussPoint* gp, FloatMatrix& answer, int lowerIndx = 1, int upperIndx = ALL_STRAINS) = 0;
        FloatMatrix *computeBasicGtoLRotationMatrix( TimeStep *tStep = 0 );
        void computeHydrodynamicLoadVector( FloatArray &answer, FloatArray flowCharacteristics, TimeStep *tStep ) override;
        virtual void calculateEquivalentLumpedNodalValues( FloatArray &answer, FloatArray vector ) = 0;
        virtual double giveTwineLength() = 0;
        virtual double giveNumberOfTwines() = 0;

    public:
        NetElement(int n, Domain* d);
        virtual ~NetElement() = default;

        void computeBodyLoadVectorAt( FloatArray &answer, Load *load, TimeStep *tStep, ValueModeType mode ) override;
        // Calculates nodal loads equivalent to the given surface loads on the element. At this moment, the method is implemented for hydrdynamic load only.
        void computeBoundarySurfaceLoadVector( FloatArray &answer, BoundaryLoad *load, int boundary, CharType type, ValueModeType mode, TimeStep *tStep, bool global = true ) override;
        void computeConstitutiveMatrixAt( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep ) override;
        int computeLoadLSToLRotationMatrix(FloatMatrix& answer, int iSurf, GaussPoint* gp) override { return 0; }
        void computeLocalNodalCoordinates(std::vector< FloatArray >& lxy);
        void computeMassMatrix( FloatMatrix &answer, TimeStep *tStep ) override { this->computeLumpedMassMatrix( answer, tStep ); }
        double computeVolumeAround( GaussPoint *gp ) override;
        MaterialMode giveMaterialMode() override { return _PlaneStress; }
        FEICellGeometry* giveCellGeometryWrapper();
        IntegrationRule* giveIntegrationRule(int i) override { return integrationRulesArray[i].get(); }
        std::vector< FloatArray > giveNodeCoordinates();
        void updateInternalState(TimeStep* tStep) override;
	};
}
#endif