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
class FEI3dLineLin;
	class NetElement : public NLStructuralElement
	{
	protected:
        // These elements don't use numerical integration. This is used just to access data associated with an integration point.
        static FEI3dLineLin interp;
        FEICellGeometry* cellGeometryWrapper;
        /**
        * Transformation Matrix form GtoL(3,3) is stored
        * at the element level for computation efficiency
        */
        FloatMatrix* GtoLRotationMatrix;
        // Local vertex coordinates
        std::vector< FloatArray > localCoords;

        virtual void computeBmatrixAt(double xi, double eta, FloatMatrix& answer) = 0;
        void computeBmatrixAt(GaussPoint* gp, FloatMatrix& answer, int lowerIndx = 1, int upperIndx = ALL_STRAINS) = 0;
        const FloatMatrix* computeGtoLRotationMatrix();

        // Numerical integration is not used in these elements. This is used to generate 1 GP to be used only in manipulating the associated element materials and cross-sections.
        void computeGaussPoints() override;
        FEInterpolation *giveInterpolation() const override;

    public:
        NetElement(int n, Domain* d);
        virtual ~NetElement() = default;

        void computeConstitutiveMatrixAt( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep ) override;
        double computeLength() override;
        int computeLoadGToLRotationMtrx(FloatMatrix& answer);
        int computeLoadLSToLRotationMatrix(FloatMatrix& answer, int iSurf, GaussPoint* gp) override { return 0; }
        void computeLocalNodalCoordinates(std::vector< FloatArray >& lxy);
        void computeMassMatrix( FloatMatrix &answer, TimeStep *tStep ) override { this->computeLumpedMassMatrix( answer, tStep ); }
        void computeSurfaceNMatrix(FloatMatrix& answer, int boundaryID, const FloatArray& lcoords) override;
        MaterialMode giveMaterialMode() override { return _PlaneStress; }
        double computeSurfaceVolumeAround(GaussPoint* gp, int iSurf) override;
        FEICellGeometry* giveCellGeometryWrapper();
        IntegrationRule* giveIntegrationRule(int i) override { return integrationRulesArray[i].get(); }
        std::vector< FloatArray > giveNodeCoordinates();
        void updateInternalState(TimeStep* tStep) override;
	};
}
#endif