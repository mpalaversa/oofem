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
 
#ifndef qdelement_h
#define qdelement_h

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
	class QdElement : public NLStructuralElement
	{
	protected:
        // Enumeration for output location of strains and stresses in element's x-y plane.
        enum class OutputLocationXY {
            GaussPoints,
            Centre,
            Corners,
            All
        };

        // Enumeration for output type of strains and stresses.
        enum class OutputType {
            Standard,
            Principal,
            VM,
            All
        };

        FEICellGeometry* cellGeometryWrapper;
        /**
        * Transformation Matrix form GtoL(3,3) is stored
        * at the element level for computation efficiency
        */
        FloatMatrix* GtoLRotationMatrix;
        // Local vertex coordinates
        std::vector< FloatArray > localCoords;

        OutputLocationXY outputAtXY;
		OutputType outputType;

        virtual void computeBmatrixAt(double xi, double eta, FloatMatrix& answer) = 0;
        void computeBmatrixAt(GaussPoint* gp, FloatMatrix& answer, int lowerIndx = 1, int upperIndx = ALL_STRAINS) override;
        const FloatMatrix* computeGtoLRotationMatrix();
        virtual void computeStrainVectorAt(FloatArray& answer, double xi, double eta, TimeStep* tStep) = 0;

    public:
        QdElement(int n, Domain* d);
        virtual ~QdElement() = default;

        int computeLoadGToLRotationMtrx(FloatMatrix& answer);
        int computeLoadLSToLRotationMatrix(FloatMatrix& answer, int iSurf, GaussPoint* gp) override { return 0; }
        void computeLocalNodalCoordinates(std::vector< FloatArray >& lxy);
        int computeNumberOfDofs() override;
        void computeStrainVector(FloatArray& answer, GaussPoint* gp, TimeStep* tStep) override;
        // When we want to get membrane or plate stresses at the output, the corresponding strains must be supplied. For shells, a method for strain calculation is invoked from within the method.
        virtual void computeStressVectorAtCentre(FloatArray& answer, TimeStep* tStep, const FloatArray& strain = 0) = 0;
        void computeSurfaceNMatrix(FloatMatrix& answer, int boundaryID, const FloatArray& lcoords) override;
        double computeSurfaceVolumeAround(GaussPoint* gp, int iSurf) override;
        OutputLocationXY getOutputLocationInXYPlane() { return outputAtXY; }
        OutputType getOutputType() { return outputType; }
        FEICellGeometry* giveCellGeometryWrapper();
        IntegrationRule* giveIntegrationRule(int i) override { return integrationRulesArray[i].get(); }
        std::vector< FloatArray > giveNodeCoordinates();
        void updateInternalState(TimeStep* tStep) override;
	};
}
#endif