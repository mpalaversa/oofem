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
 
#ifndef netrec4trla_h
#define netrec4trla_h

#include "error.h"
#include "feinterpol2d.h"
#include "femcmpnn.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "sm/Elements/netelement.h"

#include <cstdio>
#include <vector>
#include <memory>

#define _IFT_NetRec4TrLa_Name "netrec4trla"
#define _IFT_NetRec4TrLa_Mask "mask"
#define _IFT_NetRec4TrLa_dx0 "dx0"
#define _IFT_NetRec4TrLa_dy0 "dy0"
#define _IFT_NetRec4TrLa_nx "nx"
#define _IFT_NetRec4TrLa_ny "ny"

namespace oofem {
	class NetRec4TrLa : public NetElement
	{
        protected:
            int mask;
            // Mesh half length (for square meshes), mesh diagonals (for diamond meshes), no. of meshes along element side in x- and y-direction respectively
            double a0, dx0, dy0, nx, ny;
            FloatArray initialDimensions;
            /**
             * Calculates element's dimensions in the current configuration.
             * Assumes rectangular shape of the element.
             * Returns array having element's dimension in the local x-direction at position 1 and element's
             * dimension in the local y-direction at position 2.
             */
            FloatArray calculateElementDimensions( TimeStep *tStep );
            FloatArray calculateInternalDisplacements( TimeStep *tStep );
            //FloatMatrix giveShapeFunctionAtXY( double x, double y );

        public:
            NetRec4TrLa(int n, Domain* d);
            virtual ~NetRec4TrLa() = default;

		    bool computeGtoLRotationMatrix( FloatMatrix &answer ) override;
            void computeLumpedMassMatrix( FloatMatrix &answer, TimeStep *tStep ) override;
            int computeNumberOfDofs() override { return 8; }
            int computeNumberOfGlobalDofs() override { return 12; }
            void computeStiffnessMatrix( FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep ) override;
            void computeStrainVector( FloatArray &answer, GaussPoint *gp, TimeStep *tStep ) override;
            void computeStressVector( FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep ) override;
            void giveDofManDofIDMask( int inode, IntArray & ) const override;
            const char *giveClassName() const override { return "NetRec4TrLa"; }
            const char *giveInputRecordName() const override { return _IFT_NetRec4TrLa_Name; }
            void giveInternalForcesVector( FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0 ) override;
            void initializeFrom( InputRecord &ir ) override;

            void computeBmatrixAt( GaussPoint *gp, FloatMatrix &answer, int lowerIndx = 1, int upperIndx = ALL_STRAINS ) override{};
            void computeConsistentMassMatrix( FloatMatrix &answer, TimeStep *tStep, double &mass, const double *ipDensity = NULL ) override {};
            void computeBmatrixAt( double xi, double eta, FloatMatrix &answer ) override{};
	};
}
#endif