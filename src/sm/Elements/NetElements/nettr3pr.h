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
 
#ifndef nettr3pr_h
#define nettr3pr_h

#include "error.h"
#include "feinterpol2d.h"
#include "femcmpnn.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "sm/Elements/netelement.h"

#include <cstdio>
#include <vector>
#include <memory>

#define _IFT_NetTr3Pr_Name "nettr3pr"
#define _IFT_NetTr3Pr_L0 "l0"
#define _IFT_NetTr3Pr_U1 "u1"
#define _IFT_NetTr3Pr_V1 "v1"
#define _IFT_NetTr3Pr_U2 "u2"
#define _IFT_NetTr3Pr_V2 "v2"
#define _IFT_NetTr3Pr_U3 "u3"
#define _IFT_NetTr3Pr_V3 "v3"
#define _IFT_NetTr3Pr_sr "sr"

namespace oofem {
class FEI2dTrLin;
	class NetTr3Pr : public NetElement
	{
        protected:
            // These elements don't use numerical integration. This is used just to access data associated with an integration point.
            static FEI2dTrLin interp;
            int mask;
            // Undeformed length of a twine
            double L0;
            // Nodes' coordinates in the twine coord. syst.
            double U1, V1, U2, V2, U3, V3;
            // An auxilliary variable (see Priour D. A Finite Element Method for Netting Application to Fish Cages and Fishing Gear)
            double d;

            FloatArray calculateRelativeVelocity( FloatArray velocity, TimeStep *tStep ) override;
            virtual FloatArray calculateRelativeAcceleration( FloatArray velocity, TimeStep *tStep ) override;
            void calculateEquivalentLumpedNodalValues( FloatArray &answer, FloatArray vector ) override;
            // Numerical integration is not used in these elements. This is used to generate 1 GP to be used only in manipulating the associated element materials and cross-sections.
            void computeGaussPoints() override;
            void computeHydrodynamicLoadVector( FloatArray &answer, FloatArray flowCharacteristics, TimeStep *tStep ) override;
            double giveTwineLength() override { return L0; };
            // Returns the total number of twines within the element
            double giveNumberOfTwines() override { return d; };

        public:
            NetTr3Pr(int n, Domain* d);
            virtual ~NetTr3Pr() = default;

            void computeLumpedMassMatrix( FloatMatrix &answer, TimeStep *tStep ) override;
            int computeNumberOfDofs() override { return 9; }
            int computeNumberOfGlobalDofs() override { return 9; }
            void computeStiffnessMatrix( FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep ) override;
            void computeStrainVector( FloatArray &answer, GaussPoint *gp, TimeStep *tStep ) override;
            void computeStressVector( FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep ) override;
            // Returns unit vector of the twine coord. syst. in U-direction expressed in the global coord. syst.
            FloatArray computeUTwine( TimeStep *tStep );
            // Returns a unit vector of the twine coord. syst. in V-direction expressed in the global coord. syst.
            FloatArray computeVTwine( TimeStep *tStep );
            std::unique_ptr<IntegrationRule> giveBoundarySurfaceIntegrationRule( int order, int boundary ) override;
            const char *giveClassName() const override { return "NetTr3Pr"; }
            void giveDofManDofIDMask( int inode, IntArray & ) const override;
            const char *giveInputRecordName() const override { return _IFT_NetTr3Pr_Name; }
            void giveInternalForcesVector( FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0 ) override;
            FEInterpolation *giveInterpolation() const override;
            void initializeFrom( InputRecord &ir ) override;

            void computeBmatrixAt( GaussPoint *gp, FloatMatrix &answer, int lowerIndx = 1, int upperIndx = ALL_STRAINS ) override{};
            void computeConsistentMassMatrix( FloatMatrix &answer, TimeStep *tStep, double &mass, const double *ipDensity = NULL ) override {};
            void computeBmatrixAt( double xi, double eta, FloatMatrix &answer ) override{};
	};
}
#endif