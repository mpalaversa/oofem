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

#ifndef truss3dnl_h
#define truss3dnl_h

#include "../sm/Elements/Bars/truss3d.h"

#define _IFT_Truss3dnl_Name "truss3dnl"
#define _IFT_Truss3dnl_initialStretch "initstretch"

namespace oofem {
class DecoupledMaterial;
/**
 * This class implements a nonlinear two-node truss bar element for three-dimensional
 * analysis.
 */
class Truss3dnl : public Truss3d
{
protected:
    double initialStretch;

    FloatArray viscousForce;

public:
    Truss3dnl(int n, Domain * d);
    virtual ~Truss3dnl() { }

    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_Truss3dnl_Name; }
    const char *giveClassName() const override { return "Truss3dnl"; }

    void initializeFrom(InputRecord &ir) override;
    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override;
    double computeVolumeAround( GaussPoint *gp ) override;
    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0) override;

    FloatArray giveViscousForce() override { return this->viscousForce; }

protected:

    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, TimeStep *tStep, bool lin = false);
    void computeBlMatrixAt(GaussPoint *gp, FloatMatrix &answer);
    void computeBnlMatrixAt(GaussPoint *gp, FloatMatrix &answer, TimeStep *tStep, bool lin = false);
    void computeInitialStressStiffness(FloatMatrix &answer, GaussPoint *gp, TimeStep *tStep);
    void computeHydrodynamicLoadVector( FloatArray &answer, FloatArray velocity, TimeStep *tStep ) override;

    // Returns characteristic dimension of the cross-section used in calculating hydrodynamic loads on the element.
    virtual double giveCharacteristicHydrodynamicDimension();
    // Returns characteristic dimension of the cross-section used in calculating weight, buoyant, inertia and added-mass force on the element.
    virtual double giveCharacteristicWeightDimension();
};
} // end namespace oofem
#endif // truss3dnl_h
