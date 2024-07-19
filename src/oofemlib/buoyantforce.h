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

#ifndef buoyantforce_h
#define buoyantforce_h

#include "bodyload.h"
#include "bcgeomtype.h"
#include "valuemodetype.h"

#define _IFT_BuoyantForce_Name "buoyantforce"

namespace oofem {
/**
 * This class implements a body load acting on a body immersed in a fluid.
 * The inherited attribute 'componentArray' contains components of the associated
 * acceleration vector, which gives the direction of the buoyant force (it should in general
 * act in the direction opposite to the gravity force).
 * Since magnitude of the buoyant force is equivalent to the weight of the fluid displaced by
 * the body, the FE representing the body must have a fluid material associated with
 * it (e.g. see the DecoupledFluidMaterial class).
 */
class OOFEM_EXPORT BuoyantForce : public BodyLoad
{
public:
    /// Constructor
    BuoyantForce(int i, Domain * d) : BodyLoad(i, d) { }
    /**
     * Computes components values of the buoyant force vector at a given point (coordinates given in Global c.s.).
     * taking into account corresponding load time function value while respecting load response mode.
     * @param answer Component values at given point and time.
     * @param tStep Time step.
     * @param coords Global coordinates, which are used to evaluate components values.
     * @param mode Determines response mode-
     */
    void computeValueAt(FloatArray &answer, TimeStep *tStep, const FloatArray &coords, ValueModeType mode) override;

    bcValType giveBCValType() const override { return ForceLoadBVT; }
    bcGeomType giveBCGeoType() const override { return BodyLoadBGT; }

    void setBuoyantForceComponents(FloatArray newComponents);

    BodyLoadType giveBodyLoadType() override { return BodyLoadType::BuoyantForce; }
    const char *giveClassName() const override { return "BuoyantForce"; }
    const char *giveInputRecordName() const override { return _IFT_BuoyantForce_Name; }
};
} // end namespace oofem
#endif // buoyantforce_h