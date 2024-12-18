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

#ifndef consistentnetelement_h
#define consistentnetelement_h

#include "../sm/Elements/Bars/truss3dnl.h"

#define _IFT_ConsistentNetElement_Name "consistentnetelement"
#define _IFT_ConsistentNetElement_gf "gf"
#define _IFT_ConsistentNetElement_l0 "l0"

namespace oofem {
class DecoupledMaterial;
/**
 * This class implements a nonlinear two-node consistent net element for three-dimensional
 * analysis.
 */
class ConsistentNetElement : public Truss3dnl
{
protected:
    // Undeformed length of a twine and the globalisation factor
    double l0, gf;
    FloatArray viscousForce;

    void computeBodyLoadVectorAt( FloatArray &answer, Load *load, TimeStep *tStep, ValueModeType mode ) override;
    double giveCharacteristicHydrodynamicDimension() override;
    double giveCharacteristicWeightDimension() override;

public:
    ConsistentNetElement( int n, Domain *d );
    virtual ~ConsistentNetElement() {}
    const char *giveInputRecordName() const override { return _IFT_ConsistentNetElement_Name; }
    const char *giveClassName() const override { return "ConsistentNetElement"; }
    void initializeFrom( InputRecord &ir ) override;
};
} // end namespace oofem
#endif // consistentnetelement_h