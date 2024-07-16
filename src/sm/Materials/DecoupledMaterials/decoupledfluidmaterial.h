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

#ifndef decoupledfluidmaterial_h
#define decoupledfluidmaterial_h

#include "sm/Materials/DecoupledMaterials/decoupledmaterial.h"

///@name Input fields for DecoupledFluidMaterial
//@{
#define _IFT_DecoupledFluidMaterial_Name "decoupledfluidmaterial"
#define _IFT_DecoupledFluidMaterial_mu "mu"
//@}

namespace oofem {

/**
 * This is a decoupled material for fluids.
 *
 * @author Marin Palaversa
 */
class OOFEM_EXPORT DecoupledFluidMaterial : public DecoupledMaterial
{
public:

    /**
     * Constructor. Creates material with given number, belonging to given domain.
     * @param n Material number.
     * @param d Domain to which new material will belong.
     */
    DecoupledFluidMaterial(int n, Domain *d);

    const char *giveClassName() const override { return "DecoupledFluidMaterial"; }

    const char *giveInputRecordName() const override { return _IFT_DecoupledFluidMaterial_Name; }

    void initializeFrom( InputRecord &ir ) override;
    /*
    friend class CrossSection;
    friend class StructuralCrossSection;
    friend class SimpleCrossSection;
    friend class LayeredCrossSection;*/
};
} // end namespace oofem
#endif // decoupledfluidmaterial_h