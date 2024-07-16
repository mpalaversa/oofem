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

#ifndef decoupledmaterial_h
#define decoupledmaterial_h

#include "material.h"

///@name Input fields for DecoupledMaterial
//@{
#define _IFT_DecoupledMaterial_type "type"
//@}

namespace oofem {
/**
 * An abstract class representing decoupled materials used to associate secondary material properties with a FE (the primary
 * being those pertinent to the basic abstraction the FE represents). For example, a structural FE has the primary material
 * associated with it that provides material properties used in a structural analysis such as the Young's modulus,
 * density, etc. The secondary properties are associated with the FE and are used in, for example, a fluid flow
 * analysis based on a simplified model defined as loads in the structural analysis (e.g. see computeHydrodynamicLoadVector
 * method of StructuralElement class).
 * 
 * @author Marin Palaversa
 */
class OOFEM_EXPORT DecoupledMaterial : public Material
{
public:
    /// <summary>
    ///  List of decoupled material types.
    /// </summary>
    enum DecoupledMaterialType {
        Unknown,
        DecoupledFluidMaterial
    };

    /**
     * Constructor. Creates a material with given number, belonging to given domain.
     * @param n Material number.
     * @param d Domain to which new material will belong.
     */
    DecoupledMaterial(int n, Domain *d);

    const char *giveClassName() const override { return "DecoupledMaterial"; }

    /**
     * Returns a decoupled material type of the receiver.
     */
    DecoupledMaterialType giveType() { return type; }

    void initializeFrom( InputRecord &ir ) override;

    friend class DecoupledCrossSection;

protected:
    /// Type of the decoupled material (e.g. 1 denotes a DecoupledFluidMaterial). This should be replaced by an enumerator.
    DecoupledMaterialType type;
};
} // end namespace oofem
#endif // decoupledmaterial_h