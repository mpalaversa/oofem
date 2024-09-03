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

#include "sm/CrossSections/decoupledcrosssection.h"
#include "classfactory.h"

namespace oofem {
REGISTER_CrossSection( DecoupledCrossSection );

void DecoupledCrossSection ::initializeFrom( InputRecord &ir )
{
    CrossSection ::initializeFrom( ir );

    characteristicDim = 0.0;
    IR_GIVE_FIELD( ir, characteristicDim, _DecoupledCrossSection_characteristicdim );

    materialNumber = -1;
    IR_GIVE_FIELD( ir, materialNumber, _DecoupledCrossSection_material );

    IR_GIVE_OPTIONAL_FIELD( ir, addedMassCoeff, _DecoupledCrossSection_addedmasscoeff );

    userDefinedDragCoeff = 0.0;
    IR_GIVE_OPTIONAL_FIELD( ir, userDefinedDragCoeff, _DecoupledCrossSection_dragcoeff );

    sn = 0.0;
    IR_GIVE_OPTIONAL_FIELD( ir, sn, _DecoupledCrossSection_sn );
}

bool DecoupledCrossSection::isDecoupled() {
    return true;
}

Material* DecoupledCrossSection::giveMaterial()
{
    return this->giveDomain()->giveMaterial( this->giveMaterialNumber() );
}

double DecoupledCrossSection::giveMagnitudeOfMaterialProperty(int property) {
    return this->giveMaterial()->give( property );
}
} // end namespace oofem
