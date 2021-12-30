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

#include "sm/Elements/Beams/hybbeam.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Element(HybBeam);
HybBeam::HybBeam(int n, Domain* aDomain) : EccBeam(n, aDomain)
{
    numberOfDofMans = 2;
    numberOfGaussPoints = 1;
}

HybBeam :: ~HybBeam()
{

}

void
HybBeam::computeStiffnessMatrix(FloatMatrix& answer, MatResponseMode rMode, TimeStep* tStep)
{
    double L = this->computeLength();
    FloatMatrix D;
    this->computeConstitutiveMatrixAt(D, rMode, this->giveDefaultIntegrationRulePtr()->getIntegrationPoint(0), tStep);
    answer.resize(12, 12);
    answer.at(1, 1) = D.at(1, 1) / L;
    answer.at(1, 7) = -D.at(1, 1) / L;
    answer.at(2, 2) = 12 * D.at(4, 4) / (L * L * L);
    answer.at(2, 6) = 6 * D.at(4, 4) / (L * L);
    answer.at(2, 8) = -12 * D.at(4, 4) / (L * L * L);
    answer.at(2, 12) = 6 * D.at(4, 4) / (L * L);
    answer.at(3, 3) = 12 * D.at(3, 3) / (L * L * L);
    answer.at(3, 5) = -6 * D.at(3, 3) / (L * L);
    answer.at(3, 9) = -12 * D.at(3, 3) / (L * L * L);
    answer.at(3, 11) = -6 * D.at(3, 3) / (L * L);
    answer.at(4, 4) = D.at(2, 2) / L;
    answer.at(4, 10) = -D.at(2, 2) / L;
    answer.at(5, 5) = 4 * D.at(3, 3) / L;
    answer.at(5, 9) = 6 * D.at(3, 3) / (L * L);
    answer.at(5, 11) = 2 * D.at(3, 3) / L;
    answer.at(6, 6) = 4 * D.at(4, 4) / L;
    answer.at(6, 8) = -6 * D.at(4, 4) / (L * L);
    answer.at(6, 12) = 2 * D.at(4, 4) / L;
    answer.at(7, 7) = D.at(1, 1) / L;
    answer.at(8, 8) = 12 * D.at(4, 4) / (L * L * L);
    answer.at(8, 12) = -6 * D.at(4, 4) / (L * L);
    answer.at(9, 9) = 12 * D.at(3, 3) / (L * L * L);
    answer.at(9, 11) = 6 * D.at(3, 3) / (L * L);
    answer.at(10, 10) = D.at(2, 2) / L;
    answer.at(11, 11) = 4 * D.at(3, 3) / L;
    answer.at(12, 12) = 4 * D.at(4, 4) / L;
    answer.symmetrized();
}
} // end namespace oofem