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

#ifndef hybbeam_h
#define hybbeam_h

#include "sm/Elements/Beams/eccbeam.h"

///@name Input fields for HybBeam
//@{
#define _IFT_HybBeam_Name "hybbeam"
#define _IFT_HybBeam_refnode "refnode"
//@}

namespace oofem {
    /**
     * This class implements a 3D beam element based on the Euler-Bernoulli beam theory.
     * The element has 2 nodes and 6 DOFs at each node (three translations and three rotations).
     * The axial translations (compression/extension) and the rotation about the longitudinal
     * axis (torsion) is approximated by linear basis functions. The bending translations are
     * approximated by a 3rd order polynomial. The bending rotations are calulcated as their
     * first derivatives.
     * The element is not isoparametric and its stiffness matrix is given directly (no numerical
     * integration is used). The element supports only static linear analysis at this moment.
     * The difference between this element and EccBeam is that this elements tries to take into
     * account the plating attached to it and is thus suitable for structures such as stiffeners attached
     * to plating where the point of attachement is away from the stiffener's centroid. The difference
     * between this element and ModEccBeam is that the user can input a "real" moment of inertia of the
     * cross section, i.e. one based on the effective breadth.
     */
class HybBeam : public EccBeam
{
private:
    int refNode = 1;

public:
    HybBeam(int n, Domain* d);
    virtual ~HybBeam();

    const char* giveClassName() const override { return "HybBeam"; }
    const char* giveInputRecordName() const override { return _IFT_HybBeam_Name; }

    void computeStiffnessMatrix(FloatMatrix& answer, MatResponseMode rMode, TimeStep* tStep) override;
};
} // end namespace oofem
#endif // hybbeam_h