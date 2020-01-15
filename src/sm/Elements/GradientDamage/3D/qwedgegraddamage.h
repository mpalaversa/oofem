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

#ifndef qwedgegraddamage_h
#define qwedgegraddamage_h

#include "../sm/Elements/3D/qwedge.h"
#include "../sm/Elements/GradientDamage/graddamageelement.h"

#define _IFT_QWedgeGradDamage_Name "qwedgegraddamage"

namespace oofem {
class FEI3dWedgeLin;

/**
 * @author M. Horak
 */
class QWedgeGradDamage : public QWedge, public GradientDamageElement
{
protected:
    static FEI3dWedgeLin interpolation_lin;

public:
    QWedgeGradDamage(int, Domain *);
    virtual ~QWedgeGradDamage() { }

    void initializeFrom(InputRecord &ir) override;
    void giveDofManDofIDMask(int inode, IntArray &answer) const override;
    void giveDofManDofIDMask_u(IntArray &answer) const override;
    void giveDofManDofIDMask_d(IntArray &answer) const override;
    
    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_QWedgeGradDamage_Name; }
    const char *giveClassName() const override { return "QWedgeGradDamage"; }
    int computeNumberOfDofs() override { return 51; }
    MaterialMode giveMaterialMode() override { return _3dMat; }

protected:
    void computeGaussPoints() override;
    void computeNdMatrixAt(GaussPoint *gp, FloatArray &answer) override;
    void computeBdMatrixAt(GaussPoint *gp, FloatMatrix &answer) override;
    StructuralElement *giveStructuralElement() override { return this; }
    NLStructuralElement *giveNLStructuralElement() override { return this; }

    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0) override { GradientDamageElement :: giveInternalForcesVector(answer, tStep, useUpdatedGpRecord); }
    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override { GradientDamageElement :: computeStiffnessMatrix(answer, rMode, tStep); }
    void giveLocationArray_u(IntArray &answer) override { }
    void giveLocationArray_d(IntArray &answer) override { }

};
}
#endif // end namespace oofem
