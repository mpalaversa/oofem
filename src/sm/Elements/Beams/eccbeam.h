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

#ifndef eccbeam_h
#define eccbeam_h

#include "sm/Elements/Beams/beambaseelement.h"

///@name Input fields for EccBeam
//@{
#define _IFT_EccBeam_Name "eccbeam"
#define _IFT_EccBeam_refnode "refnode"
#define _IFT_EccBeam_offset_y "offset_y"
#define _IFT_EccBeam_offset_z "offset_z"
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
 */
class EccBeam : public BeamBaseElement
{
private:
    double offset_y = 0.0;
    double offset_z = 0.0;
    int refNode = 1;

protected:
    void computeBmatrixAt(double x, FloatMatrix& answer);
    void computeBmatrixAt(GaussPoint*, FloatMatrix&, int = 1, int = ALL_STRAINS) override;
    void computeConstitutiveMatrixAt(FloatMatrix& answer, MatResponseMode rMode, GaussPoint* gp, TimeStep* tStep) override;

public:
    EccBeam(int n, Domain* d);
    virtual ~EccBeam();

    void computeGaussPoints() override;
    bool computeGtoLRotationMatrix(FloatMatrix& answer) override;
    double computeLength() override;
    int computeNumberOfDofs() override { return 12; }
    virtual void computeStiffnessMatrix(FloatMatrix& answer, MatResponseMode rMode, TimeStep* tStep) override;
    void computeStressVector(FloatArray& answer, const FloatArray& strain, GaussPoint* gp, TimeStep* tStep) override;
    double computeVolumeAround(GaussPoint* gp) override;
    const char* giveClassName() const override { return "EccBeam"; }
    void giveDofManDofIDMask(int inode, IntArray&) const override;
    const char* giveInputRecordName() const override { return _IFT_EccBeam_Name; }
    integrationDomain giveIntegrationDomain() const override { return _Line; }
    void giveInternalForcesVector(FloatArray& answer, TimeStep* tStep, int useUpdatedGpRecord = 0) override;
    int giveLocalCoordinateSystem(FloatMatrix& answer) override;
    MaterialMode giveMaterialMode() override { return _3dBeam; }
    void initializeFrom(InputRecord& ir) override;
    void printOutputAt(FILE* file, TimeStep* tStep) override;
    void updateLocalNumbering(EntityRenumberingFunctor& f) override;
};
} // end namespace oofem
#endif // eccbeam_h