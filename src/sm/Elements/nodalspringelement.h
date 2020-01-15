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

#ifndef nodalspringelement_h
#define nodalspringelement_h

#include "sm/Elements/structuralelement.h"

///@name Input fields for spring element
//@{
#define _IFT_NodalSpringElement_Name "nodalspring"
#define _IFT_NodalSpringElement_dofmask "dofmask"
#define _IFT_NodalSpringElement_springConstants "k"
#define _IFT_NodalSpringElement_masses "m"

//@}

namespace oofem {
/**
 * This class implements a simple linear spring element connecting the given node and the ground.
 * Its purpose is to introduce additional spring stiffness and mass for nodal DOFs.
 * The orientation of spring is assumed to be constant during simulation.
 * 
 * Note: the extension of spring element to (material) nonlinear case would require
 * to introduce link to material model (like in any normal element) 
 * instead of using spring constants.
 *
 * The longitudinal spring constant [Force/Length], torsional spring constant [Force*Length/Radians].
 */
class NodalSpringElement : public StructuralElement
{
 protected:
    /// Spring constants 
    FloatArray springConstants;
    /// total mass of the spring; to be distributed to nodes
    FloatArray masses;
    /// Dof mask
    IntArray dofMask;

public:
    NodalSpringElement(int n, Domain * d);
    virtual ~NodalSpringElement() { }

    void computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep) override;
    void computeMassMatrix(FloatMatrix &answer, TimeStep *tStep) override
    { computeLumpedMassMatrix(answer, tStep); }
    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override;
    void computeInitialStressMatrix(FloatMatrix &answer, TimeStep *tStep) override
    { answer.clear(); }

    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0) override;

    int computeNumberOfDofs() override { return dofMask.giveSize(); }
    int computeNumberOfGlobalDofs() override;

    void giveDofManDofIDMask(int inode, IntArray &answer) const override;

    void updateInternalState(TimeStep *tStep) override { }
    void updateYourself(TimeStep *tStep) override { }
    int checkConsistency() override { return 1; }
    void printOutputAt(FILE *file, TimeStep *tStep) override;

#ifdef __OOFEG
    //void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep);
    //void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType);
    //void drawScalar(oofegGraphicContext &gc, TimeStep *tStep);
#endif

    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_NodalSpringElement_Name; }
    const char *giveClassName() const override { return "NodalSpringElement"; }
    void initializeFrom(InputRecord &ir) override;
    Element_Geometry_Type giveGeometryType() const override { return EGT_point; }
    bool isCast(TimeStep *tStep) override { return true; }

protected:
    void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) override
    { answer.clear(); }
    void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override
    { answer.clear(); }

    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer,
                          int lowerIndx = 1, int upperIndx = ALL_STRAINS) override
    { answer.clear(); }
    void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer) override { answer.clear(); }
    bool computeGtoLRotationMatrix(FloatMatrix &answer) override;
};
} // end namespace oofem
#endif // nodalspringelement_h
