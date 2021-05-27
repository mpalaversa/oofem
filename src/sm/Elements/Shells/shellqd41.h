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
 
#ifndef shellqd41_h
#define shellqd41_h

#include "sm/Elements/nlstructuralelement.h"
#include "sm/Elements/PlaneStress/linquad3d_planestress.h"
#include "sm/Elements/Plates/qdkt.h"

#define _IFT_ShellQd41_Name "shellqd41"

namespace oofem {

class ShellQd41 : public NLStructuralElement
{
    LinQuad3DPlaneStress* membrane;
    QDKTPlate* plate;

public:
    ShellQd41(int n, Domain* d);
    virtual ~ShellQd41() { }

    void computeBmatrixAt(GaussPoint* gp, FloatMatrix& answer, int lowerIndx = 1, int upperIndx = ALL_STRAINS) override;
    void computeBmatrixAt(int elementVertex, FloatMatrix& answer, int lowerIndx = 1, int upperIndx = ALL_STRAINS) { }
    void computeBmatrixPlateAt(GaussPoint* gp, FloatMatrix& answer);
    void computeConstitutiveMatrixAt(FloatMatrix& answer, MatResponseMode rMode, GaussPoint* gp, TimeStep* tStep) override;
    void computeGaussPoints() override;
    bool computeGtoLRotationMatrix(FloatMatrix& answer) override;
    int computeLoadGToLRotationMtrx(FloatMatrix& answer) override { return membrane->computeLoadGToLRotationMtrx(answer); }
    int computeNumberOfDofs() override { return 24; }
    void computeStiffnessMatrix(FloatMatrix& answer, MatResponseMode rMode, TimeStep* tStep) override;
    void computeStressVector(FloatArray& answer, const FloatArray& strain, GaussPoint* gp, TimeStep* tStep) { answer = this->giveStructuralCrossSection()->giveRealStress_PlaneStress(strain, gp, tStep); }
    void computeStressVectorTop(FloatArray& answer, const FloatArray& strain, GaussPoint* gp, TimeStep* tStep) { answer = this->giveStructuralCrossSection()->giveRealStress_PlaneStress(strain, gp, tStep); }
    void computeStressVectorBottom(FloatArray& answer, const FloatArray& strain, GaussPoint* gp, TimeStep* tStep) { answer = this->giveStructuralCrossSection()->giveRealStress_PlaneStress(strain, gp, tStep); }
    double computeVolumeAround(GaussPoint* gp) override;
    const char* giveClassName() const override { return "ShellQd41"; }
    int giveDefaultIntegrationRule() const override { return plate->giveDefaultIntegrationRule(); }
    void giveDofManDofIDMask(int inode, IntArray&) const override;
    const char* giveInputRecordName() const override { return _IFT_ShellQd41_Name; }
    integrationDomain giveIntegrationDomain() const override { return _Square; }
    IntegrationRule* giveIntegrationRule(int i) override { return plate->giveIntegrationRule(i); }
    FEInterpolation* giveInterpolation() const override { return plate->giveInterpolation(); }
    MaterialMode giveMaterialMode() override { return _PlaneStress; }
    std::vector< FloatArray > giveNodeCoordinates();
    void initializeFrom(InputRecord& ir) override;
    void setCrossSection(int csIndx) override;

    //IntegrationRule* giveDefaultIntegrationRulePtr() override { return plate->giveDefaultIntegrationRulePtr(); }
    //Element_Geometry_Type giveGeometryType() const override { return EGT_quad_1; }

private:
    float drillCoeff = 100;
};
}
#endif