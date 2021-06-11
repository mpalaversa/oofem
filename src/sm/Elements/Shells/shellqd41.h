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
#define _IFT_ShellQd41_outputAtXY "outputatxy"
#define _IFT_ShellQd41_outputType "outputtype"
#define _IFT_ShellQd41_outputAtZ "outputatz"
#define _IFT_ShellQd41_outputCategory "outputcategory"

namespace oofem {
    // Enumeration for output location of strains and stresses in element's x-y plane.
    enum class OutputLocationXY {
        GaussPoints,
        Centroid,
        Corners,
        All
    };
    // Enumeration for output category of strains and stresses.
    enum class OutputCategory {
        Membrane,
        Plate,
        Combined,
        All
    };
    // Enumeration for output type of strains and stresses.
    enum class OutputType {
        Standard,
        Principal,
        VM,
        All
    };

class ShellQd41 : public NLStructuralElement
{
    LinQuad3DPlaneStress* membrane;
    QDKTPlate* plate;

public:
    ShellQd41(int n, Domain* d);
    virtual ~ShellQd41() { }

    void computeBmatrixAt(GaussPoint* gp, FloatMatrix& answer, int lowerIndx = 1, int upperIndx = ALL_STRAINS) override;
    // Computes B matrix for the plate part of the element at natural coordinates (xi, eta).
    void computeBmatrixPlateAt(double xi, double eta, FloatMatrix& answer);
    void computeBmatrixPlateAt(GaussPoint* gp, FloatMatrix& answer);
    void computeConstitutiveMatrixAt(FloatMatrix& answer, MatResponseMode rMode, GaussPoint* gp, TimeStep* tStep) override;
    void computeGaussPoints() override;
    bool computeGtoLRotationMatrix(FloatMatrix& answer) override;
    int computeLoadGToLRotationMtrx(FloatMatrix& answer) override { return membrane->computeLoadGToLRotationMtrx(answer); }
    int computeLoadLSToLRotationMatrix(FloatMatrix& answer, int iSurf, GaussPoint* gp) override { return plate->computeLoadLSToLRotationMatrix(answer, iSurf, gp); }
    void computeMembraneStrainVectorAt(FloatArray& answer, double xi, double eta, TimeStep* tStep);
    void computePlateCurvaturesAt(FloatArray& answer, double xi, double eta, TimeStep* tStep);
    void computePlateStrainVectorAt(FloatArray& answer, double xi, double eta, TimeStep* tStep);
    int computeNumberOfDofs() override { return 24; }
    void computeStiffnessMatrix(FloatMatrix& answer, MatResponseMode rMode, TimeStep* tStep) override;
    void computeStrainVector(FloatArray& answer, GaussPoint* gp, TimeStep* tStep) override;
    void computeStrainVectorAtCentroid(FloatArray& answer, TimeStep* tStep);
    void computeStressVector(FloatArray& answer, const FloatArray& strain, GaussPoint* gp, TimeStep* tStep) override;
    void computeStressVectorAtCentroid(FloatArray& answer, TimeStep* tStep, const FloatArray& strain = 0);
    void computeSurfaceNMatrix(FloatMatrix& answer, int boundaryID, const FloatArray& lcoords) override;
    double computeSurfaceVolumeAround(GaussPoint* gp, int iSurf) override { return membrane->computeVolumeAround(gp); }
    double computeVolumeAround(GaussPoint* gp) override;
    OutputCategory getOutputCategory() { return outputCategory; }
    OutputLocationXY getOutputLocationInXYPlane() { return outputAtXY; }
    double getOutputLocationInZ() { return outputAtZ; }
    OutputType getOutputType() { return outputType; }
    const char* giveClassName() const override { return "ShellQd41"; }
    int giveDefaultIntegrationRule() const override { return plate->giveDefaultIntegrationRule(); }
    void giveDofManDofIDMask(int inode, IntArray&) const override;
    const char* giveInputRecordName() const override { return _IFT_ShellQd41_Name; }
    integrationDomain giveIntegrationDomain() const override { return _Square; }
    IntegrationRule* giveIntegrationRule(int i) override { return plate->giveIntegrationRule(i); }
    // giveInternalForcesVector is used only in non-linear analysis. This should be changed when non-linear analysis capabilities are implemented.
    void giveInternalForcesVector(FloatArray& answer, TimeStep* tStep, int useUpdatedGpRecord) override { answer.zero(); }
    FEInterpolation* giveInterpolation() const override { return plate->giveInterpolation(); }
    MaterialMode giveMaterialMode() override { return _PlaneStress; }
    std::vector< FloatArray > giveNodeCoordinates();
    void giveCharacteristicOutput( FloatArray &answer, TimeStep *tStep ) override { this->getStressesTopBottom( answer, tStep ); }
    void getStressesTopBottom(FloatArray& answer, TimeStep* tStep);
    void giveSurfaceDofMapping(IntArray& answer, int iSurf) const override;
    void initializeFrom(InputRecord& ir) override;
    void setCrossSection(int csIndx) override;
    void updateInternalState(TimeStep* tStep) override;
    void updateLocalNumbering(EntityRenumberingFunctor& f) override;

    //IntegrationRule* giveDefaultIntegrationRulePtr() override { return plate->giveDefaultIntegrationRulePtr(); }
    //Element_Geometry_Type giveGeometryType() const override { return EGT_quad_1; }

private:
    OutputLocationXY outputAtXY = OutputLocationXY::GaussPoints;
    OutputCategory outputCategory = OutputCategory::Combined;
    OutputType outputType = OutputType::Standard;
    double outputAtZ = 0.0;
   
    float drillCoeff = 100;
};
}
#endif