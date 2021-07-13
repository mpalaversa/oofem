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
 
#include "sm/Elements/Plates/pltqd4dkt.h"

#include "classfactory.h"
#include "fei2dquadlin.h"
#include "gausspoint.h"
#include "sm/CrossSections/structuralcrosssection.h"
#include "sm/Materials/structuralms.h"

namespace oofem {
REGISTER_Element(PltQd4DKT);

FEI2dQuadLin PltQd4DKT::interp_lin(1, 2);

PltQd4DKT :: PltQd4DKT(int n, Domain* aDomain) : QdPlate(n, aDomain)
{
    numberOfDofMans = 4;
    numberOfGaussPoints = 4;
}

void
PltQd4DKT::computeBmatrixAt(double ksi, double eta, FloatMatrix& answer)
// Returns the [3x12] strain-displacement matrix {B} of the receiver,
// evaluated at gp.
{
    // get node coordinates
    std::vector< FloatArray > locCoords = giveNodeCoordinates();
    double x1 = locCoords.at(0)[0], x2 = locCoords.at(1)[0], x3 = locCoords.at(2)[0], x4 = locCoords.at(3)[0], y1 = locCoords.at(0)[1], y2 = locCoords.at(1)[1], y3 = locCoords.at(2)[1], y4 = locCoords.at(3)[1], z1 = locCoords.at(0)[2], z2 = locCoords.at(1)[2], z3 = locCoords.at(2)[2], z4 = locCoords.at(3)[2];

    // geometrical characteristics of element sides
    double dx4 = x2 - x1;
    double dy4 = y2 - y1;
    double l4 = sqrt(dx4 * dx4 + dy4 * dy4);

    double dx5 = x3 - x2;
    double dy5 = y3 - y2;
    double l5 = sqrt(dx5 * dx5 + dy5 * dy5);

    double dx6 = x4 - x3;
    double dy6 = y4 - y3;
    double l6 = sqrt(dx6 * dx6 + dy6 * dy6);

    double dx7 = x1 - x4;
    double dy7 = y1 - y4;
    double l7 = sqrt(dx7 * dx7 + dy7 * dy7);

    double c4 = dy4 / l4;
    double s4 = -dx4 / l4;

    double c5 = dy5 / l5;
    double s5 = -dx5 / l5;

    double c6 = dy6 / l6;
    double s6 = -dx6 / l6;

    double c7 = dy7 / l7;
    double s7 = -dx7 / l7;

    // transformation matrix from vertex dofs (w,fi_x, fi_y) to quadratic rotation DOFs
    double T101 = -3. / 2. / l4 * c4;
    double T102 = -1. / 4. * c4 * c4 + 1. / 2. * s4 * s4;
    double T103 = -1. / 4. * c4 * s4 - 1. / 2. * c4 * s4;
    double T104 = 3. / 2. / l4 * c4;
    double T105 = -1. / 4. * c4 * c4 + 1. / 2. * s4 * s4;
    double T106 = -1. / 4. * c4 * s4 - 1. / 2. * c4 * s4;

    double T201 = -3. / 2. / l4 * s4;
    double T202 = -1. / 4. * c4 * s4 - 1. / 2. * c4 * s4;
    double T203 = -1. / 4. * s4 * s4 + 1. / 2. * c4 * c4;
    double T204 = 3. / 2. / l4 * s4;
    double T205 = -1. / 4. * c4 * s4 - 1. / 2. * c4 * s4;
    double T206 = -1. / 4. * s4 * s4 + 1. / 2. * c4 * c4;

    double T304 = -3. / 2. / l5 * c5;
    double T305 = -1. / 4. * c5 * c5 + 1. / 2. * s5 * s5;
    double T306 = -1. / 4. * c5 * s5 - 1. / 2. * c5 * s5;
    double T307 = 3. / 2. / l5 * c5;
    double T308 = -1. / 4. * c5 * c5 + 1. / 2. * s5 * s5;
    double T309 = -1. / 4. * c5 * s5 - 1. / 2. * c5 * s5;

    double T404 = -3. / 2. / l5 * s5;
    double T405 = -1. / 4. * c5 * s5 - 1. / 2. * c5 * s5;
    double T406 = -1. / 4. * s5 * s5 + 1. / 2. * c5 * c5;
    double T407 = 3. / 2. / l5 * s5;
    double T408 = -1. / 4. * c5 * s5 - 1. / 2. * c5 * s5;
    double T409 = -1. / 4. * s5 * s5 + 1. / 2. * c5 * c5;

    double T507 = -3. / 2. / l6 * c6;
    double T508 = -1. / 4. * c6 * c6 + 1. / 2. * s6 * s6;
    double T509 = -1. / 4. * c6 * s6 - 1. / 2. * c6 * s6;
    double T510 = 3. / 2. / l6 * c6;
    double T511 = -1. / 4. * c6 * c6 + 1. / 2. * s6 * s6;
    double T512 = -1. / 4. * c6 * s6 - 1. / 2. * c6 * s6;

    double T607 = -3. / 2. / l6 * s6;
    double T608 = -1. / 4. * c6 * s6 - 1. / 2. * c6 * s6;
    double T609 = -1. / 4. * s6 * s6 + 1. / 2. * c6 * c6;
    double T610 = 3. / 2. / l6 * s6;
    double T611 = -1. / 4. * c6 * s6 - 1. / 2. * c6 * s6;
    double T612 = -1. / 4. * s6 * s6 + 1. / 2. * c6 * c6;

    double T701 = 3. / 2. / l7 * c7;
    double T702 = -1. / 4. * c7 * c7 + 1. / 2. * s7 * s7;
    double T703 = -1. / 4. * c7 * s7 - 1. / 2. * c7 * s7;
    double T710 = -3. / 2. / l7 * c7;
    double T711 = -1. / 4. * c7 * c7 + 1. / 2. * s7 * s7;
    double T712 = -1. / 4. * c7 * s7 - 1. / 2. * c7 * s7;

    double T801 = 3. / 2. / l7 * s7;
    double T802 = -1. / 4. * c7 * s7 - 1. / 2. * c7 * s7;
    double T803 = -1. / 4. * s7 * s7 + 1. / 2. * c7 * c7;
    double T810 = -3. / 2. / l7 * s7;
    double T811 = -1. / 4. * c7 * s7 - 1. / 2. * c7 * s7;
    double T812 = -1. / 4. * s7 * s7 + 1. / 2. * c7 * c7;

    // derivatives of quadratic interpolation functions
    // we do not have "midside" nodes -> explicit here
    double N1dk = 0.25 * (2.0 * ksi + eta) * (1.0 + eta);
    double N2dk = 0.25 * (2.0 * ksi - eta) * (1.0 + eta);
    double N3dk = 0.25 * (2.0 * ksi + eta) * (1.0 - eta);
    double N4dk = 0.25 * (2.0 * ksi - eta) * (1.0 - eta);
    double N7dk = -ksi * (1.0 - eta);
    double N8dk = 0.5 * (1.0 - eta * eta);
    double N5dk = -ksi * (1.0 + eta);
    double N6dk = -0.5 * (1.0 - eta * eta);

    double N3de = 0.25 * (2.0 * eta + ksi) * (1.0 - ksi);
    double N4de = 0.25 * (2.0 * eta - ksi) * (1.0 + ksi);
    double N1de = 0.25 * (2.0 * eta + ksi) * (1.0 + ksi);
    double N2de = 0.25 * (2.0 * eta - ksi) * (1.0 - ksi);
    double N7de = -0.5 * (1.0 - ksi * ksi);
    double N8de = -eta * (1.0 + ksi);
    double N5de = 0.5 * (1.0 - ksi * ksi);
    double N6de = -eta * (1.0 - ksi);

    double detJ = 1. / 8. * ((y4 - y2) * (x3 - x1) - (y3 - y1) * (x4 - x2)) +
        ksi / 8 * ((y3 - y4) * (x2 - x1) - (y2 - y1) * (x3 - x4)) +
        eta / 8 * ((y4 - y1) * (x3 - x2) - (y3 - y2) * (x4 - x1));

    double dxdk = -1.0 / detJ * ((y3 - y2) + (y4 - y1 + ksi * (y1 - y2 + y3 - y4))) / 4.0;
    double dxde = 1.0 / detJ * ((y2 - y1) + (y3 - y4 + eta * (y1 - y2 + y3 - y4))) / 4.0;
    double dydk = 1.0 / detJ * ((x3 - x2) + (x4 - x1 + ksi * (x1 - x2 + x3 - x4))) / 4.0;
    double dyde = -1.0 / detJ * ((x2 - x1) + (x3 - x4 + eta * (x1 - x2 + x3 - x4))) / 4.0;

    double dN102 = N1dk * dxdk + N1de * dxde;
    double dN104 = N2dk * dxdk + N2de * dxde;
    double dN106 = N3dk * dxdk + N3de * dxde;
    double dN108 = N4dk * dxdk + N4de * dxde;
    double dN110 = N5dk * dxdk + N5de * dxde;
    double dN112 = N6dk * dxdk + N6de * dxde;
    double dN114 = N7dk * dxdk + N7de * dxde;
    double dN116 = N8dk * dxdk + N8de * dxde;

    double dN201 = -N1dk * dydk - N1de * dyde;
    double dN203 = -N2dk * dydk - N2de * dyde;
    double dN205 = -N3dk * dydk - N3de * dyde;
    double dN207 = -N4dk * dydk - N4de * dyde;
    double dN209 = -N5dk * dydk - N5de * dyde;
    double dN211 = -N6dk * dydk - N6de * dyde;
    double dN213 = -N7dk * dydk - N7de * dyde;
    double dN215 = -N8dk * dydk - N8de * dyde;

    answer.resize(3, 12);

    answer.at(1, 1) = T201 * dN110 + T801 * dN116;
    answer.at(1, 2) = T202 * dN110 + T802 * dN116;
    answer.at(1, 3) = dN102 + T203 * dN110 + T803 * dN116;
    answer.at(1, 4) = T204 * dN110 + T404 * dN112;
    answer.at(1, 5) = T205 * dN110 + T405 * dN112;
    answer.at(1, 6) = dN104 + T206 * dN110 + T406 * dN112;
    answer.at(1, 7) = T407 * dN112 + T607 * dN114;
    answer.at(1, 8) = T408 * dN112 + T608 * dN114;
    answer.at(1, 9) = dN106 + T409 * dN112 + T609 * dN114;
    answer.at(1, 10) = T610 * dN114 + T810 * dN116;
    answer.at(1, 11) = T611 * dN114 + T811 * dN116;
    answer.at(1, 12) = dN108 + T612 * dN114 + T812 * dN116;

    answer.at(2, 1) = T101 * dN209 + T701 * dN215;
    answer.at(2, 2) = dN201 + T102 * dN209 + T702 * dN215;
    answer.at(2, 3) = T103 * dN209 + T703 * dN215;
    answer.at(2, 4) = T104 * dN209 + T304 * dN211;
    answer.at(2, 5) = dN203 + T105 * dN209 + T305 * dN211;
    answer.at(2, 6) = T106 * dN209 + T306 * dN211;
    answer.at(2, 7) = T307 * dN211 + T507 * dN213;
    answer.at(2, 8) = dN205 + T308 * dN211 + T508 * dN213;
    answer.at(2, 9) = T309 * dN211 + T509 * dN213;
    answer.at(2, 10) = T510 * dN213 + T710 * dN215;
    answer.at(2, 11) = dN207 + T511 * dN213 + T711 * dN215;
    answer.at(2, 12) = T512 * dN213 + T712 * dN215;

    answer.at(3, 1) = -T101 * dN110 - T201 * dN209 - T701 * dN116 - T801 * dN215;
    answer.at(3, 2) = -dN102 - T102 * dN110 - T202 * dN209 - T702 * dN116 - T802 * dN215;
    answer.at(3, 3) = -dN201 - T103 * dN110 - T203 * dN209 - T703 * dN116 - T803 * dN215;
    answer.at(3, 4) = -T104 * dN110 - T204 * dN209 - T304 * dN112 - T404 * dN211;
    answer.at(3, 5) = -dN104 - T105 * dN110 - T205 * dN209 - T305 * dN112 - T405 * dN211;
    answer.at(3, 6) = -dN203 - T106 * dN110 - T206 * dN209 - T306 * dN112 - T406 * dN211;
    answer.at(3, 7) = -T307 * dN112 - T407 * dN211 - T507 * dN114 - T607 * dN213;
    answer.at(3, 8) = -dN106 - T308 * dN112 - T408 * dN211 - T508 * dN114 - T608 * dN213;
    answer.at(3, 9) = -dN205 - T309 * dN112 - T409 * dN211 - T509 * dN114 - T609 * dN213;
    answer.at(3, 10) = -T510 * dN114 - T610 * dN213 - T710 * dN116 - T810 * dN215;
    answer.at(3, 11) = -dN108 - T511 * dN114 - T611 * dN213 - T711 * dN116 - T811 * dN215;
    answer.at(3, 12) = -dN207 - T512 * dN114 - T612 * dN213 - T712 * dN116 - T812 * dN215;
}

void
PltQd4DKT::computeConstitutiveMatrixAt(FloatMatrix& answer, MatResponseMode rMode, GaussPoint* gp, TimeStep* tStep)
{
    answer = this->giveStructuralCrossSection()->giveKirchhoffPlateStiffMtrx(rMode, gp, tStep);
}

void
PltQd4DKT::computeCurvaturesAt(FloatArray& answer, double xi, double eta, TimeStep* tStep) {
    FloatArray u;
    this->computeVectorOf(VM_Total, tStep, u);
    /* This is to be uncommented once adapted for this element (if necessary).
    if (initialDisplacements) {
        u.subtract(*initialDisplacements);
    }*/

    FloatMatrix BMatrix;
    computeBmatrixAt(xi, eta, BMatrix);
    answer.beProductOf(BMatrix, u);
}

bool
PltQd4DKT::computeGtoLRotationMatrix(FloatMatrix& answer)
// Returns the rotation matrix of the receiver of the size [24,24]
// r(local) = T * r(global)
// for one node (r written transposed): {w,R_u,R_v} = T * {u,v,w,R_u,R_v,R_w}
{
    if (GtoLRotationMatrix == NULL)
        QdElement::computeGtoLRotationMatrix();

    answer.resize(12, 24);

    for (int i = 1; i <= 3; i++) {
        answer.at(2, i + 3) = answer.at(5, i + 9) = answer.at(8, i + 15) = answer.at(11, i + 21) = GtoLRotationMatrix->at(1, i);
        answer.at(3, i + 3) = answer.at(6, i + 9) = answer.at(9, i + 15) = answer.at(12, i + 21) = GtoLRotationMatrix->at(2, i);
        answer.at(1, i) = answer.at(4, i + 6) = answer.at(7, i + 12) = answer.at(10, i + 18) = GtoLRotationMatrix->at(3, i);
    }

    return 1;
}

void
PltQd4DKT::computeStrainVectorAt(FloatArray& answer, double xi, double eta, TimeStep* tStep) {
    FloatArray curvatures;
    switch (outputType) {
    case OutputType::Standard:
        computeCurvaturesAt(curvatures, xi, eta, tStep);
        answer.beScaled(outputAtZ, curvatures);
        break;
    case OutputType::Principal:
        OOFEM_ERROR("Not yet implemented for a %s element.", giveClassName());
        break;
    case OutputType::VM:
        OOFEM_ERROR("Not yet implemented for a %s element.", giveClassName());
        break;
    case OutputType::All:
        OOFEM_ERROR("Not yet implemented for a %s element.", giveClassName());
        break;
    default:
        OOFEM_ERROR("Something went wrong. An unknown output category requested for element %d.", giveGlobalNumber());
        break;
    }
}

void
PltQd4DKT::computeStressVector(FloatArray& answer, const FloatArray& strain, GaussPoint* gp, TimeStep* tStep) {
    answer = this->giveStructuralCrossSection()->giveRealStress_KirchhoffPlate(strain, gp, tStep);
}

double
PltQd4DKT::computeSurfaceVolumeAround(GaussPoint* gp, int iSurf)
{
    double detJ, weight;

    weight = gp->giveWeight();
    detJ = fabs(this->interp_lin.giveTransformationJacobian(gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this)));
    return detJ * weight;
}

double
PltQd4DKT::computeVolumeAround(GaussPoint* gp)
// Returns the portion of the receiver which is attached to gp. For this element, it is equal to the surface volume.
{
    return computeSurfaceVolumeAround(gp, 1);
}

FEInterpolation*
PltQd4DKT::giveInterpolation() const {
    return &interp_lin;
}

bool
PltQd4DKT::giveRotationMatrix(FloatMatrix& answer)
{
    bool is_GtoL, is_NtoG;
    FloatMatrix GtoL, NtoG;
    IntArray nodes;
    nodes.enumerate(this->giveNumberOfDofManagers());

    is_GtoL = this->computeGtoLRotationMatrix(GtoL);
    is_NtoG = this->computeDofTransformationMatrix(NtoG, nodes, true);

    if (is_GtoL && NtoG.isNotEmpty()) {
        answer.beProductOf(GtoL, NtoG);
    }
    else if (is_GtoL) {
        answer = GtoL;
    }
    else if (is_NtoG) {
        answer = NtoG;
    }
    else {
        answer.clear();
        return false;
    }
    return true;
}

void
PltQd4DKT::initializeFrom(InputRecord& ir)
{
    StructuralElement::initializeFrom(ir);

    int outputAtXYTemp, outputTypeTemp, outputAtZTemp;
    outputAtXYTemp = outputTypeTemp = outputAtZTemp = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, outputAtXYTemp, _IFT_PltQd4DKT_outputAtXY);
    IR_GIVE_OPTIONAL_FIELD(ir, outputTypeTemp, _IFT_PltQd4DKT_outputType);
    IR_GIVE_OPTIONAL_FIELD(ir, outputAtZ, _IFT_PltQd4DKT_outputAtZ);

    switch (outputAtXYTemp) {
    case 1:
        outputAtXY = OutputLocationXY::GaussPoints;
        break;
    case 2:
        outputAtXY = OutputLocationXY::Centre;
        break;
    case 3:
        outputAtXY = OutputLocationXY::Corners;
        break;
    case 4:
        outputAtXY = OutputLocationXY::All;
        break;
    default:
        outputAtXY = OutputLocationXY::GaussPoints;
        break;
    }
    switch (outputTypeTemp) {
    case 1:
        outputType = OutputType::Standard;
        break;
    case 2:
        outputType = OutputType::Principal;
        break;
    case 3:
        outputType = OutputType::VM;
        break;
    case 4:
        outputType = OutputType::All;
        break;
    default:
        outputType = OutputType::Standard;
        break;
    }
}

void
PltQd4DKT::printOutputAt(FILE* file, TimeStep* tStep)
{
    FloatArray v;
    fprintf(file, "element %d (%8d):\n", this->giveLabel(), number);
    int nPointsXY{ 0 };
    switch (outputAtXY)
    {
    case OutputLocationXY::GaussPoints:
        nPointsXY = giveIntegrationRulesArray()[0]->giveNumberOfIntegrationPoints();
        break;
    case OutputLocationXY::Centre:
        nPointsXY = 1;
        break;
    case OutputLocationXY::Corners:
        nPointsXY = giveNumberOfDofManagers();
        break;
    case OutputLocationXY::All:
        nPointsXY = giveIntegrationRulesArray()[0]->giveNumberOfIntegrationPoints() + giveNumberOfDofManagers() + 1;
        break;
    default:
        OOFEM_ERROR("Something went wrong. The following options for output location at XY plane can be selected for a membrane element: 'All', 'Centroid', 'Gauss points' and 'Corners'.");
        break;
    }

    StructuralMaterialStatus* ms;
    GaussPoint* gp;
    for (int i = 0; i < nPointsXY; i++) {
        switch (outputAtXY)
        {
        case OutputLocationXY::GaussPoints:
            fprintf(file, "  GP %d :", i + 1);
            gp = integrationRulesArray[0]->getIntegrationPoint(i);
            ms = static_cast<StructuralMaterialStatus*>(gp->giveMaterialStatus());

            /*
            fprintf(file, "  forces     ");
            for (auto& val : this->giveMidplaneIPValue(i, IST_ShellForceTensor, tStep)) {
                fprintf(file, " %.4e", val);
            }

            fprintf(file, "\n          moments    ");
            for (auto& val : this->giveMidplaneIPValue(i, IST_ShellMomentTensor, tStep)) {
                fprintf(file, " %.4e", val);
            }

            fprintf(file, "\n          strains    ");
            for (auto& val : ms->giveStrainVector()) {
                fprintf(file, " %.4e", val);
            }

            fprintf(file, "\n          curvatures ");
            for (auto& val : this->giveMidplaneIPValue(i, IST_CurvatureTensor, tStep)) {
                fprintf(file, " %.4e", val);
            }*/

            fprintf(file, "\n          at z = %f :\n", outputAtZ);


            fprintf(file, "                      strains    ");
            for (auto& val : ms->giveStrainVector()) {
                fprintf(file, " %.4e", val);
            }

            fprintf(file, "\n                      stresses   ");
            for (auto& val : ms->giveStressVector()) {
                fprintf(file, " %.4e", val);
            }
            fprintf(file, "\n");
            break;
        case OutputLocationXY::Centre:
            fprintf(file, "  Centre");
            gp = integrationRulesArray[0]->getIntegrationPoint(0);
            ms = static_cast<StructuralMaterialStatus*>(gp->giveMaterialStatus());

            /*
            fprintf(file, "  forces     ");
            for (auto& val : this->giveMidplaneIPValue(i, IST_ShellForceTensor, tStep)) {
                fprintf(file, " %.4e", val);
            }

            fprintf(file, "\n          moments    ");
            for (auto& val : this->giveMidplaneIPValue(i, IST_ShellMomentTensor, tStep)) {
                fprintf(file, " %.4e", val);
            }

            fprintf(file, "\n          strains    ");
            for (auto& val : ms->giveStrainVector()) {
                fprintf(file, " %.4e", val);
            }

            fprintf(file, "\n          curvatures ");
            for (auto& val : this->giveMidplaneIPValue(i, IST_CurvatureTensor, tStep)) {
                fprintf(file, " %.4e", val);
            }*/

            fprintf(file, "\n          at z = %f :\n", outputAtZ);


            fprintf(file, "                      strains    ");
            for (auto& val : ms->giveStrainVector()) {
                fprintf(file, " %.4e", val);
            }

            fprintf(file, "\n                      stresses   ");
            for (auto& val : ms->giveStressVector()) {
                fprintf(file, " %.4e", val);
            }
            break;
        case OutputLocationXY::Corners:
            fprintf(file, "  Node %d :", i + 1);
            break;
        case OutputLocationXY::All:
            if (i < 4)
                fprintf(file, "  GP %d :", i + 1);
            else if (i >= 4 && i < 8)
                fprintf(file, "  Node %d :", i - 3);
            else
                fprintf(file, "  Centroid :");
            break;
        default:
            OOFEM_ERROR("Something went wrong. The following options for output location at XY plane can be selected for a membrane element: 'All', 'Centroid', 'Gauss points' and 'Corners'.");
            break;
        }
    }
}
}