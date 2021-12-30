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

#include "sm/Elements/Beams/eccbeam.h"
#include "classfactory.h"
#include "engngm.h"
#include "fei3dlinelin.h"
#include "gaussintegrationrule.h"
#include "node.h"

namespace oofem {
REGISTER_Element(EccBeam);

EccBeam::EccBeam(int n, Domain* aDomain) : BeamBaseElement(n, aDomain)
{
    numberOfDofMans = 2;
    numberOfGaussPoints = 1;
}

EccBeam :: ~EccBeam()
{

}

void
EccBeam::computeBmatrixAt(double x, FloatMatrix& answer)
// eps = {\eps_x, \kappa_x, \kappa_y, \kappa_z}^T
{
    double l{ 0.0 }, rod{ 0.0 }, beam1{ 0.0 }, beam2{ 0.0 }, beam3{ 0.0 };

    l = this->computeLength();
    rod = 1 / l;
    beam1 = (6 / (l * l)) * ((2 * x) / l - 1);
    beam2 = (2 / l) * (2 - 3 * x / l);
    beam3 = (2 / l) * (1 - 3 * x / l);

    answer.resize(4, 12);

    answer.at(1, 1) = -rod;
    answer.at(1, 7) = rod;

    answer.at(2, 4) = -rod;
    answer.at(2, 10) = rod;

    answer.at(3, 3) = beam1;
    answer.at(3, 5) = beam2;
    answer.at(3, 9) = -beam1;
    answer.at(3, 11) = beam3;

    answer.at(4, 2) = beam1;
    answer.at(4, 6) = -beam2;
    answer.at(4, 8) = -beam1;
    answer.at(4, 12) = -beam3;
}

void
EccBeam::computeBmatrixAt(GaussPoint* gp, FloatMatrix& answer, int li, int ui)
// eps = {\eps_x, \kappa_x, \kappa_y, \kappa_z}^T
{
    double l{ 0.0 }, rod{ 0.0 }, beam1{ 0.0 }, beam2{ 0.0 }, beam3{ 0.0 };
    TimeStep* tStep = this->domain->giveEngngModel()->giveCurrentStep();

    l = this->computeLength();
    rod = 1 / l;
    beam1 = (6 / (l*l)) * ((2 * gp->giveNaturalCoordinate(1)) / l - 1);
    beam2 = (2 / l) * (2 - 3 * gp->giveNaturalCoordinate(1) / l);
    beam3 = (2 / l) * (1 - 3 * gp->giveNaturalCoordinate(1) / l);

    answer.resize(4, 12);

    answer.at(1, 1) = -rod;
    answer.at(1, 7) = rod;

    answer.at(2, 4) = -rod;
    answer.at(2, 10) = rod;

    answer.at(3, 3) = beam1;
    answer.at(3, 5) = beam2;
    answer.at(3, 9) = -beam1;
    answer.at(3, 11) = beam3;

    answer.at(4, 2) = beam1;
    answer.at(4, 6) = -beam2;
    answer.at(4, 8) = -beam1;
    answer.at(4, 12) = -beam3;
}

void
EccBeam::computeConstitutiveMatrixAt(FloatMatrix& answer, MatResponseMode rMode, GaussPoint* gp, TimeStep* tStep)
{
    answer = this->giveStructuralCrossSection()->giveEBBeamStiffMtrx(rMode, gp, tStep);
}

void
EccBeam::computeGaussPoints()
{
    if (integrationRulesArray.size() == 0) {
        // the gauss point is used only when methods from crosssection and/or material
        // classes are requested
        integrationRulesArray.resize(1);
        integrationRulesArray[0] = std::make_unique< GaussIntegrationRule >(1, this, 1, 2);
        this->giveCrossSection()->setupIntegrationPoints(*integrationRulesArray[0], this->numberOfGaussPoints, this);
    }
}

bool
EccBeam::computeGtoLRotationMatrix(FloatMatrix& answer)
{
    FloatMatrix lcs;
    int ndofs = computeNumberOfGlobalDofs();
    answer.resize(ndofs, ndofs);
    answer.zero();

    this->giveLocalCoordinateSystem(lcs);
    for (int i = 1; i <= 3; i++) {
        for (int j = 1; j <= 3; j++) {
            answer.at(i, j) = lcs.at(i, j);
            answer.at(i + 3, j + 3) = lcs.at(i, j);
            answer.at(i + 6, j + 6) = lcs.at(i, j);
            answer.at(i + 9, j + 9) = lcs.at(i, j);
        }
    }

    for (int i = 13; i <= ndofs; i++) {
        answer.at(i, i) = 1.0;
    }

    return true;
}

double
EccBeam::computeLength()
{
    double dx, dy, dz, length = 0.0;
    Node* nodeA, * nodeB;

    if (length == 0.) {
        nodeA = this->giveNode(1);
        nodeB = this->giveNode(2);
        dx = nodeB->giveCoordinate(1) - nodeA->giveCoordinate(1);
        dy = nodeB->giveCoordinate(2) - nodeA->giveCoordinate(2);
        dz = nodeB->giveCoordinate(3) - nodeA->giveCoordinate(3);
        length = sqrt(dx * dx + dy * dy + dz * dz);
    }

    return length;
}

void
EccBeam::computeStiffnessMatrix(FloatMatrix& answer, MatResponseMode rMode, TimeStep* tStep)
{
    double L = this->computeLength();
    FloatMatrix D;
    this->computeConstitutiveMatrixAt(D, rMode, this->giveDefaultIntegrationRulePtr()->getIntegrationPoint(0), tStep);
    answer.resize(12, 12);
    answer.at(1, 1) = D.at(1, 1) / L;
    answer.at(1, 5) = -answer.at(1, 1) * offset_z;
    answer.at(1, 6) = -answer.at(1, 1) * offset_y;
    answer.at(1, 7) = - D.at(1, 1) / L;
    answer.at(1, 5) = -answer.at(1, 11) * offset_z;
    answer.at(1, 6) = -answer.at(1, 12) * offset_y;
    answer.at(2, 2) = 12 * D.at(4, 4) / (L * L * L);
    answer.at(2, 6) = 6 * D.at(4, 4) / (L * L);
    answer.at(2, 8) = - 12 * D.at(4, 4) / (L * L * L);
    answer.at(2, 12) = 6 * D.at(4, 4) / (L * L);
    answer.at(3, 3) = 12 * D.at(3, 3) / (L * L * L);
    answer.at(3, 5) = - 6 * D.at(3, 3) / (L * L);
    answer.at(3, 9) = - 12 * D.at(3, 3) / (L * L * L);
    answer.at(3, 11) = - 6 * D.at(3, 3) / (L * L);
    answer.at(4, 4) = D.at(2, 2) / L;
    answer.at(4, 10) = - D.at(2, 2) / L;
    answer.at(5, 5) = 4 * D.at(3, 3) / L - answer.at(1, 1) * offset_z*offset_z;
    answer.at(5, 7) = - D.at(1, 1) / L * offset_z;
    answer.at(5, 9) = 6 * D.at(3, 3) / (L * L);
    answer.at(5, 11) = 2 * D.at(3, 3) / L - answer.at(1, 1) * offset_z * offset_z;
    answer.at(6, 6) = 4 * D.at(4, 4) / L - answer.at(1, 1) * offset_y * offset_y;
    answer.at(6, 7) = -D.at(1, 1) / L * offset_y;
    answer.at(6, 8) = - 6 * D.at(4, 4) / (L * L);
    answer.at(6, 12) = 2 * D.at(4, 4) / L - answer.at(1, 1) * offset_y * offset_y;
    answer.at(7, 7) = D.at(1, 1) / L;
    answer.at(7, 11) = -D.at(1, 1) / L * offset_z;
    answer.at(7, 12) = -D.at(1, 1) / L * offset_y;
    answer.at(8, 8) = 12 * D.at(4, 4) / (L * L * L);
    answer.at(8, 12) = - 6 * D.at(4, 4) / (L * L);
    answer.at(9, 9) = 12 * D.at(3, 3) / (L * L * L);
    answer.at(9, 11) = 6 * D.at(3, 3) / (L * L);
    answer.at(10, 10) = D.at(2, 2) / L;
    answer.at(11, 11) = 4 * D.at(3, 3) / L - answer.at(1, 1) * offset_z * offset_z;
    answer.at(12, 12) = 4 * D.at(4, 4) / L - answer.at(1, 1) * offset_y * offset_y;
    answer.symmetrized();
}

void
EccBeam::computeStressVector(FloatArray& answer, const FloatArray& strain, GaussPoint* gp, TimeStep* tStep)
{
    answer = this->giveStructuralCrossSection()->giveGeneralizedStress_EBBeam(strain, gp, tStep);
}

double
EccBeam::computeVolumeAround(GaussPoint* gp)
{
    return gp->giveWeight() * 0.5 * this->computeLength();
}

void
EccBeam::giveDofManDofIDMask(int inode, IntArray& answer) const
{
    answer = {
        D_u, D_v, D_w, R_u, R_v, R_w
    };
}

void
EccBeam::giveInternalForcesVector(FloatArray& answer, TimeStep* tStep, int useUpdatedGpRecord)
{
    FloatArray u, strain, stress;
    FloatMatrix B;
    double L = this->computeLength();
    double x{ 0.0 };
    this->computeVectorOf(VM_Total, tStep, u);
    answer.clear();
    while (x <= L) {
        computeBmatrixAt(x, B);
        strain.beProductOf(B, u);
        this->computeStressVector(stress, strain, this->giveDefaultIntegrationRulePtr()->getIntegrationPoint(0), tStep);
        answer.addSubVector(stress, answer.giveSize() + 1);
        x += L;
    }

    FloatMatrix D;
    computeConstitutiveMatrixAt(D, ElasticStiffness, this->giveDefaultIntegrationRulePtr()->getIntegrationPoint(0), tStep);
    answer.at(2) = answer.at(8) = D.at(4, 4) / (L * L) * (-12 * u.at(8) / L + 12 * u.at(2) / L + 6 * u.at(12) + 6 * u.at(6));
    answer.at(3) = answer.at(9) = -D.at(3, 3) / (L * L) * (-12 * u.at(9) / L + 12 * u.at(3) / L - 6 * u.at(11) - 6 * u.at(5));
}

int
EccBeam::giveLocalCoordinateSystem(FloatMatrix& answer)
//
// returns a unit vectors of local coordinate system at element
// stored rowwise (mainly used by some materials with ortho and anisotrophy)
//
{
    FloatArray lx, ly, lz, help(3);
    Node* nodeA, * nodeB;
    nodeA = this->giveNode(1);
    nodeB = this->giveNode(2);

    lx.beDifferenceOf(nodeB->giveCoordinates(), nodeA->giveCoordinates());
    lx.normalize();

    Node* refNode = this->giveDomain()->giveNode(this->refNode);
    help.beDifferenceOf(refNode->giveCoordinates(), nodeA->giveCoordinates());

    lz.beVectorProductOf(lx, help);
    lz.normalize();

    ly.beVectorProductOf(lz, lx);
    ly.normalize();

    answer.resize(3, 3);
    answer.zero();
    for (int i = 1; i <= 3; i++) {
        answer.at(1, i) = lx.at(i);
        answer.at(2, i) = ly.at(i);
        answer.at(3, i) = lz.at(i);
    }

    return 1;
}

void
EccBeam::initializeFrom(InputRecord& ir)
{
    BeamBaseElement::initializeFrom(ir);
    
    if (ir.hasField(_IFT_EccBeam_refnode)) {
        IR_GIVE_FIELD(ir, refNode, _IFT_EccBeam_refnode);
        if (refNode == 0) {
            OOFEM_WARNING("Wrong reference node specified. Using default orientation.");
        }
    }
    
    if (ir.hasField(_IFT_EccBeam_offset_y))
        IR_GIVE_FIELD(ir, offset_y, _IFT_EccBeam_offset_y);

    if (ir.hasField(_IFT_EccBeam_offset_z))
        IR_GIVE_FIELD(ir, offset_z, _IFT_EccBeam_offset_z);
}

void
EccBeam::printOutputAt(FILE* File, TimeStep* tStep)
{
    FloatArray rl, Fl;

    fprintf(File, "beam element %d (%8d) :\n", this->giveLabel(), this->giveNumber());

    // ask for global element displacement vector
    this->computeVectorOf(VM_Total, tStep, rl);
    // ask for global element end forces vector
    this->giveInternalForcesVector(Fl, tStep);

    fprintf(File, "  Element forces    ");
    for (auto& val : Fl) {
        fprintf(File, " %.4e", val);
    }

    fprintf(File, "\n");
}

void
EccBeam::updateLocalNumbering(EntityRenumberingFunctor& f)
{
    BeamBaseElement::updateLocalNumbering(f);
    if (this->refNode) {
        this->refNode = f(this->refNode, ERS_DofManager);
    }
}
} // end namespace oofem