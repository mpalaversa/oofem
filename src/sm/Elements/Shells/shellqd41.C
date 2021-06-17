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
 
#include "sm/Elements/Shells/shellqd41.h"

#include "classfactory.h"
#include "gaussintegrationrule.h"
#include "gausspoint.h"
#include "load.h"

namespace oofem {
REGISTER_Element(ShellQd41);

ShellQd41 ::ShellQd41(int n, Domain* aDomain) : NLStructuralElement(n, aDomain)
{
    numberOfDofMans = 4;
    numberOfGaussPoints = 4;

    membrane = new LinQuad3DPlaneStress(n, aDomain);
    plate = new QDKTPlate(n, aDomain);
}
/*
void
ShellQd41::computeBmatrixAt(GaussPoint* gp, FloatMatrix& answer, int li, int ui)
{
    FloatMatrix BMatrixPlate;
    computeBmatrixPlateAt(gp, BMatrixPlate);

    FloatMatrix BMatrixMembrane;
    membrane->computeBmatrixAt(gp, BMatrixMembrane, li, ui);

    // Sets order of matrix B for ShellQd41 and all terms to 0.
    answer.resize(3, 24);
    answer.zero();
    // Set values for terms in the first row.
    answer.at(1, 1) = BMatrixMembrane(0, 0);
    answer.at(1, 3) = BMatrixPlate(0, 0);
    answer.at(1, 4) = BMatrixPlate(0, 1);
    answer.at(1, 5) = BMatrixPlate(0, 2);
    answer.at(1, 7) = BMatrixMembrane(0, 2);
    answer.at(1, 9) = BMatrixPlate(0, 3);
    answer.at(1, 10) = BMatrixPlate(0, 4);
    answer.at(1, 11) = BMatrixPlate(0, 5);
    answer.at(1, 13) = BMatrixMembrane(0, 4);
    answer.at(1, 15) = BMatrixPlate(0, 6);
    answer.at(1, 16) = BMatrixPlate(0, 7);
    answer.at(1, 17) = BMatrixPlate(0, 8);
    answer.at(1, 19) = BMatrixMembrane(0, 6);
    answer.at(1, 21) = BMatrixPlate(0, 9);
    answer.at(1, 22) = BMatrixPlate(0, 10);
    answer.at(1, 23) = BMatrixPlate(0, 11);
    // Set values for terms in the second row.
    answer.at(2, 2) = BMatrixMembrane(1, 1);
    answer.at(2, 3) = BMatrixPlate(1, 0);
    answer.at(2, 4) = BMatrixPlate(1, 1);
    answer.at(2, 5) = BMatrixPlate(1, 2);
    answer.at(2, 8) = BMatrixMembrane(1, 3);
    answer.at(2, 9) = BMatrixPlate(1, 3);
    answer.at(2, 10) = BMatrixPlate(1, 4);
    answer.at(2, 11) = BMatrixPlate(1, 5);
    answer.at(2, 14) = BMatrixMembrane(1, 5);
    answer.at(2, 15) = BMatrixPlate(1, 6);
    answer.at(2, 16) = BMatrixPlate(1, 7);
    answer.at(2, 17) = BMatrixPlate(1, 8);
    answer.at(2, 20) = BMatrixMembrane(1, 7);
    answer.at(2, 21) = BMatrixPlate(1, 9);
    answer.at(2, 22) = BMatrixPlate(1, 10);
    answer.at(2, 23) = BMatrixPlate(1, 11);
    // Set values for terms in the third row.
    answer.at(3, 1) = BMatrixMembrane(2, 0);
    answer.at(3, 2) = BMatrixMembrane(2, 1);
    answer.at(3, 3) = BMatrixPlate(2, 0);
    answer.at(3, 4) = BMatrixPlate(2, 1);
    answer.at(3, 5) = BMatrixPlate(2, 2);
    answer.at(3, 7) = BMatrixMembrane(2, 2);
    answer.at(3, 8) = BMatrixMembrane(2, 3);
    answer.at(3, 9) = BMatrixPlate(2, 3);
    answer.at(3, 10) = BMatrixPlate(2, 4);
    answer.at(3, 11) = BMatrixPlate(2, 5);
    answer.at(3, 13) = BMatrixMembrane(2, 4);
    answer.at(3, 14) = BMatrixMembrane(2, 5);
    answer.at(3, 15) = BMatrixPlate(2, 6);
    answer.at(3, 16) = BMatrixPlate(2, 7);
    answer.at(3, 17) = BMatrixPlate(2, 8);
    answer.at(3, 19) = BMatrixMembrane(2, 6);
    answer.at(3, 20) = BMatrixMembrane(2, 7);
    answer.at(3, 21) = BMatrixPlate(2, 9);
    answer.at(3, 22) = BMatrixPlate(2, 10);
    answer.at(3, 23) = BMatrixPlate(2, 11);
}
*/

void
ShellQd41::computeBmatrixAt(GaussPoint* gp, FloatMatrix& answer, int li, int ui)
{
    FloatMatrix BMatrixMembrane;
    membrane->computeBmatrixAt(gp, BMatrixMembrane, li, ui);

    // Sets order of matrix B for ShellQd41 and all terms to 0.
    answer.resize(3, 24);
    answer.zero();
    // Set values for terms in the first row.
    answer.at(1, 1) = BMatrixMembrane(0, 0);
    answer.at(1, 7) = BMatrixMembrane(0, 2);
    answer.at(1, 13) = BMatrixMembrane(0, 4);
    answer.at(1, 19) = BMatrixMembrane(0, 6);
    // Set values for terms in the second row.
    answer.at(2, 2) = BMatrixMembrane(1, 1);
    answer.at(2, 8) = BMatrixMembrane(1, 3);
    answer.at(2, 14) = BMatrixMembrane(1, 5);
    answer.at(2, 20) = BMatrixMembrane(1, 7);
    // Set values for terms in the third row.
    answer.at(3, 1) = BMatrixMembrane(2, 0);
    answer.at(3, 2) = BMatrixMembrane(2, 1);
    answer.at(3, 7) = BMatrixMembrane(2, 2);
    answer.at(3, 8) = BMatrixMembrane(2, 3);
    answer.at(3, 13) = BMatrixMembrane(2, 4);
    answer.at(3, 14) = BMatrixMembrane(2, 5);
    answer.at(3, 19) = BMatrixMembrane(2, 6);
    answer.at(3, 20) = BMatrixMembrane(2, 7);
}

void
ShellQd41::computeBmatrixPlateAt(double ksi, double eta, FloatMatrix& answer)
// Returns the [5x12] strain-displacement matrix {B} of the receiver,
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

    answer.resize(5, 12);
    answer.zero();

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

    // Note: no shear strains, no shear forces => the 4th and 5th rows are zero
}

void
ShellQd41::computeBmatrixPlateAt(GaussPoint* gp, FloatMatrix& answer) {
    computeBmatrixPlateAt(gp->giveNaturalCoordinate(1), gp->giveNaturalCoordinate(2), answer);
}

void
ShellQd41::computeBodyLoadVectorAt(FloatArray& answer, Load* forLoad, TimeStep* tStep, ValueModeType mode) {
    double dens, dV;
    FloatArray force, ntf, loadVectorFromPlate;
    FloatMatrix T;

    if ((forLoad->giveBCGeoType() != BodyLoadBGT) || (forLoad->giveBCValType() != ForceLoadBVT)) {
        OOFEM_ERROR("unknown load type");
    }

    // note: force is assumed to be in global coordinate system.
    forLoad->computeComponentArrayAt(force, tStep, mode);
    // transform from global to element local c.s
    if (this->computeLoadGToLRotationMtrx(T)) {
        force.rotatedWith(T, 'n');
    }
    FloatArray loadVector;
    loadVector.resize(3);
    int j{ 1 };
    for(int i = 1; i <= force.giveSize(); i++)
        if (force.at(i) != double(0)) {
            loadVector.at(1) = force.at(i);
            j = i;
        }

    loadVectorFromPlate.clear();
    FloatArray NMatrixTemp;
    FloatMatrix NMatrix;
    if (loadVector.giveSize()) {
        for (GaussPoint* gp : *this->giveDefaultIntegrationRulePtr()) {
            giveInterpolation()->evalN(NMatrixTemp, gp->giveSubPatchCoordinates(), *membrane->giveCellGeometryWrapper());
            NMatrix.beNMatrixOf(NMatrixTemp, 3);
            dV = computePlateVolumeAround(gp) * this->giveCrossSection()->give(CS_Thickness, gp);
            dens = this->giveCrossSection()->give('d', gp);
            ntf.beTProductOf(NMatrix, loadVector);
            loadVectorFromPlate.add(dV * dens, ntf);
        }
    }
    else {
        return;
    }

    answer.resize(24);
    answer.zero();
    for (int i = 1; i <= 12; i+=3) {
        // takes into account only local w-displacement (add \theta_x and \theta_y in the future)
        answer.at(j) = loadVectorFromPlate.at(i);
        //answer.at(j+1) = loadVectorFromPlate.at(i+1);
        //answer.at(j+2) = loadVectorFromPlate.at(i+2);
        j += 6;
    }
}

void
ShellQd41::computeConstitutiveMatrixAt(FloatMatrix& answer, MatResponseMode rMode, GaussPoint* gp, TimeStep* tStep)
{   
    answer.add(this->giveStructuralCrossSection()->giveStiffnessMatrix_PlaneStress(rMode, gp, tStep));
    FloatMatrix DPlate;
    DPlate.beSubMatrixOf(this->giveStructuralCrossSection()->give2dPlateStiffMtrx(rMode, gp, tStep), 1, 3, 1, 3);
    answer.add(DPlate);
}

void
ShellQd41::computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    if (integrationRulesArray.size() == 0) {
        integrationRulesArray.resize(1);
        integrationRulesArray[0] = std::make_unique<GaussIntegrationRule>(1, this, 1, 5);
        this->giveCrossSection()->setupIntegrationPoints(*integrationRulesArray[0], numberOfGaussPoints, this);
    }
}

bool
ShellQd41::computeGtoLRotationMatrix(FloatMatrix& answer)
// Returns the rotation matrix of the receiver of the size [24,24]
// r(local) = T * r(global)
// for one node (r written transposed): {u,v,w,R_u,R_v,R_w} = T * {u,v,w,R_u,R_v,R_w}
{
    if (membrane->GtoLRotationMatrix == NULL)
        membrane->computeGtoLRotationMatrix();;

    answer.resize(24, 24);
    answer.zero();

    for (int i = 1; i <= 3; i++) {
        answer.at(1, i) = answer.at(4, i + 3) = answer.at(7, i + 6) = answer.at(10, i + 9) = answer.at(13, i + 12) = answer.at(16, i + 15) = answer.at(19, i + 18) = answer.at(22, i + 21) = membrane->GtoLRotationMatrix->at(1, i);
        answer.at(2, i) = answer.at(5, i + 3) = answer.at(8, i + 6) = answer.at(11, i + 9) = answer.at(14, i + 12) = answer.at(17, i + 15) = answer.at(20, i + 18) = answer.at(23, i + 21) = membrane->GtoLRotationMatrix->at(2, i);
        answer.at(3, i) = answer.at(6, i + 3) = answer.at(9, i + 6) = answer.at(12, i + 9) = answer.at(15, i + 12) = answer.at(18, i + 15) = answer.at(21, i + 18) = answer.at(24, i + 21) = membrane->GtoLRotationMatrix->at(3, i);
    }

    return 1;
}

void
ShellQd41::computeMembraneStrainVectorAt(FloatArray& answer, double xi, double eta, TimeStep* tStep) {
    FloatMatrix b;
    FloatArray u;
    
    FloatArray uMemb;
    int j = 1;
    switch (outputType) {
    case OutputType::Standard:
        this->computeVectorOf(VM_Total, tStep, u);
        /* This is to be uncommented once adapted for ShellQd41 (if necessary).
        if (initialDisplacements) {
            u.subtract(*initialDisplacements);
        }*/
        uMemb.resize(8);
        for (int i = 1; i <= 8; i += 2) {
            uMemb.at(i) = u.at(j);
            uMemb.at(i + 1) = u.at(j + 1);
            j += 6;
        }
        membrane->computeBmatrixAt(xi, eta, b);
        answer.beProductOf(b, uMemb);
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
ShellQd41::computePlateCurvaturesAt(FloatArray& answer, double xi, double eta, TimeStep* tStep) {
    FloatArray u;
    this->computeVectorOf(VM_Total, tStep, u);
    /* This is to be uncommented once adapted for ShellQd41 (if necessary).
    if (initialDisplacements) {
        u.subtract(*initialDisplacements);
    }*/

    FloatArray uPlate;
    uPlate.resize(12);
    int j = 3;
    for (int i = 1; i <= 12; i += 3) {
        uPlate.at(i) = u.at(j);
        uPlate.at(i + 1) = u.at(j + 1);
        uPlate.at(i + 2) = u.at(j + 2);
        j += 6;
    }

    FloatMatrix temp;
    computeBmatrixPlateAt(xi, eta, temp);
    FloatMatrix b;
    b.beSubMatrixOf(temp, 1, 3, 1, 12);
    answer.beProductOf(b, uPlate);
}

void
ShellQd41::computePlateStrainVectorAt(FloatArray& answer, double xi, double eta, TimeStep* tStep) {
    FloatArray curvatures;
    switch (outputType) {
    case OutputType::Standard:
        computePlateCurvaturesAt(curvatures, xi, eta, tStep);
        curvatures.at(3) = 2 * curvatures.at(3);
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

double
ShellQd41::computePlateVolumeAround(GaussPoint* gp)
// Returns the portion of the receiver which is attached to gp.
{
    double detJ, weight;

    weight = gp->giveWeight();
    detJ = fabs(giveInterpolation()->giveTransformationJacobian(gp->giveNaturalCoordinates(), *membrane->giveCellGeometryWrapper()));
    return detJ * weight; ///@todo What about thickness?
}

void
ShellQd41::computeStiffnessMatrix(FloatMatrix& answer, MatResponseMode rMode, TimeStep* tStep)
{
    StructuralCrossSection* cs = this->giveStructuralCrossSection();

    auto tempDrillCoeff = cs->give(CS_RelDrillingStiffness, this->giveDefaultIntegrationRulePtr()->getIntegrationPoint(0));
    if (tempDrillCoeff != 0.0)
        drillCoeff = tempDrillCoeff;
    /*
    IntArray membDofManArray = membrane->giveDofManArray();
    for (int i = 0; i <= dofManArray.giveSize(); i++) {
        membDofManArray[i] = dofManArray.at(i);
    }*/

    FloatMatrix BMatrixPlate;
    FloatMatrix BMatrixMembrane;

    FloatMatrix DMatrixPlate;
    FloatMatrix DMatrixMembrane;

    FloatMatrix DBPlate;
    FloatMatrix DBMembrane;

    FloatMatrix stiffMatPlate;
    FloatMatrix stiffMatMembrane;

    FloatMatrix temp;

    bool matStiffSymmFlag = cs->isCharacteristicMtrxSymmetric(rMode);

    answer.clear();

    if (!this->isActivated(tStep)) {
        return;
    }

    // Compute matrix from material stiffness (total stiffness for small def.) - B^T * dS/dE * B
    if (integrationRulesArray.size() == 1) {
        for (auto& gp : *this->giveDefaultIntegrationRulePtr()) {
            // Engineering (small strain) stiffness
            if (nlGeometry == 0) {
                computeBmatrixPlateAt(gp, temp);
                BMatrixPlate.beSubMatrixOf(temp, 1, 3, 1, 12);
                membrane->computeBmatrixAt(gp, BMatrixMembrane);

                DMatrixPlate = cs->giveKirchhoffPlateStiffMtrx(rMode, gp, tStep);
                DMatrixMembrane = cs->giveStiffnessMatrix_PlaneStress(rMode, gp, tStep);
            }
            else if (nlGeometry == 1)
                OOFEM_ERROR("Only linear geometry is allowed for this element.");

            DBPlate.beProductOf(DMatrixPlate, BMatrixPlate);
            DBMembrane.beProductOf(DMatrixMembrane, BMatrixMembrane);

            double dV = computeVolumeAround(gp);

            if (matStiffSymmFlag) {
                stiffMatPlate.plusProductSymmUpper(BMatrixPlate, DBPlate, this->computePlateVolumeAround(gp));
                stiffMatMembrane.plusProductSymmUpper(BMatrixMembrane, DBMembrane, dV);
            }
            else {
                OOFEM_ERROR("Stiffness matrix of this element must be symmetric.");
            }
        }

        answer.resize(24, 24);
        answer.zero();

        // Assembling element stiffness matrix composed of the "membrane" and the "plate" part.
        // 1st row/column:
        answer.at(1, 1) = stiffMatMembrane.at(1, 1);
        answer.at(1, 2) = answer.at(2, 1) = stiffMatMembrane.at(1, 2);
        answer.at(1, 7) = answer.at(7, 1) = stiffMatMembrane.at(1, 3);
        answer.at(1, 8) = answer.at(8, 1) = stiffMatMembrane.at(1, 4);
        answer.at(1, 13) = answer.at(13, 1) = stiffMatMembrane.at(1, 5);
        answer.at(1, 14) = answer.at(14, 1) = stiffMatMembrane.at(1, 6);
        answer.at(1, 19) = answer.at(19, 1) = stiffMatMembrane.at(1, 7);
        answer.at(1, 20) = answer.at(20, 1) = stiffMatMembrane.at(1, 8);
        // 2nd row/column:
        answer.at(2, 2) = stiffMatMembrane.at(2, 2);
        answer.at(2, 7) = answer.at(7, 2) = stiffMatMembrane.at(2, 3);
        answer.at(2, 8) = answer.at(8, 2) = stiffMatMembrane.at(2, 4);
        answer.at(2, 13) = answer.at(13, 2) = stiffMatMembrane.at(2, 5);
        answer.at(2, 14) = answer.at(14, 2) = stiffMatMembrane.at(2, 6);
        answer.at(2, 19) = answer.at(19, 2) = stiffMatMembrane.at(2, 7);
        answer.at(2, 20) = answer.at(20, 2) = stiffMatMembrane.at(2, 8);
        // 3rd row/column:
        answer.at(3, 3) = stiffMatPlate.at(1, 1);
        answer.at(3, 4) = answer.at(4, 3) = stiffMatPlate.at(1, 2);
        answer.at(3, 5) = answer.at(5, 3) = stiffMatPlate.at(1, 3);
        answer.at(3, 9) = answer.at(9, 3) = stiffMatPlate.at(1, 4);
        answer.at(3, 10) = answer.at(10, 3) = stiffMatPlate.at(1, 5);
        answer.at(3, 11) = answer.at(11, 3) = stiffMatPlate.at(1, 6);
        answer.at(3, 15) = answer.at(15, 3) = stiffMatPlate.at(1, 7);
        answer.at(3, 16) = answer.at(16, 3) = stiffMatPlate.at(1, 8);
        answer.at(3, 17) = answer.at(17, 3) = stiffMatPlate.at(1, 9);
        answer.at(3, 21) = answer.at(21, 3) = stiffMatPlate.at(1, 10);
        answer.at(3, 22) = answer.at(22, 3) = stiffMatPlate.at(1, 11);
        answer.at(3, 23) = answer.at(23, 3) = stiffMatPlate.at(1, 12);
        // 4th row/column:
        answer.at(4, 4) = stiffMatPlate.at(2, 2);
        answer.at(4, 5) = answer.at(5, 4) = stiffMatPlate.at(2, 3);
        answer.at(4, 9) = answer.at(9, 4) = stiffMatPlate.at(2, 4);
        answer.at(4, 10) = answer.at(10, 4) = stiffMatPlate.at(2, 5);
        answer.at(4, 11) = answer.at(11, 4) = stiffMatPlate.at(2, 6);
        answer.at(4, 15) = answer.at(15, 4) = stiffMatPlate.at(2, 7);
        answer.at(4, 16) = answer.at(16, 4) = stiffMatPlate.at(2, 8);
        answer.at(4, 17) = answer.at(17, 4) = stiffMatPlate.at(2, 9);
        answer.at(4, 21) = answer.at(21, 4) = stiffMatPlate.at(2, 10);
        answer.at(4, 22) = answer.at(22, 4) = stiffMatPlate.at(2, 11);
        answer.at(4, 23) = answer.at(23, 4) = stiffMatPlate.at(2, 12);
        // 5th row/column:
        answer.at(5, 5) = stiffMatPlate.at(3, 3);
        answer.at(5, 9) = answer.at(9, 5) = stiffMatPlate.at(3, 4);
        answer.at(5, 10) = answer.at(10, 5) = stiffMatPlate.at(3, 5);
        answer.at(5, 11) = answer.at(11, 5) = stiffMatPlate.at(3, 6);
        answer.at(5, 15) = answer.at(15, 5) = stiffMatPlate.at(3, 7);
        answer.at(5, 16) = answer.at(16, 5) = stiffMatPlate.at(3, 8);
        answer.at(5, 17) = answer.at(17, 5) = stiffMatPlate.at(3, 9);
        answer.at(5, 21) = answer.at(21, 5) = stiffMatPlate.at(3, 10);
        answer.at(5, 22) = answer.at(22, 5) = stiffMatPlate.at(3, 11);
        answer.at(5, 23) = answer.at(23, 5) = stiffMatPlate.at(3, 12);
        // 6th row/column:
        answer.at(6, 6) = drillCoeff;
        // 7th row/column:
        answer.at(7, 7) = stiffMatMembrane.at(3, 3);
        answer.at(7, 8) = answer.at(8, 7) = stiffMatMembrane.at(3, 4);
        answer.at(7, 13) = answer.at(13, 7) = stiffMatMembrane.at(3, 5);
        answer.at(7, 14) = answer.at(14, 7) = stiffMatMembrane.at(3, 6);
        answer.at(7, 19) = answer.at(19, 7) = stiffMatMembrane.at(3, 7);
        answer.at(7, 20) = answer.at(20, 7) = stiffMatMembrane.at(3, 8);
        // 8th row/column:
        answer.at(8, 8) = stiffMatMembrane.at(4, 4);
        answer.at(8, 13) = answer.at(13, 8) = stiffMatMembrane.at(4, 5);
        answer.at(8, 14) = answer.at(14, 8) = stiffMatMembrane.at(4, 6);
        answer.at(8, 19) = answer.at(19, 8) = stiffMatMembrane.at(4, 7);
        answer.at(8, 20) = answer.at(20, 8) = stiffMatMembrane.at(4, 8);
        // 9th row/column:
        answer.at(9, 9) = stiffMatPlate.at(4, 4);
        answer.at(9, 10) = answer.at(10, 9) = stiffMatPlate.at(4, 5);
        answer.at(9, 11) = answer.at(11, 9) = stiffMatPlate.at(4, 6);
        answer.at(9, 15) = answer.at(15, 9) = stiffMatPlate.at(4, 7);
        answer.at(9, 16) = answer.at(16, 9) = stiffMatPlate.at(4, 8);
        answer.at(9, 17) = answer.at(17, 9) = stiffMatPlate.at(4, 9);
        answer.at(9, 21) = answer.at(21, 9) = stiffMatPlate.at(4, 10);
        answer.at(9, 22) = answer.at(22, 9) = stiffMatPlate.at(4, 11);
        answer.at(9, 23) = answer.at(23, 9) = stiffMatPlate.at(4, 12);
        // 10th row/column:
        answer.at(10, 10) = stiffMatPlate.at(5, 5);
        answer.at(10, 11) = answer.at(11, 10) = stiffMatPlate.at(5, 6);
        answer.at(10, 15) = answer.at(15, 10) = stiffMatPlate.at(5, 7);
        answer.at(10, 16) = answer.at(16, 10) = stiffMatPlate.at(5, 8);
        answer.at(10, 17) = answer.at(17, 10) = stiffMatPlate.at(5, 9);
        answer.at(10, 21) = answer.at(21, 10) = stiffMatPlate.at(5, 10);
        answer.at(10, 22) = answer.at(22, 10) = stiffMatPlate.at(5, 11);
        answer.at(10, 23) = answer.at(23, 10) = stiffMatPlate.at(5, 12);
        // 11th row/column:
        answer.at(11, 11) = stiffMatPlate.at(6, 6);
        answer.at(11, 15) = answer.at(15, 11) = stiffMatPlate.at(6, 7);
        answer.at(11, 16) = answer.at(16, 11) = stiffMatPlate.at(6, 8);
        answer.at(11, 17) = answer.at(17, 11) = stiffMatPlate.at(6, 9);
        answer.at(11, 21) = answer.at(21, 11) = stiffMatPlate.at(6, 10);
        answer.at(11, 22) = answer.at(22, 11) = stiffMatPlate.at(6, 11);
        answer.at(11, 23) = answer.at(23, 11) = stiffMatPlate.at(6, 12);
        // 12th row/column:
        answer.at(12, 12) = drillCoeff;
        // 13th row/column:
        answer.at(13, 13) = stiffMatMembrane.at(5, 5);
        answer.at(13, 14) = answer.at(14, 13) = stiffMatMembrane.at(5, 6);
        answer.at(13, 19) = answer.at(19, 13) = stiffMatMembrane.at(5, 7);
        answer.at(13, 20) = answer.at(20, 13) = stiffMatMembrane.at(5, 8);
        // 14th row/column:
        answer.at(14, 14) = stiffMatMembrane.at(6, 6);
        answer.at(14, 19) = answer.at(19, 14) = stiffMatMembrane.at(6, 7);
        answer.at(14, 20) = answer.at(20, 14) = stiffMatMembrane.at(6, 8);
        // 15th row/column:
        answer.at(15, 15) = stiffMatPlate.at(7, 7);
        answer.at(15, 16) = answer.at(16, 15) = stiffMatPlate.at(7, 8);
        answer.at(15, 17) = answer.at(17, 15) = stiffMatPlate.at(7, 9);
        answer.at(15, 21) = answer.at(21, 15) = stiffMatPlate.at(7, 10);
        answer.at(15, 22) = answer.at(22, 15) = stiffMatPlate.at(7, 11);
        answer.at(15, 23) = answer.at(23, 15) = stiffMatPlate.at(7, 12);
        // 16th row/column:
        answer.at(16, 16) = stiffMatPlate.at(8, 8);
        answer.at(16, 17) = answer.at(17, 16) = stiffMatPlate.at(8, 9);
        answer.at(16, 21) = answer.at(21, 16) = stiffMatPlate.at(8, 10);
        answer.at(16, 22) = answer.at(22, 16) = stiffMatPlate.at(8, 11);
        answer.at(16, 23) = answer.at(23, 16) = stiffMatPlate.at(8, 12);
        // 17th row/column:
        answer.at(17, 17) = stiffMatPlate.at(9, 9);
        answer.at(17, 21) = answer.at(21, 17) = stiffMatPlate.at(9, 10);
        answer.at(17, 22) = answer.at(22, 17) = stiffMatPlate.at(9, 11);
        answer.at(17, 23) = answer.at(23, 17) = stiffMatPlate.at(9, 12);
        // 18th row/column:
        answer.at(18, 18) = drillCoeff;
        // 19th row/column:
        answer.at(19, 19) = stiffMatMembrane.at(7, 7);
        answer.at(19, 20) = answer.at(20, 19) = stiffMatMembrane.at(7, 8);
        // 20th row/column:
        answer.at(20, 20) = stiffMatMembrane.at(8, 8);
        // 21st row/column:
        answer.at(21, 21) = stiffMatPlate.at(10, 10);
        answer.at(21, 22) = answer.at(22, 21) = stiffMatPlate.at(10, 11);
        answer.at(21, 23) = answer.at(23, 21) = stiffMatPlate.at(10, 12);
        // 22nd row/column:
        answer.at(22, 22) = stiffMatPlate.at(11, 11);
        answer.at(22, 23) = answer.at(23, 22) = stiffMatPlate.at(11, 12);
        // 23rd row/column:
        answer.at(23, 23) = stiffMatPlate.at(12, 12);
        // 24th row/column:
        answer.at(24, 24) = drillCoeff;
    }
    else
        OOFEM_ERROR("It is not allowed to define more than 1 integration rule for this element.");
}

void
ShellQd41::computeStrainVector(FloatArray& answer, GaussPoint* gp, TimeStep* tStep)
// Computes the vector containing the strains at the Gauss point gp of
// the receiver, at time step tStep.
{
    if (!this->isActivated(tStep)) {
        answer.resize(StructuralMaterial::giveSizeOfVoigtSymVector(gp->giveMaterialMode()));
        answer.zero();
        return;
    }

    switch (outputCategory) {
    case OutputCategory::Membrane:
        computeMembraneStrainVectorAt(answer, gp->giveNaturalCoordinate(1), gp->giveNaturalCoordinate(2), tStep);
        break;
    case OutputCategory::Plate:
        computePlateStrainVectorAt(answer, gp->giveNaturalCoordinate(1), gp->giveNaturalCoordinate(2), tStep);
        break;
    case OutputCategory::Combined:
        OOFEM_ERROR("Not yet implemented for a %s element.", giveClassName());
        break;
    default:
        OOFEM_ERROR("Something went wrong. An unknown output category requested for element %d.", giveGlobalNumber());
        break;
    }
}

void
ShellQd41::computeStrainVectorAtCentroid(FloatArray& answer, TimeStep* tStep) {
    FloatArray membraneStrains, plateStrains;
    
    switch (outputCategory) {
    case OutputCategory::Membrane:
        computeMembraneStrainVectorAt(answer, 0.0, 0.0, tStep);
        break;
    case OutputCategory::Plate:
        computePlateStrainVectorAt(answer, 0.0, 0.0, tStep);
        break;
    case OutputCategory::Combined:
        computeMembraneStrainVectorAt(membraneStrains, 0.0, 0.0, tStep);
        computePlateStrainVectorAt(plateStrains, 0.0, 0.0, tStep);
        membraneStrains.add(plateStrains);
        answer.resize(3);
        answer = membraneStrains;
        break;
    default:
        OOFEM_ERROR("Something went wrong. An unknown output category requested for element %d.", giveGlobalNumber());
        break;
    }
}

void
ShellQd41::computeStressVector(FloatArray& answer, const FloatArray& strain, GaussPoint* gp, TimeStep* tStep) {
    switch (outputCategory) {
    case OutputCategory::Membrane:
        answer = this->giveStructuralCrossSection()->giveRealStress_PlaneStress(strain, gp, tStep);
        break;
    case OutputCategory::Plate:
        answer = this->giveStructuralCrossSection()->giveRealStress_KirchhoffPlate(strain, gp, tStep);
        break;
    case OutputCategory::Combined:
        OOFEM_ERROR("Not yet implemented for a %s element.", giveClassName());
        break;
    default:
        OOFEM_ERROR("Something went wrong. An unknown output category requested for element %d.", giveGlobalNumber());
        break;
    }
}

void
ShellQd41::computeStressVectorAtCentroid(FloatArray& answer, TimeStep* tStep, const FloatArray& strain) {
    FloatArray membraneStrains, plateStrains, membraneStresses, plateStresses;
    
    switch (outputCategory) {
    case OutputCategory::Membrane:
        answer = this->giveStructuralCrossSection()->giveRealStress_PlaneStress(strain, this->giveIntegrationRulesArray()[0]->getIntegrationPoint(0), tStep);
        break;
    case OutputCategory::Plate:
        plateStresses = this->giveStructuralCrossSection()->giveRealStress_KirchhoffPlate(strain, this->giveIntegrationRulesArray()[0]->getIntegrationPoint(0), tStep);
        answer.resize(6);
        for (int i = 1; i <= 6; i++) {
            if (i < 4)
                answer.at(i) = plateStresses.at(i);
            else
                answer.at(i) = -plateStresses.at(i - 3);
        }
        break;
    case OutputCategory::Combined:
        computeMembraneStrainVectorAt(membraneStrains, 0.0, 0.0, tStep);
        membraneStresses = this->giveStructuralCrossSection()->giveRealStress_PlaneStress(membraneStrains, this->giveIntegrationRulesArray()[0]->getIntegrationPoint(0), tStep);;
        computePlateStrainVectorAt(plateStrains, 0.0, 0.0, tStep);
        plateStresses = this->giveStructuralCrossSection()->giveRealStress_KirchhoffPlate(plateStrains, this->giveIntegrationRulesArray()[0]->getIntegrationPoint(0), tStep);

        answer.resize(6);
        for (int i = 1; i <= 6; i++) {
            if (i < 4)
                answer.at(i) = plateStresses.at(i) + membraneStresses.at(i);
            else
                answer.at(i) = -plateStresses.at(i - 3) + membraneStresses.at(i - 3);
        }
        break;
    default:
        OOFEM_ERROR("Something went wrong. An unknown output category requested for element %d.", giveGlobalNumber());
        break;
    }
}

void
ShellQd41::getStressesTopBottom(FloatArray& answer, TimeStep* tStep) {
    outputAtXY = OutputLocationXY::Centroid;
    outputAtZ = this->giveStructuralCrossSection()->give(CS_Thickness, this->giveDefaultIntegrationRulePtr()->getIntegrationPoint(0))/2;
    outputCategory = OutputCategory::Combined;
    outputType = OutputType::Standard;

    computeStressVectorAtCentroid(answer, tStep);
}

void
ShellQd41::giveSurfaceDofMapping(IntArray& answer, int iSurf) const
{
    if (iSurf == 1 || iSurf == 2) {
        answer.enumerate(24);
    }
    else {
        OOFEM_ERROR("wrong surface number");
    }
}

void
ShellQd41::computeSurfaceNMatrix(FloatMatrix& answer, int boundaryID, const FloatArray& lcoords)
{
    FloatArray n_vec;
    this->giveInterpolation()->boundarySurfaceEvalN(n_vec, boundaryID, lcoords, *membrane->giveCellGeometryWrapper());
    answer.beNMatrixOf(n_vec, 6);
}
/*
double
ShellQd41::computeSurfaceVolumeAround(GaussPoint* gp, int iSurf)
{
    FEInterpolation* fei = this->giveInterpolation();
    const FloatArray& lcoords = gp->giveNaturalCoordinates();
    double J = fei->boundarySurfaceGiveTransformationJacobian(iSurf, lcoords, *membrane->giveCellGeometryWrapper());

    return (gp->giveWeight() * J);
}
*/
std::vector< FloatArray >
ShellQd41::giveNodeCoordinates()
{
    IntArray newDofMans(this->dofManArray.giveSize());
    for (int i = 1; i <= this->dofManArray.giveSize(); i++)
        newDofMans.at(i) = this->dofManArray.at(i);
    
    membrane->setDofManagers(newDofMans);

    std::vector< FloatArray > c;
    membrane->computeLocalNodalCoordinates(c);
    return c;
}

double
ShellQd41::computeVolumeAround(GaussPoint* gp)
{
    return membrane->computeVolumeAround(gp);
}

void
ShellQd41::giveDofManDofIDMask(int inode, IntArray& answer) const
{
    answer = { D_u, D_v, D_w, R_u, R_v, R_w };
}

void
ShellQd41::initializeFrom(InputRecord& ir)
{
    StructuralElement::initializeFrom(ir);
    plate->initializeFrom(ir);
    membrane->initializeFrom(ir);

    int outputAtXYTemp, outputTypeTemp, outputAtZTemp, outputCategoryTemp;
    outputAtXYTemp = outputTypeTemp = outputAtZTemp = outputCategoryTemp = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, outputAtXYTemp, _IFT_ShellQd41_outputAtXY);
    IR_GIVE_OPTIONAL_FIELD(ir, outputTypeTemp, _IFT_ShellQd41_outputType);
    IR_GIVE_OPTIONAL_FIELD(ir, outputCategoryTemp, _IFT_ShellQd41_outputCategory);
    IR_GIVE_OPTIONAL_FIELD(ir, outputAtZ, _IFT_ShellQd41_outputAtZ);

    switch (outputAtXYTemp) {
    case 1:
        outputAtXY = OutputLocationXY::GaussPoints;
        break;
    case 2:
        outputAtXY = OutputLocationXY::Centroid;
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
    switch (outputCategoryTemp) {
    case 1:
        outputCategory = OutputCategory::Membrane;
        break;
    case 2:
        outputCategory = OutputCategory::Plate;
        break;
    case 3:
        outputCategory = OutputCategory::Combined;
        break;
    case 4:
        outputCategory = OutputCategory::All;
        break;
    default:
        outputCategory = OutputCategory::Membrane;
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
ShellQd41::setCrossSection(int csIndx)
{
    StructuralElement::setCrossSection(csIndx);
    plate->setCrossSection(csIndx);
    membrane->setCrossSection(csIndx);
}

void
ShellQd41::updateInternalState(TimeStep* tStep)
// Updates receiver at the end of the step.
{
    FloatArray stress, strain;
    // force updating strains & stresses
    switch (outputAtXY) {
    case OutputLocationXY::GaussPoints:
        for (auto& iRule : integrationRulesArray) {
            for (GaussPoint* gp : *iRule) {
                this->computeStrainVector(strain, gp, tStep);
                this->computeStressVector(stress, strain, gp, tStep);
            }
        }
        break;
    case OutputLocationXY::Centroid:
        this->computeStrainVectorAtCentroid(strain, tStep);
        this->computeStressVectorAtCentroid(stress, tStep, strain);
        break;
    case OutputLocationXY::Corners:
        OOFEM_ERROR("Not implemented for this element");
        break;
    case OutputLocationXY::All:
        OOFEM_ERROR("Not implemented for this element");
        break;
    deafult:
        OOFEM_ERROR("Something went wrong. Output requested at an unknown location in x-y plane of element %d.", giveGlobalNumber());
        break;
    }
}

void
ShellQd41::updateLocalNumbering(EntityRenumberingFunctor& f)
{
    // Update numbering of the related DOFs for the membrane and plate part.
    membrane->updateLocalNumbering(f);
    plate->updateLocalNumbering(f);
    // Update numbering of the related DOFs for ShellQd41.
    for (auto& dnum : dofManArray)
        dnum = f(dnum, ERS_DofManager);
}
}