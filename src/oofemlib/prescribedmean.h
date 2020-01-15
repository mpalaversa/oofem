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
#ifndef PRESCRIBEDMEAN_H
#define PRESCRIBEDMEAN_H

#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <memory>

#include "activebc.h"
#include "node.h"
#include "boundaryload.h"

namespace oofem
{

///@name Input fields for PrescribedMean
//@{
#define _IFT_PrescribedMean_Name "prescribedmean"
#define _IFT_PrescribedMean_DofID "dofid"
#define _IFT_PrescribedMean_Mean "mean"
#define _IFT_PrescribedMean_Edge "edge"

//@}

class OOFEM_EXPORT PrescribedMean : public ActiveBoundaryCondition
{
private:
    std::unique_ptr<Node> lambdaDman;

    double c;

    static double domainSize;

    int dofid;

    void computeDomainSize();

    bool elementEdges;

    IntArray elements;
    IntArray sides;
    IntArray lambdaIDs;

public:
    PrescribedMean(int n, Domain * d) : ActiveBoundaryCondition(n, d), lambdaDman(std::make_unique<Node>(0, this->domain) ) {}
    virtual ~PrescribedMean() {}

    void initializeFrom(InputRecord &ir) override;

    void assemble(SparseMtrx &answer, TimeStep *tStep, CharType type,
                  const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s, double scale = 1.0) override;

    void assembleVector(FloatArray &answer, TimeStep *tStep,
                        CharType type, ValueModeType mode,
                        const UnknownNumberingScheme &s, FloatArray *eNorm=nullptr) override;

    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep,
                                  CharType type, ValueModeType mode,
                                  const UnknownNumberingScheme &s, FloatArray *eNorm=nullptr);

    void giveExternalForcesVector(FloatArray &answer, TimeStep *tStep,
                                  CharType type, ValueModeType mode,
                                  const UnknownNumberingScheme &s);

    int giveNumberOfInternalDofManagers() override { return 1; }

    DofManager *giveInternalDofManager(int i) override { return lambdaDman.get(); }

    const char *giveClassName() const override { return "PrescribedMean"; }
    const char *giveInputRecordName() const override { return _IFT_PrescribedMean_Name; }
};

}

#endif // PRESCRIBEDMEAN_H
