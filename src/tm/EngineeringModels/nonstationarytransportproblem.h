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

#ifndef nonstationarytransportproblem_h
#define nonstationarytransportproblem_h

#include "sparselinsystemnm.h"
#include "sparsemtrx.h"
#include "primaryfield.h"
#include "dofdistributedprimaryfield.h"
#include "tm/Materials/transportmaterial.h"
#include "tm/EngineeringModels/stationarytransportproblem.h"
#include "linsystsolvertype.h"

///@name Input fields for NonStationaryTransportProblem
//@{
#define _IFT_NonStationaryTransportProblem_Name "nonstationaryproblem"
#define _IFT_NonStationaryTransportProblem_initt "initt"
#define _IFT_NonStationaryTransportProblem_deltat "deltat"
#define _IFT_NonStationaryTransportProblem_deltatfunction "deltatfunction"
#define _IFT_NonStationaryTransportProblem_prescribedtimes "prescribedtimes"
#define _IFT_NonStationaryTransportProblem_alpha "alpha"
#define _IFT_NonStationaryTransportProblem_lumpedcapa "lumpedcapa"
#define _IFT_NonStationaryTransportProblem_changingproblemsize "changingproblemsize"
//@}

namespace oofem {


/**
 * Callback class for assembling element external forces:
 * - edge or surface load on elements
 * - add internal source vector on elements
 * @author Mikael Öhman
 */
class TransportExternalForceAssembler : public VectorAssembler
{
public:
    void vectorFromElement(FloatArray &vec, Element &element, TimeStep *tStep, ValueModeType mode) const override;
};


/**
 * Callback class for assembling mid point effective tangents
 * @author Mikael Öhman
 */
class MidpointLhsAssembler : public MatrixAssembler
{
protected:
    double lumped;
    double alpha;

public:
    MidpointLhsAssembler(bool lumped, double alpha);
    void matrixFromElement(FloatMatrix &mat, Element &element, TimeStep *tStep) const override;
};


/**
 * Callback class for assembling CBS pressure matrices
 * @author Mikael Öhman
 */
class IntSourceLHSAssembler : public MatrixAssembler
{
public:
    void matrixFromElement(FloatMatrix &mat, Element &element, TimeStep *tStep) const override;
};


/**
 * This class represents linear nonstationary transport problem.
 */
class NonStationaryTransportProblem : public StationaryTransportProblem
{
protected:
    /**
     * Contains last time stamp of internal variable update.
     * This update is made via various services
     * (like those for computing real internal forces or updating the internal state).
     */
    StateCounterType internalVarUpdateStamp;

    LinSystSolverType solverType = ST_Direct; ///@todo Remove this and use nonlinear methods.
    std :: unique_ptr< SparseLinearSystemNM > linSolver; ///@todo Remove this and use nonlinear methods.

    /// Right hand side vector from boundary conditions.
    FloatArray bcRhs;

    /// Initial time from which the computation runs. Default is zero.
    double initT = 0.;
    /// Length of time step.
    double deltaT = 0.;
    double alpha = 0.;

    /// If set then stabilization using lumped capacity will be used.
    int lumpedCapacityStab = 0;

    /// Associated time function for time step increment.
    int dtFunction = 0;

    /// Specified times where the problem is solved
    FloatArray discreteTimes;

    /// Determines if there are change in the problem size (no application/removal of Dirichlet boundary conditions).
    bool changingProblemSize = false;

public:
    NonStationaryTransportProblem(int i, EngngModel * _master);

    void solveYourselfAt(TimeStep *tStep) override;
    void updateYourself(TimeStep *tStep) override;
    double giveUnknownComponent(ValueModeType, TimeStep *tStep, Domain *d, Dof *dof) override;
    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    void updateDomainLinks() override;

    TimeStep *giveNextStep() override;
    TimeStep *giveSolutionStepWhenIcApply(bool force = false) override;
    NumericalMethod *giveNumericalMethod(MetaStep *mStep) override;

    void initializeFrom(InputRecord &ir) override;
    int checkConsistency() override;

    // identification
    const char *giveInputRecordName() const { return _IFT_NonStationaryTransportProblem_Name; }
    const char *giveClassName() const override { return "NonStationaryTransportProblem"; }
    fMode giveFormulation() override { return TL; }

    /// Allows to change number of equations during solution.
    int requiresUnknownsDictionaryUpdate() override { return changingProblemSize; }
    bool requiresEquationRenumbering(TimeStep *tStep) override { return changingProblemSize; }
    //Store solution vector to involved DoFs
    //void updateDofUnknownsDictionary(DofManager *dman, TimeStep *tStep) override;

    int giveUnknownDictHashIndx(ValueModeType mode, TimeStep *tStep) override;

    /**
     * Returns time function for time step increment.
     * Used time function should provide step lengths as function of step number.
     * Initial step with number 0 is considered as [ -dt(0), 0 ], first step is [ 0, dt(1) ], ...
     */
    Function *giveDtFunction();

    /**
     * Returns the time step length for given step number n, initial step is number 0.
     */
    double giveDeltaT(int n);

    /**
     * Returns time for time step number n (array discreteTimes must be specified)
     */
    double giveDiscreteTime(int n);

#ifdef __CEMHYD_MODULE
    void averageOverElements(TimeStep *tStep);
#endif

protected:
    virtual void assembleAlgorithmicPartOfRhs(FloatArray &rhs,
                                              const UnknownNumberingScheme &s, TimeStep *tStep);

    /**
     * This function is normally called at the first time to project initial conditions to previous (0^th) solution vector.
     * @param tStep Previous solution step.
     */
    virtual void applyIC(TimeStep *tStep);

    /**
     * Assembles part of RHS due to Dirichlet boundary conditions.
     * @param answer Global vector where the contribution will be added.
     * @param tStep Solution step.
     * @param mode Mode of result.
     * @param s A map of non-default equation numbering if required.
     * @param d Domain.
     */
    virtual void assembleDirichletBcRhsVector(FloatArray &answer, TimeStep *tStep, ValueModeType mode,
                                              const UnknownNumberingScheme &s, Domain *d);
    /**
     * Copy unknowns in DOF's from previous to current position.
     * @param mode What the unknown describes (increment, total value etc.).
     * @param fromTime From which time step to obtain value.
     * @param toTime To which time to copy.
     */
    virtual void copyUnknownsInDictionary(ValueModeType mode, TimeStep *fromTime, TimeStep *toTime);

    /**
     * Updates IP values on elements.
     * @param tStep Solution step.
     */
    virtual void updateInternalState(TimeStep *tStep);
};
} // end namespace oofem
#endif // nonstationarytransportproblem_h
