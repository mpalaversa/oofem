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

#ifndef nodalload_h
#define nodalload_h

#include "load.h"
#include "bcgeomtype.h"
#include "bctype.h"
#include "valuemodetype.h"
#include "element.h"
#include "node.h"

#include <memory>

namespace oofem {
///@name Input fields for nodal loads
//@{
#define _IFT_NodalLoad_Name "nodalload"
#define _IFT_NodalLoad_cstype "cstype"
#define _IFT_NodalLoad_loadtype "loadtype"
#define _IFT_NodalLoad_sn "sn"
#define _IFT_NodalLoad_a3 "a3"
#define _IFT_NodalLoad_b4 "b4"
//@}
//class Element;
class TimeStep;

/**
 * Class implementing a concentrated load (force, moment,...) that acts
 * directly on a dof manager (node or element side, if it has associated DOFs).
 * This load could not be applied on an element.
 * A nodal load is usually attribute of one or more nodes or element sides.
 *
 * The component array, which size should be same as number of DOFs in particular
 * node/side is read from input.
 *
 * The attribute componentArray contains, for example for the case of a
 * plane beam structure, 2 forces and 1 moment on a right place. (6 dof per
 * node is assumed)
 *
 * @note{Load is not restricted to structural problems. For example, in thermal
 * analysis, a nodal load would be a concentrated heat source.}
 */
class OOFEM_EXPORT NodalLoad : public Load
{
protected:
    /**
     * Load coordinate system.
     * It is actually used only when local coordinate system in node is defined and load is specified in global
     * coordinate system
     */
    CoordSystType coordSystemType;
    /// Load type (its physical meaning).
    bcType lType;
    // Solidity ratio (must be specified in the input file when NodalLoad with loadType 9 is used)
    double sn;
    // Coefficients that can be specified in the input file when NodalLoad with loadType 9 is used.
    // Both are 0.05 by default.
    double a3, b4;

public:
    /**
     * Constructor. Creates nodal load object with given number, belonging to given domain.
     * @param n Load  number.
     * @param d Domain to which new object will belongs.
     */
    NodalLoad( int n, Domain *d );

    const char *giveInputRecordName() const override { return _IFT_NodalLoad_Name; }
    void computeValueAt(FloatArray &answer, TimeStep *tStep, const FloatArray &coords, ValueModeType mode) override
    { computeComponentArrayAt(answer, tStep, mode); }
    CoordSystType giveCoordSystMode() override { return coordSystemType; }

    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;
    const char *giveClassName() const override { return "NodalLoad"; }
    bcGeomType giveBCGeoType() const override { return NodalLoadBGT; }
    bcType giveType() const override { return lType; }
    double giveSolidityRatio() { return sn; }

     /**
     * A method for calculating hydrodynamic load based on "Kristiansen, T., Faltinsen, O.M., 2012. Modelling of
     * current loads on aquaculture net cages. Journal of Fluids and Structures 34 (2012), 218–235."
     * @param answer Force based on the implemented method.
     * @param associatedElements Elements associated with the node where the load is evaluated.
     * @param node Node at which the load is evaluated.
     * @param tStep Current time step.
     */
    void computeHydrodynamicKFLoad( FloatArray &answer, std::vector<Element*> associatedElements, DofManager *node, TimeStep *tStep );
    double computeDragCoefficientOfCircularCylinder( double density, double mu, double dt, double relativeVelocity, double sn );
    bool NodalLoad::isEdgePairAlreadyFound( std::vector<IntArray> panelEdges, int edge1, int edge2 );

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;
};
} // end namespace oofem
#endif // nodalload_h
