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

#ifndef decoupledcrosssection_h
#define decoupledcrosssection_h

#include "crosssection.h"

///@name Input fields for DecoupledCrossSection
//@{
#define _IFT_DecoupledCrossSection_Name "decoupledcrosssection"
#define _DecoupledCrossSection_CharacteristicDim "characteristicdim"
#define _DecoupledCrossSection_dragcoeff "dragcoeff"
#define _DecoupledCrossSection_Material "material"
//@}


namespace oofem {
class Material;
class Element;
class FloatArray;
class FloatMatrix;

/// @todo This should be an abstract class whose derived classes represent cross-sections of different geometries (see
/// the hierarchy of the DecoupledMaterial class). At this point, everything is implemented in this class supporting only
/// cross-sections that have one characteristic dimension such as circles and squares.
/**
 * Decoupled cross-sections are used to associate secondary cross-section properties with a FE (the primary being those
 * pertinent to the basic abstraction the FE represents). For example, a structural FE has the primary cross-section
 * associated with it that provides data used in a structural analysis such as area of the cross-section for a truss FE.
 * The secondary properties are associated with the FE and are used in, for example, a fluid flow
 * analysis based on a simplified model defined as loads in the structural analysis (e.g. see computeHydrodynamicLoadVector
 * method of StructuralElement class).
 */
class OOFEM_EXPORT DecoupledCrossSection : public CrossSection
{
protected:
    /// <summary>
    /// Represents a characteristic dimension of a decoupled cross-section
    /// (e.g. diameter for a circular cross-section or length of a side for
    /// a square cross-section).
    /// </summary>
    double characteristicDim, userDefinedDragCoeff;
    int materialNumber;

public:
    /**
     * Constructor. Creates cross section with given number, belonging to given domain.
     * @param n Cross section number.
     * @param d Domain to which new cross section will belong.
     */
    DecoupledCrossSection(int n, Domain *d) : CrossSection(n, d)  { }

    const char *giveClassName() const override { return "DecoupledCrossSection"; }
    
    const char *giveInputRecordName() const override { return _IFT_DecoupledCrossSection_Name; }

    double giveMagnitudeOfMaterialProperty( int property );
    
    double giveCharacteristicDimension() { return characteristicDim; }
    double giveDragCoefficient() override { return userDefinedDragCoeff; }
    
    Material *giveMaterial();

    int giveMaterialNumber() { return materialNumber; }
    
    void initializeFrom( InputRecord &ir ) override;

    bool isDecoupled() override;

    void restoreIPContext( DataStream &stream, ContextMode mode, GaussPoint *gp ) override { }
    void saveIPContext( DataStream &stream, ContextMode mode, GaussPoint *gp ) override { }
    int giveIPValue( FloatArray &answer, GaussPoint *ip, InternalStateType type, TimeStep *tStep ) override { return 0; }
    int setupIntegrationPoints( IntegrationRule &irule, int npointsXY, int npointsZ, Element *element ) override { return 0; }
    bool hasProperty( CrossSectionProperty a ) override { return false; }
    double predictRelativeComputationalCost( GaussPoint *ip ) override { return 0.0; }
    int packUnknowns( DataStream &buff, TimeStep *tStep, GaussPoint *ip ) override { return 0; }
    int unpackAndUpdateUnknowns( DataStream &buff, TimeStep *tStep, GaussPoint *ip ) override { return 0; }
    int estimatePackSize( DataStream &buff, GaussPoint *ip ) override { return 0; }
    Material *giveMaterial( IntegrationPoint *ip ) const override { return ip->giveElement()->giveMaterial(); }
};
} // end namespace oofem
#endif // decoupledcrosssection_h
