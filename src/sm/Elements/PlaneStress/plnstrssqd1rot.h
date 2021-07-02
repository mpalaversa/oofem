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
 
#ifndef plnstrssqd1rot_h
#define plnstrssqd1rot_h

#include "qdmembrane.h"
#include "femcmpnn.h"
#include "error.h"
#include "floatmatrix.h"
#include "floatarray.h"

#include <cstdio>
#include <vector>
#include <memory>

#define _IFT_PlnStrssQd1Rot_Name "plnstrssqd1rot"

namespace oofem {
	class PlnStrssQd1Rot : public QdMembrane
	{
	public:
		PlnStrssQd1Rot(int n, Domain* d);
		virtual ~PlnStrssQd1Rot() {}

		const char* giveClassName() const override { return "PlnStrssQd1Rot"; }
		const char* giveInputRecordName() const override { return _IFT_PlnStrssQd1Rot_Name; }
	};
}
#endif