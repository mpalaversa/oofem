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
 
#include "sm/Elements/Shells/qdshell.h"

#include "classfactory.h"
#include "sm/CrossSections/structuralcrosssection.h"
#include "sm/Materials/structuralms.h"

namespace oofem {
	QdShell::QdShell(int n, Domain* aDomain) : QdElement(n, aDomain)
	{
		outputCategory = OutputCategory::Combined;
		outputAtZ = 0.0;
	}

	void
	QdShell::giveDofManDofIDMask(int inode, IntArray& answer) const
	{
		answer = { D_u, D_v, D_w, R_u, R_v, R_w };
	}

	void
	QdShell::getStressesTopBottom(FloatArray& answer, TimeStep* tStep) {
		// Remove the following 3 lines of code when the method is considered generic.
		outputAtXY = OutputLocationXY::Centre;
		outputCategory = OutputCategory::Combined;
		outputType = OutputType::Standard;
		
		outputAtZ = this->giveStructuralCrossSection()->give(CS_Thickness, this->giveDefaultIntegrationRulePtr()->getIntegrationPoint(0)) / 2;

		FloatArray strains;
		computeStrainVector(strains, tStep);
		computeStressVector(answer, strains, tStep);
	}

	void
	QdShell::giveSurfaceDofMapping(IntArray& answer, int iSurf) const
	{
		if (iSurf == 1 || iSurf == 2) {
			answer.enumerate(24);
		}
		else {
			OOFEM_ERROR("wrong surface number");
		}
	}

    void
    QdShell::printOutputAt(FILE* file, TimeStep* tStep)
    {
        FloatArray v;
        fprintf(file, "element %d (%8d):\n", this->giveLabel(), number);
        int nPointsXY{ 0 }, nPointsZ{ 0 };
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
            OOFEM_ERROR("Something went wrong. The following options for output location at XY plane can be selected for ShellQd41 element: 'All', 'Centroid', 'Gauss points' and 'Corners'.");
            break;
        }

        FloatArray strains, stresses;
        switch (outputCategory)
        {
        case OutputCategory::Membrane:
            nPointsZ = 1;
            strains.resize(3);
            stresses.resize(3);
            break;
        case OutputCategory::Plate:
            nPointsZ = 1;
            strains.resize(3);
            stresses.resize(3);
            break;
        case OutputCategory::Combined:
            nPointsZ = 2;
            strains.resize(6);
            stresses.resize(6);
            break;
        case OutputCategory::All:
            nPointsZ = 2;
            break;
        default:
            OOFEM_ERROR("Something went wrong. The following options for output location at z can be selected for ShellQd41 element: 'All', 'Membrane', 'Plate' and 'Combined'.");
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

                if (nPointsZ == 1) {
                    if (outputCategory == OutputCategory::Plate) {
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
                    }
                    else {
                        fprintf(file, "\n          strains    ");
                        for (auto& val : ms->giveStrainVector())
                            fprintf(file, " %.4e", val);

                        fprintf(file, "\n          stresses    ");
                        for (auto& val : ms->giveStressVector())
                            fprintf(file, " %.4e", val);
                    }
                }
                else {
                    strains = ms->giveStrainVector();
                    stresses = ms->giveStressVector();
                    for (int j = 0; j < nPointsZ; j++) {
                        if (j == 0) {
                            fprintf(file, "\n           at z = %f :\n", outputAtZ);

                            fprintf(file, "                      strains    ");
                            for (int k = 1; k <= 6; k++) {
                                fprintf(file, " %.4e", strains.at(k));
                            }

                            fprintf(file, "\n                      stresses    ");
                            for (int k = 1; k <= 6; k++) {
                                fprintf(file, " %.4e", stresses.at(k));
                            }
                        }
                        else {
                            fprintf(file, "\n           at z = -%f :\n", outputAtZ);

                            fprintf(file, "                      strains    ");
                            for (int k = 7; k <= 12; k++) {
                                fprintf(file, " %.4e", strains.at(k));
                            }

                            fprintf(file, "\n                      stresses   ");
                            for (int k = 7; k <= 12; k++) {
                                fprintf(file, " %.4e", stresses.at(k));
                            }
                        }
                    }
                }

                fprintf(file, "\n");
                break;
            case OutputLocationXY::Centre:
                fprintf(file, "  Centroid");
                gp = integrationRulesArray[0]->getIntegrationPoint(0);
                ms = static_cast<StructuralMaterialStatus*>(gp->giveMaterialStatus());

                if (nPointsZ == 1) {
                    if (outputCategory == OutputCategory::Plate) {
                        fprintf(file, "\n          at z = %f :\n", outputAtZ);


                        fprintf(file, "                      strains    ");
                        for (auto& val : ms->giveStrainVector()) {
                            fprintf(file, " %.4e", val);
                        }

                        fprintf(file, "\n                      stresses   ");
                        for (auto& val : ms->giveStressVector()) {
                            fprintf(file, " %.4e", val);
                        }
                    }
                    else {
                        fprintf(file, "\n          strains    ");
                        for (auto& val : ms->giveStrainVector())
                            fprintf(file, " %.4e", val);

                        fprintf(file, "\n          stresses    ");
                        for (auto& val : ms->giveStressVector())
                            fprintf(file, " %.4e", val);
                    }
                }
                else {
                    strains = ms->giveStrainVector();
                    stresses = ms->giveStressVector();
                    for (int j = 0; j < nPointsZ; j++) {
                        if (j == 0) {
                            fprintf(file, "\n          at z = %f :\n", outputAtZ);

                            fprintf(file, "                      strains    ");
                            for (int k = 1; k <= 6; k++) {
                                fprintf(file, " %.4e", strains.at(k));
                            }

                            fprintf(file, "\n                      stresses    ");
                            for (int k = 1; k <= 6; k++) {
                                fprintf(file, " %.4e", stresses.at(k));
                            }
                        }
                        else {
                            fprintf(file, "\n          at z = -%f :\n", outputAtZ);

                            fprintf(file, "                      strains    ");
                            for (int k = 7; k <= 12; k++) {
                                fprintf(file, " %.4e", strains.at(k));
                            }

                            fprintf(file, "\n                      stresses   ");
                            for (int k = 7; k <= 12; k++) {
                                fprintf(file, " %.4e", stresses.at(k));
                            }
                        }
                    }
                }

                fprintf(file, "\n");
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
                OOFEM_ERROR("Something went wrong. The following options for output location at XY plane can be selected for ShellQd41 element: 'All', 'Centroid', 'Gauss points' and 'Corners'.");
                break;
            }
        }
    }
}