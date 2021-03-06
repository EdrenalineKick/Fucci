/*

Copyright (c) 2005-2017, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "LocalNeighbourhoodDensity.hpp"
#include "Exception.hpp"

template <unsigned DIM>
LocalNeighbourhoodDensity<DIM>::LocalNeighbourhoodDensity()
        : AbstractCellBasedSimulationModifier<DIM>()
{
}

template <unsigned DIM>
LocalNeighbourhoodDensity<DIM>::~LocalNeighbourhoodDensity()
{
}

template <unsigned DIM>
void LocalNeighbourhoodDensity<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM, DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template <unsigned DIM>
void LocalNeighbourhoodDensity<DIM>::SetupSolve(AbstractCellPopulation<DIM, DIM>& rCellPopulation, std::string outputDirectory)
{
    UpdateCellData(rCellPopulation);
}

template <unsigned DIM>
void LocalNeighbourhoodDensity<DIM>::UpdateCellData(AbstractCellPopulation<DIM, DIM>& rCellPopulation)
{
    rCellPopulation.Update();

    for (typename AbstractCellPopulation<DIM, DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        CellPtr pCell = *cell_iter;

        c_vector<double, DIM> cellCoords = rCellPopulation.GetLocationOfCellCentre(pCell);

        unsigned neighbourCount = 0;

        for (typename AbstractCellPopulation<DIM, DIM>::Iterator nested_cell_iter = rCellPopulation.Begin();
             nested_cell_iter != rCellPopulation.End();
             ++nested_cell_iter)
        {
            CellPtr pComparedCell = *nested_cell_iter;
            c_vector<double, DIM> comparedCoords = rCellPopulation.GetLocationOfCellCentre(pComparedCell);
            c_vector<double, DIM> squaredDistance = zero_vector<double>(DIM);

            for (unsigned i = 0; i < DIM; i++)
            {
                squaredDistance[i] = pow((cellCoords[i] - comparedCoords[i]), 2.0);
            }

            double distance = pow(sum(squaredDistance), 0.5);

            if (distance <= 3.5)
            {
                neighbourCount++;
            }
        }

        pCell->GetCellData()->SetItem("Neighbourhood Count", neighbourCount);
    }
}

template <unsigned DIM>
void LocalNeighbourhoodDensity<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(LocalNeighbourhoodDensity)