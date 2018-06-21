#ifndef PROJECTS_FUCCI_TEST_TESTMODELG2PHASECELLCYCLE_HPP_
#define PROJECTS_FUCCI_TEST_TESTMODELG2PHASECELLCYCLE_HPP_

#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"
#include "CellsGenerator.hpp"
#include "G2TysonNovakCellCycleModel.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "OffLatticeSimulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "VoronoiDataWriter.hpp"
#include "FakePetscSetup.hpp"

class TestG2PhaseCellCycle : public AbstractCellBasedTestSuite
{
public:
    void TestG2Phase()
    {
        HoneycombMeshGenerator generator(2, 2);    // Parameters are: cells across, cells up
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        CellsGenerator<G2TysonNovakCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes());

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);


        cell_population.AddPopulationWriter<VoronoiDataWriter>();

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("G2Test");
        simulator.SetEndTime(30.0);

        simulator.SetSamplingTimestepMultiple(12);
        
        simulator.Solve();
    }
};

#endif /* PROJECTS_FUCCI_TEST_TESTMODELSCENARIOSODES_HPP_ */