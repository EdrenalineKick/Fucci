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
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "CellProliferativePhasesWriter.hpp"
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
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        cell_population.AddPopulationWriter<VoronoiDataWriter>();

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("G2TestDoob");
        simulator.SetEndTime(30.0);

        simulator.SetSamplingTimestepMultiple(12);
        
        simulator.Solve();

        // //generate mesh
		// HoneycombMeshGenerator generator(2,2,2);
		// MutableMesh<2,2>* p_mesh = generator.GetMesh();

		// //Get location indices of the nodes on the mesh
		// std::vector<unsigned> location_indices = generator.GetCellLocationIndices();
		// std::vector<unsigned> real_indices;
		
		// for (unsigned i = 0 ; i < location_indices.size(); i++)
		// {
		// 	unsigned cell_index = location_indices[i];
		// 	if (i == 1)
		// 	{
		// 		real_indices.push_back(cell_index);
		// 	}
		// }

		// //Generate Cell Pointer vector and pointer to stem cell type and wild type mutation state
		// std::vector<CellPtr> cells;
		// MAKE_PTR(StemCellProliferativeType, p_stem_type);
		// boost::shared_ptr<AbstractCellProperty> p_state = CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>();


		// //Iterate over real_indices, assigning cell cycle and proliferative type to the cells
		// for (unsigned i = 0; i<real_indices.size(); i++)
		// {
        //     std::cout << "sweg";
		// 	//Set cell cycle
		// 	G2TysonNovakCellCycleModel* p_cycle_model = new G2TysonNovakCellCycleModel();
        //     std::cout << "bleh";
		// 	p_cycle_model->SetDimension(2);

		// 	//To avoid a 'pulsing' behaviour with birth events, we set each cell's initial age to be
		// 	// ~U(-12, 0) in the past, as each cell cycle duration is U(11, 13).
		// 	p_cycle_model->SetBirthTime(0);

		// 	CellPtr p_cell(new Cell(p_state, p_cycle_model));
		// 	p_cell->SetCellProliferativeType(p_stem_type); //Set the cell to be stem cell type
        //     std::cout << "bloh";
		// 	p_cell->InitialiseCellCycleModel();

		// 	cells.push_back(p_cell);
		// }

		// //Link mesh and cells together into a cell population
		// MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, real_indices);

		// cell_population.AddPopulationWriter<VoronoiDataWriter>();
		// cell_population.AddCellWriter<CellProliferativePhasesWriter>();
			
		// OffLatticeSimulation<2> simulator(cell_population);

		// MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
		// p_linear_force->SetCutOffLength(2.0);
		// p_linear_force->SetMeinekeSpringStiffness(45);
		// simulator.AddForce(p_linear_force);

		// simulator.SetOutputDirectory("Birthtime 0");
		// simulator.SetEndTime(60);
		// simulator.SetSamplingTimestepMultiple(60);

		// simulator.Solve();
    }
};

#endif /* PROJECTS_FUCCI_TEST_TESTMODELSCENARIOSODES_HPP_ */