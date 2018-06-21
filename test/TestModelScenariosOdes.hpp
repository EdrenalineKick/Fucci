#ifndef PROJECTS_FUCCI_TEST_TESTMODELSCENARIOSODES_HPP_
#define PROJECTS_FUCCI_TEST_TESTMODELSCENARIOSODES_HPP_

#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"
#include "SmartPointers.hpp"

#include "ModelScenariosOdeSystem.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"
#include "FakePetscSetup.hpp"

class TestModelScenariosOde: public AbstractCellBasedTestSuite
{
    public:
    void TestSolvingOdes()
    {
        ModelScenariosOdeSystem my_ode;

        BackwardEulerIvpOdeSolver euler_solver(5);
        std::vector<double> initial_conditions = my_ode.GetInitialConditions();
        OdeSolution solutions = euler_solver.Solve(&my_ode,initial_conditions,0,15,0.01,0.1);
        //OdeSolution solutions = euler_solver.Solve(&my_ode,)

        for (unsigned i = 0; i < solutions.rGetTimes().size(); i++)
        {
            std::cout << solutions.rGetTimes()[i] << " "
                      << solutions.rGetSolutions()[i][0] << " "
                      << solutions.rGetSolutions()[i][1] << " "
                      << solutions.rGetSolutions()[i][2] << " "
                      << solutions.rGetSolutions()[i][3] << " "
                      << solutions.rGetSolutions()[i][4] << "\n";

        }
        solutions.WriteToFile("SolvingModelScenarios","my_ode_solution", "sec");
    }
};

#endif /* PROJECTS_FUCCI_TEST_TESTMODELSCENARIOSODES_HPP_ */