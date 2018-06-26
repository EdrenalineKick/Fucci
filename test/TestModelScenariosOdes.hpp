#ifndef PROJECTS_FUCCI_TEST_TESTMODELSCENARIOSODES_HPP_
#define PROJECTS_FUCCI_TEST_TESTMODELSCENARIOSODES_HPP_

#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"
#include "SmartPointers.hpp"

#include "TysonNovak2001OdeSystem.hpp"
#include "G2TysonNovak2001OdeSystem.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"
#include "FakePetscSetup.hpp"

class TestModelScenariosOde: public AbstractCellBasedTestSuite
{
    public:
    void TestSolvingOdes()
    {
        G2TysonNovak2001OdeSystem my_ode;

        RungeKutta4IvpOdeSolver euler_solver;
        std::vector<double> initial_conditions = my_ode.GetInitialConditions();
        OdeSolution solutions = euler_solver.Solve(&my_ode,initial_conditions,0,5,0.0001,0.1);
        //OdeSolution solutions = euler_solver.Solve(&my_ode,)

        for (unsigned i = 0; i < solutions.rGetTimes().size(); i++)
        {
            std::cout << solutions.rGetTimes()[i] << " "
                      << solutions.rGetSolutions()[i][0] << " "
                      << solutions.rGetSolutions()[i][1] << " "
                      << solutions.rGetSolutions()[i][2] << " "
                      << solutions.rGetSolutions()[i][3] << " "
                      << solutions.rGetSolutions()[i][4] << " "
                      << solutions.rGetSolutions()[i][5] << "\n";

        }
        solutions.WriteToFile("G2Ode","my_ode_solution", "sec");
    }
};

#endif /* PROJECTS_FUCCI_TEST_TESTMODELSCENARIOSODES_HPP_ */