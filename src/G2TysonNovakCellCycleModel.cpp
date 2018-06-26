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
#include "G2TysonNovakCellCycleModel.hpp"

G2TysonNovakCellCycleModel::G2TysonNovakCellCycleModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
        : AbstractImprovedOdeBasedPhaseBasedCellCycleModel(SimulationTime::Instance()->GetTime(), pOdeSolver)
{
    if (!mpOdeSolver)
    {
#ifdef CHASTE_CVODE
        mpOdeSolver = CellCycleModelOdeSolver<G2TysonNovakCellCycleModel, CvodeAdaptor>::Instance();
        mpOdeSolver->Initialise();
        // Chaste solvers always check for stopping events, CVODE needs to be instructed to do so
        mpOdeSolver->CheckForStoppingEvents();
        mpOdeSolver->SetMaxSteps(10000);
#else
        mpOdeSolver = CellCycleModelOdeSolver<G2TysonNovakCellCycleModel, RungeKutta4IvpOdeSolver>::Instance();
        mpOdeSolver->Initialise();
#endif //CHASTE_CVODE
    }
}

G2TysonNovakCellCycleModel::~G2TysonNovakCellCycleModel()
{
}

G2TysonNovakCellCycleModel::G2TysonNovakCellCycleModel(const G2TysonNovakCellCycleModel& rModel)
        : AbstractImprovedOdeBasedPhaseBasedCellCycleModel(rModel)
{
    /*
     * Initialize only those member variables defined in this class.
     * Create the new cell-cycle model's ODE system and use the current
     * values of the state variables in mpOdeSystem as an initial condition.
     *
     * The member variables mDivideTime and mG2PhaseStartTime are
     * initialized in the AbstractOdeBasedPhaseBasedCellCycleModel
     * constructor.
     *
     * The member variables mCurrentCellCyclePhase, mG1Duration,
     * mMinimumGapDuration, mStemCellG1Duration, mTransitCellG1Duration,
     * mSDuration, mG2Duration and mMDuration are initialized in the
     * AbstractPhaseBasedCellCycleModel constructor.
     *
     * The member variables mBirthTime, mReadyToDivide and mDimension
     * are initialized in the AbstractCellCycleModel constructor.
     */
    assert(rModel.GetOdeSystem());
    SetOdeSystem(new G2TysonNovak2001OdeSystem);
    SetStateVariables(rModel.GetOdeSystem()->rGetStateVariables());
}

AbstractCellCycleModel* G2TysonNovakCellCycleModel::CreateCellCycleModel()
{
    // Create a new cell-cycle model
    return new G2TysonNovakCellCycleModel(*this);
}

void G2TysonNovakCellCycleModel::Initialise()
{
    assert(mpOdeSystem == nullptr);
    assert(mpCell != nullptr);

    mpOdeSystem = new G2TysonNovak2001OdeSystem;
    mpOdeSystem->SetStateVariables(mpOdeSystem->GetInitialConditions());

    AbstractImprovedOdeBasedPhaseBasedCellCycleModel::Initialise();
}

void G2TysonNovakCellCycleModel::UpdateCellCyclePhase()
{
    assert(mpOdeSystem != nullptr);

    double current_time = SimulationTime::Instance()->GetTime();


    // Update the phase from M to G1 when necessary
    if (mCurrentCellCyclePhase == M_PHASE)
    {
        double m_duration = GetMDuration();
        if (GetAge() >= m_duration)
        {
            mCurrentCellCyclePhase = G_ONE_PHASE;
            mLastTime = m_duration + mBirthTime;
        }
        else
        {
            // Still dividing; don't run ODEs
            return;
        }
    }

    if (current_time > mLastTime)
    {
        if (!this->mFinishedRunningOdes)
        {
            // Update whether a stopping event has occurred
            this->mFinishedRunningOdes = SolveOdeToTime(current_time);
            //UpdateG1Time(*this);
            
            // Check no concentrations have gone negative
            for (unsigned i=0; i<mpOdeSystem->GetNumberOfStateVariables(); i++)
            {
                if (mpOdeSystem->rGetStateVariables()[i] < -DBL_EPSILON)
                {
                    // LCOV_EXCL_START
                    EXCEPTION("A protein concentration " << i << " has gone negative (" <<
                              mpOdeSystem->rGetStateVariables()[i] << ")\n"
                              << "Chaste predicts that the CellCycleModel numerical method is probably unstable.");
                    // LCOV_EXCL_STOP
                }
            }

            if (this->mFinishedRunningOdes)
            {
                // Update durations of each phase
                mG1Duration = GetG1Time() - mBirthTime - GetMDuration();
                SetSDuration(0.1);
                mG2PhaseStartTime = GetG1Time() + GetSDuration();
                mG2Duration = GetOdeStopTime() - mG2PhaseStartTime;
                mDivideTime = GetOdeStopTime() + 0.2;

                std::cout << "Stop Time" << GetOdeStopTime();
                std::cout << "G1 " << GetG1Time() << " G2 " <<mG2PhaseStartTime <<  " Divide " << mDivideTime << "\n";
                std::cout << "Total time " << GetMDuration() + GetG1Duration() + GetSDuration() + GetG2Duration() <<"\n";
                // Update phase
                if (current_time >= mG2PhaseStartTime)
                {
                    mCurrentCellCyclePhase = G_TWO_PHASE;
                }
                else
                {
                    mCurrentCellCyclePhase = S_PHASE;
                }
            }
        }
        else
        {
            // ODE model finished, just increasing time until division...
            if (current_time >= mG2PhaseStartTime)
            {
                std::cout << "ODE finished \n";
                ReadyToDivide();
            }
        }
    }
}

void G2TysonNovakCellCycleModel::UpdateG1Time(G2TysonNovakCellCycleModel& rModel)
{
    assert(mpOdeSystem != nullptr);
    assert(mpCell != nullptr);

    bool g1_true = static_cast<G2TysonNovak2001OdeSystem*>(mpOdeSystem)->CalculateG1Event(SimulationTime::Instance()->GetTime(), mpOdeSystem->rGetStateVariables());

    if (g1_true == true)
    {
        rModel.G1Trigger();
        rModel.GetG1Time();
    }
}

void G2TysonNovakCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    // No new parameters to output, so just call method on direct parent class
    AbstractImprovedOdeBasedPhaseBasedCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(G2TysonNovakCellCycleModel)
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(G2TysonNovakCellCycleModel)
