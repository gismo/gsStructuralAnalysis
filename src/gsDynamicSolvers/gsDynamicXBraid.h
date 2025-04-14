 /** @file gsDynamicXBraid.h

    @brief This class provides an interface to XBraid for time integration methods deriving from \ref gsDynamicBase

    For more information about XBraid, please check \ref gsXBraid or

    [XBraid on GitHub] https://github.com/XBraid/xbraid/
    [XBraid cite] XBraid: Parallel multigrid in time. http://llnl.gov/casc/xbraid.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#pragma once
#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicBase.h>
#include <gsIO/gsOptionList.h>

#ifdef gsXBraid_ENABLED
#include <gsXBraid/gsXBraid.h>
#endif

namespace gismo
{

/**
    @brief Performs the arc length method to solve a nonlinear system of equations.

    \tparam T coefficient type

    \ingroup gsStructuralAnalysis
*/
template <class T>
class gsDynamicXBraid : public gsXBraid< gsMatrix<T> >
{

public:
    typedef typename std::function<void(const index_t&,const T&,const gsVector<T>&,const gsVector<T>&,const gsVector<T>&)> callback_type;

    virtual ~gsDynamicXBraid() {};

#ifdef gsXBraid_ENABLED
    /// Constructor
    ///
    /// @param[in]  solver    A dynamic solver
    /// @param[in]  comm      The MPI communication channel
    /// @param[in]  tstart    The start time of the simulation
    /// @param[in]  tend      The end time of the simulation
    /// @param[in]  numDofs   The number of degrees of freedom
    /// @param[in]  numSteps  The number time steps to perform
    ///
    gsDynamicXBraid(
                    const gsDynamicBase<T> * solver,
                    const gsMpiComm& comm,
                    const T& tstart,
                    const T& tend,
                    const index_t numDofs,
                    index_t numSteps = 10)
    :
    gsXBraid< gsMatrix<T> >::gsXBraid(comm, tstart, tend, (int)numSteps),
    m_solver(solver),
    m_numDofs(numDofs)
    {
        defaultOptions();
        m_callback = [](const index_t&,const T&,const gsVector<T>&,const gsVector<T>&,const gsVector<T>&){return;}; // set empty callback
    }
#else
    /// Constructor
    ///
    /// @param[in]  solver    A dynamic solver
    /// @param[in]  comm      The MPI communication channel
    /// @param[in]  tstart    The start time of the simulation
    /// @param[in]  tend      The end time of the simulation
    /// @param[in]  numDofs   The number of degrees of freedom
    /// @param[in]  numSteps  The number time steps to perform
    ///
    gsDynamicXBraid(
                    const gsDynamicBase<T> * solver,
                    const gsMpiComm& comm,
                    const T& tstart,
                    const T& tend,
                    const index_t numDofs,
                    index_t numSteps = 10)
    {
        GISMO_ERROR("XBraid needs to be enabled to use the gsDynamicXBraid interface")
    }
#endif

// Member functions
public:
    /**
     * @brief      Dets the default options
     */
    void defaultOptions()
    {
        // see gsXBraid/filedata/pde/heat2d_square_ibvp1.xml
        m_options.addInt("CFactor", "Coarsening factor of the parallel-in-time multigrid solver", 2);
        m_options.addInt("access", "Access level (never [=0], =after finished [=1(default)], each iteration [=2]", 1);
        m_options.addInt("maxIter", "Maximum iteration numbers  of the parallel-in-time multigrid solver", 100);
        m_options.addInt("maxLevel", "Maximum numbers of parallel-in-time multigrid levels", 30);
        m_options.addInt("minCLevel", "Minimum level of the parallel-in-time multigrid solver", 2);
        m_options.addInt("norm", "Temporal norm of the parallel-in-time multigrid solver (1-norm [=1], 2-norm [=2(default)], inf-norm [=3])", 2);
        m_options.addInt("numFMG", "Number of full multigrid steps of the parallel-in-time multigrid solver", 1);
        m_options.addInt("numFMGVcyc", "Number of full multigrid V-cycles of the parallel-in-time multigrid solver", 1);
        m_options.addInt("numMaxRef", "Maximum number of refinements of the parallel-in-time multigrid solver", 1);
        m_options.addInt("numRelax", "Number of relaxation steps of the parallel-in-time multigrid solver", 1);
        m_options.addInt("numStorage", "Number of storage of the parallel-in-time multigrid solver", -1);
        m_options.addInt("print", "Print level (no output [=0], runtime inforation [=1], run statistics [=2(default)], debug [=3])", 2);
        m_options.addReal("absTol", "Absolute tolerance of the parallel-in-time multigrid solver", 1e-6);
        m_options.addReal("relTol", "Relative tolerance of the parallel-in-time multigrid solver", 1e-3);
        m_options.addSwitch("fmg", "Perform full multigrid (default is off)", 0);
        m_options.addSwitch("incrMaxLevels", "Increase the maximum number of parallel-in-time multigrid levels after performing a refinement (default is off)", 0);
        m_options.addSwitch("periodic", "Periodic time grid (default is off)", 0);
        m_options.addSwitch("refine", "Perform refinement in time (default off)", 0);
        m_options.addSwitch("sequential", "Set the initial guess of the parallel-in-time multigrid solver as the sequential time stepping solution (default is off)", 0);
        m_options.addSwitch("skip", "Skip all work on the first down cycle of the parallel-in-time multigrid solver (default on)", 1);
        m_options.addSwitch("spatial", "Perform spatial coarsening and refinement (default is off)", 0);
        m_options.addSwitch("tol", "Tolerance type (absolute [=true], relative [=false(default)]", 0);


        m_options.addSwitch("extraVerbose", "Extra verbosity", 0);
    }

    /**
     * @brief      Initializes the class
     */
    void initialize()
    {
        // see gsXBraid/examples/xbraid_heatEquation_example.cpp
        this->SetCFactor(m_options.getInt("CFactor"));
        this->SetMaxIter(m_options.getInt("maxIter"));
        this->SetMaxLevels(m_options.getInt("maxLevel"));
        this->SetMaxRefinements(m_options.getInt("numMaxRef"));
        this->SetMinCoarse(m_options.getInt("minCLevel"));
        this->SetNFMG(m_options.getInt("numFMG"));
        this->SetNFMGVcyc(m_options.getInt("numFMGVcyc"));
        this->SetNRelax(m_options.getInt("numRelax"));
        this->SetAccessLevel(m_options.getInt("access"));
        this->SetPrintLevel(m_options.getInt("print"));
        this->SetStorage(m_options.getInt("numStorage"));
        this->SetTemporalNorm(m_options.getInt("norm"));

        if (m_options.getSwitch("fmg"))           this->SetFMG();
        if (m_options.getSwitch("incrMaxLevels")) this->SetIncrMaxLevels();
        if (m_options.getSwitch("periodic"))      this->SetPeriodic(1); else this->SetPeriodic(0);
        if (m_options.getSwitch("refine"))        this->SetRefine(1);   else this->SetRefine(0);
        if (m_options.getSwitch("sequential"))    this->SetSeqSoln(1);  else this->SetSeqSoln(0);
        if (m_options.getSwitch("skip"))          this->SetSkip(1);     else this->SetSkip(0);
        if (m_options.getSwitch("spatial"))       this->SetSpatialCoarsenAndRefine();
        if (m_options.getSwitch("tol"))           this->SetAbsTol(m_options.getReal("absTol"));
        else                                      this->SetRelTol(m_options.getReal("relTol"));
    }

#ifdef gsXBraid_ENABLED
    /// See \ref gsXBraid for the documentation
    braid_Int Init(braid_Real    t, braid_Vector *u_ptr) override
    {
        gsMatrix<T>* u = new gsMatrix<T>(3*m_numDofs, 1);

        // Does this mean zero displacements?
        u->setZero();

        if (m_solver->solutionU().rows()==m_numDofs)
            u->col(0).segment(0          ,m_numDofs) = m_solver->solutionU();
        if (m_solver->solutionV().rows()==m_numDofs)
            u->col(0).segment(m_numDofs  ,m_numDofs) = m_solver->solutionV();
        if (m_solver->solutionA().rows()==m_numDofs)
            u->col(0).segment(2*m_numDofs,m_numDofs) = m_solver->solutionA();

        *u_ptr = (braid_Vector) u;
        return braid_Int(0);
    }
#endif

    /**
     * @brief      Access the options of the class
     *
     * @return     The options stored in the class
     */
    gsOptionList & options() { return m_options; }

    /**
     * @brief      Set the class options
     *
     * @param[in]  options  The options to be set
     */
    void setOptions(gsOptionList options) {m_options.update(options,gsOptionList::addIfUnknown); };

    /**
     * @brief      Provide a callback function to be executed after the simulation
     *
     * @param[in]  callback  The callback function
     */
    void setCallback(callback_type callback) const {m_callback = callback;}

#ifdef gsXBraid_ENABLED
    /// See \ref gsXBraid for the documentation
    braid_Int Step(braid_Vector    u, braid_Vector    ustop, braid_Vector    fstop, BraidStepStatus &status) override
    {
        gsVector<T>* u_ptr = (gsVector<T>*) u;
        // gsMatrix<T>* ustop_ptr = (gsMatrix<T>*) ustop; // the guess is not used

        // XBraid forcing
        if (fstop != NULL)
        {
            gsVector<T>* fstop_ptr = (gsVector<T>*) fstop;
            *u_ptr += *fstop_ptr;
        }

        // Get time step information
        std::pair<braid_Real, braid_Real> time = static_cast<gsXBraidStepStatus&>(status).timeInterval();
        if (m_options.getSwitch("extraVerbose")) gsInfo<<"Solving interval ["<<time.first<<" , "<<time.second<<"] (level "<<static_cast<gsXBraidStepStatus&>(status).level()<<")\n";
        T t  = time.first;
        T dt = time.second - time.first;

        // Solve time step
        gsVector<T> U = (*u_ptr).segment(0          ,m_numDofs);
        gsVector<T> V = (*u_ptr).segment(m_numDofs  ,m_numDofs);
        gsVector<T> A = (*u_ptr).segment(2*m_numDofs,m_numDofs);

        gsStatus stepStatus = m_solver->step(t,dt,U,V,A);

        u_ptr->segment(0          ,m_numDofs) = U;
        u_ptr->segment(m_numDofs  ,m_numDofs) = V;
        u_ptr->segment(2*m_numDofs,m_numDofs) = A;

        // Carry out adaptive refinement in time
        if (static_cast<gsXBraidStepStatus&>(status).level() == 0)
        {
            if (stepStatus==gsStatus::Success)
            {
                braid_Real error = static_cast<gsXBraidStepStatus&>(status).error();
                if (error != braid_Real(-1.0))
                {
                    braid_Int rfactor = (braid_Int) std::ceil( std::sqrt( error / 1e-3) );
                    status.SetRFactor(rfactor);
                }
                else
                    status.SetRFactor(1);
            }
            // Refine if solution interval failed
            else
            {
                if (m_options.getSwitch("extraVerbose")) gsInfo<<"Step "<<(static_cast<gsXBraidStepStatus&>(status)).timeIndex()<<" did not converge";
                status.SetRFactor((braid_Int)2);
            }
        }
        return braid_Int(0);
    }

    /// See \ref gsXBraid for the documentation
    braid_Int SpatialNorm(  braid_Vector  u,
                            braid_Real   *norm_ptr) override
    {
        gsVector<T>* u_ptr = (gsVector<T>*) u;
        *norm_ptr = u_ptr->segment(0,m_numDofs).norm(); // Displacement-based norm
        // *norm_ptr = u_ptr->norm();
        return braid_Int(0);
    }

    /// See \ref gsXBraid for the documentation
    braid_Int BufSize(braid_Int *size_ptr, BraidBufferStatus &status) override
    {
        *size_ptr = sizeof(T)*(m_numDofs*3+2); // +2 comes from rows, cols of the solution vector u.
        return braid_Int(0);
    }

    /// See \ref gsXBraid for the documentation
    braid_Int Access(braid_Vector u, BraidAccessStatus &status) override
    {
        gsVector<T>* u_ptr = (gsVector<T>*) u;
        m_callback((index_t)    static_cast<gsXBraidAccessStatus&>(status).timeIndex(),
                   (T)          static_cast<gsXBraidAccessStatus&>(status).time(),
                   (gsVector<T>)(*u_ptr).segment(0          ,m_numDofs),
                   (gsVector<T>)(*u_ptr).segment(m_numDofs  ,m_numDofs),
                   (gsVector<T>)(*u_ptr).segment(2*m_numDofs,m_numDofs)
                   );
        return braid_Int(0);
    }

    /// See \ref gsXBraid for the documentation
    /*
        NOTE: This routine is not implemented. How to do it:
        1. Make the Coarsen and Refine routines virtual
        2. Within the example, overload this class, and define Coarsen and Refine for the problem
        3. Update the m_numDoFs within!!
    */
    braid_Int Coarsen(braid_Vector fu, braid_Vector *cu_ptr, BraidCoarsenRefStatus &status) override
    {
        gsMatrix<T> *fu_ptr = (gsMatrix<T>*) fu;
        gsMatrix<T>* cu     = new gsMatrix<T>();
        *cu = *fu_ptr;
        *cu_ptr = (braid_Vector) cu;
        return braid_Int(0);
    }

    /// See \ref gsXBraid for the documentation
    braid_Int Refine(braid_Vector cu, braid_Vector *fu_ptr, BraidCoarsenRefStatus &status) override
    {
        gsMatrix<T> *cu_ptr = (gsMatrix<T>*) cu;
        gsMatrix<T>* fu     = new gsMatrix<T>();
        *fu = *cu_ptr;
        *fu_ptr = (braid_Vector) fu;
        return braid_Int(0);
    }
#endif

// Class members
protected:

    const gsDynamicBase<T> * m_solver;
    index_t m_numDofs;
    gsOptionList m_options;
    mutable callback_type m_callback;
};

} // namespace gismo