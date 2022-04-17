/** @file gsThinShell_ArcLength.cpp

    @brief Code for the arc-length method of a shell based on loads

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#include <gismo.h>

#include <gsKLShell/gsThinShellAssembler.h>
#include <gsKLShell/getMaterialMatrix.h>
// #include <gsThinShell/gsArcLengthIterator.h>
#include <gsStructuralAnalysis/gsArcLengthIterator.h>

#include <gsHSplines/gsKdNode.h>


using namespace gismo;

template <class T>
void initStepOutput( const std::string name, const gsMatrix<T> & points);

template <class T>
void writeStepOutput(const T lambda, const gsMultiPatch<T> & deformation, const std::string name, const gsMatrix<T> & points, const index_t extreme=-1, const index_t kmax=100);

template<index_t dim, class T>
gsTensorBSpline<dim,T> gsSpaceTimeFit(const gsMatrix<T> & solutionCoefs, const gsVector<T> & times, const gsVector<T> & ptimes, gsMultiBasis<T> & spatialBasis, index_t deg = 2);


template<class T, class solution_t >
class gsSpaceTimeHierarchy
{

  typedef kdnode<1, T> node_t;

  typedef typename node_t::point point_t;

public:

  ~gsSpaceTimeHierarchy()
  {
    for (size_t l=0; l!=m_tree.size(); ++l)
      for (typename std::map<T,node_t * >::const_iterator it = m_tree[l].begin(); it != m_tree[l].end(); ++it)
      {
        gsDebug<<"delete m_tree["<<l<<"]["<<it->first<<"])\n";
        delete it->second;
      }
  }

  /**
   * @brief      { function_description }
   *
   * @param[in]  times      The times (MONOTONICALLY INCREASING!)
   * @param[in]  solutions  The solutions
   */
  gsSpaceTimeHierarchy(std::vector<T> times, std::vector<solution_t> solutions)
  :
  m_points(times.size())
  {
    GISMO_ASSERT(times.size()==solutions.size(),"Sizes must agree");

    m_solutions.push_back(std::map<T,solution_t>());
    m_solutions.push_back(std::map<T,solution_t>());
    for (size_t k=0; k!=times.size(); ++k)
    {
      m_solutions[0].insert({times.at(k),solutions.at(k)});

    }

    _defaultOptions();
    _init();
  }

  /**
   * @brief      { function_description }
   *
   * Initializes the class based on times without solutions. This constructor prepares a queue only; it does not fill the tree with solutions
   *
   * @param[in]  times      The times (MONOTONICALLY INCREASING!)
   */
  gsSpaceTimeHierarchy(std::vector<T> times)
  :
  m_points(times.size())
  {
    m_tree.resize(m_maxLevel+1);
    gsVector<T,1> low,upp;
    index_t k=0;
    for (typename std::vector<T>::const_iterator it = times.begin(); it != std::prev(times.end()); it++, k++)
    {
      /// Add parents to tree
      low<<*it;
      upp<<*std::next(it);
      m_tree[0][*it] = new node_t(low,upp);

      // Push items on level 0,
      m_queue.push({0,*it});
    }

    _defaultOptions();
  }

private:
  void _init()
  {
    m_tree.resize(m_maxLevel+1);
    gsVector<T,1> low,upp;
    index_t k=0;
    for (typename std::map<T,solution_t>::const_iterator it = m_solutions[0].begin(); it != std::prev(m_solutions[0].end()); it++, k++)
    {
      /// Add parents to tree
      low<<it->first;
      upp<<std::next(it)->first;
      m_tree[0][it->first] = new node_t(low,upp);
    }

    for (typename std::map<T,node_t * >::const_iterator it = m_tree[0].begin(); it != m_tree[0].end(); ++it)
    {
      low = it->second->lowCorner();
      upp = it->second->uppCorner();
      it->second->split(0,(upp[0]+low[0])/2);
      m_tree[1][it->second->left ->lowCorner()[0]] = it->second->left;
      m_tree[1][it->second->right->lowCorner()[0]] = it->second->right;

      /// Add childs to queue
      m_queue.push({1,it->second->left ->lowCorner()[0]});
      m_queue.push({1,it->second->right->lowCorner()[0]});
    }

    /// Add childs to queue
    // for (typename std::map<T,node_t * >::const_iterator it = m_tree[1].begin(); it != m_tree[1].end(); ++it)
    //   m_queue.push({1,it->first});
  }

  std::pair<point_t,point_t> _interval(node_t * node)
  {
    std::pair<point_t,point_t> interval;

    if (node->isLeaf())
    {
      interval.first  = node->lowCorner();
      interval.second = node->uppCorner();
      return interval;
    }
    else
    {
      node_t * left  = node->left;
      node_t * right = node->right;
      while (!left->isLeaf() || !right->isLeaf())
      {
        if (!left->isLeaf())
          left  = left ->left;
        if (!right->isLeaf())
          right = right->right;
      }
      interval.first  = left->lowCorner();
      interval.second = right->uppCorner();
      return interval;
    }
  }

  // void _initLevel(index_t level)
  // {
  //   GISMO_ASSERT(m_ptimes.size()==m_indices.size(),"times and indices have different level");
  //   if (m_ptimes.size() < level+1 && level < m_maxLevel)
  //   {
  //     gsVector<T> time;
  //     time.setLinSpaced(m_points*math::pow(2,level),0,1);
  //     m_ptime = std::vector(time.data(),time.data() + time.rows()*time.cols());
  //     gsVector<index_t> index;
  //     index.setLinSpaced(m_points*math::pow(2,level),0,m_points*math::pow(2,level)-1);
  //     m_index = std::vector(index.data(),index.data() + index.rows()*index.cols());

  //     m_ptimes.push_back(std::map<index_t,T>());
  //     m_indices.push_back(std::map<T,index_t>());
  //     GISMO_ASSERT(m_ptime.size()==m_index.size(),"Sizes are not the same!");
  //     for (index_t k=0; k!=m_ptime.size(); k++)
  //     {
  //       m_ptimes[level].insert({m_index[k],m_ptime[k]});
  //       m_indices[level].insert({m_ptime[k],m_index[k]});
  //     }

  //     // m_solutions.push_back(std::map<index_t,solution_t>());
  //     // m_times.push_back(std::map<index_t,T>());
  //   }
  // }
  //


  void _defaultOptions()
  {
    m_maxLevel = 5;
    m_options.addInt("MaxLevel","Sets the maximum level for hierarchical refinement",m_maxLevel);
  }

  void _applyOptions()
  {
    m_maxLevel = m_options.getInt("MaxLevel");
    m_tree.resize(m_maxLevel+1);

  }

public:

  gsOptionList & options() { return m_options; }

  std::tuple<T,    T,  solution_t, solution_t> getJob()
  //         time, dt, start,      guess
  {
    this->_applyOptions();

    GISMO_ASSERT(!m_queue.empty(),"The queue is empty! Something went wrong.");

    T time = m_queue.front().second;
    index_t level = m_queue.front().first;

    GISMO_ASSERT(m_tree[level].count(time)>0,"Node at level "<<level<<" with time "<<time<<"not found");

    GISMO_ASSERT(!m_tree[level][time]->isRoot(),"Node cannot be a parent!");

    gsVector<T> low, upp;
    std::tie(low,upp) = _interval(m_tree[level][time]);
    T tstart = low[0];
    T tend   = upp[0];
    T dt     = tend-tstart;

    gsDebugVar(level);

    solution_t start, guess;
    if (m_tree[level][time]->isLeftChild())
    {
      std::tie(low,upp) = _interval(m_tree[level][time]->parent);
      tend = upp(0);

      // Find best available start point
      index_t startLevel = level;
      while (startLevel > 0 && m_solutions[startLevel-1].count(tstart)==0)
        startLevel -=1;

      GISMO_ASSERT(m_solutions[startLevel-1].count(tstart)!=0,"Cannot find start point at tstart = "<<tstart<<" in level "<<startLevel-1);
      start = m_solutions[startLevel-1][tstart];

      // Find best available start point
      index_t guessLevel = level;
      while (guessLevel > 0 && m_solutions[guessLevel-1].count(tend)==0)
        guessLevel -=1;

      GISMO_ASSERT(m_solutions[guessLevel-1].count(tend)!=0,"Cannot find end point at tend = "<<tend<<" in level "<<guessLevel-1);
      guess = m_solutions[guessLevel-1][tend];
    }
    else if (m_tree[level][time]->isRightChild())
    {
      std::tie(low,upp) = _interval(m_tree[level][time]->parent);
      tend = upp(0);
      GISMO_ASSERT(m_solutions[level].count(tstart)!=0,"Cannot find start point at tstart = "<<tstart<<" in level "<<level);
      start = m_solutions[level][tstart];
      GISMO_ASSERT(m_solutions[level-1].count(tend)!=0,"Cannot find end point at tend = "<<tend<<" in level "<<level-1);
      guess = m_solutions[level-1][tend];
    }
    else
      GISMO_ERROR("Node is not a child!");

    return std::make_tuple(tstart, dt, start, guess);
  }

  bool reference_into(solution_t & result)
  {
    T time = m_queue.front().second;
    index_t level = m_queue.front().first;

    gsVector<T> low, upp;
    std::tie(low,upp) = _interval(m_tree[level][time]);

    if (m_solutions[level-1].count(upp[0])!=0)
    {
      result =  m_solutions[level-1][upp[0]];
      return true;
    }
    else
      return false;
  }

  void submitJob(solution_t solution)
  {
    T time = m_queue.front().second;
    index_t level = m_queue.front().first;
    node_t * tmp = m_tree[level][time];

    gsVector<T> low, upp;
    std::tie(low,upp) = _interval(tmp);
    T tend   = upp[0];

    if (m_solutions.size() < level + 1)
      m_solutions.resize(level+1);

    m_solutions[level][tend] = solution;
    gsDebug<<"Job submitted at t = "<<tend<<"; level "<<level<<"\n";
  }

  void addJobHere()
  {
    T time = m_queue.front().second;
    index_t level = m_queue.front().first;

    if (level+1 > m_maxLevel)
    {
      gsWarn<<"Max level reached!";
      return;
    }

    gsVector<T> low, upp;
    std::tie(low,upp) = _interval(m_tree[level][time]);

    m_tree[level][time]->split(0,(upp[0]+low[0])/2);

    m_tree[level+1][m_tree[level][time]->left ->lowCorner()[0]] = m_tree[level][time]->left;
    m_tree[level+1][m_tree[level][time]->right->lowCorner()[0]] = m_tree[level][time]->right;

    gsDebugVar(m_tree[level][time]->left ->lowCorner()[0]);
    gsDebugVar(m_tree[level][time]->right->lowCorner()[0]);

    m_queue.push({level+1,m_tree[level][time]->left ->lowCorner()[0]});
    m_queue.push({level+1,m_tree[level][time]->right->lowCorner()[0]});
  }

  index_t currentLevel()
  {
    return m_queue.front().first;
  }

  void removeJob()
  {
    m_queue.pop();
  }

  bool empty()
  {
    return m_queue.empty();
  }

  std::pair<std::vector<T>,std::vector<solution_t>> getFlatSolution()
  {
    std::vector<T> times;
    std::vector<solution_t> solutions;

    std::map<T,solution_t> map;
    for (size_t l = 0; l!=m_solutions.size(); ++l)
      for (typename std::map<T,solution_t>::const_iterator it = m_solutions[l].begin(); it!=m_solutions[l].end(); ++it)
        map[it->first] = it->second; // overwrites if duplicate

    for (typename std::map<T,solution_t>::const_iterator it = map.begin(); it!=map.end(); ++it)
    {
      times.push_back(it->first);
      solutions.push_back(it->second);
    }

    return std::make_pair(times,solutions);
  }

  void print()
  {
    for (size_t l = 0; l!=m_solutions.size(); ++l)
    {
      gsInfo<<"level = "<<l<<":\n";
      for (typename std::map<T,solution_t>::const_iterator it = m_solutions[l].begin(); it!=m_solutions[l].end(); ++it)
      {
        gsInfo<<"\t";
        gsInfo<<"time = "<<it->first<<"\tsol = "<<it->second<<"\n";
      }
    }
  }

  void printTree()
  {
    gsVector<T,1> low, upp;
    for (size_t l=0; l!=m_tree.size(); ++l)
      for (typename std::map<T,node_t * >::const_iterator it = m_tree[l].begin(); it != m_tree[l].end(); ++it)
      {
        std::tie(low,upp) = _interval(it->second);
        if (l!=0)
          gsDebug<<"interval = ["<<low[0]<<","<<upp[0]<<"], level = "<<l<<"\n";
        else
          gsDebug<<"interval = ["<<low[0]<<","<<upp[0]<<"], level = "<<l<<"\n";
      }
  }

  void printQueue()
  {
    gsVector<T,1> low, upp;
    std::queue<std::pair<index_t,T>> queue = m_queue;
    gsInfo<<"Queue has "<<queue.size()<<" elements\n";
    while (!queue.empty())
    {
      gsInfo<<"Level "<<queue.front().first<<"; time "<<queue.front().second<<"\n";
      queue.pop();
    }
    gsInfo<<"Queue is empty\n";

  }


protected:
  index_t m_points;
  index_t m_maxLevel;

  std::vector<std::map<T,solution_t>>       m_solutions;

  std::vector<std::map<T,node_t *>>         m_tree;

  std::queue<std::pair<index_t,T>>   m_queue;

  gsOptionList m_options;

};

int main (int argc, char** argv)
{
    // Input options
    int numElevate    = 1;
    int numHref       = 1;
    bool plot         = false;
    bool mesh         = false;
    bool stress       = false;
    bool membrane     = false;
    bool quasiNewton  = false;
    int quasiNewtonInt= -1;
    bool adaptive     = false;
    int step          = 10;
    int method        = 2; // (0: Load control; 1: Riks' method; 2: Crisfield's method; 3: consistent crisfield method; 4: extended iterations)
    bool deformed     = false;

    bool composite = false;

    real_t relax      = 1.0;

    int testCase      = 0;

    int result        = 0;

    bool write        = false;

    index_t maxit     = 20;
    // index_t iniLevels  = 2;
    // index_t maxLevels  = 4;
    index_t maxLevel  = 2;

    // Arc length method options
    real_t dL        = 0.5; // Ard length to find bifurcation
    real_t tol        = 1e-6;
    real_t tolU       = 1e-6;
    real_t tolF       = 1e-3;

    index_t deg_z = 2;

    // MPI Stuff
    double     number;
    double     max_num;
    double     recv_num;
    int        proc_count;
    int        my_rank;
    int        iiter;
    // !MPI Stuff

    std::string wn("data.csv");

    std::string assemberOptionsFile("options/solver_options.xml");

    gsCmdLine cmd("Arc-length analysis for thin shells.");
    cmd.addString( "f", "file", "Input XML file for assembler options", assemberOptionsFile );

    cmd.addInt("t", "testcase", "Test case: 0: clamped-clamped, 1: pinned-pinned, 2: clamped-free", testCase);

    cmd.addInt("r","hRefine", "Number of dyadic h-refinement (bisection) steps to perform before solving", numHref);
    cmd.addInt("e","degreeElevation", "Number of degree elevation steps to perform on the Geometry's basis before solving", numElevate);
    cmd.addSwitch("composite", "Composite material", composite);

    cmd.addInt("m","Method", "Arc length method; 1: Crisfield's method; 2: RIks' method.", method);
    cmd.addReal("L","dL", "arc length", dL);
    // cmd.addInt("I","inilvl", "Initial levels", iniLevels);
    // cmd.addInt("M","maxlvl", "Max levels", maxLevels);
    cmd.addInt("l","level", "Max level", maxLevel);
    cmd.addReal("A","relaxation", "Relaxation factor for arc length method", relax);

    cmd.addInt("q","QuasiNewtonInt","Use the Quasi Newton method every INT iterations",quasiNewtonInt);
    cmd.addInt("N", "maxsteps", "Maximum number of steps", step);
    cmd.addInt("z", "degz", "Degree of fitting splin", deg_z);

    cmd.addSwitch("adaptive", "Adaptive length ", adaptive);
    cmd.addSwitch("quasi", "Use the Quasi Newton method", quasiNewton);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addSwitch("mesh", "Plot mesh?", mesh);
    cmd.addSwitch("stress", "Plot stress in ParaView format", stress);
    cmd.addSwitch("membrane", "Use membrane model (no bending)", membrane);
    cmd.addSwitch("deformed", "plot on deformed shape", deformed);
    cmd.addSwitch("write", "write to file", write);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // Conditional compilation
    #ifdef GISMO_WITH_MPI
      gsInfo << "Gismo was compiled with MPI support.\n";
    #else
      gsInfo << "Gismo was compiled without MPI support.\n";
    #endif

    // GISMO_ASSERT(maxLevels>iniLevels,"Max levels must  be more than initial levels!");


    // Initialize the MPI environment
    const gsMpi & mpi = gsMpi::init(argc, argv);

    // Get current wall time
    double wtime = mpi.wallTime();

    // Get the world communicator
    gsMpiComm comm = mpi.worldComm();
    MPI_Request req;

    //Get size and rank of the processor
    proc_count = comm.size();
    my_rank = comm.rank();

    // GISMO_ASSERT(proc_count > 1,"At least two processes are required.\n");

/* /////////////////////////////////////////////////////////////////////////////////////////
    Problem initialization

    Constructs all parameters required for the simulation
    Defines all solver options

    ALL PROCESSES
///////////////////////////////////////////////////////////////////////////////////////// */

    gsFileData<> fd(assemberOptionsFile);
    gsOptionList opts;
    fd.getFirst<gsOptionList>(opts);

    gsMultiPatch<> mp;
    real_t aDim;
    real_t bDim;


    real_t thickness;
    real_t Exx, Eyy, Gxy;
    real_t PoissonRatio = 0.3;
    real_t Density    = 1e0;

    if (composite)
    {
      Exx = 3300;
      Eyy = 1100;
      Gxy = 660;
      PoissonRatio = 0.25;
    }
    else
    {
      Exx  = 3102.75;
      PoissonRatio = 0.3;
    }


    if (testCase==1)
    {
      thickness = 6.35;
    }
    else if (testCase==2)
    {
      thickness = 12.7;
    }
    else if (testCase==3)
    {
      thickness = 16.75;
    }

    gsReadFile<>("surface/scordelis_lo_roof_shallow.xml", mp);

    for(index_t i = 0; i< numElevate; ++i)
      mp.patch(0).degreeElevate();    // Elevate the degree

    // h-refine
    for(index_t i = 0; i< numHref; ++i)
      mp.patch(0).uniformRefine();

    gsMultiBasis<> dbasis(mp);
    gsInfo<<"Basis (patch 0): "<< mp.patch(0).basis() << "\n";

    // Boundary conditions
    gsBoundaryConditions<> BCs;
    gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();

    // Initiate Surface forces
    std::string tx("0");
    std::string ty("0");
    std::string tz("0");

    gsVector<> tmp(mp.targetDim());
    gsVector<> neu(mp.targetDim());
    tmp.setZero();
    neu.setZero();
    gsConstantFunction<> neuData(neu,mp.targetDim());

    // Unscaled load
    real_t Load = 0;

    gsMatrix<> writePoints(2,3);
    writePoints.col(0)<< 0.0,0.5;
    writePoints.col(1)<< 0.5,0.5;
    writePoints.col(2)<< 1.0,0.5;

    GISMO_ASSERT(mp.targetDim()==3,"Geometry must be surface (targetDim=3)!");
    // Diaphragm conditions
    BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
    BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
    BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
    // BCs.addCornerValue(boundary::southwest, 0.0, 0, 0); // (corner,value, patch, unknown)
    BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
    BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
    BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

    Load = -1e1;
    // Point loads
    gsVector<> point(2);
    gsVector<> load (3);
    point<< 0.5, 0.5 ;
    load << 0.0, 0.0, Load ;
    pLoads.addLoad(point, load, 0 );

/* /////////////////////////////////////////////////////////////////////////////////////////
    Initial file IO

    blablabla

    Master process
///////////////////////////////////////////////////////////////////////////////////////// */

    std::string output = "solution";
    std::string dirname = "ArcLengthResults";

    dirname = dirname + "/" +  "Roof_t="+ std::to_string(thickness) + "-r=" + std::to_string(numHref) + "-e" + std::to_string(numElevate) +"_solution";
    output =  "solution";
    wn = "data.txt";
    std::string line = "line.txt";

    std::string commands = "mkdir -p " + dirname;
    const char *command = commands.c_str();
    system(command);

    // plot geometry
    if (plot)
      gsWriteParaview(mp,dirname + "/" + "mp",1000,mesh);

    if (write)
    {
      initStepOutput(dirname + "/" + wn, writePoints);
      initStepOutput(dirname + "/" + line, writePoints);
    }

/* /////////////////////////////////////////////////////////////////////////////////////////
    Initialize solvers

    blablabla

    ALL PROCESSES
///////////////////////////////////////////////////////////////////////////////////////// */

    // Initialise solution object
    gsMultiPatch<> mp_def = mp;

    gsMaterialMatrixBase<real_t>* materialMatrix;

    real_t pi = math::atan(1)*4;
    index_t kmax = 3;

    std::vector<gsFunctionSet<> * > Gs(kmax);
    std::vector<gsFunctionSet<> * > Ts(kmax);
    std::vector<gsFunctionSet<> * > Phis(kmax);

    gsMatrix<> Gmat = gsCompositeMatrix(Exx,Eyy,Gxy,PoissonRatio,PoissonRatio*Eyy/Exx);
    Gmat.resize(Gmat.rows()*Gmat.cols(),1);
    gsConstantFunction<> Gfun(Gmat,3);
    Gs[0] = Gs[1] = Gs[2] = &Gfun;

    gsConstantFunction<> phi1,phi2,phi3;
    phi1.setValue(pi/2,3);
    phi2.setValue(0,3);
    phi3.setValue(pi/2,3);

    Phis[0] = &phi1;
    Phis[1] = &phi2;
    Phis[2] = &phi3;

    gsConstantFunction<> thicks(thickness/kmax,3);
    Ts[0] = Ts[1] = Ts[2] = &thicks;

    gsConstantFunction<> force(tmp,3);
    gsFunctionExpr<> t(std::to_string(thickness), 3);
    gsFunctionExpr<> E(std::to_string(Exx),3);
    gsFunctionExpr<> nu(std::to_string(PoissonRatio),3);
    gsFunctionExpr<> rho(std::to_string(Density),3);

    std::vector<gsFunction<>*> parameters;
    parameters.resize(2);
    parameters[0] = &E;
    parameters[1] = &nu;

    gsOptionList options;
    if (composite)
    {
        materialMatrix = new gsMaterialMatrixComposite<3,real_t>(mp,Ts,Gs,Phis);
    }
    else
    {
        options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
        options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
        materialMatrix = getMaterialMatrix<3,real_t>(mp,t,parameters,rho,options);
    }

    gsThinShellAssemblerBase<real_t>* assembler;
    assembler = new gsThinShellAssembler<3, real_t, true >(mp,dbasis,BCs,force,materialMatrix);

    // Construct assembler object
    assembler->setOptions(opts);
    assembler->setPointLoads(pLoads);

    typedef std::function<gsSparseMatrix<real_t> (gsVector<real_t> const &)>                                Jacobian_t;
    typedef std::function<gsVector<real_t> (gsVector<real_t> const &, real_t, gsVector<real_t> const &) >   ALResidual_t;
    // Function for the Jacobian
    Jacobian_t Jacobian = [&assembler](gsVector<real_t> const &x)
    {
      gsMultiPatch<> mp_def;
      assembler->constructSolution(x,mp_def);
      assembler->assembleMatrix(mp_def);
      gsSparseMatrix<real_t> m = assembler->matrix();
      return m;
    };
    // Function for the Residual
    ALResidual_t ALResidual = [&assembler](gsVector<real_t> const &x, real_t lam, gsVector<real_t> const &force)
    {
      gsMultiPatch<> mp_def;
      assembler->constructSolution(x,mp_def);
      assembler->assembleVector(mp_def);
      gsVector<real_t> Fint = -(assembler->rhs() - force);
      gsVector<real_t> result = Fint - lam * force;
      return result; // - lam * force;
    };
    // Assemble linear system to obtain the force vector
    assembler->assemble();
    gsVector<> Force = assembler->rhs();

    gsArcLengthIterator<real_t> arcLength(Jacobian, ALResidual, Force);

    arcLength.options().setInt("Solver",0); // LDLT solver
    arcLength.options().setInt("BifurcationMethod",0); // 0: determinant, 1: eigenvalue
    arcLength.options().setInt("Method",method);
    arcLength.options().setReal("Length",dL);
    arcLength.options().setInt("AngleMethod",0); // 0: step, 1: iteration
    arcLength.options().setSwitch("AdaptiveLength",adaptive);
    arcLength.options().setInt("AdaptiveIterations",5);
    arcLength.options().setReal("Scaling",0.0);
    arcLength.options().setReal("Tol",tol);
    arcLength.options().setReal("TolU",tolU);
    arcLength.options().setReal("TolF",tolF);
    arcLength.options().setInt("MaxIter",maxit);
    arcLength.options().setSwitch("Verbose",true);
    arcLength.options().setReal("Relaxation",relax);
    if (quasiNewtonInt>0)
    {
      quasiNewton = true;
      arcLength.options().setInt("QuasiIterations",quasiNewtonInt);
    }
    arcLength.options().setSwitch("Quasi",quasiNewton);

    gsInfo<<arcLength.options();
    arcLength.applyOptions();
    arcLength.initialize();

    gsMultiPatch<> deformation = mp;

    // Make objects for previous solutions
    real_t Lguess,Lold, L0, Lref;
    gsMatrix<> Uguess,Uold, U0, Uref;
    Uold.setZero(Force.size(),1);
    U0.setZero(Force.size(),1);
    L0 = Lold = 0.0;

    gsMatrix<> solVector;
    real_t indicator = 0.0;
    arcLength.setIndicator(indicator); // RESET INDICATOR
    real_t dL0 = dL;


    typedef std::pair<gsVector<real_t>,real_t> solution_t;

    /*
      \a solutions is a container the for each level contains the solutions per point
      \a points contains all the points across levels in the format (level, U, lambda) -------------------------------> OVERKILL? WHY NEEDED?
      \a refPoints is a container that contains (level, U, lambda) of the points from which a refinement should START in level+1
      \a errors is a container that contains the error[l][i] e_i at the ith point of level l
    */

/* /////////////////////////////////////////////////////////////////////////////////////////
    Initialize hierarchy

    blablabla

    Master process
///////////////////////////////////////////////////////////////////////////////////////// */

    std::vector<solution_t> solutions;
    std::vector<real_t> times;
    real_t timescaling;
    real_t s = 0;

    index_t level = 0;
    gsInfo<<"------------------------------------------------------------------------------------\n";
    gsInfo<<"\t\t\tLevel "<<level<<" (dL = "<<dL<<") -- Coarse grid \n";
    gsInfo<<"------------------------------------------------------------------------------------\n";

    index_t stepi = step; // number of steps for level i
    stepi = step * (math::pow(2,level));

      gsDebugVar(assembler->matrix());


    // Add the undeformed solution
    solutions.push_back({U0,L0});
    times.push_back(s);
    // Add other solutions
    for (index_t k=1; k<stepi+1; k++)
    {
      s+=dL;
      gsInfo<<"Load step "<< k<<"\t"<<"dL = "<<dL<<"; curve time = "<<s<<"\n";
      // assembler->constructSolution(solVector,solution);
      arcLength.step();

      // gsInfo<<"m_U = "<<arcLength.solutionU()<<"\n";
      if (!(arcLength.converged()))
        GISMO_ERROR("Loop terminated, arc length method did not converge.\n");

      real_t lambda = arcLength.solutionL();
      solutions.push_back({arcLength.solutionU(),lambda});
      times.push_back(s);
    }


    gsSpaceTimeHierarchy<real_t,solution_t> hierarchy(times,solutions);
    hierarchy.options().setInt("MaxLevel",maxLevel);

    gsSpaceTimeHierarchy<real_t,solution_t> hierarchy2(times);
    hierarchy2.printQueue();

    hierarchy2.getJob();


//     /////////////////////////////////////////////////////////////////////////////////////////////
//     /////////////////////////////////////////////////////////////////////////////////////////////
//     /////////////////////////////////////////////////////////////////////////////////////////////

//     // Store solution coefficients in a matrix
//     index_t blocksize = mp.patch(0).coefs().rows();
//     gsMatrix<> solutionCoefs(solutions.size(),3*blocksize);
//     gsVector<> loads(solutions.size());
//     gsMultiPatch<> mp_tmp;

//     for (index_t k=0; k!= solutions.size(); k++)
//     {
//       assembler->constructSolution(solutions[k].first,mp_tmp);
//       solutionCoefs.row(k) = mp_tmp.patch(0).coefs().reshape(1,3*blocksize);

//       loads.at(k) = solutions[k].second;
//     }

//     gsTensorBSpline<3,real_t> fit = gsSpaceTimeFit<3,real_t>(solutionCoefs,loads,gsAsVector<>(times),dbasis,deg_z);

//     typename gsTensorBSpline<3,real_t>::BoundaryGeometryType target;

//     gsParaviewCollection collection(dirname + "/" + output);
//     gsParaviewCollection datacollection(dirname + "/" + "data");

//     gsField<> solField;

//     if (plot || write)
//     {
//       gsVector<> xi;
//       xi.setLinSpaced(100,times[0],times[times.size()-1]);

//       for (index_t k = 0; k!=xi.size(); k++)
//       {
//         fit.slice(2,xi.at(k),target);
//         gsGeometry<real_t> * slice = target.clone().release();
//         real_t lambda = slice->coefs()(0,3);
//         slice->embed(3);

//         deformation.patch(0) = *slice;
//         deformation.patch(0).coefs() -= mp.patch(0).coefs();// assuming 1 patch here

//         if (plot)
//         {
//           solField = gsField<>(mp,deformation);
//           gsWriteParaview(solField,"slice");

//           std::string fileName = dirname + "/" + output + util::to_string(k);
//           gsWriteParaview<>(solField, fileName, 1000,mesh);
//           fileName = output + util::to_string(k) + "0";
//           collection.addTimestep(fileName,xi[k],".vts");
//           if (mesh) collection.addTimestep(fileName,k,"_mesh.vtp");
//         }
//         if (write)
//         {
//             writeStepOutput(lambda,deformation, dirname + "/" + line, writePoints,1, 201);
//         }
//       }

//       {
//         for (index_t k=0; k!= solutions.size(); k++)
//         {
//           assembler->constructSolution(solutions[k].first,mp_tmp);

//           real_t lambda = solutions[k].second;

//           real_t Time = times[k];

//           deformation.patch(0) = mp_tmp.patch(0);
//           deformation.patch(0).coefs() -= mp.patch(0).coefs();// assuming 1 patch here

//           if (plot)
//           {
//             solField = gsField<>(mp,deformation);
//             std::string fileName = dirname + "/" + "data" + util::to_string(k);
//             gsWriteParaview<>(solField, fileName, 1000,mesh);
//             fileName = "data" + util::to_string(k) + "0";
//             datacollection.addTimestep(fileName,Time,".vts");
//             if (mesh) datacollection.addTimestep(fileName,k,"_mesh.vtp");
//           }
//           if (write)
//           {
//               writeStepOutput(lambda,deformation, dirname + "/" + wn, writePoints,1, 201);
//           }
//         }
//       }

//       if (plot)
//       {
//         collection.save();
//         datacollection.save();
//       }
//     }
//     /////////////////////////////////////////////////////////////////////////////////////////////
//     /////////////////////////////////////////////////////////////////////////////////////////////
//     /////////////////////////////////////////////////////////////////////////////////////////////

//     gsSpaceTimeHierarchy<real_t,solution_t> hierarchy(times,solutions,maxLevel);

//     hierarchy.printTree();

//     solution_t start, guess, reference;
//     real_t ttmp = 0;
//     real_t dLtmp;
//     index_t it = 0;
//     index_t itmax = 100;
//     real_t TOL = 1e-2;
//     while (!hierarchy.empty() && it < itmax)
//     {
//       std::tie(ttmp,dLtmp,start,guess) = hierarchy.getJob();

//       gsDebugVar(ttmp);
//       gsDebugVar(dLtmp);

//       std::tie(Uold,Lold) = start;
//       std::tie(Uguess,Lguess) = guess;

//       arcLength.setLength(dLtmp);
//       arcLength.setSolution(Uold,Lold);
//       arcLength.resetStep();

//       arcLength.setInitialGuess(Uguess,Lguess);

//       gsInfo<<"Starting from (lvl,|U|,L) = ("<<hierarchy.currentLevel()<<","<<Uold.norm()<<","<<Lold<<"), curve time = "<<ttmp<<"\n";
//       arcLength.step();

//       bool success = hierarchy.reference_into(reference);
//       if (success)
//       {
//         std::tie(Uref,Lref) = reference;
//       }
//       else
//       {
//         fit.slice(2,ttmp+dLtmp,target);
//         gsGeometry<real_t> * slice = target.clone().release();
//         Lref  = slice->coefs()(0,3);
//         slice->embed(3);
//         gsMultiPatch<> mp_tmp2(*slice);
//         mp_tmp2.patch(0).coefs() -= mp.patch(0).coefs();
//         Uref = assembler->constructSolutionVector(mp_tmp2);
//       }

//       gsVector<> DeltaU = Uref - arcLength.solutionU();
//       real_t DeltaL = Lref - arcLength.solutionL();

//       real_t error = arcLength.distance(DeltaU,DeltaL) / (dLtmp);

//       hierarchy.submitJob(std::make_pair(arcLength.solutionU(),arcLength.solutionL()));

//       if (error > TOL)
//       {
//         gsDebug<<"ERROR > TOL\n";
//         hierarchy.addJobHere();
//       }

//       hierarchy.removeJob();

//       it++;
//     }

//     std::tie(times,solutions) = hierarchy.getFlatSolution();

//     gsDebugVar(gsAsVector(times));

//     /////////////////////////////////////////////////////////////////////////////////////////////
//     /////////////////////////////////////////////////////////////////////////////////////////////
//     /////////////////////////////////////////////////////////////////////////////////////////////

//     gsParaviewCollection collection2(dirname + "/" + output + "_refit");
//     gsParaviewCollection datacollection2(dirname + "/" + "data" + "_refit");


//     // Store solution coefficients in a matrix
//     blocksize = mp.patch(0).coefs().rows();
//     solutionCoefs = gsMatrix<> (solutions.size(),3*blocksize);
//     loads = gsVector<>(solutions.size());

//     for (index_t k=0; k!= solutions.size(); k++)
//     {
//       assembler->constructSolution(solutions[k].first,mp_tmp);
//       solutionCoefs.row(k) = mp_tmp.patch(0).coefs().reshape(1,3*blocksize);

//       loads.at(k) = solutions[k].second;
//     }

//     gsTensorBSpline<3,real_t> fit2 = gsSpaceTimeFit<3,real_t>(solutionCoefs,loads,gsAsVector<>(times),dbasis,deg_z);

//     if (plot || write)
//     {
//       gsVector<> xi;
//       xi.setLinSpaced(100,times[0],0.99999*times[times.size()-1]);

//       for (index_t k = 0; k!=xi.size(); k++)
//       {
//         fit2.slice(2,xi.at(k),target);
//         gsGeometry<real_t> * slice = target.clone().release();
//         real_t lambda = slice->coefs()(0,3);
//         slice->embed(3);

//         deformation.patch(0) = *slice;
//         deformation.patch(0).coefs() -= mp.patch(0).coefs();// assuming 1 patch here

//         if (plot)
//         {
//           solField = gsField<>(mp,deformation);
//           gsWriteParaview(solField,"slice");

//           std::string fileName = dirname + "/" + output + util::to_string(k);
//           gsWriteParaview<>(solField, fileName, 1000,mesh);
//           fileName = output + util::to_string(k) + "0";
//           collection2.addTimestep(fileName,xi[k],".vts");
//           if (mesh) collection2.addTimestep(fileName,k,"_mesh.vtp");
//         }
//         if (write)
//         {
//             writeStepOutput(lambda,deformation, dirname + "/" + line, writePoints,1, 201);
//         }
//       }

//       {
//         for (index_t k=0; k!= solutions.size(); k++)
//         {
//           assembler->constructSolution(solutions[k].first,mp_tmp);

//           real_t lambda = solutions[k].second;

//           real_t Time = times[k];

//           deformation.patch(0) = mp_tmp.patch(0);
//           deformation.patch(0).coefs() -= mp.patch(0).coefs();// assuming 1 patch here

//           if (plot)
//           {
//             solField = gsField<>(mp,deformation);
//             std::string fileName = dirname + "/" + "data" + util::to_string(k);
//             gsWriteParaview<>(solField, fileName, 1000,mesh);
//             fileName = "data" + util::to_string(k) + "0";
//             datacollection2.addTimestep(fileName,Time,".vts");
//             if (mesh) datacollection2.addTimestep(fileName,k,"_mesh.vtp");
//           }
//           if (write)
//           {
//               writeStepOutput(lambda,deformation, dirname + "/" + wn, writePoints,1, 201);
//           }
//         }
//       }

//       if (plot)
//       {
//         collection2.save();
//         datacollection2.save();
//       }
//     }

  return result;
}

template <class T>
void initStepOutput(const std::string name, const gsMatrix<T> & points)
{
  std::ofstream file;
  file.open(name,std::ofstream::out);
  file  << std::setprecision(20)
        << "Deformation norm" << ",";
        for (index_t k=0; k!=points.cols(); k++)
        {
          file<< "point "<<k<<" - x" << ","
              << "point "<<k<<" - y" << ","
              << "point "<<k<<" - z" << ",";
        }

  file  << "Lambda" << ","
        << "Indicator"
        << "\n";
  file.close();

  gsInfo<<"Step results will be written in file: "<<name<<"\n";
}

template <class T>
void writeStepOutput(const T lambda, const gsMultiPatch<T> & deformation, const std::string name, const gsMatrix<T> & points, const index_t extreme, const index_t kmax) // extreme: the column of point indices to compute the extreme over (default -1)
{
  gsMatrix<T> P(2,1), Q(2,1);
  gsMatrix<T> out(3,points.cols());
  gsMatrix<T> tmp;

  for (index_t p=0; p!=points.cols(); p++)
  {
    P<<points.col(p);
    deformation.patch(0).eval_into(P,tmp);
    out.col(p) = tmp;
  }

  std::ofstream file;
  file.open(name,std::ofstream::out | std::ofstream::app);
  if (extreme==-1)
  {
    file  << std::setprecision(6)
          << "NA" << ",";
          for (index_t p=0; p!=points.cols(); p++)
          {
            file<< out(0,p) << ","
                << out(1,p) << ","
                << out(2,p) << ",";
          }

    file  << lambda << ","
          << "NA" << ","
          << "\n";
  }
  else if (extreme==0 || extreme==1)
  {
    gsMatrix<T> out2(kmax,points.cols()); // evaluation points in the rows, output (per coordinate) in columns
    for (int p = 0; p != points.cols(); p ++)
    {
      Q.at(1-extreme) = points(1-extreme,p);
      for (int k = 0; k != kmax; k ++)
      {
        Q.at(extreme) = 1.0*k/(kmax-1);
        deformation.patch(0).eval_into(Q,tmp);
        out2(k,p) = tmp.at(2); // z coordinate
      }
    }

    file  << std::setprecision(6)
          << "NA" << ",";
          for (index_t p=0; p!=points.cols(); p++)
          {
            file<< out(0,p) << ","
                << out(1,p) << ","
                << std::max(abs(out2.col(p).maxCoeff()),abs(out2.col(p).minCoeff())) << ",";
          }

    file  << lambda << ","
          << "NA" << ","
          << "\n";
  }
  else
    GISMO_ERROR("Extremes setting unknown");

  file.close();
}

template<int dim, class T>
gsTensorBSpline<dim,T> gsSpaceTimeFit(const gsMatrix<T> & solutionCoefs, const gsVector<T> & times, const gsVector<T> & ptimes, gsMultiBasis<T> & spatialBasis, index_t deg)
{
  GISMO_ASSERT(solutionCoefs.rows()==times.rows(),"Number of time and solution steps should match! "<<solutionCoefs.rows()<<"!="<<times.rows());
  GISMO_ASSERT(solutionCoefs.cols() % dim == 0,"Is the dimension correct?"<<solutionCoefs.cols()<<"!="<<dim);
  index_t nsteps = times.rows();
  index_t bsize = solutionCoefs.cols()/dim;

  // Prepare fitting basis
  gsKnotVector<> kv(ptimes.minCoeff(),ptimes.maxCoeff(),nsteps-(deg+1),deg+1);
  gsBSplineBasis<T> lbasis(kv);


  //////// TO DO:
  //////// - Make multi-patch compatible
  //////// - Include dimension d

  // for (index_t p = 0; p!=dbasis.nBases(); ++p)
  // {
    gsTensorBSplineBasis<dim,T> tbasis(
                                            static_cast<gsBSplineBasis<T> *>(&spatialBasis.basis(0).component(0))->knots(),
                                            static_cast<gsBSplineBasis<T> *>(&spatialBasis.basis(0).component(1))->knots(),
                                            kv
                                            );
  // }

  gsMatrix<> rhs(times.size(),(dim+1)*bsize);
  gsVector<> ones; ones.setOnes(bsize);

  for (index_t lam = 0; lam!=nsteps; ++lam)
  {
    rhs.block(lam,0,1,dim * bsize) = solutionCoefs.row(lam);
    rhs.block(lam,dim*bsize,1,bsize) = times.at(lam) * ones.transpose();
  }

  // get the Greville Abcissae (anchors)
  gsMatrix<> anchors = lbasis.anchors();

  // Get the collocation matrix at the anchors
  gsSparseMatrix<> C;
  lbasis.collocationMatrix(anchors,C);

  gsSparseSolver<>::LU solver;
  solver.compute(C);

  gsMatrix<> sol, coefs((nsteps)*bsize,dim+1);
  sol = solver.solve(rhs);

  for (index_t lam = 0; lam!=nsteps; ++lam)
  {
    gsMatrix<> tmp = sol.block(lam,0,1,dim * bsize);
    coefs.block(lam * bsize,0,bsize,dim) = tmp.reshape(bsize,dim);
    coefs.block(lam * bsize,dim,bsize,1) = sol.block(lam,dim*bsize,1,bsize).transpose();
  }

  // gsTensorBSpline<3,T> tspline = tbasis.makeGeometry(give(coefs)).release();
  gsTensorBSpline<dim,T> tspline(tbasis,give(coefs));
  return tspline;
}

template<index_t dim, class T>
class gsSpaceTimeFitter
{
public:
  gsSpaceTimeFitter  ( const gsMatrix<T> & solutionCoefs,
                      const gsVector<T> & times,
                      const gsVector<T> & ptimes,
                      const gsMultiBasis<T> & spatialBasis,
                      const index_t deg = 2)
  :
  m_data(solutionCoefs),
  m_times(times),
  m_ptimes(ptimes),
  m_bases(spatialBasis),
  m_deg(deg)
  {

  }

  gsTensorBSpline<dim,T> compute()
  {
    GISMO_ASSERT(m_data.rows()==m_times.rows(),"Number of time and solution steps should match! "<<m_data.rows()<<"!="<<m_times.rows());
    GISMO_ASSERT(m_data.cols() % dim == 0,"Is the dimension correct?"<<m_data.cols()<<"!="<<dim);
    index_t nsteps = m_times.rows();
    index_t bsize = m_data.cols()/dim;

    // Prepare fitting basis
    gsKnotVector<> kv(m_ptimes.minCoeff(),m_ptimes.maxCoeff(),nsteps-(m_deg+1),m_deg+1);
    gsBSplineBasis<T> lbasis(kv);


    //////// TO DO:
    //////// - Make multi-patch compatible
    //////// - Include dimension d

    // for (index_t p = 0; p!=dbasis.nBases(); ++p)
    // {
      gsTensorBSplineBasis<dim,T> tbasis(
                                              static_cast<gsBSplineBasis<T> *>(&m_bases.basis(0).component(0))->knots(),
                                              static_cast<gsBSplineBasis<T> *>(&m_bases.basis(0).component(1))->knots(),
                                              kv
                                              );
    // }

    gsMatrix<> rhs(m_times.size(),(dim+1)*bsize);
    gsVector<> ones; ones.setOnes(bsize);

    for (index_t lam = 0; lam!=nsteps; ++lam)
    {
      rhs.block(lam,0,1,dim * bsize) = m_data.row(lam);
      rhs.block(lam,dim*bsize,1,bsize) = m_times.at(lam) * ones.transpose();
    }

    // get the Greville Abcissae (anchors)
    gsMatrix<> anchors = lbasis.anchors();

    // Get the collocation matrix at the anchors
    gsSparseMatrix<> C;
    lbasis.collocationMatrix(anchors,C);

    gsSparseSolver<>::LU solver;
    solver.compute(C);

    gsMatrix<> sol, coefs((nsteps)*bsize,dim+1);
    sol = solver.solve(rhs);

    for (index_t lam = 0; lam!=nsteps; ++lam)
    {
      gsMatrix<> tmp = sol.block(lam,0,1,dim * bsize);
      coefs.block(lam * bsize,0,bsize,dim) = tmp.reshape(bsize,dim);
      coefs.block(lam * bsize,dim,bsize,1) = sol.block(lam,dim*bsize,1,bsize).transpose();
    }

    // gsTensorBSpline<3,T> tspline = tbasis.makeGeometry(give(coefs)).release();
    gsTensorBSpline<dim,T> tspline(tbasis,give(coefs));
    return tspline;
  }

protected:
  gsMatrix<T> m_data;
  gsVector<T> m_times;
  gsVector<T> m_ptimes;
  gsMultiBasis<T> m_bases;
  index_t m_deg;

};
