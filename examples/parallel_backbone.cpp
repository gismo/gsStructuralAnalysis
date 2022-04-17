/** @file mpi_example.cpp

    @brief Testing MPI with G+Smo

    Execute (eg. with 10 processes):
       mpirun -np 10 ./bin/mpi_example

    or provide a hosts file on a cluster:
       mpirun -hostfile hosts.txt ./bin/mpi_example

    If your cluster is using srun:
       srun -N 10 ./bin/mpi_example

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, C. Hofer, R. Schneckenleitner
*/

#include <gismo.h>
#include <time.h>
#include <unistd.h>

real_t                    next_random

  ( void )

{
  static index_t  initialized = 0;
  index_t         next;

  if ( ! initialized )
  {
    index_t  my_rank;
    index_t  flag;

    MPI_Initialized ( &flag );

    if ( flag )
    {
      MPI_Comm_rank ( MPI_COMM_WORLD, &my_rank );

      srand ( (unsigned index_t) my_rank );
    }

    initialized = 1;
  }

  next = rand ();

  return ((real_t) next / (real_t) RAND_MAX);
}

using namespace gismo;

index_t main(index_t argc, char **argv)
{
  const index_t  N    = 100;
  const index_t  root = 0;
  const index_t  tag  = 1;

  real_t     number;
  real_t     max_num;
  real_t     recv_num;
  index_t        proc_count;
  index_t        my_rank;
  index_t        iiter;

  gsCmdLine cmd("An example for testing MPI with G+Smo.\n"
      "Execute (eg. with 10 processes):                                      "
      "  *  mpirun -np 10 ./bin/mpi_example\n"
      "or provide a hosts file on a cluster:                                 "
      "  *  mpirun -hostfile hosts.txt ./bin/mpi_example\n"
      "If your cluster is using srun:                                        "
      "  *  srun -N 10 ./bin/mpi_example"
  );
  try { cmd.getValues(argc,argv); } catch (index_t rv) { return rv; }

/////////////////////////////MPI/////////////////////////////
  // Conditional compilation
#ifdef GISMO_WITH_MPI
  gsInfo << "Gismo was compiled with MPI support.\n";
#else
  gsInfo << "Gismo was compiled without MPI support.\n";
#endif


  // Initialize the MPI environment
  const gsMpi & mpi = gsMpi::init(argc, argv);

  // Get current wall time
  real_t wtime = mpi.wallTime();

  // Get the world communicator
  gsMpiComm comm = mpi.worldComm();
  MPI_Request req;

  //Get size and rank of the processor
  proc_count = comm.size();
  my_rank = comm.rank();

  GISMO_ASSERT(proc_count > 1,"At least two processes are required.\n");

/////////////////////////////////////////////////////////////

  gsInfo<<"[MPI process "<<my_rank<<"] Problem initialization...";
  real_t maxSleep = 1;
  real_t sleepTime = maxSleep * next_random();

  sleep(sleepTime);
  gsInfo<<" ("<<sleepTime<<" seconds)\n";



  // MPI
  index_t njobs = 0;
  std::queue<index_t> m_workers;
  real_t tstart = 0.0;
  real_t tend   = 1.0;
  index_t N0        = 10;
  real_t dt0    = (tend - tstart) / N0;

  real_t x, y, z;

  gsMpiStatus status;

  index_t globalID = 0;

  /*
      The example works as follows:
      The master process has a list of xy-values for which a z-value should be known.
      x and y are functions of t and the workers will compute the z-values.
      The system is
      x(t) = cos(t)
      y(t) = sin(t)
      z(t) = x(t)*y(t)

      The guess will be computed as follows
      x(t+dt) = cos(t+dt) \approx -cos(t)+dt*sin(t)
      y(t+dt) = sin(t+dt) \approx  sin(t)+dt*cos(t)
      z(t+dt) \approx (-cos(t)+dt*sin(t))*sin(t)+dt*cos(t)

      Therefore, the data that is sent from the master process contains:

      send_typle_t <level,ID,x,stop>
      recv_typle_t <level,ID,y>

   */
  // typedef std::tuple<index_t,index_t,real_t,bool> send_tuple_t;
  // typedef std::tuple<index_t,index_t,real_t,real_t> recv_tuple_t;

  typedef std::pair<gsVector<real_t,3>,real_t> solution_t; // fix size to real_t,3
  //                 ID, start time, time step, start, guess, stop?
  typedef std::tuple<index_t,real_t,real_t,real_t,solution_t,solution_t,bool> send_tuple_t;
  // ID, solutions, times, distances
  typedef std::tuple<index_t,std::vector<solution_t>,std::vector<real_t>,std::vector<real_t>> recv_tuple_t;

  // !MPI



  ////////////////////////////////////////////////////////////////////////////////////////////////
  // INIT
  ////////////////////////////////////////////////////////////////////////////////////////////////
  //                    ID,      tstart, tend,   dt,     start     , guess
  std::queue<std::tuple<index_t, real_t, real_t, real_t, solution_t, solution_t,bool>> m_queue; // xi,exact



  bool message = true;

  if ( my_rank == 0 )
  {
    send_tuple_t send;
    recv_tuple_t recv;

    gsInfo<<"[MPI process "<<my_rank<<"] Making queue ...";
    gsVector<real_t> coords(3);
    real_t angle, nextangle, dangle;
    real_t t, tpp, dt;
    solution_t solution, refsolution;

    std::tuple<index_t, real_t, real_t, real_t, solution_t, solution_t, bool> item;
    for (index_t k = 0; k!= N0+1; k++)
    {
      dt = dt0;
      t = k*dt;

      angle = t * 2 * M_PI;
      x = math::cos(angle);
      y = math::sin(angle);
      z = x*y;
      coords<<x,y,z;
      solution = std::make_pair(coords,angle);

      tpp = t+dt;
      nextangle = tpp * 2 * M_PI;
      dangle = nextangle-angle;
      x = math::cos(nextangle) - math::sin(nextangle)*dangle;
      y = math::sin(nextangle) + math::cos(nextangle)*dangle;
      z = x * y;
      coords<<x,y,z;
      refsolution = std::make_pair(coords,nextangle);

      item = std::make_tuple(globalID,t, tpp, dt0, solution, refsolution,false);

      m_queue.push(item);
      globalID++;
    }
    gsInfo<<"finished.\n";

    gsInfo<<"[MPI process "<<my_rank<<"] Adding workers...\n";
    for (index_t w = 1; w!=proc_count; w++)
      m_workers.push(w);

    while (!m_queue.empty() && !m_workers.empty())
    {
      send = m_queue.front();
      m_queue.pop();

      gsInfo<<"[MPI process "<<my_rank<<"] Sending a job to "<<m_workers.front()<<": ID: "<<std::get<0>(send)
                                      <<"; t = "<<std::get<1>(send)
                                      <<"; t+1 = "<<std::get<2>(send)
                                      <<"; dt = "<<std::get<3>(send)
                                      <<"; x(t) = "<<std::get<4>(send).first.at(0)
                                      <<"; y(t) = "<<std::get<4>(send).first.at(1)
                                      <<"; z(t) = "<<std::get<4>(send).first.at(2)
                                      <<"; x'(t) = "<<std::get<5>(send).first.at(0)
                                      <<"; y'(t) = "<<std::get<5>(send).first.at(1)
                                      <<"; z'(t) = "<<std::get<5>(send).first.at(2)
                                      <<"\n";


      comm.isend(&send, 1, m_workers.front(),&req,tag);
        gsDebugVar(std::get<4>(send).first);
      m_workers.pop();
      njobs++;
    }

    while (njobs > 0)
    {
      gsInfo<<"[MPI process "<<my_rank<<"] "<<njobs<<" job(s) running\n";
  gsDebugVar("hi");
      comm.recv(&recv,1,MPI_ANY_SOURCE,tag,&status);
  gsDebugVar("hi");
      gsInfo<<"[MPI process "<<my_rank<<"] Received a job from "<<status.MPI_SOURCE<<": ID: "<<std::get<0>(recv)<<"; solutions.size() = "<<std::get<1>(recv).size()<<"; times.size() = "<<std::get<2>(recv).size()<<"; distances.size() = "<<std::get<3>(recv).size()<<"\n";
      njobs--;
      m_workers.push(status.MPI_SOURCE);

      // if (std::get<3>(recv) < 0.75 && std::get<3>(recv) > 0.25)
      // {
      //   gsInfo<<"[MPI process "<<my_rank<<"] special!"<<"\n";
      //   // m_queue.push({l,k,std::pow(2,-l)*dt0*k})
      // }

      while (!m_queue.empty() && !m_workers.empty())
      {
        // push
        send = m_queue.front();

        m_queue.pop();
        // Question 1a: How to send the tuple <index_t,index_t,real_t> to slave?
        // Question 2 : How to send the job to ANY slave (i.e. the first available one)

        gsInfo<<"[MPI process "<<my_rank<<"] Sending a job to "<<m_workers.front()<<": ID: "<<std::get<0>(send)
                                        <<"; t = "<<std::get<1>(send)
                                        <<"; t+1 = "<<std::get<2>(send)
                                        <<"; dt = "<<std::get<3>(send)
                                        <<"; x(t) = "<<std::get<4>(send).first.at(0)
                                        <<"; y(t) = "<<std::get<4>(send).first.at(1)
                                        <<"; z(t) = "<<std::get<4>(send).first.at(2)
                                        <<"; x'(t) = "<<std::get<5>(send).first.at(0)
                                        <<"; y'(t) = "<<std::get<5>(send).first.at(1)
                                        <<"; z'(t) = "<<std::get<5>(send).first.at(2)
                                        <<"\n";
        comm.isend(&send, 1, m_workers.front(),&req,tag);
        m_workers.pop();
        njobs++;
      }
    }

    std::get<6>(send) = true;
    for (index_t w = 1; w!=proc_count; w++)
      comm.isend(&send, 1,w,&req,tag);
  }
  else
  {
    while (true)
    {
  gsDebugVar("hi");
      recv_tuple_t send;
      send_tuple_t recv;
      index_t Nintervals = 2;
      std::vector<real_t> times(Nintervals);
      std::vector<solution_t> solutions(Nintervals);
      std::vector<real_t> distances(Nintervals+1);

      gsVector<real_t> coords(3), coordsprev(3);
      real_t x, y, z, xpp, ypp, zpp, xref, yref, zref, t, tpp, dt0, dt;
      real_t angle, nextangle;

  gsDebugVar("hi");
      comm.recv(&recv,1,0,tag,MPI_STATUS_IGNORE);
  gsDebugVar("hi");

      gsDebugVar(std::get<4>(recv).second);
      gsDebugVar(std::get<4>(recv).first);

      if (std::get<6>(recv))
      {
        gsInfo<<"[MPI process "<<my_rank<<"] I have to stop!!"<<"\n";
        break;
      }
  gsDebugVar("hi");

      t = std::get<1>(recv);
  gsDebugVar("hi");
      tpp = std::get<2>(recv);
  gsDebugVar("hi");
      dt0 = std::get<3>(recv);
  gsDebugVar("hi");
  coordsprev = std::get<4>(recv).first;
  gsDebugVar("hi");
      x = coordsprev.at(0);
  gsDebugVar("hi");
      y = coordsprev.at(1);
  gsDebugVar("hi");
      z = coordsprev.at(2);
  gsDebugVar("hi");
  coords = std::get<5>(recv).first;
  gsDebugVar("hi");
      xref = coords.at(0);
  gsDebugVar("hi");
      yref = coords.at(1);
  gsDebugVar("hi");
      zref = coords.at(2);
  gsDebugVar("hi");
      gsInfo<<"[MPI process "<<my_rank<<"] Received a job from "<<0<<": ID: "<<std::get<0>(recv)<<"; t = "<<t<<"; t+1 = "<<tpp<<"; dt = "<<dt<<"; x(t) = "<<x<<"; y(t) = "<<y<<"; z(t) = "<<z<<"; x'(t) = "<<xref<<"; y'(t) = "<<yref<<"; z'(t) = "<<zref<<"\n";
  gsDebugVar("hi");

      coordsprev<<x,y,z;

      dt = dt0/Nintervals;
      for (index_t k = 0; k!=Nintervals; k++)
      {
        gsDebugVar(t);
        gsDebugVar(tpp);
        gsDebugVar(t+dt*k);
        nextangle = (t + dt*k) * 2 * M_PI;
        xpp = math::cos(nextangle);
        ypp = math::sin(nextangle);
        zpp = xpp*ypp;
        coords<<xpp,ypp,zpp;

        times.at(k) = t+k*dt;
        solutions.at(k) = std::make_pair(coords,nextangle);
        distances.at(k) = (coords-coordsprev).norm();

        coordsprev = coords;
      }
      coords<<xref,yref,zref;
      distances.at(Nintervals) = (coords-coordsprev).norm();


      std::get<0>(send) = std::get<0>(recv); // ID
      std::get<1>(send) = solutions; // ID
      std::get<2>(send) = times; // ID
      std::get<3>(send) = distances; // ID

      real_t maxSleep = 1;
      real_t sleepTime = maxSleep * next_random();

      sleep(sleepTime);
      gsInfo<<"[MPI process "<<my_rank<<"] slept for "<<sleepTime<<" seconds\n";

      gsInfo<<"[MPI process "<<my_rank<<"] Sending a job to    "<<0<<": ID: "<<std::get<0>(send)<<"; solutions.size() = "<<std::get<1>(send).size()<<"; times.size() = "<<std::get<2>(send).size()<<"; distances.size() = "<<std::get<3>(send).size()<<"\n";
      comm.send(&send,1,0,tag);
    }
  }


  return 0;
}