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

  typedef std::pair<gsVector<real_t>,real_t> solution_t; // fix size to real_t,3
  //                 ID, start time, time step, start, guess, stop?
  typedef std::tuple<index_t,real_t,real_t,real_t,gsVector<real_t>,real_t,gsVector<real_t>,real_t,bool> send_tuple_t;
  // ID, solutions, times, distances
  typedef std::tuple<index_t,std::vector<solution_t>,std::vector<real_t>,std::vector<real_t>> recv_tuple_t;

  // !MPI



  ////////////////////////////////////////////////////////////////////////////////////////////////
  // INIT
  ////////////////////////////////////////////////////////////////////////////////////////////////
  //                    ID,      tstart, tend,   dt,     start                   , guess
  std::queue<std::tuple<index_t, real_t, real_t, real_t, gsVector<real_t>, real_t, gsVector<real_t>, real_t,bool>> m_queue; // xi,exact



  bool message = true;
  gsVector<real_t> coords(3), nextcoords(3);
  real_t angle, nextangle, dangle;
  real_t t, tpp, dt;
  bool stop;
  index_t ID;

  std::vector<real_t>           times;
  std::vector<gsVector<real_t>> solutions;
  std::vector<real_t>           angles;
  std::vector<real_t>           distances;

  if ( my_rank == 0 )
  {
    send_tuple_t send;
    recv_tuple_t recv;

    gsInfo<<"[MPI process "<<my_rank<<"] Making queue ...";

    std::tuple<index_t, real_t, real_t, real_t, gsVector<real_t>, real_t, gsVector<real_t>, real_t, bool> item;
    for (index_t k = 0; k!= N0+1; k++)
    {
      dt = dt0;
      t = k*dt;

      angle = t * 2 * M_PI;
      x = math::cos(angle);
      y = math::sin(angle);
      z = x*y;
      coords<<x,y,z;

      tpp = t+dt;
      nextangle = tpp * 2 * M_PI;
      dangle = nextangle-angle;
      x = math::cos(nextangle) - math::sin(nextangle)*dangle;
      y = math::sin(nextangle) + math::cos(nextangle)*dangle;
      z = x * y;
      nextcoords<<x,y,z;

      item = std::make_tuple(globalID,t, tpp, dt0, coords, angle, nextcoords, nextangle,false);

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

      std::tie(ID,t, tpp, dt0, coords, angle, nextcoords, nextangle,stop) = send;

      gsInfo<<"[MPI process "<<my_rank<<"] Sending a job to "<<m_workers.front()<<": ID: "<<ID
                                      <<"; t = "<<t
                                      <<"; t+1 = "<<tpp
                                      <<"; dt = "<<dt0
                                      <<"; x(t) = "<<coords.at(0)
                                      <<"; y(t) = "<<coords.at(1)
                                      <<"; z(t) = "<<coords.at(2)
                                      <<"; x'(t) = "<<nextcoords.at(0)
                                      <<"; y'(t) = "<<nextcoords.at(1)
                                      <<"; z'(t) = "<<nextcoords.at(2)
                                      <<"\n";


      comm.isend(&ID, 1, m_workers.front(),&req,0);
      comm.isend(&t, 1, m_workers.front(),&req,1);
      comm.isend(&tpp, 1, m_workers.front(),&req,2);
      comm.isend(&dt0, 1, m_workers.front(),&req,3);
      comm.isend(coords.data(), coords.size(), m_workers.front(),&req,4);
      comm.isend(&angle, 1, m_workers.front(),&req,5);
      comm.isend(nextcoords.data(), nextcoords.size(), m_workers.front(),&req,6);
      comm.isend(&nextangle, 1, m_workers.front(),&req,7);
      comm.isend(&stop, 1, m_workers.front(),&req,8);
      m_workers.pop();
      njobs++;
    }

    index_t size;
    while (njobs > 0)
    {
      gsInfo<<"[MPI process "<<my_rank<<"] "<<njobs<<" job(s) running\n";

      index_t SAME_SOURCE = MPI_ANY_SOURCE;
      // ID
      comm.recv(&ID,    1,SAME_SOURCE,0,&status);
      SAME_SOURCE = status.MPI_SOURCE;

      // Size of solutions
      comm.recv(&size,  1,SAME_SOURCE,1,&status);

      times.resize(size);
      solutions.resize(size);
      angles.resize(size);
      distances.resize(size+1);

      for (index_t k=0; k!=size; k++)
      {
        comm.recv(&(times.at(k))    ,1,SAME_SOURCE,4*k+2+0,&status);
        index_t solsize;
        comm.recv(&solsize,  1,SAME_SOURCE,4*k+2+1);
        gsVector<> solution(solsize);
        comm.recv(solution.data(),solsize,SAME_SOURCE,4*k+2+2,&status);
        solutions.at(k) = solution;
        comm.recv(&(angles.at(k))   ,1,SAME_SOURCE,4*k+2+3,&status);
        comm.recv(&(distances.at(k)),1,SAME_SOURCE,4*k+2+4,&status);
      }
      comm.recv(&(distances.at(size)),1,SAME_SOURCE,4*size+2+0,&status);

      gsInfo<<"[MPI process "<<my_rank<<"] Received a job from "<<status.MPI_SOURCE<<": ID: "<<ID<<"; size = "<<size<<"\n";
      gsInfo<<"\t\t"<<"time\tsolution\tangle\tdistance\n";
      for (index_t k=0; k!=solutions.size(); k++)
        gsInfo<<"\t\t"<<times.at(k)<<"\t"<<solutions.at(k).transpose()<<"\t"<<angles.at(k)<<"\t"<<distances.at(k)<<"\n";
      gsInfo<<"\t\t\t\t\t\t\t\t"<<distances.back()<<"\n";

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

        index_t ID;
        bool stop;
        std::tie(ID,t, tpp, dt0, coords, angle, nextcoords, nextangle,stop) = send;

        gsInfo<<"[MPI process "<<my_rank<<"] Sending a job to "<<m_workers.front()<<": ID: "<<ID
                                        <<"; t = "<<t
                                        <<"; t+1 = "<<tpp
                                        <<"; dt = "<<dt0
                                        <<"; x(t) = "<<coords.at(0)
                                        <<"; y(t) = "<<coords.at(1)
                                        <<"; z(t) = "<<coords.at(2)
                                        <<"; x'(t) = "<<nextcoords.at(0)
                                        <<"; y'(t) = "<<nextcoords.at(1)
                                        <<"; z'(t) = "<<nextcoords.at(2)
                                        <<"\n";



        comm.isend(&ID, 1, m_workers.front(),&req,0);
        comm.isend(&t, 1, m_workers.front(),&req,1);
        comm.isend(&tpp, 1, m_workers.front(),&req,2);
        comm.isend(&dt0, 1, m_workers.front(),&req,3);
        comm.isend(coords.data(), coords.size(), m_workers.front(),&req,4);
        comm.isend(&angle, 1, m_workers.front(),&req,5);
        comm.isend(nextcoords.data(), nextcoords.size(), m_workers.front(),&req,6);
        comm.isend(&nextangle, 1, m_workers.front(),&req,7);
        comm.isend(&stop, 1, m_workers.front(),&req,8);
        m_workers.pop();
        njobs++;
      }
    }

    stop = true;
    for (index_t w = 1; w!=proc_count; w++)
      comm.isend(&stop, 1,w,&req,8);
  }
  else
  {
    while (true)
    {
      recv_tuple_t send;
      send_tuple_t recv;
      index_t Nintervals = 2;

      index_t ID;
      gsVector<real_t> coords(3), coordsprev(3);
      real_t t, tpp, dt0, dt;
      real_t x, y, z, xpp, ypp, zpp, xref, yref, zref;
      real_t angle, nextangle;

      comm.recv(&stop,        1, 0, 8, MPI_STATUS_IGNORE);
      gsInfo<<"[MPI process "<<my_rank<<"] stop = "<<stop<<"\n";

      if (stop)
      {
        gsInfo<<"[MPI process "<<my_rank<<"] I have to stop!!"<<"\n";
        break;
      }

      times.resize(Nintervals);
      solutions.resize(Nintervals);
      angles.resize(Nintervals);
      distances.resize(Nintervals+1);

      comm.recv(&ID,          1, 0, 0, MPI_STATUS_IGNORE);
      comm.recv(&t,           1, 0, 1, MPI_STATUS_IGNORE);
      comm.recv(&tpp,         1, 0, 2, MPI_STATUS_IGNORE);
      comm.recv(&dt0,         1, 0, 3, MPI_STATUS_IGNORE);
      comm.recv(coordsprev.data(),  coordsprev.size(), 0, 4, MPI_STATUS_IGNORE);
      comm.recv(&angle,       1, 0, 5, MPI_STATUS_IGNORE);
      comm.recv(coords.data()    ,  coords.size()    , 0, 6, MPI_STATUS_IGNORE);
      comm.recv(&nextangle,   1, 0, 7, MPI_STATUS_IGNORE);

      x = coordsprev.at(0);
      y = coordsprev.at(1);
      z = coordsprev.at(2);

      xref = coords.at(0);
      yref = coords.at(1);
      zref = coords.at(2);
      gsInfo<<"[MPI process "<<my_rank<<"] Received a job from "<<0<<": ID: "<<ID<<"; t = "<<t<<"; t+1 = "<<tpp<<"; dt = "<<dt<<"; x(t) = "<<x<<"; y(t) = "<<y<<"; z(t) = "<<z<<"; x'(t) = "<<xref<<"; y'(t) = "<<yref<<"; z'(t) = "<<zref<<"\n";

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
        solutions.at(k) = coords;
        angles.at(k) = nextangle;
        distances.at(k) = (coords-coordsprev).norm();

        coordsprev = coords;
      }
      coords<<xref,yref,zref;
      distances.at(Nintervals) = (coords-coordsprev).norm();

      real_t maxSleep = 1;
      real_t sleepTime = maxSleep * next_random();

      sleep(sleepTime);
      gsInfo<<"[MPI process "<<my_rank<<"] slept for "<<sleepTime<<" seconds\n";

      gsInfo<<"[MPI process "<<my_rank<<"] Sending a job to    "<<0<<": ID: "<<ID<<"; size = "<<solutions.size()<<"\n";
      gsInfo<<"\t\t"<<"time\tsolution\tangle\tdistance\n";
      for (index_t k=0; k!=solutions.size(); k++)
        gsInfo<<"\t\t"<<times.at(k)<<"\t"<<solutions.at(k).transpose()<<"\t"<<angles.at(k)<<"\t"<<distances.at(k)<<"\n";
      gsInfo<<"\t\t\t\t\t\t\t\t"<<distances.back()<<"\n";

      index_t size = solutions.size();

      // ID
      comm.send(&ID,1,0,0);
      // Number of objects to send
      comm.send(&size,1,0,1);
      for (index_t k=0; k!=solutions.size(); k++)
      {
        comm.send(&(times.at(k))    ,1,0,4*k+2+0);
        index_t solsize = solutions.at(k).size();
        comm.send(&solsize,1,0,4*k+2+1);
        comm.send(solutions.at(k).data(),solsize,0,4*k+2+2);
        comm.send(&(angles.at(k))   ,1,0,4*k+2+3);
        comm.send(&(distances.at(k)),1,0,4*k+2+4);
      }
      comm.send(&(distances.at(solutions.size())),1,0,4*solutions.size()+2+0);

    }
  }


  return 0;
}