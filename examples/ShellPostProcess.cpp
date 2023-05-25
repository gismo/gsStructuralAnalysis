/** @file matplotlib_example.cpp

    @brief Testing file reading and writing

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <iostream>
#include <fstream>
#include <filesystem>

#include <gismo.h>

using namespace gismo;

int exists(std::string fn)
{
  std::ifstream ifile;
  ifile.open(fn);
  if(ifile)
  {
     return 1;
  }
  else
  {
     return 0;
  }
}

int main(int argc, char *argv[])
{


  std::string path;
  std::string fn;

  std::string geoname = "geometry";
  std::string solname = "solution";

  std::string output = "solution";

  bool paraview = false;
  bool write = false;
  bool mesh = false;

  index_t N = 1000;

  //fn = "basis_thbs_01.xml"; //default example

  gsCmdLine cmd("Hi, give me the path of your output");
  cmd.addString("g","geo", "Geometry name", geoname);
  cmd.addString("s","sol", "Solution base name", solname);
  cmd.addString("o","out", "Output", output);

  cmd.addSwitch("paraview", "write to paraview", paraview);
  cmd.addSwitch("write", "write to files", write);
  cmd.addSwitch("mesh", "Plot mesh", mesh);

  cmd.addInt("N","pts", "number of points", N);

  cmd.addPlainString("filename", "Path to output; should at least contain geometry.xml and solution.xml", path);

  try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

  if (path.empty() )
  {
      gsInfo<< cmd.getMessage();
      gsInfo<<"\nType "<< argv[0]<< " -h, to get the list of command line options.\n";
      return 0;
  }

  std::vector<std::string> filenames;


  gsInfo<<"Initializing...\n";

  gsMultiPatch<> mp, sol;
  fn = geoname + ".xml";
  gsDebugVar(exists(path + fn));
  GISMO_ENSURE(exists(path + fn),"File "<<path+fn<<" does not exist!");
  gsReadFile<>(path + fn,mp);

  gsInfo<<"\tRead file "<<path + fn<<"\n";

  if (exists(path + solname + ".xml"))
  {
    fn = solname + ".xml";
    gsInfo<<"\tFound file "<<path + fn<<"\n";
    filenames.push_back(fn);
  }

  bool exist = true;
  index_t k = 0;
  fn = solname + std::to_string(0) + ".xml";
  exist = exists(path+fn);
  while (exist)
  {
    gsInfo<<"\tFound file "<<path + fn<<"\n";
    filenames.push_back(fn);

    k++;
    fn = solname + std::to_string(k) + ".xml";
    exist = exists(path+fn);
  }

  gsInfo<<"Done\n";

  gsInfo<<"Writing...\n";

  gsParaviewCollection collection(path + output);

  k = 0;

  if (paraview || write)
  {
    for (std::vector<std::string>::iterator it = filenames.begin(); it!=filenames.end(); it++, k++ )
    {
      gsInfo<<"\t"<<*it<<std::flush;

      gsReadFile<>(path + *it,sol);
      gsField<> solField(mp,sol);

      if (paraview)
      {
        std::string fileName = path + output + util::to_string(k);
        gsWriteParaview<>(solField, fileName, N, mesh);
        fileName = output + util::to_string(k) + "0";
        collection.addPart(fileName + ".vts",k);
        if (mesh) collection.addPart(fileName + "_mesh.vtp",k);
        gsInfo<<"---> "<<fileName<<".vts";
      }
      if (write)
      {
        std::ofstream file(path + output + util::to_string(k) + ".csv");

        for (size_t p = 0; p!=sol.nPatches(); p++)
        {
          gsMatrix<> ab = sol.patch(p).support();
          gsVector<> a = ab.col(0);
          gsVector<> b = ab.col(1);

          gsVector<unsigned> np = uniformSampleCount(a, b, N);
          gsMatrix<> pts = gsPointGrid(a, b, np);

          gsMatrix<> eval_geo   = mp.patch(p).eval(pts);//pts
          gsMatrix<> eval_field = solField.isParametric() ? solField.function(p).eval(pts) : solField.function(p).eval(eval_geo);

          eval_geo.transposeInPlace();
          eval_field.transposeInPlace();

          gsMatrix<> res(eval_geo.rows(),eval_geo.cols() + eval_field.cols());
          res << eval_geo,eval_field;

          for(index_t  i = 0; i < res.rows(); i++)
          {
            for(index_t j = 0; j < res.cols(); j++)
            {
              std::string str = std::to_string(res(i,j));
              if(j+1 == res.cols())
              {
                file<<str;
              }
              else
              {
                file<<str<<',';
              }
            }
            file<<'\n';
          }

        }
        gsInfo<<"---> "<<path + output + util::to_string(k) + ".csv";

      }
      gsInfo<<"\n";
    }

    gsInfo<<"Done\n";


    if (paraview)
        collection.save();
  }


// #ifdef GISMO_WITH_MATPLOTLIB


// #endif



  return 0;
}
