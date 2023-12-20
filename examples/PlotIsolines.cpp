/** @file PlotIsolines.cpp

    @brief Testing file reading and writing

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <gismo.h>

using namespace gismo;

void writeToCSVfile(std::string name, gsMatrix<> matrix)
{
  std::ofstream file(name.c_str());
  for(int  i = 0; i < matrix.rows(); i++)
  {
    for(int j = 0; j < matrix.cols(); j++)
    {
       std::string str = std::to_string(matrix(i,j));
       if(j+1 == matrix.cols())
       {
           file<<std::setprecision(10)<<str;
       }
       else
       {
           file<<std::setprecision(10)<<str<<',';
       }
    }
    file<<'\n';
  }
}

int main(int argc, char *argv[])
{
  std::string input;
  std::string output = "data";
  index_t N = 1000;
  index_t patch = -1;

  std::vector<real_t> u;
  std::vector<index_t> dirs;

  gsCmdLine cmd("Hi, give me the path of your output");
  cmd.addString("i","in", "Input", input);
  cmd.addString("o","outDir", "Output directory", output);
  cmd.addInt("p","patch", "Patch", patch);
  cmd.addInt("N","npts", "number of points", N);

  cmd.addMultiReal("u","ucoord","u-coordinates of the isoline, all coordinates are scaled from 0 to 1",u);
  cmd.addMultiInt("d","dirs","directions of the isolines",dirs);

  try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

  gsMultiPatch<> mp;
  GISMO_ASSERT(!input.empty(),"Input cannot be empty!");
  GISMO_ASSERT(u.size()==dirs.size(),"Number of points and number of directions shouold be equal!");

  gsReadFile<>(input,mp);

  index_t ncols = (patch==-1) ? mp.nPatches() : 1;
  gsMatrix<> tmp;
  gsMatrix<> result(N,ncols*mp.geoDim());

  gsMatrix<> supp, coords(mp.parDim(),N);
  index_t dir;

  gsFileManager::mkdir(output);
  char sep = gsFileManager::getNativePathSeparator();
  for (index_t k = 0; k!=u.size(); k++)
  {
    dir = dirs[k];
    GISMO_ASSERT(u[k] <=1 && u[k] >= 0,"u coordinate must be in [0,1]");
    if (patch==-1)
    {
      for (size_t p = 0; p!=mp.nPatches(); p++)
      {
        supp = mp.patch(p).parameterRange();
        coords.row(dir).setConstant( u[k]*( supp(dir,1)-supp(dir,0) )+supp(dir,0) );
        coords.row(1-dir).setLinSpaced(N,0,1);
        coords.row(1-dir).array() *= ( supp(1-dir,1)-supp(1-dir,0) );
        coords.row(1-dir).array() += supp(1-dir,0);
        mp.patch(p).eval_into(coords,tmp);
        tmp.transposeInPlace();
        result.block(0,p*mp.geoDim(),tmp.rows(),tmp.cols()) = tmp;
      }
    }
    else
    {
      supp = mp.patch(patch).parameterRange();
      coords.row(dir).setConstant( u[k]*( supp(dir,1)-supp(dir,0) )+supp(dir,0) );
      coords.row(1-dir).setLinSpaced(N,0,1);
      coords.row(1-dir).array() *= ( supp(1-dir,1)-supp(1-dir,0) );
      coords.row(1-dir).array() += supp(1-dir,0);
      mp.patch(patch).eval_into(coords,result);
      result.transposeInPlace();
    }
    writeToCSVfile(output + sep + "line_u=" + std::to_string(u[k]) + "_dir=" + std::to_string(dir) + ".csv",result);
  }


  return 0;
}
