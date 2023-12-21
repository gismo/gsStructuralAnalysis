#include <gsCore/gsTemplateTools.h>

#include <gsStructuralAnalysis/src/gsStructuralAnalysisTools/gsPanelCreator.h>

namespace gismo
{

#ifdef GISMO_WITH_PYBIND11

  namespace py = pybind11;

  void pybind11_init_gsStructuralAnalysis(py::module &m)
  {
    gismo::pybind11_init_gsPanelCreator( m );
    /*Bindings for gsStructuralAnalysis go here*/
  }

#endif
}

