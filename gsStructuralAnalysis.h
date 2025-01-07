#include <gsCore/gsTemplateTools.h>

namespace gismo
{

#ifdef GISMO_WITH_PYBIND11

  namespace py = pybind11;

  void pybind11_init_gsStructuralAnalysis(py::module &m);

#endif
}

