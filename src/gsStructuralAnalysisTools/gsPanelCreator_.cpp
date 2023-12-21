#include <gsCore/gsTemplateTools.h>

#include <gsStructuralAnalysis/src/gsStructuralAnalysisTools/gsPanelCreator.h>
#include <gsStructuralAnalysis/src/gsStructuralAnalysisTools/gsPanelCreator.hpp>

namespace gismo
{

STRUCT_TEMPLATE_INST gsPanelCreator<real_t>;

#ifdef GISMO_WITH_PYBIND11

namespace py = pybind11;

void pybind11_init_gsPanelCreator(py::module &m)
{
  	using Class = gsPanelCreator<real_t>;
	py::class_<Class>(m, "gsPanelCreator")
        .def(py::init<>())
        .def_static("Plate", &Class::Plate,
    				"Creates a plate of length Lp and with Wp",
	    			py::arg("Lp"),py::arg("Wp"),py::arg("x")=1, py::arg("y")=0, py::arg("z")=0)
        .def_static("Strip", &Class::Strip,
    				"Creates a strip of length Lb and height Hw",
	    			py::arg("Lb"),py::arg("Hw"),py::arg("x")=1, py::arg("y")=0, py::arg("z")=0)
        .def_static("IBeam", &Class::IBeam,
    				"Creates an I-beam of length Lb, with web height Hw and with flange width Wf",
	    			py::arg("Lb"),py::arg("Hw"),py::arg("Wf"),py::arg("x")=1, py::arg("y")=0, py::arg("z")=0)
        .def_static("TBeam", &Class::TBeam,
    				"Creates a T-beam of length Lb, with web height Hw and with flange width Wf",
	    			py::arg("Lb"),py::arg("Hw"),py::arg("Wf"),py::arg("x")=1, py::arg("y")=0, py::arg("z")=0)
        .def_static("LBeam", &Class::LBeam,
    				"Creates an L-beam of length Lb, with web height Hw and with flange width Wf",
	    			py::arg("Lb"),py::arg("Hw"),py::arg("Wf"),py::arg("x")=1, py::arg("y")=0, py::arg("z")=0)
        .def_static("PanelT", &Class::PanelT,
    				"Creates a panel of length Lp and width Wp with a T-stiffener with web height Hw and flange width Wf",
	    			py::arg("Lp"),py::arg("Wp"),py::arg("Hw"),py::arg("Wf"),py::arg("x")=1, py::arg("y")=0, py::arg("z")=0)
        .def_static("PanelL", &Class::PanelL,
    				"Creates a panel of length Lp and width Wp with an L-stiffener with web height Hw and flange width Wf",
	    			py::arg("Lp"),py::arg("Wp"),py::arg("Hw"),py::arg("Wf"),py::arg("x")=1, py::arg("y")=0, py::arg("z")=0)
        .def_static("PlateGirderL", &Class::PlateGirderL,
    				"Creates a panel of length Lp and witdt Wp with a T-girder with web height Hwg and flange width Wfg and with an L-stiffener web height Hws and flange width Wfs",
	    			py::arg("Lp"),py::arg("Wp"),py::arg("Hwg"),py::arg("Wfg"),py::arg("Hws"),py::arg("Wfs"),py::arg("x")=1, py::arg("y")=0, py::arg("z")=0)

    ;
}
#endif
}
