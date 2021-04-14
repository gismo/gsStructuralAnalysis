#include <gismo.h>

#include <gsCore/gsTemplateTools.h>

#include <gsStructuralAnalysis/gsModalSolver.h>
#include <gsStructuralAnalysis/gsModalSolver.hpp>

namespace gismo
{
		CLASS_TEMPLATE_INST gsModalSolver<real_t,Spectra::GEigsMode::Cholesky>;
		CLASS_TEMPLATE_INST gsModalSolver<real_t,Spectra::GEigsMode::RegularInverse>;
		CLASS_TEMPLATE_INST gsModalSolver<real_t,Spectra::GEigsMode::ShiftInvert>;
		CLASS_TEMPLATE_INST gsModalSolver<real_t,Spectra::GEigsMode::Buckling>;
		CLASS_TEMPLATE_INST gsModalSolver<real_t,Spectra::GEigsMode::Cayley>;
}
