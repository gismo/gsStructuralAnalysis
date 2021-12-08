#include <gismo.h>

#include <gsCore/gsTemplateTools.h>
#include <gsStructuralAnalysis/gsEigenProblemBase.h>
#include <gsStructuralAnalysis/gsEigenProblemBase.hpp>

#include <gsStructuralAnalysis/gsBucklingSolver.h>

#include <gsStructuralAnalysis/gsModalSolver.h>

namespace gismo
{
		CLASS_TEMPLATE_INST gsEigenProblemBase<real_t,Spectra::GEigsMode::Cholesky>;
		CLASS_TEMPLATE_INST gsEigenProblemBase<real_t,Spectra::GEigsMode::RegularInverse>;
		CLASS_TEMPLATE_INST gsEigenProblemBase<real_t,Spectra::GEigsMode::ShiftInvert>;
		CLASS_TEMPLATE_INST gsEigenProblemBase<real_t,Spectra::GEigsMode::Buckling>;
		CLASS_TEMPLATE_INST gsEigenProblemBase<real_t,Spectra::GEigsMode::Cayley>;

		CLASS_TEMPLATE_INST gsBucklingSolver<real_t,Spectra::GEigsMode::Cholesky>;
		CLASS_TEMPLATE_INST gsBucklingSolver<real_t,Spectra::GEigsMode::RegularInverse>;
		CLASS_TEMPLATE_INST gsBucklingSolver<real_t,Spectra::GEigsMode::ShiftInvert>;
		CLASS_TEMPLATE_INST gsBucklingSolver<real_t,Spectra::GEigsMode::Buckling>;
		CLASS_TEMPLATE_INST gsBucklingSolver<real_t,Spectra::GEigsMode::Cayley>;

		CLASS_TEMPLATE_INST gsModalSolver<real_t,Spectra::GEigsMode::Cholesky>;
		CLASS_TEMPLATE_INST gsModalSolver<real_t,Spectra::GEigsMode::RegularInverse>;
		CLASS_TEMPLATE_INST gsModalSolver<real_t,Spectra::GEigsMode::ShiftInvert>;
		CLASS_TEMPLATE_INST gsModalSolver<real_t,Spectra::GEigsMode::Buckling>;
		CLASS_TEMPLATE_INST gsModalSolver<real_t,Spectra::GEigsMode::Cayley>;
}
