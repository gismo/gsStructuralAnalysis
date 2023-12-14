
#include <gsCore/gsTemplateTools.h>
#include <gsStructuralAnalysis/src/gsEigenSolvers/gsEigenProblemBase.h>
#include <gsStructuralAnalysis/src/gsEigenSolvers/gsEigenProblemBase.hpp>

#include <gsStructuralAnalysis/src/gsEigenSolvers/gsBucklingSolver.h>

#include <gsStructuralAnalysis/src/gsEigenSolvers/gsModalSolver.h>

namespace gismo
{
		CLASS_TEMPLATE_INST gsEigenProblemBase<real_t>;

		CLASS_TEMPLATE_INST gsBucklingSolver<real_t>;

		CLASS_TEMPLATE_INST gsModalSolver<real_t>;
}
