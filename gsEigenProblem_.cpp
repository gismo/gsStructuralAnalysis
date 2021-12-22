#include <gismo.h>

#include <gsCore/gsTemplateTools.h>
#include <gsStructuralAnalysis/gsEigenProblemBase.h>
#include <gsStructuralAnalysis/gsEigenProblemBase.hpp>

#include <gsStructuralAnalysis/gsBucklingSolver.h>

#include <gsStructuralAnalysis/gsModalSolver.h>

namespace gismo
{
		CLASS_TEMPLATE_INST gsEigenProblemBase<real_t>;

		CLASS_TEMPLATE_INST gsBucklingSolver<real_t>;

		CLASS_TEMPLATE_INST gsModalSolver<real_t>;
}
