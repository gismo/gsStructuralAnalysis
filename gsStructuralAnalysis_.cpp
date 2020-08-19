
#include <gismo.h>

#include <gsCore/gsTemplateTools.h>

#include <gsStructuralAnalysis/gsArcLengthIterator.h>
#include <gsStructuralAnalysis/gsArcLengthIterator.hpp>

#include <gsStructuralAnalysis/gsTimeIntegrator.h>

#include <gsStructuralAnalysis/gsBucklingSolver.h>
#include <gsStructuralAnalysis/gsBucklingSolver.hpp>

#include <gsStructuralAnalysis/gsModalSolver.h>
#include <gsStructuralAnalysis/gsModalSolver.hpp>

namespace gismo
{
		CLASS_TEMPLATE_INST gsArcLengthIterator<real_t>;
		CLASS_TEMPLATE_INST gsTimeIntegrator<real_t>;
		CLASS_TEMPLATE_INST gsBucklingSolver<real_t>;
		CLASS_TEMPLATE_INST gsModalSolver<real_t>;
}
