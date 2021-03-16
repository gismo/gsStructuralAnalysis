#include <gismo.h>

#include <gsCore/gsTemplateTools.h>

#include <gsStructuralAnalysis/gsTimeIntegrator.h>
#include <gsStructuralAnalysis/gsTimeIntegrator.hpp>

namespace gismo
{
		CLASS_TEMPLATE_INST gsTimeIntegrator<real_t>;
}
