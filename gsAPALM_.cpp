#include <gsCore/gsTemplateTools.h>

#include <gsStructuralAnalysis/gsAPALMData.h>
#include <gsStructuralAnalysis/gsAPALMData.hpp>

#include <gsStructuralAnalysis/gsAPALM.h>
#include <gsStructuralAnalysis/gsAPALM.hpp>

namespace gismo
{
	CLASS_TEMPLATE_INST gsAPALMData<real_t>;
	CLASS_TEMPLATE_INST gsAPALM<real_t>;
}
