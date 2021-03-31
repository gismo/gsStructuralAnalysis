#include <gismo.h>

#include <gsCore/gsTemplateTools.h>

#include <gsStructuralAnalysis/gsDynamicRelaxationLC.h>
#include <gsStructuralAnalysis/gsDynamicRelaxationALM.h>

namespace gismo
{
	CLASS_TEMPLATE_INST gsDynamicRelaxationLC<real_t>;
	CLASS_TEMPLATE_INST gsDynamicRelaxationALM<real_t>;
}
