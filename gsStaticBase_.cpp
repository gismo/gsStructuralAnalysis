
#include <gsCore/gsTemplateTools.h>

#include <gsStructuralAnalysis/gsStaticBase.h>

#include <gsStructuralAnalysis/gsStaticDR.h>
#include <gsStructuralAnalysis/gsStaticDR.hpp>

#include <gsStructuralAnalysis/gsStaticNewton.h>
#include <gsStructuralAnalysis/gsStaticNewton.hpp>

namespace gismo
{
		CLASS_TEMPLATE_INST gsStaticBase<real_t>;
		CLASS_TEMPLATE_INST gsStaticDR<real_t>;
		CLASS_TEMPLATE_INST gsStaticNewton<real_t>;
}
