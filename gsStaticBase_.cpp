
#include <gsCore/gsTemplateTools.h>

#include <gsStructuralAnalysis/gsStaticBase.h>

#include <gsStructuralAnalysis/gsStaticDR.h>
#include <gsStructuralAnalysis/gsStaticDR.hpp>

#include <gsStructuralAnalysis/gsStaticNewton.h>
#include <gsStructuralAnalysis/gsStaticNewton.hpp>

#include <gsStructuralAnalysis/gsStaticComposite.h>
#include <gsStructuralAnalysis/gsStaticComposite.hpp>

namespace gismo
{
		CLASS_TEMPLATE_INST gsStaticBase<real_t>;
		CLASS_TEMPLATE_INST gsStaticDR<real_t>;
		CLASS_TEMPLATE_INST gsStaticNewton<real_t>;
		CLASS_TEMPLATE_INST gsStaticComposite<real_t>;
}
