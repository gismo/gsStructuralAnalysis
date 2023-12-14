
#include <gsCore/gsTemplateTools.h>

#include <gsStructuralAnalysis/src/gsStaticSolvers/gsStaticBase.h>

#include <gsStructuralAnalysis/src/gsStaticSolvers/gsStaticDR.h>
#include <gsStructuralAnalysis/src/gsStaticSolvers/gsStaticDR.hpp>

#include <gsStructuralAnalysis/src/gsStaticSolvers/gsStaticNewton.h>
#include <gsStructuralAnalysis/src/gsStaticSolvers/gsStaticNewton.hpp>

#include <gsStructuralAnalysis/src/gsStaticSolvers/gsStaticComposite.h>
#include <gsStructuralAnalysis/src/gsStaticSolvers/gsStaticComposite.hpp>

namespace gismo
{
		CLASS_TEMPLATE_INST gsStaticBase<real_t>;
		CLASS_TEMPLATE_INST gsStaticDR<real_t>;
		CLASS_TEMPLATE_INST gsStaticNewton<real_t>;
		CLASS_TEMPLATE_INST gsStaticComposite<real_t>;
}
