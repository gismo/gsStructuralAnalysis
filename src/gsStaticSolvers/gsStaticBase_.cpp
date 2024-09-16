
#include <gsCore/gsTemplateTools.h>

#include <gsStructuralAnalysis/src/gsStaticSolvers/gsStaticBase.h>

#include <gsStructuralAnalysis/src/gsStaticSolvers/gsStaticDR.h>
#include <gsStructuralAnalysis/src/gsStaticSolvers/gsStaticDR.hpp>

#include <gsStructuralAnalysis/src/gsStaticSolvers/gsStaticNewton.h>
#include <gsStructuralAnalysis/src/gsStaticSolvers/gsStaticNewton.hpp>

#include <gsStructuralAnalysis/src/gsStaticSolvers/gsStaticOpt.h>
#include <gsStructuralAnalysis/src/gsStaticSolvers/gsStaticOpt.hpp>

#include <gsOptimizer/gsGradientDescent.h>
#ifdef gsHLBFGS_ENABLED
#include <gsHLBFGS/gsHLBFGS.h>
#endif

#include <gsStructuralAnalysis/src/gsStaticSolvers/gsStaticComposite.h>
#include <gsStructuralAnalysis/src/gsStaticSolvers/gsStaticComposite.hpp>

namespace gismo
{
		CLASS_TEMPLATE_INST gsStaticBase<real_t>;
		CLASS_TEMPLATE_INST gsStaticDR<real_t>;
		CLASS_TEMPLATE_INST gsStaticNewton<real_t>;
		CLASS_TEMPLATE_INST gsStaticComposite<real_t>;

		CLASS_TEMPLATE_INST gsOptProblemStatic<real_t>;
		CLASS_TEMPLATE_INST gsStaticOpt<real_t, gsGradientDescent<real_t>>;
#ifdef gsHLBFGS_ENABLED
		CLASS_TEMPLATE_INST gsStaticOpt<real_t, gsHLBFGS<real_t>>;
#endif

}
