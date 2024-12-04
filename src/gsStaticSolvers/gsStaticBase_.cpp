
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
#ifdef gsOptim_ENABLED
#include <gsOptim/gsOptim.h>
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
#ifdef gsOptim_ENABLED
		CLASS_TEMPLATE_INST gsStaticOpt<real_t, gsOptimBFGS<real_t>>;
		CLASS_TEMPLATE_INST gsStaticOpt<real_t, gsOptimLBFGS<real_t>>;
		CLASS_TEMPLATE_INST gsStaticOpt<real_t, gsOptimCG<real_t>>;
		CLASS_TEMPLATE_INST gsStaticOpt<real_t, gsOptimGD<real_t>>;
		// CLASS_TEMPLATE_INST gsStaticOpt<real_t, gsOptimNewton<real_t>>;
		CLASS_TEMPLATE_INST gsStaticOpt<real_t, gsOptimNM<real_t>>;
		CLASS_TEMPLATE_INST gsStaticOpt<real_t, gsOptimDE<real_t>>;
		CLASS_TEMPLATE_INST gsStaticOpt<real_t, gsOptimDEPRMM<real_t>>;
		CLASS_TEMPLATE_INST gsStaticOpt<real_t, gsOptimPSO<real_t>>;
		CLASS_TEMPLATE_INST gsStaticOpt<real_t, gsOptimPSODV<real_t>>;
		CLASS_TEMPLATE_INST gsStaticOpt<real_t, gsOptimSUMT<real_t>>;
#endif

}
