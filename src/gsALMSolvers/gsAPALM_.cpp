#include <gsCore/gsTemplateTools.h>

#include <gsStructuralAnalysis/src/gsALMSolvers/gsAPALMData.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsAPALMData.hpp>

#include <gsStructuralAnalysis/src/gsALMSolvers/gsAPALMDataContainer.h>

#include <gsStructuralAnalysis/src/gsALMSolvers/gsAPALM.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsAPALM.hpp>

// #include <gsStructuralAnalysis/src/gsALMSolvers/gsAPALMApp.h>
// #include <gsStructuralAnalysis/src/gsALMSolvers/gsAPALMApp.hpp>

namespace gismo
{
	CLASS_TEMPLATE_INST gsAPALMData<real_t,std::pair<real_t,real_t>>;
	CLASS_TEMPLATE_INST gsAPALMData<real_t,std::pair<gsVector<real_t>,real_t>>;
	CLASS_TEMPLATE_INST gsAPALMData<real_t,std::pair<gsMatrix<real_t>,real_t>>;

	CLASS_TEMPLATE_INST gsAPALMDataContainer<real_t,std::pair<real_t,real_t>>;
	CLASS_TEMPLATE_INST gsAPALMDataContainer<real_t,std::pair<gsVector<real_t>,real_t>>;
	CLASS_TEMPLATE_INST gsAPALMDataContainer<real_t,std::pair<gsMatrix<real_t>,real_t>>;

	CLASS_TEMPLATE_INST gsAPALM<real_t>;
}
