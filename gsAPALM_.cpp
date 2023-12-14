#include <gsCore/gsTemplateTools.h>

#include <gsStructuralAnalysis/gsAPALMData.h>
#include <gsStructuralAnalysis/gsAPALMData.hpp>

#include <gsStructuralAnalysis/gsAPALMDataContainer.h>

#include <gsStructuralAnalysis/gsAPALM.h>
#include <gsStructuralAnalysis/gsAPALM.hpp>

// #include <gsStructuralAnalysis/gsAPALMApp.h>
// #include <gsStructuralAnalysis/gsAPALMApp.hpp>

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
