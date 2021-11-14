#include <gismo.h>

#include <gsCore/gsTemplateTools.h>

#include <gsStructuralAnalysis/gsSpaceTimeHierarchy.h>
// #include <gsStructuralAnalysis/gsSpaceTimeHierarchy.hpp>

namespace gismo
{
		CLASS_TEMPLATE_INST gsSpaceTimeHierarchy<real_t,std::pair<gsVector<real_t>,real_t>>;
}
