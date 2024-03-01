#include <gsCore/gsTemplateTools.h>

#ifdef gsXbraid_ENABLED
#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicXBraid.h>
#endif

namespace gismo
{
#ifdef gsXbraid_ENABLED
	CLASS_TEMPLATE_INST gsDynamicXBraid<real_t>;
#endif
}
