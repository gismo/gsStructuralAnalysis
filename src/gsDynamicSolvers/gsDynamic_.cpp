#include <gsCore/gsTemplateTools.h>

#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicBase.h>
#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicBase.hpp>

#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicExplicitEuler.h>
#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicExplicitEuler.hpp>

#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicImplicitEuler.h>
#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicImplicitEuler.hpp>

// #include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicNewmark.h>
// #include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicNewmark.hpp>

// #include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicBathe.h>
// #include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicBathe.hpp>


namespace gismo
{
	CLASS_TEMPLATE_INST gsDynamicBase<real_t>;
	CLASS_TEMPLATE_INST gsDynamicExplicitEuler<real_t>;
	CLASS_TEMPLATE_INST gsDynamicImplicitEuler<real_t>;
	// CLASS_TEMPLATE_INST gsDynamicNewmark<real_t>;
	// CLASS_TEMPLATE_INST gsDynamicBathe<real_t>;

}
