#include <gsCore/gsTemplateTools.h>

#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicBase.h>
#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicBase.hpp>

#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicExplicitEuler.h>
#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicExplicitEuler.hpp>

#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicImplicitEuler.h>
#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicImplicitEuler.hpp>

#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicRK4.h>
#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicRK4.hpp>

// #include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicCentralDifference.h>
// #include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicCentralDifference.hpp>

#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicNewmark.h>
#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicNewmark.hpp>

#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicWilson.h>
#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicWilson.hpp>

#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicBathe.h>
#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicBathe.hpp>

namespace gismo
{

	enum struct gsStatus;
	template <class T> struct gsStructuralAnalysisOps;

	CLASS_TEMPLATE_INST gsDynamicBase<real_t>;
	
	CLASS_TEMPLATE_INST gsDynamicExplicitEuler<real_t,false>;
	CLASS_TEMPLATE_INST gsDynamicExplicitEuler<real_t,true>;

	CLASS_TEMPLATE_INST gsDynamicImplicitEuler<real_t,false>;
	CLASS_TEMPLATE_INST gsDynamicImplicitEuler<real_t,true>;

	CLASS_TEMPLATE_INST gsDynamicRK4<real_t,false>;
	CLASS_TEMPLATE_INST gsDynamicRK4<real_t,true>;

	// CLASS_TEMPLATE_INST gsDynamicCentralDifference<real_t,false>;
	// CLASS_TEMPLATE_INST gsDynamicCentralDifference<real_t,true>;

	CLASS_TEMPLATE_INST gsDynamicNewmark<real_t,false>;
	CLASS_TEMPLATE_INST gsDynamicNewmark<real_t,true>;

    TEMPLATE_INST gsStatus
    gsDynamicNewmark<real_t,true>::_step_impl<true>(const real_t t, const real_t dt, gsVector<real_t> & U, gsVector<real_t> & V, gsVector<real_t> & A) const;

	CLASS_TEMPLATE_INST gsDynamicWilson<real_t,false>;
	CLASS_TEMPLATE_INST gsDynamicWilson<real_t,true>;

	CLASS_TEMPLATE_INST gsDynamicBathe<real_t,false>;
	CLASS_TEMPLATE_INST gsDynamicBathe<real_t,true>;

}
