#include <gsCore/gsTemplateTools.h>

#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMBase.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMBase.hpp>

#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMLoadControl.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMLoadControl.hpp>

#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMRiks.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMRiks.hpp>

#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMCrisfield.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMCrisfield.hpp>

#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMConsistentCrisfield.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMConsistentCrisfield.hpp>

namespace gismo
{
	CLASS_TEMPLATE_INST gsALMBase<real_t>;
	CLASS_TEMPLATE_INST gsALMLoadControl<real_t>;
	CLASS_TEMPLATE_INST gsALMRiks<real_t>;
	CLASS_TEMPLATE_INST gsALMCrisfield<real_t>;
	CLASS_TEMPLATE_INST gsALMConsistentCrisfield<real_t>;
}
