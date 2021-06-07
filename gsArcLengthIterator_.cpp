#include <gismo.h>

#include <gsCore/gsTemplateTools.h>

#include <gsStructuralAnalysis/gsALMBase.h>
#include <gsStructuralAnalysis/gsALMBase.hpp>

#include <gsStructuralAnalysis/gsArcLengthIterator.h>
#include <gsStructuralAnalysis/gsArcLengthIterator.hpp>

#include <gsStructuralAnalysis/gsALMLoadControl.h>
#include <gsStructuralAnalysis/gsALMLoadControl.hpp>

#include <gsStructuralAnalysis/gsALMRiks.h>
#include <gsStructuralAnalysis/gsALMRiks.hpp>

#include <gsStructuralAnalysis/gsALMCrisfield.h>
#include <gsStructuralAnalysis/gsALMCrisfield.hpp>

#include <gsStructuralAnalysis/gsALMConsistentCrisfield.h>
#include <gsStructuralAnalysis/gsALMConsistentCrisfield.hpp>

#include <gsStructuralAnalysis/gsALMExplicitIterations.h>
#include <gsStructuralAnalysis/gsALMExplicitIterations.hpp>

namespace gismo
{
		CLASS_TEMPLATE_INST gsArcLengthIterator<real_t>;

		CLASS_TEMPLATE_INST gsALMBase<real_t>;
		CLASS_TEMPLATE_INST gsALMLoadControl<real_t>;
		CLASS_TEMPLATE_INST gsALMRiks<real_t>;
		CLASS_TEMPLATE_INST gsALMCrisfield<real_t>;
		CLASS_TEMPLATE_INST gsALMConsistentCrisfield<real_t>;
		CLASS_TEMPLATE_INST gsALMExplicitIterations<real_t>;
}
