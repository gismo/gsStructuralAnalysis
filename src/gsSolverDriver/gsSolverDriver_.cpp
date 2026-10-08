 /** @file gsSolverDriver_.cpp

    @brief Explicit instantiations of the solver-driver templates.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#include <gsCore/gsTemplateTools.h>

#include <gsStructuralAnalysis/src/gsSolverDriver/gsStateSpace.h>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsSolverDriverTypes.h>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsContinuationProblem.h>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsConvergenceTest.h>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsTangentSolver.h>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsTangentSolver.hpp>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsPathConstraint.h>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsPathConstraint.hpp>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsPathConstraintCrisfield.h>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsPathConstraintCrisfield.hpp>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsStepController.h>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsNewtonCorrector.h>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsContinuationDriver.h>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsContinuationDriver.hpp>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsStabilityAnalysis.h>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsStabilityAnalysis.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST  gsStateSpace<real_t>;
    CLASS_TEMPLATE_INST  gsVectorStateSpace<real_t>;
    CLASS_TEMPLATE_INST  gsFunctionalStateSpace<real_t>;
    STRUCT_TEMPLATE_INST gsStepHistory<real_t>;
    STRUCT_TEMPLATE_INST gsPathPoint<real_t>;
    STRUCT_TEMPLATE_INST gsStepContext<real_t>;
    STRUCT_TEMPLATE_INST gsContinuationProblem<real_t>;
    CLASS_TEMPLATE_INST gsConvergenceTest<real_t>;
    CLASS_TEMPLATE_INST gsRelResidualTest<real_t>;
    CLASS_TEMPLATE_INST gsAbsResidualTest<real_t>;
    CLASS_TEMPLATE_INST gsRelUpdateTest<real_t>;
    CLASS_TEMPLATE_INST gsPathConstraintTest<real_t>;
    CLASS_TEMPLATE_INST gsMinIterTest<real_t>;
    CLASS_TEMPLATE_INST gsMaxIterTest<real_t>;
    CLASS_TEMPLATE_INST gsDivergedTest<real_t>;
    CLASS_TEMPLATE_INST gsRoundoffFloorTest<real_t>;
    CLASS_TEMPLATE_INST gsAllTest<real_t>;
    CLASS_TEMPLATE_INST gsAnyTest<real_t>;
    CLASS_TEMPLATE_INST gsTangentSolver<real_t>;
    CLASS_TEMPLATE_INST gsPathConstraint<real_t>;
    CLASS_TEMPLATE_INST gsFixedConstraint<real_t>;
    CLASS_TEMPLATE_INST gsLoadConstraint<real_t>;
    CLASS_TEMPLATE_INST gsDisplacementConstraint<real_t>;
    CLASS_TEMPLATE_INST gsRiksConstraint<real_t>;
    CLASS_TEMPLATE_INST gsCrisfieldConstraint<real_t>;
    CLASS_TEMPLATE_INST gsConsistentCrisfieldConstraint<real_t>;
    CLASS_TEMPLATE_INST gsStepController<real_t>;
    CLASS_TEMPLATE_INST gsNewtonCorrector<real_t>;
    CLASS_TEMPLATE_INST gsContinuationDriver<real_t>;
    CLASS_TEMPLATE_INST gsStabilityAnalysis<real_t>;
}
