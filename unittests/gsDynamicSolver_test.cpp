/** @file gsDynamicSolver_test.cpp

    @brief Regression recordings of gsDynamicNewmark<real_t,true> on a nonlinear 2-dof oscillator

    * Recording: integrates a fixed-dt trajectory and compares (t, |U|, |V|, |A|, Newton passes)
                 per step against reference arrays embedded in this file. Newton passes are
                 counted through the Jacobian op, since numIterations() carries no count.
                 Steps whose pass count changes under a scaling of the tolerances are marked
                 "boundary" and their pass count is not compared.
    * Stateful: the const step(t,dt,U,V,A) and the stateful step(dt) must agree bitwise.

    == GENERATOR ==
         Setting the environment variable GSSA_STAGE0_PRINT prints the reference arrays as
         C++ initialisers to gsInfo and skips the comparison against the embedded references.

    == BASIC REFERENCE ==
         - TEST(NAME_OF_TEST) { body_of_test }

    == CHECK MACRO REFERENCE ==
         - CHECK(EXPR);
         - CHECK_EQUAL(EXPECTED,ACTUAL);
         - CHECK_CLOSE(EXPECTED,ACTUAL,EPSILON);

    Author(s): H.M.Verhelst
 **/

#include "gismo_unittest.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicNewmark.h>
#include "gsStage0Oscillator.h"

SUITE(gsDynamicSolver_test)
{

namespace
{

const index_t N = 50;
const real_t  dt = real_t(0.02);
const real_t  baseTol = real_t(1e-10);

struct Recording
{
    std::vector<double> t, nU, nV, nA;
    std::vector<int>    passes;
    bool                allSuccess;
};

gsVector<real_t> initialU() { gsVector<real_t> U(2); U << real_t(0.5), real_t(1.0); return U; }

template <class Solver>
void setOptions(Solver & solver, const real_t alpha, const real_t tolScale)
{
    gsOptionList & o = solver.options();
    o.setReal("alpha", alpha);
    o.setReal("delta", real_t(0.5));
    o.setReal("TolF", baseTol*tolScale);
    o.setReal("TolU", baseTol*tolScale);
    o.setInt("MaxIter", 100);
    o.setString("Solver", "SimplicialLDLT");
    o.setSwitch("Quasi", false);
    solver.setTimeStep(dt);
}

// Integrates N steps with the stateful step(dt). Tolerances are scaled by \a tolScale.
// Newton passes per step are the Jacobian-call increase; the residual-call increase must be passes+1.
Recording runTrajectory(const real_t alpha, const real_t tolScale)
{
    gsStage0Oscillator<real_t> osc;
    gsDynamicNewmark<real_t,true> solver(osc.massOp(), osc.dampingOp(), osc.jacobianOp(), osc.residualOp());
    setOptions(solver, alpha, tolScale);

    const gsVector<real_t> U0 = initialU();
    const gsVector<real_t> V0 = gsVector<real_t>::Zero(2);
    solver.setU(U0);
    solver.setV(V0);
    solver.setA(osc.initialAcceleration(U0, V0, 0));

    Recording rec;
    rec.allSuccess = true;
    for (index_t k = 0; k < N; ++k)
    {
        osc.resetCounters();
        const gsStatus status = solver.step();
        CHECK(status == gsStatus::Success);
        if (status != gsStatus::Success)
        {
            rec.allSuccess = false;
            break;
        }
        CHECK_EQUAL(osc.jacobianCalls + 1, osc.residualCalls);
        rec.t.push_back(double(solver.time()));
        rec.nU.push_back(double(solver.solutionU().norm()));
        rec.nV.push_back(double(solver.solutionV().norm()));
        rec.nA.push_back(double(solver.solutionA().norm()));
        rec.passes.push_back(int(osc.jacobianCalls));
    }
    return rec;
}

// Boundary flag per step: pass count differs between any two of the runs.
std::vector<int> boundaryFlags(const Recording & a, const Recording & b, const Recording & c)
{
    std::vector<int> flag(a.passes.size(), 0);
    for (size_t k = 0; k < flag.size(); ++k)
        flag[k] = (a.passes[k] != b.passes[k]) || (a.passes[k] != c.passes[k]) || (b.passes[k] != c.passes[k]);
    return flag;
}

template <class V>
void printArray(const char * type, const char * name, const V & v, const char * fmt)
{
    gsInfo << "static const " << type << " " << name << "[] = {\n";
    char buf[64];
    for (size_t k = 0; k < v.size(); ++k)
    {
        std::snprintf(buf, sizeof(buf), fmt, v[k]);
        gsInfo << buf << (k+1 < v.size() ? "," : "") << ((k+1)%5==0 || k+1==v.size() ? "\n" : " ");
    }
    gsInfo << "};\n";
}

// Reference recordings of gsDynamicNewmark<double,true> on gsStage0Oscillator, nominal run
// (alpha 0.25, delta 0.5, TolF = TolU = 1e-10, SimplicialLDLT, dt 0.02, N 50 steps).
// Oscillator: m1 = m2 = 1, k = 1, k3 = 10, c = 0.05, F0 = 1, omega = 2; U0 = (0.5, 1), V0 = 0,
// A0 = initialAcceleration(U0, V0, 0). Boundary flags come from runs with both tolerances scaled by 0.7 and 1.4.
// Generated on 2026-10-03 (gsStructuralAnalysis b5ec853) by
//   GSSA_STAGE0_PRINT=1 OMP_NUM_THREADS=1 bin/unittests gsDynamicSolver_test
static const double refT[] = {
0.02, 0.040000000000000001, 0.059999999999999998, 0.080000000000000002, 0.10000000000000001,
0.12000000000000001, 0.14000000000000001, 0.16, 0.17999999999999999, 0.19999999999999998,
0.21999999999999997, 0.23999999999999996, 0.25999999999999995, 0.27999999999999997, 0.29999999999999999,
0.32000000000000001, 0.34000000000000002, 0.36000000000000004, 0.38000000000000006, 0.40000000000000008,
0.4200000000000001, 0.44000000000000011, 0.46000000000000013, 0.48000000000000015, 0.50000000000000011,
0.52000000000000013, 0.54000000000000015, 0.56000000000000016, 0.58000000000000018, 0.6000000000000002,
0.62000000000000022, 0.64000000000000024, 0.66000000000000025, 0.68000000000000027, 0.70000000000000029,
0.72000000000000031, 0.74000000000000032, 0.76000000000000034, 0.78000000000000036, 0.80000000000000038,
0.8200000000000004, 0.84000000000000041, 0.86000000000000043, 0.88000000000000045, 0.90000000000000047,
0.92000000000000048, 0.9400000000000005, 0.96000000000000052, 0.98000000000000054, 1.0000000000000004
};
static const double refNormU[] = {
1.1177248124958392, 1.116805382231566, 1.1152924843025553, 1.1132040302392541, 1.1105589400191227,
1.1073770119318247, 1.1036787828791961, 1.099485383158842, 1.0948183898007295, 1.0896996823601732,
1.0841513047391289, 1.0781953361399619, 1.0718537736867157, 1.0651484286158386, 1.0581008372786564,
1.0507321875463864, 1.0430632605954337, 1.035114387500607, 1.0269054195947158, 1.0184557111761512,
1.0097841128665699, 1.0009089737381685, 0.99184815023911344, 0.98261901993771561, 0.97323849816969299,
0.9637230557957277, 0.9540887364452717, 0.94435117182432049, 0.93452559388756906, 0.92462684290822172,
0.91466937071252741, 0.90466723857333009, 0.89463410947172328, 0.88458323463409949, 0.87452743443076186,
0.86447907388041578, 0.85445003314193591, 0.84445167349132555, 0.8344947993788856, 0.82458961724081681,
0.8147456918025322, 0.80497190065961277, 0.79527638795823208, 0.78566651802137788, 0.776148829781387,
0.76672899288381802, 0.75741176632280993, 0.74820096045366202, 0.73909940320399758, 0.7301089112697875
};
static const double refNormV[] = {
0.034553514612398292, 0.068157806625341405, 0.10070417561973762, 0.13209249051264163, 0.16223304848756984,
0.19104814161270156, 0.21847329057199125, 0.24445811927486294, 0.26896685857180952, 0.29197848116020797,
0.31348648240883353, 0.33349833279809737, 0.35203463666025864, 0.36912803875777012, 0.38482192494390799,
0.39916896580516242, 0.41222955296800212, 0.4240701768916938, 0.43476179271301413, 0.44437821730525645,
0.45299459639324169, 0.46068597553440321, 0.46752600321621102, 0.47358578839757393, 0.47893292869345938,
0.48363071922574302, 0.48773754610451353, 0.49130646273694428, 0.49438494186440579, 0.49701479157687867,
0.49923221970280357, 0.50106802804507033, 0.50254791600644588, 0.50369287224151948, 0.50451963304978065,
0.50504118719288915, 0.50526730853945789, 0.50520510024073073, 0.50485953683113005, 0.50423399353778187,
0.50333075499377011, 0.50215149832442363, 0.50069774808673395, 0.49897130269315537, 0.49697463367771083,
0.49471126042765901, 0.49218610379416156, 0.48940582231406754, 0.48637913463745192, 0.48311713118375166
};
static const double refNormA[] = {
1.705352734995889, 1.6550984755958096, 1.5996589414201623, 1.5395869110111549, 1.4755552434594557,
1.4083437244540742, 1.3388243860682882, 1.2679458047556618, 1.1967166126364515, 1.1261880333726946,
1.0574346548256166, 0.99153189662372387, 0.929527849249284, 0.87240668861167825, 0.82104134326768918,
0.77613537833138668, 0.73815876941327985, 0.70728869271969663, 0.68337172834974391, 0.66592347957574483,
0.65417286324261126, 0.64714406163122817, 0.64375703991816746, 0.64292424419582572, 0.64362691305220787,
0.64496428128063521, 0.64617732939984318, 0.64665314590469369, 0.64591683253438503, 0.64361676348204511,
0.63950728419725633, 0.63343136583473936, 0.62530457801292716, 0.61510100032848714, 0.60284126397130178,
0.58858269422535614, 0.57241143171614772, 0.554436391392786, 0.53478494186274017, 0.51360023790183684,
0.4910402105212423, 0.46727831429990691, 0.44250625906411695, 0.41693912490697177, 0.39082348851837284,
0.36444947417061224, 0.33816793205604484, 0.31241403438407861, 0.28773786955330061, 0.26483966587883945
};
static const int refPasses[] = {
3, 3, 3, 3, 3,
3, 3, 3, 3, 3,
3, 3, 3, 3, 3,
3, 3, 3, 3, 3,
3, 3, 3, 3, 3,
3, 3, 3, 3, 3,
3, 3, 3, 3, 3,
3, 3, 3, 3, 3,
3, 3, 3, 3, 3,
3, 3, 3, 3, 3
};
static const int refBoundary[] = {
0, 0, 0, 0, 0,
0, 0, 0, 0, 0,
0, 0, 0, 0, 0,
0, 0, 0, 0, 0,
0, 0, 0, 0, 0,
0, 0, 0, 0, 0,
0, 0, 0, 0, 0,
0, 0, 0, 0, 0,
0, 0, 0, 0, 0,
0, 0, 0, 0, 0
};

template <class T, size_t n>
size_t length(const T (&)[n]) { return n; }

} // anonymous namespace

TEST(newmark_nonlinear_oscillator_recording)
{
    UNITTEST_TIME_CONSTRAINT(1000);

    const bool generate = (std::getenv("GSSA_STAGE0_PRINT") != NULL);
    const real_t alpha = real_t(0.25);

    const Recording nominal = runTrajectory(alpha, real_t(1));
    CHECK(nominal.allSuccess);
    CHECK_EQUAL(N, index_t(nominal.passes.size()));
    if (nominal.passes.size() != size_t(N))
        return;

    std::vector<int> flags;
    if (generate)
    {
        const Recording lo = runTrajectory(alpha, real_t(0.7));
        const Recording hi = runTrajectory(alpha, real_t(1.4));
        CHECK(lo.allSuccess && hi.allSuccess);
        CHECK_EQUAL(N, index_t(lo.passes.size()));
        CHECK_EQUAL(N, index_t(hi.passes.size()));
        if (lo.passes.size() != size_t(N) || hi.passes.size() != size_t(N))
            return;
        flags = boundaryFlags(lo, nominal, hi);

        gsInfo << "// REFERENCES-GENERATED-BEGIN\n";
        printArray("double", "refT", nominal.t, "%.17g");
        printArray("double", "refNormU", nominal.nU, "%.17g");
        printArray("double", "refNormV", nominal.nV, "%.17g");
        printArray("double", "refNormA", nominal.nA, "%.17g");
        printArray("int", "refPasses", nominal.passes, "%d");
        printArray("int", "refBoundary", flags, "%d");
        gsInfo << "// REFERENCES-GENERATED-END\n";
    }
    else
    {
        CHECK_EQUAL(size_t(N), length(refT));
        CHECK_EQUAL(size_t(N), length(refNormU));
        CHECK_EQUAL(size_t(N), length(refNormV));
        CHECK_EQUAL(size_t(N), length(refNormA));
        CHECK_EQUAL(size_t(N), length(refPasses));
        CHECK_EQUAL(size_t(N), length(refBoundary));
        if (length(refT) != size_t(N) || length(refNormU) != size_t(N) || length(refNormV) != size_t(N) ||
            length(refNormA) != size_t(N) || length(refPasses) != size_t(N) || length(refBoundary) != size_t(N))
            return;

        flags.assign(refBoundary, refBoundary + N);
        for (index_t k = 0; k < N; ++k)
        {
            CHECK(std::abs(nominal.t[k]  - refT[k])     <= 1e-14*std::abs(refT[k]));
            CHECK(std::abs(nominal.nU[k] - refNormU[k]) <= 1e-8*std::abs(refNormU[k]));
            CHECK(std::abs(nominal.nV[k] - refNormV[k]) <= 1e-8*std::abs(refNormV[k]));
            CHECK(std::abs(nominal.nA[k] - refNormA[k]) <= 1e-8*std::abs(refNormA[k]));
            if (!refBoundary[k])
                CHECK_EQUAL(refPasses[k], nominal.passes[k]);
        }
    }

    index_t nBoundary = 0, nNonlinear = 0;
    for (index_t k = 0; k < N; ++k)
    {
        nBoundary  += flags[k];
        nNonlinear += (nominal.passes[k] >= 3);
    }
    CHECK(5*(N - nBoundary) >= 4*N);
    CHECK(5*nNonlinear >= 4*N);
}

TEST(newmark_const_step_matches_stateful_step_bitwise)
{
    UNITTEST_TIME_CONSTRAINT(1000);

    gsStage0Oscillator<real_t> oscA, oscB;
    gsDynamicNewmark<real_t,true> stateful(oscA.massOp(), oscA.dampingOp(), oscA.jacobianOp(), oscA.residualOp());
    gsDynamicNewmark<real_t,true> constant(oscB.massOp(), oscB.dampingOp(), oscB.jacobianOp(), oscB.residualOp());
    setOptions(stateful, real_t(0.25), real_t(1));
    setOptions(constant, real_t(0.25), real_t(1));

    const gsVector<real_t> U0 = initialU();
    const gsVector<real_t> V0 = gsVector<real_t>::Zero(2);
    const gsVector<real_t> A0 = oscA.initialAcceleration(U0, V0, 0);
    stateful.setU(U0);
    stateful.setV(V0);
    stateful.setA(A0);

    gsVector<real_t> U = U0, V = V0, A = A0;
    real_t t = 0;
    for (index_t k = 0; k < N; ++k)
    {
        const gsStatus s1 = stateful.step(dt);
        const gsStatus s2 = static_cast<const gsDynamicNewmark<real_t,true> &>(constant).step(t, dt, U, V, A);
        t += dt;
        CHECK(s1 == gsStatus::Success);
        CHECK(s2 == gsStatus::Success);
        if (s1 != gsStatus::Success || s2 != gsStatus::Success)
            break;
        bool equal = true;
        for (index_t i = 0; i < 2; ++i)
            equal = equal && (stateful.solutionU()(i) == U(i)) && (stateful.solutionV()(i) == V(i))
                          && (stateful.solutionA()(i) == A(i));
        CHECK(equal);
        if (!equal)
            break;
    }
    CHECK(stateful.time() == t);
}

}
