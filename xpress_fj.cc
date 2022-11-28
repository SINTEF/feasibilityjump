/*
 * XPress integration for the Feasibility Jump heuristic.
 */

extern "C"
{
#include "xprs.h"
}

#include <cstdio>
#include <string>
#include <chrono>
#include <vector>
#include <cassert>
#include <mutex>
#include <thread>
#include <functional>
#include <atomic>
#include <cmath>
#include <climits>
#include <chrono>

#include "feasibilityjump.hh"

const int NUM_THREADS = 2;

std::atomic_size_t totalNumSolutionsFound(0);
std::atomic_size_t totalNumSolutionsAdded(0);
std::atomic_bool presolveFinished(false);
std::atomic_bool heuristicFinished(false);
std::chrono::steady_clock::time_point startTime;

struct Solution
{
    std::vector<double> assignment;
    bool includesContinuous;
};

std::vector<Solution> heuristicSolutions;
std::mutex heuristicSolutions_mutex;

std::mutex presolvedProblem_mutex;
std::mutex nonPresolvedProblem_mutex;

struct ProblemInstance
{
    int numCols;
    std::vector<char> varTypes;
    std::vector<double> lb;
    std::vector<double> ub;
    std::vector<double> objCoeffs;

    int numRows;
    int numNonZeros;
    std::vector<char> rowtypes;
    std::vector<double> rhs;
    std::vector<double> rhsrange;
    std::vector<int> rowStart;
    std::vector<int> colIdxs;
    std::vector<double> colCoeffs;
};

struct FJData
{
    std::vector<int> originalIntegerCols;
    XPRSprob originalProblemCopy = nullptr;
    XPRSprob presolvedProblemCopy = nullptr;
    ProblemInstance originalData;
    ProblemInstance presolvedData;
};

FJData gFJData;

// A function that receives the result of any solution added with XPRSaddmipsol.
void XPRS_CC userSolNotify(XPRSprob problem, void *_data, const char *solName, int status)
{
    if (status == 0)
        printf(FJ_LOG_PREFIX "XPress received solution: An error occurred while processing the solution.\n");
    if (status == 1)
        printf(FJ_LOG_PREFIX "XPress received solution: Solution is feasible.\n");
    if (status == 2)
        printf(FJ_LOG_PREFIX "XPress received solution: Solution is feasible after reoptimizing with fixed globals.\n");
    if (status == 3)
        printf(FJ_LOG_PREFIX "XPress received solution: A local search heuristic was applied and a feasible solution discovered.\n");
    if (status == 4)
        printf(FJ_LOG_PREFIX "XPress received solution: A local search heuristic was applied but a feasible solution was not found.\n");
    if (status == 5)
        printf(FJ_LOG_PREFIX "XPress received solution: Solution is infeasible and a local search could not be applied.\n");
    if (status == 6)
        printf(FJ_LOG_PREFIX "XPress received solution: Solution is partial and a local search could not be applied.\n");
    if (status == 7)
        printf(FJ_LOG_PREFIX "XPress received solution: Failed to reoptimize the problem with globals fixed to the provided solution. Likely because a time or iteration limit was reached.\n");
    if (status == 8)
        printf(FJ_LOG_PREFIX "XPress received solution: Solution is dropped. This can happen if the MIP problem is changed or solved to completion before the solution could be processed.\n");
}

std::string inputFilename;
std::string outDir;

void XPRS_CC intsol(XPRSprob problem, void *data)
{
    double time = std::chrono::duration_cast<std::chrono::milliseconds>(
                      std::chrono::steady_clock::now() - startTime)
                      .count() /
                  1000.0;

    double objVal;
    XPRSgetdblattrib(problem, XPRS_MIPOBJVAL, &objVal);
    printf("XPRESS %g %g\n", time, objVal);

    if (outDir.size() > 0)
    {
        char path[512];
        sprintf(path, "%s/%.2f_%g_%s.sol", outDir.c_str(), time, objVal, inputFilename.c_str());

        printf("Writing to filename '%s'\n", path);
        XPRSwriteslxsol(problem, path, "");
    }
}

int XPRS_CC checktime(XPRSprob problem, void *_data)
{
    if (presolveFinished && totalNumSolutionsFound != totalNumSolutionsAdded)
    {
        std::lock_guard<std::mutex> guard(heuristicSolutions_mutex);

        for (auto &sol : heuristicSolutions)
        {

            std::vector<double> values;
            values.reserve(gFJData.originalIntegerCols.size());
            for (auto &idx : gFJData.originalIntegerCols)
                values.push_back(sol.assignment[idx]);
            assert(values.size() == gFJData.originalIntegerCols.size());

            // Add only the integer values.
            XPRSaddmipsol(problem, values.size(), values.data(),
                          gFJData.originalIntegerCols.data(), nullptr);

            // We could also have added all values and let the solver handle it, like this:
            //   auto data = new std::vector<double>(sol.assignment);
            //   XPRSaddmipsol(problem, data->size(), data->data(), nullptr, nullptr);
        }
        totalNumSolutionsAdded += heuristicSolutions.size();
        heuristicSolutions.clear();
    }
    if (heuristicFinished)
    {
        printf(FJ_LOG_PREFIX "all threads terminated. Removing MIP solver callback.\n");
        XPRSremovecbchecktime(problem, nullptr, nullptr);
    }

    return 0;
}

void XPRS_CC presolve_callback(XPRSprob problem, void *_data)
{
    presolveFinished = true;
    checktime(problem, nullptr);
}

ProblemInstance getXPRSProblemData(XPRSprob problem)
{
    ProblemInstance data;

    XPRSgetintattrib(problem, XPRS_COLS, &data.numCols);
    data.varTypes = std::vector<char>(data.numCols);
    XPRSgetcoltype(problem, data.varTypes.data(), 0, data.numCols - 1);
    data.lb = std::vector<double>(data.numCols);
    XPRSgetlb(problem, data.lb.data(), 0, data.numCols - 1);
    data.ub = std::vector<double>(data.numCols);
    XPRSgetub(problem, data.ub.data(), 0, data.numCols - 1);
    data.objCoeffs = std::vector<double>(data.numCols);
    XPRSgetobj(problem, data.objCoeffs.data(), 0, data.numCols - 1);

    XPRSgetintattrib(problem, XPRS_ROWS, &data.numRows);
    data.rowtypes = std::vector<char>(data.numRows);
    XPRSgetrowtype(problem, data.rowtypes.data(), 0, data.numRows - 1);
    data.rhs = std::vector<double>(data.numRows);
    XPRSgetrhs(problem, data.rhs.data(), 0, data.numRows - 1);
    data.rhsrange = std::vector<double>(data.numRows);
    XPRSgetrhsrange(problem, data.rhsrange.data(), 0, data.numRows - 1);

    XPRSgetrows(problem, nullptr, nullptr, nullptr, 0, &data.numNonZeros, 0, data.numRows - 1);
    printf(FJ_LOG_PREFIX "copying %d x %d matrix with %d nonzeros.\n",
           data.numCols, data.numRows, data.numNonZeros);
    data.rowStart = std::vector<int>(data.numRows + 1);
    data.colIdxs = std::vector<int>(data.numNonZeros);
    data.colCoeffs = std::vector<double>(data.numNonZeros);

    XPRSgetrows(problem,
                data.rowStart.data(),
                data.colIdxs.data(),
                data.colCoeffs.data(),
                data.numNonZeros,
                &data.numNonZeros,
                0,
                data.numRows - 1);

    return data;
}

bool copyDataToHeuristicSolver(FeasibilityJumpSolver &solver, ProblemInstance &data, int relaxContinuous)
{
    printf("initializing FJ with %d vars %d constraints\n", data.numCols, data.numRows);
    for (int colIdx = 0; colIdx < data.numCols; colIdx += 1)
    {
        VarType vartype = VarType::Continuous;
        if (data.varTypes[colIdx] == 'C')
        {
            vartype = VarType::Continuous;
        }
        else if (data.varTypes[colIdx] == 'I')
        {
            vartype = VarType::Integer;
        }
        else if (data.varTypes[colIdx] == 'B')
        {
            vartype = VarType::Integer;
        }
        else
        {
            printf(FJ_LOG_PREFIX "unsupported variable type '%c' (%d).\n",
                   data.varTypes[colIdx], data.varTypes[colIdx]);
            return false;
        }

        solver.addVar(vartype, data.lb[colIdx], data.ub[colIdx], data.objCoeffs[colIdx]);
    }

    for (int rowIdx = 0; rowIdx < data.numRows; rowIdx += 1)
    {
        RowType rowtype;
        if (data.rowtypes[rowIdx] == 'N')
        {
            continue;
        }
        else if (data.rowtypes[rowIdx] == 'L')
        {
            rowtype = RowType::Lte;
        }
        else if (data.rowtypes[rowIdx] == 'G')
        {
            rowtype = RowType::Gte;
        }
        else if (data.rowtypes[rowIdx] == 'E')
        {
            rowtype = RowType::Equal;
        }
        else if (data.rowtypes[rowIdx] == 'R')
        {
            // For the range constraint, we need two linear inequalities:
            // rhs - range <= lhs <= rhs

            if (data.rhsrange[rowIdx] < 0.0)
            {
                printf(FJ_LOG_PREFIX "unsupported negative range value '%g'.\n",
                       data.rhsrange[rowIdx]);
                return false;
            }

            solver.addConstraint(RowType::Gte,
                                 data.rhs[rowIdx] - data.rhsrange[rowIdx],
                                 data.rowStart[rowIdx + 1] - data.rowStart[rowIdx],
                                 &data.colIdxs[data.rowStart[rowIdx]],
                                 &data.colCoeffs[data.rowStart[rowIdx]],
                                 relaxContinuous);
            solver.addConstraint(RowType::Lte,
                                 data.rhs[rowIdx],
                                 data.rowStart[rowIdx + 1] - data.rowStart[rowIdx],
                                 &data.colIdxs[data.rowStart[rowIdx]],
                                 &data.colCoeffs[data.rowStart[rowIdx]],
                                 relaxContinuous);
            continue;
        }
        else
        {
            printf(FJ_LOG_PREFIX "unsupported constraint type '%c'. Ignoring constraint.\n", data.rowtypes[rowIdx]);
            return false;
        }

        solver.addConstraint(rowtype,
                             data.rhs[rowIdx],
                             data.rowStart[rowIdx + 1] - data.rowStart[rowIdx],
                             &data.colIdxs[data.rowStart[rowIdx]],
                             &data.colCoeffs[data.rowStart[rowIdx]],
                             relaxContinuous);
    }

    return true;
}

// An object containing a function to be executed when the object is destructed.
struct Defer
{
    std::function<void(void)> func;
    Defer(std::function<void(void)> pFunc) : func(pFunc){};
    ~Defer() { func(); }
};

void mapHeuristicSolution(FJStatus &status, bool usePresolved)
{

    Solution s;
    bool conversionOk = false;
    if (usePresolved)
    {
        assert(status.numVars == gFJData.presolvedData.numCols);

        printf(FJ_LOG_PREFIX "received a solution from presolved instance.\n");
        XPRSprob copy;
        XPRScreateprob(&copy);
        XPRScopyprob(copy, gFJData.presolvedProblemCopy, "");
        auto &data = gFJData.presolvedData;

        for (int i = 0; i < data.numCols; i += 1)
            if (data.varTypes[i] != 'C')
            {
                XPRSchgbounds(copy, 1, &i, "B", &status.solution[i]);
            }

        XPRSsetintcontrol(copy, XPRS_LPITERLIMIT, INT_MAX - 2);
        XPRSmipoptimize(copy, "");
        int status;
        XPRSgetintattrib(copy, XPRS_MIPSTATUS, &status);
        switch (status)
        {
        case XPRS_MIP_OPTIMAL:
            s.assignment.resize(gFJData.originalData.numCols);
            s.includesContinuous = true;

            XPRSgetmipsol(copy, s.assignment.data(), nullptr);
            conversionOk = true;
            break;
        case XPRS_MIP_LP_NOT_OPTIMAL:
            printf(FJ_LOG_PREFIX "Unexpected status: global search incomplete.\n");
            break;
        default:
            printf(FJ_LOG_PREFIX "Unexpected solution status (%i).\n", status);
        }
    }
    else
    {
        printf(FJ_LOG_PREFIX "received a solution from non-presolved instance.\n");
        s.assignment = std::vector<double>(status.solution, status.solution + status.numVars);
        s.includesContinuous = true;
        conversionOk = true;
    }

    if (conversionOk)
    {
        {
            std::lock_guard<std::mutex> guard(heuristicSolutions_mutex);
            heuristicSolutions.push_back(s);
            totalNumSolutionsFound += 1;
        }

        XPRSprob copy;
        XPRScreateprob(&copy);
        XPRScopyprob(copy, gFJData.originalProblemCopy, "");
        XPRSaddcbintsol(copy, intsol, nullptr, 0);

        for (int i = 0; i < gFJData.originalData.numCols; i += 1)
            if (gFJData.originalData.varTypes[i] != 'C')
                XPRSchgbounds(copy, 1, &i, "B", &s.assignment[i]);
        XPRSmipoptimize(copy, "");
        int status;
        XPRSgetintattrib(copy, XPRS_MIPSTATUS, &status);
        switch (status)
        {
        case XPRS_MIP_OPTIMAL:
            printf(FJ_LOG_PREFIX "MIP extension successful.\n");
            break;
        case XPRS_MIP_LP_NOT_OPTIMAL:
            printf(FJ_LOG_PREFIX "checktime Unexpected status: global search incomplete.\n");
            break;
        default:
            printf(FJ_LOG_PREFIX "checktime Unexpected solution status (%i).\n", status);
        }
    }
}

const int maxEffort = 100000000;

// Starts background threads running the Feasibility Jump heuristic.
// Also installs the check-time callback to report any feasible solutions
// back to the MIP solver.
void start_feasibility_jump_heuristic(XPRSprob problem, size_t maxTotalSolutions, bool heuristicOnly, bool relaxContinuous = false, bool exponentialDecay = false, int verbose = 0)
{

    // Copy the problem to the heuristic.
    XPRScreateprob(&gFJData.originalProblemCopy);
    XPRScopyprob(gFJData.originalProblemCopy, problem, "");

    {
        auto allThreadsTerminated = std::make_shared<Defer>([]()
                                                            { heuristicFinished = true; });

        for (int thread_idx = 0; thread_idx < NUM_THREADS; thread_idx += 1)
        {
            auto seed = thread_idx;
            bool usePresolved = thread_idx % 2 == 1;
            double decayFactor = (!exponentialDecay) ? 1.0 : 0.9999;

            std::thread(
                [verbose, maxTotalSolutions, usePresolved, seed,
                 relaxContinuous, decayFactor, allThreadsTerminated]()
                {
                    // Prepare data for the non-presolved version.
                    {
                        std::lock_guard<std::mutex> guard(nonPresolvedProblem_mutex);
                        if (gFJData.originalData.numCols == 0)
                        {
                            gFJData.originalData = getXPRSProblemData(gFJData.originalProblemCopy);
                        }
                    }

                    // Produce the presolved solution
                    if (usePresolved)
                    {
                        std::lock_guard<std::mutex> guard(presolvedProblem_mutex);
                        if (gFJData.presolvedProblemCopy == nullptr)
                        {
                            XPRScreateprob(&gFJData.presolvedProblemCopy);
                            XPRSsetlogfile(gFJData.presolvedProblemCopy, "presolve.log");
                            XPRScopyprob(gFJData.presolvedProblemCopy, gFJData.originalProblemCopy, "");
                            XPRSsetintcontrol(gFJData.presolvedProblemCopy, XPRS_LPITERLIMIT, 0);

                            XPRSmipoptimize(gFJData.presolvedProblemCopy, "");
                            gFJData.presolvedData = getXPRSProblemData(gFJData.presolvedProblemCopy);
                        }
                    }

                    ProblemInstance &data = usePresolved ? gFJData.presolvedData : gFJData.originalData;
                    FeasibilityJumpSolver solver(seed, verbose, decayFactor);
                    bool copyOk = copyDataToHeuristicSolver(solver, data, relaxContinuous);
                    if (!copyOk)
                        return;

                    solver.solve(
                        nullptr, [maxTotalSolutions, usePresolved](FJStatus status) -> CallbackControlFlow
                        {
                            
    
                            double time = std::chrono::duration_cast<std::chrono::milliseconds>(
                                std::chrono::steady_clock::now() - startTime).count() /1000.0;

                            // If we received a solution, put it on the queue.
                            if (status.solution != nullptr)
                            {
                                printf("FJSOL %g %g\n", time, status.solutionObjectiveValue);
    
                                mapHeuristicSolution(status, usePresolved);
                            }
    
                            // If we have enough solutions or spent enough time, quit.
                            auto quitNumSol = totalNumSolutionsFound >= maxTotalSolutions;
                            if(quitNumSol) printf(FJ_LOG_PREFIX "quitting because number of solutions %zd >= %zd.\n", totalNumSolutionsFound.load(), maxTotalSolutions);
                            auto quitEffort = status.effortSinceLastImprovement > maxEffort;
                            if(quitEffort) printf(FJ_LOG_PREFIX "quitting because effort %d > %d.\n", status.effortSinceLastImprovement , maxEffort);
                            
                            auto quit = quitNumSol || quitEffort || heuristicFinished;
                            if (quit)
                                printf(FJ_LOG_PREFIX "effort rate: %g Mops/sec\n", status.totalEffort / time / 1.0e6);
                            return quit ? CallbackControlFlow::Terminate : CallbackControlFlow::Continue; });
                })
                .detach();
        }
    }

    if (heuristicOnly)
    {
        while (!heuristicFinished)
        {
            std::this_thread::sleep_for(std::chrono::milliseconds(50));
        }
        printf(FJ_LOG_PREFIX "all threads exited.\n");
    }
}

#define CHECK_RETURN(call)                                  \
    do                                                      \
    {                                                       \
        int result_ = call;                                 \
        if (result_ != 0)                                   \
        {                                                   \
            fprintf(stderr, "Line %d: %s failed with %d\n", \
                    __LINE__, #call, result_);              \
            returnCode = result_;                           \
            goto cleanup;                                   \
        }                                                   \
    } while (0)

int printUsage()
{
    printf("Usage: xpress_fj [--save-solutions|-s OUTDIR] [--verbose|-v] [--heuristic-only|-h] [--exponential-decay|-e] [--relax-continuous|-r] INFILE\n");
    return 1;
}

int main(int argc, char *argv[])
{
    int verbose = 0;
    bool heuristicOnly = false;
    bool relaxContinuous = false;
    bool exponentialDecay = false;

    std::string inputPath;
    for (int i = 1; i < argc; i += 1)
    {
        std::string argvi(argv[i]);
        if (argvi == "--save-solutions" || argvi == "-s")
        {
            if (i + 1 < argc)
                outDir = std::string(argv[i + 1]);
            else
                return printUsage();
            i += 1;
        }
        else if (argvi == "--verbose" || argvi == "-v")
            verbose += 1;
        else if (argvi == "--heuristic-only" || argvi == "-h")
            heuristicOnly = true;
        else if (argvi == "--relax-continuous" || argvi == "-r")
            relaxContinuous = true;
        else if (argvi == "--exponential-decay" || argvi == "-e")
            exponentialDecay = true;
        else if (!inputPath.empty())
            return printUsage();
        else
            inputPath = argvi;
    }

    if (inputPath.empty())
        return printUsage();

    inputFilename = inputPath.substr(inputPath.find_last_of("/\\") + 1);

    int returnCode = 0;
    XPRSprob problem = nullptr;
    if (XPRSinit("") != 0)
    {
        char message[512];
        XPRSgetlicerrmsg(message, sizeof(message));
        fprintf(stderr, "Licensing error: %s\n", message);
        return 1;
    }

    startTime = std::chrono::steady_clock::now();

    CHECK_RETURN(XPRScreateprob(&problem));
    CHECK_RETURN(XPRSsetlogfile(problem, "xpress.log"));
    CHECK_RETURN(XPRSreadprob(problem, inputPath.c_str(), ""));

    // Install the solution callback to report when
    // a solution was found by the solver.
    CHECK_RETURN(XPRSaddcbintsol(problem, intsol, nullptr, 0));
    CHECK_RETURN(XPRSaddcbusersolnotify(problem, userSolNotify, nullptr, 0));

    start_feasibility_jump_heuristic(problem, 1, heuristicOnly, relaxContinuous, exponentialDecay, verbose);

    if (!heuristicOnly)
    {

        int numColsAll;
        int numColsOrig;
        XPRSgetintattrib(problem, XPRS_COLS, &numColsAll);
        XPRSgetintattrib(problem, XPRS_ORIGINALCOLS, &numColsOrig);
        assert(numColsAll == numColsOrig);

        // Prepare the list of integer variables to be used in
        // `checktime` for adding mip solutions.
        auto varTypes = std::vector<char>(numColsAll);
        XPRSgetcoltype(problem, varTypes.data(), 0, numColsAll - 1);
        for (int i = 0; i < numColsAll; i += 1)
            if (varTypes[i] != 'C')
                gFJData.originalIntegerCols.push_back(i);

        // An error in XPress causes solutions to be rejected sometimes
        // if they are added using `addmipsol` while presolve is running.
        // So, we need a callback to detect that presolve has finished.
        // We wait until this has happened before adding any heuristic solutions
        // to XPress.
        CHECK_RETURN(XPRSaddcbpresolve(problem, presolve_callback, nullptr, 0));

        // Install a callback that converts the heuristic solutions
        // into XPRSaddmipsol calls...
        CHECK_RETURN(XPRSaddcbchecktime(problem, checktime, nullptr, 0));
        // ...and then solve normally.
        XPRSsetintcontrol(problem, XPRS_THREADS, 1);

        double time = std::chrono::duration_cast<std::chrono::milliseconds>(
                          std::chrono::steady_clock::now() - startTime)
                          .count() /
                      1000.0;

        int timeout = std::ceil(60.0 - time);
        XPRSsetintcontrol(problem, XPRS_MAXTIME, -timeout);
        CHECK_RETURN(XPRSmipoptimize(problem, ""));

        double ub, lb;
        CHECK_RETURN(XPRSgetdblattrib(problem, XPRS_MIPBESTOBJVAL, &ub));
        CHECK_RETURN(XPRSgetdblattrib(problem, XPRS_BESTBOUND, &lb));

        time = std::chrono::duration_cast<std::chrono::milliseconds>(
                   std::chrono::steady_clock::now() - startTime)
                   .count() /
               1000.0;

        printf("{\"exit_time\":%g, \"lb\": %g, \"ub\": %g}\n", time, lb, ub);

        heuristicFinished = true;

        int status;
        XPRSgetintattrib(problem, XPRS_MIPSTATUS, &status);
        switch (status)
        {
        case XPRS_MIP_OPTIMAL:
            printf("XPRESS Solved, optimal.\n");

            break;
        default:
            printf(FJ_LOG_PREFIX "Unexpected solution status (%i).\n", status);
        }
    }

cleanup:
    if (returnCode > 0)
    {
        /* There was an error with the solver. Get the error code and error message.
         * If prob is still NULL then the error was in XPRScreateprob() and
         * we cannot find more detailed error information.
         */
        if (problem != NULL)
        {
            int errorCode = -1;
            char errorMessage[512] = {0};
            XPRSgetintattrib(problem, XPRS_ERRORCODE, &errorCode);
            XPRSgetlasterror(problem, errorMessage);
            fprintf(stderr, "Error %d: %s\n", errorCode, errorMessage);
        }
    }

    // Normally, we would clean up XPress at this point, if it should be used
    // as part of a larger program.  However, for this benchmarking program, we
    // avoid having to deal with waiting for the threads to shut down if we
    // just skip de-allocating the global Xpress structs.  If we do XPRSfree
    // here, we risk that a thread running the FJ heuristic will run some
    // XPress function between the XPRSfree call and the program shutting down,
    // which is a use-after-free error.
    //
    // XPRSdestroyprob(problem);
    // XPRSfree();

    return returnCode;
}
