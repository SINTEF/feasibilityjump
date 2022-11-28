/*
 * Solve MPS using XPress Optimizer.
 * Should be the same as the command-line XPress.
 * We are testing the interface to compare with the
 * heuristic-enhanced solver.
 */

extern "C"
{
#include "xprs.h"
}

#include <cstdio>
#include <string>
#include <chrono>
#include <vector>
#include <cmath>

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

std::string mpsFileName;
std::string outDir;
std::chrono::steady_clock::time_point startTime;

void XPRS_CC intsol(XPRSprob problem, void *data)
{
    double time = std::chrono::duration_cast<std::chrono::milliseconds>(
                      std::chrono::steady_clock::now() - startTime)
                      .count() /
                  1000.0;

    double objVal;
    XPRSgetdblattrib(problem, XPRS_MIPOBJVAL, &objVal);
    printf("Got MIP solution with objective %g\n", objVal);

    char path[512];
    sprintf(path, "%s/%.2f_%g_%s.sol", outDir.c_str(), time, objVal, mpsFileName.c_str());

    printf("Writing to filename '%s'\n", path);
    XPRSwriteslxsol(problem, path, "");
}

std::vector<int> originalIntegerCols;

int printUsage()
{
    printf("Usage: xpress_baseline [--timeout|-t TIMEOUT] [--save-solutions|-s OUTDIR] INFILE\n");
    return 1;
}


int main(int argc, char *argv[])
{

    int timeout = INT32_MAX/2;
    double ub, lb;
    double time;

    std::string mpsFile;
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
        else if (argvi == "--timeout" || argvi == "-t")
        {
            if (i + 1 < argc)
                timeout = std::stoi(argv[i + 1]);
            else
                return printUsage();
            i += 1;
        }
        else if (!mpsFile.empty())
            return printUsage();
        else
            mpsFile = argvi;
    }

    mpsFileName = mpsFile.substr(mpsFile.find_last_of("/\\") + 1);

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

    int numColsAll;
    std::vector<char> varTypes;

    CHECK_RETURN(XPRScreateprob(&problem));
    CHECK_RETURN(XPRSsetlogfile(problem, "xpress.log"));
    CHECK_RETURN(XPRSreadprob(problem, mpsFile.c_str(), ""));

    // The following operations are done only for fairness of comparison with
    // `xpress_fj`, the integration of the Feasibilty Jump heuristic with
    // XPress. Because `xpress_fj` uses the external interface of XPress for
    // integration, we cannot avoid making a copy of the problem data before
    // launching the heuristic in a background thread. This would probably not
    // be neccessary in an integration using the internal interface of XPress.
    //
    XPRSprob problemClone;
    XPRScreateprob(&problemClone);
    XPRScopyprob(problemClone, problem, "");
    
    XPRSgetintattrib(problem, XPRS_COLS, &numColsAll);
    varTypes.resize(numColsAll);
    XPRSgetcoltype(problem, varTypes.data(), 0, numColsAll - 1);
    for (int i = 0; i < numColsAll; i += 1)
        if (varTypes[i] != 'C')
            originalIntegerCols.push_back(i);


    // Install the solution callback to report when
    // a solution was found by the solver.
    XPRSaddcbintsol(problem, intsol, nullptr, 0);
    XPRSsetintcontrol(problem, XPRS_THREADS, 1);

    /* Search for an integer solution */

    time = std::chrono::duration_cast<std::chrono::milliseconds>(
               std::chrono::steady_clock::now() - startTime)
               .count() /
           1000.0;

    timeout = std::ceil(60.0 - time);
    XPRSsetintcontrol(problem, XPRS_MAXTIME, -timeout);
    CHECK_RETURN(XPRSmipoptimize(problem, ""));

    CHECK_RETURN(XPRSgetdblattrib(problem, XPRS_MIPBESTOBJVAL, &ub));
    CHECK_RETURN(XPRSgetdblattrib(problem, XPRS_BESTBOUND, &lb));

    time = std::chrono::duration_cast<std::chrono::milliseconds>(
               std::chrono::steady_clock::now() - startTime)
               .count() /
           1000.0;

    printf("{\"exit_time\":%g, \"lb\": %g, \"ub\": %g}\n", time, lb, ub);

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
    XPRSdestroyprob(problem);
    XPRSfree();

    return returnCode;
}
