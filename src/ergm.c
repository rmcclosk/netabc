#include <R.h>
#include <Rinternals.h>
#include <Rembedded.h>
#include "ergm.h"

void R_library(const char *libname);

void start_R(void)
{
    char *argv[] = {"R", "--vanilla", "--silent"};
    int argc = sizeof(argv)/sizeof(argv[0]);
    Rf_initEmbeddedR(argc, argv);
}

void stop_R(void)
{
    Rf_endEmbeddedR(0);
}

void create_ergm(int nnode)
{
    SEXP call, n, directed;
    SEXP net;
    SEXP formula, model;

    R_library("network");
    R_library("ergm");

    // create an empty network
    PROTECT(n = ScalarInteger(nnode));
    PROTECT(directed = ScalarLogical(0));
    PROTECT(call = lang3(install("network.initialize"), n, directed));
    PROTECT(net = eval(call, R_GlobalEnv));
    defineVar(install("net"), net, R_GlobalEnv);

    // create the model structure
    PROTECT(formula = mkString("net ~ edges + triangles"));
    PROTECT(call = lang2(install("as.formula"), formula));
    PROTECT(formula = eval(call, R_GlobalEnv));
    PROTECT(call = lang2(install("ergm"), formula));
    PROTECT(model = eval(call, R_GlobalEnv));
    
    PROTECT(call = lang2(install("print"), model));
    eval(call, R_GlobalEnv);
}

void R_library(const char *libname)
{
    SEXP arg = PROTECT(mkString(libname));
    SEXP call = PROTECT(lang2(install("library"), arg));
    SEXP call2 = PROTECT(lang2(install("suppressPackageStartupMessages"), call));
    eval(call2, R_GlobalEnv);
    UNPROTECT(3);
}

