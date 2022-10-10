/*
** This file causes the entry points of my .C routines to be preloaded.
**
** Added on 24/09/2017 at the request of R CMD check.
**
** It adds one more layer of protection by declaring the number of arguments,
** and perhaps a tiny bit of speed.
*/
#include <R.h>
// #include <Rinternals.h>   // Causes some conflict of its macro LENGTH(x) with
                             // another macro of the same name which is present
                             // in one of the C++ headers perhaps included
                             // by some of *.h files. This issue was discovered
                             // with the bayesSurv package
                             // But it seems that Rinternals.h are not needed
                             // for smoothSurv compilation.
#include <R_ext/Rdynload.h>

#include "smoothSurvReg84.h"

static const R_CMethodDef Centries[] = {
    {"C_smoothSurvReg84", (DL_FUNC) &smoothSurvReg84, 43},    // smoothSurvReg.fit, minPenalty
    {NULL, NULL, 0}    
};

extern "C" void R_init_smoothSurv(DllInfo *dll){
    R_registerRoutines(dll, Centries, NULL,  NULL,     NULL);
    /*                      .C        .Call  .Fortran  .External  */
    
    /* The following line makes only those routines defined above
       available to outside packages, i.e., internal C++ things
       are now invisible.
    */
    R_useDynamicSymbols(dll, FALSE);
 
    /*
    ** This line makes them only be available via the symbols above
    **  i.e., .C("bayesBisurvreg", ) won't work but .C(C_bayesBisurvreg, )  will
    */
    R_forceSymbols(dll, TRUE);
}
    
