/*
 * Copyright (c) 2019-2020 Szabolcs Horvát.
 *
 * See the file LICENSE.txt for copying permission.
 */

#include "IGRNG.h"
#include "WolframCompileLibrary.h"
#include <random>


// We use a hard-coded value of 2^31 - 1 because:
//  - 'mint' is 32-bit on some platforms and this is th max supported value.
//  - 'unsigned long int', used by igraph, is 32-bit even on 64-bit Windows. The value must fit it in.
//  - We want reproducible results across platforms with a given random seed.
static const mint rng_Mma_max_int = 0x7fffffff;

static LibraryFunctionPointer randomInteger;
static LibraryFunctionPointer randomReal;
static LibraryFunctionPointer randomNormal;


void rng_Mma_get_function_pointers() {
    randomInteger = mma::libData->compileLibraryFunctions->getFunctionCallPointer("RandomInteger");
    randomReal    = mma::libData->compileLibraryFunctions->getFunctionCallPointer("RandomReal");
    randomNormal  = mma::libData->compileLibraryFunctions->getFunctionCallPointer("RandomNormal");
}


int igraph_rng_Mma_init(void **state) {
    IGRAPH_ERROR("Mathematica RNG error, unsupported function 'init' called.", IGRAPH_EINTERNAL);
    return IGRAPH_SUCCESS;
}

void igraph_rng_Mma_destroy(void *state) {
    // This function does not return an error code, so cannot use IGRAPH_ERROR.
    IGRAPH_FATAL("Mathematica RNG error, unsupported function 'destroy' called.");
}

int igraph_rng_Mma_seed(void *state, unsigned long int seed) {
    IGRAPH_ERROR("Mathematica RNG error, unsupported function 'seed' called.", IGRAPH_EINTERNAL);
    return IGRAPH_SUCCESS;
}

unsigned long int igraph_rng_Mma_get(void *state) {
    MArgument FPA[3];

    mint lo = 0, hi = rng_Mma_max_int;
    mint res;

    MArgument_getIntegerAddress(FPA[0]) = &lo;
    MArgument_getIntegerAddress(FPA[1]) = &hi;
    MArgument_getIntegerAddress(FPA[2]) = &res;

    int err = randomInteger(mma::libData, 2, FPA, FPA[2]);
    if (err)
        throw mma::LibraryError("RNG: Error calling RandomInteger", err);

    mma::libData->compileLibraryFunctions->WolframLibraryData_cleanUp(mma::libData, 1);

    return res;
}

igraph_real_t igraph_rng_Mma_get_real(void *state) {
    MArgument FPA[3];

    mreal lo = 0.0, hi = 1.0;
    mreal res;

    MArgument_getRealAddress(FPA[0]) = &lo;
    MArgument_getRealAddress(FPA[1]) = &hi;
    MArgument_getRealAddress(FPA[2]) = &res;

    int err = randomReal(mma::libData, 2, FPA, FPA[2]);
    if (err)
        throw mma::LibraryError("RNG: Error calling RandomReal.", err);

    mma::libData->compileLibraryFunctions->WolframLibraryData_cleanUp(mma::libData, 1);

    return res;
}

igraph_real_t igraph_rng_Mma_get_norm(void *state) {
    MArgument FPA[3];

    mreal mean = 0.0, stddev = 1.0;
    mreal res;

    MArgument_getRealAddress(FPA[0]) = &mean;
    MArgument_getRealAddress(FPA[1]) = &stddev;
    MArgument_getRealAddress(FPA[2]) = &res;

    int err = randomNormal(mma::libData, 2, FPA, FPA[2]);
    if (err)
        throw mma::LibraryError("RNG: Error calling RandomNormal.", err);

    mma::libData->compileLibraryFunctions->WolframLibraryData_cleanUp(mma::libData, 1);

    return res;
}


static const igraph_rng_type_t igraph_rngtype_Mma = {
    /* name= */      "Mathematica",
    /* min=  */      0,
    /* max=  */      rng_Mma_max_int,
    /* init= */      igraph_rng_Mma_init,
    /* destroy= */   igraph_rng_Mma_destroy,
    /* seed= */      igraph_rng_Mma_seed,
    /* get= */       igraph_rng_Mma_get,
    /* get_real= */  igraph_rng_Mma_get_real,
    /* get_norm= */  igraph_rng_Mma_get_norm,
    /* get_geom= */  0,
    /* get_binom= */ 0,
    /* get_exp= */   0,
    /* get_gamma= */ 0
};


static igraph_rng_t *rng_array[2];


/* Public functions */

void rngInit() {
    // This must run only once at startup (not on subsequent package loads).
    // We achieve this by using 'static'
    static igraph_rng_t rng_Default = *igraph_rng_default();
    static igraph_rng_t rng_Mma = {
        &igraph_rngtype_Mma,
        0,
        // prevent re-seeding after setting new generator; see RNG_BEGIN()
        // 0 = non-default
        // 1 = default, not yet seeded
        // 2 = default, already seeded
        /* def= */ 0
    };

    rng_array[0] = &rng_Mma;         // Mathematica's RNG
    rng_array[1] = &rng_Default;     // igraph's default RNG

    rng_Mma_get_function_pointers();

    // Seeding from the time is no longer needed because currenty there is no way to change the RNG method
    // without also seeding the generator in IGraph/M. See the definition of IGSeedRandom[].
    // igCheck(igraph_rng_seed(&rng_Default, std::random_device()()));

    rngSet(0); // set Mathematica's generator as default
}

void rngSet(mint id) {
    if (id < 0 || id >= sizeof(rng_array) / sizeof(igraph_rng_t *))
        throw mma::LibraryError("setRng: invalid random number generator ID");
    igraph_rng_set_default(rng_array[id]);
}

void rngSeed(mint s) {
    igCheck(igraph_rng_seed(igraph_rng_default(), s));
}

const char *rngName() {
    return igraph_rng_name(igraph_rng_default());
}
