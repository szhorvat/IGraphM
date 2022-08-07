/*
 * Copyright (c) 2019-2022 Szabolcs Horvát.
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


igraph_error_t igraph_rng_Mma_init(void **state) {
    IGRAPH_ERROR("Mathematica RNG error, unsupported function 'init' called.", IGRAPH_EINTERNAL);
    return IGRAPH_SUCCESS;
}

void igraph_rng_Mma_destroy(void *state) {
    // This function does not return an error code, so cannot use IGRAPH_ERROR.
    IGRAPH_FATAL("Mathematica RNG error, unsupported function 'destroy' called.");
}

igraph_error_t igraph_rng_Mma_seed(void *state, igraph_uint_t seed) {
    IGRAPH_ERROR("Mathematica RNG error, unsupported function 'seed' called.", IGRAPH_EINTERNAL);
    return IGRAPH_SUCCESS;
}

igraph_uint_t igraph_rng_Mma_get(void *state) {
    MArgument FPA[3];

    // Use a full 32-bit range regardless of whether mint/igraph_integer_t
    // are 32-bit or 64-bit. This is to ensure the same random sequence on
    // all platforms. We assume that int is 32-bit.
    mint lo = INT_MIN, hi = INT_MAX;
    mint res;

    MArgument_getIntegerAddress(FPA[0]) = &lo;
    MArgument_getIntegerAddress(FPA[1]) = &hi;
    MArgument_getIntegerAddress(FPA[2]) = &res;

    int err = randomInteger(mma::libData, 2, FPA, FPA[2]);
    if (err)
        throw mma::LibraryError("RNG: Error calling RandomInteger", err);

    mma::libData->compileLibraryFunctions->WolframLibraryData_cleanUp(mma::libData, 1);

    // Transfer 32 bits stored in a signed type into an unsigned one.
    unsigned int u = (int) res;

    return res;
}

igraph_integer_t igraph_rng_Mma_get_int(void *state, igraph_integer_t l, igraph_integer_t h) {
    MArgument FPA[3];

    mint lo = l, hi = h;
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

    // Ensure sampling from half-open interval [0, 1).
    // It is not documented whether RandomReal[] may return an exact 1.0,
    // so we guard against this possibility.
    do {
        int err = randomReal(mma::libData, 2, FPA, FPA[2]);
        if (err)
            throw mma::LibraryError("RNG: Error calling RandomReal.", err);
    } while (res == 1.0);

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
    /* bits= */      32,
    /* init= */      igraph_rng_Mma_init,
    /* destroy= */   igraph_rng_Mma_destroy,
    /* seed= */      igraph_rng_Mma_seed,
    /* get= */       igraph_rng_Mma_get,
    /* get_int= */   igraph_rng_Mma_get_int,
    /* get_real= */  igraph_rng_Mma_get_real,
    /* get_norm= */  igraph_rng_Mma_get_norm,
    /* get_geom= */  nullptr,
    /* get_binom= */ nullptr,
    /* get_exp= */   nullptr,
    /* get_gamma= */ nullptr,
    /* get_pois= */  nullptr
};

/* Pointers to each available RNG */
static igraph_rng_t *rng_array[5];

/* Pointers to the RNG types for not-yet-initialized generators
 * or NULL for already initialized ones. */
static const igraph_rng_type_t *rng_type[5];


/* Public functions */

void rngInit() {
    // This must run only once at startup (not on subsequent package loads).
    // We achieve this by using 'static'
    static igraph_rng_t rng_Default = *igraph_rng_default();
    static igraph_rng_t rng_Mma = {
        &igraph_rngtype_Mma,
        nullptr,
        true
    };
    static igraph_rng_t rng_MT19937;
    static igraph_rng_t rng_PCG32;
    static igraph_rng_t rng_PCG64;

    rng_array[0] = &rng_Mma;     rng_type[0] = nullptr; // Mathematica's RNG
    rng_array[1] = &rng_Default; rng_type[1] = nullptr; // igraph's default RNG
    rng_array[2] = &rng_MT19937; rng_type[2] = &igraph_rngtype_mt19937;
    rng_array[3] = &rng_PCG32;   rng_type[3] = &igraph_rngtype_pcg32;
    rng_array[4] = &rng_PCG64;   rng_type[4] = &igraph_rngtype_pcg64;

    rng_Mma_get_function_pointers();

    // Seeding from the time is no longer needed because currenty there is no way to change the RNG method
    // without also seeding the generator in IGraph/M. See the definition of IGSeedRandom[].
    // igCheck(igraph_rng_seed(&rng_Default, std::random_device()()));

    rngSet(0); // set Mathematica's generator as default
}

void rngSet(mint id) {
    if (id < 0 || id >= sizeof(rng_array) / sizeof(rng_array[0]))
        throw mma::LibraryError("setRng: invalid random number generator ID");
    if (rng_type[id] != nullptr) {
        igCheck(igraph_rng_init(rng_array[id], rng_type[id])); // expected to fail for PCG64 on 32-bit systems
        rng_type[id] = nullptr; // mark as already initialized
    }
    igraph_rng_set_default(rng_array[id]);
}

void rngSeed(mint s) {
    igCheck(igraph_rng_seed(igraph_rng_default(), s));
}

const char *rngName() {
    return igraph_rng_name(igraph_rng_default());
}
