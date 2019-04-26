/*
 * Copyright (c) 2019 Szabolcs Horvát.
 *
 * See the file LICENSE.txt for copying permission.
 */

#ifndef IG_RNG_H
#define IG_RNG_H

#include "IGCommon.h"

void rngInit();         // Initialize RNGs for IGraph/M. Called once from IGlobal init().
const char *rngName();  // Get RNG name. Used to determine how to seed from WL. The WL generator can't be seeded from C.
void rngSeed(mint s);   // Seed the RNG. Works only for igraph's RNGs, will error for WL RNG.
void rngSet(mint id);   // Choose which RNG to use. See Method option of IGSeedRandom[].

#endif // IG_RNG_H
