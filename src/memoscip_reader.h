#ifndef MEMOSCIP_READER_H
#define MEMOSCIP_READER_H

#include <string>
#include "scip/scip.h"

// Read a KPCG instance from `path` and add variables/constraints/objective
// to the given `scip` instance. On error, returns false and `err` contains
// a short message.
bool readKPCGInstance(SCIP* scip, const std::string& path, int encoding, std::string& err);

// Read a binary CSP instance in XCSP 2.1 format (extensional constraints, conflicts or supports)
// from `path`. One binary variable x<i>_<v> is created per value v of each variable i.
// encoding 0: each forbidden pair (a, b) of a constraint on (i, j) gives x<i>_<a> + x<j>_<b> <= 1.
// otherwise: one binary variable per pair of values of each constraint, forbidden pairs fixed to 0.
bool readXCSPInstance(SCIP* scip, const std::string& path, int encoding, std::string& err);

#endif // MEMOSCIP_READER_H
