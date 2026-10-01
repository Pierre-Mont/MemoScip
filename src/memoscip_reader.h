#ifndef MEMOSCIP_READER_H
#define MEMOSCIP_READER_H

#include <string>
#include "scip/scip.h"

// Read a KPCG instance from `path` and add variables/constraints/objective
// to the given `scip` instance. On error, returns false and `err` contains
// a short message.
bool readKPCGInstance(SCIP* scip, const std::string& path, int encoding, std::string& err);

#endif // MEMOSCIP_READER_H
