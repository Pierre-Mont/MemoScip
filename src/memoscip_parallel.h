#pragma once
/** Option -parallel: rejects a learned cut nearly parallel to an earlier kept learned cut.
 *
 *  The cuts are compared in a canonical space where the cosine does not depend on how a cut is written:
 *  - every slack is replaced by its exact definition (model slacks slack_<i> from their equality row, learned slacks
 *    sc_<k> from their cut), so only the original variables remain;
 *  - each domain block (equality row sum_v x_v = 1 over binaries with coefficient 1) is centred: the all-ones
 *    direction of the block is projected out, so the result does not depend on which value is eliminated.
 *  A new cut is rejected if |cos| with a kept learned cut is > 1 - 1e-4. The model rows are not compared. */
#include <string>
#include <vector>
#include "scip/scip.h"

/** Returns true if the cut consvals . consvars >= rhs must be rejected. If it is accepted, its canonical vector is kept
 *  until parallelRegister is called with the same name. */
bool parallelCheck(SCIP* scip, const std::string& cutname, const std::vector<SCIP_VAR*>& consvars,
                   const std::vector<SCIP_Real>& consvals);

/** Called when the cut is added with its slack: consvals . consvars = rhs, where consvars contains the new slack.
 *  Records the definition of the slack and keeps the cut for the next comparisons. */
void parallelRegister(SCIP* scip, const std::string& cutname, SCIP_VAR* slack,
                      const std::vector<SCIP_VAR*>& consvars, const std::vector<SCIP_Real>& consvals);

/** Prints the statistics of the check. */
void parallelPrintStatistics();
