#include "memoscip_types.h"

#include "scip/pub_lp.h"
#include "scip/scip_lp.h"

#include <cassert>
#include <cmath>
#include <iostream>
void SolCollector::init(SCIP* scip)
{
    int ncols = SCIPgetNLPCols(scip);
    int nrows = SCIPgetNLPRows(scip);
    dualsol.resize(nrows);
    redcost.resize(ncols);
}

SCIP_RETCODE SolCollector::collect(SCIP* scip)
{
    SCIP_LPI* lpi = nullptr;
    SCIP_CALL(SCIPgetLPI(scip, &lpi));
    SCIP_CALL(SCIPlpiGetSol(lpi, &objval,
        NULL, dualsol.data(),
        NULL, redcost.data()));
    
    // Only remove floating point noise: a small dual times a large row coefficient is not negligible
    for (size_t i = 0; i < dualsol.size(); ++i)
    {
        if(fabs(dualsol[i]) < 1e-12)
            dualsol[i] = 0.0;
    }
    // Reduced costs recomputed from the duals above (d = c - y^T A), so that the proofs and the reduced costs match exactly
    SCIP_COL** cols = SCIPgetLPCols(scip);
    for (size_t i = 0; i < redcost.size(); ++i)
    {
        SCIP_COL* col = cols[i];
        SCIP_ROW** colrows = SCIPcolGetRows(col);
        SCIP_Real* colvals = SCIPcolGetVals(col);
        SCIP_Real d = SCIPcolGetObj(col);
        for (int k = 0; k < SCIPcolGetNLPNonz(col); ++k)
        {
            int lppos = SCIProwGetLPPos(colrows[k]);
            if(lppos >= 0)
                d -= dualsol[lppos] * colvals[k];
        }
        redcost[i] = fabs(d) < 1e-12 ? 0.0 : d;
    }
    return SCIP_OKAY;
}
