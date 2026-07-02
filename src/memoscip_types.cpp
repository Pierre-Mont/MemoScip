#include "memoscip_types.h"

#include "scip/scip_lp.h"

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
    for (size_t i = 0; i < dualsol.size(); ++i)
    {
        if(SCIPisDualfeasZero(scip, dualsol[i]))
            dualsol[i] = 0.0;
    }
    for (size_t i = 0; i < redcost.size(); ++i)
    {
        if(SCIPisDualfeasZero(scip, redcost[i]))
            redcost[i] = 0.0;
    }
    return SCIP_OKAY;
}
