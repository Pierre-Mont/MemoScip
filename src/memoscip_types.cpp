#include "memoscip_types.h"

#include "scip/scip_lp.h"

void SolCollector::init(SCIP* scip)
{
    int ncols = SCIPgetNLPCols(scip);
    int nrows = SCIPgetNLPRows(scip);
    primsol.resize(ncols);
    dualsol.resize(nrows);
    activity.resize(nrows);
    redcost.resize(ncols);
}

SCIP_RETCODE SolCollector::collect(SCIP* scip)
{
    SCIP_LPI* lpi = nullptr;
    SCIP_CALL(SCIPgetLPI(scip, &lpi));
    SCIP_CALL(SCIPlpiGetSol(lpi, &objval,
        primsol.data(), dualsol.data(),
        activity.data(), redcost.data()));
    return SCIP_OKAY;
}
