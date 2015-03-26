#include "byLocation.h"
#include "byAlignment.h"
#include <Rcpp.h>


// annotateByLocation
//SEXP annotateByLocation(SEXP inputData, SEXP dfGENETXData, SEXP dfGENETXGroupData, SEXP dfEXONData, SEXP stackEXONData, SEXP dfORGDBData, SEXP vecORGDBColumnsData, SEXP vecSEPData, SEXP vecRangesData, SEXP listFunctionParamData);
RcppExport SEXP ProbeAnnotator_annotateByLocation(SEXP inputDataSEXP, SEXP dfGENETXDataSEXP, SEXP dfGENETXGroupDataSEXP, SEXP dfEXONDataSEXP, SEXP stackEXONDataSEXP, SEXP dfORGDBDataSEXP, SEXP vecORGDBColumnsDataSEXP, SEXP vecSEPDataSEXP, SEXP vecRangesDataSEXP, SEXP listFunctionParamDataSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< SEXP >::type inputData(inputDataSEXP );
        Rcpp::traits::input_parameter< SEXP >::type dfGENETXData(dfGENETXDataSEXP );
        Rcpp::traits::input_parameter< SEXP >::type dfGENETXGroupData(dfGENETXGroupDataSEXP );
        Rcpp::traits::input_parameter< SEXP >::type dfEXONData(dfEXONDataSEXP );
        Rcpp::traits::input_parameter< SEXP >::type stackEXONData(stackEXONDataSEXP );
        Rcpp::traits::input_parameter< SEXP >::type dfORGDBData(dfORGDBDataSEXP );
        Rcpp::traits::input_parameter< SEXP >::type vecORGDBColumnsData(vecORGDBColumnsDataSEXP );
        Rcpp::traits::input_parameter< SEXP >::type vecSEPData(vecSEPDataSEXP );
        Rcpp::traits::input_parameter< SEXP >::type vecRangesData(vecRangesDataSEXP );
        Rcpp::traits::input_parameter< SEXP >::type listFunctionParamData(listFunctionParamDataSEXP );
        //SEXP __result = R_NilValue;
	SEXP __result = annotateByLocation(inputData, dfGENETXData, dfGENETXGroupData, dfEXONData, stackEXONData, dfORGDBData, 
		      vecORGDBColumnsData, vecSEPData, vecRangesData, listFunctionParamData);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}



// annotateByAlignment
//SEXP annotateByAlignment(SEXP inputData, SEXP dfGENETXData, SEXP dfGENETXGroupData, SEXP dfEXONData, SEXP stackEXONData, SEXP dfORGDBData, SEXP vecORGDBColumnsData, SEXP vecSEPData, SEXP vecRangesData, SEXP listFunctionParamData);
RcppExport SEXP ProbeAnnotator_annotateByAlignment(SEXP inputDataSEXP, SEXP dfGENETXDataSEXP, SEXP dfGENETXGroupDataSEXP, SEXP dfEXONDataSEXP, SEXP stackEXONDataSEXP, SEXP dfORGDBDataSEXP, SEXP vecORGDBColumnsDataSEXP, SEXP vecSEPDataSEXP, SEXP vecRangesDataSEXP, SEXP listFunctionParamDataSEXP) {
    BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< SEXP >::type inputData(inputDataSEXP);
        Rcpp::traits::input_parameter< SEXP >::type dfGENETXData(dfGENETXDataSEXP);
        Rcpp::traits::input_parameter< SEXP >::type dfGENETXGroupData(dfGENETXGroupDataSEXP);
        Rcpp::traits::input_parameter< SEXP >::type dfEXONData(dfEXONDataSEXP);
        Rcpp::traits::input_parameter< SEXP >::type stackEXONData(stackEXONDataSEXP);
        Rcpp::traits::input_parameter< SEXP >::type dfORGDBData(dfORGDBDataSEXP);
        Rcpp::traits::input_parameter< SEXP >::type vecORGDBColumnsData(vecORGDBColumnsDataSEXP);
        Rcpp::traits::input_parameter< SEXP >::type vecSEPData(vecSEPDataSEXP);
        Rcpp::traits::input_parameter< SEXP >::type vecRangesData(vecRangesDataSEXP);
        Rcpp::traits::input_parameter< SEXP >::type listFunctionParamData(listFunctionParamDataSEXP);
        //SEXP __result = R_NilValue;
        //annotateByLocation(inputData, txDbData, txGroupData, txDbKeyData, txAnnotRangeData, orgDbData, orgDbKeyData, orgDbColumnsData, sepInterData, dfEXONData, stackEXONData);
        SEXP __result = annotateByAlignment(inputData, dfGENETXData, dfGENETXGroupData, dfEXONData, stackEXONData, dfORGDBData,
                vecORGDBColumnsData, vecSEPData, vecRangesData, listFunctionParamData);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
    END_RCPP
}
