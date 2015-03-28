
#include "annotatorClass.h"

class ByLocation : public Annotator {
    //own data
private:
    Rcpp::CharacterVector inputGroup;
    Rcpp::CharacterVector inputChrom;
    Rcpp::CharacterVector inputStrand;
    std::vector<int> inputStart;
    std::vector<int> inputEnd;

    //construction
public:
    ByLocation(SEXP listFunctionParamData) : Annotator(listFunctionParamData) { }

protected:
    //methods
    void SetInput(SEXP inputData) {
        Rcpp::DataFrame input = Rcpp::DataFrame(inputData);
        inputStrand = input["strand"];
        inputName = Rcpp::as<std::vector<std::string> >(input["ID"]);
        inputChrom = input["chr"];
        inputStart = Rcpp::as<std::vector<int> >(input["start"]);
        inputEnd = Rcpp::as<std::vector<int> >(input["end"]);
        inputGroup = input["txGroup"];

        inputLength = inputStart.size();
    }

    void GetMapping_INTRONEXON() {
        for (int i = 0; i < inputGroup.length(); i++) {
            std::tr1::unordered_map<std::string, int>::const_iterator got = hashTxGroup.find(std::string(inputGroup[i]));
            if (got != hashTxGroup.end()) {
                int from = dfGENETXGroupIndex1[got->second];
                int to = dfGENETXGroupIndex2[got->second];

                std::vector<Mapping> mappings;
                for (int j = from; j <= to; ++j) {
                    DoGeneMapping_EXONINTRON(inputStart[i], inputEnd[i], txDbFullStart[j], txDbFullEnd[j],
                            txDbExtendedStart[j], txDbExtendedEnd[j],
                            txDbStart[j], txDbEnd[j],
                            txDbPromotorStart[j], txDbPromotorEnd[j],
                            mappings,
                            txDbName[j], txDbTxName[j],
                            exonNumber,
                            exonStart, exonEnd,
                            txDbExonID[j]);
                }

                MappingToResult(mappings, i);
            } else {
                //error -group not found
                return;
            }
        }
    }

    void GetMapping_NO_INTRONEXON() {
        for (int i = 0; i < inputGroup.length(); i++) {
            std::tr1::unordered_map<std::string, int>::const_iterator got = hashTxGroup.find(std::string(inputGroup[i]));
            if (got != hashTxGroup.end()) {
                int from = dfGENETXGroupIndex1[got->second];
                int to = dfGENETXGroupIndex2[got->second];

                std::vector<Mapping> mappings;
                for (int j = from; j <= to; ++j) {
                    DoGeneMapping_NO_INTRONEXON(inputStart[i], inputEnd[i], txDbFullStart[j], txDbFullEnd[j],
                            txDbExtendedStart[j], txDbExtendedEnd[j],
                            txDbStart[j], txDbEnd[j],
                            txDbPromotorStart[j], txDbPromotorEnd[j],
                            mappings,
                            txDbName[j], txDbTxName[j]);
                }

                MappingToResult(mappings, i);
            } else {
                //error -group not found
                return;
            }
        }
    }

    void GetMapping_EXON_NO_INTRON() {
        for (int i = 0; i < inputGroup.length(); i++) {
            std::tr1::unordered_map<std::string, int>::const_iterator got = hashTxGroup.find(std::string(inputGroup[i]));
            if (got != hashTxGroup.end()) {
                int from = dfGENETXGroupIndex1[got->second];
                int to = dfGENETXGroupIndex2[got->second];

                std::vector<Mapping> mappings;
                for (int j = from; j <= to; ++j) {
                    DoGeneMapping_EXON_NO_INTRON(inputStart[i], inputEnd[i], txDbFullStart[j], txDbFullEnd[j],
                            txDbExtendedStart[j], txDbExtendedEnd[j],
                            txDbStart[j], txDbEnd[j],
                            txDbPromotorStart[j], txDbPromotorEnd[j],
                            mappings,
                            txDbName[j], txDbTxName[j],
                            exonNumber,
                            exonStart, exonEnd,
                            txDbExonID[j]);
								}

                MappingToResult(mappings, i);
            } else {
	        			Rprintf("Error: group %s not found", std::string(inputGroup[i]).c_str());
                return;
            }
        }
    }

};

SEXP annotateByLocation(SEXP inputData, SEXP dfGENETXData, SEXP dfGENETXGroupData, SEXP dfEXONData, SEXP stackEXONData, SEXP dfORGDBData,
        SEXP vecORGDBColumnsData, SEXP vecSEPData, SEXP vecRangesData, SEXP listFunctionParamData) {

    ByLocation annot = ByLocation(listFunctionParamData);
    return annot.Annotate(inputData, dfGENETXData, dfGENETXGroupData, dfEXONData, stackEXONData, dfORGDBData,
            vecORGDBColumnsData, vecSEPData, vecRangesData);
}

