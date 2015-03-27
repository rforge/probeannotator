class Annotator {
protected:
    // variables
    // from dfGENETX (TxDb)
    Rcpp::CharacterVector txDbStrand;
    std::vector<std::string> txDbName;
    std::vector<std::string> txDbTxName;
    Rcpp::CharacterVector txDbChrom;
    Rcpp::CharacterVector txDbGroup;
    std::vector<int> txDbStart;
    std::vector<int> txDbEnd;
    std::vector<int> txDbPromotorStart;
    std::vector<int> txDbPromotorEnd;
    std::vector<int> txDbExtendedStart;
    std::vector<int> txDbExtendedEnd;
    std::vector<int> txDbFullStart;
    std::vector<int> txDbFullEnd;
    //from dfEXON 
    std::vector<std::string> exonNumber;
    Rcpp::IntegerVector exonStart;
    Rcpp::IntegerVector exonEnd;
    //from stackEXON 
    std::vector<std::vector<int> > txDbExonID;
    //from dfGENETXGroup
    Rcpp::CharacterVector dfGENETXGroupChrom;
    Rcpp::IntegerVector dfGENETXGroupIndex1;
    Rcpp::IntegerVector dfGENETXGroupIndex2;
    std::tr1::unordered_map<std::string, int> hashTxGroup;
    std::tr1::unordered_map<std::string, int> hashTxName;
    //from dfORGDB
    std::tr1::unordered_map<std::string, int> hashOrgDb;
    Rcpp::CharacterVector vecORGDBColumns;
    Rcpp::CharacterVector orgDbName;
    std::vector<Rcpp::CharacterVector> orgDbVector;
    //result
    Rcpp::CharacterVector vec_annot;
    Rcpp::CharacterVector vec_txname;
    Rcpp::CharacterVector vec_loctype;
    Rcpp::CharacterVector vec_locnum;
    Rcpp::CharacterVector vec_start;
    Rcpp::CharacterVector vec_end;
    std::vector<Rcpp::CharacterVector> res_orgDb;
    //parameters
    std::vector<std::string> vecSEPARATORS;
    int mappingStyle;
    bool verbose;
    int inputLength;
    std::vector<std::string> inputName;

public:

    Annotator(SEXP listFunctionParamData) {
        Rcpp::List listFunctionParam = Rcpp::List(listFunctionParamData);

        //0=mappingStyle
        //1=verbose
        Rcpp::NumericVector NV(listFunctionParam[0]);
        mappingStyle = NV[0];
        verbose = (NV[1] == 1) ? true : false;

    }

    //-mappingStyle: 0=INTRONEXON, 1=NO_INTROEXON, 2=EXON_NO_INTRON

    SEXP Annotate(SEXP inputData, SEXP dfGENETXData, SEXP dfGENETXGroupData, SEXP dfEXONData, SEXP stackEXONData, SEXP dfORGDBData,
            SEXP vecORGDBColumnsData, SEXP vecSEPData, SEXP vecRangesData) {
        //set background data
        if (verbose) Rprintf("SetData()\n");
        this->SetData(dfGENETXData, dfGENETXGroupData, dfEXONData, stackEXONData, dfORGDBData,
                vecORGDBColumnsData, vecSEPData, vecRangesData);

        //get probes
        if (verbose) Rprintf("SetInput()\n");
        this->SetInput(inputData);
        //inputData.delete();

        //set result data
        if (verbose) Rprintf("SetResult(length = %d)\n", inputLength);
        this->SetResult(inputLength);

        //do mapping & group results
        if (verbose) Rprintf("GetMapping(mappingStyle = %d)\n", mappingStyle);
        if (mappingStyle == 0)
            this->GetMapping_INTRONEXON();
        else if (mappingStyle == 1)
            this->GetMapping_NO_INTRONEXON();
        else if (mappingStyle == 2)
            this->GetMapping_EXON_NO_INTRON();

        if (verbose) Rprintf("GetResult()\n");
        return this->GetResult();
    }
protected:

    virtual void GetMapping_INTRONEXON() { }

    virtual void GetMapping_NO_INTRONEXON() { }

    virtual void GetMapping_EXON_NO_INTRON() { }

    virtual void SetInput(SEXP inputData) { }

    void MappingToResult(std::vector<Mapping> mappings, int i) {
        
        if (mappings.size() > 0) {
            //group
            if (mappings.size() > 1) {
                sort(mappings.begin(), mappings.end());
                //mappings.erase(unique(mappings.begin(), mappings.end()), mappings.end());
                std::string lastname = mappings[0].fullName;
                for (unsigned int i = 1; i < mappings.size(); ++i) {
                    lastname = mappings[i - 1].fullName;
                    while (i < mappings.size() && lastname == mappings[i].fullName) {
                        mappings[i - 1].txname += vecSEPARATORS[0] + mappings[i].txname;
                        mappings.erase(mappings.begin() + i);
                    }
                }
            }
						
    				std::string temp_gene = mappings[0].gene;
            std::string temp_loctype = mappings[0].loctype;
            std::string temp_locnum = mappings[0].loctypeNum;
            std::string temp_start = mappings[0].Start();
            std::string temp_end = mappings[0].End();
            std::string temp_txname = mappings[0].txname;

            std::vector<std::string> orgDbColumnTemp;

            std::tr1::unordered_map<std::string, int>::const_iterator got_gene = hashOrgDb.find(mappings[0].gene);
            if (got_gene != hashOrgDb.end()) {
                for (unsigned int j = 0; j < orgDbVector.size(); ++j) {
                    orgDbColumnTemp.push_back(std::string(orgDbVector[j][got_gene->second]));
                }
            } else {
                for (unsigned int j = 0; j < orgDbVector.size(); ++j) {
                    orgDbColumnTemp.push_back("NA");
                }
            }
            mappings.erase(mappings.begin() + 0);
            //for (mapping& s : mappings) { 
            for (unsigned int k = 0; k < mappings.size(); ++k) {
                Mapping s = mappings[k];
                temp_gene += vecSEPARATORS[1] + s.gene;
                temp_txname += vecSEPARATORS[1] + s.txname;
                temp_loctype += vecSEPARATORS[1] + s.loctype;
                temp_locnum += vecSEPARATORS[1] + s.loctypeNum;
                temp_start += vecSEPARATORS[1] + s.Start();
                temp_end += vecSEPARATORS[1] + s.End();
                std::tr1::unordered_map<std::string, int>::const_iterator got_gene2 = hashOrgDb.find(s.gene);
                if (got_gene2 != hashOrgDb.end()) {
                    for (unsigned int j = 0; j < orgDbVector.size(); ++j)
                        orgDbColumnTemp[j] += vecSEPARATORS[1] + orgDbVector[j][got_gene2->second];
                } else {
                    for (unsigned int j = 0; j < orgDbVector.size(); ++j)
                        orgDbColumnTemp[j] += vecSEPARATORS[1] + "NA";
                }
            }

            vec_start[i] = temp_start;
            vec_end[i] = temp_end;
            for (unsigned int j = 0; j < orgDbVector.size(); ++j)
                res_orgDb[j][i] = orgDbColumnTemp[j];
            vec_annot[i] = temp_gene;
            vec_txname[i] = temp_txname;
            vec_loctype[i] = temp_loctype;
            vec_locnum[i] = temp_locnum;
        } else {
            vec_start[i] = NA_INTEGER;
            vec_end[i] = NA_INTEGER;
            vec_locnum[i] = R_NaString;
            //for(int j = 0; j < orgDbVector.size(); ++j)
            //	res_orgDb[j][i] =  R_NaString;
            vec_annot[i] = R_NaString;
            vec_loctype[i] = "intragenic";
            vec_txname[i] = R_NaString;
        }
    }

    void SetResult(int length) {
        if (length >= 0) {
            //result...
            vec_annot = Rcpp::CharacterVector(length);
            vec_txname = Rcpp::CharacterVector(length);
            vec_loctype = Rcpp::CharacterVector(length);
            vec_locnum = Rcpp::CharacterVector(length);
            vec_start = Rcpp::CharacterVector(length);
            vec_end = Rcpp::CharacterVector(length);
            for (unsigned int i = 0; i < orgDbVector.size(); ++i) {
                Rcpp::CharacterVector vec = Rcpp::CharacterVector(length);
                res_orgDb.push_back(vec);
            }
        } else {
            //must have size...
        }
    }

    SEXP GetResult() {
        Rcpp::DataFrame df = Rcpp::DataFrame::create();
        df[ "probename"] = inputName;
        df[ "entrezid" ] = vec_annot;
        df[ "txname" ] = vec_txname;
        df[ "loctype" ] = vec_loctype;
        df[ "locnum" ] = vec_locnum;
        df[ "start" ] = vec_start;
        df[ "end" ] = vec_end;
        for (unsigned int i = 0; i < orgDbVector.size(); ++i) {
            df[ std::string(vecORGDBColumns[i + 1]) ] = res_orgDb[i];
        }

        return df;
    }

    void SetData(SEXP dfGENETXData, SEXP dfGENETXGroupData, SEXP dfEXONData, SEXP stackEXONData, SEXP dfORGDBData, SEXP vecORGDBColumnsData, SEXP vecSEPData, SEXP vecRangesData) {
        if(verbose) Rprintf("Get Data\n");
        //is.data.frame
        Rcpp::DataFrame dfGENETX = Rcpp::DataFrame(dfGENETXData);
        Rcpp::DataFrame dfGENETXGroup = Rcpp::DataFrame(dfGENETXGroupData);
        Rcpp::DataFrame dfEXON = Rcpp::DataFrame(dfEXONData);
        Rcpp::DataFrame dfORGDB = Rcpp::DataFrame(dfORGDBData);
        //is.list
        Rcpp::List stackEXON = Rcpp::List(stackEXONData);
        //is.character
        vecORGDBColumns = Rcpp::CharacterVector(vecORGDBColumnsData);
        Rcpp::CharacterVector vecSEP = Rcpp::CharacterVector(vecSEPData);
        std::vector<int> vecRanges = Rcpp::as<std::vector<int> >(vecRangesData);
        //Rprintf("vecRanges:%d\t%d\t%d\t%d\n", vecRanges[0],vecRanges[1],vecRanges[2],vecRanges[3]);


        //parameters
        for (int i = 0; i < vecSEP.length(); ++i)
            vecSEPARATORS.push_back(std::string(vecSEP[i]));

        //from dfGENETX (TxDb)
        //'GENEID', 'TXNAME', 'chr', 'strand', 'probeEnd', 'start' ,'end', 'txGroup'
        txDbStrand = dfGENETX["strand"];
        txDbName = Rcpp::as<std::vector<std::string> >(dfGENETX[ "GENEID" ]);
        txDbTxName = Rcpp::as<std::vector<std::string> >(dfGENETX[ "TXNAME" ]);
        txDbChrom = dfGENETX["chr"];
        txDbStart = Rcpp::as<std::vector<int> >(dfGENETX["start"]);
        txDbEnd = Rcpp::as<std::vector<int> >(dfGENETX["end"]);
        txDbGroup = dfGENETX["txGroup"];
        txDbPromotorStart = std::vector<int>(txDbStrand.length());
        txDbPromotorEnd = std::vector<int>(txDbStrand.length());
        txDbExtendedStart = std::vector<int>(txDbStrand.length());
        txDbExtendedEnd = std::vector<int>(txDbStrand.length());
        txDbFullStart = std::vector<int>(txDbStrand.length());
        txDbFullEnd = std::vector<int>(txDbStrand.length());
       for (int i = 0; i < txDbStrand.length(); ++i) {
            hashTxName.insert(std::pair<std::string, int> (std::string(txDbTxName[i]), i));
            txDbExtendedStart[i] = txDbStart[i] - vecRanges[1];
            txDbExtendedEnd[i] = txDbEnd[i] + vecRanges[1];
            if (txDbStrand[i] == "+") {
                txDbPromotorStart[i] = txDbStart[i] - vecRanges[0];
                txDbPromotorEnd[i] = txDbStart[i] + vecRanges[0];
                txDbFullStart[i] = txDbStart[i] - vecRanges[2];
                txDbFullEnd[i] = txDbEnd[i] + vecRanges[3];
            } else {
                txDbPromotorStart[i] = txDbEnd[i] - vecRanges[0];
                txDbPromotorEnd[i] = txDbEnd[i] + vecRanges[0];
                txDbFullEnd[i] = txDbEnd[i] + vecRanges[2];
                txDbFullStart[i] = txDbStart[i] - vecRanges[3];
            }
        }

        //from dfEXON 
        exonNumber = Rcpp::as<std::vector<std::string> >(dfEXON["EXONID"]);
        std::tr1::unordered_map<int, int> hashTemp;
        for (unsigned int i = 0; i < exonNumber.size(); ++i) {
            hashTemp.insert(std::pair<int, int > (atoi(exonNumber[i].c_str()), i));
        }
        exonStart = dfEXON["start"];
        exonEnd = dfEXON["end"];
        //from stackEXON 
        txDbExonID = std::vector<std::vector<int> >(stackEXON.length());
        for (int i = 0; i < stackEXON.length(); ++i) {
            std::vector<int> tempo = Rcpp::as<std::vector<int> >(stackEXON[i]);
            txDbExonID[i] = tempo;
        }


        //from dfGENETXGroup
        //'group', 'index1', 'index2'
        dfGENETXGroupChrom = dfGENETXGroup["group"];
        dfGENETXGroupIndex1 = dfGENETXGroup["index1"];
        dfGENETXGroupIndex2 = dfGENETXGroup["index2"];
        for (int i = 0; i < dfGENETXGroupChrom.length(); ++i) {
            hashTxGroup.insert(std::pair<std::string, int> (std::string(dfGENETXGroupChrom[i]), i));
        }


        //from dfORGDB
        //'ENTREZID', orgDbColumns
        orgDbName = dfORGDB[ "ENTREZID" ];
        for (int i = 1; i < vecORGDBColumns.length(); ++i) {
            orgDbVector.push_back(dfORGDB[ std::string(vecORGDBColumns[i]) ]);
        }
        for (int i = 0; i < orgDbName.length(); ++i) {
            hashOrgDb.insert(std::pair<std::string, int> (std::string(orgDbName[i]), i));
        }

    }



};

