#include <fstream> 
#include "annotatorClass.h"

class ByAlignment : public Annotator {
    //own data
private:
    std::vector<Probe> probes;
    std::tr1::unordered_map<std::string, int> hash_probes;
    std::string file2;
    std::string column_separator2;
    std::string probe_separator;
    int score_column_index, probe_column_index, ref_column_index, offset_column_index;
    int min_score;
    int hasProbe;
    int upstream;
    int downstream;

    //construction
public:

    ByAlignment(SEXP listFunctionParamData) : Annotator(listFunctionParamData) {
    }

    //methods
protected:

    void SetInput(SEXP inputData) {
        Rcpp::List input = Rcpp::List(inputData);
        //0 = file1
        //1 = file2
        //2 = separator in file2
        //3=  separator in probeset
        Rcpp::CharacterVector CV(input[0]);
        std::string f = std::string(CV[0]);
        file2 = std::string(CV[1]);
        column_separator2 = std::string(CV[3]);
        probe_separator = std::string(CV[4]);
        if(probe_separator.length() > 0) hasProbe = 1;
        
        //0:3=score_column_index, probe_column_index, ref_column_index,offset_column_index;
        //4=min_score
        //5=upstream
        //6=downstream
        Rcpp::NumericVector NV(input[1]);
        score_column_index = NV[0];
        probe_column_index = NV[1];
        ref_column_index = NV[2];
        offset_column_index = NV[3];
        min_score = NV[4];
        downstream = NV[5];
        upstream = NV[6];



        //file stream
        std::ifstream ifs(f.c_str());
        //string for lines in ifs 
        std::string line;

        //result
        //integer to track the number of probes (distinct)
        int count = 0;

        bool add = true;
        //if header TRUE, skip first line
        //if (header == TRUE)
        //    std::getline(ifs, line);


        //read lines
        while (!ifs.eof()) {
            //clear line
            line.clear();
            //get line
            std::getline(ifs, line);
            if (line[0] == '>') {
                //get probe_name
                std::string probe_name = line.substr(1, line.find_first_of(' ') - 1);
                std::string probe_setname;
                add = true;
                if(hasProbe == 1)
                {
                    // probe set used
                    std::vector<std::string> probe_split;
                    split(probe_split, probe_name, probe_separator);
                    if(probe_split.size() != 2) 
                    {
                        //can't split probe name
                        //-> skip probe
                        add = false;
                    } else {
                        probe_setname = probe_split[0];
                    }
                } else {
                    // probe set not used
                    probe_setname = "";
                }
                //check if already exists in probes
                //if TRUE -> do nothing
                //if FALSE -> add
                std::tr1::unordered_map<std::string, int>::const_iterator got = hash_probes.find(probe_name);
                if (got == hash_probes.end()) {
                    //create probe
                    Probe probe = Probe(probe_name, probe_setname);
                    //add to hashtable
                    hash_probes.insert(std::pair<std::string, int> (probe_name, count));
                    //add to vector
                    probes.push_back(probe);
                    inputName.push_back(probe_name);
                    //increase probe count*/
                    count++;
                }
            }
        }
        ifs.close();
        inputLength = count;
    }

    void GetMapping_INTRONEXON() {
        //file stream
        std::ifstream ifs(file2.c_str());
        //string for lines in ifs 
        std::string line;
        //monitoring
        int conv_notok = 0, conv_ok = 0, score_ok = 0, count_map = 0;

        while (!ifs.eof()) {

            //clear line
            line.clear();
            //get line
            std::getline(ifs, line);
            if (!line.empty()) {
                //we want to extract reference_name, probe_name, score
                std::vector<std::string> line_split;
                split(line_split, line, column_separator2);
                //-get reference_name, probe_name, score
                int score = atoi(line_split[score_column_index].c_str());
                std::string reference_name = line_split[ref_column_index];
                std::string probe_name = line_split[probe_column_index];
                int offset = atoi(line_split[offset_column_index].c_str());
                //check if reference_name already exists in dataset
                //if TRUE -> get reference_name (converted)
                //if FALSE -> error
                std::tr1::unordered_map<std::string, int>::const_iterator got_refconversion = hashTxName.find(reference_name);
                int ref_index = -1;
                if (got_refconversion != hashTxName.end()) {
                    //transcript found...
                    conv_ok++;
                    ref_index = got_refconversion->second;
                } else {
                    //transcript not found...
                    conv_notok++;
                }


                if (score < min_score && ref_index >= 0) {
                    score_ok++;
                    //get probe
                    //if TRUE -> take index
                    //if FALSE -> send error 	
                    int probe_index;
                    std::tr1::unordered_map<std::string, int>::const_iterator got_probe = hash_probes.find(probe_name);
                    if (got_probe == hash_probes.end()) {
                        Rprintf("%s: probe not found!\n", probe_name.c_str());
                        //error
                        return;
                    } else {
                        probe_index = got_probe->second;
                    }


                    Annotator::DoGeneMapping_EXONINTRON(offset + txDbStart[ref_index] - downstream, offset + txDbStart[ref_index] - downstream + 20, 
                            txDbFullStart[ref_index], txDbFullEnd[ref_index],
                            txDbExtendedStart[ref_index], txDbExtendedEnd[ref_index],
                            txDbStart[ref_index], txDbEnd[ref_index],
                            txDbPromotorStart[ref_index], txDbPromotorEnd[ref_index],
                            probes[probe_index].mappings,
                            txDbName[ref_index], txDbTxName[ref_index],
                            exonNumber,
                            exonStart, exonEnd,
                            txDbExonID[ref_index]);


                    //increase maps count
                    count_map++;

                }
            }
        }

        if (verbose) Rprintf("MappingToResult()\n");
        for (unsigned int i = 0; i < probes.size(); ++i) {
            MappingToResult(probes[i].mappings, i);
        }
    }

    void GetMapping_NO_INTRONEXON() {
        //file stream
        std::ifstream ifs(file2.c_str());
        //string for lines in ifs 
        std::string line;
        //monitoring
        int conv_notok = 0, conv_ok = 0, score_ok = 0, count_map = 0;

        while (!ifs.eof()) {

            //clear line
            line.clear();
            //get line
            std::getline(ifs, line);
            if (!line.empty()) {
                //we want to extract reference_name, probe_name, score
                std::vector<std::string> line_split;
                split(line_split, line, column_separator2);
                //-get reference_name, probe_name, score
                int score = atoi(line_split[score_column_index].c_str());
                std::string reference_name = line_split[ref_column_index];
                std::string probe_name = line_split[probe_column_index];
                int offset = atoi(line_split[offset_column_index].c_str());
                //check if reference_name already exists in dataset
                //if TRUE -> get reference_name (converted)
                //if FALSE -> error
                std::tr1::unordered_map<std::string, int>::const_iterator got_refconversion = hashTxName.find(reference_name);
                int ref_index = -1;
                if (got_refconversion != hashTxName.end()) {
                    //transcript found...
                    conv_ok++;
                    ref_index = got_refconversion->second;
                } else {
                    //transcript not found...
                    conv_notok++;
                }


                if (score < min_score && ref_index >= 0) {
                    score_ok++;
                    //get probe
                    //if TRUE -> take index
                    //if FALSE -> send error 	
                    int probe_index;
                    std::tr1::unordered_map<std::string, int>::const_iterator got_probe = hash_probes.find(probe_name);
                    if (got_probe == hash_probes.end()) {
                        //Rcpp::Rcout << probe_name << ": probe not found!\n";
                        //error
                        return;
                    } else {
                        probe_index = got_probe->second;
                    }


                    //do mapping & add to vector
                    Annotator::DoGeneMapping_NO_INTRONEXON(offset + txDbStart[ref_index] - downstream, offset + txDbStart[ref_index] - downstream + 20, 
                            txDbFullStart[ref_index], txDbFullEnd[ref_index],
                            txDbExtendedStart[ref_index], txDbExtendedEnd[ref_index],
                            txDbStart[ref_index], txDbEnd[ref_index],
                            txDbPromotorStart[ref_index], txDbPromotorEnd[ref_index], probes[probe_index].mappings,
                            txDbName[ref_index], txDbTxName[ref_index]);

                    //increase maps count
                    count_map++;

                }
            }
        }

        //Rcpp::Rcout << "MappingToResult()\n";
        for (unsigned int i = 0; i < probes.size(); ++i) {
            MappingToResult(probes[i].mappings, i);
        }
    }

    void GetMapping_EXON_NO_INTRON() {
        //file stream
        std::ifstream ifs(file2.c_str());
        //string for lines in ifs 
        std::string line;
        //monitoring
        int conv_notok = 0, conv_ok = 0, score_ok = 0, count_map = 0;

        while (!ifs.eof()) {

            //clear line
            line.clear();
            //get line
            std::getline(ifs, line);
            if (!line.empty()) {

                //we want to extract reference_name, probe_name, score
                std::vector<std::string> line_split;
                split(line_split, line, column_separator2);
                //-get reference_name, probe_name, score
                int score = atoi(line_split[score_column_index].c_str());
                std::string reference_name = line_split[ref_column_index];
                std::string probe_name = line_split[probe_column_index];
                int offset = atoi(line_split[offset_column_index].c_str());
                //check if reference_name already exists in dataset
                //if TRUE -> get reference_name (converted)
                //if FALSE -> error
                std::tr1::unordered_map<std::string, int>::const_iterator got_refconversion = hashTxName.find(reference_name);
                int ref_index = -1;
                if (got_refconversion != hashTxName.end()) {
                    //transcript found...
                    conv_ok++;
                    ref_index = got_refconversion->second;
                } else {
                    //transcript not found...
                    //Rcpp::Rcout << "reference not found \t" << reference_name << std::endl;
                    conv_notok++;
                }


                if (score < min_score && ref_index >= 0) {
                    score_ok++;
                    //get probe
                    //if TRUE -> take index
                    //if FALSE -> send error 	
                    int probe_index;
                    std::tr1::unordered_map<std::string, int>::const_iterator got_probe = hash_probes.find(probe_name);
                    if (got_probe == hash_probes.end()) {
                        //Rcpp::Rcout << probe_name << ": probe not found!\n";
                        //error
                        return;
                    } else {
                        probe_index = got_probe->second;
                    }

                    Annotator::DoGeneMapping_EXON_NO_INTRON(offset + txDbStart[ref_index] - downstream, offset + txDbStart[ref_index] - downstream + 20, 
                            txDbFullStart[ref_index], txDbFullEnd[ref_index],
                            txDbExtendedStart[ref_index], txDbExtendedEnd[ref_index],
                            txDbStart[ref_index], txDbEnd[ref_index],
                            txDbPromotorStart[ref_index], txDbPromotorEnd[ref_index],
                            probes[probe_index].mappings,
                            txDbName[ref_index], txDbTxName[ref_index],
                            exonNumber,
                            exonStart, exonEnd,
                            txDbExonID[ref_index]);

                    //increase maps count
                    count_map++;

                }
            }
        }


        //Rcpp::Rcout << "\t" << score_ok << "\t" << conv_ok << "\t" << conv_notok << std::endl;
        for (unsigned int i = 0; i < probes.size(); ++i) {
            MappingToResult(probes[i].mappings, i);
        }
    }

private:

    template <typename Container>
    Container& split(Container& result, const typename Container::value_type& s, const typename Container::value_type& delimiters) {
        result.clear();
        size_t current;
        size_t next = -1;
        do {
            current = next + 1;
            next = s.find_first_of(delimiters, current);
            result.push_back(s.substr(current, next - current));
        } while (next != Container::value_type::npos);
        return result;
    }
};

SEXP annotateByAlignment(SEXP inputData, SEXP dfGENETXData, SEXP dfGENETXGroupData, SEXP dfEXONData, SEXP stackEXONData, SEXP dfORGDBData,
        SEXP vecORGDBColumnsData, SEXP vecSEPData, SEXP vecRangesData, SEXP listFunctionParamData) {

    ByAlignment annot = ByAlignment(listFunctionParamData);
    return annot.Annotate(inputData, dfGENETXData, dfGENETXGroupData, dfEXONData, stackEXONData, dfORGDBData,
            vecORGDBColumnsData, vecSEPData, vecRangesData);
}
