#include <iostream>
#include <sstream>  // Required for stringstreams
#include <fstream> 
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <Rcpp.h>

std::string IntToString ( int number )
{
  std::ostringstream oss;

  // Works just like cout
  oss<< number;

  // Return the underlying string
  return oss.str();
}

template <typename T>
    bool IsInBound(const T& value, const T& low, const T& high) {
    return !(value < low) && (value < high);
}

struct mapping {
	mapping(std::string const& f, std::string const&l, std::string const& s, std::string const& e) {
		gene = f;
		loctype = l;
		start = s;
		end = e;
	}
  bool operator<( mapping const& rhs ) const
     { return gene < rhs.gene; }
  bool operator==( mapping const& rhs ) const
     { return gene == rhs.gene; }


  std::string gene;
  std::string loctype;
  std::string start;
  std::string end;
};


SEXP annotateByLocation(SEXP inputData, SEXP txDbData, SEXP txGroupData, SEXP txDbKeyData, SEXP txAnnotRangeData, 
																	SEXP orgDbData, SEXP orgDbKeyData, SEXP orgDbColumnsData)  {
	Rcpp::DataFrame input = Rcpp::DataFrame(inputData); 
	Rcpp::DataFrame txDb = Rcpp::DataFrame(txDbData);
	Rcpp::DataFrame txGroup = Rcpp::DataFrame(txGroupData);
	Rcpp::CharacterVector txDbKey = Rcpp::CharacterVector(txDbKeyData);
	Rcpp::IntegerVector txAnnotRange = Rcpp::IntegerVector(txAnnotRangeData);
	Rcpp::DataFrame orgDb = Rcpp::DataFrame(orgDbData);
	Rcpp::CharacterVector orgDbKey = Rcpp::CharacterVector(orgDbKeyData);
	Rcpp::CharacterVector orgDbColumns = Rcpp::CharacterVector(orgDbColumnsData);


	std::string key = std::string(txDbKey[0]);
	//from input
	//'id', 'chr', 'strand', 'probeEnd', 'start' ,'end', 'txGroup'
	Rcpp::CharacterVector inputStrand = input["strand"];
	Rcpp::CharacterVector inputName = input["ID"];
	Rcpp::CharacterVector inputChrom = input["chr"];
	Rcpp::IntegerVector inputStart = input["start"];
	Rcpp::IntegerVector inputEnd = input["end"];
	Rcpp::CharacterVector inputGroup = input["txGroup"];
	
	//from txDb
	//txDbKey[0], 'chr', 'strand', 'probeEnd', 'start' ,'end', 'txGroup'
	Rcpp::CharacterVector txDbStrand = txDb["strand"]; 
	Rcpp::CharacterVector txDbName = txDb[ key.c_str() ]; 
	Rcpp::CharacterVector txDbChrom = txDb["chr"]; 
	Rcpp::IntegerVector txDbStart = txDb["start"];
	Rcpp::IntegerVector txDbEnd = txDb["end"];
	Rcpp::CharacterVector txDbGroup = txDb["txGroup"];
	Rcpp::IntegerVector txDbPromotorStart = Rcpp::IntegerVector(txDbName.length());
	Rcpp::IntegerVector txDbPromotorEnd = Rcpp::IntegerVector(txDbName.length());
	Rcpp::IntegerVector txDbExtendedStart = Rcpp::IntegerVector(txDbName.length());
	Rcpp::IntegerVector txDbExtendedEnd = Rcpp::IntegerVector(txDbName.length());
	for(int i=0; i < txDbName.length(); ++i)
	{
		txDbExtendedStart[i] = txDbStart[i]-txAnnotRange[1];
		txDbExtendedEnd[i] = txDbEnd[i]+txAnnotRange[1];
		if(txDbStrand[i] == "+") {
			txDbPromotorStart[i] = txDbStart[i]-txAnnotRange[0];
			txDbPromotorEnd[i] = txDbStart[i]+txAnnotRange[0];
		} else {
			txDbPromotorStart[i] = txDbEnd[i]-txAnnotRange[0];
			txDbPromotorEnd[i] = txDbEnd[i]+txAnnotRange[0];
		}
	}
	//from txGroup
	//'group', 'index1', 'index2'
	Rcpp::CharacterVector txGroupChrom = txGroup["group"];
	Rcpp::IntegerVector txGroupIndex1 = txGroup["index1"];
	Rcpp::IntegerVector txGroupIndex2 = txGroup["index2"];
	std::tr1::unordered_map<std::string, int> hashTxGroup;
	for(int i = 0; i < txGroupChrom.length(); ++i)
	{
		hashTxGroup.insert( std::pair<std::string,int> (std::string(txGroupChrom[i]), i) );
	}	


	//from orgDb
	//orgDbKey[0], orgDbColumns
	std::tr1::unordered_map<std::string, int> hashOrgDb;
	Rcpp::CharacterVector orgDbName = orgDb[ std::string(orgDbKey[0]) ];
	std::vector<Rcpp::CharacterVector> orgDbVector;
	for(int i = 1; i < orgDbColumns.length(); ++i)
	{
		orgDbVector.push_back( orgDb[ std::string(orgDbColumns[i]) ] );
	}
	for(int i = 0; i < orgDbName.length(); ++i)
	{
		hashOrgDb.insert( std::pair<std::string,int> (std::string(orgDbName[i]), i) );
	}
	//result...
	Rcpp::CharacterVector vec_annot = Rcpp::CharacterVector(inputName.length());
	Rcpp::CharacterVector vec_loctype = Rcpp::CharacterVector(inputName.length());
	Rcpp::CharacterVector vec_start = Rcpp::CharacterVector(inputName.length());
	Rcpp::CharacterVector vec_end = Rcpp::CharacterVector(inputName.length());
	std::vector<Rcpp::CharacterVector> res_orgDb;
	for(int i = 1; i < orgDbColumns.length(); ++i)
	{
		Rcpp::CharacterVector vec = Rcpp::CharacterVector(inputName.length());
		res_orgDb.push_back( vec );
	}	

	
	for (int i = 0; i < inputGroup.length(); i++) {

		std::tr1::unordered_map<std::string,int>::const_iterator got = hashTxGroup.find( std::string(inputGroup[i]) );				
		if ( got != hashTxGroup.end() )  {
			int from = txGroupIndex1[got->second];
			int to = txGroupIndex2[got->second];

			std::vector<mapping> temp1;
		
			for(int j = from; j <= to; ++j)
			{
				if(IsInBound(inputStart[i], txDbExtendedStart[j],txDbExtendedEnd[j]) || IsInBound(inputEnd[i], txDbExtendedStart[j],txDbExtendedEnd[j])) {

					if(IsInBound(inputStart[i], txDbStart[j] ,txDbEnd[j]) ||  IsInBound(inputEnd[i], txDbStart[j] ,txDbEnd[j]))
					{
							temp1.push_back(mapping(std::string(txDbName[j]), "gene", IntToString(txDbStart[j]), IntToString(txDbEnd[j]) ));
					}
					else if(IsInBound(inputStart[i], txDbPromotorStart[j],txDbPromotorEnd[j]) || IsInBound(inputEnd[i], txDbPromotorStart[j],txDbPromotorEnd[j]))
					{
							temp1.push_back(mapping(std::string(txDbName[j]), "promoter", IntToString(txDbStart[j]), IntToString(txDbEnd[j]) ));
					} else  
					{
							temp1.push_back(mapping(std::string(txDbName[j]), "extended", IntToString(txDbStart[j]), IntToString(txDbEnd[j]) ));
					}
				}
			}


			sort(temp1.begin(), temp1.end());
			temp1.erase(unique(temp1.begin(), temp1.end()), temp1.end());

			//std::string temp_gene = temp1[0].gene;
			//std::string temp_loctype = temp1[0].loctype;
			//std::string temp_start = temp1[0].start;
			//std::string temp_end  = temp1[0].end;

			if(temp1.size() > 0) {



				std::string temp_gene = temp1[0].gene;
				std::string temp_loctype = temp1[0].loctype;
				std::string temp_start = temp1[0].start;
				std::string temp_end  = temp1[0].end;

				std::vector<std::string> orgDbColumnTemp;

				std::tr1::unordered_map<std::string,int>::const_iterator got_gene = hashOrgDb.find( temp1[0].gene );				
				if ( got_gene != hashOrgDb.end() )  {
					for(int j = 1; j < orgDbColumns.length(); ++j)
					{
						orgDbColumnTemp.push_back( std::string( orgDbVector[j-1][got_gene->second] ) );
					}
				} else {
					for(int j = 1; j < orgDbColumns.length(); ++j)
					{
						orgDbColumnTemp.push_back( "NA" );
					}
				}
				temp1.erase(temp1.begin() + 0);

				//for (mapping& s : temp1) { 
				for (int k = 0; k < temp1.size(); ++k) {
					mapping s= temp1[k]; 
					temp_gene += "|" + s.gene; 
					temp_loctype += "|" + s.loctype; 
					temp_start += "|" + s.start; 
					temp_end += "|" + s.end; 	
					std::tr1::unordered_map<std::string,int>::const_iterator got_gene2 = hashOrgDb.find( s.gene );				
					if ( got_gene2 != hashOrgDb.end() )  {
						for(int j = 1; j < orgDbColumns.length(); ++j)
						{
							orgDbColumnTemp[j-1] += "|" + orgDbVector[j-1][got_gene2->second];
						}
					} else {
						for(int j = 1; j < orgDbColumns.length(); ++j)
						{
							orgDbColumnTemp[j-1] += "|NA";
						}
					}
				}
				
				vec_start[i] = temp_start;
				vec_end[i] =  temp_end;
				for(int j = 1; j < orgDbColumns.length(); ++j)
				{
					res_orgDb[j-1][i] = orgDbColumnTemp[j-1];
				}
				vec_annot[i] = temp_gene;
				vec_loctype[i] = temp_loctype;
			} else {
				vec_start[i] = NA_INTEGER;
				vec_end[i] =  NA_INTEGER;
				for(int j = 1; j < orgDbColumns.length(); ++j)
				{
					res_orgDb[j-1][i] =  R_NaString;
				}
				vec_annot[i] =  R_NaString;
				vec_loctype[i] =  "intrageneic";
			}
			
		} else {
			//error -group not found
			return(R_NilValue);
		}
	}


	Rcpp::DataFrame df = Rcpp::DataFrame::create();
	df[ std::string(txDbKey[0]) ] = vec_annot;
	df[ "loctype" ] = vec_loctype;
	df[ "start" ] = vec_start;
	df[ "end" ] = vec_end;
	for(int i = 1; i < orgDbColumns.length(); ++i)
	{
		df[ std::string(orgDbColumns[i]) ] = res_orgDb[i-1];
	}	

	return df;
}



// annotateByLocation
SEXP annotateByLocation(SEXP inputData, SEXP txDbData, SEXP txGroupData, SEXP txDbKeyData, SEXP txAnnotRangeData, SEXP orgDbData, SEXP orgDbKeyData, SEXP orgDbColumnsData);
RcppExport SEXP ProbeAnnotator_annotateByLocation(SEXP inputDataSEXP, SEXP txDbDataSEXP, SEXP txGroupDataSEXP, SEXP txDbKeyDataSEXP, SEXP txAnnotRangeDataSEXP, SEXP orgDbDataSEXP, SEXP orgDbKeyDataSEXP, SEXP orgDbColumnsDataSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< SEXP >::type inputData(inputDataSEXP );
        Rcpp::traits::input_parameter< SEXP >::type txDbData(txDbDataSEXP );
        Rcpp::traits::input_parameter< SEXP >::type txGroupData(txGroupDataSEXP );
        Rcpp::traits::input_parameter< SEXP >::type txDbKeyData(txDbKeyDataSEXP );
        Rcpp::traits::input_parameter< SEXP >::type txAnnotRangeData(txAnnotRangeDataSEXP );
        Rcpp::traits::input_parameter< SEXP >::type orgDbData(orgDbDataSEXP );
        Rcpp::traits::input_parameter< SEXP >::type orgDbKeyData(orgDbKeyDataSEXP );
        Rcpp::traits::input_parameter< SEXP >::type orgDbColumnsData(orgDbColumnsDataSEXP );
        SEXP __result = annotateByLocation(inputData, txDbData, txGroupData, txDbKeyData, txAnnotRangeData, orgDbData, orgDbKeyData, orgDbColumnsData);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
