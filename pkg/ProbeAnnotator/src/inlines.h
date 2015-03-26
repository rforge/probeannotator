#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include "structs.h"
#include <Rcpp.h>
/*
is value within bounds
*/
inline bool IsInBounds(const int& value, const int& low, const int& high) {
		return !(value < low) && (value < high);
}




/*
do gene mapping WITHOUT exons & introns
*/
inline void DoGeneMapping_NO_INTRONEXON(int inputStart, int inputEnd, int txDbFullStart, int txDbFullEnd, 
					int txDbExtendedStart, int txDbExtendedEnd, int txDbStart, int txDbEnd, int txDbPromotorStart, int txDbPromotorEnd, 
					std::vector<Mapping> mappings,
					std::string geneName, std::string txName
 				      )
{
	if(IsInBounds(inputStart, txDbFullStart,txDbFullEnd) || IsInBounds(inputEnd, txDbFullStart,txDbFullEnd)) 
	{
		if(IsInBounds(inputStart, txDbStart ,txDbEnd) ||  IsInBounds(inputEnd, txDbStart ,txDbEnd))
		{
			mappings.push_back( Mapping(geneName, "gene", txDbStart, txDbEnd, "", txName) );
		}
		else if(IsInBounds(inputStart, txDbPromotorStart,txDbPromotorEnd) || IsInBounds(inputEnd, txDbPromotorStart,txDbPromotorEnd))
		{
			mappings.push_back( Mapping(geneName, "promotor", txDbStart, txDbEnd, "", txName) );
		} 
		else if(IsInBounds(inputStart, txDbExtendedStart,txDbExtendedEnd) || IsInBounds(inputEnd, txDbExtendedStart,txDbExtendedEnd)) 
	  
		{
			mappings.push_back( Mapping(geneName, "extended", txDbStart, txDbEnd, "", txName) );
		}
	}
}

/*
do gene mapping WITH exons & introns
*/
inline void __attribute__((__always_inline__)) DoGeneMapping_EXONINTRON(int inputStart, int inputEnd, int txDbFullStart, int txDbFullEnd,  
					int txDbExtendedStart, int txDbExtendedEnd, int txDbStart, int txDbEnd, int txDbPromotorStart, int txDbPromotorEnd, 
					std::vector<Mapping> &mappings,
					std::string geneName, std::string txName, 
					std::vector<std::string> &exonNum,
					Rcpp::IntegerVector &exonStart, Rcpp::IntegerVector &exonEnd, 
					std::vector<int> &exonID)
					//std::tr1::unordered_map<std::string,std::vector<int> > &hashEXON)
	{
	if(IsInBounds(inputStart, txDbFullStart,txDbFullEnd) || IsInBounds(inputEnd, txDbFullStart,txDbFullEnd)) 
	{

	//--------IF GENE
	if(IsInBounds(inputStart, txDbStart ,txDbEnd) ||  IsInBounds(inputEnd, txDbStart ,txDbEnd))
	{
		//std::tr1::unordered_map<std::string,std::vector<int> >::const_iterator gotEXON = hashEXON.find( txName );
		//if(gotEXON != hashEXON.end()) {
		if(exonID.size() > 0) 
		{
		//std::vector<int> temp_vec1 = gotEXON->second;
		unsigned int k1;
		if(IsInBounds(inputStart, txDbStart ,exonStart[exonID[0]]) ||  IsInBounds(inputEnd, txDbStart ,exonStart[exonID[0]]))
		{
			mappings.push_back( Mapping(geneName, "intron", txDbStart, txDbEnd, "0" , txName) );
					return;
		} else
		{
			for(k1 = 0; k1 < exonID.size()-1; ++k1)
			{
				int k = exonID[k1];
				if(IsInBounds(inputStart, exonStart[k] ,exonEnd[k]) ||  IsInBounds(inputEnd, exonStart[k] ,exonEnd[k]))
				{
					//EXON
					mappings.push_back( Mapping(geneName, "exon", txDbStart, txDbEnd, exonNum[k] , txName) );
					return;
				} 
				int kp = exonID[k1+1];
				if(IsInBounds(inputStart, exonEnd[k] ,exonStart[kp]) ||  IsInBounds(inputEnd, exonEnd[k] ,exonStart[kp])) 
				{
					//INTRON;
					mappings.push_back( Mapping(geneName, "intron", txDbStart, txDbEnd, "k" , txName) );
					return;
				}
			}	  
		}
		if(k1 == exonID.size()-1)
		{
			//EXON
			mappings.push_back( Mapping(geneName, "exon", txDbStart, txDbEnd, exonNum[ exonID[k1-1] ] , txName) );
			return;
		}
		} 
		else {		
			//GENE
			mappings.push_back( Mapping(geneName, "gene", txDbStart, txDbEnd, "", txName) );
			return;
		}
	} //-------IF PROMOTOR
	else if(IsInBounds(inputStart, txDbPromotorStart,txDbPromotorEnd) || IsInBounds(inputEnd, txDbPromotorStart,txDbPromotorEnd))
	{
		mappings.push_back( Mapping(geneName, "promotor", txDbStart, txDbEnd, "", txName) );
		return;
	} 
	else if(IsInBounds(inputStart, txDbExtendedStart,txDbExtendedEnd) || IsInBounds(inputEnd, txDbExtendedStart,txDbExtendedEnd)) 
	{
		mappings.push_back( Mapping(geneName, "extended", txDbStart, txDbEnd, "", txName) );
		return;
	}
	} 
}



/*
do gene mapping WITH exons & introns
*/
inline void DoGeneMapping_EXON_NO_INTRON(int inputStart, int inputEnd, int txDbFullStart, int txDbFullEnd,
					int txDbExtendedStart, int txDbExtendedEnd, int txDbStart, int txDbEnd, int txDbPromotorStart, int txDbPromotorEnd, 
					std::vector<Mapping> &mappings,
					std::string geneName, std::string txName, 
					std::vector<std::string> &exonNum,
					Rcpp::IntegerVector &exonStart, Rcpp::IntegerVector &exonEnd, 
					std::vector<int> &exonID)
					//std::tr1::unordered_map<std::string,std::vector<int> > &hashEXON)
	{
	if(IsInBounds(inputStart, txDbFullStart,txDbFullEnd) || IsInBounds(inputEnd, txDbFullStart,txDbFullEnd)) {

	//--------IF GENE
	if(IsInBounds(inputStart, txDbStart ,txDbEnd) ||  IsInBounds(inputEnd, txDbStart ,txDbEnd))
	{
	//std::tr1::unordered_map<std::string,std::vector<int> >::const_iterator gotEXON = hashEXON.find( txName );
	//if(gotEXON != hashEXON.end()) 
	if(exonID.size() > 0)
	{
		//Rcpp::Rcout << "exon in\n";
		//std::vector<int> temp_vec1 = gotEXON->second;
		unsigned int k1;
		
		int temp_start = txDbStart;
		int l;
		for(k1 = 0; k1 < exonID.size()-1; ++k1)
		{
			int k = exonID[k1];
			l = exonEnd[k] - exonStart[k];
			if(IsInBounds(inputStart, temp_start, temp_start + l) ||  IsInBounds(inputEnd, temp_start , temp_start + l))
			{
				//EXON
				mappings.push_back( Mapping(geneName, "exon", txDbStart, txDbEnd, exonNum[k] , txName) );
				break;
			} 
			
			temp_start += l;
		}
		if(k1 == exonID.size())
		{
			//EXON
			mappings.push_back( Mapping(geneName, "exon", txDbStart, txDbEnd, exonNum[ exonID[k1-1] ] , txName) );
		}		
	} else {			
		//GENE
		mappings.push_back( Mapping(geneName, "gene", txDbStart, txDbEnd, "", txName) );
	}
	} //-------IF PROMOTOR
	else if(IsInBounds(inputStart, txDbPromotorStart,txDbPromotorEnd) || IsInBounds(inputEnd, txDbPromotorStart,txDbPromotorEnd))
	{
		mappings.push_back( Mapping(geneName, "promotor", txDbStart, txDbEnd, "", txName) );
	} 
       	else if(IsInBounds(inputStart, txDbExtendedStart,txDbExtendedEnd) || IsInBounds(inputEnd, txDbExtendedStart,txDbExtendedEnd)) 
	{
		mappings.push_back( Mapping(geneName, "extended", txDbStart, txDbEnd, "", txName) );
	}
	}
}
