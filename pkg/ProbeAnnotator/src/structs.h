// structs.h
#ifndef INCLUDED_STRUCTS_H
#define INCLUDED_STRUCTS_H
//code...
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>

/*
Mapping Structure
 */
struct Mapping {
    //variables
    std::string fullName;
    std::string gene;
    std::string loctype;
    std::string loctypeNum;
    int start;
    int end;
    std::string txname;

    //constructor

    Mapping(
            std::string const& gene,
            std::string const& loctype,
            int const& start,
            int const& end,
            std::string const& loctypeNum,
            std::string const& txname) {
        this->gene = gene;
        this->loctype = loctype;
        this->start = start;
        this->end = end;
        this->loctypeNum = loctypeNum;
        this->txname = txname;
        //fullName ~ unique ID
        this->fullName = gene + loctype + loctypeNum;
    }

    //operator '<' and '=='

    bool operator<(Mapping const& rhs) const {
        return fullName < rhs.fullName;
    }

    bool operator==(Mapping const& rhs) const {
        return fullName == rhs.fullName;
    }

    //methods

    std::string IntToString(int number) {
        std::ostringstream oss;

        // Works just like cout
        oss << number;

        // Return the underlying string
        return oss.str();
    }

    std::string Start() {
        return IntToString(start);
    }

    std::string End() {
        return IntToString(end);
    }

};

struct Probe {
    //variables
    std::string name;
    std::string setname;
    std::vector<Mapping> mappings;

    //constructor

    Probe(std::string const& name, std::string const& setname) {
        this->name = name;
        this->setname = setname;
    }
};
#endif 
