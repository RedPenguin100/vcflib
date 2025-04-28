/*
    vcflib C++ library for parsing and manipulating VCF files

    Copyright © 2010-2020 Erik Garrison
    Copyright © 2020      Pjotr Prins

    This software is published under the MIT License. See the LICENSE file.
*/

#ifndef BEDREADER_H
#define BEDREADER_H

#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <iterator>
#include "IntervalTree.h"

using namespace std;


// stores the positional information of a bed target entry
class BedTarget {

public:

    string seq;  // sequence name
    int left;    // left position
    int right;   // right position, adjusted to 0-base
    string desc; // descriptive information, target name typically

    explicit BedTarget(const string& s);

    BedTarget(string s, int l, int r, string d = "")
        : seq(std::move(s))
        , left(l)
        , right(r)
        , desc(std::move(d))
    { }

};



class BedReader {

    bool _isOpen;
    ifstream file;

public:

    bool isOpen(void) { return _isOpen; }

    vector<BedTarget> targets;
    map<string, IntervalTree<size_t, BedTarget*> > intervals; // intervals by reference sequence

    vector<BedTarget> entries(void);

    vector<BedTarget*> targetsContained(const BedTarget& target);

    vector<BedTarget*> targetsOverlapping(const BedTarget& target);

    BedReader(void)
        : _isOpen(false)
    {
    }

    BedReader(string& fname)
        : _isOpen(false) {
        open(fname);
    }

    void addTargets(vector<BedTarget>& targets);

    void open(const string& fname);
};



#endif

