#include "cigar.hpp"

#include <iostream>

namespace vcflib {

// generates cigar from allele parsed by parsedAlternates
// Note: this function is not used in vcflib
string varCigar(const vector<VariantAllele>& vav, bool xForMismatch) {
    string cigar;
    pair<int, string> element;
    for (const auto& va : vav) {
        if (va.ref != va.alt) {
            if (element.second == "M") {
                cigar += convert(element.first) + element.second;
                element.second = ""; element.first = 0;
            }
            if (va.ref.size() == va.alt.size()) {
                cigar += convert(va.ref.size()) + (xForMismatch ? "X" : "M");
            } else if (va.ref.size() > va.alt.size()) {
                cigar += convert(va.ref.size() - va.alt.size()) + "D";
            } else {
                cigar += convert(va.alt.size() - va.ref.size()) + "I";
            }
        } else {
            if (element.second == "M") {
                element.first += va.ref.size();
            } else {
                element = make_pair(va.ref.size(), "M");
            }
        }
    }
    if (element.second == "M") {
        cigar += convert(element.first) + element.second;
    }
    element.second = ""; element.first = 0;
    return cigar;
}

string mergeCigar(const string& c1, const string& c2) {
    vector<pair<int, char> > cigar1 = splitCigar(c1);
    vector<pair<int, char> > cigar2 = splitCigar(c2);
    // check if the middle elements are the same
    if (cigar1.back().second == cigar2.front().second) {
        cigar1.back().first += cigar2.front().first;
        cigar2.erase(cigar2.begin());
    }
    for (const auto& c : cigar2) {
        cigar1.push_back(c);
    }
    return joinCigar(cigar1);
}

vector<pair<int, char> > splitUnpackedCigar(const string& cigarStr) {
    vector<pair<int, char> > cigar;
    int num = 0;
    char type = cigarStr[0];
    // cerr << "[" << cigarStr << "]" << endl; // 18,12,14
    for (char c: cigarStr) {
        // cerr << "[" << c << "]";
        if (isdigit(c)) {
          cerr << "Is this a valid unpacked CIGAR? <" << cigarStr << ">?" << endl;
          exit(1);
        }
        if (c != type) {
          cigar.push_back(make_pair(num, type));
          //cerr << num << ":" << type << ", ";
          type = c;
          num = 0;
        }
        num += 1;
    }
    cigar.push_back(make_pair(num, type));
    //cerr << num << ":" << type << ", ";
    return cigar;
}

vector<pair<int, char> > splitCigar(const string& cigarStr) {
    vector<pair<int, char> > cigar;
    string number;
    char type = '\0';
    // strings go [Number][Type] ...
    for (const auto c : cigarStr) {
        if (isdigit(c)) {
            if (type == '\0') {
                number += c;
            } else {
                // signal for next token, push back the last pair, clean up
                cigar.push_back(make_pair(atoi(number.c_str()), type));
                number.clear();
                type = '\0';
                number += c;
            }
        } else {
            type = c;
        }
    }
    if (!number.empty() && type != '\0') {
        cigar.push_back(make_pair(atoi(number.c_str()), type));
    }
    return cigar;
}

list<pair<int, char> > splitCigarList(const string& cigarStr) {
    list<pair<int, char> > cigar;
    string number;
    char type = '\0';
    // strings go [Number][Type] ...
    for (const auto c : cigarStr) {
        if (isdigit(c)) {
            if (type == '\0') {
                number += c;
            } else {
                // signal for next token, push back the last pair, clean up
                cigar.push_back(make_pair(atoi(number.c_str()), type));
                number.clear();
                type = '\0';
                number += c;
            }
        } else {
            type = c;
        }
    }
    if (!number.empty() && type != '\0') {
        cigar.push_back(make_pair(atoi(number.c_str()), type));
    }
    return cigar;
}

vector<pair<int, char> > cleanCigar(const vector<pair<int, char> >& cigar) {
    vector<pair<int, char> > cigarClean;
    for (const auto& c : cigar) {
        if (c.first > 0) {
            cigarClean.push_back(c);
        }
    }
    return cigarClean;
}

string joinCigar(const vector<pair<int, char> >& cigar) {
    string cigarStr;
    bool has_error = false;
    for (auto c: cigar) {
        auto len = c.first;
        if (len < 0) has_error = true;
        if (len != 0) {
            cigarStr += convert(len) + c.second;
        }
    }
    if (has_error) {
        cerr << "ERROR: joinCigar creates illegal cigar " << cigarStr << endl;
        exit(1);
    }
    return cigarStr;
}

string joinCigarList(const list<pair<int, char> >& cigar) {
    string cigarStr;
    for (const auto& c : cigar) {
        cigarStr += convert(c.first) + c.second;
    }
    return cigarStr;
}

int cigarRefLen(const vector<pair<int, char> >& cigar) {
    int len = 0;
    for (const auto& c : cigar) {
        if (c.second == 'M' || c.second == 'D' || c.second == 'X') {
            len += c.first;
        }
    }
    return len;
}

bool isEmptyCigarElement(const pair<int, char>& elem) {
    return elem.first == 0;
}

vector<pair<int, string> > old_splitCigar(const string& cigarStr) {
    vector<pair<int, string> > cigar;
    string number;
    string type;
    // strings go [Number][Type] ...
    for (const auto c : cigarStr) {
        if (isdigit(c)) {
            if (type.empty()) {
                number += c;
            } else {
                // signal for next token, push back the last pair, clean up
                cigar.push_back(make_pair(atoi(number.c_str()), type));
                number.clear();
                type.clear();
                number += c;
            }
        } else {
            type += c;
        }
    }
    if (!number.empty() && !type.empty()) {
        cigar.push_back(make_pair(atoi(number.c_str()), type));
    }
    return cigar;
}

string old_joinCigar(const vector<pair<int, string> >& cigar) {
    string cigarStr;
    for (const auto& c : cigar) {
        if (c.first) {
            cigarStr += convert(c.first) + c.second;
        }
    }
    return cigarStr;
}

string old_joinCigar(const vector<pair<int, char> >& cigar) {
    string cigarStr;
    for (const auto& c : cigar) {
        if (c.first) {
            cigarStr += convert(c.first) + string(1, c.second);
        }
    }
    return cigarStr;
}


}
