#include <stdlib.h>
#include <vector>

#ifndef VARIANT_CPP
#define VARIANT_CPP

using namespace std;
typedef vector<string> Row;

class Variant {
public:
  Variant() {}
  Variant(Row &row) {
    _chr = atoi(row[0].c_str());
    _bit = atoi(row[1].c_str());
    _pos = atoi(row[3].c_str()) - 1;
    _ref = row[4];
    _alt = row[5];
    _type = row[2];
  }

  int GetChrom()   { return _chr;  }
  int GetBit()     { return _bit;  }
  int GetPos()     { return _pos;  }
  string GetRef()  { return _ref;  }
  string GetAlt()  { return _alt;  }
  string GetType() { return _type; }

  bool subsumes(Variant &var) {
    return (var.GetBit() == (_bit + 1))
      && (_ref.size() > var.GetRef().size())
      && (_ref.find(var.GetRef()) != string::npos);
  }

  bool modifiesRef(Variant &var) {
    if (_type == "SNP" && (var.GetType() == "SV" || var.GetType() == "INDEL")) {
      return var.GetPos() == _pos;
    }
    return false;
  }

private:
  int _chr;
  int _bit;
  int _pos;
  string _ref;
  string _alt;
  string _type;

};

ostream& operator<<(ostream &strm, Variant &var) {
  strm << var.GetChrom() << '\t';
  strm << var.GetBit()   << '\t';
  strm << var.GetPos()   << '\t';
  strm << var.GetRef()   << '\t';
  strm << var.GetAlt()   << '\t';
  strm << var.GetType();
  return strm;
}

#endif