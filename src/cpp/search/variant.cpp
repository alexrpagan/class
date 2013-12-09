#include <stdlib.h>
#include <vector>
#include <cassert>
#include <cmath>

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
    _sv_len = (_type == "SV") ? atoi(row[6].c_str()) : 0;
  }

  int GetChrom()   { return _chr;  }
  int GetBit()     { return _bit;  }
  int GetPos()     { return _pos;  }
  string GetRef()  { return _ref;  }
  string GetAlt()  { return _alt;  }
  string GetType() { return _type; }

  bool IsDeletion() { return (_type == "INDEL" || _type == "SV") && _ref.size() > 1; }
  bool IsInsert()   { return _type == "INDEL" && _ref.size() == 1; }

  int GetLengthMod() {
    if ( _type == "SNP" ) {
      return 0;
    } else if( _type == "INDEL" ) {
      return _alt.size() - _ref.size();
    } else if( _type == "SV" ) {
      if (_alt == "<DEL>") {
        return 0;
      }
      return _alt.size() - _ref.size();
    } else {
      assert(false);
    }
    return 0;
  }

  bool subsumes(Variant &var) {
    return var.GetPos() >= _pos && var.GetPos() <= _pos + abs(GetLengthMod());
  }

private:
  int _chr;
  int _bit;
  int _pos;
  int _sv_len;
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