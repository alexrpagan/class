#include "Variant.h"
#include "split.h"
#include "ewah.h"
#include <string>
#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

using namespace std;
using namespace vcf;

namespace po = boost::program_options;
namespace fs = boost::filesystem;

namespace {
  const int SUCCESS = 0;
  const int ERROR = 1;
  const int BUFFER_SIZE   = 1024 * 1024;
  const string VT_DB_FILE = "vt.db";
  string VT_KEY = "VT";
}

int err_out(const string& err_msg, const int err_code) {
  cerr << err_msg << endl;
  return err_code;
}

int main(int argc, char** argv) {

  VariantCallFile variantFile;
  string vcf_file_name, outdir;

  po::options_description desc("Allowed options");
  desc.add_options()
    ("outdir", po::value<string>(&outdir),        "Name of output dir")
    ("vcf",    po::value<string>(&vcf_file_name), "Name of vcf file")
  ;

  po::variables_map vm;

  try
  {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if(vm.count("help")) {
      cerr << desc << endl;
      return SUCCESS;
    }

    const string outdir_err = "You must supply a valid output directory.";
    if (vm.count("outdir")) {
      fs::path outdir_path(outdir);
      if(fs::create_directory(outdir_path)) {
        return err_out(outdir_err, ERROR);
      }
    } else {
      return err_out(outdir_err, ERROR);
    }

    if (vm.count("vcf")) {
      variantFile.open(vcf_file_name);
    } else {
      variantFile.open(std::cin);
    }
    if (!variantFile.is_open()) return ERROR;

  }
  catch(po::error& e)
  {
    cerr << "ERROR: " << e.what() << endl << endl;
    cerr << desc << endl;
    return ERROR;
  }

  ofstream vdb_file;

  Variant var(variantFile);

  int sampleSize = -1;

  while (variantFile.getNextVariant(var)) {

    map<string, map<string, vector<string> > >::iterator s     = var.samples.begin();
    map<string, map<string, vector<string> > >::iterator sEnd  = var.samples.end();

    if (sampleSize < 0) sampleSize = var.samples.size();

    // TODO: dump this information into variant db streamm
    cout << var.sequenceName << "\t"
       << var.position     << "\t"
       << var.ref          << "\t"
       << var.getInfoValueString(VT_KEY, 0) << "\t";

    var.printAlt(cout);     cout << "\t";
    var.printAlleles(cout); cout << "\t";

    for (; s != sEnd; ++s) {
      map<string, vector<string> >& sample = s->second;
      string& genotype = sample["GT"].front(); // XXX assumes we can only have one GT value
      cout << genotype << "\t";
    }
    cout << endl;

  }
  return SUCCESS;

}
