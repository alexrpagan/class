#include "Variant.h"
#include "split.h"
#include "ewah.h"
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/shared_ptr.hpp>

using namespace std;
using namespace vcf;

namespace po = boost::program_options;
namespace fs = boost::filesystem;

namespace {

  const int SUCCESS = 0;
  const int ERROR = 1;

  const char* HELP_FLAG = "help";
  const char* OUT_FLAG  = "out";
  const char* VCF_FLAG  = "vcf";

  const string VT_DB_FILE = "vt.db";
  string VT_KEY = "VT";

}

typedef EWAHBoolArray<uint32_t> Bitmap;
typedef map<string, boost::shared_ptr<Bitmap> > Genotypes;

int
err_out(const string& err_msg, const int err_code) {
  cerr << err_msg << endl;
  return err_code;
}

void
dump_bitmap_part(
  fs::path outdir
  , const string& sample_name
  , const int cpy_idx
  , boost::shared_ptr<Bitmap> bitmap_part
){
  ofstream bitmap_part_out;
  stringstream fname("");
  fname << sample_name << "-" << cpy_idx;
  fs::path part_path(fname.str());
  bitmap_part_out.open((outdir / part_path).string().c_str());
  bitmap_part->write(bitmap_part_out);
  bitmap_part_out.close();
}

boost::shared_ptr<Bitmap>
get_genotype(Genotypes &genotypes, string sample) {
  if (genotypes.find(sample) == genotypes.end()) {
    boost::shared_ptr<Bitmap> bitmap(new Bitmap());
    genotypes[sample] = bitmap;
  }
  return genotypes[sample];
}

int
main(int argc, char** argv) {

  VariantCallFile variantFile;
  string vcf_file_name, outdir;

  po::options_description desc("Allowed options");
  desc.add_options()
    (HELP_FLAG, "Get some help")
    (OUT_FLAG,  po::value<string>(&outdir),        "Name of output dir")
    (VCF_FLAG,  po::value<string>(&vcf_file_name), "Name of vcf file")
  ;

  po::variables_map vm;

  try
  {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if(vm.count(HELP_FLAG)) {
      cerr << desc << endl;
      return SUCCESS;
    }

    const string outdir_err = "You must supply a valid output directory.";
    if (vm.count(OUT_FLAG)) {
      fs::path outdir_path(outdir);
      if( !fs::exists(outdir_path) && !fs::create_directory(outdir_path) ) {
        return err_out(outdir_err, ERROR);
      }
    } else {
      return err_out(outdir_err, ERROR);
    }

    if (vm.count(VCF_FLAG)) {
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

  ofstream vdb;
  fs::path outdir_path(outdir);
  fs::path vdb_path(VT_DB_FILE);
  vdb.open((outdir_path / vdb_path).string().c_str());

  if (!vdb.is_open()) return ERROR;

  Variant var(variantFile);

  // TODO: make this work for haploid organisms too
  Genotypes cpy1;
  Genotypes cpy2;

  int sampleSize = -1, variantCount = 0;

  while (variantFile.getNextVariant(var)) {

    vdb << variantCount << "\t"
        << var.getInfoValueString(VT_KEY, 0) << "\t"
        << var.position     << "\t"
        << var.ref          << "\t";

    var.printAlt(vdb);
    vdb << "\t" << endl;

    map<string, map<string, vector<string> > >::iterator s     = var.samples.begin();
    map<string, map<string, vector<string> > >::iterator sEnd  = var.samples.end();

    for (; s != sEnd; ++s) {
      boost::shared_ptr<Bitmap> bitmap1 = get_genotype(cpy1, s->first);
      boost::shared_ptr<Bitmap> bitmap2 = get_genotype(cpy2, s->first);

      map<string, vector<string> >& sample = s->second;
      // XXX assumes we can only have one GT value
      string& genotype = sample["GT"].front();
      vector<string> gt = split(genotype, "|/");
      int cpy_idx = 0;
      for (vector<string>::iterator g = gt.begin(); g != gt.end(); ++g, ++cpy_idx) {
        int idx = atoi(g->c_str());
        if (idx != 0 && idx != 1) {
          cerr << "Invalid genotype entry: " << idx << endl;
          cerr << var;
          return ERROR;
        }
        if (idx == 1) {
          switch(cpy_idx) {
            case 0: { bitmap1->set(variantCount); break; }
            case 1: { bitmap2->set(variantCount); break; }
            default: {
              cerr << "Invalid ploidy: " << cpy_idx << endl;
              cerr << var;
              return ERROR;
            }
          }
        }
      }
    }
    variantCount++;
  }

  for(Genotypes::iterator iter = cpy1.begin(); iter != cpy1.end(); ++iter) {
    dump_bitmap_part(outdir_path, iter->first, 1, iter->second);
    dump_bitmap_part(outdir_path, iter->first, 2, cpy2[iter->first]);
  }

  // clean up
  vdb.close();

  return SUCCESS;

}
