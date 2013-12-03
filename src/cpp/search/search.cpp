// TODO: remove unused INCLUDES

#include <fstream>
#include <sstream>

#include "ewah.h"
#include "tabix.hpp"

#include <ncbi_pch.hpp>
#include <corelib/ncbiapp.hpp>
#include <corelib/ncbienv.hpp>
#include <corelib/ncbiargs.hpp>

#include <objmgr/object_manager.hpp>

#include <objects/seqalign/Seq_align_set.hpp>

#include <algo/blast/api/sseqloc.hpp>
#include <algo/blast/api/local_blast.hpp>
#include <algo/blast/api/bl2seq.hpp>
#include <algo/blast/api/uniform_search.hpp>
#include <algo/blast/api/blast_types.hpp>
#include <algo/blast/api/blast_aux.hpp>
#include <algo/blast/api/objmgr_query_data.hpp>
#include <algo/blast/api/blast_options_handle.hpp>
#include <algo/blast/api/blast_nucl_options.hpp>
#include <algo/blast/api/blast_prot_options.hpp>

#include <algo/blast/blastinput/blast_input.hpp>
#include <algo/blast/blastinput/blast_fasta_input.hpp>

#include <boost/filesystem.hpp>
#include <boost/shared_ptr.hpp>

USING_NCBI_SCOPE;
USING_SCOPE(blast);
USING_SCOPE(objects);

namespace fs = boost::filesystem;

namespace {

  const int SUCCESS = 0;
  const int ERROR = 1;

  const bool IS_PROTEIN = false;

  const string COARSE_EVALUE_FLAG = "coarse_evalue";
  const string BLAST_PROGRAM      = "blastn";
  const string EVALUE_FLAG        = "evalue";
  const string INPUT_FLAG         = "in";
  const string CHR_FLAG           = "chr";
  const string VDB_FLAG           = "vdb";
  const string DB_FLAG            = "db";

  const string VDB_BASE_FILENAME = "vt";
  const string VDB_FILENAME = "vt.db.gz";

}

typedef EWAHBoolArray<uint32_t> Bitmap;
typedef map<string, boost::shared_ptr<Bitmap> > Genotypes;

typedef vector<string> Row;
typedef list<Row> Rows;

class CSearchApplication : public CNcbiApplication
{
private:

  virtual void Init(void);
  virtual int  Run(void);
  virtual void Exit(void);

  CSearchResultSet RunBlast(string dbname, TSeqLocVector query_loc, CRef<CBlastOptionsHandle> opts);
  fs::path getVDBFilePath(const string vdb_dir);
  Genotypes SlurpGenotypes(const string vdb_dir);
  void PrintErrorMessages(CSearchResults &queryResult);
  void PrintQueryResult(CSearchResults &queryResult);
  string getRegion(string chr, int start, int end);
  Bitmap bitmapFromArray(vector<size_t> array);

};

void
CSearchApplication::Init(void)
{

  auto_ptr<CArgDescriptions> arg_desc(new CArgDescriptions);

  arg_desc->SetUsageContext(GetArguments().GetProgramBasename(),
                "Blast over bitvector-compressed genomes");

  arg_desc->AddKey(DB_FLAG, "Database",
    "Location of reference blast database", CArgDescriptions::eString);

  arg_desc->AddKey(VDB_FLAG, "VDB",
    "Location of variant database", CArgDescriptions::eString);

  arg_desc->AddKey(CHR_FLAG, "chr",
    "The chromosome to blast over", CArgDescriptions::eString);

  arg_desc->AddKey(INPUT_FLAG, "QueryFile",
    "FASTA file containing queries", CArgDescriptions::eInputFile);

  arg_desc->AddDefaultKey(EVALUE_FLAG, "evalue",
    "E-value threshold for final hits", CArgDescriptions::eDouble, "1e-30");

  arg_desc->AddDefaultKey(COARSE_EVALUE_FLAG, "CoarseEvalue",
    "E-value threshold for coarse search against reference sequence", CArgDescriptions::eDouble, "1e-20");

  SetupArgDescriptions(arg_desc.release());

}

int
CSearchApplication::Run(void)
{
  const CArgs& args = GetArgs();

  EProgram program = ProgramNameToEnum(BLAST_PROGRAM);

  CRef<CBlastOptionsHandle> opts(CBlastOptionsFactory::Create(program));

  if(args[COARSE_EVALUE_FLAG].AsDouble()) {
    opts->SetEvalueThreshold(args[COARSE_EVALUE_FLAG].AsDouble());
  }
  opts->Validate();

  CRef<CObjectManager> objmgr = CObjectManager::GetInstance();
  if (!objmgr) {
    throw std::runtime_error("Could not initialize object manager");
  }

  // this chunk of madness prepares the query file.
  SDataLoaderConfig dlconfig(IS_PROTEIN);
  CBlastInputSourceConfig iconfig(dlconfig);
  CBlastFastaInputSource fasta_input(args[INPUT_FLAG].AsInputFile(), iconfig);
  CScope scope(*objmgr);
  CBlastInput blast_input(&fasta_input);
  TSeqLocVector query_loc = blast_input.GetAllSeqLocs(scope);

  string dbname  = args[DB_FLAG].AsString();
  string vdbpath = args[VDB_FLAG].AsString();
  string chr     = args[CHR_FLAG].AsString();

  Genotypes genotypes = SlurpGenotypes(vdbpath);

  string vdb_file_name = getVDBFilePath(vdbpath).string();
  Tabix vdb(vdb_file_name);

  // coarse blast.
  CSearchResultSet results = RunBlast(dbname, query_loc, opts);

  for (unsigned int i = 0; i < results.GetNumResults(); i++) {
    CSearchResults &queryResult = results[i];
    PrintErrorMessages(queryResult);
    const list<CRef<CSeq_align> > &seqAlignList = queryResult.GetSeqAlign()->Get();

    if (seqAlignList.size() > 1000) {
      cout << "punted query " << i << endl;
      continue;
    }

    ITERATE(list <CRef<CSeq_align> >, seqAlignIter, seqAlignList) {

      int start = (*seqAlignIter)->GetSeqStart(1);
      int stop  = (*seqAlignIter)->GetSeqStop(1);
      string region = getRegion(chr, start, stop);

      cout << region << endl;

      Rows variants;
      string line;
      vdb.setRegion(region);
      while(vdb.getNextLine(line)) {
        Row row;
        stringstream ss(line);
        string cell;
        while(getline(ss, cell, '\t')) {
          row.push_back(cell);
        }
        variants.push_back(row);
      }

      Bitmap variant_filter;
      for (Rows::iterator it = variants.begin(); it != variants.end(); ++it) {
        int idx = atoi((*it)[1].c_str());
        variant_filter.set(idx);
      }

      map<vector<size_t>, boost::shared_ptr<vector<string> > > variant_genotypes;

      // TODO: parallelize this step
      for (Genotypes::iterator it = genotypes.begin(); it != genotypes.end(); ++it) {
        Bitmap and_result;
        it->second->logicaland(variant_filter, and_result);
        vector<size_t> bits_set = and_result.toArray();

        if (variant_genotypes.find(bits_set) == variant_genotypes.end()) {
          boost::shared_ptr<vector<string> > vec(new vector<string>());
          variant_genotypes.insert(make_pair(bits_set, vec));
        }

        variant_genotypes[bits_set]->push_back(it->first);
      }

      map<vector<size_t>, boost::shared_ptr<vector<string> > >::iterator it   = variant_genotypes.begin();

      cout << variant_filter.numberOfOnes() << " variants (total) " << endl;
      for(; it != variant_genotypes.end(); ++it) {
        cout << bitmapFromArray(it->first) << " variants : " << it->second->size() << " genomes " << endl;
      }

      // Construct FASTA strings from unique variant sets

      // Run fine blast

      // Report results.

    }
  }

  return SUCCESS;
}

Bitmap
CSearchApplication::bitmapFromArray(vector<size_t> array) {
  Bitmap bits;
  for (unsigned int i = 0; i < array.size(); ++i) {
    bits.set(array[i]);
  }
  return bits;
}

string
CSearchApplication::getRegion(string chr, int start, int end) {
  ostringstream ss;
  ss << chr << ":" << start << "-" << end;
  return ss.str();
}

fs::path
CSearchApplication::getVDBFilePath(const string vdb_dir) {
  fs::path db_path(vdb_dir);
  fs::path vdb_name(VDB_FILENAME);
  return (db_path / vdb_name);
}

Genotypes
CSearchApplication::SlurpGenotypes(const string vdb_dir) {
  Genotypes genotypes;
  fs::directory_iterator end_iter;
  fs::path vdb_path(vdb_dir);
  fs::path vt_path(VDB_FILENAME);
  if (fs::exists(vdb_path) && fs::is_directory(vdb_path)) {
    for(fs::directory_iterator dir_iter(vdb_path) ; dir_iter != end_iter ; ++dir_iter) {
      fs::path base = dir_iter->path().stem();
      if (base.string().find(VDB_BASE_FILENAME) == std::string::npos) {
        boost::shared_ptr<Bitmap> bitmap(new Bitmap());
        ifstream in;
        in.open(dir_iter->path().string().c_str());
        bitmap->read(in);
        in.close();
        genotypes.insert(make_pair(base.string(), bitmap));
      }
    }
  }
  return genotypes;
}

CSearchResultSet
CSearchApplication::RunBlast(string dbname, TSeqLocVector query_loc, CRef<CBlastOptionsHandle> opts)
{
  const CSearchDatabase target_db(dbname, CSearchDatabase::eBlastDbIsNucleotide);
  CRef<IQueryFactory> query_factory(new CObjMgr_QueryFactory(query_loc));
  CLocalBlast blaster(query_factory, opts, target_db);
  return *blaster.Run();
}

void
CSearchApplication::PrintQueryResult(CSearchResults &queryResult)
{
  const list<CRef<CSeq_align> > &seqAlignList = queryResult.GetSeqAlign()->Get();
  ITERATE(list <CRef<CSeq_align> >, seqAlignIter, seqAlignList) {
    cout << "start-stop (target): id "
      << (*seqAlignIter)->GetSeq_id(1).GetSeqIdString()
      << " " << (*seqAlignIter)->GetSeqStart(1)
      << "-" << (*seqAlignIter)->GetSeqStop(1) << endl;
  }
}

void
CSearchApplication::PrintErrorMessages(CSearchResults &queryResult) {
    // print any error messages
    TQueryMessages messages = queryResult.GetErrors(eBlastSevWarning);
    if (messages.size() > 0) {
      CConstRef<CSeq_id> seq_id = queryResult.GetSeqId();
      if (seq_id.NotEmpty()) {
        cerr << "ID: " << seq_id->AsFastaString() << endl;
      } else {
        cerr << "ID: " << "Unknown" << endl;
      }
      ITERATE(vector<CRef<CSearchMessage> >, it, messages) {
        cerr << (*it)->GetMessage() << endl;
      }
    }
}

void
CSearchApplication::Exit(void)
{
  SetDiagStream(0);
}

int
main(int argc, const char* argv[])
{
  return CSearchApplication().AppMain(argc, argv, 0, eDS_Default, 0);
}
