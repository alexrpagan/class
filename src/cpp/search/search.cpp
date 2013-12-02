// TODO: remove unused INCLUDES

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

USING_NCBI_SCOPE;
USING_SCOPE(blast);
USING_SCOPE(objects);

namespace {

  const int SUCCESS = 0;
  const int ERROR = 1;

  const bool IS_PROTEIN = false;

  const string COARSE_EVALUE_FLAG = "coarse_evalue";
  const string BLAST_PROGRAM      = "blastn";
  const string EVALUE_FLAG        = "evalue";
  const string INPUT_FLAG         = "in";
  const string VDB_FLAG           = "vdb";
  const string DB_FLAG            = "db";

  const string VDB_FILENAME = "vt.db";

}

class CSearchApplication : public CNcbiApplication
{
private:
  virtual void Init(void);
  virtual int  Run(void);
  virtual void Exit(void);

  CSearchResultSet RunBlast(string dbname, TSeqLocVector query_loc, CRef<CBlastOptionsHandle> opts);
  void PrintErrorMessages(CSearchResults &queryResult);
  void PrintQueryResult(CSearchResults &queryResult);
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

  string dbname = args[DB_FLAG].AsString();

  // TODO: Pre-load bitvectors for chromosome of interest.

  CSearchResultSet results = RunBlast(dbname, query_loc, opts);
  for (unsigned int i = 0; i < results.GetNumResults(); i++) {
    CSearchResults &queryResult = results[i];
    PrintErrorMessages(queryResult);
    PrintQueryResult(queryResult);

    // TODO: search tabix DB for variants

  }

  return SUCCESS;
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
