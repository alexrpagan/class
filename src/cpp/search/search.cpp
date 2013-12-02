#include <cstring>

#include <ncbi_pch.hpp>
#include <corelib/ncbiapp.hpp>
#include <corelib/ncbienv.hpp>
#include <corelib/ncbiargs.hpp>

#include <objmgr/object_manager.hpp>

#include <objects/seqalign/Seq_align_set.hpp>

#include <algo/blast/api/sseqloc.hpp>
#include <algo/blast/api/local_blast.hpp>
#include <algo/blast/api/bl2seq.hpp> // added
#include <algo/blast/api/uniform_search.hpp>
#include <algo/blast/api/blast_types.hpp>
#include <algo/blast/api/blast_aux.hpp>
#include <algo/blast/api/objmgr_query_data.hpp>
#include <algo/blast/api/blast_options_handle.hpp>
#include <algo/blast/api/blast_nucl_options.hpp>
#include <algo/blast/api/blast_prot_options.hpp>

#include <algo/blast/blastinput/blast_input.hpp>
#include <algo/blast/blastinput/blast_fasta_input.hpp>

#ifndef SKIP_DOXYGEN_PROCESSING
USING_NCBI_SCOPE;
USING_SCOPE(blast);
USING_SCOPE(objects);
#endif /* SKIP_DOXYGEN_PROCESSING */

namespace {

  const int SUCCESS = 0;
  const int ERROR = 1;

  const bool IS_PROTEIN = false;

  const string BLAST_PROGRAM = "blastn";
  const string INPUT_FLAG    = "in";
  const string EVALUE_FLAG   = "evalue";
  const string DB_FLAG       = "db";
}

class CSearchApplication : public CNcbiApplication
{
private:
  virtual void Init(void);
  virtual int  Run(void);
  virtual void Exit(void);
};

void CSearchApplication::Init(void)
{
  auto_ptr<CArgDescriptions> arg_desc(new CArgDescriptions);

  arg_desc->SetUsageContext(GetArguments().GetProgramBasename(),
                "Parallel blast over compressed genomes");

  arg_desc->AddKey
    (DB_FLAG, "Database",
     "Location of reference blast database",
     CArgDescriptions::eString);

  arg_desc->AddKey(INPUT_FLAG, "QueryFile",
    "FASTA file containing queries", CArgDescriptions::eInputFile);

  arg_desc->AddDefaultKey(EVALUE_FLAG, "evalue",
    "E-value threshold for final hits", CArgDescriptions::eDouble, "1e-30");

  // arg_desc->AddDefaultKey("coarse_evalue", "coarse_evalue",
  // "E-value threshold for coarse search against uniques database", CArgDescriptions::eDouble, "1e-20");

  SetupArgDescriptions(arg_desc.release());

}

int CSearchApplication::Run(void)
{

  const CArgs& args = GetArgs();

  EProgram program = ProgramNameToEnum(BLAST_PROGRAM);

  CRef<CBlastOptionsHandle> opts(CBlastOptionsFactory::Create(program));

  if(args[EVALUE_FLAG].AsDouble()) {
    opts->SetEvalueThreshold(args[EVALUE_FLAG].AsDouble());
  }

  opts->Validate();

  CRef<CObjectManager> objmgr = CObjectManager::GetInstance();
  if (!objmgr) {
    throw std::runtime_error("Could not initialize object manager");
  }

  SDataLoaderConfig dlconfig(IS_PROTEIN);
  CBlastInputSourceConfig iconfig(dlconfig);
  CBlastFastaInputSource fasta_input(args[INPUT_FLAG].AsInputFile(), iconfig);
  CScope scope(*objmgr);

  CBlastInput blast_input(&fasta_input);
  TSeqLocVector query_loc = blast_input.GetAllSeqLocs(scope);

  const CSearchDatabase target_db(args[DB_FLAG].AsString(), CSearchDatabase::eBlastDbIsNucleotide);
  CRef<IQueryFactory> query_factory(new CObjMgr_QueryFactory(query_loc));
  CLocalBlast blaster(query_factory, opts, target_db);
  CSearchResultSet results = *blaster.Run();

  int num_results = 0;
  for (int i = 0; i < results.GetNumResults(); i++) {
    const list < CRef <CSeq_align> > &seqAlignList = results[i].GetSeqAlign()->Get();
    num_results += seqAlignList.size();
  }
  cout << "num_results: " << num_results << endl;
  cout << endl;

  return SUCCESS;

}

void CSearchApplication::Exit(void)
{
  SetDiagStream(0);
}

#ifndef SKIP_DOXYGEN_PROCESSING
int main(int argc, const char* argv[])
{
  return CSearchApplication().AppMain(argc, argv, 0, eDS_Default, 0);
}
#endif /* SKIP_DOXYGEN_PROCESSING */
