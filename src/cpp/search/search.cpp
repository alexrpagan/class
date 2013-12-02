#include <ncbi_pch.hpp>
#include <corelib/ncbiapp.hpp>
#include <corelib/ncbienv.hpp>
#include <corelib/ncbiargs.hpp>

USING_NCBI_SCOPE;

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

  SetupArgDescriptions(arg_desc.release());
}

int CSearchApplication::Run(void)
{
  // Get arguments
  const CArgs& args = GetArgs();

  return 0;
}

void CSearchApplication::Exit(void)
{
  // TODO: cleanup goes here.
}

int main(int argc, const char* argv[])
{
  return CSearchApplication().AppMain(argc, argv);
}
