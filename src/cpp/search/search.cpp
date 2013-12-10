// TODO: remove unused INCLUDES

#include <fstream>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <stdlib.h>
#include <sys/sysinfo.h>

#include "variant.cpp"

#include "timer.cpp"
#include "ewah.h"
#include "tabix.hpp"
#include "Fasta.h"

#include <ncbi_pch.hpp>
#include <corelib/ncbiapp.hpp>
#include <corelib/ncbistre.hpp>
#include <corelib/ncbienv.hpp>
#include <corelib/ncbiargs.hpp>
#include <corelib/ncbidbg.hpp>

#include <objmgr/object_manager.hpp>

#include <objects/seqalign/Seq_align_set.hpp>

#include <objtools/blast/seqdb_reader/seqdb.hpp>
#include <objtools/blast/seqdb_reader/seqdbexpert.hpp>
#include <objtools/blast/blastdb_format/seq_writer.hpp>
#include <objtools/blast/blastdb_format/blastdb_formatter.hpp>
#include <objtools/blast/blastdb_format/blastdb_seqid.hpp>

#include <algo/blast/api/sseqloc.hpp>
#include <algo/blast/api/local_blast.hpp>
#include <algo/blast/api/bl2seq.hpp>
#include <algo/blast/api/uniform_search.hpp>
#include <algo/blast/api/blast_types.hpp>
#include <algo/blast/api/blast_aux.hpp>
#include <algo/blast/api/objmgr_query_data.hpp>
#include <algo/blast/api/blast_options_handle.hpp>
#include <algo/blast/api/blast_nucl_options.hpp>

#include <algo/blast/blastinput/blast_input.hpp>
#include <algo/blast/blastinput/blast_fasta_input.hpp>

#include <boost/filesystem.hpp>
#include <boost/shared_ptr.hpp>

USING_NCBI_SCOPE;
USING_SCOPE(blast);
USING_SCOPE(objects);

namespace fs = boost::filesystem;

// TODO: move these to class definition?
namespace {

  const int SUCCESS = 0;
  const int ERROR = 1;

  const unsigned int PUNT_THRESH = 1000;

  const int FUZZ        = 1500;
  const int FRINGE      = 50;
  const int FASTA_WIDTH = 80;

  const bool IS_PROTEIN = false;

  const string COARSE_EVALUE_FLAG = "coarse_evalue";
  const string BLAST_PROGRAM      = "blastn";
  const string EVALUE_FLAG        = "evalue";
  const string INPUT_FLAG         = "in";
  const string FULL_FLAG          = "full";
  const string CHR_FLAG           = "chr";
  const string REF_FLAG           = "ref";
  const string VDB_FLAG           = "vdb";
  const string DB_FLAG            = "db";
  const string PERF_FLAG          = "p";
  const string VERBOSE_FLAG       = "v";

  const string VDB_BASE_FILENAME  = "vt";
  const string VDB_FILENAME       = "vt.db.gz";

}

inline bool isValidNuc(char c)
{
  switch(c) {
  case 'A':
  case 'C':
  case 'G':
  case 'T':
  case 'U':
  case 'R':
  case 'Y':
  case 'S':
  case 'W':
  case 'K':
  case 'M':
  case 'B':
  case 'D':
  case 'H':
  case 'V':
  case 'N':
    return true;
  default:
    return false;
  }
}

// TODO: move these typedefs into class def.
typedef EWAHBoolArray<uint32_t> Bitmap;
typedef map<string, boost::shared_ptr<Bitmap> > Genotypes;
typedef vector< CRef<CBlastDBSeqId> > TQueries;

typedef vector<string> Row;
typedef list<Row> Rows;

class SearchApp : public CNcbiApplication
{
private:

  string _chr_name;

  bool _verbose;

  int _curr_query;

  // perf
  bool _report_perf;
  Timer _timer;
  vector<double> _ref_seq_retrieval_times;
  vector<double> _variant_fetch_times;
  vector<double> _variant_filter_times;
  vector<double> _fine_blast_times;

  // readers
  FastaReference _fastaRef;

  // bitvectors
  Genotypes _genotypes;

  virtual void Init(void);
  virtual int  Run(void);
  virtual void Exit(void);

  CSearchResultSet RunBlast(string dbname, TSeqLocVector query_loc, CRef<CBlastOptionsHandle> opts);
  fs::path getVDBFilePath(const string vdb_dir);
  string getRegion(int start, int end);
  Bitmap bitmapFromArray(vector<size_t> array);
  void slurpGenotypes(const string vdb_dir);
  void PrintErrorMessages(CSearchResults &queryResult);
  void getVariantSequence(string &outstr, vector<Variant> &variants, int start, int end);
  void getReferenceSequence(string &outstr, int start, int end);
  void reportPerf(string label, vector<double>& samples);
  void printVariant(Row &variant);

  CSearchResultSet fullBlast(TSeqLocVector query_loc, CRef<CBlastOptionsHandle> opts, Uint8 *full_db_size);

};

// TODO: are nucs equivalent?

void
SearchApp::Init(void)
{

  auto_ptr<CArgDescriptions> arg_desc(new CArgDescriptions);

  arg_desc->SetUsageContext(GetArguments().GetProgramBasename(),
                "Blast over bitvector-compressed genomes");

  arg_desc->AddDefaultKey(FULL_FLAG, "Full",
    "Full blastdb", CArgDescriptions::eString, "");

  arg_desc->AddKey(DB_FLAG, "Database",
    "Location of reference blast database", CArgDescriptions::eString);

  arg_desc->AddKey(VDB_FLAG, "VDB",
    "Location of variant database", CArgDescriptions::eString);

  arg_desc->AddKey(CHR_FLAG, "chr",
    "The chromosome to blast over", CArgDescriptions::eString);

  arg_desc->AddKey(REF_FLAG, "reference",
    "Raw fasta of reference seqences. Must be fasta indexed.", CArgDescriptions::eString);

  arg_desc->AddKey(INPUT_FLAG, "QueryFile",
    "FASTA file containing queries", CArgDescriptions::eInputFile);

  arg_desc->AddDefaultKey(EVALUE_FLAG, "evalue",
    "E-value threshold for final hits", CArgDescriptions::eDouble, "1e-30");

  arg_desc->AddDefaultKey(COARSE_EVALUE_FLAG, "CoarseEvalue",
    "E-value threshold for coarse search against reference sequence", CArgDescriptions::eDouble, "1e-20");

  arg_desc->AddFlag(PERF_FLAG, "Report performance data", true);

  arg_desc->AddFlag(VERBOSE_FLAG, "Be loquacious.", true);

  SetupArgDescriptions(arg_desc.release());

}

int
SearchApp::Run(void)
{
  const CArgs& args = GetArgs();

  string dbname  = args[DB_FLAG].AsString();
  string vdbpath = args[VDB_FLAG].AsString();
  string refpath = args[REF_FLAG].AsString();

  double coarse_evalue = args[COARSE_EVALUE_FLAG].AsDouble();
  double fine_evalue   = args[EVALUE_FLAG].AsDouble();

  // init members.
  _chr_name    = args[CHR_FLAG].AsString();
  _report_perf = args[PERF_FLAG].AsBoolean();
  _verbose     = args[VERBOSE_FLAG].AsBoolean();

  _timer.update_time();
  slurpGenotypes(vdbpath);
  if( _report_perf ) {
    cerr << "Reading bitvectors: " << _timer.update_time() << endl ;
  }

  string vdb_file_name = getVDBFilePath(vdbpath).string();
  Tabix vdb(vdb_file_name);

  set<int> variants_to_ignore;

  _timer.update_time();
  vector<Variant> variantDB;
  string region = getRegion(0, std::numeric_limits<int>::max());
  vdb.setRegion(region);
  string line;
  bool first_variant = true;
  Variant last_var;
  Bitmap length_modifying_variants;
  while(vdb.getNextLine(line)) {
    Row row;
    stringstream ss(line);
    string cell;
    while(getline(ss, cell, '\t')) row.push_back(cell);
    Variant var(row);
    variantDB.push_back(var);
    if(first_variant) {
      first_variant = false;
      last_var = var;
    } else {
      if (var.GetPos() + 1 < last_var.GetPos() + last_var.GetRef().size()) {
        variants_to_ignore.insert(var.GetBit());
      } else {
        if (var.GetLengthMod() != 0) {
          length_modifying_variants.set(var.GetBit());
        }
        last_var = var;
      }
    }
  }
  if( _report_perf ) {
    cerr << "Reading variant DB: " << _timer.update_time() << endl ;
  }

  if (_verbose) {
    cerr << "Skipping " << variants_to_ignore.size() << " variants "<< endl;
  }

  _fastaRef.open(refpath);

  // Prep coarse blast.
  EProgram program = ProgramNameToEnum(BLAST_PROGRAM);
  CRef<CBlastOptionsHandle> opts(CBlastOptionsFactory::Create(program));

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

  bool compare_with_full = args[FULL_FLAG].AsString().size() > 0;

  Uint8 full_db_size;
  CSearchResultSet full_results;
  if ( compare_with_full ) {
    opts->SetEvalueThreshold(fine_evalue);
    opts->Validate();
    full_results = fullBlast(query_loc, opts, &full_db_size);
  }

  int full_ctr = 0, missed_ctr = 0, coarse_hits = 0;

  // set evalue-back to coarse
  opts->SetEvalueThreshold(coarse_evalue);
  opts->SetDbLength(full_db_size);
  opts->Validate();
  // coarse blast.
  CSearchResultSet results = RunBlast(dbname, query_loc, opts);

  for (unsigned int query_idx = 0; query_idx < results.GetNumResults(); query_idx++) {
    _curr_query = query_idx;
    if ( _verbose ) {
      cout << "# Query: " << query_idx << endl;
      cerr << "# Query: " << query_idx << endl;
    }

    CSearchResults &queryResult = results[query_idx];
    PrintErrorMessages(queryResult);

    const list<CRef<CSeq_align> > &seqAlignList = queryResult.GetSeqAlign()->Get();

    coarse_hits = seqAlignList.size();

    if ( seqAlignList.size() > PUNT_THRESH ) {
      cerr << "Too many hits. Skipping! " << query_idx << " (" << seqAlignList.size() << ")" << endl;
      continue;
    }

    set<pair<int,int> > class_hits;
    map<string, vector<pair<int, int> > > indiv_to_hits;

    ITERATE(list <CRef<CSeq_align> >, seqAlignIter, seqAlignList) {
      int start = (*seqAlignIter)->GetSeqStart(1);
      int stop  = (*seqAlignIter)->GetSeqStop(1);

      // TODO: use fasta index to get sequence length and to do bounds check here.
      int adj_start = start - FRINGE;
      int adj_stop  = stop  + FRINGE;
      start = adj_start;
      stop  = adj_stop;

      vector<Variant> variants;
      string line;
      // NB: variant DB positions are 1-indexed.
      string region = getRegion(start + 1, stop);
      _timer.update_time();
      vdb.setRegion(region);
      while(vdb.getNextLine(line)) {
        Row row;
        stringstream ss(line);
        string cell;
        while(getline(ss, cell, '\t')) {
          row.push_back(cell);
        }
        Variant var(row);
        variants.push_back(var);
      }
      _variant_fetch_times.push_back(_timer.update_time());

      map<size_t, Variant> variants_by_bit;
      Bitmap variant_filter;
      for (vector<Variant>::iterator it = variants.begin(); it != variants.end(); ++it) {
        variant_filter.set(it->GetBit());
        variants_by_bit.insert(make_pair(it->GetBit(), *it));
      }

      _timer.update_time();
      map<vector<size_t>, boost::shared_ptr<vector<string> > > variant_genotypes;
      for (Genotypes::iterator it = _genotypes.begin(); it != _genotypes.end(); ++it) {
        Bitmap and_result;
        it->second->logicaland(variant_filter, and_result);
        // Convert to vector of bits to key variant_genotypes.
        // As the bitset class does not implement cmp.
        vector<size_t> bits_set = and_result.toArray();
        if (variant_genotypes.find(bits_set) == variant_genotypes.end()) {
          boost::shared_ptr<vector<string> > vec(new vector<string>());
          variant_genotypes.insert(make_pair(bits_set, vec));
        }
        variant_genotypes[bits_set]->push_back(it->first);
      }
      _variant_filter_times.push_back(_timer.update_time());

      map<vector<size_t>, boost::shared_ptr<vector<string> > >::iterator it = variant_genotypes.begin();

      vector<pair<vector<string>, string> > fine_blast_targets;
      for(; it != variant_genotypes.end(); ++it) {
        vector<Variant> variants;
        for(unsigned int variant_idx = 0; variant_idx < (it->first).size(); ++variant_idx) {
          variants.push_back(variants_by_bit[(it->first)[variant_idx]]);
        }
        string var_seq;
        getVariantSequence(
            var_seq
          , variants
          , adj_start
          , adj_stop
        );
        fine_blast_targets.push_back(make_pair(*(it->second), var_seq));
      }

      // TODO: experiment with batching for fine blast!
      stringstream fine_blast_input_ss;
      for (size_t seq_idx = 0; seq_idx < fine_blast_targets.size(); ++seq_idx) {
        fine_blast_input_ss << ">SEQ" << seq_idx << endl;
        fine_blast_input_ss << fine_blast_targets[seq_idx].second << endl;
      }

      string fine_blast_input_str = fine_blast_input_ss.str();

      //trim trailing endl to make fasta parser happy.
      fine_blast_input_str.erase(fine_blast_input_str.length() - 1, 1);

      opts->SetEvalueThreshold(fine_evalue);
      opts->Validate();

      CRef<CObjectManager> fine_blast_objmgr = CObjectManager::GetInstance();
      if (!fine_blast_objmgr) {
        throw std::runtime_error("Could not initialize fine blast object manager");
      }

      CBlastFastaInputSource fine_blast_input(fine_blast_input_str, iconfig);
      CBlastInput input_from_str(&fine_blast_input);
      CScope fine_blast_scope(*fine_blast_objmgr);
      TSeqLocVector target_from_str = input_from_str.GetAllSeqLocs(fine_blast_scope);

      CBl2Seq fine_blaster(query_loc[query_idx], target_from_str, *opts);
      TSeqAlignVector fine_results(fine_blaster.Run());
      _fine_blast_times.push_back(_timer.update_time());

      for (size_t fine_hit_idx = 0; fine_hit_idx < fine_results.size(); ++fine_hit_idx) {
        CConstRef<CSeq_align_set> fine_align_set = fine_results[fine_hit_idx];
        const list <CRef<CSeq_align> > &fine_align_list = fine_align_set->Get();
        ITERATE(list<CRef<CSeq_align> >, fine_iter, fine_align_list) {

          int sequence_idx = atoi((*fine_iter)->GetSeq_id(1).GetSeqIdString().c_str()) - 1;
          int fine_start   = (*fine_iter)->GetSeqStart(1);
          int fine_stop    = (*fine_iter)->GetSeqStop(1);

          pair<int, int> hit = make_pair(adj_start + fine_start, adj_start + fine_stop);
          vector<string> affected_indivs = fine_blast_targets[sequence_idx].first;
          vector<string>::iterator indiv_it = affected_indivs.begin();
          for ( ; indiv_it != affected_indivs.end(); ++indiv_it ) {
            if (indiv_to_hits.find(*indiv_it) == indiv_to_hits.end()) {
              vector<pair<int, int> > indiv_hits;
              indiv_to_hits.insert(make_pair(*indiv_it, indiv_hits));
            }
            indiv_to_hits[*indiv_it].push_back(hit);
          }

          if ( _verbose ) {
            cerr << "Hit before translation" << endl;
            for (int i = 0; i < fine_blast_targets[sequence_idx].first.size(); i++) {
              cerr << fine_blast_targets[sequence_idx].first[i] << ", ";
            }
            cerr << endl;
            //cout << fine_blast_targets[sequence_idx].second << endl;
            //cout << fine_start << "-" << fine_stop << endl;
            cerr << adj_start + fine_start << "-" << adj_start + fine_stop << endl;
          }
        }
      }
    }

    map<string, vector<pair<int, int> > >::iterator hit_trans_it = indiv_to_hits.begin();
    for (; hit_trans_it != indiv_to_hits.end(); ++hit_trans_it) {

      string indiv = hit_trans_it->first;
      vector<pair<int, int> > hits = hit_trans_it->second;
      sort(hits.begin(), hits.end());

      int num_hits = hits.size();
      int reported = 0;

      pair<int, int> first_hit = hits.front();
      pair<int, int> last_hit  = hits.back();
      assert(last_hit.first >= first_hit.first);

      // Only loop over length-modifying variants
      Bitmap and_res;
      boost::shared_ptr<Bitmap> bits = _genotypes[indiv];
      bits->logicaland(length_modifying_variants, and_res);
      vector<size_t> set_bits = and_res.toArray();

      if (_verbose) {
        cerr << "Looping over " << set_bits.size()
             << " length modifying variants for " << indiv << endl;
      }

      int offset = 0;
      vector<pair<int, int> >::iterator hit_iter = hits.begin();
      vector<size_t>::iterator it = set_bits.begin();
      for (; it != set_bits.end(); ++it) {
        // get the next set bit
        int bit = *it;
        assert(!variants_to_ignore.count(bit));
        Variant var = variantDB[bit];
        assert(var.GetBit() == bit);
        assert(var.GetLengthMod() != 0);
        // loop through any variants we have passed, reporting hits with the offset up to this point
        for(; hit_iter != hits.end() && (*hit_iter).first < var.GetPos(); ++hit_iter) {
          int start_with_offset = offset + hit_iter->first;
          int stop_with_offset  = offset + hit_iter->second;
          cout << indiv
               << "\t" << start_with_offset
               << "\t" << stop_with_offset << endl;
          if ( compare_with_full ) {
            class_hits.insert(make_pair(start_with_offset - FUZZ, stop_with_offset + FUZZ));
          }
          reported++;
        }
        if (hit_iter == hits.end()) {
          break;
        }
        offset += var.GetLengthMod();
      }

      // individual has no more variants -- just print the rest of the hits with the current offset.
      for (; hit_iter != hits.end(); ++hit_iter) {
        int start_with_offset = offset + hit_iter->first;
        int stop_with_offset  = offset + hit_iter->second;
        cout << indiv
             << "\t" << start_with_offset
             << "\t" << stop_with_offset << endl;
        if ( compare_with_full ) {
          class_hits.insert(make_pair(start_with_offset - FUZZ, stop_with_offset + FUZZ));
        }
        reported++;
      }
      assert(num_hits == reported);
    }

    /*
      Iterate through all full hits for this query.
    */
    if ( compare_with_full ) {
      int full_query = 0, missed_query = 0;
      const list <CRef<CSeq_align> > &full_seqAlignList = full_results[query_idx].GetSeqAlign()->Get();
      ITERATE(list < CRef <CSeq_align> >, seqAlign_it, full_seqAlignList) {
        full_ctr++, full_query++;
        int start = (*seqAlign_it)->GetSeqStart(1);
        int end   = (*seqAlign_it)->GetSeqStop(1);
        if (!class_hits.count(make_pair(start, end))) {
          bool found_imperfect = false;
          if (_verbose) {
            // look for hits that aren't quite right
            cerr << "full blast: hit at pos " << start << "-" << end << " ";
            cerr << "Missed. ";
            double evalue;
            (*seqAlign_it)->GetNamedScore(CSeq_align::eScore_EValue, evalue);
            cerr << evalue << endl;
          }
          for (set<pair<int, int> >::iterator it = class_hits.begin(); it != class_hits.end(); it++) {
            int hit_start = it->first;
            int hit_end   = it->second;
            if (hit_start <= end && hit_end >= start) {
              if (((hit_end - FUZZ) - (hit_start + FUZZ)) == end - start) {
                found_imperfect = true;
              }
              break;
            }
          }
          if (!found_imperfect) {
            missed_ctr++, missed_query++;
          }
        }
      }
      if (full_query > 0) {
        double query_accuracy = 1 - missed_query / (double) full_query;
        if (query_accuracy < .8) {
          if(_verbose) {
            cerr << "Bad query: " << query_accuracy << ": " << query_idx
                 << ", " << full_query << " full hits. "
                 << coarse_hits << " coarse hits." << endl;
          }
        }
      }
    }
  }

  if ( _report_perf ) {
    reportPerf("Retrieving reference subseqence",     _ref_seq_retrieval_times);
    reportPerf("Fetching overlapping variants",       _variant_fetch_times);
    reportPerf("Bitwise-and over variant BVs",        _variant_filter_times);
    reportPerf("Initializing and running fine blast", _fine_blast_times);
    if ( compare_with_full ) {
      cerr << endl;
      cerr << "Overall accuracy: " << 1-missed_ctr/(double) full_ctr << endl;
      cerr << "Missed " << missed_ctr << "/" << full_ctr << endl;
    }
  }

  return SUCCESS;
}

// bool pairCompareByStartPos(const std::pair<int, int>& firstElem, const std::pair<int, int>& secondElem) {
//   return firstElem.first < secondElem.first;
// }


CSearchResultSet
SearchApp::fullBlast(
    TSeqLocVector query_loc
    , CRef<CBlastOptionsHandle> opts
    , Uint8 *full_db_size
) {
  const CArgs& args = GetArgs();
  const CSearchDatabase target_db(args[FULL_FLAG].AsString(), CSearchDatabase::eBlastDbIsNucleotide);

  *full_db_size = target_db.GetSeqDb()->GetTotalLength();
  cerr << "Full db size: " << *full_db_size << endl;

  CRef<IQueryFactory> query_factory(new CObjMgr_QueryFactory(query_loc));

  _timer.update_time();
  CLocalBlast blaster(query_factory, opts, target_db);
  cerr << "Time spent initializing full BLAST: " << _timer.update_time() << endl;

  _timer.update_time();
  CSearchResultSet results = *blaster.Run();
  cerr << "Full BLAST search time: " << _timer.update_time() << endl;

  return results;
}

void
SearchApp::reportPerf(string label, vector<double> &samples) {
  double sum = accumulate(samples.begin(), samples.end(), 0.0);
  double mean = sum / samples.size();
  double sq_sum = inner_product(samples.begin(), samples.end(), samples.begin(), 0.0);
  double stdev  = sqrt(sq_sum / samples.size() - mean * mean);
  cerr << label << ":" << endl;
  cerr << "samples: "  << samples.size() << endl;
  cerr << "sum:     "  << sum    << endl;
  cerr << "mean:    "  << mean   << endl;
  cerr << "stdev:   "  << stdev  << endl;
  cerr << endl;
}

void
SearchApp::getVariantSequence(
    string &outbuf
  , vector<Variant> &variants
  , int ref_start_pos
  , int ref_end_pos
) {

  // positions of last variant
  int last_start = 0, last_end = ref_start_pos;
  bool first = true;

  Variant last_variant;
  for (vector<Variant>::iterator it = variants.begin(); it != variants.end(); ++it) {
    Variant variant = *it;
    if (first) {
      last_variant = variant;
      last_start   = variant.GetPos();
      first = false;
    }
    if (variant.GetPos() < last_end) {
      if (_verbose) {
        // cerr << "Unusual overlapping or out-of-order variants:" << endl;
        // cerr << variant.GetPos() << " - " << last_end << endl;
        // for (vector<Variant>::iterator inner_it = variants.begin(); inner_it != variants.end(); ++inner_it) {
        //   cerr << *inner_it << endl;
        // }
      }
      continue;
    }
    // get all of the reference sequence up to this variant.
    string ref_before_var;
    if (variant.GetPos() - last_end > 0) {
      getReferenceSequence(ref_before_var, last_end, variant.GetPos());
      int ref_size = ref_before_var.size();
      // sanity check
      assert(ref_before_var[ref_size-1] == variant.GetRef()[0]);
      ref_before_var.resize(ref_size-1);
    }
    outbuf.append(ref_before_var);
    string alt = variant.GetAlt();
    if (alt != "<DEL>") {
      outbuf.append(variant.GetAlt());
    } else {
      if (_verbose) {
        // cerr << "Impreci se SV; ignoring." << endl;
      }
      outbuf.append(variant.GetRef());
    }
    last_start = variant.GetPos();
    last_end   = last_start + variant.GetRef().size();
    last_variant = variant;
  }
  if (last_end < ref_end_pos) {
    string endref;
    getReferenceSequence(endref, last_end, ref_end_pos);
    outbuf += endref;
  }
}


void
SearchApp::getReferenceSequence(string& outstr, int start, int end){
  _timer.update_time();
  const int seq_len = end - start + 1;
  string seq = _fastaRef.getSubSequence(_chr_name, start, seq_len);
  assert (seq.length() == seq_len);
  outstr = seq;
  _ref_seq_retrieval_times.push_back(_timer.update_time());
}

Bitmap
SearchApp::bitmapFromArray(vector<size_t> array) {
  Bitmap bits;
  for (unsigned int i = 0; i < array.size(); ++i) {
    bits.set(array[i]);
  }
  return bits;
}

string
SearchApp::getRegion(int start, int end) {
  ostringstream ss;
  ss << _chr_name << ":" << start << "-" << end;
  return ss.str();
}

fs::path
SearchApp::getVDBFilePath(const string vdb_dir) {
  fs::path db_path(vdb_dir);
  fs::path vdb_name(VDB_FILENAME);
  return (db_path / vdb_name);
}

void
SearchApp::slurpGenotypes(const string vdb_dir) {
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
        _genotypes.insert(make_pair(base.string(), bitmap));
      }
    }
  }
}

CSearchResultSet
SearchApp::RunBlast(string dbname, TSeqLocVector query_loc, CRef<CBlastOptionsHandle> opts)
{
  const CSearchDatabase target_db(dbname, CSearchDatabase::eBlastDbIsNucleotide);
  CRef<IQueryFactory> query_factory(new CObjMgr_QueryFactory(query_loc));

  _timer.update_time();
  CLocalBlast blaster(query_factory, opts, target_db);
  if (_report_perf) {
    cerr << "Coarse BLAST set-up: " << _timer.update_time() << endl;
  }

  _timer.update_time();
  CSearchResultSet res = *blaster.Run();
  if (_report_perf) {
    cerr << "Coarse BLAST search: " << _timer.update_time() << endl;
  }

  return res;
}

void
SearchApp::PrintErrorMessages(CSearchResults &queryResult) {
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
SearchApp::Exit(void)
{
  SetDiagStream(0);
}

int
main(int argc, const char* argv[])
{
  return SearchApp().AppMain(argc, argv, 0, eDS_Default, 0);
}
