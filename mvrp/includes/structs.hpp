#ifndef CRED_STRUCTS_HPP
#define CRED_STRUCTS_HPP
#include <vector>
#include <map>
#include <tuple>
#include "libff/common/profiling.hpp"
#include "libff/algebra/curves/public_params.hpp"
#include "libff/common/default_types/ec_pp.hpp"
#include "basic_types.hpp"
#include "math_operations.hpp"
#include "ipa.hpp"
#include "kzg.hpp"
#include "range.hpp"
#include "circuit.hpp"
#include <xassert/XAssert.h>

#ifndef NO_PROCPS
#include <proc/readproc.h>
#endif

using namespace std;

namespace cred {

class AgendaItem {
  string label;
  long long enter_time;
  long long enter_cpu_time;
  long long leave_time;
  long long leave_cpu_time;

  size_t times_enter;
  size_t times_leave;
public:
  AgendaItem() {};
  AgendaItem(const string& label);
  void mark_enter();
  void mark_leave();
  void print(long long& start_time, long long& start_cpu_time);
  inline double duration() {
#ifdef CRED_DEBUG
    assert(times_enter == times_leave);
#endif
    long long time_from_last = (leave_time - enter_time) / this->times_enter;
    return time_from_last * 1e-9;
  }
};

class Agenda {
  string mode;  // base or improved
  size_t N;     // max bits
  size_t D;     // total bits
  size_t num;   // 
  size_t bit_len;

  size_t ipa_size;
  size_t proof_size;

  std::vector<string> label_records;
  std::map<string, AgendaItem> items;
  std::map<string, size_t> v_sizes;

  long long start_time;
  long long start_cpu_time;
public:
  Agenda() {};
  Agenda(const string mode, const size_t N, const size_t D, const size_t num, const size_t bit_len);
  void create_item(const string& label);
  void mark_item_end(const string& label);
  void mark_mem_usage(const string& label) {
    struct proc_t usage;
    look_up_our_self(&usage);
    if(v_sizes.find(label) == v_sizes.end()) {
      v_sizes[label] = usage.vsize;
    } else {
      v_sizes[label] = v_sizes[label] > usage.vsize ? v_sizes[label] : usage.vsize;
    }
  }
  void mark_proof_size(const CRangeProof& pi) {
    if(mode == "base") {
      this->proof_size = pi.size_base();
    } else {
      this->proof_size = pi.size_improved();
    }
  }
  void mark_ipa_size(const CRangeProof& pi) {
    assert(mode == "base");
    this->ipa_size = pi._pi_ipa.size();
  }

  string to_string_rough();
  string ipa_record();
  string bullet_record();
  void print();
  void write_file(const char* filename);
};



// Agenda agenda;

};
#endif