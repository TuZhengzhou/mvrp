#include <vector>
#include <ostream>
#include <iterator>
#include <sstream>
#include <openssl/sha.h>
#include "structs.hpp"

using cred::G1;
using cred::G2;
using cred::GT;
using cred::Fr;

/************************** CCredSRS ***************************/
namespace cred {

AgendaItem::AgendaItem(const string& label) {
  this->label = label;
  this->enter_time      = libff::get_nsec_time();
  this->enter_cpu_time  = libff::get_nsec_cpu_time();
  this->leave_time        = 0;
  this->leave_cpu_time    = 0;

  this->times_enter = 1;
  this->times_leave = 0;
};

void AgendaItem::mark_enter() {
  this->enter_time     += libff::get_nsec_time();
  this->enter_cpu_time += libff::get_nsec_cpu_time();

  this->times_enter += 1;
}

void AgendaItem::mark_leave() {
  this->leave_time     += libff::get_nsec_time();
  this->leave_cpu_time += libff::get_nsec_cpu_time();

  this->times_leave += 1;
}

void AgendaItem::print(long long& start_time, long long& start_cpu_time) {
  printf("%-35s\t", this->label.c_str());
  assert(this->times_enter == this->times_leave);
  
  // long long time_from_start = leave_time - start_time;
  long long time_from_last = (leave_time - enter_time) / this->times_enter;

  // long long cpu_time_from_start = leave_cpu_time - start_cpu_time;
  long long cpu_time_from_last = (leave_cpu_time - enter_cpu_time) / this->times_enter;

  if (time_from_last != 0) {
      double parallelism_from_last = 1.0 * cpu_time_from_last / time_from_last;
      printf("[%0.8fs x%0.2f]\n", time_from_last * 1e-9, parallelism_from_last);
  } else {
      printf("[             ]\n");
  }
}

Agenda::Agenda(const string mode, const size_t N, const size_t D, const size_t num, const size_t bit_len) {
  start_time = libff::get_nsec_time();
  start_cpu_time = libff::get_nsec_cpu_time();

  this->mode = mode;
  this->N = N;
  this->D = D;
  this->num = num;
  this->bit_len = bit_len;
  items.clear();
};

void Agenda::create_item(const string& label) {
  if(this->items.find(label) == this->items.end()) {
    label_records.push_back(label);

    // long long now = libff::get_nsec_time();
    AgendaItem item = AgendaItem(label);
    items[label] = item;
  } else {
    items[label].mark_enter();
  }


}

void Agenda::mark_item_end(const string& label) {
  assert(this->items.find(label) != this->items.end());

  items[label].mark_leave();
}

string Agenda::to_string_rough() {
  char ret[256];
  if(mode == "base") {
    sprintf(ret, "%2lu bit \\times\\thinspace %4lu & %4lu & %0.1f & %0.1f\n", bit_len, num, proof_size, items["Prove_base"].duration() * 1000, items["Verify_base"].duration() * 1000);
  } else {
    sprintf(ret, "%2lu bit \\times\\thinspace %4lu & %4lu & %0.1f & %0.1f\n", bit_len, num, proof_size, items["Prove_improved"].duration() * 1000, items["Verify_improved"].duration() * 1000);
  }
  return string(ret);
}

string Agenda::ipa_record() {
  char ret[256];
  assert(mode == "base");
  sprintf(ret, "%2lu bit \\times\\thinspace %4lu & %4lu & %0.1f & %0.1f\n", bit_len, num, ipa_size, items["Prove::IpaProve"].duration() * 1000, items["Verify_base::IpaVerify"].duration() * 1000);
  return string(ret);
}

string Agenda::bullet_record() {
  char ret[256];
  assert(mode == "base");
  float bullet_prove_t = 1000*(items["bullet::Round1"].duration() + items["bullet::Round2_plus_pi_tilde"].duration() \
                               - items["Prove::_pi_tilde"].duration() + items["bullet::Round3"].duration());
  sprintf(ret, "%2lu bit \\times\\thinspace %4lu & %4lu & %0.1f & %0.1f\n", bit_len, num, ipa_size, bullet_prove_t, items["Verify_base::IpaVerify"].duration() * 1000);
  return string(ret);
}

void Agenda::print() {
  printf("Mode: %s\n", mode.c_str());
  printf("Args: N_%lu, D_%lu, num_%lu, bitlen_%lu\n", N, D, num, bit_len);
  printf("ProofSize: %lu\n", proof_size);
  for(auto label : label_records) {
    items[label].print(this->start_time, this->start_cpu_time);
  }

  printf("\n\n");
  for(auto item: v_sizes) {
    printf("* Peak vsize (physical memory+swap) in mebibytes (%s): %lu\n", item.first.c_str(), item.second >> 20);
  }
}

void Agenda::write_file(const char* filename) {
  FILE* fp = freopen(filename,"w", stdout);
  if(fp == NULL) cout<<"error"<<endl;

  this->print();
  fclose(fp);

  fp = freopen("/dev/tty","w",stdout); // 恢复默认输出
  // fclose(fp);
}


};
