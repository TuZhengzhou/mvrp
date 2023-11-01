#include <iostream>
#include <vector>
#include <map>
#include <assert.h>
#include <iterator>
#include "libff/common/profiling.hpp"
#include "libff/algebra/curves/mnt/mnt6/mnt6_pp.hpp"
#include "libff/algebra/curves/bn128/bn128_pp.hpp"
#include "libff/algebra/curves/public_params.hpp"
#include "structs.hpp"
#include "prover.hpp"
#include "verifier.hpp"
#include "ipa.hpp"
#include "kzg.hpp"
#include "circuit_generator.hpp"

// #include "prover.cpp"
// #include "verifier.cpp"

using namespace std;
using namespace cred;
using cred::G1;
using cred::G2;
using cred::GT;
using cred::Fr;

int main(int argc,char *argv[]) {
  size_t N_bit;
  size_t D_bit_max;
  size_t repeat;

  map<string, size_t> args;
  for (int i = 1; i < argc; i++) {
    char* token = strtok(argv[i], "="); // 使用'='分割参数
    if (token != nullptr) {
      std::string str = token;
      token = strtok(nullptr, "="); // 获取等号后的部分
      if (token != nullptr) {
        std::string num = token;
        // std::cout << "str: " << str << ", num: " << num << std::endl;
        args[str] = stoul(num);
      }
    }
  }
  N_bit      = (args.find("-N") != args.end()) ? min(args["-N"], (size_t)16)    : (size_t)8; // 设定值 : 默认值
  D_bit_max  = (args.find("-D") != args.end()) ? min(args["-D"], (size_t)N_bit) : (size_t)6;
  repeat     = (args.find("-R") != args.end()) ? min(args["-R"], (size_t)10)    : (size_t)2;

  libff::start_profiling();
  libff::default_ec_pp::init_public_params();

  // auto gen = CircuitParaGenerator(4, 2);
  // bool check = gen.check();

  // return 0;

  size_t N, D, num, bit_len;
  string mode;
  Agenda agenda;
  CredProof pi;
  bool result;

  N = (1 << N_bit) + 1;
  CredSRS srs = CredSRS(N);

  const string base_format = "/root/mvrp/our_cred/static/base/%s-N_%lu-D_%lu-num_%lu-bitlen_%lu.txt";
  const string improved_format = "/root/mvrp/our_cred/static/improved/%s-N_%lu-D_%lu-num_%lu-bitlen_%lu.txt";

  /** N = 1 << 20, D = 2 << 19, bit_len from 2 << 1 to 2 << 8, num from 2 << 18 to 2 << 11 */
  size_t D_bit, i, j;
  string base_record, improved_record, ipa_record;
  for(D_bit = 6; D_bit <= D_bit_max; D_bit++) {
    // N = (1 << (D_bit+D_bit));
    D = 1 << D_bit;
    char file_name[256];
  
    i = D_bit == 6 ? 3 : 6;
    for(; i <= 6; i++) {
      bit_len = 1 << i;
      num = D / bit_len;
      
      std::vector<Fr>     set_values_vec(num);
      std::vector<size_t> indexs(num);
      std::vector<Range>  range_vec(num);
      for(j = 0; j < num; j++) {
          set_values_vec[j] = Fr::one() + Fr::one();
          indexs[j] = j+1;
          range_vec[j] = Range(Fr::zero(), -Fr::one(), bit_len);
      }
      auto set      = SET(set_values_vec);
      auto ranges   = Ranges(indexs, range_vec);

      auto u = Fr::random_element() * G1::one();
      auto ipa_sys = IPAProveSystem(u);
      auto prover   = Prover(srs, set, ranges, ipa_sys);
      auto C        = prover._C;
      auto verifier = Verifier(srs, C, num, ranges, ipa_sys);

      mode = "base";
      sprintf(file_name, base_format.c_str(), mode.c_str(), N, D, num, bit_len);
      agenda = Agenda(mode, N, D, num, bit_len);
      for(j = 0; j < repeat; j++) {
        pi = prover.prove(agenda, false);
        result = verifier.verify(pi, agenda, false);
      }
      agenda.print();
      agenda.write_file(file_name);
      base_record += agenda.to_string_rough();
      ipa_record  += agenda.ipa_record();

      mode = "improved";
      sprintf(file_name, improved_format.c_str(), mode.c_str(), N, D, num, bit_len);
      agenda = Agenda(mode, N, D, num, bit_len);
      for(size_t i = 0; i < repeat; i++) {
        pi = prover.prove(agenda, true);
        result = verifier.verify(pi, agenda, true);
      }
      agenda.print();
      agenda.write_file(file_name);
      improved_record += agenda.to_string_rough();
    }
  }
  printf("%s\n", base_record.c_str());
  printf("%s\n", ipa_record.c_str());
  printf("%s\n", improved_record.c_str());

  string base_file = "/root/mvrp/our_cred/static/base.txt";
  string ipa_file = "/root/mvrp/our_cred/static/ipa.txt";
  string improved_file = "/root/mvrp/our_cred/static/improved.txt";
  write_file(base_file, base_record);
  write_file(ipa_file, ipa_record);
  write_file(improved_file, improved_record);

  std::cout << "output OK:>" << endl;
}

