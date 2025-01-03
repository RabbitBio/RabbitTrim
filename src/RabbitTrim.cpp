#include <sys/sysinfo.h>
#include <string>
#include <vector>
#include "CLI11.hpp"
#include "param.h"
#include "handler/HandlerPE.h"
#include "handler/HandlerSE.h"

// void __attribute__ ((constructor)) init_cpu_flag
// {
//   __builtin_cpu_init();
// }
// // define some macros related to cpuid
// #if __builtin_cpu_supports("sse2")
// #define TRIM_SSE2_ENABLED 1
// #if __builtin_cpu_supports("avx2")
// #define TRIM_AVX2_ENABLED 1
// #if __builtin_cpu_supports("avx512bw")
// #define TRIM_AVX512BW_ENABLED 1

int main( int argc, char **argv) {

  // process the command line parameters
  CLI::App app{"RabbitTrim : Optimised versions of Trimmomatic and Ktrim tools"};
  app.require_subcommand(1);
  CLI::App* ktrim = app.add_subcommand("ktrim", "Based Ktrim");
  CLI::App* trimmomatic = app.add_subcommand("trimmomatic", "Based Trimmomatic");
  // std::string mode;
  bool isPE = false; 
  bool isSE = false; 
  int threads = 1;
  int phred= 0;
  std::vector<std::string> forwardFiles;
  std::vector<std::string> reverseFiles;
  std::string output;
  std::string trimLog;
  std::string stats;
  std::vector<std::string> steps;
  bool quiet = false;
  bool validatePairing = false;
  int compressLevel = 4;
  // bool useIgzip = false;
  bool usePigz = true;
  int pigzThreadsNum = 2;


  // auto trimm_option_m = trimmomatic->add_option("-m,--mode", mode, "specify mode of data handled: single-end(SE) or pair-end data(PE)");
  auto trimm_flag_pe = trimmomatic->add_flag("-P,--PE,!--no-PE", isPE, "specify whether to be pair end data ");
  auto trimm_flag_se = trimmomatic->add_flag("-S,--SE,!--no-SE", isSE, "specify whether to be single end data");
  auto trimm_option_t = trimmomatic->add_option("-t,--threads", threads, "specify the number of threads");
  auto trimm_option_p = trimmomatic->add_option("-p,--phred", phred, "specify the baseline of the phred score");
  auto trimm_option_f = trimmomatic->add_option("-f,--forward", forwardFiles, "specify the path to the forward read file"); 
  auto trimm_option_r = trimmomatic->add_option("-r,--reverse", reverseFiles, "specify the path to the reverse read file"); 
  auto trimm_option_o = trimmomatic->add_option("-o,--output", output, "specify the path to the output file"); 
  auto trimm_option_l = trimmomatic->add_option("--log", trimLog, "specify the path to the trim log file"); 
  auto trimm_option_stat = trimmomatic->add_option("--stats", stats, "specify the path to the trim statistical data file"); 
  auto trimm_option_s = trimmomatic->add_option("-s,--steps", steps, "specify the steps to be performed");
  auto trimm_option_c = trimmomatic->add_option("-c,--compressLevel", compressLevel, "specify the compression level for the output file");
  auto trimm_flag_q = trimmomatic->add_flag("--quiet,!--no-quiet", quiet, "specify whether to print program runtime information");
  auto trimm_flag_v = trimmomatic->add_flag("--validatePair,!--no-validatePair", validatePairing, "specify whether to validate pair data");
  // auto trimm_flag_igzip = trimmomatic->add_flag("--igzip,!--no-igzip", useIgzip, "specify whether to use igzip");
  // auto trimm_flag_pigz = trimmomatic->add_flag("--pigz,!--no-pigz", usePigz, "specify whether to use pigz");
  // auto trimm_option_pigz_th = trimmomatic->add_option("-g,--pigzThreadsNum", pigzThreadsNum, "specify the max thread number of pigz");


  // trimm_option_m->required()->check(CLI::IsMember({"PE", "SE"}));
  trimm_flag_pe->excludes(trimm_flag_se)->needs(trimm_option_f)->needs(trimm_option_r);
  trimm_flag_se->excludes(trimm_flag_pe)->excludes(trimm_option_r)-> needs(trimm_option_f);
  trimm_option_t->check(CLI::PositiveNumber);
  trimm_option_p->check(CLI::IsMember(std::set<int> {0, 33, 64}));
  trimm_option_f->check(CLI::ExistingPath);
  trimm_option_r->check(CLI::ExistingPath);
  trimm_option_o->required();
  // trimm_option_l->required();
  trimm_option_stat->required();
  trimm_option_s->required();
  trimm_option_c->check(CLI::Range(0, 9).description("compression level is limited to be between 0 and 9"));
  // trimm_option_pigz_th->needs(trimm_flag_pigz) -> check(CLI::PositiveNumber) -> check(CLI::Bound(2, 64));


  std::string seqA;
  std::string seqB;
  std::string seqKit;
  int minLen = 36;
  int minQual = 20;
  int window = 5;
  double mismatch = 0.125;

  auto kt_flag_pe  = ktrim->add_flag("-P,--PE,!--no-PE", isPE, "specify whether to be pair end data ");
  auto kt_flag_se  = ktrim->add_flag("-S,--SE,!--no-SE", isSE, "specify whether to be single end data");
  auto kt_option_f = ktrim->add_option("-f,--forward", forwardFiles, "specify the path to the forward read file"); 
  auto kt_option_r = ktrim->add_option("-r,--reverse", reverseFiles, "specify the path to the reverse read file"); 
  auto kt_option_o = ktrim->add_option("-o,--output", output, "specify the path to the output file"); 
  auto kt_option_t = ktrim->add_option("-t,--threads", threads, "specify the number of threads");
  auto kt_option_p = ktrim->add_option("-p,--phred", phred, "specify the baseline of the phred score");
  auto kt_option_m = ktrim->add_option("-m,--mismatch", mismatch, "Set the proportion of mismatches allowed during index and sequence comparison, Default: 0.125");
  auto kt_option_l = ktrim->add_option("-l,--minlen", minLen, "specify the minimum read size to be kept");
  auto kt_option_q = ktrim->add_option("-q,--quality", minQual, "specify the minimum quality score");
  auto kt_option_w = ktrim->add_option("-w,--window", window, "specify the window size for quality check");
  auto kt_option_a = ktrim->add_option("-a,--seqA", seqA, "specify the adapter sequence of read1");
  auto kt_option_b = ktrim->add_option("-b,--seqB", seqB, "specify the adapter sequence of read2");
  auto kt_option_k = ktrim->add_option("-k,--seqKit", seqKit, "specify the sequencing kit to use built-in adapters");
  // auto kt_option_s = ktrim->add_option("-s,--steps", steps, "specify the steps to be performed");
  // auto kt_option_stat = ktrim->add_option("--stats", stats, "specify the path to the trim statistical data file"); 
  auto kt_option_c = ktrim->add_option("-c,--compressLevel", compressLevel, "specify the compression level for the output file");
  // auto kt_flag_igzip = ktrim->add_flag("--igzip,!--no-igzip", useIgzip, "specify whether to use igzip");
  // auto kt_flag_pigz  = ktrim->add_flag("--pigz, !--no-pigz" , usePigz , "specify whether to use pigz"); // reserve for debug
  // auto kt_option_pigz_th = ktrim->add_option("-g,--pigzThreadsNum", pigzThreadsNum, "specify the max thread number of pigz"); // reserve for debug 

  kt_flag_pe->excludes(kt_flag_se)->needs(kt_option_f)->needs(kt_option_r);
  kt_flag_se->excludes(kt_flag_pe)-> excludes(kt_option_r)->needs(kt_option_f);
  kt_option_f->check(CLI::ExistingPath);
  kt_option_r->check(CLI::ExistingPath);
  kt_option_o->required();
  kt_option_t->check(CLI::PositiveNumber);
  kt_option_p->check(CLI::IsMember(std::set<int> {0, 33, 64}));
  kt_option_m->check(CLI::Range(0.0, 1.0));
  kt_option_l->check(CLI::PositiveNumber);
  kt_option_q->check(CLI::PositiveNumber);
  kt_option_w->check(CLI::PositiveNumber);
  kt_option_k->check(CLI::IsMember(std::set<std::string>{"Illumina", "Nextera", "BGI", "Transposase"}));
  kt_option_c->check(CLI::Range(0, 9).description("compression level is limited to be between 0 and 9"));
  // kt_option_s->required();
  // kt_option_stat->required();
  // kt_option_pigz_th -> needs(kt_flag_pigz) -> check(CLI::PositiveNumber) -> check(CLI::Bound(2, 64));

  // ktrim params
  CLI11_PARSE(app, argc, argv);

  rabbit::trim::RabbitTrimParam rp;
  if(threads > 63) threads = 63;
  if(threads > get_nprocs_conf()) threads = get_nprocs_conf();
  rp.threads = threads; 
  rp.phred = phred;
  rp.forwardFiles = forwardFiles;
  rp.reverseFiles = reverseFiles;
  rp.output = output;
  rp.trimLog = trimLog;
  rp.stats = stats;
  rp.steps = steps;
  rp.validatePairing = validatePairing;
  rp.compressLevel = compressLevel;
  // rp.useIgzip = useIgzip;
  rp.usePigz = usePigz;
  rp.pigzThreadsNum = pigzThreadsNum;

  rp.seqA = seqA;
  rp.seqB = seqB;
  rp.seqKit = seqKit;
  rp.minLen = minLen;
  rp.minQual = minQual;
  rp.window = window;
  rp.mismatch = mismatch;
  rp.use_default_mismatch = (*kt_option_m) ? false : true;

  double start,finish;
  rabbit::Logger logger(true, true, quiet);
  start = rabbit::trim::util::getTime();
  if(app.got_subcommand(ktrim)){
    rp.prepare();
    logger.infoln("run the subcommand : ktrim");
  }else{
    logger.infoln("run the subcommand : trimmomatic");
  }
  logger.infoln("using thread nums : " + std::to_string(threads));

  logger.infoln("RabbitTrim: started with arguments: ");
  for(int i = 1; i < argc; i++)
  {
    std::cout << argv[i] << " "; 
  }
  std::cout << std::endl;

  if(isPE){
    rabbit::trim::process_pe(rp, logger);
  }
  else{
    rabbit::trim::process_se(rp, logger);
  }

  finish = rabbit::trim::util::getTime();
  std::cout<<"Time: "<<finish-start<<" s"<< std::endl;
  return 0;
}
