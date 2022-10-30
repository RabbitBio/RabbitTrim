#include <sys/sysinfo.h>
#include <string>
#include <vector>
#include "CLI11.hpp"
#include "param.h"
#include "core/HandlerPE.h"
#include "core/HandlerSE.h"

int main( int argc, char **argv) {
    
	// process the command line parameters
    CLI::App app{"RabbitTrim : Optimised versions of Trimmomatic and Ktrim tools"};
    app.require_subcommand(1);
    CLI::App* ktrim = app.add_subcommand("ktrim", "Based Ktrim");
    CLI::App* trimmomatic = app.add_subcommand("trimmomatic", "Based Trimmomatic");
    // std::string mode;
    bool isPE = false; 
    bool isSE = false; 
    int threads = get_nprocs_conf();
    int phred = 0;
    std::vector<std::string> forwardFiles;
    std::vector<std::string> reverseFiles;
    std::string output;
    std::string trimLog;
    std::string stats;
    std::vector<std::string> steps;
    bool quiet = false;
    bool validatePairing = false;
    int compressLevel = 4;

    // auto trimm_option_m = trimmomatic->add_option("-m,--mode", mode, "specify mode of data handled: single-end(SE) or pair-end data(PE)");
    auto trimm_flag_pe = trimmomatic->add_flag("--PE,!--no-PE", isPE, "specify whether to be pair end data ");
    auto trimm_flag_se = trimmomatic->add_flag("--SE,!--no-SE", isSE, "specify whether to be single end data");
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
    
    // trimm_option_m->required()->check(CLI::IsMember({"PE", "SE"}));
    trimm_flag_pe->excludes(trimm_flag_se)->needs(trimm_option_f)->needs(trimm_option_r);
    trimm_flag_se->excludes(trimm_flag_pe);
    trimm_option_t->check(CLI::PositiveNumber);
    trimm_option_p->check(CLI::IsMember(std::set<int> {0, 33, 64}));
    trimm_option_f->check(CLI::ExistingPath);
    trimm_option_r->check(CLI::ExistingPath);
    trimm_option_o->required()->check(CLI::ExistingPath);
    trimm_option_l->required()->check(CLI::ExistingPath);
    trimm_option_stat->required()->check(CLI::ExistingPath);
    trimm_option_s->required();
    trimm_option_c->check(CLI::Range(1, 9).description("compression level is limited to be between 1 and 9"));
    CLI11_PARSE(app, argc, argv);
    
    rabbit::trim::RabbitTrimParam rp;
    rp.threads = threads;
    rp.phred = phred;
    rp.forwardFiles = forwardFiles;
    rp.reverseFiles = reverseFiles;
    rp.output = output;
    rp.trimLog = trimLog;
    rp.stats = stats;
    rp.steps = steps;
    rp.validatePairing = validatePairing;

	double start,finish;
	rabbit::Logger logger(true, true, quiet);
	start = getTime();
    if(app.got_subcommand(ktrim)){
        logger.infoln("run the subcommand : ktrim");
    }else{
        if(app.got_subcommand(trimmomatic)){
            logger.infoln("run the subcommand : trimmomatic");
        }
    }
    
    if(isPE){
        rabbit::trim::process_pe(rp, logger);
    }
    else{
        rabbit::trim::process_se(rp, logger);
    }
    
	finish = getTime();
	cout<<"Time: "<<finish-start<<" s"<<endl;
	return 0;
}
