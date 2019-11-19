#include "opts.h"

#include <iostream>
#include <getopt.h>
#include <cstdlib>
#include <stdlib.h>

#include <TEnv.h>

#include "BARS.h"

int gNEventsProcessed = -1; 
int gNCut = -1;
int gQCut = -1;
int gQCutHits = -1;
int gQCutChi2 = -1;
int gTCutTimeWindowNs = -1;
int gTCutHits = -1;
int gTCutChi2 = -1;
int gZCut = -1;
int gTDelayCut = -1;
int gQRatioCut = -1;
int gBranchCut = -1;
int gVisEventID = -1;
double gScatteringCorrection = 10; //10 ns per 15 meters 

using namespace BARS;

// At the end of list you should add 'App::opt_NULL' element
static const struct App::ProgramOption_t options_list[]{
  // --help and --config options will be added automaticly
	{App::opt_Verbose, NOT_REQUIRED},
	{App::opt_Input,   NOT_REQUIRED},
	{App::opt_Output,  NOT_REQUIRED},
	{App::opt_Cluster, NOT_REQUIRED},
	{App::opt_Season,  NOT_REQUIRED},
	{App::opt_Run,     NOT_REQUIRED},
	{
		{
			"number", 'n',
			required_argument,
			"set number of events you want to process",
			[](char* argv) {gNEventsProcessed = atoi(argv);},
			[]() {;}
		},
		NOT_REQUIRED
	},
	{
		{
			"view", 'w',
			required_argument,
			"set ID of the event that should be visualized",
			[](char* argv) {gVisEventID = atoi(argv);},
			[]() {;}
		},
		NOT_REQUIRED
	},
  // !!! must be at the end
	{App::opt_NULL, NOT_REQUIRED}
};

// Read resource file, fill parameters
//_____________________________________________________________________
void readRC(const char* rcpath)
{
	TEnv env;    
	if (-1 == env.ReadFile(App::RCPath, kEnvAll)) {
		std::cerr << "ERROR: failed to read configuration file at '" << App::RCPath << "'" << std::endl;
		exit(1);
	}
	gNCut = env.GetValue("NCut", 70);
	gQCut = env.GetValue("QCut", 2000);
	gQCutHits = env.GetValue("QCutHits", 6);
	gQCutChi2 = env.GetValue("QCutChi2", 50);
	gTCutTimeWindowNs = env.GetValue("TCutTimeWindowNs", 50);
	gTCutHits = env.GetValue("TCutHits", 20);
	gTCutChi2 = env.GetValue("TCutChi2", 50);
	gZCut = env.GetValue("ZCut",600);
	gTDelayCut = env.GetValue("TDelayCut",400);
	gQRatioCut = env.GetValue("QRatioCut",80);
	gBranchCut = env.GetValue("BranchCut",0);
}

// Parse options passed to the application.
//_____________________________________________________________________
void parseOpts(int argc, char **argv)
{
	App::ParseProgramOptions(argc, argv, options_list);
}

//_____________________________________________________________________
void checkParams()
{
	App::CheckProgramOptions(options_list);
}