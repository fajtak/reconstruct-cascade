#include "opts.h"

#include <iostream>
#include <getopt.h>
#include <cstdlib>
#include <stdlib.h>

#include <TEnv.h>

#include "BARS.h"

int gNEventsProcessed = -1; 
int gNCut = -1;
double gQCut = -1;
int gQCutHits = -1;
double gQCutChi2 = -1;
double gTCutTimeWindowNs = -1;
int gTCutHits = -1;
double gTCutChi2 = -1;
double gZCut = -1;
double gTDelayCut = -1;
double gQRatioCut = -1;
int gBranchCut = -1;
int gVisEventID = -1;
double gScatteringCorrection = 10; //10 ns per 15 meters 
double gLikelihoodCut = -1;
bool gMC = false;
std::string gProductionID = "";

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
	{
		{
			"tag", 't',
			required_argument,
			"set production ID of the data that will be processed",
			[](char* argv) {gProductionID = argv;},
			[]() {;}
		},
		NOT_REQUIRED
	},
	{
		{
			"mc", 'm',
			no_argument,
			"use mc data as the input for the program",
			[](char* argv) {gMC = true;},
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
	gQCut = env.GetValue("QCut", 15.0);
	gQCutHits = env.GetValue("QCutHits", 6);
	gQCutChi2 = env.GetValue("QCutChi2", 50);
	gTCutTimeWindowNs = env.GetValue("TCutTimeWindowNs", 50);
	gTCutHits = env.GetValue("TCutHits", 20);
	gTCutChi2 = env.GetValue("TCutChi2", 50);
	gZCut = env.GetValue("ZCut",600);
	gTDelayCut = env.GetValue("TDelayCut",400);
	gQRatioCut = env.GetValue("QRatioCut",80);
	gBranchCut = env.GetValue("BranchCut",0);
	gLikelihoodCut = env.GetValue("LikelihoodCut",2.5);
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