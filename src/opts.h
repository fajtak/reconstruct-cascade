#pragma once
#ifndef PARSEOPTS
#define PARSEOPTS

#include <string>

extern int gNEventsProcessed;
extern int gVisEventID;
extern int gNCut;
extern double gQCut;
extern int gQCutHits;
extern double gQCutChi2;
extern double gTCutTimeWindowNs;
extern int gTCutHits;
extern double gTCutChi2;
extern double gZCut;
extern double gTDelayCut;
extern double gQRatioCut;
extern int gBranchCut;
extern double gScatteringCorrection;
extern double gLikelihoodCut;
extern bool gMC;
extern std::string gProductionID;

void parseOpts(int argc, char** argv);
void readRC(const char* rcpath);
void checkParams();

#endif /*PARSEOPTS*/
