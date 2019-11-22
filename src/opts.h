#pragma once
#ifndef PARSEOPTS
#define PARSEOPTS

#include <string>

extern int gNEventsProcessed;
extern int gVisEventID;
extern int gNCut;
extern int gQCut;
extern int gQCutHits;
extern int gQCutChi2;
extern int gTCutTimeWindowNs;
extern int gTCutHits;
extern int gTCutChi2;
extern int gZCut;
extern int gTDelayCut;
extern int gQRatioCut;
extern int gBranchCut;
extern double gScatteringCorrection;
extern std::string gProductionID;

void parseOpts(int argc, char** argv);
void readRC(const char* rcpath);
void checkParams();

#endif /*PARSEOPTS*/
