#include "TChain.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TH1F.h"

#include <iostream>

std::vector<double> stringXPositions = {0,0,0,0,0,0,0,0};
std::vector<double> stringYPositions = {0,0,0,0,0,0,0,0};

void setStringPositions(int cluster)
{
	switch (cluster) {
		case 0: 
			/* code here */
			break;
		case 1:
			/* code here */
			break;
		case 2:
			stringXPositions = {-270.462,-225.154,-221.309,-257.28,-310.08,-337.212,-319.92,-277.709};
			stringYPositions = {-37.7285,-64.0684,-117.803,-154.07,-146.06,-101.399,-55.170,-98.5479};
			break;
		case 3:
			break;
		case 4:
			stringXPositions = {-163.912,-119.26,-113.895,-152.28,-202.59,-230.829,-213.25,-170.297};
			stringYPositions = {-628.257,-656.49,-707.52,-744.24,-738.58,-694.133,-645.06,-685.35};
			break;
		default:
			break;
	}
}

int studyReconstructedCascades(int year, int cluster)
{
	TChain reconstructedCascades("nt_cascades");
	TString filesDir = Form("/media/fajtak/Data/BaikalData/exp%d_barsv051/cluster%d/",year,cluster);

	setStringPositions(cluster);

	auto dir = gSystem->OpenDirectory(filesDir.Data());
	while (auto f = gSystem->GetDirEntry(dir)) {
	  	if (!strcmp(f,".") || !strcmp(f,"..")) continue;
	  	TString fullFilePath = filesDir + f + "/recCasc_nTuple.root";
	  	if (!gSystem->AccessPathName(fullFilePath))
	  		reconstructedCascades.Add(TString(filesDir) + f + "/recCasc_nTuple.root");
	}
	gSystem->FreeDirectory(dir);

	// reconstructedCascades.Add(filesDir.Data()); // add files,

	float runID, eventID, nPulses, nPulsesT, qRatio, closeHits, X, Y, Z, time;
	reconstructedCascades.SetBranchAddress("runID", &runID);
	reconstructedCascades.SetBranchAddress("EventID", &eventID);
	reconstructedCascades.SetBranchAddress("NPulses", &nPulses);
	reconstructedCascades.SetBranchAddress("NPulsesT", &nPulsesT);
	reconstructedCascades.SetBranchAddress("QRatio", &qRatio);
	reconstructedCascades.SetBranchAddress("CloseHits", &closeHits);
	reconstructedCascades.SetBranchAddress("X", &X);
	reconstructedCascades.SetBranchAddress("Y", &Y);
	reconstructedCascades.SetBranchAddress("Z", &Z);
	reconstructedCascades.SetBranchAddress("Time", &time);

	TGraph* cascadePositions = new TGraph(reconstructedCascades.GetEntries());
	TGraph* stringPositions = new TGraph(8,&stringXPositions[0],&stringYPositions[0]);

	TH1F* h_nHits = new TH1F("h_nHits","Number of hits created by cascade;NHits [#];NoE [#]",100,0,100);

	for (int i = 0; i < reconstructedCascades.GetEntries(); ++i)
	{
		reconstructedCascades.GetEntry(i);
		std::cout << i << "\t" << runID << "\t" << eventID << "\t\t" << nPulsesT << "\t" << qRatio << "\t" << closeHits << std::endl;
		cascadePositions->SetPoint(i,X,Y);
		h_nHits->Fill(nPulsesT);
	}

	TCanvas* c_cascadePositions = new TCanvas("c_cascadePositions","Results",800,600);

	cascadePositions->SetMarkerStyle(21);
	cascadePositions->SetTitle("Positions of reconstructed cascades;X [m];Y [m]");
	cascadePositions->Draw("AP");
	stringPositions->SetMarkerStyle(20);
	stringPositions->SetMarkerColor(kRed);
	stringPositions->Draw("PSame");

	TCanvas* c_nHits = new TCanvas("c_nHits","Results",800,600);
	h_nHits->Draw();

	return 1;
}


// // wildcards work
// // define variables and assign them to the corresponding branches
// float pot, cur, temp, pres;
// my_tuple->SetBranchAddress("Potential", &pot);
// my_tuple->SetBranchAddress("Current", &cur);
// my_tuple->SetBranchAddress("Temperature", &temp);
// my_tuple->SetBranchAddress("Pressure", &pres);
// cout << "Potential\tCurrent\tTemperature\tPressure\n";
// for (size_t irow=0; irow<in_chain.GetEntries(); ++irow){
// in_chain.GetEntry(irow); // loads all variables that have
// // been connected to branches
// cout << pot << "\t" << cur << "\t" << temp <<
// "\t" << pres << endl;
//}