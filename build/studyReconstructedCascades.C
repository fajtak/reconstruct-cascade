#include "TChain.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TMath.h"

#include <iostream>
#include <cstdlib>

std::vector<double> stringXPositions = {0,0,0,0,0,0,0,0};
std::vector<double> stringYPositions = {0,0,0,0,0,0,0,0};

void setStringPositions(int cluster)
{
	switch (cluster) {
		case 0: 
			stringXPositions = {-13.76,32.14,45.06,5.13,-45.034,-76.205,-59.85,-14.465};
			stringYPositions = {-211.35,-235.88,-285.45,-325.828,-319.815,-281.633,-231.37,-270.172};
			break;
		case 1:
			stringXPositions = {-195.186,-164.79,-180.08,-227.51,-276.243,-279.59,-248.17,-222.698};
			stringYPositions = {-340.621,-384.09,-435.125,-450.127,-424.314,-372.59,-337.03,-391.09};
			break;
		case 2:
			stringXPositions = {-270.25,-228.578,-220.89,-261.89,-309.856,-337.48,-319.744,-282.265};
			stringYPositions = {-37.36,-65.264,-117.78,-153.57,-146.267,-101.438,-55.245,-96.822};
			break;
		case 3:
			stringXPositions = {65.849,108.73,113.87,74.189,25.1,-2.48,16.08,58.372};
			stringYPositions = {-435.471,-462.39,-514.683,-549.898,-544.25,-500.533,-453,-491.965};
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

	TString filesDir;
	if(const char* env_p = std::getenv("CEPH_MNT"))
	{
    	// filesDir = Form("%s/exp%d_barsv051/cluster%d/",env_p,year,cluster);
    	// filesDir = Form("/Data/BaikalData/dataLog/exp20%d/cluster%d/",year,cluster);
    	// filesDir = Form("/Data/BaikalData/dataLog_HQ/exp20%d/cluster%d/",year,cluster);
    	// filesDir = Form("/Data/BaikalData/dataLog_grid/exp20%d/cluster%d/",year,cluster);
    	filesDir = Form("/Data/BaikalData/dataLog_gridII/exp20%d/cluster%d/",year,cluster);
		cout << env_p << endl;
	}else{	
		cout << "SET ENVIROMENT VARIABLE CEPH_MNT" << endl;
		return -1;
	}

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

	float runID, eventID, nPulses, nPulsesT, qRatio, closeHits, X, Y, Z, time, likelihood,energy,theta,phi;
	reconstructedCascades.SetBranchAddress("runID", &runID);
	reconstructedCascades.SetBranchAddress("EventID", &eventID);
	reconstructedCascades.SetBranchAddress("NPulses", &nPulses);
	reconstructedCascades.SetBranchAddress("NPulsesT", &nPulsesT);
	reconstructedCascades.SetBranchAddress("QRatio", &qRatio);
	reconstructedCascades.SetBranchAddress("CloseHits", &closeHits);
	reconstructedCascades.SetBranchAddress("Likelihood", &likelihood);
	reconstructedCascades.SetBranchAddress("X", &X);
	reconstructedCascades.SetBranchAddress("Y", &Y);
	reconstructedCascades.SetBranchAddress("Z", &Z);
	reconstructedCascades.SetBranchAddress("Time", &time);
	reconstructedCascades.SetBranchAddress("Energy", &energy);
	reconstructedCascades.SetBranchAddress("Theta", &theta);
	reconstructedCascades.SetBranchAddress("Phi", &phi);

	TGraph* cascadePositions = new TGraph(reconstructedCascades.GetEntries());
	TGraph* zRhoPositions = new TGraph(reconstructedCascades.GetEntries());
	TGraph* stringPositions = new TGraph(8,&stringXPositions[0],&stringYPositions[0]);
	TGraph* NvsDistance = new TGraph(reconstructedCascades.GetEntries());

	TGraph* NvsEnergy = new TGraph(reconstructedCascades.GetEntries());

	TH1F* h_nHits = new TH1F("h_nHits","Number of hits created by cascade;NHits [#];NoE [#]",100,0,100);
	TH1F* h_X = new TH1F("h_X","X positions of reconstructed cascades",2000,-1000,1000);
	TH1F* h_Y = new TH1F("h_Y","Y positions of reconstructed cascades",2000,-1000,1000);
	TH1F* h_energy = new TH1F("h_energy","Cascades Energy;Energy [TeV];NoE [#]",200,0,2000);
	TH1F* h_zenith = new TH1F("h_zenith","Zenith angle;Theta [degree];NoE [#]",180,0,180);

	for (int i = 0; i < reconstructedCascades.GetEntries(); ++i)
	{
		reconstructedCascades.GetEntry(i);
		if (nPulsesT < 25 || likelihood < 0)
			continue;
		// if (TMath::Abs(X-stringXPositions[7]) > 60 || TMath::Abs(Y-stringYPositions[7]) > 60)
			// continue; 
		if (likelihood > 2 || energy < 100)
			continue;
		std::cout << i << "\t" << runID << "\t" << eventID << "\t" << nPulsesT << "\\" << nPulses << "\t" << qRatio << "\t" << closeHits << "\t" << X-stringXPositions[7] << "\t" << Y-stringYPositions[7] << "\t" << Z << std::endl;
		std::cout << i << "\t" << likelihood << "\t" << energy << "\t\t" << theta << "\t" << phi << std::endl;
		if (Z < 600)
		{
			cascadePositions->SetPoint(i,X,Y);
			zRhoPositions->SetPoint(i,TMath::Sqrt(TMath::Power(X-stringXPositions[7],2)+TMath::Power(Y-stringYPositions[7],2)),Z);
			NvsDistance->SetPoint(i,TMath::Sqrt(TMath::Power(X-stringXPositions[7],2)+TMath::Power(Y-stringYPositions[7],2)),nPulsesT);
			NvsEnergy->SetPoint(i,energy,nPulsesT);
			h_X->Fill(X);
			h_Y->Fill(Y);			
			h_nHits->Fill(nPulsesT);		
			h_energy->Fill(energy);	
			h_zenith->Fill((theta/TMath::Pi()*180));	
			// if ()
		}
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

	TCanvas* c_XandY = new TCanvas("c_XandY","Results",800,600);		
	c_XandY->Divide(1,2);
	c_XandY->cd(1);
	h_X->Draw();
	c_XandY->cd(2);
	h_Y->Draw();

	TCanvas* c_ZRho = new TCanvas("c_ZRho","Results",800,600);
	zRhoPositions->SetTitle("Positions of reconstructed cascades;Rho [m];Z [m]");
	zRhoPositions->SetMarkerStyle(21);
	zRhoPositions->Draw("AP");

	TCanvas* c_NvsD = new TCanvas("c_NvsD","Results",800,600);
	NvsDistance->SetTitle("Number of hits vs distance from the center string;Rho [m];Nhits [#]");
	NvsDistance->SetMarkerStyle(21);
	NvsDistance->Draw("AP");

	TCanvas* c_NvsE = new TCanvas("c_NvsE","Results",800,600);
	NvsEnergy->SetTitle("Number of hits vs reconstructed energy;Energy [TeV];Nhits [#]");
	NvsEnergy->SetMarkerStyle(21);
	NvsEnergy->Draw("AP");	

	TCanvas* c_energy = new TCanvas("c_energy","Results",800,600);
	h_energy->Draw();

	TCanvas* c_theta = new TCanvas("c_theta","Results",800,600);
	h_zenith->Draw();

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
