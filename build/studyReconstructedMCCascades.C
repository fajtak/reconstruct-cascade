#include "TChain.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TVector3.h"

#include <iostream>
#include <cstdlib>

std::vector<double> stringXPositions = {0,0,0,0,0,0,0,0};
std::vector<double> stringYPositions = {0,0,0,0,0,0,0,0};

void setStringPositions(int cluster)
{
	switch (cluster) {
		case 0: 
			stringXPositions = {57.13,35.66,-16.00,-53.33,-52.76,-9.57,39.28,-0.42};
			stringYPositions = {-3.48,-47.96,-53.44,-22.91,31.04,53.26,41.88,1.60};
			break;
		default:
			break;
	}
}

int studyReconstructedMCCascades(int year = 2016, int cluster = 0)
{
	TChain reconstructedCascades("nt_cascades");

	TString filesDir = "/Data/BaikalData/mc/2018may/recCasc_nTuple.root";
	// TString filesDir = "/Data/BaikalData/mc/nuatm_feb19/recCasc_nTuple.root";


	setStringPositions(cluster);
	reconstructedCascades.Add(filesDir);

	// reconstructedCascades.Add(filesDir.Data()); // add files,

	float runID, eventID, nPulses, nPulsesT, qRatio, closeHits, X, Y, Z, time, likelihood,energy,theta,phi,trueX,trueY,trueZ,trueEnergy,trueTheta,truePhi,efficiency,eventPurity,cascadePurity;
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
	reconstructedCascades.SetBranchAddress("TrueX", &trueX);
	reconstructedCascades.SetBranchAddress("TrueY", &trueY);
	reconstructedCascades.SetBranchAddress("TrueZ", &trueZ);
	reconstructedCascades.SetBranchAddress("TrueEnergy", &trueEnergy);
	reconstructedCascades.SetBranchAddress("TrueTheta", &trueTheta);
	reconstructedCascades.SetBranchAddress("TruePhi", &truePhi);
	reconstructedCascades.SetBranchAddress("Efficiency", &efficiency);
	reconstructedCascades.SetBranchAddress("EventPurity", &eventPurity);
	reconstructedCascades.SetBranchAddress("CascadePurity", &cascadePurity);

	TGraph* cascadePositions = new TGraph(reconstructedCascades.GetEntries());
	TGraph* zRhoPositions = new TGraph(reconstructedCascades.GetEntries());
	TGraph* stringPositions = new TGraph(8,&stringXPositions[0],&stringYPositions[0]);
	TGraph* NvsDistance = new TGraph(reconstructedCascades.GetEntries());

	TGraph* NvsEnergy = new TGraph(reconstructedCascades.GetEntries());

	TH1F* h_nHits = new TH1F("h_nHits","Number of hits created by cascade;NHits [#];NoE [#]",100,0,100);
	TH1F* h_X = new TH1F("h_X","X positions of reconstructed cascades",2000,-1000,1000);
	TH1F* h_Y = new TH1F("h_Y","Y positions of reconstructed cascades",2000,-1000,1000);
	TH1F* h_energy = new TH1F("h_energy","Cascades Energy;Energy [TeV];NoE [#]",2000,0,2000);
	TH1F* h_zenith = new TH1F("h_zenith","Zenith angle;Theta [degree];NoE [#]",180,0,180);
	TH1F* h_mismatchPosition = new TH1F("h_mismatchPosition","Mismatch Position;Distance [m];NoE [#]",100,0,100);
	TH1F* h_mismatchPosLong = new TH1F("h_mismatchPosLong","Mismatch Position Longitudinal Direction;Distance [m];NoE [#]",100,0,100);
	TH1F* h_mismatchPosPerp = new TH1F("h_mismatchPosPerp","Mismatch Position Perpendicular Direction;Distance [m];NoE [#]",100,0,100);
	TH1F* h_mismatchAngle = new TH1F("h_mismatchAngle","Mismatch angle;Mismatch angle [deg];NoE [#]",360,0,360);
	TH1F* h_mismatchEnergy = new TH1F("h_mismatchEnergy","Mismatch energy; E_{rec}/E_{true} [#];NoE [#]",500,0,5);
	TH1F* h_mismatchEnergyLog = new TH1F("h_mismatchEnergyLog","Mismatch energy; log10(E_{rec}/E_{true}) [#];NoE [#]",100,-1,1);

	TH2F* h_mismatchAngleEnergy = new TH2F("h_,h_mismatchAngleEnergy","Mismatch angle vs. cascade Energy;Cascade energy [TeV];Mismatch angle [deg]",1000,0,1000,360,0,360);

	TH1F* h_efficiency = new TH1F("h_efficiency","Hit selection efficiency;Efficiency [%];NoE [#]",120,0,1.2);
	TH1F* h_eventPurity = new TH1F("h_eventPurity","Event purity [%]; NoE [#]",120,0,1.2);
	TH1F* h_cascadePurity = new TH1F("h_cascadePurity","Cascade purity [%]; NoE [#]",120,0,1.2);

	for (int i = 0; i < reconstructedCascades.GetEntries(); ++i)
	{
		reconstructedCascades.GetEntry(i);
		double mismatchPosition = TMath::Sqrt(TMath::Power(trueX-X,2)+TMath::Power(trueY-Y,2)+TMath::Power(trueZ-Z,2));
		if (likelihood != -1 && cascadePurity == 1)
		{
			h_mismatchPosition->Fill(mismatchPosition);			
			TVector3 posTrue(trueX,trueY,trueZ);
			TVector3 posRec(X,Y,Z);
			TVector3 posDiff = posRec-posTrue;
			TVector3 cascDirTrue(0,0,1);
			cascDirTrue.SetTheta(trueTheta);
			cascDirTrue.SetPhi(truePhi);
			TVector3 cascDirRec(0,0,1);
			cascDirRec.SetTheta(theta);
			cascDirRec.SetPhi(phi);

			h_mismatchPosPerp->Fill(posDiff.Perp(cascDirTrue));
			h_mismatchPosLong->Fill(posDiff*cascDirTrue);
			h_mismatchAngle->Fill(cascDirRec.Angle(cascDirTrue)/TMath::Pi()*180);
			h_mismatchAngleEnergy->Fill(trueEnergy,cascDirRec.Angle(cascDirTrue)/TMath::Pi()*180);

			h_mismatchEnergyLog->Fill(TMath::Log10(energy/trueEnergy));
			h_mismatchEnergy->Fill(energy/trueEnergy);
			h_efficiency->Fill(efficiency);
			h_eventPurity->Fill(eventPurity);
			h_cascadePurity->Fill(cascadePurity);

			h_energy->Fill(trueEnergy);

			cout << theta << " " << trueTheta << " " << (theta-trueTheta)/TMath::Pi()*180 << " " << phi << " " << truePhi << " " << (phi-truePhi)/TMath::Pi()*180 << endl;

			// if (mismatchPosition > 10)
			// {
			// 	std::cout << i << "\t" << runID << "\t" << eventID << "\t\t" << nPulsesT << "\t" << qRatio << "\t" << closeHits << "\t" << X << "\t" << Y << "\t" << Z << std::endl;
			// 	std::cout << i << "\t" << likelihood << "\t" << energy << "\t\t" << theta << "\t" << phi << std::endl;
			// }
		}
		
		if (0 && Z < 560 && likelihood != -1)
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

	// TCanvas* c_cascadePositions = new TCanvas("c_cascadePositions","Results",800,600);

	// cascadePositions->SetMarkerStyle(21);
	// cascadePositions->SetTitle("Positions of reconstructed cascades;X [m];Y [m]");
	// cascadePositions->Draw("AP");
	// stringPositions->SetMarkerStyle(20);
	// stringPositions->SetMarkerColor(kRed);
	// stringPositions->Draw("PSame");

	// TCanvas* c_nHits = new TCanvas("c_nHits","Results",800,600);
	// h_nHits->Draw();

	// TCanvas* c_XandY = new TCanvas("c_XandY","Results",800,600);		
	// c_XandY->Divide(1,2);
	// c_XandY->cd(1);
	// h_X->Draw();
	// c_XandY->cd(2);
	// h_Y->Draw();

	// TCanvas* c_ZRho = new TCanvas("c_ZRho","Results",800,600);
	// zRhoPositions->SetTitle("Positions of reconstructed cascades;Rho [m];Z [m]");
	// zRhoPositions->SetMarkerStyle(21);
	// zRhoPositions->Draw("AP");

	// TCanvas* c_NvsD = new TCanvas("c_NvsD","Results",800,600);
	// NvsDistance->SetTitle("Number of hits vs distance from the center string;Rho [m];Nhits [#]");
	// NvsDistance->SetMarkerStyle(21);
	// NvsDistance->Draw("AP");

	// TCanvas* c_NvsE = new TCanvas("c_NvsE","Results",800,600);
	// NvsEnergy->SetTitle("Number of hits vs reconstructed energy;Energy [TeV];Nhits [#]");
	// NvsEnergy->SetMarkerStyle(21);
	// NvsEnergy->Draw("AP");	

	TCanvas* c_energy = new TCanvas("c_energy","Results",800,600);
	h_energy->Draw();

	// TCanvas* c_theta = new TCanvas("c_theta","Results",800,600);
	// h_zenith->Draw();

	TCanvas* c_mismatchPosition = new TCanvas("c_mismatchPosition","Results",800,600);
	h_mismatchPosition->Draw();

	Double_t x, q;
	q = 0.5; // 0.5 for "median"
	h_mismatchPosition->ComputeIntegral(); // just a precaution
	h_mismatchPosition->GetQuantiles(1, &x, &q);
	std::cout << "Median mismatch position = " << x << std::endl;

	TCanvas* c_mismatchPosPerp = new TCanvas("c_mismatchPosPerp","Results",800,600);
	h_mismatchPosPerp->Draw();

	h_mismatchPosPerp->ComputeIntegral(); // just a precaution
	h_mismatchPosPerp->GetQuantiles(1, &x, &q);
	std::cout << "Median mismatch position perpendicular = " << x << std::endl;

	TCanvas* c_mismatchPosLong = new TCanvas("c_mismatchPosLong","Results",800,600);
	h_mismatchPosLong->Draw();

	h_mismatchPosLong->ComputeIntegral(); // just a precaution
	h_mismatchPosLong->GetQuantiles(1, &x, &q);
	std::cout << "Median mismatch position longitudinal = " << x << std::endl;

	TCanvas* c_mismatchAngle = new TCanvas("c_mismatchAngle","Results",800,600);
	h_mismatchAngle->Draw();

	h_mismatchAngle->ComputeIntegral(); // just a precaution
	h_mismatchAngle->GetQuantiles(1, &x, &q);
	std::cout << "Median mismatch angle = " << x << std::endl;

	TCanvas* c_mismatchAngleEnergy = new TCanvas("c_mismatchAngleEnergy","Results",800,600);
	h_mismatchAngleEnergy->Draw("colz");

	TCanvas* c_mismatchEnergyLog = new TCanvas("c_mismatchEnergyLog","Results",800,600);
	h_mismatchEnergyLog->Draw();

	TCanvas* c_mismatchEnergy = new TCanvas("c_mismatchEnergy","Results",800,600);
	h_mismatchEnergy->Draw();

	TCanvas* c_effPur = new TCanvas("c_effPur","Results",800,600);
	c_effPur->Divide(2,2);
	c_effPur->cd(1);
	h_efficiency->Draw();
	c_effPur->cd(2);
	h_eventPurity->Draw();
	c_effPur->cd(3);
	h_cascadePurity->Draw();
	c_effPur->cd(4);

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
