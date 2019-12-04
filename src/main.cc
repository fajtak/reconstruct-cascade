//program headers
#include "opts.h"

//system headers
#include <iostream>
#include <vector>

//ROOT dependencies
#include "TTree.h"
#include "TFile.h"
#include "TVirtualFitter.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TNtuple.h"

//BARS dependencies
#include "BARS.h"
#include "BExtractedImpulseTel.h"
#include "BExtractedHeader.h"
#include "BGeomTel.h"
#include "BEvoGeomApply.h"

TVirtualFitter* fMinuit;

const int gNOMs = 288;
const int gNStrings = 8;
const double ReciprocalSpeedOfLightinWater = 4.57;

vector<TVector3> gOMpositions(288);
vector<double> gOMtimeCal(288);
vector<double> gOMQCal(288);

double gLogTable4D[200][351][21][21]{0};

struct PulsesVariables
{
  int OMID;
  double time;
  double charge;
  double expectedCharge;
};

vector<PulsesVariables> g_pulses; //global vector of PulsesVariables

void PrintGPulses()
{
	for (int i = 0; i < g_pulses.size(); ++i)
	{
		cout << i << " " << g_pulses[i].OMID << " " << g_pulses[i].charge << " " << g_pulses[i].time << endl;
	}
}

using namespace BARS;

TH1F* h_nHits = new TH1F("h_nHits","Number of hits;Hits [#]; NoE [#]",200,0,200);
TH1F* h_chargePerHit = new TH1F("h_chargePerHit","Charge per hit;Q [FADC];NoE [#]",10000,0,100000);
TH1F* h_nHitsAfterQ = new TH1F("h_nHitsAfterQ","Number of hits after Q Filter; Nhits [#];NoE [#]",20,0,20);
TH1F* h_chi2AfterQ = new TH1F("h_chi2AfterQ","Chi2 after Q Filter; chi^2 [#];NoE [#]",100,0,100);
TH1F* h_nHitsAfterT = new TH1F("h_nHitsAfterT","Number of hits after T Filter; NHits [#];NoE [#]",100,0,100);
TH1F* h_chi2AfterT= new TH1F("h_chi2AfterT","Chi2 after T Filter; chi^2 [#];NoE [#]",100,0,100);
TGraph* g_cascadePosition = new TGraph();
TH1F* h_timeRes = new TH1F("h_timeRes","Time residuals;DeltaT [ns];NoE [#]",400,-200,200);
TH1F* h_nHitOM = new TH1F("h_nHitOM","Number of hits per OM per run",288,0,288);
TH2F* h_XYposition = new TH2F("h_XYposition","XvsY cascade positions;X [m];Y [m]",500,-500,500,500,-500,500);

// Function that checks if all arguments have been set
bool CheckParams()
{
	if (BARS::App::Season == -1)
    {
    	std::cout << "Set season with -s" << std::endl;
    	return false;
    }

    if (BARS::App::Cluster == -1)
    {
    	std::cout << "Set cluster with -c" << std::endl;
    	return false;
    }

    if (BARS::App::Run == -1)
    {
    	std::cout << "Set run with -r" << std::endl;
    	return false;
    }
    return true;
}

// Input: calculated parameters R,Z,phi,cosTheta Output: given lower indexes in 4D array 
int GetIndexes(double* inputs, int* outputs)
{
	if (inputs[0] < 0)
		return -1;
	if (inputs[2] < -1 || inputs[2] > 1)
		return -1;
	if (inputs[3] < 0 || inputs[3] > TMath::Pi())
		return -1;

	if (inputs[0] < 199 )
		outputs[0] = (int)floor(inputs[0]);
	else
	{
		outputs[0] = 198;
		inputs[0] = 199;
	}
	if (inputs[1] < 150 && inputs[1] >= -200)
		outputs[1] = (int)floor(inputs[1]+200);
	else
	{
		if (inputs[1] < -200)
		{
			inputs[1] = -200;
			outputs[1] = 0;
		}else{
			outputs[1] = 349;
			inputs[1] = 150;
		}
	}
	outputs[2] = (int)floor((inputs[2]+1)/0.1);
	outputs[3] = (int)floor(inputs[3]/(TMath::Pi()/20));

	return 0;
}

double interpolate(double x, double x1, double x2, double f1, double f2)
{
	return (x-x1)/(x2-x1)*f2 + (x2-x)/(x2-x1)*f1;
}

double interpolate4D(double* inputs, int* indexes)
{
	double R0Z0Phi0 = interpolate(inputs[3],(indexes[3])*(TMath::Pi()/20),(indexes[3]+1)*(TMath::Pi()/20),gLogTable4D[indexes[0]][indexes[1]][indexes[2]][indexes[3]],gLogTable4D[indexes[0]][indexes[1]][indexes[2]][indexes[3]+1]);
	double R0Z0Phi1 = interpolate(inputs[3],(indexes[3])*(TMath::Pi()/20),(indexes[3]+1)*(TMath::Pi()/20),gLogTable4D[indexes[0]][indexes[1]][indexes[2]+1][indexes[3]],gLogTable4D[indexes[0]][indexes[1]][indexes[2]+1][indexes[3]+1]); 

	double R0Z1Phi0 = interpolate(inputs[3],(indexes[3])*(TMath::Pi()/20),(indexes[3]+1)*(TMath::Pi()/20),gLogTable4D[indexes[0]][indexes[1]+1][indexes[2]][indexes[3]],gLogTable4D[indexes[0]][indexes[1]][indexes[2]][indexes[3]+1]);
	double R0Z1Phi1 = interpolate(inputs[3],(indexes[3])*(TMath::Pi()/20),(indexes[3]+1)*(TMath::Pi()/20),gLogTable4D[indexes[0]][indexes[1]+1][indexes[2]+1][indexes[3]],gLogTable4D[indexes[0]][indexes[1]][indexes[2]][indexes[3]+1]);

	double R1Z0Phi0 = interpolate(inputs[3],(indexes[3])*(TMath::Pi()/20),(indexes[3]+1)*(TMath::Pi()/20),gLogTable4D[indexes[0]+1][indexes[1]][indexes[2]][indexes[3]],gLogTable4D[indexes[0]][indexes[1]][indexes[2]][indexes[3]+1]);
	double R1Z0Phi1 = interpolate(inputs[3],(indexes[3])*(TMath::Pi()/20),(indexes[3]+1)*(TMath::Pi()/20),gLogTable4D[indexes[0]+1][indexes[1]][indexes[2]+1][indexes[3]],gLogTable4D[indexes[0]][indexes[1]][indexes[2]][indexes[3]+1]);

	double R1Z1Phi0 = interpolate(inputs[3],(indexes[3])*(TMath::Pi()/20),(indexes[3]+1)*(TMath::Pi()/20),gLogTable4D[indexes[0]+1][indexes[1]+1][indexes[2]][indexes[3]],gLogTable4D[indexes[0]][indexes[1]][indexes[2]][indexes[3]+1]);
	double R1Z1Phi1 = interpolate(inputs[3],(indexes[3])*(TMath::Pi()/20),(indexes[3]+1)*(TMath::Pi()/20),gLogTable4D[indexes[0]+1][indexes[1]+1][indexes[2]+1][indexes[3]],gLogTable4D[indexes[0]][indexes[1]][indexes[2]][indexes[3]+1]);

	double R0Z0 = interpolate(inputs[2],(indexes[2])*0.1-1,(indexes[2]+1)*0.1-1,R0Z0Phi0,R0Z0Phi1);
	double R0Z1 = interpolate(inputs[2],(indexes[2])*0.1-1,(indexes[2]+1)*0.1-1,R0Z1Phi0,R0Z1Phi1);

	double R1Z0 = interpolate(inputs[2],(indexes[2])*0.1-1,(indexes[2]+1)*0.1-1,R1Z0Phi0,R1Z0Phi1);
	double R1Z1 = interpolate(inputs[2],(indexes[2])*0.1-1,(indexes[2]+1)*0.1-1,R1Z1Phi0,R1Z1Phi1);

	double R0 = interpolate(inputs[1],(indexes[1]-200),(indexes[1]+1-200),R0Z0,R0Z1);
	double R1 = interpolate(inputs[1],(indexes[1]-200),(indexes[1]+1-200),R1Z0,R1Z1);

	return interpolate(inputs[0],(indexes[0]),(indexes[0]+1),R0,R1);
}

// interpolation of the values from 4D log likelihood table. Input: 4D array, R,Z,phi,cosTheta
double GetInterpolatedValue(double* in)
{
	int out[4]{0};

	int returnValue = GetIndexes(in,out);
	if (returnValue == -1)
	{
		cerr << "Table parameters ERROR" << endl;
		exit(1);
	}
	// cout << "IN values " << in[0] << " " << in[1] << " " << in[2] << " " << in[3] << endl;
	// cout << "Out values " << out[0] << " " << out[1] << " " << out[2] << " " << out[3] << endl;
	// cout << interpolate4D(in,out) << endl;
	return interpolate4D(in,out);
}

double ExpectedTime(double matrixX, double matrixY, double matrixZ, double matrixTime,int OMID)
{
  	double distanceFromLEDmatrix = TMath::Sqrt(TMath::Power(gOMpositions[OMID].X()-matrixX,2)+TMath::Power(gOMpositions[OMID].Y()-matrixY,2)+TMath::Power(gOMpositions[OMID].Z()-matrixZ,2));
  	double scattering_correction = (gOMpositions[OMID].Z() < matrixZ) ? (gScatteringCorrection/15.0)*(matrixZ - gOMpositions[OMID].Z()) : 0;
  	// double scattering_correction = (gOMpositions[OMID].Z < LEDmatrixZ) ? (scatteringCorrection/15.0)*(LEDmatrixZ - myOMs[n].Z) : 0;

  	double expectedTime = matrixTime + distanceFromLEDmatrix*ReciprocalSpeedOfLightinWater + scattering_correction; //v nanosekundach

  	return  expectedTime;
}

// minimization function
void chi2(Int_t &npar, Double_t* gin, Double_t &f, Double_t* par, Int_t iflag) //keep all these nonsense parameters
{
	double theChi2 = 0;
	double constant = 1.0/(g_pulses.size() - 4);
  	double error = 1.0/3; //error is 3 ns - photomultiplier
  	double theChi;

	for (int i = 0; i < g_pulses.size(); ++i) 
	{       
		// chi calculation (measured - theoretical)
		theChi = constant*error*(g_pulses[i].time - ExpectedTime(par[0], par[1], par[2], par[3], g_pulses[i].OMID))*TMath::Log10(g_pulses[i].charge/150+1);
	  	// theChi = constant*error*(g_pulses[i].Times - ExpectedTime(myOMs,g_pulses[i].OMID, par[0], par[1], par[2], par[3]))*log(g_pulses[i].Charge);
		theChi2 += theChi*theChi;  
	}
	f = theChi2; // function returns calculated chi2, it is global function which is used in SetFCN()
}

void MEstimator(Int_t &npar, Double_t* gin, Double_t &f, Double_t* par, Int_t iflag)
{	
	double sum = 0;
	for (int i = 0; i < g_pulses.size(); ++i) 
	{       
		double Tres = (g_pulses[i].time - ExpectedTime(par[0], par[1], par[2], par[3], g_pulses[i].OMID));
		double element = g_pulses[i].charge/150*TMath::Sqrt(1+TMath::Power(Tres,2)/2);
		sum += element;
	}
	f = sum; // function returns calculated chi2, it is global function which is used in SetFCN()

}

// cascadeParameters[7]: X,Y,Z,Time,Energy,Theta,Phi
// tableParameters[4]: R,Z,Phi',CosTheta'
void GetParameters(const Double_t* cascadeParameters,const int OMID, double* tableParameters)
{
	TVector3 OMpos(gOMpositions[OMID].X()-cascadeParameters[0],gOMpositions[OMID].Y()-cascadeParameters[1],gOMpositions[OMID].Z()-cascadeParameters[2]);

	TVector3 showerRef(0,0,1);
	showerRef.SetTheta(cascadeParameters[5] - TMath::Pi());
	showerRef.SetPhi(cascadeParameters[6] + TMath::Pi());
	// showerRef.Print();
	// OMpos.Print();
	
	tableParameters[0] = OMpos.Perp(showerRef);
	tableParameters[1] = OMpos*showerRef;

	TVector3 OMorien(0,0,-1);
	tableParameters[2] = showerRef*OMorien;

	TVector3 x2 = showerRef.Orthogonal();
	x2.SetMag(1.0);
	// x2.Print();
	TVector3 x3 = showerRef.Cross(x2);
	x3.SetMag(1.0);
	// x3.Print();
	TVector3 rhoProjection(OMpos*x2,OMpos*x3,0);
	TVector3 omegaProjection(OMorien*x2,OMorien*x3,0);

	tableParameters[3] = rhoProjection.Angle(omegaProjection);
}

void logLikelihood(Int_t &npar, Double_t* gin, Double_t &f, Double_t* par, Int_t iflag)
{
	double logLike = 0;
	double tableParameters[4]{0};
	// cout << "Calculating logLike" << endl;
	// cout << par[0] << " " << par[1] << " " << par[2] << " " << par[3] << " " << par[4] << " " << par[5] << " " << par[6] << endl;
	for (int i = 0; i < g_pulses.size(); ++i)
	{
		GetParameters(par,g_pulses[i].OMID,tableParameters);
		// gOMpositions[g_pulses[i].OMID].Print();
		// cout << i << " " << tableParameters[0] << " " << tableParameters[1] << " " << tableParameters[2] << " " << tableParameters[3] << endl;
		double expectedNPE = GetInterpolatedValue(tableParameters);
		// cout << expectedNPE << " " << expectedNPE*100000000*par[4] << " " << g_pulses[i].charge/120 << " " << TMath::PoissonI(g_pulses[i].charge/120,expectedNPE*100000000*par[4]) << endl;
		
		if (TMath::Poisson(g_pulses[i].charge/gOMQCal[g_pulses[i].OMID],expectedNPE*100000000*par[4]) > 10e-350)
		{
			logLike -= TMath::Log10(TMath::Poisson(g_pulses[i].charge/gOMQCal[g_pulses[i].OMID],expectedNPE*100000000*par[4]));
			// cout << TMath::Log10(TMath::PoissonI(g_pulses[i].charge/120,expectedNPE*100000000*par[4])) << endl;
		}
		else
		{
			logLike -= -350;
			// cout << -350 << endl;
		}
	}
	// cout << "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF: " << logLike << endl;
	f = logLike;
}

void FillExpectedCharges(double X, double Y, double Z, double energy, double theta, double phi)
{
	double tableParameters[4]{0};
	double par[7]{X,Y,Z,0,energy,theta,phi};
	for (int i = 0; i < g_pulses.size(); ++i)
	{
		GetParameters(par,g_pulses[i].OMID,tableParameters);
		g_pulses[i].expectedCharge = GetInterpolatedValue(tableParameters)*100000000*energy;
	}
}

double FitMatrixPosition(TVector3 &R2, double &T2)
{ 
	fMinuit->ReleaseParameter(0);
	fMinuit->ReleaseParameter(1);
	fMinuit->ReleaseParameter(2);
	fMinuit->ReleaseParameter(3);

	//cout<<"fitEVent "<<endl;
	fMinuit->SetParameter(0,"LED X",R2.X(),0.01,-1000,1000); //estimation of start point
	fMinuit->SetParameter(1,"LED Y",R2.Y(),0.01,-1000,1000);
	fMinuit->SetParameter(2,"LED Z",R2.Z(),0.01,0,1000);
	fMinuit->SetParameter(3,"Time",T2,1,0,20000);

	fMinuit->FixParameter(4);
	fMinuit->FixParameter(5);
	fMinuit->FixParameter(6);

	double arglist[10] = {500,0.01}; //500 means 500 iterations and then MIGRAD stops
	fMinuit->ExecuteCommand("MIGRAD",arglist,2); //this line performs minimalisation

	double chi2, edm, errdef;
	int nvpar, nparx;

	fMinuit->GetStats(chi2,edm,errdef,nvpar,nparx);

	R2.SetX(fMinuit->GetParameter(0));
	R2.SetY(fMinuit->GetParameter(1));
	R2.SetZ(fMinuit->GetParameter(2));
	T2 = fMinuit->GetParameter(3);

	 // cout<<"Chi2  "<<chi2/(g_pulses.size()-4)<<endl;
	/*cout<<"Pozicia kaskady po Fite "<<R2.X()<<endl;
	cout<<"Pozicia kaskady po Fite "<<R2.Y()<<endl;
	cout<<"Pozicia kaskady po Fite "<<R2.Z()<<endl;
	cout<<"Cas                     "<<T2<<endl;
	cout<<" "<<endl;*/

	return chi2/(g_pulses.size()-4);  
}

// Fitter setting
void SetFitter(void)
{
	fMinuit = TVirtualFitter::Fitter(0,7); // the second number is number of parameters
	double arg = -1;
	fMinuit->ExecuteCommand("SET PRINTOUT",&arg,1); // these two lines means that it wont be able to print results on the screen
	fMinuit->ExecuteCommand("SET NOW", &arg ,1);
	// fMinuit->SetFCN(chi2);
	fMinuit->SetFCN(MEstimator);
}

double FitMatrixDirection(TVector3 &R2, double &T2, double &energy, double &theta, double &phi)
{ 
	fMinuit->ReleaseParameter(4);
	fMinuit->ReleaseParameter(5);
	fMinuit->ReleaseParameter(6);

	// cout << "StartFitting" << endl;
	fMinuit->SetFCN(logLikelihood);
	fMinuit->SetParameter(0,"LED X",R2.X(),0.01,-1000,1000); //estimation of start point
	fMinuit->SetParameter(1,"LED Y",R2.Y(),0.01,-1000,1000);
	fMinuit->SetParameter(2,"LED Z",R2.Z(),0.01,0,1000);
	fMinuit->SetParameter(3,"Time",T2,1,0,20000);
	fMinuit->SetParameter(4,"Energy",10,10,0,10000);
	fMinuit->SetParameter(5,"Theta",theta,0.1,0,TMath::Pi());
	fMinuit->SetParameter(6,"Phi",phi,0.1,0,2*TMath::Pi());

	fMinuit->FixParameter(0);
	fMinuit->FixParameter(1);
	fMinuit->FixParameter(2);
	fMinuit->FixParameter(3);

	double arglist[10] = {500,0.01}; //500 means 500 iterations and then MIGRAD stops
	fMinuit->ExecuteCommand("SIMPLEX",arglist,2); //this line performs minimalisation

	double chi2, edm, errdef;
	int nvpar, nparx;

	fMinuit->GetStats(chi2,edm,errdef,nvpar,nparx);
	energy = fMinuit->GetParameter(4);
	theta = fMinuit->GetParameter(5);
	phi = fMinuit->GetParameter(6);

	return chi2/(g_pulses.size()-4);  
}

void PrintConfig(void)
{
	std::cout << std::string(80,'*') << std::endl;
	std::cout << "Reconstruction configuration: " << std::endl;
	std::cout << std::string(80,'*') << std::endl;
	std::cout << "NFilter: " << gNCut << endl;
	std::cout << "QFilter: " << gQCut << endl;
	std::cout << "QFilterHits: " << gQCutHits << endl;
	std::cout << "QFilterChi2: " << gQCutChi2 << endl;
	std::cout << "TFilterTimeWindow: " << gTCutTimeWindowNs << endl;
	std::cout << "TFilterHits: " << gTCutHits << endl;
	std::cout << "TFilterChi2: " << gTCutChi2 << endl;
	std::cout << "ZFilter: " << gZCut << endl;
	std::cout << "TDelayFilter: " << gTDelayCut << endl;
	std::cout << "QRatioFilter: " << gQRatioCut << endl;
	std::cout << "BranchFilter: " << gBranchCut << endl;
	std::cout << "LikelihoodFilter: " << gLikelihoodCut << endl;
	std::cout << std::string(80,'*') << std::endl;
}

bool NFilterPassed(BExtractedImpulseTel* impulseTel, int &nPulses)
{
	nPulses = impulseTel->GetNimpulse();
	h_nHits->Fill(nPulses);
	if (nPulses >= gNCut)
	{
		return true;
	}
	return false;
}

bool OMIDAlreadyInGPulses(BExtractedImpulse* impulse)
{
	bool OMIDAlreadyIn = false;
	for (int k = 0; k < g_pulses.size(); ++k)
	{
		if (g_pulses[k].OMID == impulse->GetNch())
		{
			if (g_pulses[k].charge < impulse->GetQ())
			{
				g_pulses[k].charge = impulse->GetQ();
				g_pulses[k].time = impulse->GetT()*5 - gOMtimeCal[impulse->GetNch()];
			}
			OMIDAlreadyIn = true;
			break;
		}
	}
	return OMIDAlreadyIn;
}

bool QFilterPassed(BExtractedImpulseTel* impulseTel)
{
	g_pulses.clear();
	int nPulsesWithHigherCharge = 0;
	for (int i = 0; i < impulseTel->GetNimpulse(); ++i)
	{
		h_chargePerHit->Fill(impulseTel->GetQ(i));
		// cout << i << " " << impulseTel->GetNch(i) << " " << impulseTel->GetQ(i) << " " << impulseTel->GetT(i) << endl;
		// h_nHitOM->Fill(impulseTel->GetNch(i));
		if (impulseTel->GetQ(i) > gQCut && gOMtimeCal[impulseTel->GetNch(i)] != 0)
		{
			if (OMIDAlreadyInGPulses(impulseTel->At(i)))
				continue;
			nPulsesWithHigherCharge++;
			g_pulses.push_back(PulsesVariables{impulseTel->GetNch(i),5*(impulseTel->GetT(i))-gOMtimeCal[impulseTel->GetNch(i)],impulseTel->GetQ(i)});
		}
	}
	h_nHitsAfterQ->Fill(nPulsesWithHigherCharge);
	if (nPulsesWithHigherCharge >= gQCutHits)
		return true;
	else
		return false;
}

bool TFilterPassed(BExtractedImpulseTel* impulseTel, TVector3& matrixPosition, double& matrixTime, int &nPulses)
{
	nPulses = 0;
	g_pulses.clear();

	for(int i = 0; i < impulseTel->GetNimpulse(); i++)
	{   
		if(impulseTel->GetQ(i) > 0 && gOMtimeCal[impulseTel->GetNch(i)] != 0)
		{
			double distanceFromLEDmatrix = (matrixPosition - gOMpositions[impulseTel->GetNch(i)]).Mag();
			double scattering_correction = (gOMpositions[impulseTel->GetNch(i)].Z() < matrixPosition.Z()) ? (gScatteringCorrection/15.0)*(matrixPosition.Z() - gOMpositions[impulseTel->GetNch(i)].Z()) : 0;
	        double expectedTime = matrixTime + distanceFromLEDmatrix*ReciprocalSpeedOfLightinWater + scattering_correction;

	        h_timeRes->Fill(5*(impulseTel->GetT(i)) - gOMtimeCal[impulseTel->GetNch(i)] - expectedTime);
	        if(5*(impulseTel->GetT(i)) - gOMtimeCal[impulseTel->GetNch(i)] >= expectedTime-gTCutTimeWindowNs && 5*(impulseTel->GetT(i)) - gOMtimeCal[impulseTel->GetNch(i)] <= expectedTime + gTCutTimeWindowNs)
	        {
	        	if (OMIDAlreadyInGPulses(impulseTel->At(i)))
	        		continue;
	        	g_pulses.push_back(PulsesVariables{impulseTel->GetNch(i),impulseTel->GetT(i)*5 - gOMtimeCal[impulseTel->GetNch(i)],impulseTel->GetQ(i)});
	        	nPulses++;
	        }
	    }
	}
	h_nHitsAfterT->Fill(nPulses);
	if (nPulses >= gTCutHits)
		return true;
	else
		return false;
}

void EstimateInitialMatrixPosition(TVector3& position, double& matrixTime)
{
	double maxCharge = numeric_limits<double>::min();
	int maxChargePulseID = -1;
	for (int i = 0; i < g_pulses.size(); ++i)
	{
		if (g_pulses[i].charge > maxCharge)
		{
			maxCharge = g_pulses[i].charge;
			maxChargePulseID = i;
		}
	}
	position = gOMpositions[g_pulses[maxChargePulseID].OMID]; 
	matrixTime = g_pulses[maxChargePulseID].time;
	// position.Print();
	// cout << matrixTime << endl;
}

int EventVisualization(BExtractedImpulseTel* impulseTel, TVector3& position, double& matrixTime, int eventID)
{
	int nHitsPerString[gNStrings]{0};
	TGraph* g_hits[gNStrings];
	TGraph* g_lowerLimit[gNStrings];
	TGraph* g_upperLimit[gNStrings];
	TMultiGraph* mg_hitsMatrix[gNStrings];
	TGraph* g_ledMatrix[gNStrings];
	TGraph* g_QvsL = new TGraph(g_pulses.size());
	g_QvsL->SetMarkerStyle(20);

	int nOMsPerString = gNOMs/gNStrings;

	for(int k = 0; k < impulseTel->GetNimpulse(); k++)
	{
		nHitsPerString[impulseTel->GetNch(k)/36]++;
	}
	for (int i = 0; i < gNStrings; ++i)
	{
		g_hits[i] = new TGraph(nHitsPerString[i]);
		g_lowerLimit[i] = new TGraph(nOMsPerString);
		g_upperLimit[i] = new TGraph(nOMsPerString);
		mg_hitsMatrix[i] = new TMultiGraph(Form("mg_%d",i),Form("String_%d;Calibrated time [ns]; OM Z position [m]",i+1));
		nHitsPerString[i] = 0;
		g_ledMatrix[i] = new TGraph(1);
		g_ledMatrix[i]->SetPoint(0,matrixTime,position.Z());
	}
	for(int k = 0; k < impulseTel->GetNimpulse(); k++)
	{
		int stringID = impulseTel->GetNch(k)/36;
		g_hits[stringID]->SetPoint(nHitsPerString[stringID],impulseTel->GetT(k)*5- gOMtimeCal[impulseTel->GetNch(k)],gOMpositions[impulseTel->GetNch(k)].Z());
		nHitsPerString[stringID]++;
	}
	for (int j = 0; j < gNOMs; ++j)
	{
		double distanceFromLEDmatrix = (position - gOMpositions[j]).Mag();
	    double expectedTime = matrixTime + distanceFromLEDmatrix*ReciprocalSpeedOfLightinWater;
		g_lowerLimit[j/nOMsPerString]->SetPoint(j%nOMsPerString,expectedTime-50,gOMpositions[j].Z());
		g_upperLimit[j/nOMsPerString]->SetPoint(j%nOMsPerString,expectedTime+50,gOMpositions[j].Z());
	}
	for (int i = 0; i < g_pulses.size(); ++i)
	{
		double distanceFromLEDmatrix = (position - gOMpositions[g_pulses[i].OMID]).Mag();
		g_QvsL->SetPoint(i,distanceFromLEDmatrix,g_pulses[i].charge);
	}

	TCanvas *cEvent = new TCanvas(Form("cEvent_%d",eventID),Form("cEvent_%d",eventID),200,10,600,400);  
	cEvent->Divide(3,3);
	for (int i = 0; i < gNStrings; ++i)
	{
		cEvent->cd(i+1);
		mg_hitsMatrix[i]->Add(g_hits[i],"P");
		mg_hitsMatrix[i]->Add(g_ledMatrix[i],"P");
		mg_hitsMatrix[i]->Add(g_lowerLimit[i],"L");
		mg_hitsMatrix[i]->Add(g_upperLimit[i],"L");
		mg_hitsMatrix[i]->Draw("AP");
		// mg_hitsMatrix[i]->SetBit(kCanDelete);
		g_hits[i]->SetMarkerStyle(20);
		g_ledMatrix[i]->SetMarkerStyle(20);
		g_ledMatrix[i]->SetMarkerColor(kRed);
		g_lowerLimit[i]->SetLineColor(kGreen);
		g_upperLimit[i]->SetLineColor(kGreen);
	}
	cEvent->cd(9);
	g_QvsL->Draw("AP");
	g_QvsL->SetTitle("Charge vs. Distance;Distance from cascade [m]; Charge [FADC channels]");
	cEvent->Write();

	delete cEvent;
	delete g_QvsL;
	
	for (int i = 0; i < gNStrings; ++i)
	{
		delete mg_hitsMatrix[i];
	}
	return 0;
}

int ChargeVisualization(int eventID)
{
	int nHitsPerString[gNStrings]{0};
	TGraph* g_MeasQ[gNStrings];
	TGraph* g_ExpQ[gNStrings];
	TMultiGraph* mg_QSum[gNStrings];

	int nOMsPerString = gNOMs/gNStrings;

	for(int k = 0; k < g_pulses.size(); k++)
	{
		nHitsPerString[g_pulses[k].OMID/36]++;
	}
	for (int i = 0; i < gNStrings; ++i)
	{
		g_MeasQ[i] = new TGraph(nHitsPerString[i]);
		g_ExpQ[i] = new TGraph(nHitsPerString[i]);
		mg_QSum[i] = new TMultiGraph(Form("mg_%d",i),Form("String_%d; OM ID [#]; Charge [p.e.]",i+1));
		nHitsPerString[i] = 0;
	}
	for(int k = 0; k < g_pulses.size(); k++)
	{
		int stringID = g_pulses[k].OMID/36;
		g_MeasQ[stringID]->SetPoint(nHitsPerString[stringID],g_pulses[k].OMID,g_pulses[k].charge/gOMQCal[g_pulses[k].OMID]);
		g_ExpQ[stringID]->SetPoint(nHitsPerString[stringID],g_pulses[k].OMID,g_pulses[k].expectedCharge);
		nHitsPerString[stringID]++;
	}

	TCanvas *cCharge = new TCanvas(Form("cCharge_%d",eventID),Form("cCharge_%d",eventID),200,10,600,400);  
	cCharge->Divide(3,3);
	for (int i = 0; i < gNStrings; ++i)
	{
		cCharge->cd(i+1);
		mg_QSum[i]->Add(g_MeasQ[i],"P");
		mg_QSum[i]->Add(g_ExpQ[i],"P");
		mg_QSum[i]->Draw("AP");
		g_MeasQ[i]->SetMarkerStyle(20);
		g_MeasQ[i]->SetMarkerColor(kBlue);
		g_ExpQ[i]->SetMarkerStyle(20);
		g_ExpQ[i]->SetMarkerColor(kGreen);
	}
	cCharge->Write();

	delete cCharge;
	
	for (int i = 0; i < gNStrings; ++i)
	{
		delete mg_QSum[i];
	}
	return 0;
}

bool ZFilterPassed(TVector3& position)
{
	if (position.Z() < gZCut)
		return true;
	else
		return false;
}

bool TDelayFilterPassed(double matrixTime)
{
	double smallestDelay = numeric_limits<double>::max();;
	for (int i = 0; i < g_pulses.size(); ++i)
	{
		if ((g_pulses[i].time - matrixTime) < gTDelayCut)
			return true;
	}
	return false;
}

bool QRatioFilterPassed(BExtractedImpulseTel* impulseTel, double &qRatio)
{
	double allCharge = 0;
	for (int i = 0; i < impulseTel->GetNimpulse(); ++i)
	{
		allCharge += impulseTel->GetQ(i)>0?impulseTel->GetQ(i):0;
	}
	double insideCharge = 0;
	for (int i = 0; i < g_pulses.size(); ++i)
	{
		insideCharge += g_pulses[i].charge>0?g_pulses[i].charge:0;
	}
	// PrintGPulses();
	// cout << "allCharge: " << allCharge << " insideCharge: " << insideCharge << " Ratio: " << (double)insideCharge/allCharge*100 << endl;
	qRatio = (double)insideCharge/allCharge*100;
	if (qRatio > gQRatioCut)
		return true;
	else
		return false;
}

bool BranchFilterPassed(TVector3& position)
{
	int upper = 0;
	int lower = 0;
	for (int i = 0; i < g_pulses.size(); ++i)
	{
		if (gOMpositions[g_pulses[i].OMID].Z() > position.Z())
			upper++;
		else
			lower++;
	}
	// cout << upper << " " << lower << endl;
	if (upper+gBranchCut > lower)
		return true;
	else
		return false;
}

struct OMIDdistance
{
	int OMID;
	double distance;
};

bool compareEntries(OMIDdistance first, OMIDdistance second)
{
	return (first.distance <= second.distance);
}

bool CloseHitsFilterPassed(TVector3& matrixPosition, int &closeHits)
{
	closeHits = 0;

	vector<OMIDdistance> vOMIDdistance;
	// matrixPosition.Print();

	for (int i = 0; i < gNOMs; ++i)
	{
		if (gOMQCal[i] != -1 && gOMpositions[i].Z() > matrixPosition.Z()-20 && gOMtimeCal[i] != 0)
		{
			// cout << i << endl;
			// gOMpositions[i].Print();
			double OMdistance = (matrixPosition - gOMpositions[i]).Mag();
			vOMIDdistance.push_back(OMIDdistance{i,OMdistance});
			// cout << OMdistance <<endl;
		}
	}
	sort(vOMIDdistance.begin(),vOMIDdistance.end(),compareEntries);

	for (int i = 0; i < 15; ++i)
	{
		// cout << i << " " << vOMIDdistance[i].OMID << " " << vOMIDdistance[i].distance << endl;
		for (int j = 0; j < g_pulses.size(); ++j)
		{
			if (g_pulses[j].OMID == vOMIDdistance[i].OMID)
			{
				closeHits++;
				break;
			}
		}
	}

	// PrintGPulses();
	// cout << "closeHits: " << closeHits << endl;

	if (closeHits > 8)
		return true;
	else
		return false;
}

bool LikelihoodFilterPassed(TVector3 &matrixPosition, double &matrixTime, double &energy, double &theta, double &phi, double &likelihood)
{
	double lowestLog = 10000;
	double cascadeEnergy = 0;
	for (int k = 0; k < 12; ++k)
	{
		for (int l = 0; l < 8; ++l)
		{
			double cascadeTheta = TMath::Pi()/12*k;
			double cascadePhi = TMath::Pi()*2/8*l;
			double recentLog = FitMatrixDirection(matrixPosition,matrixTime,cascadeEnergy,cascadeTheta,cascadePhi);
			// double recentLog = 100000;
			if (recentLog < lowestLog)
			{
				lowestLog = recentLog;
				energy = cascadeEnergy;
				theta = cascadeTheta;
				phi = cascadePhi;
			}
			// cout << k << " " << l << " " <<  << " " << cascadeEnergy << " " << cascadeTheta << " " << cascadePhi << endl;
		}
	}
	likelihood = lowestLog;
	if (lowestLog < gLikelihoodCut)
		return true;
		// cout << i << " " << energy << " " << theta << " " << phi << endl;
	return false;
}

void FillCascPos(TVector3 matrixPosition)
{
	g_cascadePosition->Set(g_cascadePosition->GetN()+1);
	double rho = sqrt(pow(matrixPosition.X()-gOMpositions[287].X(),2)+pow(matrixPosition.Y()-gOMpositions[287].Y(),2));
	g_cascadePosition->SetPoint(g_cascadePosition->GetN()-1,rho,matrixPosition.Z());
	// cout << "X: " << matrixPosition.X()-gOMpositions[287].X() << " Y: " << matrixPosition.Y()-gOMpositions[287].Y() << " Z: " << matrixPosition.Z()<< endl;
	// cout << "Rho: " << rho << endl;
	h_XYposition->Fill(matrixPosition.X()-gOMpositions[287].X(),matrixPosition.Y()-gOMpositions[287].Y());
}

int DoTheMagic(TTree* tree, BExtractedImpulseTel* impulseTel)
{
	// if --view (-w) switch is used, only specified eventID is visualized and program stops
	if (gVisEventID != -1)
	{
		TString outputFileName = BARS::Data::Directory(BARS::App::Cluster, BARS::App::Season, BARS::App::Run, BARS::Data::JOINT, gProductionID.c_str());
		outputFileName += Form("recCasc_vis_%d.root",gVisEventID);
		TFile* outputFile = new TFile(outputFileName,"RECREATE");
		tree->GetEntry(gVisEventID);
		TVector3 matrixPosition(0,0,0); 
		double matrixTime = 0;
		EventVisualization(impulseTel,matrixPosition,matrixTime,gVisEventID);
		cout << "Visualization of event: " << gVisEventID << " has been produced!" << endl;
		return 0;
	}

	TString outputFileName2 = BARS::Data::Directory(BARS::App::Cluster, BARS::App::Season, BARS::App::Run, BARS::Data::JOINT, gProductionID.c_str());
	outputFileName2 += "recCasc_nTuple.root";
	TFile* outputFile2 = new TFile(outputFileName2,"RECREATE");

   	TNtuple* nt_cascades = new TNtuple("nt_cascades","NTuple of reconstructed cascades","runID:EventID:NPulses:NPulsesT:QRatio:CloseHits:Likelihood:X:Y:Z:Time:Energy:Theta:Phi");

	TString outputFileName = BARS::Data::Directory(BARS::App::Cluster, BARS::App::Season, BARS::App::Run, BARS::Data::JOINT, gProductionID.c_str());
	outputFileName += "recCasc_vis.root";
	TFile* outputFile = new TFile(outputFileName,"RECREATE");

	int nEntries = (gNEventsProcessed == -1)?tree->GetEntries():gNEventsProcessed;
	int nNFilter = 0;
	int nQFilter = 0;
	int nQFilterChi2 = 0;
	int nTFilter = 0;
	int nTFilterChi2 = 0;
	int nZFilter = 0;
	int nTDelayFilter = 0;
	int nQRatioFilter = 0;
	int nBranchFilter = 0;
	int nCloseHitsFilter = 0;
	int nLikelihoodFilter = 0;
	
	int nPulses = 0;
	int nPulsesT = 0;
	double qRatio = 0;
	int closeHits = 0;

	// Loop through all the events
	for (int i = 0; i < nEntries; ++i)
	{
		if (i%(nEntries/10) == 0)
		{
			cout << round((double)(i)/nEntries*100) << "% ";
			cout << std::flush;
		}
		tree->GetEntry(i);
		if (!NFilterPassed(impulseTel,nPulses))
			continue;
		nNFilter++;
		if (!QFilterPassed(impulseTel))
			continue;
		nQFilter++;

		//initialization of R and T for LED matrix
		TVector3 matrixPosition(0,0,0); 
		double matrixTime = 0;
		EstimateInitialMatrixPosition(matrixPosition,matrixTime);
		fMinuit->SetFCN(MEstimator);
		double chi2QResult = FitMatrixPosition(matrixPosition,matrixTime);
		h_chi2AfterQ->Fill(chi2QResult);
		if (chi2QResult > gQCutChi2)
			continue;
		nQFilterChi2++;
		if (!TFilterPassed(impulseTel,matrixPosition,matrixTime,nPulsesT))
			continue;
		nTFilter++;
		fMinuit->SetFCN(chi2);
		double chi2TResult = FitMatrixPosition(matrixPosition,matrixTime);
		h_chi2AfterT->Fill(chi2TResult);
		if (chi2TResult > gTCutChi2)
			continue;
		nTFilterChi2++;
		TFilterPassed(impulseTel,matrixPosition,matrixTime,nPulsesT);
		if (!ZFilterPassed(matrixPosition))
			continue;
		nZFilter++;
		FillCascPos(matrixPosition);
		if (!TDelayFilterPassed(matrixTime))
			continue;
		nTDelayFilter++;
		if (!QRatioFilterPassed(impulseTel,qRatio))
			continue;
		nQRatioFilter++;
		if (!BranchFilterPassed(matrixPosition))
			continue;
		nBranchFilter++;
		if (!CloseHitsFilterPassed(matrixPosition,closeHits))
			continue;
		nCloseHitsFilter++;
		EventVisualization(impulseTel,matrixPosition,matrixTime,i);
		if (nCloseHitsFilter > 100)
		{
			cout << "File was identified as a matrixRun and terminated" << endl;
			break;
		}
		double cascadeEnergy = 0;
		double cascadeTheta = 0;
		double cascadePhi = 0;
		double likelihood = 0;
		if (!LikelihoodFilterPassed(matrixPosition,matrixTime,cascadeEnergy,cascadeTheta,cascadePhi,likelihood))
		{
			nt_cascades->Fill((double)BARS::App::Run,(double)i,(double)nPulses,(double)nPulsesT,qRatio,(double)closeHits,-1,matrixPosition.X(),matrixPosition.Y(),matrixPosition.Z(),matrixTime,0,0,0);
			continue;
		}
		nLikelihoodFilter++;
		nt_cascades->Fill((double)BARS::App::Run,(double)i,(double)nPulses,(double)nPulsesT,qRatio,(double)closeHits,likelihood,matrixPosition.X(),matrixPosition.Y(),matrixPosition.Z(),matrixTime,cascadeEnergy,cascadeTheta,cascadePhi);
		FillExpectedCharges(matrixPosition.X(),matrixPosition.Y(),matrixPosition.Z(),cascadeEnergy,cascadeTheta,cascadePhi);
		ChargeVisualization(i);
	}

	cout << "\nNentries: \t" << nEntries << endl;
	cout << "After NFilter: \t" << nNFilter << endl;
	cout << "After QFilter: \t" << nQFilter << endl;
	cout << "After QFilterChi2: \t" << nQFilterChi2 << endl;
	cout << "After TFilter: \t" << nTFilter << endl;
	cout << "After TFilterChi2: \t" << nTFilterChi2 << endl;
	cout << "After ZFilter: \t" << nZFilter << endl;
	cout << "After TDelayFilter: \t" << nTDelayFilter << endl;
	cout << "After QRatioFilter: \t" << nQRatioFilter << endl;
	cout << "After BranchFilter: \t" << nBranchFilter << endl;
	cout << "After CloseHitsFilter: \t" << nCloseHitsFilter << endl;
	cout << "After LikelihoodFiter: \t" << nLikelihoodFilter << endl;

	outputFile->Close();
	outputFile2->cd();

	if (nCloseHitsFilter < 100)
		nt_cascades->Write();

	outputFile2->Close();

	// delete outputFile;

	return 0;
}


int SetOMsDynamic(BGeomTel* bgeom) //dynamic posiions
{
	int nOKOMs = 0;

	for (int i = 0; i < bgeom->GetNumOMs(); ++i)
	{
		if (bgeom->At(i)->HasData())
		{   
			gOMpositions[i] = bgeom->At(i)->GetPosition();
			// cout << i << " " << gOMpositions[i].X() << " " << gOMpositions[i].Y() << endl;
			nOKOMs++;
		}
	}
	return nOKOMs;
}

int ReadGeometry(BExtractedHeader* header) // read dynamic geometry
{
	const char* geometryFileName = BARS::Geom::File(BARS::App::Cluster, BARS::App::Season, BARS::Geom::OM_EXPORT_LINEAR);

	TTree* geometryTree = nullptr;
	TFile* geomFile = new TFile(geometryFileName,"READ");
	geometryTree = (TTree*)geomFile->Get("Events");

 	// check if the file is opened
	if (!geomFile || !geometryTree)
	{
		cout<<"No TTree called Event in the geometry was found!"<<endl;
		return 99;
	}

	long double startTime = header->GetTime().GetSec();

	BDynamicGeometry* telGeometry = NULL;
	geometryTree->SetBranchAddress("BGeomTel.", &telGeometry);

    //extract time of the first and the last geometry record
	geometryTree->GetEntry(0);
	double geometryStartTime = telGeometry->GetTime().GetSec();
	geometryTree->GetEntry(geometryTree->GetEntries()-1);
	double geometryEndTime = telGeometry->GetTime().GetSec();


  	//check if the geometry file covers the whole run
  	if (geometryStartTime > startTime)
  	{
  		geometryTree->GetEntry(0);
  		int nOKOMs = SetOMsDynamic(telGeometry);
    	cerr<< "The precise dynamic geometry for this run was not available (startGeometry > startRun)" << endl;
    	cerr << "StartGeometry: " << geometryStartTime << " startRun: " << startTime << endl;
    	cerr<< "The first accessible detector geometry is used. The time difference: " << (geometryStartTime-startTime)/3600.0/24.0 << " days." << endl;
    	return 0;
  	}
  	if (geometryEndTime < startTime)
  	{
  		geometryTree->GetEntry(geometryTree->GetEntries()-1);
  		int nOKOMs = SetOMsDynamic(telGeometry);
    	cerr<< "The precise dynamic geometry for this run was not available (endGeometry < startRun)" << endl;
    	cerr<< "The last accessible detector geometry is used. The time difference: " << (startTime-geometryEndTime)/3600.0/24.0 << " days." << endl;
    	return 0;
  	}
  	
  	// iterate through all the geometry entries
	for (int i = 0; i < geometryTree->GetEntries(); ++i)
	{
		geometryTree->GetEntry(i);   

	    if(telGeometry->GetTime().AsDouble() >= startTime) //to druhe cislo je zaciatok nasho subrunu
	    {
	    	// cout << "SubRun start: " << telGeometry->GetTime().AsDouble()<<endl;
	    	int nOKOMs = SetOMsDynamic(telGeometry);
	    	return 0;
	    }  
	} 
	return -1;
}

int ReadCalFile(const char* fileName, double multConst)
{
	ifstream inputFile;
    inputFile.open(fileName);

    if (!inputFile)
    {
    	cerr << "Calibration file: " << fileName << " was NOT found. Program termination!" << endl;
    	return -1;
  	}

    string dummyLine;

    getline(inputFile,dummyLine);
    getline(inputFile,dummyLine);

    double readValue = 0;

    for (int i = 0; i < gNOMs; ++i)
    {
		if (i%36 == 0)
		{
		    getline(inputFile,dummyLine);
		}
		inputFile >> readValue;
		if (readValue != -10000 && readValue != -1000)
			gOMtimeCal[i] += readValue*multConst;
		else
			gOMtimeCal[i] = 0;
		// cout << i << " " << gOMtimeCal[i] << endl;
		if ((i+1)%36 == 0)
		{
		    getline(inputFile,dummyLine);
		    getline(inputFile,dummyLine);
		}
    }
    inputFile.close();
    return 0;
}

int ReadTimeCal(void)
{
	const char* offsetFileName = BARS::Calib::File(BARS::App::Cluster, BARS::App::Season, BARS::Calib::OFFSET, "offsets");
	if (ReadCalFile(offsetFileName,-2.5) == -1)
		return -1;
	const char* timeCalName = BARS::Calib::File(BARS::App::Cluster, BARS::App::Season, BARS::Calib::TIME, "timecalib_dzh");
	if (ReadCalFile(timeCalName,1) == -1)
		return -1;
	// for (int i = 0; i < gNOMs; ++i)
	// {
	// 	cout << gOMtimeCal[i] << endl;
	// }
    return 0;
}

// The charge calibration parameters are read directly from qcalib file created by Zhenya's DQM in the same folder
// It also let us know about operating OMs. If there is -1 in the qcalib for OM we assume it is dead.
int ReadQCal(void)
{
	const char* filePath = BARS::Data::File(BARS::App::Cluster, BARS::App::Season, BARS::App::Run, BARS::Data::JOINT,gProductionID.c_str());
	TString jointFileName(filePath);
	TString chargeFileName(jointFileName(0,jointFileName.Length()-17));
	chargeFileName += "qcalib";

	ifstream inputFile;
    inputFile.open(chargeFileName);

    if (!inputFile)
    {
    	cerr << "Calibration file: " << chargeFileName << " was NOT found. Program termination!" << endl;
    	return -1;
  	}

    string dummyLine;
    double readValue = 0;

    for (int i = 0; i < gNOMs; ++i)
    {
		inputFile >> readValue;
		gOMQCal[i] = readValue;
		if ((i+1)%36 == 0)
		{
		    getline(inputFile,dummyLine);
		    getline(inputFile,dummyLine);
		}
    }
    inputFile.close();
}

void SaveHistograms()
{
	TString outputFileName = BARS::Data::Directory(BARS::App::Cluster, BARS::App::Season, BARS::App::Run, BARS::Data::JOINT, gProductionID.c_str());
	outputFileName += "recCasc_hist.root";
	TFile* outputFile = new TFile(outputFileName,"RECREATE");
	h_nHits->Write();
	h_chargePerHit->Write();
	h_nHitsAfterQ->Write();
	h_chi2AfterQ->Write();
	h_nHitsAfterT->Write();
	h_chi2AfterT->Write();
	g_cascadePosition->Write();
	h_timeRes->Write();
	h_nHitOM->Write();
	h_XYposition->Write();
	delete outputFile;
}

int ReadLogTable()
{
	cout << "4D LogTable reading starts" << endl;
	ifstream fTab ("/Data/BaikalData/showerTable/hq001200_double.dqho2011", ios::in|ios::binary|ios::ate);

	if (!fTab.is_open())
		return -1;

	streampos size = 8;
	char * memblock = new char[size];

	fTab.seekg (0, ios::beg);
	fTab.read(memblock,size);
	int dummy;

	for(int i = 0; i < 200; i++){           // step of R in meters
	  for(int j = 0; j < 351; j++){         // step of Z in meters
	      for(int k = 0; k < 21; k++){      // step of phi
	        for(int m = 0; m < 21; m++){     // step of cos(theta)
				fTab.read(memblock,size);
				double* value = (double*)memblock;
				gLogTable4D[i][j][k][m] = *value;
				// cout << i << " " << j << " " << k << " " << m << " " << gLogTable4D[i][j][k][m] << endl;
	        }
	        // cin >> dummy;
	      }
	  }
	}
	fTab.close();
	cout << "LogTable ends" << endl;
	return 0;
}

int main(int argc, char** argv) 
{
    // Init should be called at the beggining of all BARS programms
    App::Init(argc, argv, 0, parseOpts, readRC, checkParams);
    SetFitter();

    if (!CheckParams)
    	return -1;

    if (gProductionID == "")
    	gProductionID = "barsv051";

    // const char* filePath = BARS::Data::File(BARS::App::Cluster, BARS::App::Season, BARS::App::Run, BARS::Data::JOINT,"r01_i01_j01_t01");
    const char* filePath = BARS::Data::File(BARS::App::Cluster, BARS::App::Season, BARS::App::Run, BARS::Data::JOINT,gProductionID.c_str());

    if (ReadLogTable() == -1)
    {
    	std::cout << "Problem with 4D LogLikelihood file!" << std::endl;
    	return -1;
    }
    // ReadLogTable();

    if (!BARS::App::FileExists(filePath))
    {
    	std::cout << "File: " << filePath << " was not found!" << endl;
    	return -2;
    }

    TFile* file = new TFile(filePath,"READ");
    TTree* tree = (TTree*)file->Get("Events");

	BExtractedImpulseTel* impulseTel = NULL;
	BExtractedHeader* header = NULL;
	BGeomTel* bgeom = NULL;
	tree->SetBranchAddress("BJointImpulseTel.",&impulseTel);
	tree->SetBranchAddress("BJointHeader.",&header);  

	PrintConfig();
	cout << "Season: " << BARS::App::Season << " Cluster: " << BARS::App::Cluster << " Run: " <<  BARS::App::Run << endl;
	tree->GetEntry(5);
	if (ReadGeometry(header) == -1)
	{
		std::cout << "Problem with geometry file!" << std::endl;
		return -1;
	}
	if (ReadTimeCal() == -1)
	{
		std::cout << "Problem with time calibration files!" << std::endl;
		return -1;
	}
	if (ReadQCal() == -1)
	{
		std::cout << "Problem with charge calibration files!" << std::endl;
		return -1;	
	}

	DoTheMagic(tree,impulseTel);
	SaveHistograms();

    return 0;
}
