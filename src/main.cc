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
#include "TChain.h"
#include "TRandom2.h"

//BARS dependencies
#include "BARS.h"
#include "BExtractedImpulseTel.h"
#include "BExtractedHeader.h"
#include "BGeomTel.h"
#include "BEvoGeomApply.h"
#include "BEvent.h"
#include "BImpulse.h"
#include "BGeomTelMC.h"
#include "BEventMask.h"
#include "BSource.h"

using namespace std;

TVirtualFitter* fMinuit;

const int gNOMs = 288;
const int gNStrings = 8;
const double ReciprocalSpeedOfLightinWater = 4.57;

vector<TVector3> gOMpositions(288);
vector<double> gOMtimeCal(288,-1);
vector<double> gOMQCal(288,0);

double gLogTable4D[200][351][21][21]{0};

struct PulsesVariables
{
  int OMID;
  double time;
  double charge;
  double expectedCharge;
  int MCflag;
};

struct mcCascade
{
	int eventID, nHits, nNoiseHits;
	float showerEnergy, cosTheta, phi, position[3], charge[288], time[288], noiseCharge[1000], noiseTime[1000];
	unsigned short chID[288], noiseChID[1000];
};

vector<PulsesVariables> g_pulses; //global vector of PulsesVariables

void PrintGPulses()
{
	for (int i = 0; i < g_pulses.size(); ++i)
	{
		cout << i << " " << g_pulses[i].OMID << " " << g_pulses[i].charge << " " << g_pulses[i].time << " " << g_pulses[i].MCflag << " " << gOMpositions[g_pulses[i].OMID].X() << " " << gOMpositions[g_pulses[i].OMID].Y() << " " << gOMpositions[g_pulses[i].OMID].Z() << endl;
	}
}

using namespace BARS;

TH1F* h_nHits = new TH1F("h_nHits","Number of hits;Hits [#]; NoE [#]",200,0,200);
TH1F* h_chargePerHit = new TH1F("h_chargePerHit","Charge per hit;Q [p.e.];NoE [#]",1000,0,1000);
TH1F* h_chargePerEvent = new TH1F("h_chargePerEvent","Charge per event;Q [p.e.];NoE [#]",1000,0,10000);
TH1F* h_nHitsAfterQ = new TH1F("h_nHitsAfterQ","Number of hits after Q Filter; Nhits [#];NoE [#]",20,0,20);
TH1F* h_chi2AfterQ = new TH1F("h_chi2AfterQ","Chi2 after Q Filter; chi^2 [#];NoE [#]",1000,0,10000);
TH1F* h_nHitsAfterT = new TH1F("h_nHitsAfterT","Number of hits after T Filter; NHits [#];NoE [#]",100,0,100);
TH1F* h_chi2AfterT= new TH1F("h_chi2AfterT","Chi2 after T Filter; chi^2 [#];NoE [#]",100,0,100);
TGraph* g_cascadePosition = new TGraph();
TH1F* h_timeRes = new TH1F("h_timeRes","Time residuals;DeltaT [ns];NoE [#]",400,-200,200);
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
	if (inputs[2] < 0 || inputs[2] > TMath::Pi())
		return -1;
	if (inputs[3] < -1 || inputs[3] > 1)
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
	outputs[3] = (int)floor((1-inputs[3])/0.1);
	outputs[2] = (int)floor(inputs[2]/(TMath::Pi()/20));

	return 0;
}

double interpolate(double x, double x1, double x2, double f1, double f2)
{
	return (x-x1)/(x2-x1)*f2 + (x2-x)/(x2-x1)*f1;
}

double interpolate4D(double* inputs, int* indexes)
{
	double R0Z0Phi0 = interpolate(inputs[3],1-(indexes[3])*0.1,1-(indexes[3]+1)*0.1,gLogTable4D[indexes[0]][indexes[1]][indexes[2]][indexes[3]],gLogTable4D[indexes[0]][indexes[1]][indexes[2]][indexes[3]+1]);
	double R0Z0Phi1 = interpolate(inputs[3],1-(indexes[3])*0.1,1-(indexes[3]+1)*0.1,gLogTable4D[indexes[0]][indexes[1]][indexes[2]+1][indexes[3]],gLogTable4D[indexes[0]][indexes[1]][indexes[2]+1][indexes[3]+1]); 

	double R0Z1Phi0 = interpolate(inputs[3],1-(indexes[3])*0.1,1-(indexes[3]+1)*0.1,gLogTable4D[indexes[0]][indexes[1]+1][indexes[2]][indexes[3]],gLogTable4D[indexes[0]][indexes[1]][indexes[2]][indexes[3]+1]);
	double R0Z1Phi1 = interpolate(inputs[3],1-(indexes[3])*0.1,1-(indexes[3]+1)*0.1,gLogTable4D[indexes[0]][indexes[1]+1][indexes[2]+1][indexes[3]],gLogTable4D[indexes[0]][indexes[1]][indexes[2]][indexes[3]+1]);

	double R1Z0Phi0 = interpolate(inputs[3],1-(indexes[3])*0.1,1-(indexes[3]+1)*0.1,gLogTable4D[indexes[0]+1][indexes[1]][indexes[2]][indexes[3]],gLogTable4D[indexes[0]][indexes[1]][indexes[2]][indexes[3]+1]);
	double R1Z0Phi1 = interpolate(inputs[3],1-(indexes[3])*0.1,1-(indexes[3]+1)*0.1,gLogTable4D[indexes[0]+1][indexes[1]][indexes[2]+1][indexes[3]],gLogTable4D[indexes[0]][indexes[1]][indexes[2]][indexes[3]+1]);

	double R1Z1Phi0 = interpolate(inputs[3],1-(indexes[3])*0.1,1-(indexes[3]+1)*0.1,gLogTable4D[indexes[0]+1][indexes[1]+1][indexes[2]][indexes[3]],gLogTable4D[indexes[0]][indexes[1]][indexes[2]][indexes[3]+1]);
	double R1Z1Phi1 = interpolate(inputs[3],1-(indexes[3])*0.1,1-(indexes[3]+1)*0.1,gLogTable4D[indexes[0]+1][indexes[1]+1][indexes[2]+1][indexes[3]],gLogTable4D[indexes[0]][indexes[1]][indexes[2]][indexes[3]+1]);

	double R0Z0 = interpolate(inputs[2],(indexes[2])*(TMath::Pi()/20),(indexes[2]+1)*(TMath::Pi()/20),R0Z0Phi0,R0Z0Phi1);
	double R0Z1 = interpolate(inputs[2],(indexes[2])*(TMath::Pi()/20),(indexes[2]+1)*(TMath::Pi()/20),R0Z1Phi0,R0Z1Phi1);

	double R1Z0 = interpolate(inputs[2],(indexes[2])*(TMath::Pi()/20),(indexes[2]+1)*(TMath::Pi()/20),R1Z0Phi0,R1Z0Phi1);
	double R1Z1 = interpolate(inputs[2],(indexes[2])*(TMath::Pi()/20),(indexes[2]+1)*(TMath::Pi()/20),R1Z1Phi0,R1Z1Phi1);

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

  	double expectedTime = matrixTime + distanceFromLEDmatrix*ReciprocalSpeedOfLightinWater;// + scattering_correction; //v nanosekundach

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
		// theChi = constant*error*(g_pulses[i].time - ExpectedTime(par[0], par[1], par[2], par[3], g_pulses[i].OMID))*TMath::Log10(g_pulses[i].charge+1);
		theChi = constant*error*(g_pulses[i].time - ExpectedTime(par[0], par[1], par[2], par[3], g_pulses[i].OMID));
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
		double element = g_pulses[i].charge*TMath::Sqrt(1+TMath::Power(Tres,2)/2);
		sum += element;
	}
	f = sum; // function returns calculated chi2, it is global function which is used in SetFCN()

}

// cascadeParameters[7]: X,Y,Z,Time,Energy,Theta,Phi
// tableParameters[4]: R,Z,Phi',CosTheta'
void GetParameters(const Double_t* cascadeParameters,const int OMID, double* tableParameters)
{
	TVector3 OMpos(gOMpositions[OMID].X()-cascadeParameters[0],gOMpositions[OMID].Y()-cascadeParameters[1],gOMpositions[OMID].Z()-cascadeParameters[2]);
	// cout << gOMpositions[OMID].X() << " " << gOMpositions[OMID].Y() << " " << gOMpositions[OMID].Z() << endl;
	// cout << cascadeParameters[0] << " " << cascadeParameters[1] << " " << cascadeParameters[2] << endl;
	// OMpos.Print();

	TVector3 showerRef(0,0,1);
	showerRef.SetTheta(cascadeParameters[5] - TMath::Pi());
	showerRef.SetPhi((cascadeParameters[6] + TMath::Pi()));
	// showerRef.Print();
	// OMpos.Print();
	
	tableParameters[0] = OMpos.Perp(showerRef);
	tableParameters[1] = OMpos*showerRef;

	// cout << "New" << endl;
	TVector3 y = OMpos.Cross(showerRef);
	// y.Print();
	y.SetMag(1.0);
	TVector3 x = OMpos.Cross(y);
	// x.Print();
	x.SetMag(1.0);

	// tableParameters[3] = rhoProjection.Angle(omegaProjection);
	TVector3 OMorien(0,0,-1);
	tableParameters[3] = OMpos*OMorien/OMpos.Mag();

	TVector3 par = OMpos;
	// par.Print();
	par.SetMag(OMpos*OMorien/OMpos.Mag());
	TVector3 phi = OMorien - par;

	tableParameters[2] = phi.Angle(x);

}

bool NotInGPulses(int OMID)
{
	bool inGPulses = false;
	for (int i = 0; i < g_pulses.size(); ++i)
	{
		if (g_pulses[i].OMID == OMID)
		{
			inGPulses = true;
			break;
		}
	}
	return !inGPulses;
}

void logLikelihood(Int_t &npar, Double_t* gin, Double_t &f, Double_t* par, Int_t iflag)
{
	// cout << "In log" << endl;
	double logLike = 0;
	double tableParameters[4]{0};
	// cout << "Calculating logLike" << endl;
	// cout << par[0] << " " << par[1] << " " << par[2] << " " << par[3] << " " << par[4] << " " << par[5] << " " << par[6] << endl;
	for (int i = 0; i < g_pulses.size(); ++i)
	{
		// cout << "GPulses: " << g_pulses[i].OMID << endl;
		GetParameters(par,g_pulses[i].OMID,tableParameters);
		// gOMpositions[g_pulses[i].OMID].Print();
		// cout << i << " " << tableParameters[0] << " " << tableParameters[1] << " " << tableParameters[2] << " " << tableParameters[3] << endl;
		double expectedNPE = GetInterpolatedValue(tableParameters);
		// cout << expectedNPE << " " << expectedNPE*100000000*par[4] << " " << g_pulses[i].charge << " " << TMath::PoissonI(g_pulses[i].charge,expectedNPE*100000000*par[4]) << endl;
		
		if (TMath::Poisson(g_pulses[i].charge,expectedNPE*100000000*par[4]) > 10e-350)
		{
			logLike -= TMath::Log10(TMath::Poisson(g_pulses[i].charge,expectedNPE*100000000*par[4]));
			// cout << TMath::Log10(TMath::PoissonI(g_pulses[i].charge,expectedNPE*100000000*par[4])) << endl;
		}
		else
		{
			logLike -= -350;
			// cout << -350 << endl;
		}
	}
	// cout << "After g pulses" << endl;
	for (int i = 0; i < gNOMs; ++i)
	{
		if (NotInGPulses(i) && gOMQCal[i] != -1)
		{
			GetParameters(par,i,tableParameters);
			double expectedNPE = GetInterpolatedValue(tableParameters);
			// cout << i << " " << TMath::Poisson(0,expectedNPE*100000000*par[4]) << endl;

			if (TMath::Poisson(0,expectedNPE*100000000*par[4]) > 10e-350)
			{
				logLike -= TMath::Log10(TMath::Poisson(0,expectedNPE*100000000*par[4]));
				// cout << TMath::Log10(TMath::PoissonI(g_pulses[i].charge,expectedNPE*100000000*par[4])) << endl;
				// cout << i << " " << TMath::Log10(TMath::PoissonI(0,expectedNPE*100000000*par[4])) << endl;
			}
			else
			{
				logLike -= -350;
				// cout << i << " " << -350 << endl;
			}
		}
	}
	// cout << "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF: " << logLike << endl;
	f = logLike;
}

double FitMatrixPosition(TVector3 &R2, double &T2)
{ 
	fMinuit->ReleaseParameter(0);
	fMinuit->ReleaseParameter(1);
	fMinuit->ReleaseParameter(2);
	fMinuit->ReleaseParameter(3);

	//cout<<"fitEVent "<<endl;
	fMinuit->SetParameter(0,"LED X",R2.X(),1,-750,750); //estimation of start point
	fMinuit->SetParameter(1,"LED Y",R2.Y(),1,-750,750);
	fMinuit->SetParameter(2,"LED Z",R2.Z(),1,-750,750);
	if (!gMCMu && !gMCNu && !gMCCas)
		fMinuit->SetParameter(3,"Time",T2,1,14000,17000);
	else
		fMinuit->SetParameter(3,"Time",T2,1,-100,1600);

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
	fMinuit->SetFCN(chi2);
	// fMinuit->SetFCN(MEstimator);
}

double FitMatrixDirection(TVector3 &R2, double &T2, double &energy, double &theta, double &phi)
{ 
	// cout << "Fit Matrix direction" << endl;
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

	// cout << "FIT" << endl;
	double arglist[10] = {500,0.01}; //500 means 500 iterations and then MIGRAD stops
	fMinuit->ExecuteCommand("SIMPLEX",arglist,2); //this line performs minimalisation

	double chi2, edm, errdef;
	int nvpar, nparx;

	fMinuit->GetStats(chi2,edm,errdef,nvpar,nparx);
	energy = fMinuit->GetParameter(4);
	theta = fMinuit->GetParameter(5);
	phi = fMinuit->GetParameter(6);

	return chi2/(g_pulses.size()-4);  
	// return chi2/gNOMs;  
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

bool NFilterPassed(BEvent* event, int &nPulses)
{
	nPulses = event->NHits();
	h_nHits->Fill(nPulses);
	if (event->NHits() >= gNCut)
	{
		return true;
	}
	return false;
}

bool NFilterPassed(mcCascade* cascade, int &nPulses)
{
	nPulses = cascade->nHits+cascade->nNoiseHits;
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
			if (g_pulses[k].charge < impulse->GetQ()/gOMQCal[impulse->GetNch()])
			{
				g_pulses[k].charge = impulse->GetQ()/gOMQCal[impulse->GetNch()];
				g_pulses[k].time = impulse->GetT()*5 - gOMtimeCal[impulse->GetNch()];
			}
			OMIDAlreadyIn = true;
			break;
		}
	}
	return OMIDAlreadyIn;
}

bool OMIDAlreadyInGPulses(BImpulse* impulse)
{
	bool OMIDAlreadyIn = false;
	for (int k = 0; k < g_pulses.size(); ++k)
	{
		if (g_pulses[k].OMID == impulse->Channel())
		{
			if (g_pulses[k].charge < impulse->Q())
			{
				g_pulses[k].charge = impulse->Q();
				g_pulses[k].time = impulse->T();
			}
			OMIDAlreadyIn = true;
			break;
		}
	}
	return OMIDAlreadyIn;
}

bool QFilterPassed(BExtractedImpulseTel* impulseTel, double &overallCharge)
{
	g_pulses.clear();
	int nPulsesWithHigherCharge = 0;
	overallCharge = 0;
	for (int i = 0; i < impulseTel->GetNimpulse(); ++i)
	{
		if (gOMQCal[impulseTel->GetNch(i)] < 0)
			cout << "ERROR" << endl;
		h_chargePerHit->Fill(impulseTel->GetQ(i)/gOMQCal[impulseTel->GetNch(i)]);
		overallCharge += impulseTel->GetQ(i)/gOMQCal[impulseTel->GetNch(i)];
		// cout << i << " " << impulseTel->GetNch(i) << " " << impulseTel->GetQ(i) << " " << impulseTel->GetT(i) << endl;
		if (impulseTel->GetQ(i)/gOMQCal[impulseTel->GetNch(i)] > gQCut && gOMtimeCal[impulseTel->GetNch(i)] != 0)
		{
			if (OMIDAlreadyInGPulses(impulseTel->At(i)))
				continue;
			nPulsesWithHigherCharge++;
			g_pulses.push_back(PulsesVariables{impulseTel->GetNch(i),5*(impulseTel->GetT(i))-gOMtimeCal[impulseTel->GetNch(i)],impulseTel->GetQ(i)/gOMQCal[impulseTel->GetNch(i)]});
		}
	}
	h_nHitsAfterQ->Fill(nPulsesWithHigherCharge);
	h_chargePerEvent->Fill(overallCharge);
	if (nPulsesWithHigherCharge >= gQCutHits && overallCharge > gQCutOverall)
		return true;
	else
		return false;
}

bool QFilterPassed(BEvent* event)
{
	g_pulses.clear();
	int nPulsesWithHigherCharge = 0;
	for (int i = 0; i < event->NHits(); ++i)
	{
		h_chargePerHit->Fill(event->Q(i));
		// cout << i << " " << impulseTel->GetNch(i) << " " << impulseTel->GetQ(i) << " " << impulseTel->GetT(i) << endl;
		if (event->Q(i) > gQCut)
		{
			if (OMIDAlreadyInGPulses(event->GetImpulse(i)))
				continue;
			nPulsesWithHigherCharge++;
			g_pulses.push_back(PulsesVariables{event->HitChannel(i),event->T(i),event->Q(i)});
		}
	}
	h_nHitsAfterQ->Fill(nPulsesWithHigherCharge);
	if (nPulsesWithHigherCharge >= gQCutHits)
		return true;
	else
		return false;
}

bool QFilterPassed(mcCascade* cascade)
{
	g_pulses.clear();
	int nPulsesWithHigherCharge = 0;
	for (int i = 0; i < cascade->nHits; ++i)
	{
		h_chargePerHit->Fill(cascade->charge[cascade->chID[i]-1]);
		// cout << i << " " << impulseTel->GetNch(i) << " " << impulseTel->GetQ(i) << " " << impulseTel->GetT(i) << endl;
		if (cascade->charge[cascade->chID[i]-1] > gQCut)
		{
			// if (OMIDAlreadyInGPulses(event->GetImpulse(i)))
				// continue;
			nPulsesWithHigherCharge++;
			g_pulses.push_back(PulsesVariables{cascade->chID[i]-1,cascade->time[cascade->chID[i]-1],cascade->charge[cascade->chID[i]-1]});
		}
	}
	for (int i = 0; i < cascade->nNoiseHits; ++i)
	{
		h_chargePerHit->Fill(cascade->noiseCharge[i]);
		// cout << i << " " << impulseTel->GetNch(i) << " " << impulseTel->GetQ(i) << " " << impulseTel->GetT(i) << endl;
		if (cascade->noiseCharge[i] > gQCut)
		{
			// if (OMIDAlreadyInGPulses(event->GetImpulse(i)))
				// continue;
			nPulsesWithHigherCharge++;
			g_pulses.push_back(PulsesVariables{cascade->noiseChID[i],cascade->noiseTime[i],cascade->noiseCharge[i]});
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
	        double expectedTime = matrixTime + distanceFromLEDmatrix*ReciprocalSpeedOfLightinWater;// + scattering_correction;

	        h_timeRes->Fill(5*(impulseTel->GetT(i)) - gOMtimeCal[impulseTel->GetNch(i)] - expectedTime);
	        if(5*(impulseTel->GetT(i)) - gOMtimeCal[impulseTel->GetNch(i)] >= expectedTime-gTCutTimeWindowNs && 5*(impulseTel->GetT(i)) - gOMtimeCal[impulseTel->GetNch(i)] <= expectedTime + gTCutTimeWindowNs)
	        {
	        	if (OMIDAlreadyInGPulses(impulseTel->At(i)))
	        		continue;
	        	g_pulses.push_back(PulsesVariables{impulseTel->GetNch(i),impulseTel->GetT(i)*5 - gOMtimeCal[impulseTel->GetNch(i)],impulseTel->GetQ(i)/gOMQCal[impulseTel->GetNch(i)]});
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

bool TFilterPassed(BEvent* event, BEventMaskMC* eventMask, TVector3& matrixPosition, double& matrixTime, int &nPulses)
{
	nPulses = 0;
	g_pulses.clear();

	for(int i = 0; i < event->NHits(); i++)
	{   
		if(event->Q(i) > 0)
		{
			double distanceFromLEDmatrix = (matrixPosition - gOMpositions[event->HitChannel(i)]).Mag();
			// double scattering_correction = (gOMpositions[event->HitChannel(i)].Z() < matrixPosition.Z()) ? (gScatteringCorrection/15.0)*(matrixPosition.Z() - gOMpositions[event->HitChannel(i)].Z()) : 0;
	        double expectedTime = matrixTime + distanceFromLEDmatrix*ReciprocalSpeedOfLightinWater;// + scattering_correction;

	        h_timeRes->Fill(event->T(i) - expectedTime);
	        if(event->T(i) >= expectedTime-gTCutTimeWindowNs && event->T(i) <= expectedTime + gTCutTimeWindowNs)
	        {
	        	if (OMIDAlreadyInGPulses(event->GetImpulse(i)))
	        		continue;
	        	g_pulses.push_back(PulsesVariables{event->HitChannel(i),event->T(i),event->Q(i),0,eventMask->GetFlag(i)});
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

bool TFilterPassed(mcCascade* cascade, TVector3& matrixPosition, double& matrixTime, int &nPulses)
{
	nPulses = 0;
	g_pulses.clear();

	for(int i = 0; i < cascade->nHits; i++)
	{   
		if(cascade->charge[cascade->chID[i]-1] > 0)
		{
			double distanceFromLEDmatrix = (matrixPosition - gOMpositions[cascade->chID[i]-1]).Mag();
			// double scattering_correction = (gOMpositions[event->HitChannel(i)].Z() < matrixPosition.Z()) ? (gScatteringCorrection/15.0)*(matrixPosition.Z() - gOMpositions[event->HitChannel(i)].Z()) : 0;
	        double expectedTime = matrixTime + distanceFromLEDmatrix*ReciprocalSpeedOfLightinWater;// + scattering_correction;

	        h_timeRes->Fill(cascade->time[cascade->chID[i]-1] - expectedTime);
	        if(cascade->time[cascade->chID[i]-1] >= expectedTime-gTCutTimeWindowNs && cascade->time[cascade->chID[i]-1] <= expectedTime + gTCutTimeWindowNs)
	        {
	        	// if (OMIDAlreadyInGPulses(event->GetImpulse(i)))
	        		// continue;
	        	g_pulses.push_back(PulsesVariables{cascade->chID[i]-1,cascade->time[cascade->chID[i]-1],cascade->charge[cascade->chID[i]-1],0,0});
	        	nPulses++;
	        }
	    }
	}
	for(int i = 0; i < cascade->nNoiseHits; i++)
	{   
		if(cascade->noiseCharge[i] > 0)
		{
			double distanceFromLEDmatrix = (matrixPosition - gOMpositions[cascade->noiseChID[i]]).Mag();
			// double scattering_correction = (gOMpositions[event->HitChannel(i)].Z() < matrixPosition.Z()) ? (gScatteringCorrection/15.0)*(matrixPosition.Z() - gOMpositions[event->HitChannel(i)].Z()) : 0;
	        double expectedTime = matrixTime + distanceFromLEDmatrix*ReciprocalSpeedOfLightinWater;// + scattering_correction;

	        h_timeRes->Fill(cascade->noiseTime[i] - expectedTime);
	        if(cascade->noiseTime[i] >= expectedTime-gTCutTimeWindowNs && cascade->noiseTime[i] <= expectedTime + gTCutTimeWindowNs)
	        {
	        	// if (OMIDAlreadyInGPulses(event->GetImpulse(i)))
	        		// continue;
	        	g_pulses.push_back(PulsesVariables{cascade->noiseChID[i],cascade->noiseTime[i],cascade->noiseCharge[i],0,0});
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
	// cout << "High charge based initial estimation" << endl;
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

void EstimateInitialMatrixPositionMC(TVector3& position, double& matrixTime, mcCascade* cascade)
{
	position = TVector3(cascade->position[0],cascade->position[1],cascade->position[2]);
	// position = TVector3(15,136,-100);
	// position = TVector3(-57,-39,-170);
	matrixTime = 0;
}


void EstimateInitialMatrixPositionMC(TVector3& position, double& matrixTime)//, mcCascade* cascade)
{
	// cout << "Grid based initial estimation" << endl;
	int nPar = 0;
	double* gin = new double(0);
	double likelihoodValue = 0;
	double cascadeParameters[7];
	int iflag = 0;

	double lowestLog = 1000000;

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

	for (int i = 0; i < 30; ++i)
	{
		for (int j = 0; j < 30; ++j)
		{
			for (int k = 0; k < 60; ++k)
			{
				// Int_t &npar, Double_t* gin, Double_t &f, Double_t* par, Int_t iflag
				cascadeParameters[0] = -150+10*i+gOMpositions[269].X();
				cascadeParameters[1] = -150+10*j+gOMpositions[269].Y();
				cascadeParameters[2] = -300+10*k+gOMpositions[269].Z();
				double distance = TMath::Sqrt(TMath::Power(cascadeParameters[0]-gOMpositions[g_pulses[maxChargePulseID].OMID].X(),2)+TMath::Power(cascadeParameters[1]-gOMpositions[g_pulses[maxChargePulseID].OMID].Y(),2)+TMath::Power(cascadeParameters[2]-gOMpositions[g_pulses[maxChargePulseID].OMID].Z(),2));
				cascadeParameters[3] = g_pulses[maxChargePulseID].time - distance*ReciprocalSpeedOfLightinWater;
				chi2(nPar,gin,likelihoodValue,cascadeParameters,iflag);
				// cout << likelihoodValue << " " << cascadeParameters[0] << " " << cascadeParameters[1] << " " << cascadeParameters[2] << endl;
				if (likelihoodValue < lowestLog)
				{
					lowestLog = likelihoodValue;
					position[0] = cascadeParameters[0];
					position[1] = cascadeParameters[1];
					position[2] = cascadeParameters[2];
					matrixTime = cascadeParameters[3];
				}
			}
		}
	}
	// cout << lowestLog << " " << position[0] << " " << position[1] << " " << position[2] << " " << matrixTime << endl;
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
	g_QvsL->SetTitle("Charge vs. Distance;Distance from cascade [m]; Charge [p.e.]");
	cEvent->Write();

	delete cEvent;
	delete g_QvsL;
	
	for (int i = 0; i < gNStrings; ++i)
	{
		delete mg_hitsMatrix[i];
	}
	return 0;
}

int EventVisualization(BEvent* event, TVector3& position, double& matrixTime, int eventID)
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

	for(int k = 0; k < event->NHits(); k++)
	{
		nHitsPerString[event->GetImpulseChannelID(k)/36]++;
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
	for(int k = 0; k < event->NHits(); k++)
	{
		int stringID = event->GetImpulseChannelID(k)/36;
		g_hits[stringID]->SetPoint(nHitsPerString[stringID],event->T(k),gOMpositions[event->GetImpulseChannelID(k)].Z());
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
	g_QvsL->SetTitle("Charge vs. Distance;Distance from cascade [m]; Charge [p.e.]");
	cEvent->Write();

	delete cEvent;
	delete g_QvsL;
	
	for (int i = 0; i < gNStrings; ++i)
	{
		delete mg_hitsMatrix[i];
	}
	return 0;
}

int EventVisualization(mcCascade* cascade, TVector3& position, double& matrixTime, int eventID)
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

	for (int k = 0; k < cascade->nHits; k++)
	{
		nHitsPerString[(cascade->chID[k]-1)/36]++;
	}
	for (int i = 0; i < cascade->nNoiseHits; ++i)
	{
		nHitsPerString[(cascade->noiseChID[i])/36]++;
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
	for(int k = 0; k < cascade->nHits; k++)
	{
		int stringID = (cascade->chID[k]-1)/36;
		g_hits[stringID]->SetPoint(nHitsPerString[stringID],cascade->time[cascade->chID[k]-1],gOMpositions[cascade->chID[k]-1].Z());
		nHitsPerString[stringID]++;
	}
	for(int k = 0; k < cascade->nNoiseHits; k++)
	{
		int stringID = (cascade->noiseChID[k])/36;
		g_hits[stringID]->SetPoint(nHitsPerString[stringID],cascade->noiseTime[k],gOMpositions[cascade->noiseChID[k]].Z());
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
	// for (int i = 0; i < cascade->nNoiseHits; ++i)
	// {
	// 	double distanceFromLEDmatrix = (position - gOMpositions[cascade->noiseChID[i]]).Mag();
	// 	g_QvsL->SetPoint(i,distanceFromLEDmatrix,cascade->charge[cascade->chID[i]]-1);
	// }

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
	g_QvsL->SetTitle("Charge vs. Distance;Distance from cascade [m]; Charge [p.e.]");
	cEvent->Write();

	delete cEvent;
	delete g_QvsL;
	
	for (int i = 0; i < gNStrings; ++i)
	{
		delete mg_hitsMatrix[i];
	}
	return 0;
}

int ChargeVisualization(int eventID, double X, double Y, double Z, double energy, double theta, double phi)
{
	double tableParameters[4]{0};
	double par[7]{X,Y,Z,0,energy,theta,phi};

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
		g_MeasQ[stringID]->SetPoint(nHitsPerString[stringID],g_pulses[k].OMID,g_pulses[k].charge);
		nHitsPerString[stringID]++;
	}

	for (int i = 0; i < gNOMs; ++i)
	{
		int stringID = i/36;
		GetParameters(par,i,tableParameters);
		g_ExpQ[stringID]->SetPoint(i%36,i,GetInterpolatedValue(tableParameters)*100000000*energy);		
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

int PrintCascadeJSON(int eventID, double X, double Y, double Z, double time, double theta, double phi)
{
	if (gMCMu || gMCNu || gMCCas)
		return 1;

	TString outputFileName;
	if (App::Output == "")
		outputFileName = BARS::Data::Directory(BARS::Data::JOINT, BARS::App::Season, BARS::App::Cluster, BARS::App::Run, gProductionID.c_str());
	else
		outputFileName = Form("%s/data/exp%d/cluster%d/%04d/",App::Output,BARS::App::Season,BARS::App::Cluster,BARS::App::Run);
	char fname[100];
	std::sprintf(fname,"cascade_season%d_cluster%d_run%d_evt%d.json",BARS::App::Season, BARS::App::Cluster, BARS::App::Run, eventID);
	outputFileName += fname;
	std::ofstream fOutputFile;
	fOutputFile.open(outputFileName);
	fOutputFile<<"{"<<std::endl;
	fOutputFile<<"\t\"eventID\": "<<eventID<<","<<std::endl;
	fOutputFile<<"\t\"season\": "<<BARS::App::Season<<","<<std::endl;
	fOutputFile<<"\t\"cluster\": "<<BARS::App::Cluster<<","<<std::endl;
	fOutputFile<<"\t\"run\": "<<BARS::App::Run<<","<<std::endl;
	fOutputFile<<"\t\"pulses\": {"<<std::endl;

	Int_t nImpulse = g_pulses.size();
	for (int i = 0; i < nImpulse-1; i++)
	{
		fOutputFile<<"\t\t\"" << i << "\": \{ \"mask\": "<< 1 <<
		  ", \"amplitude\": "<< g_pulses[i].charge <<
		  ", \"charge\": "<< g_pulses[i].charge <<
		  ", \"time\": "<< g_pulses[i].time <<
		  ", \"channelID\": "<< g_pulses[i].OMID <<
		  " },"<<std::endl;
	}
	fOutputFile<<"\t\t\"" << nImpulse-1 << "\": \{ \"mask\": "<< 1 <<
	  ", \"amplitude\": "<< g_pulses[nImpulse-1].charge <<
	  ", \"charge\": "<< g_pulses[nImpulse-1].charge <<
	  ", \"time\": "<< g_pulses[nImpulse-1].time <<
	  ", \"channelID\": "<< g_pulses[nImpulse-1].OMID <<
	  " }"<<std::endl;
	fOutputFile<<"\t},"<<std::endl;
	fOutputFile<<"\t\"geometry\": ["<<std::endl;
	for (int i=0; i<gNOMs-1; i++){
	 fOutputFile<<"\t\t{\"channelID\": "<<i<<
	            ", \"x\": "<<gOMpositions[i].X()<<
	            ", \"y\": "<<gOMpositions[i].Y()<<
	            ", \"z\": "<<gOMpositions[i].Z()<<"},"<<std::endl;
	}
	fOutputFile<<"\t\t{\"channelID\": "<<gNOMs-1<<
	            ", \"x\": "<<gOMpositions[gNOMs-1].X()<<
	            ", \"y\": "<<gOMpositions[gNOMs-1].Y()<<
	            ", \"z\": "<<gOMpositions[gNOMs-1].Z()<<"}"<<std::endl;
	fOutputFile<<"\t],"<<std::endl;
	fOutputFile<<"\t\"origins\": {"<<std::endl;
	fOutputFile<<"\t\t\"cascades\": [{"<<std::endl;
	fOutputFile<<"\t\t\t\"mc\": false,"<<std::endl;
	fOutputFile<<"\t\t\t\"title\": \"cascadeFit\","<<std::endl;
	fOutputFile<<"\t\t\t\"direction\": {"<<std::endl;
	fOutputFile<<"\t\t\t\t\"theta\": "<<std::right<<theta<<","<<std::endl;
	fOutputFile<<"\t\t\t\t\"phi\": "<<std::right<<phi<<","<<std::endl;
	fOutputFile<<"\t\t\t\t\"x\": "<<std::right<<X<<","<<std::endl;
	fOutputFile<<"\t\t\t\t\"y\": "<<std::right<<Y<<","<<std::endl;
	fOutputFile<<"\t\t\t\t\"z\": "<<std::right<<Z<<""<<std::endl;
	// fOutputFile<<"\t\t\t\t\"time\": "<<std::right<<time<<std::endl;
	fOutputFile<<"\t\t\t}"<<std::endl;
	fOutputFile<<"\t\t}]"<<std::endl;
	fOutputFile<<"\t}"<<std::endl;
	fOutputFile<<"}"<<std::endl;
	fOutputFile.close();

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
		allCharge += impulseTel->GetQ(i)>0?impulseTel->GetQ(i)/gOMQCal[impulseTel->GetNch(i)]:0;
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

bool QRatioFilterPassed(BEvent* event, double &qRatio)
{
	double allCharge = 0;
	for (int i = 0; i < event->NHits(); ++i)
	{
		allCharge += event->Q(i)>0?event->Q(i):0;
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

bool QRatioFilterPassed(mcCascade* cascade, double &qRatio)
{
	double allCharge = 0;
	for (int i = 0; i < cascade->nHits; ++i)
	{
		allCharge += cascade->charge[cascade->chID[i]-1];
	}
	for (int i = 0; i < cascade->nNoiseHits; ++i)
	{
		allCharge += cascade->noiseCharge[i];
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
	// cout << "In likelihood" << endl;
	double lowestLog = 10000;
	double cascadeEnergy = 10;
	for (int k = 0; k < 8; ++k)
	{
		for (int l = 0; l < 12; ++l)
		{
			// cout << k << " " << l << endl;
			double cascadeTheta = TMath::Pi()/8*k;
			double cascadePhi = TMath::Pi()*2/12*l;
			double recentLog = FitMatrixDirection(matrixPosition,matrixTime,cascadeEnergy,cascadeTheta,cascadePhi);
			// double recentLog = 100000;
			if (recentLog < lowestLog)
			{
				lowestLog = recentLog;
				energy = cascadeEnergy;
				theta = cascadeTheta;
				phi = cascadePhi;
			}
			// cout << k << " " << l << " " << recentLog  << " " << cascadeEnergy << " " << cascadeTheta << " " << cascadePhi << endl;
		}
	}
	likelihood = lowestLog;
	// cout << "End likelihood" << endl;
	// cout << lowestLog << " " << energy << " " << theta << " " << phi << endl;
	if (lowestLog < gLikelihoodCut)
		return true;
		// cout << i << " " << energy << " " << theta << " " << phi << endl;
	return false;
}

bool LikelihoodFilterPassedMC(TVector3 &matrixPosition, double &matrixTime, double &energy, double &theta, double &phi, double &likelihood)
{
	double lowestLog = 10000;
	double cascadeEnergy = 0;

	lowestLog = FitMatrixDirection(matrixPosition,matrixTime,energy,theta,phi);

	// for (int k = 0; k < 12; ++k)
	// {
	// 	for (int l = 0; l < 8; ++l)
	// 	{
	// 		double cascadeTheta = TMath::Pi()/12*k;
	// 		double cascadePhi = TMath::Pi()*2/8*l;
	// 		double recentLog = FitMatrixDirection(matrixPosition,matrixTime,cascadeEnergy,cascadeTheta,cascadePhi);
	// 		// double recentLog = 100000;
	// 		if (recentLog < lowestLog)
	// 		{
	// 			lowestLog = recentLog;
	// 			energy = cascadeEnergy;
	// 			theta = cascadeTheta;
	// 			phi = cascadePhi;
	// 		}
	// 		// cout << k << " " << l << " " << recentLog  << " " << cascadeEnergy << " " << cascadeTheta << " " << cascadePhi << endl;
	// 	}
	// }
	likelihood = lowestLog;
	// cout << lowestLog << " " << energy << " " << theta << " " << phi << endl;
	if (lowestLog < gLikelihoodCut)
		return true;
		// cout << i << " " << energy << " " << theta << " " << phi << endl;
	return false;
}

// bool LikelihoodFilter2Passed(TVector3 &matrixPosition, double &matrixTime, double &energy, double &theta, double &phi, double &likelihood)
// {
// 	int nThetaSteps = 50;
// 	int nPhiSteps = 50;
// 	const int nEnergySteps = 50;

// 	double EnergyValues[nEnergySteps]{1,1.5,2,2.5,3,3.5,4,4.5,5,7.5,10,15,25,50,75,100,125,150,175,200,225,250,275,300,325,350,375,400,425,450,475,500,550,600,650,700,750,800,850,900,950,1000,1.1,1.2,1.3,1.4,1.5,1.6,1.8,2};

// 	int nPar = 0;
// 	double* gin = new double(0);
// 	double likelihoodValue = 0;
// 	double cascadeParameters[7];
// 	int iflag = 0;

// 	cascadeParameters[0] = matrixPosition.X();
// 	cascadeParameters[1] = matrixPosition.Y();
// 	cascadeParameters[2] = matrixPosition.Z();
// 	cascadeParameters[3] = matrixTime;

// 	double lowestLog = 1000000;

// 	for (int i = 0; i < nThetaSteps; ++i)
// 	{
// 		for (int j = 0; j < nPhiSteps; ++j)
// 		{
// 			for (int k = 0; k < nEnergySteps; ++k)
// 			{
// 				// Int_t &npar, Double_t* gin, Double_t &f, Double_t* par, Int_t iflag
// 				cascadeParameters[4] = EnergyValues[i];
// 				cascadeParameters[5] = j*TMath::Pi()/nThetaSteps;
// 				cascadeParameters[6] = k*2*TMath::Pi()/nPhiSteps;
// 				logLikelihood(nPar,gin,likelihoodValue,cascadeParameters,iflag);
// 				if (likelihoodValue < lowestLog)
// 				{
// 					lowestLog = likelihoodValue;
// 					energy = cascadeParameters[4];
// 					theta = cascadeParameters[5];
// 					phi = cascadeParameters[6];
// 				}
// 			}
// 		}
// 	}
// 	cout << lowestLog << " " << energy << " " << theta << " " << phi << endl;
// }

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
		TString outputFileName = BARS::Data::Directory(BARS::Data::JOINT, BARS::App::Season, BARS::App::Cluster, BARS::App::Run, gProductionID.c_str());
		outputFileName += Form("recCasc_vis_%d.root",gVisEventID);
		TFile* outputFile = new TFile(outputFileName,"RECREATE");
		tree->GetEntry(gVisEventID);
		TVector3 matrixPosition(0,0,0); 
		double matrixTime = 0;
		EventVisualization(impulseTel,matrixPosition,matrixTime,gVisEventID);
		cout << "Visualization of event: " << gVisEventID << " has been produced!" << endl;
		return 0;
	}

	TString outputFileName2 = BARS::Data::Directory(BARS::Data::JOINT, BARS::App::Season, BARS::App::Cluster, BARS::App::Run, gProductionID.c_str());
	outputFileName2 += "recCasc_nTuple.root";
	TFile* outputFile2;
	if (App::Output == "")
		outputFile2 = new TFile(outputFileName2,"RECREATE");
	else
		outputFile2 = new TFile(Form("%s/exp%d/cluster%d/%04d/recCasc_nTuple.root",App::Output,BARS::App::Season,BARS::App::Cluster,BARS::App::Run),"RECREATE");

   	TNtuple* nt_cascades = new TNtuple("nt_cascades","NTuple of reconstructed cascades","runID:EventID:NPulses:NPulsesT:QTotal:QRatio:CloseHits:Likelihood:X:Y:Z:Time:Energy:Theta:Phi");

	TString outputFileName = BARS::Data::Directory(BARS::Data::JOINT, BARS::App::Season, BARS::App::Cluster, BARS::App::Run, gProductionID.c_str());
	outputFileName += "recCasc_vis.root";
	TFile* outputFile;
	if (App::Output == "")
		outputFile = new TFile(outputFileName,"RECREATE");
	else
		outputFile = new TFile(Form("%s/exp%d/cluster%d/%04d/recCasc_vis.root",App::Output,BARS::App::Season,BARS::App::Cluster,BARS::App::Run),"RECREATE");

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
	double overallCharge = 0;

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
		if (!QFilterPassed(impulseTel,overallCharge))
			continue;
		nQFilter++;

		// cout << i << endl;
		//initialization of R and T for LED matrix
		TVector3 matrixPosition(0,0,0); 
		double matrixTime = 0;
		EstimateInitialMatrixPositionMC(matrixPosition,matrixTime);
		// matrixPosition.Print();
		// fMinuit->SetFCN(chi2);
		fMinuit->SetFCN(MEstimator);
		double chi2QResult = FitMatrixPosition(matrixPosition,matrixTime);
		// matrixPosition.Print();
		h_chi2AfterQ->Fill(chi2QResult);
		if (chi2QResult > gQCutChi2)
			continue;
		nQFilterChi2++;
		if (!TFilterPassed(impulseTel,matrixPosition,matrixTime,nPulsesT))
			continue;
		nTFilter++;
		fMinuit->SetFCN(chi2);
		double chi2TResult = FitMatrixPosition(matrixPosition,matrixTime);
		// matrixPosition.Print();
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
		cout << i << endl;
		if (!LikelihoodFilterPassed(matrixPosition,matrixTime,cascadeEnergy,cascadeTheta,cascadePhi,likelihood))
		{
			nt_cascades->Fill((double)BARS::App::Run,(double)i,(double)nPulses,(double)nPulsesT,overallCharge,qRatio,(double)closeHits,-1,matrixPosition.X(),matrixPosition.Y(),matrixPosition.Z(),matrixTime,0,0,0);
			continue;
		}
		nLikelihoodFilter++;
		nt_cascades->Fill((double)BARS::App::Run,(double)i,(double)nPulses,(double)nPulsesT,overallCharge,qRatio,(double)closeHits,likelihood,matrixPosition.X(),matrixPosition.Y(),matrixPosition.Z(),matrixTime,cascadeEnergy,cascadeTheta,cascadePhi);
		ChargeVisualization(i,matrixPosition.X(),matrixPosition.Y(),matrixPosition.Z(),cascadeEnergy,cascadeTheta,cascadePhi);
		PrintCascadeJSON(i,matrixPosition.X(),matrixPosition.Y(),matrixPosition.Z(),matrixTime,cascadeTheta,cascadePhi);
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

int ReadGeometry(TTree* tree, BExtractedHeader* header) // read dynamic geometry
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

	int entryID = 0;
	long double startTime = 0;

	do
	{
		tree->GetEntry(entryID);
		startTime = header->GetTime().GetSec();
		entryID++;
	}while (startTime == 0);

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
	const char* filePath = BARS::Data::File(BARS::Data::JOINT, BARS::App::Season, BARS::App::Cluster, BARS::App::Run, gProductionID.c_str());
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
	TString outputFileName = "";

	if (!gMCMu && !gMCNu && !gMCCas)
		outputFileName = BARS::Data::Directory(BARS::Data::JOINT, BARS::App::Season, BARS::App::Cluster, BARS::App::Run, gProductionID.c_str());

	if (gMCMu)
		outputFileName = "/Data/BaikalData/mc/2018may/";

	if (gMCNu)
		outputFileName = "/Data/BaikalData/mc/nuatm_feb19/";

	if (gMCCas)
		outputFileName = "/Data/BaikalData/mc/DZH_cascades/";


	outputFileName += "recCasc_hist.root";
	TFile* outputFile;
	if (App::Output == "")
		outputFile = new TFile(outputFileName,"RECREATE");
	else
		outputFile = new TFile(Form("%s/data/exp%d/cluster%d/%04d/recCasc_hist.root",App::Output,BARS::App::Season,BARS::App::Cluster,BARS::App::Run),"RECREATE");

	h_nHits->Write();
	h_chargePerHit->Write();
	h_chargePerEvent->Write();
	h_nHitsAfterQ->Write();
	h_chi2AfterQ->Write();
	h_nHitsAfterT->Write();
	h_chi2AfterT->Write();
	g_cascadePosition->Write();
	h_timeRes->Write();
	h_XYposition->Write();
	delete outputFile;
}

// gLogTable4D structure: Rho, Z, Phi, cosTheta
int ReadLogTable()
{
	cout << "4D LogTable reading starts" << endl;
	ifstream fTab;
	if (App::Output == "")
		fTab.open("/Data/BaikalData/showerTable/hq001200_double.dqho2011", ios::in|ios::binary|ios::ate);
	else
		fTab.open("./hq001200_double.dqho2011", ios::in|ios::binary|ios::ate);

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

void PrintRunInfo(TTree* tree, BExtractedHeader* header)
{
	tree->GetEntry(0);
	long startTime = header->GetTime().GetSec();
	tree->GetEntry(tree->GetEntries()-1);
	long endTime = header->GetTime().GetSec();

	cout << "RunInfo = Number of entries & runTime [hours] & runTime [days]" << endl;
	cout << "! " << tree->GetEntries() << " " << (endTime-startTime)/3600.0 << " " << (endTime-startTime)/3600.0/24.0 << endl; 
	std::cout << std::string(80,'*') << std::endl;
}

int processExperimentalData()
{
	if (!CheckParams())
    	return -1;

    if (gProductionID == "")
    	gProductionID = "barsv051";

    // const char* filePath = BARS::Data::File(BARS::App::Cluster, BARS::App::Season, BARS::App::Run, BARS::Data::JOINT,"r01_i01_j01_t01");
    const char* filePath = BARS::Data::File(BARS::Data::JOINT, BARS::App::Season, BARS::App::Cluster, BARS::App::Run, gProductionID.c_str());

    if (ReadLogTable() == -1)
    {
    	std::cout << "Problem with 4D LogLikelihood file!" << std::endl;
    	return -1;
    }

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
	PrintRunInfo(tree,header);
	cout << "Season: " << BARS::App::Season << " Cluster: " << BARS::App::Cluster << " Run: " <<  BARS::App::Run << endl;
	if (ReadGeometry(tree,header) == -1)
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

int CountHitsFromCascades(BEventMaskMC* eventMask)
{
	int nHits = 0;
	// cout << i << " " << eventMask->GetEntries() << " " << event->NHits() << endl;
	for (int k = 0; k < eventMask->GetEntries(); ++k)
	{
		if (eventMask->GetFlag(k) > 0)
			nHits++;
		// cout << "\t" << eventMask->GetFlag(k) << endl;
	}
	return nHits;
}

int CountHitsFromMuons(BEventMaskMC* eventMask)
{
	int nHits = 0;
	// cout << i << " " << eventMask->GetEntries() << " " << event->NHits() << endl;
	for (int k = 0; k < eventMask->GetEntries(); ++k)
	{
		if (eventMask->GetFlag(k) < 0)
			nHits++;
		// cout << "\t" << eventMask->GetFlag(k) << endl;
	}
	return nHits;
}

double CountChargeFromCascades(BEvent* event, BEventMaskMC* eventMask)
{
	double charge = 0;
	// cout << i << " " << eventMask->GetEntries() << " " << event->NHits() << endl;
	for (int k = 0; k < eventMask->GetEntries(); ++k)
	{
		if (eventMask->GetFlag(k) > 0)
			charge += event->Q(k);
		// cout << "\t" << eventMask->GetFlag(k) << endl;
	}
	return charge;
}

double CountChargeFromMuons(BEvent* event, BEventMaskMC* eventMask)
{
	double charge = 0;
	// cout << i << " " << eventMask->GetEntries() << " " << event->NHits() << endl;
	for (int k = 0; k < eventMask->GetEntries(); ++k)
	{
		if (eventMask->GetFlag(k) < 0)
			charge += event->Q(k);
		// cout << "\t" << eventMask->GetFlag(k) << endl;
	}
	return charge;
}

int TrueCascadeParameters(BSourceEAS* mcEvent, TVector3 &cascadePositionTrue, double &cascadeEnergyTrue, double &cascadeThetaTrue, double &cascadePhiTrue, int &flagID)
{
	cascadeThetaTrue = mcEvent->GetMuonTrack(0)->GetPolarAngle();
	cascadePhiTrue = mcEvent->GetMuonTrack(0)->GetAzimuthAngle();

	double highestEnergy = 0;
	int muonID = -1;
	int cascadeID = -1;

	for (int k = 0; k < mcEvent->GetNumTracks(); ++k)
	{
		for (int j = 0; j < mcEvent->GetMuonTrack(k)->GetNumShowers(); ++j)
		{
			if (mcEvent->GetMuonTrack(k)->GetShower(j)->GetE() > highestEnergy)
			{
				highestEnergy = mcEvent->GetMuonTrack(k)->GetShower(j)->GetE();
				muonID = k;
				cascadeID = j;
				flagID = (cascadeID+1)*1000+(muonID+1);
			}
			// cout << k << " " << j << " " << mcEvent->GetMuonTrack(k)->GetShower(j)->GetX() << " " << mcEvent->GetMuonTrack(k)->GetShower(j)->GetY() << " " << mcEvent->GetMuonTrack(k)->GetShower(j)->GetZ() << " " << mcEvent->GetMuonTrack(k)->GetShower(j)->GetE() << endl;
		}
	}		
	cascadeEnergyTrue = mcEvent->GetMuonTrack(muonID)->GetShower(cascadeID)->GetE();
	cascadePositionTrue.SetX(mcEvent->GetMuonTrack(muonID)->GetShower(cascadeID)->GetX());
	cascadePositionTrue.SetY(mcEvent->GetMuonTrack(muonID)->GetShower(cascadeID)->GetY());
	cascadePositionTrue.SetZ(mcEvent->GetMuonTrack(muonID)->GetShower(cascadeID)->GetZ());

	return 0;
}

int EvalPurityEfficiency(BEventMaskMC* eventMask, int flagID, double &efficiency, double &eventPurity, double &cascadePurity)
{
	int nHitsFromCascade = 0;
	for (int i = 0; i < eventMask->GetEntries(); ++i)
	{
		if (eventMask->GetFlag(i) == flagID)
			nHitsFromCascade++;
	}

	int nFoundHitsFromCascade = 0;
	int nNoiseHits = 0;
	int nHitsOtherPhysics = 0;
	for (int i = 0; i < g_pulses.size(); ++i)
	{
		if (g_pulses[i].MCflag == flagID)
			nFoundHitsFromCascade++;
		if (g_pulses[i].MCflag == 0)
			nNoiseHits++;
		if (g_pulses[i].MCflag != flagID && g_pulses[i].MCflag != 0)
			nHitsOtherPhysics++;
	}

	efficiency = nFoundHitsFromCascade/(double)nHitsFromCascade;
	eventPurity = (g_pulses.size()-nNoiseHits)/(double)g_pulses.size();
	cascadePurity = (g_pulses.size()-nNoiseHits-nHitsOtherPhysics)/(double)g_pulses.size();
	// cout << "Eff: " << efficiency << " eventPur: " << eventPurity << " cascadePur: " << cascadePurity << endl; 

	return 0;	
}

int DoTheMagicMC(TChain* tree, BEvent* event, BEventMaskMC* eventMask, BSourceEAS* mcEvent)
{

	TString outputFileName2 = "";
	TString outputFileName = "";
	if (gMCMu)
	{
		outputFileName2 = "/Data/BaikalData/mc/2018may/";
		outputFileName = "/Data/BaikalData/mc/2018may/";
	}
	if (gMCNu)
	{
		outputFileName2 = "/Data/BaikalData/mc/nuatm_feb19/";
		outputFileName = "/Data/BaikalData/mc/nuatm_feb19/";
	}

	outputFileName2 += "recCasc_nTuple.root";
	TFile* outputFile2 = new TFile(outputFileName2,"RECREATE");

   	TNtuple* nt_cascades = new TNtuple("nt_cascades","NTuple of reconstructed cascades",TString("runID:EventID:NPulses:NPulsesT:QRatio:CloseHits:Likelihood:X:Y:Z:Time:InitDist:FinalDist:Energy:Theta:Phi:TrueX:TrueY:TrueZ:TrueEnergy:TrueTheta:TruePhi:Efficiency:EventPurity:CascadePurity"));

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
	TVector3 initialPosEstim;

	for (int i = 0; i < nEntries; ++i)
	{
		if (i%(nEntries/10) == 0)
		{
			cout << round((double)(i)/nEntries*100) << "% ";
			cout << std::flush;
		}
		tree->GetEntry(i);
		if (!NFilterPassed(event,nPulses))
			continue;
		nNFilter++;
		if (!QFilterPassed(event))
			continue;
		nQFilter++;

		//initialization of R and T for LED matrix
		TVector3 matrixPosition(0,0,0); 
		double matrixTime = 0;
		EstimateInitialMatrixPosition(matrixPosition,matrixTime);
		initialPosEstim = matrixPosition;
		fMinuit->SetFCN(chi2);
		double chi2QResult = FitMatrixPosition(matrixPosition,matrixTime);
		h_chi2AfterQ->Fill(chi2QResult);
		if (chi2QResult > gQCutChi2)
			continue;
		nQFilterChi2++;
		if (!TFilterPassed(event,eventMask,matrixPosition,matrixTime,nPulsesT))
			continue;
		nTFilter++;
		fMinuit->SetFCN(chi2);
		double chi2TResult = FitMatrixPosition(matrixPosition,matrixTime);
		h_chi2AfterT->Fill(chi2TResult);
		if (chi2TResult > gTCutChi2)
			continue;
		nTFilterChi2++;
		TFilterPassed(event,eventMask,matrixPosition,matrixTime,nPulsesT);
		if (!ZFilterPassed(matrixPosition))
			continue;
		nZFilter++;
		FillCascPos(matrixPosition);
		if (!TDelayFilterPassed(matrixTime))
			continue;
		nTDelayFilter++;
		if (!QRatioFilterPassed(event,qRatio))
			continue;
		nQRatioFilter++;
		if (!BranchFilterPassed(matrixPosition))
			continue;
		nBranchFilter++;
		if (!CloseHitsFilterPassed(matrixPosition,closeHits))
			continue;
		nCloseHitsFilter++;
		EventVisualization(event,matrixPosition,matrixTime,i);

		double cascadeEnergyTrue = 0;
		double cascadeThetaTrue = 0;
		double cascadePhiTrue = 0;
		TVector3 cascadePositionTrue(0,0,0);
		int flagID = 0;
		double efficiency = 0;
		double eventPurity = 0;
		double cascadePurity = 0;
		TrueCascadeParameters(mcEvent,cascadePositionTrue,cascadeEnergyTrue,cascadeThetaTrue,cascadePhiTrue,flagID);
		EvalPurityEfficiency(eventMask,flagID,efficiency,eventPurity,cascadePurity);

		double cascadeEnergy = cascadeEnergyTrue;
		double cascadeTheta = cascadeThetaTrue;
		double cascadePhi = cascadePhiTrue;
		double likelihood = 0;
		double initDist = TMath::Sqrt(TMath::Power(cascadePositionTrue.X()-initialPosEstim.X(),2)+TMath::Power(cascadePositionTrue.Y()-initialPosEstim.Y(),2)+TMath::Power(cascadePositionTrue.Z()-initialPosEstim.Z(),2));
		double finalDist = TMath::Sqrt(TMath::Power(cascadePositionTrue.X()-matrixPosition.X(),2)+TMath::Power(cascadePositionTrue.Y()-matrixPosition.Y(),2)+TMath::Power(cascadePositionTrue.Z()-matrixPosition.Z(),2));
		// cout << i << endl;
		cout << i << " " << cascadeEnergyTrue << " " << cascadeThetaTrue << " " << cascadePhiTrue << endl; 
		// LikelihoodFilter2Passed(matrixPosition,matrixTime,cascadeEnergy,cascadeTheta,cascadePhi,likelihood);
		if (!LikelihoodFilterPassed(matrixPosition,matrixTime,cascadeEnergy,cascadeTheta,cascadePhi,likelihood))
		{
			float parameters[25]{(float)BARS::App::Run,(float)i,(float)nPulses,(float)nPulsesT,qRatio,(float)closeHits,-1.0,matrixPosition.X(),matrixPosition.Y(),matrixPosition.Z(),matrixTime,initDist,finalDist,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0}; 
			nt_cascades->Fill(parameters);
			continue;
		}
		cout << i << " " << cascadeEnergy << " " << cascadeTheta << " " << cascadePhi << endl; 	
		nLikelihoodFilter++;
		ChargeVisualization(i,matrixPosition.X(),matrixPosition.Y(),matrixPosition.Z(),cascadeEnergy,cascadeTheta,cascadePhi);

		float parameters[25]{(float)BARS::App::Run,(double)i,(double)nPulses,(double)nPulsesT,qRatio,(double)closeHits,likelihood,matrixPosition.X(),matrixPosition.Y(),matrixPosition.Z(),matrixTime,initDist,finalDist,cascadeEnergy,cascadeTheta,cascadePhi,cascadePositionTrue.X(),cascadePositionTrue.Y(),cascadePositionTrue.Z(),cascadeEnergyTrue,cascadeThetaTrue,cascadePhiTrue,efficiency,eventPurity,cascadePurity};
		nt_cascades->Fill(parameters);
		// cout << string(80,'*') << endl;
		// cout << i << " " << matrixPosition.X() << " " << matrixPosition.Y() << " " << matrixPosition.Z() << " " << cascadeEnergy << " " << cascadeTheta << " " << cascadePhi << endl;
		// cout << mcEvent->GetMuonTrack(0)->GetPolarAngle() << " " << mcEvent->GetMuonTrack(0)->GetAzimuthAngle() << endl;
		// cout << "NumberTracks: " << mcEvent->GetNumTracks() << endl;
		// for (int k = 0; k < mcEvent->GetNumTracks(); ++k)
		// {
		// 	cout << "	number of cascades: " << mcEvent->GetMuonTrack(k)->GetNumShowers() << endl;
		// 	for (int j = 0; j < mcEvent->GetMuonTrack(k)->GetNumShowers(); ++j)
		// 	{
		// 		cout << k << " " << j << " " << mcEvent->GetMuonTrack(k)->GetShower(j)->GetX() << " " << mcEvent->GetMuonTrack(k)->GetShower(j)->GetY() << " " << mcEvent->GetMuonTrack(k)->GetShower(j)->GetZ() << " " << mcEvent->GetMuonTrack(k)->GetShower(j)->GetE() << endl;
		// 	}
		// }
		
		// PrintGPulses();
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

	nt_cascades->Write();

	outputFile2->Close();

	return 0;
}

int ReadGeometryMC(TFile* file)
{
	TTree* tree = (TTree*)file->Get("ArrayConfig");

	BGeomTelMC* geomMC = NULL;
	tree->SetBranchAddress("BGeomTel.",&geomMC);

	int nOKOMs = 0;

	tree->GetEntry(0);
	for (int i = 0; i < geomMC->GetNumOMs(); ++i)
	{
		gOMpositions[i] = geomMC->At(i)->GetPosition();
		nOKOMs++;
		if ((i > 35 && i < 60) || (i > 71 && i < 84) || i == 130 || i == 131 || i == 245 || i == 246 || i == 247 || i == 256)
			gOMQCal[i] = -1;
	}
	return nOKOMs;

}

int processMCData()
{
	TChain* mcFiles = new TChain("Events");
	const char* filePath = " ";

	if (gMCNu && gMCMu)
	{
		cout << "ERROR: You can't process both MC neutrinos and MC muons at once!" << endl;
		return -1;
	}

	if (gMCMu)
	{
		mcFiles->Add("/Data/BaikalData/mc/2018may/n_cors_n2m_cl2016_x*.root");
		filePath = "/Data/BaikalData/mc/2018may/n_cors_n2m_cl2016_x1001.root";
	}
	if (gMCNu)
	{
		mcFiles->Add("/Data/BaikalData/mc/nuatm_feb19/n_nuatm_gs_n2m_cl2016_x*.root");	
		filePath = "/Data/BaikalData/mc/nuatm_feb19/n_nuatm_gs_n2m_cl2016_x1001.root";
	}

	
	TFile* file = new TFile(filePath,"READ");
    // TTree* tree = (TTree*)file->Get("Events");

    BEvent* event = NULL;
    mcFiles->SetBranchAddress("BEvent.",&event);
    BEventMaskMC* eventMask = NULL;
    mcFiles->SetBranchAddress("MCEventMask.",&eventMask);
    BSourceEAS* mcEvent = NULL;
    mcFiles->SetBranchAddress("MCEventSource.",&mcEvent);

    if (ReadLogTable() == -1)
    {
    	std::cout << "Problem with 4D LogLikelihood file!" << std::endl;
    	return -1;
    }

    PrintConfig();
    cout << "RunInfo = Number of entries & runTime [hours] & runTime [days]" << endl;
    cout << "! " << mcFiles->GetEntries() << " " << (mcFiles->GetListOfFiles()->GetEntries()*3900)/3600.0 << "  " << (mcFiles->GetListOfFiles()->GetEntries()*3900)/3600.0/24.0 << endl; 
    std::cout << std::string(80,'*') << std::endl;

    if (ReadGeometryMC(file) == -1)
	{
		std::cout << "Problem with geometry file!" << std::endl;
		return -1;
	}

    DoTheMagicMC(mcFiles,event,eventMask,mcEvent);
	SaveHistograms();

    return 0;
}

int GenerateNoise(mcCascade* cascade)
{
	ifstream inputFile("/Data/BaikalData/showerTable/mc-noise-charge-2016.dat", ios::in);

   	if (!inputFile)
    {
    	cerr << "Calibration file: " << "/Data/BaikalData/showerTable/mc-noise-charge-2016.dat" << " was NOT found. Program termination!" << endl;
    	return -1;
  	}

    string dummyLine;
    double readValue = 0;

    getline(inputFile,dummyLine);
    vector<double> noiseTable;

    for (int i = 0; i < 500; ++i)
    {
		inputFile >> readValue;
		noiseTable.push_back(readValue);
		// cout << readValue << endl;
    }
    inputFile.close();

	TRandom2 ranGen;
	cascade->nNoiseHits = 0;
	for (int i = 0; i < gNOMs; ++i)
	{
		int noisePulses = ranGen.Poisson(0.25);
		// cout << "OMID: " << i << " nNoisePulses: " << noisePulses << endl;
		for (int j = 0; j < noisePulses; ++j)
		{
			double noiseCharge = 0;//ranGen.Landau(1,2);
			double ranNumber = ranGen.Uniform(1);
			for (int j = 0; j < noiseTable.size(); ++j)
			{
				if (noiseTable[j] > ranNumber)
				{
					noiseCharge = j*0.05;
					break;
				}
			}
			double noiseTime = ranGen.Uniform(5120)-2560;
			cascade->noiseCharge[cascade->nNoiseHits] = noiseCharge;
			cascade->noiseTime[cascade->nNoiseHits] = noiseTime;
			cascade->noiseChID[cascade->nNoiseHits] = i;
			cascade->nNoiseHits++;
			// cout << "\t" << j << " " << noiseTime << " " << noiseCharge << endl; 
		}
	}
	return 0;
}

int DoTheMagicMCCascades(TChain* tree, mcCascade* cascade)
{
	TString outputFileName2 = "/Data/BaikalData/mc/DZH_cascades/";
	TString outputFileName = "/Data/BaikalData/mc/DZH_cascades/";

	outputFileName2 += "recCasc_nTuple.root";
	TFile* outputFile2 = new TFile(outputFileName2,"RECREATE");

   	TNtuple* nt_cascades = new TNtuple("nt_cascades","NTuple of reconstructed cascades",TString("runID:EventID:NPulses:NPulsesT:QRatio:CloseHits:Likelihood:X:Y:Z:Time:InitDist:FinalDist:Energy:Theta:Phi:TrueX:TrueY:TrueZ:TrueEnergy:TrueTheta:TruePhi:Efficiency:EventPurity:CascadePurity"));

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
	int nEvents = 0;
	TVector3 initialPosEstim;

	for (int i = 0; i < nEntries; ++i)
	{
		if (i%(nEntries/10) == 0)
		{
			cout << round((double)(i)/nEntries*100) << "% ";
			cout << std::flush;
		}
		tree->GetEntry(i);
		if (cascade->showerEnergy > 1000 || TMath::Sqrt(TMath::Power(cascade->position[0],2)+TMath::Power(cascade->position[1],2)) > 60 || i % 1000 != 0)
			continue;
		nEvents++;
		cascade->nNoiseHits = 0;
		GenerateNoise(cascade);
		// cout << i << endl;
		if (!NFilterPassed(cascade,nPulses))
			continue;
		nNFilter++;
		if (!QFilterPassed(cascade))
			continue;
		nQFilter++;
		// PrintGPulses();
		// initialization of R and T for LED matrix
		TVector3 matrixPosition(0,0,0); 
		double matrixTime = 0;
		// EstimateInitialMatrixPositionMC(matrixPosition,matrixTime,cascade);
		EstimateInitialMatrixPositionMC(matrixPosition,matrixTime);
		// cout << "TruePosition: " << cascade->position[0] << " " << cascade->position[1] << " " << cascade->position[2] << endl;
		initialPosEstim = matrixPosition;
		// matrixPosition.Print();
		// cout << matrixTime << endl;
		// fMinuit->SetFCN(chi2);
		fMinuit->SetFCN(MEstimator);
		double chi2QResult = FitMatrixPosition(matrixPosition,matrixTime);
		h_chi2AfterQ->Fill(chi2QResult);
		// matrixPosition.Print();
		// cout << matrixTime << endl;
		// cout << chi2QResult << endl;
		if (chi2QResult > gQCutChi2)
			continue;
		nQFilterChi2++;
		if (!TFilterPassed(cascade,matrixPosition,matrixTime,nPulsesT))
			continue;
		// PrintGPulses();
		nTFilter++;
		fMinuit->SetFCN(chi2);
		double chi2TResult = FitMatrixPosition(matrixPosition,matrixTime);
		h_chi2AfterT->Fill(chi2TResult);
		// matrixPosition.Print();
		if (chi2TResult > gTCutChi2)
			continue;
		nTFilterChi2++;
		TFilterPassed(cascade,matrixPosition,matrixTime,nPulsesT);
		if (!ZFilterPassed(matrixPosition))
			continue;
		nZFilter++;
		FillCascPos(matrixPosition);
		if (!TDelayFilterPassed(matrixTime))
			continue;
		nTDelayFilter++;
		if (!QRatioFilterPassed(cascade,qRatio))
			continue;
		nQRatioFilter++;
		if (!BranchFilterPassed(matrixPosition))
			continue;
		nBranchFilter++;
		if (!CloseHitsFilterPassed(matrixPosition,closeHits))
			continue;
		nCloseHitsFilter++;
		EventVisualization(cascade,matrixPosition,matrixTime,i);

		double cascadeEnergyTrue = cascade->showerEnergy;
		double cascadeThetaTrue = TMath::Pi()-TMath::ACos(cascade->cosTheta);
		double cascadePhiTrue = cascade->phi+TMath::Pi();
		TVector3 cascadePositionTrue(cascade->position[0],cascade->position[1],cascade->position[2]);
		double efficiency = 1;
		double eventPurity = 1;
		double cascadePurity = 1;

		double cascadeEnergy = cascadeEnergyTrue;
		double cascadeTheta = cascadeThetaTrue;
		double cascadePhi = cascadePhiTrue;
		double likelihood = 0;
		double initDist = TMath::Sqrt(TMath::Power(cascade->position[0]-initialPosEstim.X(),2)+TMath::Power(cascade->position[1]-initialPosEstim.Y(),2)+TMath::Power(cascade->position[2]-initialPosEstim.Z(),2));
		double finalDist = TMath::Sqrt(TMath::Power(cascade->position[0]-matrixPosition.X(),2)+TMath::Power(cascade->position[1]-matrixPosition.Y(),2)+TMath::Power(cascade->position[2]-matrixPosition.Z(),2));
		// cout << i << endl;
		// cout << i << " " << cascadeEnergyTrue << " " << cascadeThetaTrue << " " << cascadePhiTrue << endl; 
		// LikelihoodFilter2Passed(matrixPosition,matrixTime,cascadeEnergy,cascadeTheta,cascadePhi,likelihood);
		if (!LikelihoodFilterPassed(matrixPosition,matrixTime,cascadeEnergy,cascadeTheta,cascadePhi,likelihood))
		{
			float parameters[25]{(float)BARS::App::Run,(float)i,(float)nPulses,(float)nPulsesT,qRatio,(float)closeHits,-1.0,matrixPosition.X(),matrixPosition.Y(),matrixPosition.Z(),matrixTime,initDist,finalDist,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0}; 
			nt_cascades->Fill(parameters);
			// cout << "CONTINUE" << endl;
			continue;
		}
		// cout << i << " " << cascadeEnergy << " " << cascadeTheta << " " << cascadePhi << endl; 	
		nLikelihoodFilter++;
		ChargeVisualization(i,matrixPosition.X(),matrixPosition.Y(),matrixPosition.Z(),cascadeEnergy,cascadeTheta,cascadePhi);

		float parameters[25]{(float)BARS::App::Run,(double)i,(double)nPulses,(double)nPulsesT,qRatio,(double)closeHits,likelihood,matrixPosition.X(),matrixPosition.Y(),matrixPosition.Z(),matrixTime,initDist,finalDist,cascadeEnergy,cascadeTheta,cascadePhi,cascadePositionTrue.X(),cascadePositionTrue.Y(),cascadePositionTrue.Z(),cascadeEnergyTrue,cascadeThetaTrue,cascadePhiTrue,efficiency,eventPurity,cascadePurity};
		nt_cascades->Fill(parameters);
	}

	cout << "\nNentries: \t" << nEvents << endl;
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

	nt_cascades->Write();

	outputFile2->Close();

	return 0;
}



int ReadGeometryMCCascades()
{
	TString fileName = "/Data/BaikalData/mc/DZH_cascades/array2016_phys_kebkal.dat";

	ifstream inputFile;
    inputFile.open(fileName);

    if (!inputFile)
    {
    	cerr << "Geometry file for MC cascades: " << fileName << " was NOT found. Program termination!" << endl;
    	return -1;
  	}

    int OMID,dummyI;
    double x,y,z,dummyD;

    for (int i = 0; i < gNOMs; ++i)
    {
		inputFile >> OMID >> x >> y >> z >> dummyD >> dummyI;
		gOMpositions[i] = TVector3(x,y,z);
    }
    inputFile.close();	

	return 0;
}

int processMCCascades()
{
	TChain* mcFiles = new TChain("h11");
	mcFiles->Add("/Data/BaikalData/mc/DZH_cascades/ne16_tin_c*_00*.root");

	mcCascade* cascade = new mcCascade;
    
	mcFiles->SetBranchAddress("L0",&cascade->eventID);
	mcFiles->SetBranchAddress("Esh",&cascade->showerEnergy);
	mcFiles->SetBranchAddress("cost",&cascade->cosTheta);
	mcFiles->SetBranchAddress("fj",&cascade->phi);
	mcFiles->SetBranchAddress("jch",&cascade->nHits);
	mcFiles->SetBranchAddress("xtr",cascade->position);
	mcFiles->SetBranchAddress("Npmt",cascade->chID);
	mcFiles->SetBranchAddress("tre",cascade->time);
	mcFiles->SetBranchAddress("are",cascade->charge);

	if (ReadLogTable() == -1)
    {
    	std::cout << "Problem with 4D LogLikelihood file!" << std::endl;
    	return -1;
    }

    PrintConfig();
    cout << "RunInfo = Number of entries & runTime [hours] & runTime [days]" << endl;
    cout << "! " << mcFiles->GetEntries() << " " << "------" << "  " << "------" << endl; 
    std::cout << std::string(80,'*') << std::endl;

    if (ReadGeometryMCCascades() == -1)
	{
		std::cout << "Problem with geometry file!" << std::endl;
		return -1;
	}

    DoTheMagicMCCascades(mcFiles,cascade);
    SaveHistograms();

    return 0;
}

int main(int argc, char** argv) 
{
	clock_t begin = clock();
    // Init should be called at the beggining of all BARS programms
    App::Init(argc, argv, 0, parseOpts, readRC, checkParams);
    SetFitter();

    // processMCCascades();
    if (!gMCMu && !gMCNu && !gMCCas)
    {
    	processExperimentalData();
    	clock_t end = clock();
    	cout << endl << "Elapsed time : " << double(end - begin)/CLOCKS_PER_SEC << endl;
		return 0;    	
    }
    if (gMCMu || gMCNu)
    {
	    processMCData();
	    clock_t end = clock();
	    cout << endl << "Elapsed time : " << double(end - begin)/CLOCKS_PER_SEC << endl;
	    return 0;
    }
    if (gMCCas)
    {
    	processMCCascades();
    	clock_t end = clock();
    	cout << endl << "Elapsed time : " << double(end - begin)/CLOCKS_PER_SEC << endl;
    	return 0;
    }
    return -1;
}
