//
// read MuE MC events from root file format 
// and make some simple analysis
//
//   Input parameters: 
//   1: MC root input filename (input events and parameters)
//   2: MC root output filename (analysis results)
//

#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TTree.h"
#include "MuEtree.h"

#include "TH1D.h"
#include "TLorentzVector.h"
#include "TMath.h"

#include <string>
#include <fstream>
#include <iostream>

//// to debug numerical precision
int counter = 0;
double sum = 0.;
double maxt13 = 0.;
////

using namespace std;

namespace MuE {

  class MyHistos {
    
  public:
    MyHistos(TFile *outFile);
    ~MyHistos();
    void Fill(const MuE::Event & event, const MuE::MCpars & params);
    void Normalize(const Long64_t& nevtot, const MuE::MCpars & params);
    
    TH1D* h_t13;
    TH1D* h_t24;
    TH1D* h_egamma;
    TH1D* h_eth;    
  };
  
  MyHistos::MyHistos(TFile *outFile) {
    
    outFile->cd();
    
    // store the sum of squares of weights for all the histos to be created
    TH1::SetDefaultSumw2();

    h_t13 = new TH1D("h_t13","(pmu_in - pmu_out)^2", 100, -0.143, 0.);
    h_t24 = new TH1D("h_t24","(pe_in - pe_out)^2", 100, -0.143, 0.);
    h_egamma = new TH1D("h_egamma","photon energy", 100, 0., 150.);
  }
  
  MyHistos::~MyHistos() {
    delete h_t13, h_t24, h_egamma, h_eth;
  }
  
  void MyHistos::Fill(const MuE::Event & event, const MuE::MCpars & params) {
    Double_t evwgt = event.wgt_full;
    if (evwgt != 0.) {
      //CARLO      evwgt = 1.;
      Double_t m_e = params.mass_e;
      TLorentzVector p_e_in(0.,0.,0.,m_e);
      Double_t m_mu = params.mass_mu;
      Double_t ebeam = event.EmuBeam;
      Double_t pbeam = sqrt(ebeam*ebeam-m_mu*m_mu);
      TLorentzVector p_mu_in(0.,0.,pbeam,ebeam);
      TLorentzVector p_mu_out(event.P_mu_out.px, event.P_mu_out.py, event.P_mu_out.pz, event.P_mu_out.E);
      TLorentzVector p_e_out(event.P_e_out.px, event.P_e_out.py, event.P_e_out.pz, event.P_e_out.E);
      int n_photons = event.photons.size();
      Double_t egamma = n_photons >0 ? event.photons[0].E : 0.;
      Double_t t13 = (p_mu_in - p_mu_out)*(p_mu_in - p_mu_out);
      //////
      if (t13 >= 0.) {
	counter++;
	cout << endl << "t13 = "<< t13;
	if (t13 > maxt13) maxt13 = t13;
	sum+=t13;
      }
      //////
      Double_t t24 = (p_e_in - p_e_out)*(p_e_in - p_e_out);
      
      h_t13->Fill(t13, evwgt);
      h_t24->Fill(t24, evwgt);
      if (n_photons > 0) h_egamma->Fill(egamma, evwgt);
    }
  }

  void MyHistos::Normalize(const Long64_t& nevtot, const MuE::MCpars & params) {
    Double_t wnorm = params.XsecNorm;
    ////
    cout << "nevtot = " << nevtot << endl;
    Double_t fnevtot = (Double_t)nevtot;
    cout << "fnevtot= " << fnevtot << endl;
    ////
    Double_t norm = wnorm / fnevtot;
    h_t13->Scale(norm, "width");
    h_t24->Scale(norm, "width");
    h_egamma->Scale(norm, "width");
  }   
}

int main(int argc, char* argv[]) {

  //using namespace std;
  using namespace MuE;

  if (argc != 3) {
    cout << "Usage : "
	 << argv[0] << " PATH_TO_INPUT_FILE PATH_TO_OUTPUT_FILE \n";
    exit(100);
  }

  TString ifile(argv[1]);
  TString ofile(argv[2]);
  TFile *output_file = new TFile(ofile,"RECREATE");

  cout<<endl<<"Reading MuE MC events from input file   : "<< argv[1] <<endl;
  cout      <<"Writing analysis results to output file : "<< argv[2] <<endl;

  // Set it to kTRUE if you do not run interactively
  gROOT->SetBatch(kTRUE); 

  // Initialize Root application
  TRint* app = new TRint("Root Application", &argc, argv);

  MyHistos* histos = new MyHistos(output_file);

  MuE::MCpars* params = new MuE::MCpars(); 
  MuE::Event* event = new MuE::Event();
  
  TFile* inputFile = TFile::Open(ifile,"READONLY");

  // get branch with constant parameters
  TTree* partree = (TTree*)inputFile->Get("MuE_params");
  TBranch* parbranch = partree->GetBranch("MuEpars");
  parbranch->SetAddress(&params);
  parbranch->GetEntry();
  partree->Show();

  TTree* tree = (TTree*)inputFile->Get("MuEtree");
  TBranch* branch = tree->GetBranch("MuE");
  branch->SetAddress(&event);

  Long64_t nevtot = tree->GetEntries();
  cout << endl << "Total number of events in tree = " << nevtot << endl << endl;

  /////
  maxt13 = 0.;
  /////
  for (Long64_t iEvent=0; iEvent<nevtot; ++iEvent) {

    if (tree->LoadTree(iEvent)<0) {
      cout << "*** Event "<< iEvent <<" not found in tree." <<endl;
      break;
    }

    if (iEvent % 100000 == 0 ) cout << " processing event : " << iEvent << endl;
    
    branch->GetEntry(iEvent);

    if (iEvent<5) tree->Show();
    histos->Fill(*event, *params);
  }
  
  histos->Normalize(nevtot, *params);
  //////
  cout << "counter = " << counter << endl;
  cout << "average t13 for the positive fraction = " << sum/counter << endl;
  cout << "max t13 for the positive fraction = " << maxt13 << endl;
  //////
  delete event;
  delete params;
  inputFile->Close();

  output_file->Write();
  
  if (!gROOT->IsBatch()) app->Run();
  
  return 0;
}
