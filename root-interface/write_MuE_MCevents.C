/////////////////////////////////////////////////////////
//
// Read MuE MC events from ASCII file format 
// and write them in root format
//
//   Input parameters: filenames for
//   1: ASCII input (header, generated events, footer) 
//   2: root output
//
// G.Abbiendi  5/Jul/2018
/////////////////////////////////////////////////////////

#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TTree.h"
#include "MuEtree.h"

#include <cstdlib>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

// set it to true to debug
bool DEBUG = false;

istringstream input_line (ifstream & input_file) {
  string line;
  getline(input_file, line); 
  istringstream stream(line);
  if (DEBUG) cout << "\t" << stream.str() << endl;
  return stream;
}


int main(int argc, char* argv[]) {

  if (argc != 3) {
    cout << "Usage : "
	 << argv[0] << "  input_ASCII_filename  output_ROOT_filename \n";
    exit(100);
  }

  ifstream input_file(argv[1]);
  TString ofile(argv[2]);
  TFile *output_file = new TFile(ofile,"RECREATE");

  cout<<endl<<"Reading MuE MC events from input file: "<< argv[1] <<endl;
  cout      <<"Writing them to output root file     : "<< argv[2] <<endl;
  
  // Set it to kTRUE if you do not run interactively
  gROOT->SetBatch(kTRUE); 

  // Initialize Root application
  TRint* app = new TRint("Root Application", &argc, argv);

  Int_t splitBranches = 2;

  MuE::Setup mc_params;
  TTree *partree = new TTree("MuEsetup","MuE MC parameters");
  partree->Branch("MuEparams",&mc_params,64000,splitBranches);

  MuE::Event event;
  TTree *tree = new TTree("MuEtree","MuE MC tree");
  tree->Branch("MuE",&event,64000,splitBranches);
  tree->SetAutoSave(500000);

  // Read the header section
  string line, key, dump, str1, str2;
  istringstream stream;

  cout << "Start reading header section" <<endl;
  stream = input_line(input_file);
  stream >> key;
  if (key != "<header>") {
    cout << "*** ERROR: unexpected format for header section." << endl;
    return 100;
  }
  
  stream = input_line(input_file);
  stream >> dump >> dump >> str1 >> str2;
  mc_params.program_version = str1+" "+str2;
  
  stream = input_line(input_file);
  stream >> dump >> dump >> mc_params.process_ID >> dump >> dump >> mc_params.running_on;

  stream = input_line(input_file);
  mc_params.start_time = stream.str();

  // Number of requested events
  stream = input_line(input_file);
  stream >> dump >> mc_params.Nevreq;

  // Unweighted (wgt=1) events ?  
  stream = input_line(input_file);
  stream = input_line(input_file);
  stream >> mc_params.UNWGT;

  // generator Mode
  stream = input_line(input_file);
  stream >> dump >> dump >> mc_params.Mode;

  // Initial seed for Random numbers 
  stream = input_line(input_file);
  stream = input_line(input_file);
  stream >> mc_params.rnd_ext >> mc_params.rnd_int;

  stream = input_line(input_file);
  
  // Nominal Muon Beam energy
  stream = input_line(input_file);
  stream = input_line(input_file);  
  stream >> mc_params.Ebeam; 

  // Beam energy Gaussian spread 
  stream = input_line(input_file);
  stream = input_line(input_file);
  stream >> mc_params.EbeamRMS;

  // muon charge
  stream = input_line(input_file);
  stream = input_line(input_file);
  stream >> mc_params.charge_mu;

  // muon mass
  stream = input_line(input_file);
  stream = input_line(input_file);
  stream >> mc_params.mass_mu;

  // electron mass
  stream = input_line(input_file);
  stream = input_line(input_file);
  stream >> mc_params.mass_e;

  // (fine structure constant)^-1  
  stream = input_line(input_file);
  stream = input_line(input_file);
  stream >> mc_params.invalfa0;  

  // photon cutoff (technical parameter)
  stream = input_line(input_file);
  stream >> dump >> dump >> dump >> dump >> mc_params.k0cut;

  // Minimum Energy of outgoing electron
  stream = input_line(input_file);
  stream >> dump >> dump >> dump >> dump >> dump >> dump >> mc_params.Emin_e;

  stream = input_line(input_file);
  
  // Cross section normalization factor 
  stream = input_line(input_file);
  stream = input_line(input_file);  
  stream >> mc_params.Wnorm; 

  // Assumed maximum weight (for unweighted generation)
  stream = input_line(input_file);
  stream = input_line(input_file);  
  stream >> mc_params.Wmax; 
  
  stream = input_line(input_file);
  stream >> key;
  if (key == "</header>") {
      cout << "End reading header section." <<endl;
  } else {
    cout << "*** ERROR: unexpected format for the header section." << endl;
    return 200;
  }

  // Read mu-e events 
  UInt_t nrows, nphot;
  cout << "Start reading MuE events" <<endl;

  while (getline(input_file, line)) {
    stream = istringstream(line);
    stream >> key;
    if (key == "<footer>") break;
    else if (key != "<event>") {
      cout << "*** ERROR: unexpected format for input event (begin)" << endl;
      return 300;
    }

    // read run number
    stream = input_line(input_file);
    stream >> event.RunNr;

    // read event number
    stream = input_line(input_file);
    stream >> event.EventNr;

    // read number of final state particles (-> nr of photons)
    stream = input_line(input_file);
    stream >> nrows;
    nphot = nrows!=0 ? nrows-2 : 0;

    // read MC weights (with full/no/leptonic running)
    stream = input_line(input_file);
    stream  >> event.wgt_full >> event.wgt_norun >> event.wgt_lep;

    // read LO weight (with full running)
    stream = input_line(input_file);
    stream >> event.wgt_LO;

    // read incoming muon energy
    stream = input_line(input_file);
    stream >> event.E_mu_in;

    // read scattered muon
    stream = input_line(input_file);
    stream >> event.P_mu_out.E >> event.P_mu_out.px >> event.P_mu_out.py >> event.P_mu_out.pz ;

    // read scattered electron
    stream = input_line(input_file);
    stream >> event.P_e_out.E  >> event.P_e_out.px  >> event.P_e_out.py  >> event.P_e_out.pz ;
    
    // read photons (when nphot>0)
    event.photons.clear();
    for (int i=0; i<nphot; ++i) {
      MuE::P4 P_photon;
      stream = input_line(input_file);
      stream >> P_photon.E >> P_photon.px >> P_photon.py >> P_photon.pz;
      event.photons.push_back(P_photon);
    }

    stream = input_line(input_file);
    stream >> key;
    if (key != "</event>") {
      cout << "*** ERROR: unexpected format for input event (end)" << endl;
      return 400;
    }
    
    if (input_file.good()) {   
      // fill the event tree
      tree->Fill();
      if (event.EventNr % 100000 == 0 ) cout << " processed event : " << event.EventNr << endl;
    } 
    else {
      cout<<"*** ERROR on reading: EventNr = "<<event.EventNr<<", n photons = "<<event.photons.size() <<endl;
      return 500;
    }
  }

  cout << "End reading MuE events." <<endl;
  if (DEBUG) tree->Print();

  // save the number of generated events (including zero-weight events) among the final parameters in the parameter tree 
  Long64_t NevTot = tree->GetEntries();
  mc_params.MCsums.Nevgen = NevTot;

  cout << "Start reading the footer section" <<endl;

  // estimated cross section for the generated process
  stream = input_line(input_file);
  stream = input_line(input_file);
  stream >> mc_params.MCsums.Xsec >> dump >> mc_params.MCsums.XsecErr;

  // total number of weights
  stream = input_line(input_file);
  stream = input_line(input_file);
  stream >> mc_params.MCsums.Nwgt;

  // true maximum weight at the end 
  stream = input_line(input_file);
  stream = input_line(input_file);
  stream >> mc_params.MCsums.WmaxTrue;

  // number of weights greater than the assumed maximum (Wmax)
  stream = input_line(input_file);
  stream = input_line(input_file);
  stream >> mc_params.MCsums.Nwgt_OverMax;

  // number of negative weights
  stream = input_line(input_file);
  stream = input_line(input_file);
  stream >> mc_params.MCsums.Nwgt_Negative;

  // estimated bias on cross section due to weights greater than Wmax (in unweighted generation)
  stream = input_line(input_file);
  stream = input_line(input_file);
  stream >> mc_params.MCsums.Xsec_OverMax >> dump >> mc_params.MCsums.Xsec_OverMax_Err;

  // estimated cross section contribution from negative weights
  stream = input_line(input_file);
  stream = input_line(input_file);
  stream >> mc_params.MCsums.Xsec_Negative >> dump >> mc_params.MCsums.Xsec_Negative_Err;

  // final sum of weights and squared weights
  stream = input_line(input_file);
  stream = input_line(input_file);
  stream >> mc_params.MCsums.Swgt >> mc_params.MCsums.SQwgt;

  // final sum of negative weights and squared weights
  stream = input_line(input_file);
  stream = input_line(input_file);
  stream >> mc_params.MCsums.Swgt_Negative >> mc_params.MCsums.SQwgt_Negative;

  stream = input_line(input_file);
  stream >> key;
  if (key == "</footer>") {
      cout << "End reading the footer section." <<endl;
  } else {
    cout << "*** ERROR: unexpected format for the footer section." << endl;
    return 600;
  }

  // fill the parameter tree (1 entry per mc run)
  if (!input_file.fail()) partree->Fill();
  else {
    cout<<"*** ERROR on reading the footer section."<<endl;
    return 700;
  }
  if (DEBUG) partree->Print();

  output_file->Write();
  output_file->Close();
  cout << "Generated " << NevTot << " events, written to file: "<< ofile << endl;
  cout << "Done." << endl;

  if (!gROOT->IsBatch()) app->Run();

  return 0;
}

