#ifndef MuEtree_H
#define MuEtree_H

///////////////////////////////////////////////
// Data Formats for MuE MC events
//
// G.Abbiendi  5/Jul/2018 
///////////////////////////////////////////////

#include "TROOT.h"
#include <vector>
#include <string>

namespace MuE {

  class MCstat {

  public:
    Long64_t Nevgen;
    Long64_t Nwgt, Nwgt_Negative;
    Double_t Swgt, Swgt_Negative;
    Double_t SQwgt, SQwgt_Negative;
    Long64_t Nwgt_OverMax;
    Double_t WmaxTrue;
    Double_t Xsec, XsecErr;
    Double_t Xsec_Negative, Xsec_Negative_Err;
    Double_t Xsec_OverMax, Xsec_OverMax_Err;
   
    MCstat(){};
    virtual ~MCstat(){};
    
    ClassDef(MCstat,2)
  };

  class Setup {

  public:
    std::string program_version;
    Int_t process_ID;
    std::string running_on;
    std::string start_time;
    Long64_t Nevreq;
    Bool_t UNWGT;
    std::string Mode;
    UInt_t rnd_ext, rnd_int;
    Double_t Ebeam, EbeamRMS;
    Double_t charge_mu;
    Double_t mass_mu, mass_e;
    Double_t invalfa0;
    Double_t k0cut; 
    Double_t Emin_e; 
    Double_t Wnorm; 
    Double_t Wmax;

    MCstat MCsums;

    Setup(){};
    virtual ~Setup(){};

    ClassDef(Setup,4)
  };

  class P4 {

  public:
    Double_t E;
    Double_t px;
    Double_t py;
    Double_t pz;
    
    P4(){};
    virtual ~P4(){};
    
    ClassDef(P4,2)
  };

  class Event {

  public:    
    UInt_t RunNr;
    Long64_t EventNr;
    Double_t wgt_full, wgt_norun, wgt_lep, wgt_LO;
    Double_t E_mu_in; 
    P4 P_mu_out;
    P4 P_e_out;
    std::vector<P4> photons;
    
    Event(){};
    virtual ~Event(){};
    
    ClassDef(Event,4)
  };

}

#endif
