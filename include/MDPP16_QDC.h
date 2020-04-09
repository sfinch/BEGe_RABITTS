//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jan 30 13:46:18 2019 by ROOT version 5.34/38
// from TTree MDPP16_SCP/MDPP16 data
// found on file: data_root/RABBIT_DEC18_010.root
// Then heavily modified by S. Finch
// Class for raw data from MDPP16_QDC tree, produced by mvme2root
//////////////////////////////////////////////////////////

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TVectorT.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class MDPP16_QDC {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   TFile          *file;

   //file variables
   TVectorT<double> *m;
   TVectorT<double> *b;

   // Declaration of leaf types
   Int_t           ADC_long[16];
   Int_t           ADC_short[16];
   Int_t           TDC[16];
   Bool_t          overflow[16];
   Int_t           time_stamp;
   Int_t           extendedtime;
   Double_t        seconds;

   // List of branches
   TBranch        *b_ADC_long;   //!
   TBranch        *b_ADC_short;   //!
   TBranch        *b_TDC;   //!
   TBranch        *b_time_stamp;   //!
   TBranch        *b_extendedtime;   //!
   TBranch        *b_overflow;   //!
   TBranch        *b_pileup;   //!
   TBranch        *b_seconds;   //!

   //functions
   MDPP16_QDC(int run_num);
   virtual ~MDPP16_QDC();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual void     Init();
   virtual void     Show(Long64_t entry = -1);
};

MDPP16_QDC::MDPP16_QDC(int run_num)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
    if (run_num<10){
        file = new TFile(Form("data_root/RABITTS_00%i.root", run_num));
    }
    else if (run_num<100){
        file = new TFile(Form("data_root/RABITTS_0%i.root", run_num));
    }
    else{
        file = new TFile(Form("data_root/RABITTS_%i.root", run_num));
    }
    file->GetObject("MDPP16_QDC",fChain);

    Init();
}

MDPP16_QDC::~MDPP16_QDC()
{
    if (!fChain) return;
    delete fChain->GetCurrentFile();
}

Int_t MDPP16_QDC::GetEntry(Long64_t entry)
{
    // Read contents of entry.
    if (!fChain) return 0;
    return fChain->GetEntry(entry);
}

void MDPP16_QDC::Init()
{
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
 
    // Set branch addresses and branch pointers
    if (!fChain) return;
    fChain->SetMakeClass(1);
 
    fChain->SetBranchAddress("ADC_short[16]", ADC_short, &b_ADC_short);
    fChain->SetBranchAddress("ADC_long[16]", ADC_long, &b_ADC_long);
    fChain->SetBranchAddress("TDC[16]", TDC, &b_TDC);
    fChain->SetBranchAddress("overflow[16]", overflow, &b_overflow);
    fChain->SetBranchAddress("time_stamp", &time_stamp, &b_time_stamp);
    fChain->SetBranchAddress("extendedtime", &extendedtime, &b_extendedtime);
    fChain->SetBranchAddress("seconds", &seconds, &b_seconds);
}

void MDPP16_QDC::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
    if (!fChain) return;
    fChain->Show(entry);
}

