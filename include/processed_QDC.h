//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Feb  6 13:16:29 2019 by ROOT version 5.34/38
// from TTree processed/Processed rabbit data
// found on file: data_processed/RABBIT_DEC18_12.root
//////////////////////////////////////////////////////////

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class processed_QDC {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   TFile          *file;
   TFile          *rawfile;

   // Declaration of leaf types
   Float_t         cycle_time;
   Float_t         nmon_PSD;

   Int_t           ADC_short[16];
   Int_t           ADC_long[16];
   Int_t           TDC[16];
   Bool_t          overflow[16];
   Int_t           time_stamp;
   Double_t        seconds;

   // List of branches
   TBranch        *b_cycle_time;   //!
   TBranch        *b_nmon_PSD;   //!

   TBranch        *b_ADC_long;   //!
   TBranch        *b_ADC_short;   //!
   TBranch        *b_TDC;   //!
   TBranch        *b_overflow;   //!
   TBranch        *b_seconds;   //!
   TBranch        *b_time_stamp;   //!

   processed_QDC(int run_num);
   virtual ~processed_QDC();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual void     Init();
   virtual void     Show(Long64_t entry = -1);
};

processed_QDC::processed_QDC(int run_num)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   file = new TFile(Form("data_processed/RABITTS_%i.root", run_num));
   file->GetObject("processed_QDC", fChain);
   if (run_num<10){
      fChain->AddFriend("MDPP16_QDC", Form("data_root/RABITTS_00%i.root", run_num));
      rawfile = new TFile(Form("data_root/RABITTS_00%i.root", run_num));
   }
   else if (run_num<100){
      fChain->AddFriend("MDPP16_QDC", Form("data_root/RABITTS_0%i.root", run_num));
      rawfile = new TFile(Form("data_root/RABITTS_0%i.root", run_num));
   }
   else{
      fChain->AddFriend("MDPP16_QDC", Form("data_root/RABITTS_%i.root", run_num));
      rawfile = new TFile(Form("data_root/RABITTS_%i.root", run_num));
   }

   Init();
}

processed_QDC::~processed_QDC()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t processed_QDC::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

void processed_QDC::Init()
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.

   // Set branch addresses and branch pointers
   if (!fChain) return;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("cycle_time", &cycle_time, &b_cycle_time);
   fChain->SetBranchAddress("nmon_PSD", &nmon_PSD, &b_nmon_PSD);

   fChain->SetBranchAddress("ADC_short[16]", ADC_short, &b_ADC_short);
   fChain->SetBranchAddress("ADC_long[16]", ADC_long, &b_ADC_long);
   fChain->SetBranchAddress("TDC[16]", TDC, &b_TDC);
   fChain->SetBranchAddress("overflow[16]", overflow, &b_overflow);
   fChain->SetBranchAddress("time_stamp", &time_stamp, &b_time_stamp);
   fChain->SetBranchAddress("seconds", &seconds, &b_seconds);

}

void processed_QDC::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
