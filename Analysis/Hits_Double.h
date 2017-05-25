#ifndef Hits_h
#define Hits_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TVector3.h>
#include <TMath.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <iostream>

using namespace std;
using namespace TMath;

// Header file for the classes stored in the TTree if any.

class Hits {
 public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain
  
  // New Variables
  TString  fileName;
  TString  dataDirectory;
  TString  filePath;
  
  Int_t    verbose  = 0;
  
  Long64_t nEvents  = 0;
  Long64_t maxEntry = 0;
  
  Int_t    eventID_prevEntry = -1;
  Int_t    eventEntry = 0;
  
  Int_t    photon1Entry = 0;
  Int_t    photon2Entry = 0;

  Float_t  meanEntriesPerEvent = 0.;
  
  Int_t    nPhot = 0;
  Int_t    nRay  = 0;
  Int_t    nComp  = 0;
  
  Float_t  fracRayOnly  = 0.;
  Float_t  fracPhotOnly = 0.;
  
  Float_t  fracPhot = 0.;
  Float_t  fracComp = 0.;
  
  TString  processStr;
  TString  firstProcess;
  
  Bool_t   fstHitGood = kFALSE;
  Bool_t   sndHitGood = kFALSE;
  
  Float_t  eTotMod1 = 0.;
  Float_t  eTotMod2 = 0.;
  
  Bool_t   goodPhi;
  
  Float_t  phi,theta;
  
  // vectors for calculating phi
  
  TVector3 vBeam_g1, vBeam_g2;
  TVector3 p0_g1,p0_g2;
  TVector3 p1_g1,p1_g2;
  TVector3 p2_g1,p2_g2;
  TVector3 vScat_g1,vScat_g2;

  TH1F * hPhi  = new TH1F("hPhi","hPhi;#phi (deg);Events",
			   36,-180,180);
  TH1F * hTheta = new TH1F("hTheta","hTheta;#theta (deg);Events",
			   36,0,180);
  
  // leaf types
  Int_t           PDGEncoding;
  Int_t           trackID;
  Int_t           parentID;
  Double_t        time;
  Float_t         edep;
  Float_t         stepLength;
  Float_t         posX;
  Float_t         posY;
  Float_t         posZ;
  Float_t         localPosX;
  Float_t         localPosY;
  Float_t         localPosZ;
  Int_t           gantryID;
  Int_t           moduleID;
  Int_t           clusterID;
  Int_t           pixelID;
  Int_t           unused4ID;
  Int_t           unused5ID;
  Int_t           photonID;
  Int_t           nPhantomCompton;
  Int_t           nCrystalCompton;
  Int_t           nPhantomRayleigh;
  Int_t           nCrystalRayleigh;
  Int_t           primaryID;
  Float_t         sourcePosX;
  Float_t         sourcePosY;
  Float_t         sourcePosZ;
  Int_t           sourceID;
  Int_t           eventID;
  Int_t           runID;
  Float_t         axialPos;
  Float_t         rotationAngle;
  Int_t           volumeID[10];
  Char_t          processName[15];
  Char_t          comptVolName[5];
  Char_t          RayleighVolName[5];
  
  // List of branches
  TBranch        *b_PDGEncoding;   //!
  TBranch        *b_trackID;   //!
  TBranch        *b_parentID;   //!
  TBranch        *b_time;   //!
  TBranch        *b_edep;   //!
  TBranch        *b_stepLength;   //!
  TBranch        *b_posX;   //!
  TBranch        *b_posY;   //!
  TBranch        *b_posZ;   //!
  TBranch        *b_localPosX;   //!
  TBranch        *b_localPosY;   //!
  TBranch        *b_localPosZ;   //!
  TBranch        *b_gantryID;   //!
  TBranch        *b_moduleID;   //!
  TBranch        *b_clusterID;   //!
  TBranch        *b_pixelID;   //!
  TBranch        *b_unused4ID;   //!
  TBranch        *b_unused5ID;   //!
  TBranch        *b_photonID;   //!
  TBranch        *b_nPhantomCompton;   //!
  TBranch        *b_nCrystalCompton;   //!
  TBranch        *b_nPhantomRayleigh;   //!
  TBranch        *b_nCrystalRayleigh;   //!
  TBranch        *b_primaryID;   //!
  TBranch        *b_sourcePosX;   //!
  TBranch        *b_sourcePosY;   //!
  TBranch        *b_sourcePosZ;   //!
  TBranch        *b_sourceID;   //!
  TBranch        *b_eventID;   //!
  TBranch        *b_runID;   //!
  TBranch        *b_axialPos;   //!
  TBranch        *b_rotationAngle;   //!
  TBranch        *b_volumeID;   //!
  TBranch        *b_processName;   //!
  TBranch        *b_comptVolName;   //!
  TBranch        *b_RayleighVolName;   //!
  
  Hits(TTree *tree=0);
  virtual ~Hits();
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop(Int_t userVerbose = 0,
			Long64_t maxEvents = -1);
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);
  
  void InitNewEvent();
  void EndEvent();
  void EndPreviousEvent();
  void PrintEntry(Long64_t entryNumber);
  void PrintEventStats();
  void PrintTVector3Elements(TString,
			     TVector3);
  void FinalStats();
  void SumEnergyDep();
  void CountProcesses(TString);
};

#endif

#ifdef Hits_cxx
Hits::Hits(TTree *tree) : fChain(0) 
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  
  fileName      = "CZT_Single_No_Pol_output.root";
  fileName      = "CZT_Single_Pol_output.root";
  
  fileName      = "CZT_Single_LivePol_output.root";
  //fileName      = "CZT_Single_StdOpt3_output.root";
  
  fileName      = "CZT_Double_LivePol_output.root";

  dataDirectory = "../Data/";
  
  filePath = dataDirectory + fileName;
  
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(filePath);
    if (!f || !f->IsOpen()) {
      f = new TFile(filePath);
    }
    f->GetObject("Hits",tree);

  }
  Init(tree);
}

Hits::~Hits()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Hits::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Hits::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Hits::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("PDGEncoding", &PDGEncoding, &b_PDGEncoding);
   fChain->SetBranchAddress("trackID", &trackID, &b_trackID);
   fChain->SetBranchAddress("parentID", &parentID, &b_parentID);
   fChain->SetBranchAddress("time", &time, &b_time);
   fChain->SetBranchAddress("edep", &edep, &b_edep);
   fChain->SetBranchAddress("stepLength", &stepLength, &b_stepLength);
   fChain->SetBranchAddress("posX", &posX, &b_posX);
   fChain->SetBranchAddress("posY", &posY, &b_posY);
   fChain->SetBranchAddress("posZ", &posZ, &b_posZ);
   fChain->SetBranchAddress("localPosX", &localPosX, &b_localPosX);
   fChain->SetBranchAddress("localPosY", &localPosY, &b_localPosY);
   fChain->SetBranchAddress("localPosZ", &localPosZ, &b_localPosZ);
   fChain->SetBranchAddress("gantryID", &gantryID, &b_gantryID);
   fChain->SetBranchAddress("moduleID", &moduleID, &b_moduleID);
   fChain->SetBranchAddress("clusterID", &clusterID, &b_clusterID);
   fChain->SetBranchAddress("pixelID", &pixelID, &b_pixelID);
   fChain->SetBranchAddress("unused4ID", &unused4ID, &b_unused4ID);
   fChain->SetBranchAddress("unused5ID", &unused5ID, &b_unused5ID);
   fChain->SetBranchAddress("photonID", &photonID, &b_photonID);
   fChain->SetBranchAddress("nPhantomCompton", &nPhantomCompton, &b_nPhantomCompton);
   fChain->SetBranchAddress("nCrystalCompton", &nCrystalCompton, &b_nCrystalCompton);
   fChain->SetBranchAddress("nPhantomRayleigh", &nPhantomRayleigh, &b_nPhantomRayleigh);
   fChain->SetBranchAddress("nCrystalRayleigh", &nCrystalRayleigh, &b_nCrystalRayleigh);
   fChain->SetBranchAddress("primaryID", &primaryID, &b_primaryID);
   fChain->SetBranchAddress("sourcePosX", &sourcePosX, &b_sourcePosX);
   fChain->SetBranchAddress("sourcePosY", &sourcePosY, &b_sourcePosY);
   fChain->SetBranchAddress("sourcePosZ", &sourcePosZ, &b_sourcePosZ);
   fChain->SetBranchAddress("sourceID", &sourceID, &b_sourceID);
   fChain->SetBranchAddress("eventID", &eventID, &b_eventID);
   fChain->SetBranchAddress("runID", &runID, &b_runID);
   fChain->SetBranchAddress("axialPos", &axialPos, &b_axialPos);
   fChain->SetBranchAddress("rotationAngle", &rotationAngle, &b_rotationAngle);
   fChain->SetBranchAddress("volumeID", volumeID, &b_volumeID);
   fChain->SetBranchAddress("processName", processName, &b_processName);
   fChain->SetBranchAddress("comptVolName", comptVolName, &b_comptVolName);
   fChain->SetBranchAddress("RayleighVolName", RayleighVolName, &b_RayleighVolName);
   Notify();
}

Bool_t Hits::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Hits::PrintEventStats(){
  
  if(verbose > 0 ){
    cout << " Event " << (nEvents-1) << " Stats "     << endl;
    cout << " eTotMod1        = " << eTotMod1     << endl;
    cout << " n event entries = " << eventEntry   << endl;
    
    if(sndHitGood){
      cout << " phi           = " << phi         << endl;
      cout << " theta         = " << theta       << endl;
    }
    
    cout << " --------------------------------------- " <<  endl;
  }
  
}

void Hits::PrintTVector3Elements(TString vName,
				 TVector3 v3){
  
  cout << endl;
  cout << vName << " =  (" 
       << v3.X() << ", " 
       << v3.Y() << ", "
       << v3.Z() << ")" 
       << endl;
  
}

void Hits::EndPreviousEvent(){
  EndEvent();
}

void Hits::EndEvent(){

  // print event statistics 
  if(eventID_prevEntry == -1)
    return;
  
  PrintEventStats();
  
  // iterate / accumulate event variables
  meanEntriesPerEvent += eventEntry;

  if(eventEntry == nPhot)
    fracPhotOnly = fracPhotOnly + 1.;
  else if(eventEntry == nRay)
    fracRayOnly += (Float_t) nRay;
  
  if(nPhot>0 && nComp == 0)
    fracPhot = fracPhot + 1.0;
  
  if(nComp > 0)
    fracComp = fracComp + 1.0;

  if(sndHitGood){
    hPhi->Fill(phi);
    hTheta->Fill(theta);
  }

}

void Hits::InitNewEvent(){
  
  eventID_prevEntry = eventID;
  
  nEvents++;
  
  //reset event by event variables
  eventEntry = 0;

  eTotMod1 = 0.0;
  eTotMod2 = 0.0;
  
  nRay  = 0;
  nPhot = 0;
  nComp = 0;

  fstHitGood  = kFALSE;
  sndHitGood = kFALSE;

  if( verbose > 0 ) {
    cout << " --------------------------------------- " <<  endl;
    cout << " Processing event " << eventID <<  endl;
  }
  
}

void Hits::FinalStats(){

  // Last iteration of loop
  EndEvent();

  // Calculate final variables
  meanEntriesPerEvent = (Float_t)meanEntriesPerEvent/nEvents;
  
  fracPhotOnly = fracPhotOnly/(Float_t)nEvents;
  fracRayOnly  = fracRayOnly/(Float_t)nEvents;
  fracPhot     = fracPhot/(Float_t)nEvents;
  fracComp     = fracComp/(Float_t)nEvents;
               
  // print full data set details
  cout <<  endl;
  cout << " --------------------------------------- " <<  endl;
  cout << " --------------------------------------- " <<  endl;
  cout << " ------------- Final Stats  ------------ " <<  endl;
  cout << "    "        << fileName                   <<  endl;  
  cout << " --------------------------------------- " <<  endl;
  cout << " nEvents        = " << nEvents             <<  endl;
  cout << " nentries       = " << maxEntry            <<  endl;
  cout << " hits per event = " << meanEntriesPerEvent <<  endl;
  cout << " fracPhot       = " << fracPhot            <<  endl; 
  cout << " fracComp       = " << fracComp            <<  endl; 
  cout << " fracPhotOnly   = " << fracPhotOnly        <<  endl; 
  cout << " fracRayOnly    = " << fracRayOnly         <<  endl; 
  cout << " --------------------------------------- " <<  endl;
  cout << " --------------------------------------- " <<  endl;

  meanEntriesPerEvent = 0.;
  fracPhotOnly = 0;
  fracRayOnly  = 0;
  fracPhot     = 0;
  fracComp     = 0;
  nEvents      = 0;
  eventID_prevEntry = -1;
  
  hPhi->SetMinimum(0);
  
  TCanvas * canvas =  new TCanvas("canvas","canvas",800,400);
  canvas->Divide(2,1);

  canvas->cd(1);
  hPhi->Draw();
  canvas->cd(2);
  hTheta->Draw();
}

void Hits::PrintEntry(Long64_t entryNum){

  
  cout << " -------------------------- "   << endl;
  cout << " eventID     = " << eventID     << endl;
  cout << " entry       = " << entryNum    << endl;
  cout << " eventEntry  = " << eventEntry  << endl;
  cout << " process     = " << processName << endl;
  cout << " trackID     = " << trackID     << endl;
  cout << " moduleID    = " << moduleID    << endl;
  cout << " photonID    = " << photonID    << endl;
  cout << " parentID    = " << parentID    << endl;
  cout << " edep        = " << edep        << endl;
  cout << " -------------------------- "   << endl;

}

void Hits::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

void Hits::SumEnergyDep(){
  
  if     (moduleID==0){
    eTotMod1 += edep;
  }
  else if(moduleID==1){
    eTotMod2 += edep;
  }
}

void Hits::CountProcesses(TString processStr){
  
  
  if     (processStr=="Rayl")
    nRay++;
  else if(processStr=="phot")
    nPhot++;
  else if(processStr=="compt")
    nComp++;
  else {
    if(verbose > 0)
      cout << " Process = " << processStr << endl; 
  }
}

#endif // #ifdef Hits_cxx
