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
#include <iomanip> 

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
  
  Int_t    verbose   = 0;
  Float_t  thickness = 0.;

  Long64_t nEvents  = 0;

  Int_t    nPrimaries = 0;

  Int_t    eventID_prevEntry = -1;
  Int_t    eventEntry   = -1;

  Int_t    photon0Entry = 0;
  Int_t    photon1Entry = 0;

  Long64_t nHits = 0;
  Int_t    entriesInEvent;
  Double_t meanHitsPerEvent = 0.;
  
  Int_t    nPhot  = 0;
  Int_t    nRay   = 0;
  Int_t    nComp  = 0;
  Int_t    nOther = 0;
  
  Long64_t nPhotTot  = 0;
  Long64_t nRayTot   = 0;
  Long64_t nCompTot  = 0;
  Long64_t nOtherTot = 0;

  Float_t  fracPhotFirst  = 0.;
  Float_t  fracRayFirst   = 0.;
  Float_t  fracCompFirst  = 0.;
  Float_t  fracOtherFirst = 0.;

  Float_t  fracPhotOnly  = 0.;
  Float_t  fracRayOnly   = 0.;
  Float_t  fracCompOnly  = 0.;
  
  Float_t  fracPhot   = 0.;
  Float_t  fracRay    = 0.;
  Float_t  fracComp   = 0.;
  Float_t  fracOther  = 0.;
  
  Float_t  fracPET = 0.;
  Float_t  fracQET = 0.;
  Float_t  fracE_LT_511 = 0.;

  TString  processStr;
  
  Bool_t   edepTot511 = kFALSE;
  Bool_t   goodPET = kFALSE;
  Bool_t   goodQET = kFALSE;

  Bool_t   fstHitGood = kFALSE;
  Bool_t   sndHitGood = kFALSE;
  
  Bool_t   photFirst = kFALSE;
  Bool_t   compFirst = kFALSE;
  
  Float_t  eTotMod1 = 0.;
  Float_t  eTotMod2 = 0.;
  
  static const Int_t nModules = 1;
  
  
  Float_t  phi[nModules],theta[nModules];
  Float_t  dR[nModules];
  Float_t  dXY[nModules];
  Float_t  dZ[nModules];
  
  // vectors for calculating theta, phi
  TVector3 vBeam[nModules];
  TVector3 p0[nModules];
  TVector3 p1[nModules];
  TVector3 p2[nModules];
  TVector3 vScat[nModules];

  TH1F * hPhi[nModules];
  TH1F * hTheta[nModules];
  TH1F * hdR[nModules];
  
  TH1F * hZ1Phot[nModules];
  TH1F * hZ1Comp[nModules];

  TH2F * hdZdXY[nModules];
  TH2F * hThetadXY[nModules];
  TH2F * hThetadR[nModules];

  TString hName;
  TString hTitle;
  
  Int_t   nBins  = 36;
  //Float_t xRange[2] = {-180., 180.};
  
  
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
  
  //Hits(TTree *tree=0);
  Hits();
  virtual ~Hits();
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop(Long64_t maxEvents = -1);
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);

  void InitHistos();
  void InitNewEvent();
  void EndEvent(Bool_t);
  void EndPreviousEvent();
  void PrintEntry(Long64_t entryNumber);
  void PrintEventStats(Bool_t);
  void PrintTVector3Elements(TString,
			     TVector3);
  void SetQETVariables();

  Bool_t HitDet(Int_t);

  TString GetFilePath();
  void SetFile(TString);
  void Run(Float_t  userThickness,
	   Int_t    userVerbose,
	   Long64_t maxEvents);
  void FinalStats();
  void SumEnergyDep();
  void CountProcesses(TString);
  void SetVerbose(Int_t);
  void SetThickness(Float_t);
  void PlotHistograms();
  void EndRun();
  void ResetRunVariables();
  
};

#endif

#ifdef Hits_cxx
Hits::Hits() : fChain(0) 
{

  SetFile("emstd_opt3");

}

Hits::~Hits()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}


Bool_t Hits::HitDet(Int_t hit){
  
  if     (hit==1)
    return p1[moduleID].Z() < (thickness + 50.);
  else if(hit==2)
    return p2[moduleID].Z() < (thickness + 50.);
  else 
    return false;
}

void Hits::SetQETVariables(){

  
  p1[moduleID].SetXYZ(posX,posY,posZ);

  // first hit of photon
  if     (eventEntry == 0 &&
	  processStr=="compt" ){
    
    p0[moduleID].SetXYZ(sourcePosX,
			sourcePosY,
			sourcePosZ);
    
    
    
    
    if(verbose > 2){
      PrintTVector3Elements(" p0 ",p0[moduleID]);
      PrintTVector3Elements(" p1 ",p1[moduleID]);
    }
    
    fstHitGood = kTRUE;
  }
  // second hit of photon
  else if(eventEntry == 1      && 
	  (processStr=="compt" ||
	   processStr=="phot") &&
	  fstHitGood ){
    
    p2[moduleID].SetXYZ(posX,posY,posZ);
    
    vBeam[moduleID] = p1[moduleID] - p0[moduleID];
    vScat[moduleID] = p2[moduleID] - p1[moduleID];
    
    phi[moduleID] = vScat[moduleID].DeltaPhi(vBeam[moduleID]);
    phi[moduleID] = RadToDeg()*phi[moduleID];
    theta[moduleID] = vBeam[moduleID].Angle(vScat[moduleID]);
    theta[moduleID] = RadToDeg()*theta[moduleID];
    
    dR[moduleID] = vScat[moduleID].Mag();
    dXY[moduleID] = TMath::Sqrt(vScat[moduleID].X()*
				vScat[moduleID].X()+
				vScat[moduleID].Y()*
				vScat[moduleID].Y());
    dZ[moduleID]  = vScat[moduleID].Z();
    
    sndHitGood = kTRUE;
    
    if(verbose>2){
      PrintTVector3Elements(" p2    ",p2[moduleID]);
      PrintTVector3Elements(" vBeam ",vBeam[moduleID]);
      PrintTVector3Elements(" vScat ",vScat[moduleID]);
    }
    
  }
  
}

void Hits::SetFile(TString userFileName){
  
  fileName      = userFileName + ".root"; 
  
  dataDirectory = "../Data/";
  
  filePath = dataDirectory + fileName;
  
  TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(filePath);
  if (!f || !f->IsOpen()) {
    f = new TFile(filePath);
  }
  TTree * tree;
  f->GetObject("Hits",tree);
  
  Init(tree);
  
}

TString Hits::GetFilePath(){

  cout << endl;
  cout << " filePath = " << filePath << endl;
  cout << endl;

  return filePath;
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

   //
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

void Hits::PrintEventStats(Bool_t prevEvent){

  Float_t eventID2Print = eventID;

  if(prevEvent)
    eventID2Print = eventID_prevEntry;

  if(verbose > 1 ){
    cout << endl;
    cout << " --------------------------------------- " <<  endl;
    cout << " Event " << eventID2Print << " Stats "    << endl;
    cout << endl;
    cout << " eTotMod1        = " << eTotMod1        << endl;
    cout << " n event entries = " << entriesInEvent  << endl;
    
    if(sndHitGood){
      cout << " phi["<<   moduleID << "]          = " << phi[moduleID]   << endl;
      cout << " theta["<< moduleID << "]        = " << theta[moduleID] << endl;
    }
    
    cout << " --------------------------------------- " <<  endl;
    cout << endl;
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
  EndEvent( kTRUE );
}

void Hits::EndEvent(Bool_t prevEvnt = kFALSE){
  
  if(verbose > 3){
    cout << endl; 
    cout << " Start of EndEvent " << endl;
    cout << " prevEvnt = " << prevEvnt << endl;
    cout << endl; 
  }
  
  // skip if first event
  if(eventID_prevEntry == -1){
    if(verbose > 3){
      cout << endl; 
      cout << " First Entry, returning " << endl;
      cout << endl; 
    }
    return;
  }
  
  entriesInEvent = eventEntry+1;
  
  PrintEventStats(prevEvnt);
  
  eventEntry = -1;
  
  //----------------
  // iterate / accumulate event variables
  
  meanHitsPerEvent += entriesInEvent;
  
  // single process events
  if(entriesInEvent == nPhot)
    fracPhotOnly = fracPhotOnly + 1.;
  else if(entriesInEvent == nRay)
    fracRayOnly  = fracRayOnly + 1;
  else if(entriesInEvent == nComp)
    fracCompOnly = fracCompOnly + 1;
  
  // process frequency per event
  //  (once or more per event)
  if(nPhot > 0)
    fracPhot = fracPhot + 1.0;
  
  if(nComp > 0)
    fracComp = fracComp + 1.0;

  if(nRay > 0)
    fracRay  = fracRay + 1.0;
  
  if(nOther > 0)
    fracOther = fracOther + 1.0;
    
  if(eTotMod1 < 0.51)
    fracE_LT_511 = fracE_LT_511 + 1.0;

  if(goodQET)
    fracQET = fracQET + 1.0;
  
  if(goodPET)
    fracPET = fracPET + 1.0;

  Float_t z1 = p1[moduleID].Z()-50.;
  
  //-------------------
  // Fill Histograms
  
  // First hit within desired thickness
  if(HitDet(1)){
    if     (photFirst)
      hZ1Phot[moduleID]->Fill(z1);
    else if(compFirst)
      hZ1Comp[moduleID]->Fill(z1);
    
  }
    
  // (first and )second hit(s) fulfill QET
  if(sndHitGood){
    hPhi[moduleID]->Fill(phi[moduleID]);
    hTheta[moduleID]->Fill(theta[moduleID]);
    hdR[moduleID]->Fill(dR[moduleID]);
    
    hdZdXY[moduleID]->Fill(z1,
			   dXY[moduleID]);
    
    //if(p1[moduleID].Z() < 60.)
    hThetadXY[moduleID]->Fill(theta[moduleID],dXY[moduleID]);
    hThetadR[moduleID]->Fill(theta[moduleID],dR[moduleID]);
    
  }
  
  if(verbose > 3){
    cout << endl; 
    cout << " End of EndEvent " << endl;
    cout << endl; 
  }
}

void Hits::InitHistos(){
  
  if( verbose > 2 ){
    cout << endl;
    cout << " --------------------------------------- " <<  endl;
    cout << "  InitHistos() " << endl;
    cout << " --------------------------------------- " <<  endl;
    cout << endl;
  }
  
  for( Int_t i = 0 ; i < nModules ; i++){
    
    hName.Form("hPhi_%d",i);
    hTitle.Form("hPhi_%d;#phi (deg);Events",i);
    hPhi[i] = new TH1F(hName,hTitle,
		       nBins,-180,180);
    
    hName.Form("hTheta_%d",i);
    hTitle.Form("hTheta_%d;#theta (deg);Events",i);
    hTheta[i] = new TH1F(hName,hTitle,
			 nBins,0,180);
    
    hName.Form("hdR_%d",i);
    hTitle.Form("hdR_%d; dR (mm) ;Events",i);
    hdR[i] = new TH1F(hName,hTitle,
		      nBins,0,20);
  
    hName.Form("hZ1Phot_%d",i);
    hTitle.Form("hZ1Phot_%d; Z1 (mm) ;Events",i);
    hZ1Phot[i] = new TH1F(hName,hTitle,
			  nBins,0,20);
    
    hName.Form("hZ1Comp_%d",i);
    hTitle.Form("hZ1Comp_%d; Z1 (mm) ;Events",i);
    hZ1Comp[i]= new TH1F(hName,hTitle,
			 nBins,0,20);
    
    hName.Form("hdZdXY_%d",i);
    hTitle.Form("hdZdXY_%d; Z_{1} (mm); dXY (mm)",i);
    hdZdXY[i] = new TH2F(hName,hTitle,
			 nBins,0,40,
			 nBins,0,20);
    
    hName.Form("hThetadXY_%d",i);
    hTitle.Form("hThetadXY_%d;#theta (deg); dXY (mm)",i);
    hThetadXY[i] = new TH2F(hName,hTitle,
			    nBins,0,130,
			    nBins,0,10);
  
    hName.Form("hThetadR_%d",i);
    hTitle.Form("hRTheta_%d; #theta (deg) ; dR (mm)",i);
    hThetadR[i] = new TH2F(hName,hTitle,
			   nBins,0,130,
			   nBins,0,10);
  
    
  }
}

void Hits::InitNewEvent(){
  
  eventID_prevEntry = eventID;
  
  nEvents++;
  
  //reset event by event variables
  entriesInEvent = 0;

  eTotMod1 = 0.0;
  eTotMod2 = 0.0;
  
  nRay   = 0;
  nPhot  = 0;
  nComp  = 0;
  nOther = 0;
  
  photFirst = kFALSE;
  compFirst = kFALSE;
  
  fstHitGood  = kFALSE;
  sndHitGood = kFALSE;
  
  goodPET = kFALSE;
  goodQET = kFALSE;
  
  edepTot511 = kFALSE;

  p0[moduleID].SetXYZ(-999.9,-999.9,-999.9);
  p1[moduleID].SetXYZ(0,0,0);
  p2[moduleID].SetXYZ(0,0,0);

  if( verbose > 1 ) {
    cout << " --------------------------------------- " <<  endl;
    cout << endl;
    cout << " Processing event " << eventID <<  endl;
  }
  
}

void Hits::EndRun(){
  
  FinalStats();
  
  PlotHistograms();  
  
  ResetRunVariables();
  
}

void Hits::ResetRunVariables(){ 

  nHits = 0.;
  
  meanHitsPerEvent  = 0.;
  eventID_prevEntry = -1;
  
  fracPhotOnly = 0;
  fracRayOnly  = 0;
  fracPhot     = 0;
  fracComp     = 0;
  
  nEvents      = 0;
}

void Hits::FinalStats(){

  // Last iteration of loop
  EndEvent();

  // print full data set details
  cout <<  endl;
  cout << " --------------------------------------- "   << endl;
  cout << " --------------------------------------- "   << endl;
  cout << " ------------- Final Stats  ------------ "   << endl;
  cout <<  endl;
  cout << "    File: "        << fileName               << endl;  
  cout <<  endl;
  cout << "    Analysis of single CZT detector      "   << endl;
  cout << "    GATE simulation                      "   << endl;
  cout << " --------------------------------------- "   << endl;
  
  // hits
  nPrimaries = eventID + 1;
  meanHitsPerEvent = (Double_t)meanHitsPerEvent/nEvents;
  
  cout << " --   Hits  "  <<  endl;
  cout <<  endl;
  cout << " Primaries        = " << nPrimaries          <<  endl;
  cout << " Events with hits = " << nEvents             <<  endl;
  cout << " Hits             = " << nHits               <<  endl;
  cout << " Mean hits/event  = " << meanHitsPerEvent    <<  endl;
  cout <<  endl;

  // processes

  //!!!!!!
  //nPrimaries = nEvents;
  
  fracPhotOnly  = fracPhotOnly/(Float_t)nPrimaries;
  fracRayOnly   = fracRayOnly/(Float_t)nPrimaries;
  fracCompOnly  = fracCompOnly/(Float_t)nPrimaries;
  
  fracPhot     = fracPhot/(Float_t)nPrimaries;
  fracRay      = fracRay/(Float_t)nPrimaries;
  fracComp     = fracComp/(Float_t)nPrimaries;
  fracOther    = fracOther/(Float_t)nPrimaries;

  fracPhotFirst  = fracPhotFirst/(Float_t)nPrimaries;
  fracRayFirst   = fracRayFirst/(Float_t)nPrimaries;
  fracCompFirst  = fracCompFirst/(Float_t)nPrimaries;
  fracOtherFirst = fracOtherFirst/(Float_t)nPrimaries;
  
  fracPET      = fracPET/(Float_t)nPrimaries;
  fracQET      = fracQET/(Float_t)nPrimaries;
  fracE_LT_511 = fracE_LT_511/(Float_t)nPrimaries;
  
  Float_t QETfraction = fracQET/(fracPET+fracQET);

  Float_t meanPhotPerPrimary = nPhotTot/nPrimaries;
  Float_t meanCompPerPrimary = nCompTot/nPrimaries;
  Float_t meanRayPerPrimary  = nRayTot/nPrimaries;
  
  cout << fixed << setprecision(3);


  cout << " ---------------------------------------- "   <<  endl;
  cout << " ----- Processes [ / primaries ] ------- "   <<  endl;
  
  if(verbose > 0 ){
    cout <<  endl;
    cout << "  Photoelec only      = " << fracPhotOnly  <<  endl; 
    cout << "  Compton   only      = " << fracCompOnly  <<  endl; 
    cout << "  Rayleigh  only      = " << fracRayOnly   <<  endl; 
    cout <<  endl;
    cout << "  Photoelectric       = " << fracPhot      <<  endl; 
    cout << "  Compton             = " << fracComp      <<  endl; 
    cout << "  Rayleigh            = " << fracRay       <<  endl; 
    cout << "  Other               = " << fracOther     <<  endl; 
    cout <<  endl;
    cout << "  Photoelec first     = " << fracPhotFirst  <<  endl; 
    cout << "  Compton   first     = " << fracCompFirst  <<  endl; 
    cout << "  Rayleigh  first     = " << fracRayFirst   <<  endl; 
    cout << "  Other     first     = " << fracOtherFirst <<  endl; 
    cout <<  endl;
  }
  
  cout << " Photo             511 = " << fracPET      <<  endl; 
  cout << " Compton/s + Photo 511 = " << fracQET      <<  endl; 
  cout << " Total Energy  <   511 = " << fracE_LT_511 <<  endl; 
  cout <<  endl;
  
  Float_t doublePET = fracPET*fracPET;
  Float_t doubleQET = fracQET*fracQET;
  Float_t pETQET    = 2*fracPET*fracQET;
  
  Float_t totalPET  = doublePET + pETQET;
  
  Float_t ratioQET  = doubleQET/(doubleQET+totalPET);
    
  cout << " Efficiency for detector pair            "   << endl;
  cout << " QET (Q) = " << fracQET*fracQET        << endl;
  cout << " PET (P) = " << totalPET               << endl;
  cout << " Q/(P+Q) = " << ratioQET               << endl;
  cout << " --------------------------------------- "   << endl;
}

void Hits::PlotHistograms(){
  
  if(verbose >= 1){
    cout << endl;
    cout << " Plotting Histograms " << endl;
  }
 
  for (Int_t i = 0 ; i < 2 ; i++)
    hPhi[i]->SetMinimum(0);
  
  Int_t width  = 800;
  Int_t height = 800;
  
  if(nModules==2)
    height = 800;

  TCanvas * canvas =  new TCanvas("canvas","canvas",width,height);
  
  //canvas->Divide(2,nModules);
  canvas->Divide(2,2);
  canvas->cd(1);
  hZ1Phot[0]->Draw();
  canvas->cd(2);
  hZ1Comp[0]->Draw();

  /* canvas->cd(1); */
/*   hdR[0]->Draw(); */
/*   canvas->cd(2); */
/*   hdZdXY[0]->Draw("colz"); */
/*   canvas->cd(3); */
/*   hThetadXY[0]->Draw("colz"); */
/*   canvas->cd(4); */
/*   hThetadR[0]->Draw("colz"); */
  
  for (Int_t i = 0 ; i < nModules ; i++){
    //canvas->cd(1+2*i); */
/*     hPhi[i]->Draw(); */
/*     canvas->cd(2+2*i); */
    //hTheta[i]->Draw();
    //hdR[i]->Draw();
    //hdZdXY[i]->Draw("colz");
    //hThetadXY[i]->Draw("colz");
  }
  
/*   Char_t cont = 'n'; */
/*   cout << " continue ? " << endl; */
/*   cin >> cont; */
  
}

void Hits::PrintEntry(Long64_t entryNum){
  
  cout << " -------------------------- "       << endl;
  cout << " primaryID   = " << primaryID       << endl;
  cout << " sourceID    = " << sourceID        << endl;
  cout << " runID       = " << runID           << endl;
  cout << " eventID     = " << eventID         << endl;
  cout << " entry       = " << entryNum        << endl;
  cout << " process     = " << processName     << endl;
  cout << " trackID     = " << trackID         << endl;
  cout << " moduleID    = " << moduleID        << endl;
  cout << " photonID    = " << photonID        << endl;
  cout << " time        = " << time            << endl;
  cout << " parentID    = " << parentID        << endl;
  cout << " edep        = " << edep            << endl;
  cout << " posZ        = " << posZ            << endl;
  cout << " -------------------------- "       << endl;

}

void Hits::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

void Hits::SumEnergyDep(){
  
  if( verbose > 1 ){
    cout << endl;
    cout << "  SumEnergyDep " << endl;
    cout << " -------------------- " << endl;
  }
  
  if     (moduleID==0){
    
    //if(hit is in thickness)
       eTotMod1 += edep;
       
    }
  else if(moduleID==1){
    eTotMod2 += edep;
  }
}

void Hits::CountProcesses(TString processStr){
  
  if( verbose > 1 ){
    cout << endl;
    cout << "  Count Processes " << endl;
    cout << " -------------------- " << endl;
  }
  
  if( eTotMod1 > 0.510 ){
    edepTot511 = kTRUE;
    if(verbose > 2 )
      cout << " edepTot511 = " << edepTot511 << endl;
  }
  
  // Rayleigh Scattering
  if     (processStr=="Rayl"){
    nRay++;
    nRayTot++;
    
    if(eventEntry==0)
      fracRayFirst = fracRayFirst + 1.0;
  }
  // Photoelectric
  else if(processStr=="phot"){
    nPhot++; nPhotTot++;
    
    if     (eventEntry==0){
      fracPhotFirst = fracPhotFirst + 1.0;
      photFirst = kTRUE;
      
      if(edepTot511)
	goodPET = kTRUE;
    }
    else if (edepTot511){
      if(nRay == 0 && nOther==0)
	goodQET = kTRUE;
    } 
  }
  // Compton
  else if(processStr=="compt"){
    nComp++;
    nCompTot++;
    
    if(eventEntry==0){
      fracCompFirst = fracCompFirst + 1.0;
      compFirst = kTRUE;
    }
  }
  // Other processes
  else {
    nOther++;
    nOtherTot++;
    if(eventEntry==1)
      fracOtherFirst = fracOtherFirst + 1.0;
    if(verbose > 2)
      cout << " Process = " << processStr << endl; 
  }
}

void Hits::SetVerbose(Int_t v){
  
  verbose = v;

}

void Hits::SetThickness(Float_t t){
  
  thickness = t;

}


#endif // #ifdef Hits_cxx
