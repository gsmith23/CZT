#define Hits_cxx
#include "Hits.h"

void Hits::Run(Float_t userThickness = 100.0,
	       Int_t userVerbose = 0,
	       Long64_t maxEvents = -1){
  
  SetThickness(userThickness);
  SetVerbose(userVerbose);
  Loop(maxEvents);
  
}

void Hits::Loop(Long64_t maxEvents)
{
  if (fChain == 0) return;

  if(verbose>0){
    cout << "  " << endl;
    cout << " --------------------------------------- " <<  endl;
    cout << " --------------------------------------- " <<  endl;
    cout << "          Looping over entries           " <<  endl;
    cout << " --------------------------------------- " <<  endl;
    cout << " --------------------------------------- " <<  endl;
    cout << "  " << endl;
  }
  
  Long64_t nentries = fChain->GetEntriesFast();
  
  Long64_t nbytes = 0, nb = 0;
  
  InitHistos();
  
  // jentry is entry in the chain
  for (Long64_t jentry = 0 ;
       jentry < nentries   ;
       jentry++) {
    
    // ientry is entry in the tree
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    
    nb = fChain->GetEntry(jentry);   
    nbytes += nb;
    
    if(verbose > 3){
      cout << endl; 
      cout << " eventID = " << eventID << endl; 
      cout << endl; 
    }
    
    //----------------
    // Check if the entry is for new event
    if(eventID != eventID_prevEntry){
      EndPreviousEvent();
      InitNewEvent();
    }
    
    eventEntry++;
    nHits++;
    
    // TStrings are more convenient
    processStr = processName;
    
    // phi, theta, hit position vectors
    SetQETVariables();

    SumEnergyDep();
    
    // Analyse only events
    // with two hits within
    // assigned depth value
    // (variable is thickness)
    if( HitDet(1) &&
	HitDet(2)    ) {
      CountProcesses(processName);
    
    }
    if(verbose > 1)
      PrintEntry(jentry);
    
    if(photonID > 1 && parentID == 0)
      PrintEntry(jentry);
    
    //------------------
    // End of Looping
    
    // truncate loop if user requested
    if(maxEvents >= 0 &&
       nEvents   > (maxEvents-1) ) {
      EndRun();
      break;
    }
    
    // last event     
    if( jentry == (nentries-1))
      EndRun();
    
  }
}
