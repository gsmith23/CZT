#define Hits_cxx
#include "Hits.h"

void Hits::Loop(Int_t userVerbose,
		Long64_t maxEvents)
{
  if (fChain == 0) return;
  
  SetVerbose(userVerbose);
  
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  
  maxEntry = nentries;
  
  InitHistos();
  
  // jentry is entry in the chain
  for (Long64_t jentry = 0; jentry < maxEntry ; jentry++) {
    
    // ientry is entry in the tree
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    
    nb = fChain->GetEntry(jentry);   
    nbytes += nb;
    
    //----------------
    // Check if the entry is for new event
    if(eventID != eventID_prevEntry){
      EndPreviousEvent();
      InitNewEvent();
    }
    
    eventEntry++;
    
    processStr = processName;
    
    // calculate phi
    // first hit of photon
    if     (eventEntry == 0 &&
	    processStr=="compt" ){
      
      p0[moduleID].SetXYZ(sourcePosX,
			  sourcePosY,
			  sourcePosZ);
      
      p1[moduleID].SetXYZ(posX,posY,posZ);
      
      
      if(verbose > 2){
	PrintTVector3Elements(" p0 = ",p0[moduleID]);
	PrintTVector3Elements(" p1 = ",p1[moduleID]);
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
      
      sndHitGood = kTRUE;
      
      if(verbose>2){
	PrintTVector3Elements(" p2    = ",p2[moduleID]);
	PrintTVector3Elements(" vBeam = ",vBeam[moduleID]);
	PrintTVector3Elements(" vScat = ",vScat[moduleID]);
      }
      
    }
    
    SumEnergyDep();
    
    CountProcesses(processName);
    
    if(verbose > 1)
      PrintEntry(jentry);
    
    if(photonID > 1 && parentID == 0)
      PrintEntry(jentry);
    
    //------------------
    // End of Looping
    
    // truncate loop
    if(maxEvents >= 0 &&
       nEvents   > maxEvents) {
      FinalStats();
      break;
    }
    
    // last event     
    if( jentry == (nentries-1))
      FinalStats();
    
  }
}
