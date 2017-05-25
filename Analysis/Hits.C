#define Hits_cxx
#include "Hits_Single.h"

void Hits::Loop(Int_t userVerbose,
		Long64_t maxEvents)
{
  if (fChain == 0) return;
  
  verbose = userVerbose;
  
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  
  maxEntry = nentries;
  
  // jentry is entry in the chain
   for (Long64_t jentry = 0; jentry < maxEntry ; jentry++) {
     
     // ientry is entry in the tree
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     
      nb = fChain->GetEntry(jentry);   
      nbytes += nb;
      
      //----------------
      // New event?
      if(eventID != eventID_prevEntry){
	
	if(eventID != (eventID_prevEntry+1) &&
	   verbose > 0){
	  cout << endl;
	  cout << " eventID           = " 
	       << eventID << endl;
	  cout << " eventID_prevEntry = " 
	       << eventID_prevEntry << endl;
	}

	EndPreviousEvent();
	InitNewEvent();
      }
      
      // Event Processing
      eventEntry++;
      
      //photonEntry++;
      
      processStr = processName;
      
      // calculate phi
      if     (eventEntry == 1 && 
	      processStr=="compt"){
	
	p0_g1.SetXYZ(sourcePosX,
		  sourcePosY,
		  sourcePosZ);
	
	p1_g1.SetXYZ(posX,posY,posZ);
	
	
	if(verbose > 2){
	  PrintTVector3Elements("p0_g1",p0_g1);
	  PrintTVector3Elements("p1_g1",p1_g1);
	}
	
	fstHitGood = kTRUE;
      }
      else if(eventEntry == 2      && 
	      (processStr=="compt" ||
	       processStr=="phot") &&
	      fstHitGood){
	
	p2_g1.SetXYZ(posX,posY,posZ);
	
	vBeam_g1 = p1_g1 - p0_g1;
	vScat_g1 = p2_g1 - p1_g1;
	
	phi = vScat_g1.DeltaPhi(vBeam_g1);
	phi = RadToDeg()*phi;
	theta = vBeam_g1.Angle(vScat_g1);
	theta = RadToDeg()*theta;
	
	sndHitGood = kTRUE;

	if(verbose>2){
	  PrintTVector3Elements("vBeam_g1",vBeam_g1);
	  PrintTVector3Elements("vScat_g1",vScat_g1);
	}
	
      }
	

      SumEnergyDep();
      
      CountProcesses(processName);
      
      if(verbose > 1)
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
