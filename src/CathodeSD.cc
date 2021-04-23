//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************

#include "CathodeSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4RunManager.hh"
#include "G4Run.hh"


#include <sys/stat.h>
#include <unistd.h>
#include <string>
#include <fstream>
#include <iostream>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CathodeSD::CathodeSD(const G4String& name, const G4String& hitsCollectionName) 
 : G4VSensitiveDetector(name),
   fHitsCollection(NULL)
{
  collectionName.insert(hitsCollectionName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CathodeSD::~CathodeSD() 
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CathodeSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection
  fHitsCollection = new CathodeHitsCollection(SensitiveDetectorName, collectionName[0]); 

  // Add this collection in hce
  G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fHitsCollection ); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool CathodeSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{  
	G4String processName = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
  	G4double edep = aStep->GetTotalEnergyDeposit();
  	if(edep==0.) return false;

  	CathodeHit* newHit = new CathodeHit();
	newHit->SetPos(aStep->GetPostStepPoint()->GetPosition());
	newHit->SetTime(aStep->GetTrack()->GetGlobalTime());
	newHit->SetEnergyDeposit(edep);
	
	fHitsCollection->insert( newHit );

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CathodeSD::EndOfEvent(G4HCofThisEvent*)
{
    // *eventManager=G4EventManager::GetEventManager ();
	//const G4Event * thisEvent=eventManager->GetConstCurrentEvent();
	
	const G4Run* run = G4RunManager::GetRunManager()->GetCurrentRun();
	
    G4int nofHits = fHitsCollection->entries();
             
    for ( G4int i=0; i<nofHits; i++ ) {
		(*fHitsCollection)[i]->Print();
		(*fHitsCollection)[i]->Draw();
		}

	
		std::ofstream file;
		G4String name="output_"+ std::to_string(run->GetRunID()) +".dat";
          
    	for ( G4int i=0; i<nofHits; i++ ) {
		 	  		 
  			std::ifstream ifile(name.c_str());
		
			if (ifile) {
  				// The file exists, and is open for input
  				ifile.close();
  				file.open(name.c_str(), std::ofstream::out | std::ofstream::app );
				file << (*fHitsCollection)[i]->GetPos().getX()<<", "<<(*fHitsCollection)[i]->GetPos().getY()<<", "<<(*fHitsCollection)[i]->GetPos().getZ() <<", "<<
				(*fHitsCollection)[i]->GetEnergyDeposit() <<", "<< (*fHitsCollection)[i]->GetTime() <<" "<<"\n";

				file.close();
				}
  
			else{
				file.open(name.c_str(), std::ofstream::out ); 
				file << "x [mm], y [mm], z [mm], Edep [MeV], Time [ns]\n";
				file << (*fHitsCollection)[i]->GetPos().getX()<<", "<<(*fHitsCollection)[i]->GetPos().getY()<<", "<<(*fHitsCollection)[i]->GetPos().getZ() <<", "<<
				(*fHitsCollection)[i]->GetEnergyDeposit() <<", "<< (*fHitsCollection)[i]->GetTime()<<" "<<"\n";

    			file.close();
				}
			}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......