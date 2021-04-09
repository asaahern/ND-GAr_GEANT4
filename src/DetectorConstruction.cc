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
// *******************************************************************
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4MaterialTable.hh"
#include "G4NistManager.hh"
#include "CathodeSD.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4ThreeVector.hh"

#include "G4PhysicalConstants.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4SystemOfUnits.hh"
#include "Materials.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4UnitsTable.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include <G4NistMaterialBuilder.hh>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
	: G4VUserDetectorConstruction()
{	
	world_half_side = 5*m;
	cylinder_length = 2.5*m;
	base_radius     = 2.5*m;
	steel_thickness = 10.*mm;
	cath_thickness  = 10.*mm;

	disable_gas_pure = true;
	disable_gas_mix = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction(){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
	// Option to switch on/off checking of volumes overlaps
	G4bool checkOverlaps = true;

	G4NistManager* nist = G4NistManager::Instance();

	 // definition of colours
	G4VisAttributes * blue = new G4VisAttributes(G4Colour(0. ,0.8 ,0.9, 0.6));
	blue -> SetVisibility(true);
	blue -> SetForceSolid(true);

	G4VisAttributes * red = new G4VisAttributes(G4Colour(0.3 ,0.1 ,0.));
	red -> SetVisibility(true);
	red -> SetForceSolid(true);

	G4VisAttributes * yellow = new G4VisAttributes(G4Colour(1. ,1. ,0., 1.0));
	yellow -> SetVisibility(true);
	yellow -> SetForceSolid(true);

	G4VisAttributes * lightgray = new G4VisAttributes(G4Colour(0.7 ,0.7 ,0.7));
	lightgray -> SetVisibility(true);
	lightgray -> SetForceSolid(true);
	

// ------------ Materials ------------

	G4Material* galactic = nist->FindOrBuildMaterial("G4_Galactic");

	G4Material* cath_mat = nist->FindOrBuildMaterial("G4_Cu");
	const G4int num_Rcath = 2;
	G4double ephoton_cath[num_Rcath] = {1.0*eV, 10.0*eV};
	G4double efficiency_cath[num_Rcath] = {1.0 , 1.0};
	G4double reflectivity_cath[num_Rcath] = {0.0, 0.0};


	G4Material* steel_mat = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
	const G4int num_Rsteel = 2;
	G4double ephoton_steel[num_Rsteel] = {1.0*eV, 10.0*eV};
	G4double reflectivity_steel[num_Rsteel] = {0., 0.}; 


	// gas
	G4int ncomponents;
	G4double fractionmass;

	G4double pressure = 10.*bar;
	G4double temperature = 298.15*kelvin;
	G4double Ar_dens  = 1.7836; // [mg/cm3]
	G4double CH4_dens = 0.7174; // [mg/cm3]
	G4double Ar_mass  = 39.948;
	G4double CH4_mass = 16.04; 


	G4cout << "TPC is full of pure Ar" << G4endl;
	G4Material* gas_mat = nist->ConstructNewGasMaterial("Ar_pure", "G4_Ar", temperature, pressure);
	// include refractive index of the gas as a function of the photon energy
	const G4int num_gas = 2;
	G4double ephoton_gas[num_gas] = { 0.1*eV, 10*eV};
	G4double refractiveIndex_gas[num_gas] = { 1., 1.};
	G4MaterialPropertiesTable* gas_MPT = new G4MaterialPropertiesTable();
	gas_MPT->AddProperty("RINDEX", ephoton_gas, refractiveIndex_gas, num_gas)->SetSpline(true);
	gas_mat->SetMaterialPropertiesTable(gas_MPT);

/*
	if (!disable_gas_mix) { // Argon 90% / CH4 10%
		G4cout << "TPC is full of Argon/CH4 90/10" << G4endl;
		G4double Ar_frac = 0.9;
		G4double CH4_frac = 0.1;
		G4double gas_dens = (Ar_dens*Ar_frac + CH4_dens*CH4_frac)*mg/cm3;


 		G4Material* Argon = nist->FindOrBuildMaterial("G4_Ar");
		G4Material* Methane= nist->FindOrBuildMaterial("G4_METHANE");

		G4Material* gas_mat = new G4Material("Ar_CH4", gas_dens, ncomponents=2);
		gas_mat->AddMaterial(Argon,   fractionmass=(Ar_mass*Ar_frac/(Ar_mass*Ar_frac+CH4_mass*CH4_frac)));
		gas_mat->AddMaterial(Methane, fractionmass=(CH4_mass*CH4_frac/(Ar_mass*Ar_frac+CH4_mass*CH4_frac)));		}
	}
	else{
		G4Material* gas_mat = galactic;
		G4cout << "TPC is empty!" << G4endl;
	}
*/
  // ------------- Volumes --------------

    // The experimental Hall
	G4Box* expHall_box = new G4Box("World", world_half_side, world_half_side, world_half_side);
	G4LogicalVolume* expHall_log = new G4LogicalVolume(expHall_box, galactic, "World", 0, 0, 0);
	G4VPhysicalVolume* expHall_phys = new G4PVPlacement(0, G4ThreeVector(), expHall_log, "World", 0, false, 0, checkOverlaps);

	// gas volume
	G4Tubs* gas_volume = new G4Tubs("gas_volume", 0., base_radius, cylinder_length, 0.*deg, 360.*deg);
	G4LogicalVolume* gas_volume_log = new G4LogicalVolume(gas_volume, gas_mat, "gas_volume", 0, 0, 0);
	G4VPhysicalVolume* gas_volume_phys = new G4PVPlacement(0, G4ThreeVector(), gas_volume_log, "gas_volume", expHall_log, false, 0, checkOverlaps);
	gas_volume_log -> SetVisAttributes(blue);	



	// Stainless-steel vessel

	// base volume
	G4Tubs* base = new G4Tubs("base", 0.0, base_radius, steel_thickness, 0. * deg, 360. * deg);
	G4LogicalVolume* base_log = new G4LogicalVolume(base, steel_mat, "base", 0, 0, 0);
	G4VPhysicalVolume* base_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, cylinder_length+steel_thickness), base_log, "base", expHall_log, false, 0, checkOverlaps);
	base_log -> SetVisAttributes(lightgray);
	// optical surface
	G4OpticalSurface* op_base = new G4OpticalSurface("base_Surface");
	G4LogicalSkinSurface* base_Surface = new G4LogicalSkinSurface("base_Surface", base_log, op_base);
	G4MaterialPropertiesTable* baseST2 = new G4MaterialPropertiesTable();
	baseST2->AddProperty("REFLECTIVITY", ephoton_steel, reflectivity_steel, num_Rsteel);
	op_base->SetType(dielectric_metal);
	op_base->SetFinish(ground);
	op_base->SetModel(glisur);
	op_base->SetPolish(0.4);
	op_base->SetMaterialPropertiesTable(baseST2);

	// walls volume
	G4Tubs* wall = new G4Tubs("wall", base_radius, base_radius+steel_thickness, cylinder_length, 0. * deg, 360. * deg);
	G4LogicalVolume* wall_log	= new G4LogicalVolume(wall, steel_mat, "wall", 0, 0, 0);
	G4VPhysicalVolume* wall_phys = new G4PVPlacement(0, G4ThreeVector(), wall_log, "wall", expHall_log, false, 0, checkOverlaps);
	wall_log  -> SetVisAttributes(lightgray);
	// optical surface
	G4OpticalSurface* op_wall = new G4OpticalSurface("wall_Surface");
	G4LogicalSkinSurface* wall_Surface = new G4LogicalSkinSurface("wall_Surface", wall_log, op_wall);
	G4MaterialPropertiesTable* wallST2 = new G4MaterialPropertiesTable();
	wallST2->AddProperty("REFLECTIVITY",  ephoton_steel, reflectivity_steel, num_Rsteel);
	op_wall->SetMaterialPropertiesTable(wallST2);
	op_wall->SetType(dielectric_metal);
	op_wall->SetFinish(ground);
	op_wall->SetPolish(0.2);
	op_wall->SetModel(glisur);

	
	// Cathode
	G4Tubs* Cath = new G4Tubs("Cath", 0., base_radius, cath_thickness, 0.*deg, 360.*deg);
	G4LogicalVolume* Cath_log = new G4LogicalVolume(Cath, cath_mat, "Cath", 0, 0, 0);
	G4VPhysicalVolume* Cath_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, -cylinder_length - cath_thickness), Cath_log, "Cath", expHall_log, false, 0, checkOverlaps);
	Cath_log -> SetVisAttributes(red);
	// optical surface
	G4OpticalSurface* opPhCa = new G4OpticalSurface("Cath_Surface");
	G4LogicalSkinSurface* Cath_Surface = new G4LogicalSkinSurface("PhCa_SkinSurface", Cath_log, opPhCa);
	G4MaterialPropertiesTable* PhCaST2 = new G4MaterialPropertiesTable();
	PhCaST2->AddProperty("EFFICIENCY", ephoton_cath, efficiency_cath, num_Rcath);
	PhCaST2->AddProperty("REFLECTIVITY", ephoton_cath, reflectivity_cath, num_Rcath);
	opPhCa->SetMaterialPropertiesTable(PhCaST2);
	opPhCa->SetType(dielectric_metal);
	opPhCa->SetFinish(polished);
	opPhCa->SetModel(glisur);

	
	//always return the physical World
	return expHall_phys;

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
  // ------------- Sensitive detectors ------------- 
  CathodeSD* CaSD = new CathodeSD("/CathodeSD","CathodeHitsCollection");  
  G4SDManager::GetSDMpointer()->AddNewDetector(CaSD);
  SetSensitiveDetector("Cath", CaSD, true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......