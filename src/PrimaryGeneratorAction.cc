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
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorAction.hh"

#include "Randomize.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"

#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction(),fParticleGun(0)
{
  G4int n_particle = 1;
  fParticleGun = new G4ParticleGun(n_particle);

  G4double wavelength = 600 * nm;   // Ar scintillation (VUV)
  G4double energy = h_Planck * c_light / wavelength;
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle("opticalphoton");

  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleEnergy(energy);
  fParticleGun->SetParticlePolarization(G4ThreeVector(1., 0., 0.));
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  fParticleGun->GeneratePrimaryVertex(anEvent);

  // define initial position
  G4double cylinder_length = 2.5*m;
  G4double base_radius     = 2.5*m;
  // to uniformly fill the gas volume (cylindrical) 
  G4double theta = G4UniformRand() * 360.0*deg;
  G4double rho   = std::sqrt(G4UniformRand() * base_radius * base_radius);
  G4double posX  = rho * std::cos(theta);
  G4double posY  = rho * std::sin(theta);
  G4double posZ  = (1-2* G4UniformRand()) * cylinder_length;
  // source point (horizontal trace crossing the gas volume):
  //G4double posX  = 0; //(1-2* G4UniformRand()) * 2291;
  //G4double posY  = 0; //1000;
  //G4double posZ  = 0; 
  fParticleGun->SetParticlePosition(G4ThreeVector(posX, posY, posZ));
  //G4cout << posX << ", " << posY << ", " << posZ << G4endl;

  // define momentum direction: isotropic and random
  theta = std::acos(1-2* G4UniformRand())*rad;    // x-z plane (0-2pi)
  G4double phi   = G4UniformRand() * 360.0*deg;   // x-y plane
  G4double mx    = std::sin(theta) * std::cos(phi);
  G4double my    = std::sin(theta) * std::sin(phi);
  G4double mz    = std::cos(theta);
  G4ThreeVector momentumDirection = G4ThreeVector(mx, my, mz);
  fParticleGun->SetParticleMomentumDirection(momentumDirection);
  //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0,0,1)); // launch all photons directed towards base

  // define polarization direction: random but contained in the plane perpendicular to the momentum direction of each photon
  G4ThreeVector normal(1., 0., 0.);
  G4ThreeVector product = normal.cross(momentumDirection);
  G4double modul2 = product * product;

  G4ThreeVector e_perpend(0., 0., 1.);
  if (modul2 > 0.) e_perpend = (1. / std::sqrt(modul2)) * product;
  G4ThreeVector e_paralle = e_perpend.cross(momentumDirection);

  G4double angle = G4UniformRand() * 360.0 * deg;
  G4ThreeVector polar = std::cos(angle) * e_paralle + std::sin(angle) * e_perpend;
  fParticleGun->SetParticlePolarization(polar);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
