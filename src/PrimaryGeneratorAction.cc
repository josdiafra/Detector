

#include "PrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "DetectorConstruction.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4GeneralParticleSource.hh"
#include "G4SingleParticleSource.hh"
#include "G4RadioactiveDecay.hh"
#include "PrimaryGeneratorMessenger.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction(),
   fParticleGun(0)
{
       
  fParticleGun = new G4GeneralParticleSource();
 
       
// Creamos un messenger para esta clase.
  fGunMessenger = new PrimaryGeneratorMessenger(this);

 
}


PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
  delete fGunMessenger;    
}


void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{  
    fParticleGun->GeneratePrimaryVertex(anEvent);          
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

