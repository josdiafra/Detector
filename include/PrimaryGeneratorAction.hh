
#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4GeneralParticleSource.hh"
#include "G4SingleParticleSource.hh"
#include "G4RadioactiveDecay.hh"


class G4ParticleGun;
class G4Event;

/// The primary generator action class with particle gum.
///
/// It defines a single particle which hits the calorimeter 
/// perpendicular to the input face. The type of the particle
/// can be changed via the G4 build-in commands of G4ParticleGun class 
/// (see the macros provided with this example).
class PrimaryGeneratorMessenger;
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{

public:
  
  PrimaryGeneratorAction();    
  virtual ~PrimaryGeneratorAction();

  virtual void GeneratePrimaries(G4Event* event);
  
  
private:
  
  
  G4GeneralParticleSource*  fParticleGun; // G4 particle source
  G4GeneralParticleSource*  distInhomogenea;
  
   PrimaryGeneratorMessenger* fGunMessenger; //messenger of this class
   
 };


#endif


