
#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"


SteppingAction::SteppingAction(
                      const DetectorConstruction* detectorConstruction,
                      EventAction* eventAction)
  : G4UserSteppingAction(),
    fDetConstruction(detectorConstruction),
    fEventAction(eventAction)
{
}

SteppingAction::~SteppingAction()
{ 
}

void SteppingAction::UserSteppingAction(const G4Step* step)
{
// Collect energy and track length step by step

  // get volume of the current step
  G4VPhysicalVolume* volume 
    = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  
  // energy deposit
  G4double edep = step->GetTotalEnergyDeposit();
  
  // step length
  G4double stepLength = 0.;
    
    
  if ( step->GetTrack()->GetDefinition()->GetPDGCharge() != 0. ) {
    stepLength = step->GetStepLength();
  }
      
  if ( volume == fDetConstruction->GetAbsorberPV() ) {
    fEventAction->AddAbs(edep,stepLength);
  }
    
 
}

