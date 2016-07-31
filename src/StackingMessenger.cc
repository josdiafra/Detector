//MACG 2014-08-06: Created by M.A. Cortes-Giraldo (Univ. Sevilla, Spain)


#include "globals.hh"

#include "StackingMessenger.hh"
#include "StackingAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"


//-------1---------2---------3---------4---------5---------6---------7---------8
// ----------------------------------------------------------------------------

StackingMessenger::StackingMessenger( StackingAction* myStack )
  : theStacking( myStack )
{ 
  theStackingDir = new G4UIdirectory("/stacking/");
  theStackingDir->SetGuidance("Stacking action control.");

  theTrackedSecParticleCmd = new G4UIcmdWithAString("/stacking/trackSecondary",
						    this);
  theTrackedSecParticleCmd
    ->SetGuidance("Add kind of secondary particle accepted in fUrgent stack.");
  theTrackedSecParticleCmd
    ->SetGuidance("Each call is able to issue ONE particle type, but");
  theTrackedSecParticleCmd
    ->SetGuidance(" \"GenericIon\" encompasses all charged (anti-)ions.");
  theTrackedSecParticleCmd->SetParameterName("particleName", false);
  theTrackedSecParticleCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  theClearSecondariesCmd =
    new G4UIcmdWithoutParameter("/stacking/clearSecondaries", this);
  theClearSecondariesCmd
    ->SetGuidance("This command removes ANY kind of secondary particles from ");
  theClearSecondariesCmd
    ->SetGuidance("acceptance into stack (i.e., they are not tracked.)");
  theClearSecondariesCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  theTrackAllCmd = new G4UIcmdWithABool("/stacking/trackAll", this);
  theTrackAllCmd->SetGuidance(" - 'true' = All particles are tracked.");
  theTrackAllCmd
    ->SetGuidance(" - 'false' = Only secondaries registered by means of");
  theTrackAllCmd
    ->SetGuidance("             \"/stacking/trackSecondary\" UI command");
  theTrackAllCmd
    ->SetGuidance("             are put into the stack for tracking.");
  theTrackAllCmd
    ->SetGuidance("NOTE: If 'true', other stacking commands are NOT applied.");
  theTrackAllCmd->SetParameterName("choice", false);
  theTrackAllCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
      
  
      

}


// ----------------------------------------------------------------------------

StackingMessenger::~StackingMessenger()
{
  delete theStackingDir;
  delete theTrackedSecParticleCmd;
  delete theClearSecondariesCmd;
  delete theTrackAllCmd;
}

// ----------------------------------------------------------------------------

void StackingMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == theTrackedSecParticleCmd) {
    theStacking->AddTrackedSecondary(newValue);
  }
  
  else if (command == theClearSecondariesCmd) {
    theStacking->ClearTrackedSecondaries();
  }

  else if (command == theTrackAllCmd) {
    theStacking->TrackAll(theTrackAllCmd->GetNewBoolValue(newValue));
  }
}
