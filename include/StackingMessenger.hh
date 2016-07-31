//MACG 2014-08-11: Created by M.A. Cortes-Giraldo (Univ. Sevilla, Spain)

#ifndef StackingMessenger_h
#define StackingMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class StackingAction;
class G4UIdirectory;
class G4UIcmdWithABool;
class G4UIcmdWithAString;
class G4UIcmdWithoutParameter;

// ----------------------------------------------------------------------------

class StackingMessenger: public G4UImessenger
{

public:

  StackingMessenger( StackingAction* );
  ~StackingMessenger();

  void SetNewValue( G4UIcommand*, G4String );

private:

  StackingAction*      theStacking;

  G4UIdirectory*             theStackingDir;
  G4UIcmdWithAString*        theTrackedSecParticleCmd;
  G4UIcmdWithoutParameter*   theClearSecondariesCmd;
  G4UIcmdWithABool*   theTrackAllCmd;

  // MTP - G4UIcmdWithAString is a concrete class of G4UIcommand. 
  // The command defined by this class takes a string. 
  // Incase the parameter string contains space(s), it must be enclosed 
  // by double-quotations (").
};

// ----------------------------------------------------------------------------

#endif
