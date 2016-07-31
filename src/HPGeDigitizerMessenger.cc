//HPGe detector simulation (Point) - Santiago Hurtado

#include "HPGeDigitizerMessenger.hh"

#include "HPGeDigitizer.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

HPGeDigitizerMessenger::HPGeDigitizerMessenger
(HPGeDigitizer* HPGeDigit)
  :HPGeAction(HPGeDigit)
{ 
  ThresholdCmd = new G4UIcmdWithADoubleAndUnit("/digitizer/Threshold",this);
  ThresholdCmd->SetGuidance("Energy deposition threshold for digi generation");
  ThresholdCmd->SetParameterName("choice",true);
  ThresholdCmd->SetDefaultValue((G4double)20.*keV);
  ThresholdCmd->SetRange("Threshold >=0.");
  ThresholdCmd->SetUnitCategory("Energy");  
  ThresholdCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

HPGeDigitizerMessenger::~HPGeDigitizerMessenger()
{
  delete ThresholdCmd;
}

void HPGeDigitizerMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{ 
  if( command == ThresholdCmd )
    { 
      HPGeAction->SetThreshold
	(ThresholdCmd->GetNewDoubleValue(newValue));
    }
}











 
