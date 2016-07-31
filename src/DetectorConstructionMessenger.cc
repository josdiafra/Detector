
#include "DetectorConstructionMessenger.hh"

#include "DetectorConstruction.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4SystemOfUnits.hh"

DetectorConstructionMessenger::DetectorConstructionMessenger(
                                           DetectorConstruction* action)
 : G4UImessenger(),
   fAction(action),
   setEspesorSi(0)
{
  setEspesorSi = new G4UIcmdWithADoubleAndUnit("/detector/capaSi/setEspesorSi",this);
  setEspesorSi->SetGuidance("Set thickness Si layer");
  setEspesorSi->SetParameterName("length",false);
  setEspesorSi->SetRange("length>0.");
  setEspesorSi->SetUnitCategory("Length");
  setEspesorSi->AvailableForStates(G4State_PreInit,G4State_Idle);
  setEspesorSi->SetToBeBroadcasted(false);
       
       
  setPosicionMuestra = new G4UIcmdWithADoubleAndUnit("/detector/posicionMuestra/setPosicionMuestra",this);
  setPosicionMuestra->SetGuidance("Set position source");
  setPosicionMuestra->SetParameterName("length",false);
  setPosicionMuestra->SetRange("length>0.");
  setPosicionMuestra->SetUnitCategory("Length");
  setPosicionMuestra->AvailableForStates(G4State_PreInit,G4State_Idle);
  setPosicionMuestra->SetToBeBroadcasted(false);
       
  setPosicion = new G4UIcmdWithAnInteger("/detector/posicionMuestra/setPosicionBandeja",this);
  setPosicion->SetGuidance("Set position source");
       
  setPresion = new G4UIcmdWithADoubleAndUnit("/detector/presionCamara/setPresionCamara",this);
  setPresion->SetGuidance("Set pressure chamber");
  setPresion->SetParameterName("pressure",false);
  setPresion->SetRange("pressure>0.");
  setPresion->SetUnitCategory("Pressure");
  setPresion->AvailableForStates(G4State_PreInit,G4State_Idle);
  setPresion->SetToBeBroadcasted(false);       
  

}


DetectorConstructionMessenger::~DetectorConstructionMessenger()
{
  delete setEspesorSi;
  delete setPosicionMuestra;
  delete setPosicion;
  delete setPresion;    
  
}

void DetectorConstructionMessenger::SetNewValue(G4UIcommand * command,
                                               G4String newValue)
{
  if( command == setEspesorSi )
   { fAction->setEspesorCapaSi(setEspesorSi->GetNewDoubleValue(newValue));}
    
  if( command == setPosicionMuestra )
   { fAction->setPosicionMuestraRadioactiva(setPosicionMuestra->GetNewDoubleValue(newValue));}   
    
  if( command == setPosicion )
   { fAction->setPosicionBandeja(setPosicion->GetNewIntValue(newValue));}   
    
  if( command == setPresion )
   { fAction->setPresionCamara(setPresion->GetNewDoubleValue(newValue));}   
}

