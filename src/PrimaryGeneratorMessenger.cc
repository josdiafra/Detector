
#include "PrimaryGeneratorMessenger.hh"

#include "PrimaryGeneratorAction.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4SystemOfUnits.hh"

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(
                                           PrimaryGeneratorAction* action)
 : G4UImessenger(),
   fAction(action),
   setInhomogenea(0)
{
  setInhomogenea = new G4UIcmdWithAnInteger("/gps/source/setInhomogenea",this);
  setInhomogenea->SetGuidance("Seleccionar una distribucion inhomogenea de la fuente");
  setInhomogenea->SetParameterName("setInhomogenea",false,false);
  setInhomogenea->SetRange("setInhomogenea >= 0");

}

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
  delete setInhomogenea;
  
}

void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand * command,
                                               G4String newValue)
{
  if( command == setInhomogenea )
   { fAction->setOpcionDist(setInhomogenea->GetNewIntValue(newValue));}
}
