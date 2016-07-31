

#ifndef PrimaryGeneratorMessenger_h
#define PrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "G4UIcmdWithAnInteger.hh"


class PrimaryGeneratorAction;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger;


class PrimaryGeneratorMessenger: public G4UImessenger
{
  public:
    PrimaryGeneratorMessenger(PrimaryGeneratorAction*);
    virtual ~PrimaryGeneratorMessenger();
    
    virtual void SetNewValue(G4UIcommand*, G4String);

  private:
    PrimaryGeneratorAction* fAction;
    G4UIcmdWithAnInteger* setInhomogenea;
};

#endif
