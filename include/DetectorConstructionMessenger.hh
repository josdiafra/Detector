
#ifndef DetectorConstructionMessenger_h
#define DetectorConstructionMessenger_h 1

#include "G4UImessenger.hh"
#include "G4UIcmdWithAnInteger.hh"


class DetectorConstruction;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger;


class DetectorConstructionMessenger:public G4UImessenger
{
public:
    DetectorConstructionMessenger(DetectorConstruction*);
    virtual ~DetectorConstructionMessenger();
    
    virtual void SetNewValue(G4UIcommand*, G4String);

private:
    DetectorConstruction* fAction;
    G4UIcmdWithADoubleAndUnit* setEspesorSi;
    G4UIcmdWithADoubleAndUnit* setPosicionMuestra;
    G4UIcmdWithAnInteger* setPosicion;
    
    G4UIcmdWithADoubleAndUnit* setPresion;
    
    
};

#endif
