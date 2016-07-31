//HPGe detector simulation (Point) - Santiago Hurtado


#ifndef HPGeDigitizerMessenger_h
#define HPGeDigitizerMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class HPGeDigitizer;
class G4UIcmdWithADoubleAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class HPGeDigitizerMessenger: public G4UImessenger
{
public:
  HPGeDigitizerMessenger(HPGeDigitizer*);
  ~HPGeDigitizerMessenger();
  
  void SetNewValue(G4UIcommand*, G4String);
  
private:
  HPGeDigitizer* HPGeAction; 
  G4UIcmdWithADoubleAndUnit*  ThresholdCmd;
};

#endif

 
