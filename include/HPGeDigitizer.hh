//HPGe detector simulation (Point) - Santiago Hurtado


#ifndef HPGeDigitizer_h
#define HPGeDigitizer_h 1

#include "G4VDigitizerModule.hh"
#include "HPGeDigi.hh"
#include "globals.hh"
#include "G4Event.hh"

class HPGeDigitizerMessenger;


class HPGeDigitizer : public G4VDigitizerModule
{
public:
  
  HPGeDigitizer(G4String name);
  ~HPGeDigitizer();
  
  void Digitize();
  void SetThreshold(G4double val) { Energy_Threshold = val;}    
  G4double ADC(G4double);
  G4double ADC_hole_pair(G4double);
  G4double ADC_electronic_noise(G4double);    
  
private:
  HPGeDigitsCollection*  DigitsCollection;
  G4double Energy_Threshold;
  HPGeDigitizerMessenger* digiMessenger;
  G4double totalEgauss;

};

#endif







 
