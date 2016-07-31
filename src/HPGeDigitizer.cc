//HPGe detector simulation Á(Point) - Santiago Hurtado

#include "HPGeDigitizer.hh"
#include "HPGeDigi.hh"
#include "HPGeDigitizerMessenger.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"


//#include "HPGeDetectorHit.hh"
#include "EventAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4DigiManager.hh"
#include "G4ios.hh"

#include "G4UnitsTable.hh"
#include "Randomize.hh"

#include <vector>


HPGeDigitizer::HPGeDigitizer(G4String name)
  :G4VDigitizerModule(name)
{
  G4String colName = "DigitsCollection";
  collectionName.push_back(colName);
  Energy_Threshold = 1.*keV;
  totalEgauss = 0.*keV;
      
//create a messenger for this class
  digiMessenger = new HPGeDigitizerMessenger(this);

}


HPGeDigitizer::~HPGeDigitizer()
{
  delete digiMessenger;
}

G4double HPGeDigitizer::ADC(G4double totalenergy)
{
  if (totalenergy>0.) {
    //G4cout<<"totalenergy="<<totalenergy/keV<<G4endl;
    //G4double fwhm = pow(2.29+(0.00233*totalenergy/keV)+(0.000000453*totalenergy/keV*totalenergy/keV),0.5);
    //G4cout<<"fwhm="<<fwhm<<G4endl;
    //G4double sigma = fwhm/2.355;
      G4double sigma = 0.00749*MeV;
    //G4cout<<"sigma="<<sigma<<G4endl;
    do {
      totalEgauss = G4RandGauss::shoot(totalenergy*MeV,sigma);
      //G4cout<<"totalEgauss="<<totalEgauss<<G4endl;
      //G4cout<<"L="<<totalenergy/keV-3*sigma<<"  R="<<totalenergy/keV+3*sigma<<G4endl;
    } while ((totalEgauss>(totalenergy*MeV+2.96*sigma))||(totalEgauss<(totalenergy*MeV-2.96*sigma)));
  }	  
  return totalEgauss;
}

G4double HPGeDigitizer::ADC_hole_pair(G4double totalenergy)
{
  if (totalenergy>0.) {
    //G4cout<<"totalenergy="<<totalenergy/keV<<G4endl;
    //G4double fwhm = pow(2.29+(0.00233*totalenergy/keV)+(0.000000453*totalenergy/keV*totalenergy/keV),0.5);
    //G4cout<<"fwhm="<<fwhm<<G4endl;
    //G4double sigma = fwhm/2.355;
      G4double sigma = 0.00860*MeV;
    //G4cout<<"sigma="<<sigma<<G4endl;
    do {
      totalEgauss = G4RandGauss::shoot(totalenergy*MeV,sigma);
      //G4cout<<"totalEgauss="<<totalEgauss<<G4endl;
      //G4cout<<"L="<<totalenergy/keV-3*sigma<<"  R="<<totalenergy/keV+3*sigma<<G4endl;
    } while ((totalEgauss>(totalenergy*MeV+2.96*sigma))||(totalEgauss<(totalenergy*MeV-2.96*sigma)));
  }	  
  return totalEgauss;
}

G4double HPGeDigitizer::ADC_electronic_noise(G4double totalenergy)
{
  if (totalenergy>0.) {
    //G4cout<<"totalenergy="<<totalenergy/keV<<G4endl;
    //G4double fwhm = pow(2.29+(0.00233*totalenergy/keV)+(0.000000453*totalenergy/keV*totalenergy/keV),0.5);
    //G4cout<<"fwhm="<<fwhm<<G4endl;
    //G4double sigma = fwhm/2.355;
      G4double sigma = 0.00781*MeV;
    //G4cout<<"sigma="<<sigma<<G4endl;
    do {
      totalEgauss = G4RandGauss::shoot(totalenergy*MeV,sigma);
      //G4cout<<"totalEgauss="<<totalEgauss<<G4endl;
      //G4cout<<"L="<<totalenergy/keV-3*sigma<<"  R="<<totalenergy/keV+3*sigma<<G4endl;
    } while ((totalEgauss>(totalenergy*MeV+2.96*sigma))||(totalEgauss<(totalenergy*MeV-2.96*sigma)));
  }	  
  return totalEgauss;
}











 

