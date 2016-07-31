
#include "EventAction.hh"
#include "RunAction.hh"
#include "Analysis.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>

#include "HPGeDigitizer.hh"


EventAction::EventAction()
 : G4UserEventAction(),
   fEnergyAbs(0.),
   fEnergyGap(0.),
   fTrackLAbs(0.),
   fTrackLGap(0.)
{}


EventAction::~EventAction()
{}


void EventAction::BeginOfEventAction(const G4Event* /*event*/)
{  
  // initialisation per event
  fEnergyAbs = 0.;
  fEnergyGap = 0.;
  fTrackLAbs = 0.;
  fTrackLGap = 0.;
}

void EventAction::EndOfEventAction(const G4Event* event)
{
    
    HPGeDigitizer * digi = new HPGeDigitizer("");
    
    totalE = digi->ADC(fEnergyAbs);
    
    totalE2 = digi->ADC_hole_pair(fEnergyAbs);
    
    totalE3 = digi->ADC_electronic_noise(fEnergyAbs);
    
  // Accumulate statistics
  //

  // get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
   
  // fill histograms
  analysisManager->FillH1(1, fEnergyAbs);
  analysisManager->FillH1(2, totalE);
  analysisManager->FillH1(3, totalE2);
  analysisManager->FillH1(4, totalE3);    
  
}  



