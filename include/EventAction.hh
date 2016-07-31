#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

#include "HPGeDigitizer.hh"

/// Event action class
///
/// It defines data members to hold the energy deposit and track lengths
/// of charged particles in Absober layers:
/// - fEnergyAbs, fEnergyGap
/// which are collected step by step via the functions
/// - AddAbs(), AddGap()

class HPGeDigitizer;

class EventAction : public G4UserEventAction
{
public:
  EventAction();
  virtual ~EventAction();

  virtual void  BeginOfEventAction(const G4Event* event);
  virtual void  EndOfEventAction(const G4Event* event);

  void AddAbs(G4double de, G4double dl);
  void AddGap(G4double de, G4double dl);

  G4double  fEnergyAbs;
    G4double  totalE;   //Energia convolucionada hp + e_noise
    G4double  totalE2;  // ''       ''      hp
    G4double  totalE3;  // ''       ''      e_noise
    HPGeDigitizer * digi;
    
  private:

    G4double  fEnergyGap;
    G4double  fTrackLAbs; 
    G4double  fTrackLGap;
    
  };

// inline functions

  inline void EventAction::AddAbs(G4double de, G4double dl) {
    fEnergyAbs += de; 
    fTrackLAbs += dl;
  }

  inline void EventAction::AddGap(G4double de, G4double dl) {
    fEnergyGap += de; 
    fTrackLGap += dl;
  }


#endif


