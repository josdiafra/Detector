
#ifndef StackingAction_h
#define StackingAction_h 1

#include "G4UserStackingAction.hh"
#include "globals.hh"
#include "G4VProcess.hh"
#include <set>

class G4ParticleDefinition;
class G4Track;
class StackingMessenger;
class G4VProcess;


class StackingAction : public G4UserStackingAction
{

public:

  StackingAction();
  ~StackingAction();

  G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track*);
  void NewStage();
  void PrepareNewEvent();

  // Accessors
  void AddTrackedSecondary(G4String pName);
  void ClearTrackedSecondaries();
  void TrackAll(G4bool choice);
  
  
  G4int alphaID;

private:  // Data members

  // Flag to tell whether all particle types are tracked.
  G4bool theAllTracked;

  // Flag to tell whether GenericIon is considered for tracking.
  // Here, "GenericIon" applies either for Ions and Anti-ions.
  G4bool theIonsTracked;

  // Vector storing which type of secondary particle is tracked
  std::set<G4ParticleDefinition*>* theTrackedSecondaries;

  // Messenger class
  StackingMessenger* theMessenger;
  
    
};

#endif
