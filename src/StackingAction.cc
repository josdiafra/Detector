
#include "StackingAction.hh"

#include "globals.hh"
#include "G4StackManager.hh"
#include "G4Track.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4Electron.hh"
#include "G4GenericIon.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4VProcess.hh"
#include "StackingMessenger.hh"

#include <vector>


//-----------------------------------------------------------------------------
StackingAction::StackingAction()
: G4UserStackingAction()
{
  alphaID = 0;
  theAllTracked = true;
  theIonsTracked = false;
  theTrackedSecondaries = new std::set<G4ParticleDefinition*>;
  theMessenger = new StackingMessenger(this);
  
}


//-----------------------------------------------------------------------------
StackingAction::~StackingAction()
{
  delete theTrackedSecondaries;
  delete theMessenger;
}


//-----------------------------------------------------------------------------
G4ClassificationOfNewTrack 
StackingAction::ClassifyNewTrack(const G4Track* aTrack)
{
  G4ClassificationOfNewTrack classification = fUrgent;

  G4ParticleDefinition* pDef = aTrack->GetDefinition();
  
  G4String CPName;

  return classification;
}


//-----------------------------------------------------------------------------
// Called by G4StackManager when the urgentStack becomes empty and contents 
// of the waitingStack are transfered to the urgentStack: 
void StackingAction::NewStage()
{;}


//-----------------------------------------------------------------------------
void StackingAction::PrepareNewEvent()
{;}


//-----------------------------------------------------------------------------
void StackingAction::AddTrackedSecondary(G4String pName)
{
  G4ParticleDefinition* pDef = G4ParticleTable::GetParticleTable()
  ->FindParticle(pName);

  theTrackedSecondaries->insert(pDef);
  G4cout << "Secondary \"" << pDef->GetParticleName() 
    //	 << "\" of type \"" << pDef->GetParticleType()
  << "\" are now accepted into 'fUrgent' stack."
  << G4endl;

  if (pName == "GenericIon") {
    theIonsTracked = true;    // Reset flag value
    G4cout << "NOTE: \"" << pName << "\" already includes ALL charged "
    << "(anti-)ions." 
    << G4endl;
  }
}


//-----------------------------------------------------------------------------
void StackingAction::ClearTrackedSecondaries()
{
  theTrackedSecondaries->clear();
  theIonsTracked = false;
  G4cout << "All kind of secondary particles have been dropped for tracking, "
  << "i.e., no secondary particles are being stored in any stack "
  << "at this point."
  << G4endl;
}


//-----------------------------------------------------------------------------
void StackingAction::TrackAll(G4bool choice)
{
  theAllTracked = choice;
  if (theAllTracked)
    G4cout << "ALL particles are now accepted (for tracking) "
  << "into 'fUrgent' stack."
  << G4endl;
  else
    G4cout << "Only selected secondary particles are accepted for tracking."
  << G4endl;
}

//-----------------------------------------------------------------------------



