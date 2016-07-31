//HPGe detector simulation - Santiago Hurtado

#ifndef HPGeDigi_h
#define HPGeDigi_h 1

#include "G4VDigi.hh"
#include "G4TDigiCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class HPGeDigi : public G4VDigi
{

public:

  HPGeDigi();
  ~HPGeDigi();
  HPGeDigi(const HPGeDigi&);
  const HPGeDigi& operator=(const HPGeDigi&);
  int operator==(const HPGeDigi&) const;
  
  inline void* operator new(size_t);
  inline void  operator delete(void*);
  
  void Draw();
  void Print();

private:

  G4double edep;
  G4double weight;

public:

  void SetEdep     (G4double de)      { edep = de; };
  void SetTrackWeight (G4double w) {weight = w;};
  
  G4double GetEdep()    { return edep; };   
  G4double GetTrackWeight() {return weight;};

  
  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

typedef G4TDigiCollection<HPGeDigi> HPGeDigitsCollection;

extern G4Allocator<HPGeDigi> HPGeDigiAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void* HPGeDigi::operator new(size_t)
{
  void* aDigi;
  aDigi = (void*) HPGeDigiAllocator.MallocSingle();
  return aDigi;
}


inline void HPGeDigi::operator delete(void* aDigi)
{
  HPGeDigiAllocator.FreeSingle((HPGeDigi*) aDigi);
}

#endif

