//HPGe detector simulation (Point) - Santiago Hurtado


#include "HPGeDigi.hh"

G4Allocator<HPGeDigi> HPGeDigiAllocator;


HPGeDigi::HPGeDigi()
{
  edep = 0.;
  weight=0.;
}

HPGeDigi::~HPGeDigi()
{
}

HPGeDigi::HPGeDigi(const HPGeDigi& right)
  :G4VDigi()
{
 edep      = right.edep;
 weight = right.weight;
}

const HPGeDigi& HPGeDigi::operator=(const HPGeDigi& right)
{
  edep      = right.edep;
  weight = right.weight;
  return *this;
}

int HPGeDigi::operator==(const HPGeDigi& right) const
{ 
  return((edep==right.edep)&&(weight==right.weight));
}

void HPGeDigi::Draw()
{;}

void HPGeDigi::Print()
{;}










 
