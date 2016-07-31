
#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "DetectorConstruction.hh"
#include "globals.hh"


class G4VPhysicalVolume;
class DetectorConstructionMessenger;
class G4GlobalMagFieldMessenger;

/// Detector construction class to define materials and geometry.
///
/// Four parameters define the geometry of the calorimeter :
///
/// - the thickness of an absorber plate,
/// - the thickness of a gap,
/// - the number of layers,
/// - the transverse size of the calorimeter (the input face is a square).


class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    DetectorConstruction();
    virtual ~DetectorConstruction();
    
    double espesorSi;
    double posicionMuestraRadioactiva;
    int posicion;
    
    double presionCamara;

public:
    virtual G4VPhysicalVolume* Construct();

    // get methods
    const G4VPhysicalVolume* GetAbsorberPV() const;
    double getEspesorCapaSi(){return this->espesorSi;};
    double getPosicionMuestraRadioactiva(){return this->posicionMuestraRadioactiva;};
    int getPosicionBandeja(){return this->posicion;};
    
    double getPresionCamara(){return this->presionCamara;};
    
    //
    
    //set methods
    void setEspesorCapaSi(double espesor);
    void setPosicionMuestraRadioactiva(double posicionMuestra);
    void setPosicionBandeja(int posicion);
    
    void setPresionCamara(double presion);
    
private:
    // methods
    //
    G4VPhysicalVolume* DefineVolumes();

    // data members
    //

    G4VPhysicalVolume*   fAbsorberPV; // the absorber physical volume
    G4bool  fCheckOverlaps; // option to activate checking of volumes overlaps
    
    DetectorConstructionMessenger* detectorMessenger;

};

// inline functions
inline const G4VPhysicalVolume* DetectorConstruction::GetAbsorberPV() const { 
  return fAbsorberPV; 
}

#endif

