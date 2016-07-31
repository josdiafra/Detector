
#include "DetectorConstruction.hh"
#include "G4SubtractionSolid.hh" 
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4RunManager.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "DetectorConstructionMessenger.hh"

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(),fAbsorberPV(0),fCheckOverlaps(true)
{
   G4GeometryManager::GetInstance()->OpenGeometry();
   G4PhysicalVolumeStore::GetInstance()->Clean();
   G4LogicalVolumeStore::GetInstance()->Clean();
   G4SolidStore::GetInstance()->Clean();

   this->espesorSi = (350/2)*um;
   this->posicionMuestraRadioactiva = -(75.36000/2)*mm + 2*(6.40000/2)*mm + 41*(2.00000/2)*mm + 2*(1.00000/2)*mm;
   
   this->presionCamara = (0.022/760)*atmosphere;
   
   
   
     // Creamos un messenger para esta clase.
   detectorMessenger = new DetectorConstructionMessenger(this);
}


DetectorConstruction::~DetectorConstruction()
{ 
    delete detectorMessenger;    
}


G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Define volumes
  return DefineVolumes();
}

G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
    //Instanciamos G4NistManager para definir los materiales
    G4NistManager* nist = G4NistManager::Instance();
    
    G4Material* aire = nist->FindOrBuildMaterial("G4_AIR");
    
    //Material de la capa de sicilicio del detector
    G4Material* Si = nist->FindOrBuildMaterial("G4_Si");
    
    //Vacio de la camara a 0.022 torr de presion

    G4double pressure = this->presionCamara;
    G4double density;
    G4int ncomponents = 2;
    G4double temperature = 293.15 * kelvin;
    G4double fractionmass;
    
    
    G4Element* N  = nist->FindOrBuildElement("N");
    G4Element* O  = nist->FindOrBuildElement("O");
    
    
    G4Material* vacuum = new G4Material("Vacuum",density=3.e-5*g/cm3,ncomponents,kStateGas,temperature,pressure);
    vacuum->AddElement(N, fractionmass=0.7);
    vacuum->AddElement(O, fractionmass=0.3);    
    
    //Material de la bandeja porta-muestra
    G4Material* Al = nist->FindOrBuildMaterial("G4_Al");
    
    //Vamos  definir ahora nuestro propio acero para el detector
    G4Element* SiElement  = nist->FindOrBuildElement("Si");
    G4Element* C  = nist->FindOrBuildElement("C");
    G4Element* Cr = nist->FindOrBuildElement("Cr");
    G4Element* Mn = nist->FindOrBuildElement("Mn");
    G4Element* Fe = nist->FindOrBuildElement("Fe");
    G4Element* Ni = nist->FindOrBuildElement("Ni");
    
    G4Material* StainlessSteel = new G4Material("StainlessSteel", density= 8.06*g/cm3, ncomponents=6);
    StainlessSteel->AddElement(C, fractionmass=0.001);
    StainlessSteel->AddElement(SiElement, fractionmass=0.007);
    StainlessSteel->AddElement(Cr, fractionmass=0.18);
    StainlessSteel->AddElement(Mn, fractionmass=0.01);
    StainlessSteel->AddElement(Fe, fractionmass=0.712);
    StainlessSteel->AddElement(Ni, fractionmass=0.09);

    //Material detector (caucho (rubber))
    G4Material* rubber = nist->FindOrBuildMaterial("G4_RUBBER_BUTYL");
    
    //Material detector (laton (brass))
    G4Material* laton = nist->FindOrBuildMaterial("G4_BRASS");
    
    
    //Material detector (delrin (polietileno))
    G4Material* delrin = nist->FindOrBuildMaterial("G4_POLYETHYLENE");
    
    // Muestra Am Am(OH)3
    G4Element* Am  = nist->FindOrBuildElement("Am");
    G4Element* H  = nist->FindOrBuildElement("H");
    //G4Material* dioxidoUranioMaterial = nist->FindOrBuildMaterial("G4_URANIUM_OXIDE");
    G4Material* muestraAm = new G4Material("muestraAm", density= 13.67*g/cm3, ncomponents=3);
    muestraAm->AddElement(O, fractionmass=0.165);
    muestraAm->AddElement(H, fractionmass=0.010);
    muestraAm->AddElement(Am, fractionmass=0.825);
    
    //Word Volume   
    G4double worldSizeX = (68/2)*mm;
    G4double worldSizeY = (75.36/2)*mm;
    G4double worldSizeZ =  (56.25/2)*mm;
    
    G4VSolid* worldS 
    = new G4Box("World",           // its name
                 worldSizeX, worldSizeY, worldSizeZ); // its size
    
    G4LogicalVolume* worldLV
    = new G4LogicalVolume(
                 worldS,           // its solid
                 vacuum,  // its material
                 "World");         // its name
    
    G4VPhysicalVolume* worldPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 worldLV,          // its logical volume                         
                 "World",          // its name
                 0,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
    
    
    G4double pRminAcero2 = 0.*cm;
    G4double pRmaxAcero2 = (16.0/2)*mm;
    G4double pDzAcero2 = (7.31000/2)*mm;
    G4double pSPhiAcero2 = 0;
    G4double pDPhiAcero2 = 2 * pi;
    
    
    G4Tubs* cilindroAcero2 = new G4Tubs("cilindroAcero2",pRminAcero2,pRmaxAcero2,pDzAcero2,pSPhiAcero2,pDPhiAcero2);
    
    
    G4LogicalVolume* cilindroAcero2LV = new G4LogicalVolume(cilindroAcero2,StainlessSteel,"cilindroAcero2");
    
    G4RotationMatrix * xRot = new G4RotationMatrix;  // Rotates X 
    xRot->rotateX(-90.*deg);                     // 
    
    G4double traslacionCilindroAcero2 = worldSizeY - pDzAcero2;
    
    G4VPhysicalVolume* cilindroAcero2PV
    = new G4PVPlacement(
                        xRot,                // no rotation
                        G4ThreeVector(0,traslacionCilindroAcero2,0),  // at (0,0,0)
                        cilindroAcero2LV,          // its logical volume
                        "cilindroAcero2",          // its name
                        worldLV,                // its mother  volume
                        false,            // no boolean operation
                        0,                // copy number
                        fCheckOverlaps);  // checking overlaps
    
    
    
    
  // Cilindro acero1 (macizoDetector) 
    G4double pRmin1 = 0.*cm;
    G4double pRmax1 = (32.00/2)*mm;
    G4double pDz1 = (12.3500/2)*mm;
    G4double pSPhi1 = 0;
    G4double pDPhi1 = 2 * pi;
    
    
    G4Tubs* cilindroAcero = new G4Tubs("cilindroAcero",pRmin1,pRmax1,pDz1,pSPhi1,pDPhi1);
    
    
    G4LogicalVolume* cilindroAceroLV = new G4LogicalVolume(cilindroAcero,StainlessSteel,"cilindroAcero");
                    // 
    
    G4double traslacionDetector = worldSizeY - pDz1-2*pDzAcero2;
    
    G4VPhysicalVolume* cilindroAceroPV
    = new G4PVPlacement(
                        xRot,                // no rotation
                        G4ThreeVector(0,traslacionDetector,0),  // at (0,0,0)
                        cilindroAceroLV,          // its logical volume
                        "cilindroAcero",          // its name
                        worldLV,                // its mother  volume
                        false,            // no boolean operation
                        0,                // copy number
                        fCheckOverlaps);  // checking overlaps
    
    
    //Construimos la capa muerta
    
    G4double pRminDeadLayer = 0.*cm;
    G4double pRmaxDeadLayer = (23.9/2)*mm;
    G4double pDzDeadLayer = (5.0/2)*nm;
    G4double pSPhiDeadLayer = 0;
    G4double pDPhiDeadLayer = 2 * pi;

    G4double trasSiDeadLayer = pDz1 - pDzDeadLayer; //Translacion del deadLayer
    
    G4Tubs* capaSilicioDeadLayer = new G4Tubs("capaSilicioDeadLayer",pRminDeadLayer,pRmaxDeadLayer,pDzDeadLayer,pSPhiDeadLayer,pDPhiDeadLayer);
    
    G4LogicalVolume* capaSilicioDeadLayerLV = new G4LogicalVolume(capaSilicioDeadLayer,Si,"capaSilicioDeadLayer");
    
    G4VPhysicalVolume* capaSilicioDeadLayerPV
    = new G4PVPlacement(
                        0,                // no rotation
                        (G4ThreeVector(0,0,trasSiDeadLayer)),  // at (0,0,0)
                        capaSilicioDeadLayerLV,          // its logical volume
                        "capaSilicioDeadLayer",          // its name
                        cilindroAceroLV,      // its mother  volume
                        false,            // no boolean operation
                        0,                // copy number
                        fCheckOverlaps);  // checking overlaps
    
    //Construimos la capa de silicio
    
    G4double pRmin3 = 0.*cm;
    G4double pRmax3 = (23.9000/2)*mm;
    G4double pDz3 = (this->espesorSi);
    G4double pSPhi3 = 0;
    G4double pDPhi3 = 2 * pi;

    G4double trasSi = pDz1 - 2*pDzDeadLayer - pDz3; //Translacion de la capa de silicio
    
    G4Tubs* capaSilicio = new G4Tubs("capaSilicio",pRmin3,pRmax3,pDz3,pSPhi3,pDPhi3);
    
    G4LogicalVolume* capaSilicioLV = new G4LogicalVolume(capaSilicio,Si,"capaSilicio");
    
    fAbsorberPV
    = new G4PVPlacement(
                        0,                // no rotation
                        (G4ThreeVector(0,0,trasSi)),  // at (0,0,0)
                        capaSilicioLV,          // its logical volume
                        "capaSilicio",          // its name
                        cilindroAceroLV,      // its mother  volume
                        false,            // no boolean operation
                        0,                // copy number
                        fCheckOverlaps);  // checking overlaps

    
    
    //Construimos la capa de caucho
    G4double pRmin4 = 0.*cm;
    G4double pRmax4 = (23.9/2)*mm;
    G4double pDz4 = ((0.5000/2))*mm;
    G4double pSPhi4 = 0;
    G4double pDPhi4 = 2 * pi;

    G4double trasRubber = pDz1- 2*pDzDeadLayer-2*pDz3-pDz4; //Translacion de la capa de caucho
    
    G4Tubs* capaRubber = new G4Tubs("capaRubber",pRmin4,pRmax4,pDz4,pSPhi4,pDPhi4);
    
    G4LogicalVolume* capaRubberLV = new G4LogicalVolume(capaRubber,rubber,"capaRubber");
    
    G4VPhysicalVolume* capaRubberPV
    = new G4PVPlacement(
                        0,                // no rotation
                        (G4ThreeVector(0,0,trasRubber)),  // at (0,0,0)
                        capaRubberLV,          // its logical volume
                        "capaRubber",          // its name
                        cilindroAceroLV,      // its mother  volume
                        false,            // no boolean operation
                        0,                // copy number
                        fCheckOverlaps);  // checking overlaps
    
    
    //Construimos la capa de laton
    G4double pRmin5 = 0.*cm;
    G4double pRmax5 = (23.9/2)*mm;
    G4double pDz5 = ((0.5000))*mm;
    G4double pSPhi5 = 0;
    G4double pDPhi5 = 2 * pi;

    G4double trasLaton = pDz1- 2*pDzDeadLayer-2*pDz3-2*pDz4-pDz5; //Translacion de la capa de caucho
    
    G4Tubs* capaLaton = new G4Tubs("capaLaton",pRmin5,pRmax5,pDz5,pSPhi5,pDPhi5);
    
    G4LogicalVolume* capaLatonLV = new G4LogicalVolume(capaLaton,rubber,"capaLaton");
    
    G4VPhysicalVolume* capaLatonPV
    = new G4PVPlacement(
                        0,                // no rotation
                        (G4ThreeVector(0,0,trasLaton)),  // at (0,0,0)
                        capaLatonLV,          // its logical volume
                        "capaLaton",          // its name
                        cilindroAceroLV,      // its mother  volume
                        false,            // no boolean operation
                        0,                // copy number
                        fCheckOverlaps);  // checking overlaps
    
    //Construimos la capa de delrin
    G4double pRmin6 = 0.*cm;
    G4double pRmax6 = (23.9/2)*mm;
    G4double pDz6 = ((4/2))*mm;
    G4double pSPhi6 = 0;
    G4double pDPhi6 = 2 * pi;

    G4double trasDelrin = pDz1- 2*pDzDeadLayer-2*pDz3-2*pDz4-2*pDz5-pDz6; //Translacion de la capa de caucho
    
    G4Tubs* capaDelrin = new G4Tubs("capaDelrin",pRmin6,pRmax6,pDz6,pSPhi6,pDPhi6);
    
    G4LogicalVolume* capaDelrinLV = new G4LogicalVolume(capaDelrin,delrin,"capaDelrin");
    
    G4VPhysicalVolume* capaDelrinPV
    = new G4PVPlacement(
                        0,                // no rotation
                        (G4ThreeVector(0,0,trasDelrin)),  // at (0,0,0)
                        capaDelrinLV,          // its logical volume
                        "capaDelrin",          // its name
                        cilindroAceroLV,      // its mother  volume
                        false,            // no boolean operation
                        0,                // copy number
                        fCheckOverlaps);  // checking overlaps
    
    
 //////////////////////////////////////////////////**MUESTRA**//////////////////////////////////////////////////////////	
    
    
    // Bandeja
    
    G4double bandejaX1 = worldSizeX;
    G4double bandejaY1 = (2.00000/2)*mm;
    G4double bandejaZ1 =  worldSizeZ-(2.30)*mm;
    
    G4VSolid* bandeja 
    = new G4Box("Bandeja",           // its name
                 bandejaX1, bandejaY1, bandejaZ1); // its size
    
    G4LogicalVolume* bandejaLV
    = new G4LogicalVolume(
                 bandeja,           // its solid
                 Al,  // its material
                 "Bandeja");         // its name
    
    
    G4double traslacionBandeja = this->posicionMuestraRadioactiva ;
    
    G4VPhysicalVolume* bandejaPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0,traslacionBandeja,0),  // at (0,0,0)
                 bandejaLV,          // its logical volume                         
                 "Bandeja",          // its name
                 worldLV,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
    
    //Plancheta
    
    G4double pRminPlancheta = 0.*cm;
    G4double pRmaxPlancheta = (1.2)*cm;
    G4double pDzPlancheta = (0.50000/2)*mm;
    G4double pSPhiPlancheta = 0;
    G4double pDPhiPlancheta = 2 * pi;
    
    
    G4Tubs* plancheta = new G4Tubs("plancheta",pRminPlancheta,pRmaxPlancheta,pDzPlancheta,pSPhiPlancheta,pDPhiPlancheta);
    
    
    G4LogicalVolume* planchetaLV = new G4LogicalVolume(plancheta,StainlessSteel,"plancheta");
    
    G4RotationMatrix * xRotplancheta = new G4RotationMatrix;  // Rotates X 
    xRotplancheta->rotateX(90.*deg); 
    
    G4double distanciaPlancheta = this->posicionMuestraRadioactiva + pDzPlancheta + bandejaY1;    
    
    G4VPhysicalVolume* planchetaPV
    = new G4PVPlacement(
                        xRotplancheta,                // no rotation
                        G4ThreeVector(0,distanciaPlancheta,0),  // at (0,0,0)
                        planchetaLV,          // its logical volume
                        "plancheta",          // its name
                        worldLV,      // its mother  volume
                        false,            // no boolean operation
                        0,                // copy number
                        fCheckOverlaps);  // checking overlaps
    
    
    //muestra Am
    
    G4double pRmin2 = 0.*cm;
    G4double pRmax2 = (0.9)*cm;
    G4double pDz2 = (40/2)*nm;
    G4double pSPhi2 = 0;
    G4double pDPhi2 = 2 * pi;
    
    
    G4Tubs* muestraAm241 = new G4Tubs("muestraAm241",pRmin2,pRmax2,pDz2,pSPhi2,pDPhi2);
    
    
    G4LogicalVolume* muestraAm241LV = new G4LogicalVolume(muestraAm241,muestraAm,"muestraAm241");
    
    G4RotationMatrix * xRot2 = new G4RotationMatrix;  // Rotates X 
    xRot2->rotateX(90.*deg); 
    
    G4double distanciaDetectorMuestra = distanciaPlancheta + pDz2 + pDzPlancheta;    
    
    G4VPhysicalVolume* muestraAm241PV
    = new G4PVPlacement(
                        xRot2,                // no rotation
                        G4ThreeVector(0,distanciaDetectorMuestra,0),  // at (0,0,0)
                        muestraAm241LV,          // its logical volume
                        "muestraAm241",          // its name
                        worldLV,      // its mother  volume
                        false,            // no boolean operation
                        0,                // copy number
                        fCheckOverlaps);  // checking overlaps
    
    //Soportes
    
    //Suplemento bandeja Left
    G4double suplementobandejaX1l = worldSizeX;
    G4double suplementobandejaY1l = (1.00/2)*mm;
    G4double suplementobandejaZ1l =  (4.0/2)*mm;
    
    G4VSolid* suplementobandejal 
    = new G4Box("suplementoBandejal",           // its name
                 suplementobandejaX1l, suplementobandejaY1l, suplementobandejaZ1l); // its size
    
    G4LogicalVolume* suplementobandejalLV
    = new G4LogicalVolume(
                 suplementobandejal,           // its solid
                 Al,  // its material
                 "suplementoBandejal");         // its name
    
    
    G4double traslsuplementoBandejaY1l = this->posicionMuestraRadioactiva - suplementobandejaY1l - bandejaY1;
    G4double traslsuplementoBandejaZ1l = worldSizeZ-suplementobandejaZ1l;    
    
    G4VPhysicalVolume* suplementobandejalPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0,traslsuplementoBandejaY1l,-traslsuplementoBandejaZ1l),  // at (0,0,0)
                 suplementobandejalLV,          // its logical volume                         
                 "suplementoBandejal",          // its name
                 worldLV,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
    
    
    //Suplemento bandeja Right
    G4double suplementobandejaX1r = worldSizeX;
    G4double suplementobandejaY1r = (1.00/2)*mm;
    G4double suplementobandejaZ1r =  (4.0/2)*mm;
    
    G4VSolid* suplementobandejar 
    = new G4Box("suplementoBandejar",           // its name
                 suplementobandejaX1r, suplementobandejaY1r, suplementobandejaZ1r); // its size
    
    G4LogicalVolume* suplementobandejarLV
    = new G4LogicalVolume(
                 suplementobandejar,           // its solid
                 Al,  // its material
                 "suplementoBandejar");         // its name
    
    
    G4double traslsuplementoBandejaY1r = this->posicionMuestraRadioactiva - suplementobandejaY1r - bandejaY1;
    G4double traslsuplementoBandejaZ1r = worldSizeZ-suplementobandejaZ1r;    
    
    G4VPhysicalVolume* suplementobandejarPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0,traslsuplementoBandejaY1r,traslsuplementoBandejaZ1r),  // at (0,0,0)
                 suplementobandejarLV,          // its logical volume                         
                 "suplementoBandejar",          // its name
                 worldLV,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
    
    
    
    //Soporte 1left
    
    G4double soporteX1l = worldSizeX;
    G4double soporteY1l = ( 6.40/2)*mm;
    G4double soporteZ1l =  (2.3/2)*mm;
    
    G4VSolid* soporte1l 
    = new G4Box("Soporte1l",           // its name
                 soporteX1l, soporteY1l, soporteZ1l); // its size
    
    G4LogicalVolume* soporte1lLV
    = new G4LogicalVolume(
                 soporte1l,           // its solid
                 Al,  // its material
                 "Soporte1l");         // its name
    
    
    G4double traslY1l = worldSizeY-soporteY1l;
    G4double traslZ1l = worldSizeZ-soporteZ1l;
    
    
    G4ThreeVector traslSoporte1l = G4ThreeVector(0,-traslY1l,-traslZ1l);
    
    
    G4VPhysicalVolume* soporte1lPV
    = new G4PVPlacement(
                 0,                // no rotation
                 traslSoporte1l,  // at (0,0,0)
                 soporte1lLV,          // its logical volume                         
                 "Soporte1l",          // its name
                 worldLV,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
    
    //Soporte 1rigth
    
    G4double soporteX1r = worldSizeX;
    G4double soporteY1r = (6.40/2)*mm;
    G4double soporteZ1r =  (2.3/2)*mm;
    
    G4VSolid* soporte1r 
    = new G4Box("Soporte1r",           // its name
                 soporteX1r, soporteY1r, soporteZ1r); // its size
    
    G4LogicalVolume* soporte1rLV
    = new G4LogicalVolume(
                 soporte1r,           // its solid
                 Al,  // its material
                 "Soporte1r");         // its name
    
    
    G4double traslY1r = worldSizeY-soporteY1r;
    G4double traslZ1r = worldSizeZ-soporteZ1r;
    
    
    G4ThreeVector traslSoporte1r = G4ThreeVector(0,-traslY1r,traslZ1r);
    
    
    G4VPhysicalVolume* soporte1rPV
    = new G4PVPlacement(
                 0,                // no rotation
                 traslSoporte1r,  // at (0,0,0)
                 soporte1rLV,          // its logical volume                         
                 "Soporte1r",          // its name
                 worldLV,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
    
    //Soporte 2left
    
    G4double soporteX2l = worldSizeX;
    G4double soporteY2l = (2.00/2)*mm;
    G4double soporteZ2l =  (2.3/2)*mm;
    
    G4VSolid* soporte2l 
    = new G4Box("Soporte2l",           // its name
                 soporteX2l, soporteY2l, soporteZ2l); // its size
    
    G4LogicalVolume* soporte2lLV
    = new G4LogicalVolume(
                 soporte2l,           // its solid
                 Al,  // its material
                 "Soporte2l");         // its name
    
    
    G4double traslY2l = traslY1l-soporteY1l-3*soporteY2l;
    G4double traslZ2l = worldSizeZ-soporteZ1l;
    
    
    G4ThreeVector traslSoporte2l = G4ThreeVector(0,-traslY2l,-traslZ2l);
    
    
    G4VPhysicalVolume* soporte2lPV
    = new G4PVPlacement(
                 0,                // no rotation
                 traslSoporte2l,  // at (0,0,0)
                 soporte2lLV,          // its logical volume                         
                 "Soporte2l",          // its name
                 worldLV,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
    
    
    //Soporte 2rigth
    
    G4double soporteX2r = worldSizeX;
    G4double soporteY2r = (2.00/2)*mm;
    G4double soporteZ2r =  (2.3/2)*mm;
    
    G4VSolid* soporte2r 
    = new G4Box("Soporte2r",           // its name
                 soporteX2r, soporteY2r, soporteZ2r); // its size
    
    G4LogicalVolume* soporte2rLV
    = new G4LogicalVolume(
                 soporte2r,           // its solid
                 Al,  // its material
                 "Soporte2r");         // its name
    
    
    G4double traslY2r = traslY1r-soporteY1r-3*soporteY2r;
    G4double traslZ2r = worldSizeZ-soporteZ1l;
    
    
    G4ThreeVector traslSoporte2r = G4ThreeVector(0,-traslY2r,traslZ2r);
    
    
    G4VPhysicalVolume* soporte2rPV
    = new G4PVPlacement(
                 0,                // no rotation
                 traslSoporte2r,  // at (0,0,0)
                 soporte2rLV,          // its logical volume                         
                 "Soporte2r",          // its name
                 worldLV,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
    
    
    //Soporte 3left
    
    G4double soporteX3l = worldSizeX;
    G4double soporteY3l = (2.00/2)*mm;
    G4double soporteZ3l =  (2.3/2)*mm;
    
    G4VSolid* soporte3l 
    = new G4Box("Soporte3l",           // its name
                 soporteX3l, soporteY3l, soporteZ3l); // its size
    
    G4LogicalVolume* soporte3lLV
    = new G4LogicalVolume(
                 soporte3l,           // its solid
                 Al,  // its material
                 "Soporte3l");         // its name
    
    
    G4double traslY3l = traslY2l-2*soporteY2l-2*bandejaY1;
    G4double traslZ3l = worldSizeZ-soporteZ2l;
    
    
    G4ThreeVector traslSoporte3l = G4ThreeVector(0,-traslY3l,-traslZ3l);
    
    
    G4VPhysicalVolume* soporte3lPV
    = new G4PVPlacement(
                 0,                // no rotation
                 traslSoporte3l,  // at (0,0,0)
                 soporte3lLV,          // its logical volume                         
                 "Soporte3l",          // its name
                 worldLV,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
    
    
   //Soporte 3 rigth
    
    G4double soporteX3r = worldSizeX;
    G4double soporteY3r = (2.00/2)*mm;
    G4double soporteZ3r =  (2.3/2)*mm;
    
    G4VSolid* soporte3r 
    = new G4Box("Soporte3r",           // its name
                 soporteX3r, soporteY3r, soporteZ3r); // its size
    
    G4LogicalVolume* soporte3rLV
    = new G4LogicalVolume(
                 soporte3r,           // its solid
                 Al,  // its material
                 "Soporte3r");         // its name
    
    
    G4double traslY3r = traslY2r-2*soporteY2r-2*bandejaY1;
    G4double traslZ3r = worldSizeZ-soporteZ2r;
    
    
    G4ThreeVector traslSoporte3r = G4ThreeVector(0,-traslY3r,traslZ3r);
    
    
    G4VPhysicalVolume* soporte3rPV
    = new G4PVPlacement(
                 0,                // no rotation
                 traslSoporte3r,  // at (0,0,0)
                 soporte3rLV,          // its logical volume                         
                 "Soporte3r",          // its name
                 worldLV,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
    
    
    //Soporte 4left
    
    G4double soporteX4l = worldSizeX;
    G4double soporteY4l = (2.00/2)*mm;
    G4double soporteZ4l =  (2.3/2)*mm;
    
    G4VSolid* soporte4l 
    = new G4Box("Soporte4l",           // its name
                 soporteX4l, soporteY4l, soporteZ4l); // its size
    
    G4LogicalVolume* soporte4lLV
    = new G4LogicalVolume(
                 soporte4l,           // its solid
                 Al,  // its material
                 "Soporte4l");         // its name
    
    
    G4double traslY4l = traslY3l-2*soporteY3l-2*bandejaY1;
    G4double traslZ4l = worldSizeZ-soporteZ3l;
    
    
    G4ThreeVector traslSoporte4l = G4ThreeVector(0,-traslY4l,-traslZ4l);
    
    
    G4VPhysicalVolume* soporte4lPV
    = new G4PVPlacement(
                 0,                // no rotation
                 traslSoporte4l,  // at (0,0,0)
                 soporte4lLV,          // its logical volume                         
                 "Soporte4l",          // its name
                 worldLV,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
    
    //Soporte 4 rigth
    
    G4double soporteX4r = worldSizeX;
    G4double soporteY4r = (2.00/2)*mm;
    G4double soporteZ4r =  (2.3/2)*mm;
    
    G4VSolid* soporte4r 
    = new G4Box("Soporte4r",           // its name
                 soporteX4r, soporteY4r, soporteZ4r); // its size
    
    G4LogicalVolume* soporte4rLV
    = new G4LogicalVolume(
                 soporte4r,           // its solid
                 Al,  // its material
                 "Soporte4r");         // its name
    
    
    G4double traslY4r = traslY3r-2*soporteY3r-2*bandejaY1;
    G4double traslZ4r = worldSizeZ-soporteZ3r;
    
    
    G4ThreeVector traslSoporte4r = G4ThreeVector(0,-traslY4r,traslZ4r);
    
    
    G4VPhysicalVolume* soporte4rPV
    = new G4PVPlacement(
                 0,                // no rotation
                 traslSoporte4r,  // at (0,0,0)
                 soporte4rLV,          // its logical volume                         
                 "Soporte4r",          // its name
                 worldLV,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
    
    //Soporte 5 left
    
    G4double soporteX5l = worldSizeX;
    G4double soporteY5l = (2.00/2)*mm;
    G4double soporteZ5l =  (2.3/2)*mm;
    
    G4VSolid* soporte5l 
    = new G4Box("Soporte5l",           // its name
                 soporteX5l, soporteY5l, soporteZ5l); // its size
    
    G4LogicalVolume* soporte5lLV
    = new G4LogicalVolume(
                 soporte5l,           // its solid
                 Al,  // its material
                 "Soporte4l");         // its name
    
    
    G4double traslY5l = traslY4l-2*soporteY4l-2*bandejaY1;
    G4double traslZ5l = worldSizeZ-soporteZ4l;
    
    
    G4ThreeVector traslSoporte5l = G4ThreeVector(0,-traslY5l,-traslZ5l);
    
    
    G4VPhysicalVolume* soporte5lPV
    = new G4PVPlacement(
                 0,                // no rotation
                 traslSoporte5l,  // at (0,0,0)
                 soporte5lLV,          // its logical volume                         
                 "Soporte5l",          // its name
                 worldLV,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
    
    
    //Soporte 5 rigth
    
    G4double soporteX5r = worldSizeX;
    G4double soporteY5r = (2.00/2)*mm;
    G4double soporteZ5r =  (2.3/2)*mm;
    
    G4VSolid* soporte5r 
    = new G4Box("Soporte5r",           // its name
                 soporteX5r, soporteY5r, soporteZ5r); // its size
    
    G4LogicalVolume* soporte5rLV
    = new G4LogicalVolume(
                 soporte5r,           // its solid
                 Al,  // its material
                 "Soporte5r");         // its name
    
    
    G4double traslY5r = traslY4r-2*soporteY4r-2*bandejaY1;
    G4double traslZ5r = worldSizeZ-soporteZ4r;
    
    
    G4ThreeVector traslSoporte5r = G4ThreeVector(0,-traslY5r,traslZ5r);
    
    
    G4VPhysicalVolume* soporte5rPV
    = new G4PVPlacement(
                 0,                // no rotation
                 traslSoporte5r,  // at (0,0,0)
                 soporte5rLV,          // its logical volume                         
                 "Soporte5r",          // its name
                 worldLV,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
    
    //Soporte 6 left
    
    G4double soporteX6l = worldSizeX;
    G4double soporteY6l = (2.00/2)*mm;
    G4double soporteZ6l =  (2.3/2)*mm;
    
    G4VSolid* soporte6l 
    = new G4Box("Soporte6l",           // its name
                 soporteX6l, soporteY6l, soporteZ6l); // its size
    
    G4LogicalVolume* soporte6lLV
    = new G4LogicalVolume(
                 soporte6l,           // its solid
                 Al,  // its material
                 "Soporte6l");         // its name
    
    
    G4double traslY6l = traslY5l-2*soporteY5l-2*bandejaY1;
    G4double traslZ6l = worldSizeZ-soporteZ5l;
    
    
    G4ThreeVector traslSoporte6l = G4ThreeVector(0,-traslY6l,-traslZ6l);
    
    
    G4VPhysicalVolume* soporte6lPV
    = new G4PVPlacement(
                 0,                // no rotation
                 traslSoporte6l,  // at (0,0,0)
                 soporte6lLV,          // its logical volume                         
                 "Soporte6l",          // its name
                 worldLV,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
    
    
     //Soporte 6 rigth
    
    G4double soporteX6r = worldSizeX;
    G4double soporteY6r = (2.00/2)*mm;
    G4double soporteZ6r =  (2.3/2)*mm;
    
    G4VSolid* soporte6r 
    = new G4Box("Soporte6r",           // its name
                 soporteX6r, soporteY6r, soporteZ6r); // its size
    
    G4LogicalVolume* soporte6rLV
    = new G4LogicalVolume(
                 soporte6r,           // its solid
                 Al,  // its material
                 "Soporte6r");         // its name
    
    
    G4double traslY6r = traslY5r-2*soporteY5r-2*bandejaY1;
    G4double traslZ6r = worldSizeZ-soporteZ5r;
    
    
    G4ThreeVector traslSoporte6r = G4ThreeVector(0,-traslY6r,traslZ6r);
    
    
    G4VPhysicalVolume* soporte6rPV
    = new G4PVPlacement(
                 0,                // no rotation
                 traslSoporte6r,  // at (0,0,0)
                 soporte6rLV,          // its logical volume                         
                 "Soporte6r",          // its name
                 worldLV,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
    
    
    //Soporte 7 left
    
    G4double soporteX7l = worldSizeX;
    G4double soporteY7l = (2.00/2)*mm;
    G4double soporteZ7l =  (2.3/2)*mm;
    
    G4VSolid* soporte7l 
    = new G4Box("Soporte7l",           // its name
                 soporteX7l, soporteY7l, soporteZ7l); // its size
    
    G4LogicalVolume* soporte7lLV
    = new G4LogicalVolume(
                 soporte7l,           // its solid
                 Al,  // its material
                 "Soporte7l");         // its name
    
    
    G4double traslY7l = traslY6l-2*soporteY6l-2*bandejaY1;
    G4double traslZ7l = worldSizeZ-soporteZ6l;
    
    
    G4ThreeVector traslSoporte7l = G4ThreeVector(0,-traslY7l,-traslZ7l);
    
    
    G4VPhysicalVolume* soporte7lPV
    = new G4PVPlacement(
                 0,                // no rotation
                 traslSoporte7l,  // at (0,0,0)
                 soporte7lLV,          // its logical volume                         
                 "Soporte7l",          // its name
                 worldLV,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
    
    
    //Soporte 7 rigth
    
    G4double soporteX7r = worldSizeX;
    G4double soporteY7r = (2.00/2)*mm;
    G4double soporteZ7r =  (2.3/2)*mm;
    
    G4VSolid* soporte7r 
    = new G4Box("Soporte7r",           // its name
                 soporteX7r, soporteY7r, soporteZ7r); // its size
    
    G4LogicalVolume* soporte7rLV
    = new G4LogicalVolume(
                 soporte7r,           // its solid
                 Al,  // its material
                 "Soporte7r");         // its name
    
    
    G4double traslY7r = traslY6r-2*soporteY6r-2*bandejaY1;
    G4double traslZ7r = worldSizeZ-soporteZ6r;
    
    
    G4ThreeVector traslSoporte7r = G4ThreeVector(0,-traslY7r,traslZ7r);
    
    
    G4VPhysicalVolume* soporte7rPV
    = new G4PVPlacement(
                 0,                // no rotation
                 traslSoporte7r,  // at (0,0,0)
                 soporte7rLV,          // its logical volume                         
                 "Soporte7r",          // its name
                 worldLV,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
    
    
    
    //Soporte 8 left
    
    G4double soporteX8l = worldSizeX;
    G4double soporteY8l = (2.00/2)*mm;
    G4double soporteZ8l =  (2.3/2)*mm;
    
    G4VSolid* soporte8l 
    = new G4Box("Soporte8l",           // its name
                 soporteX8l, soporteY8l, soporteZ8l); // its size
    
    G4LogicalVolume* soporte8lLV
    = new G4LogicalVolume(
                 soporte8l,           // its solid
                 Al,  // its material
                 "Soporte8l");         // its name
    
    
    G4double traslY8l = traslY7l-2*soporteY7l-2*bandejaY1;
    G4double traslZ8l = worldSizeZ-soporteZ7l;
    
    
    G4ThreeVector traslSoporte8l = G4ThreeVector(0,-traslY8l,-traslZ8l);
    
    
    G4VPhysicalVolume* soporte8lPV
    = new G4PVPlacement(
                 0,                // no rotation
                 traslSoporte8l,  // at (0,0,0)
                 soporte8lLV,          // its logical volume                         
                 "Soporte8l",          // its name
                 worldLV,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
    
    
    
     //Soporte 8 rigth
    
    G4double soporteX8r = worldSizeX;
    G4double soporteY8r = (2.00/2)*mm;
    G4double soporteZ8r =  (2.3/2)*mm;
    
    G4VSolid* soporte8r 
    = new G4Box("Soporte8r",           // its name
                 soporteX8r, soporteY8r, soporteZ8r); // its size
    
    G4LogicalVolume* soporte8rLV
    = new G4LogicalVolume(
                 soporte8r,           // its solid
                 Al,  // its material
                 "Soporte8r");         // its name
    
    
    G4double traslY8r = traslY7r-2*soporteY7r-2*bandejaY1;
    G4double traslZ8r = worldSizeZ-soporteZ7r;
    
    
    G4ThreeVector traslSoporte8r = G4ThreeVector(0,-traslY8r,traslZ8r);
    
    
    G4VPhysicalVolume* soporte8rPV
    = new G4PVPlacement(
                 0,                // no rotation
                 traslSoporte8r,  // at (0,0,0)
                 soporte8rLV,          // its logical volume                         
                 "Soporte8r",          // its name
                 worldLV,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
    
    
    //Soporte 9 left
    
    G4double soporteX9l = worldSizeX;
    G4double soporteY9l = (2.00/2)*mm;
    G4double soporteZ9l =  (2.3/2)*mm;
    
    G4VSolid* soporte9l 
    = new G4Box("Soporte9l",           // its name
                 soporteX9l, soporteY9l, soporteZ9l); // its size
    
    G4LogicalVolume* soporte9lLV
    = new G4LogicalVolume(
                 soporte9l,           // its solid
                 Al,  // its material
                 "Soporte9l");         // its name
    
    
    G4double traslY9l = traslY8l-2*soporteY8l-2*bandejaY1;
    G4double traslZ9l = worldSizeZ-soporteZ8l;
    
    
    G4ThreeVector traslSoporte9l = G4ThreeVector(0,-traslY9l,-traslZ9l);
    
    
    G4VPhysicalVolume* soporte9lPV
    = new G4PVPlacement(
                 0,                // no rotation
                 traslSoporte9l,  // at (0,0,0)
                 soporte9lLV,          // its logical volume                         
                 "Soporte9l",          // its name
                 worldLV,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
    
    
    //Soporte 9 rigth
    
    G4double soporteX9r = worldSizeX;
    G4double soporteY9r = (2.00/2)*mm;
    G4double soporteZ9r =  (2.3/2)*mm;
    
    G4VSolid* soporte9r 
    = new G4Box("Soporte9r",           // its name
                 soporteX9r, soporteY9r, soporteZ9r); // its size
    
    G4LogicalVolume* soporte9rLV
    = new G4LogicalVolume(
                 soporte9r,           // its solid
                 Al,  // its material
                 "Soporte9r");         // its name
    
    
    G4double traslY9r = traslY8r-2*soporteY8r-2*bandejaY1;
    G4double traslZ9r = worldSizeZ-soporteZ8r;
    
    
    G4ThreeVector traslSoporte9r = G4ThreeVector(0,-traslY9r,traslZ9r);
    
    
    G4VPhysicalVolume* soporte9rPV
    = new G4PVPlacement(
                 0,                // no rotation
                 traslSoporte9r,  // at (0,0,0)
                 soporte9rLV,          // its logical volume                         
                 "Soporte9r",          // its name
                 worldLV,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
    
    
    //Soporte 10 left
    
    G4double soporteX10l = worldSizeX;
    G4double soporteY10l = (2.00/2)*mm;
    G4double soporteZ10l =  (2.3/2)*mm;
    
    G4VSolid* soporte10l 
    = new G4Box("Soporte10l",           // its name
                 soporteX10l, soporteY10l, soporteZ10l); // its size
    
    G4LogicalVolume* soporte10lLV
    = new G4LogicalVolume(
                 soporte10l,           // its solid
                 Al,  // its material
                 "Soporte10l");         // its name
    
    
    G4double traslY10l = traslY9l-2*soporteY9l-2*bandejaY1;
    G4double traslZ10l = worldSizeZ-soporteZ9l;
    
    
    G4ThreeVector traslSoporte10l = G4ThreeVector(0,-traslY10l,-traslZ10l);
    
    
    G4VPhysicalVolume* soporte10lPV
    = new G4PVPlacement(
                 0,                // no rotation
                 traslSoporte10l,  // at (0,0,0)
                 soporte10lLV,          // its logical volume                         
                 "Soporte10l",          // its name
                 worldLV,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
    
    
    //Soporte 10 rigth
    
    G4double soporteX10r = worldSizeX;
    G4double soporteY10r = (2.00/2)*mm;
    G4double soporteZ10r =  (2.3/2)*mm;
    
    G4VSolid* soporte10r 
    = new G4Box("Soporte10r",           // its name
                 soporteX10r, soporteY10r, soporteZ10r); // its size
    
    G4LogicalVolume* soporte10rLV
    = new G4LogicalVolume(
                 soporte10r,           // its solid
                 Al,  // its material
                 "Soporte10r");         // its name
    
    
    G4double traslY10r = traslY9r-2*soporteY9r-2*bandejaY1;
    G4double traslZ10r = worldSizeZ-soporteZ9r;
    
    
    G4ThreeVector traslSoporte10r = G4ThreeVector(0,-traslY10r,traslZ10r);
    
    
    G4VPhysicalVolume* soporte10rPV
    = new G4PVPlacement(
                 0,                // no rotation
                 traslSoporte10r,  // at (0,0,0)
                 soporte10rLV,          // its logical volume                         
                 "Soporte10r",          // its name
                 worldLV,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
    
    
    //Soporte 11 left
    
    G4double soporteX11l = worldSizeX;
    G4double soporteY11l = (2.00/2)*mm;
    G4double soporteZ11l =  (2.3/2)*mm;
    
    G4VSolid* soporte11l 
    = new G4Box("Soporte11l",           // its name
                 soporteX11l, soporteY11l, soporteZ11l); // its size
    
    G4LogicalVolume* soporte11lLV
    = new G4LogicalVolume(
                 soporte11l,           // its solid
                 Al,  // its material
                 "Soporte11l");         // its name
    
    
    G4double traslY11l = traslY10l-2*soporteY10l-2*bandejaY1;
    G4double traslZ11l = worldSizeZ-soporteZ10l;
    
    
    G4ThreeVector traslSoporte11l = G4ThreeVector(0,-traslY11l,-traslZ11l);
    
    
    G4VPhysicalVolume* soporte11lPV
    = new G4PVPlacement(
                 0,                // no rotation
                 traslSoporte11l,  // at (0,0,0)
                 soporte11lLV,          // its logical volume                         
                 "Soporte11l",          // its name
                 worldLV,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
    
    
    //Soporte 11 rigth
    
    G4double soporteX11r = worldSizeX;
    G4double soporteY11r = (2.00/2)*mm;
    G4double soporteZ11r =  (2.3/2)*mm;
    
    G4VSolid* soporte11r 
    = new G4Box("Soporte11r",           // its name
                 soporteX11r, soporteY11r, soporteZ11r); // its size
    
    G4LogicalVolume* soporte11rLV
    = new G4LogicalVolume(
                 soporte11r,           // its solid
                 Al,  // its material
                 "Soporte11r");         // its name
    
    
    G4double traslY11r = traslY10r-2*soporteY10r-2*bandejaY1;
    G4double traslZ11r = worldSizeZ-soporteZ10r;
    
    
    G4ThreeVector traslSoporte11r = G4ThreeVector(0,-traslY11r,traslZ11r);
    
    
    G4VPhysicalVolume* soporte11rPV
    = new G4PVPlacement(
                 0,                // no rotation
                 traslSoporte11r,  // at (0,0,0)
                 soporte11rLV,          // its logical volume                         
                 "Soporte11r",          // its name
                 worldLV,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
    
    
     //Soporte 12 left
    
    G4double soporteX12l = worldSizeX;
    G4double soporteY12l = (2.00/2)*mm;
    G4double soporteZ12l =  (2.3/2)*mm;
    
    G4VSolid* soporte12l 
    = new G4Box("Soporte12l",           // its name
                 soporteX12l, soporteY12l, soporteZ12l); // its size
    
    G4LogicalVolume* soporte12lLV
    = new G4LogicalVolume(
                 soporte12l,           // its solid
                 Al,  // its material
                 "Soporte12l");         // its name
    
    
    G4double traslY12l = traslY11l-2*soporteY11l-2*bandejaY1;
    G4double traslZ12l = worldSizeZ-soporteZ11l;
    
    
    G4ThreeVector traslSoporte12l = G4ThreeVector(0,-traslY12l,-traslZ12l);
    
    
    G4VPhysicalVolume* soporte12lPV
    = new G4PVPlacement(
                 0,                // no rotation
                 traslSoporte12l,  // at (0,0,0)
                 soporte12lLV,          // its logical volume                         
                 "Soporte12l",          // its name
                 worldLV,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
    
    
    //Soporte 12 rigth
    
    G4double soporteX12r = worldSizeX;
    G4double soporteY12r = (2.00/2)*mm;
    G4double soporteZ12r =  (2.3/2)*mm;
    
    G4VSolid* soporte12r 
    = new G4Box("Soporte12r",           // its name
                 soporteX12r, soporteY12r, soporteZ12r); // its size
    
    G4LogicalVolume* soporte12rLV
    = new G4LogicalVolume(
                 soporte12r,           // its solid
                 Al,  // its material
                 "Soporte12r");         // its name
    
    
    G4double traslY12r = traslY11r-2*soporteY11r-2*bandejaY1;
    G4double traslZ12r = worldSizeZ-soporteZ11r;
    
    
    G4ThreeVector traslSoporte12r = G4ThreeVector(0,-traslY12r,traslZ12r);
    
    
    G4VPhysicalVolume* soporte12rPV
    = new G4PVPlacement(
                 0,                // no rotation
                 traslSoporte12r,  // at (0,0,0)
                 soporte12rLV,          // its logical volume                         
                 "Soporte12r",          // its name
                 worldLV,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
    
    
    //Soporte 13 left
    
    G4double soporteX13l = worldSizeX;
    G4double soporteY13l = (2.00/2)*mm;
    G4double soporteZ13l =  (2.3/2)*mm;
    
    G4VSolid* soporte13l 
    = new G4Box("Soporte13l",           // its name
                 soporteX13l, soporteY13l, soporteZ13l); // its size
    
    G4LogicalVolume* soporte13lLV
    = new G4LogicalVolume(
                 soporte13l,           // its solid
                 Al,  // its material
                 "Soporte13l");         // its name
    
    
    G4double traslY13l = traslY12l-2*soporteY12l-2*bandejaY1;
    G4double traslZ13l = worldSizeZ-soporteZ12l;
    
    
    G4ThreeVector traslSoporte13l = G4ThreeVector(0,-traslY13l,-traslZ13l);
    
    
    G4VPhysicalVolume* soporte13lPV
    = new G4PVPlacement(
                 0,                // no rotation
                 traslSoporte13l,  // at (0,0,0)
                 soporte13lLV,          // its logical volume                         
                 "Soporte13l",          // its name
                 worldLV,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
    
    
    //Soporte 13 rigth
    
    G4double soporteX13r = worldSizeX;
    G4double soporteY13r = (2.00/2)*mm;
    G4double soporteZ13r =  (2.3/2)*mm;
    
    G4VSolid* soporte13r 
    = new G4Box("Soporte13r",           // its name
                 soporteX13r, soporteY13r, soporteZ13r); // its size
    
    G4LogicalVolume* soporte13rLV
    = new G4LogicalVolume(
                 soporte13r,           // its solid
                 Al,  // its material
                 "Soporte13r");         // its name
    
    
    G4double traslY13r = traslY12r-2*soporteY12r-2*bandejaY1;
    G4double traslZ13r = worldSizeZ-soporteZ12r;
    
    
    G4ThreeVector traslSoporte13r = G4ThreeVector(0,-traslY13r,traslZ13r);
    
    
    G4VPhysicalVolume* soporte13rPV
    = new G4PVPlacement(
                 0,                // no rotation
                 traslSoporte13r,  // at (0,0,0)
                 soporte13rLV,          // its logical volume                         
                 "Soporte13r",          // its name
                 worldLV,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
    
    
    
    
    
    //
    // Visualization attributes
    //
    
    
    

    G4VisAttributes * detectorVista = new G4VisAttributes(G4Colour(0.,0.,1.));
    cilindroAceroLV->SetVisAttributes(detectorVista);
    
    G4VisAttributes * capaSiVista = new G4VisAttributes(G4Colour(0.,1.,1.));  //cyan
    capaSilicioLV->SetVisAttributes(capaSiVista);
    
    G4VisAttributes * capaRubberVista = new G4VisAttributes(G4Colour(1.,1.,0.));  //yellow
    capaRubberLV->SetVisAttributes(capaRubberVista);
    
    G4VisAttributes * capaDelrinVista = new G4VisAttributes(G4Colour(0.,1.,0.));  //red
    capaDelrinLV->SetVisAttributes(capaDelrinVista);
    
    G4VisAttributes * capaLatonVista = new G4VisAttributes(G4Colour(1.,0.,0.));  //red
    capaLatonLV->SetVisAttributes(capaLatonVista);
    
    G4VisAttributes * muestraAm241Vista = new G4VisAttributes(G4Colour(0.5,0.5,0.5));  //gray
    muestraAm241LV->SetVisAttributes(muestraAm241Vista);
    
    
    //
    // Always return the physical World
    return worldPV;
}

void DetectorConstruction::setPosicionMuestraRadioactiva(double posicionMuestra)
{
  this->posicionMuestraRadioactiva = posicionMuestra;
  G4RunManager::GetRunManager()->ReinitializeGeometry();  
}

void DetectorConstruction::setEspesorCapaSi(G4double espesor)
{
  this->espesorSi = espesor;
  G4RunManager::GetRunManager()->ReinitializeGeometry();  
}

void DetectorConstruction::setPresionCamara(G4double presion)
{
  this->presionCamara = presion;
  G4RunManager::GetRunManager()->ReinitializeGeometry();  
}


void DetectorConstruction::setPosicionBandeja(G4int posicion)
{
  
    G4double soporteY1l = ( 6.40/2)*mm;
    G4double worldSizeY = (75.36/2)*mm;   
    G4double soporteY2l = (2/2)*mm;
    G4double suplementobandejaY1l = (1.0/2)*mm;
    
    switch(posicion){
        
        case 1:
        this->posicionMuestraRadioactiva = -worldSizeY+2*soporteY1l+soporteY2l + 2*suplementobandejaY1l;
        G4RunManager::GetRunManager()->ReinitializeGeometry();  
        break;
        
        case 2:
        this->posicionMuestraRadioactiva = -worldSizeY+2*soporteY1l+5*soporteY2l + 2*suplementobandejaY1l;
        G4RunManager::GetRunManager()->ReinitializeGeometry();  
        break;
        
        case 3:
        this->posicionMuestraRadioactiva = -worldSizeY+2*soporteY1l+9*soporteY2l + 2*suplementobandejaY1l;
        G4RunManager::GetRunManager()->ReinitializeGeometry();  
        break;
        
        case 4:
        this->posicionMuestraRadioactiva = -worldSizeY+2*soporteY1l+13*soporteY2l + 2*suplementobandejaY1l;
        G4RunManager::GetRunManager()->ReinitializeGeometry();  
        break;    
        
        case 5:
        this->posicionMuestraRadioactiva = -worldSizeY+2*soporteY1l+17*soporteY2l + 2*suplementobandejaY1l;
        G4RunManager::GetRunManager()->ReinitializeGeometry();  
        break;    
        
        case 6:
        this->posicionMuestraRadioactiva = -worldSizeY+2*soporteY1l+21*soporteY2l + 2*suplementobandejaY1l;
        G4RunManager::GetRunManager()->ReinitializeGeometry();  
        break; 
        
        case 7:
        this->posicionMuestraRadioactiva = -worldSizeY+2*soporteY1l+25*soporteY2l + 2*suplementobandejaY1l;
        G4RunManager::GetRunManager()->ReinitializeGeometry();  
        break;
        
        case 8:
        this->posicionMuestraRadioactiva = -worldSizeY+2*soporteY1l+29*soporteY2l + 2*suplementobandejaY1l;
        G4RunManager::GetRunManager()->ReinitializeGeometry();  
        break;   
        
        case 9:
        this->posicionMuestraRadioactiva = -worldSizeY+2*soporteY1l+33*soporteY2l+ 2*suplementobandejaY1l;
        G4RunManager::GetRunManager()->ReinitializeGeometry();  
        break;       
        
        case 10:
        this->posicionMuestraRadioactiva = -worldSizeY+2*soporteY1l+37*soporteY2l + 2*suplementobandejaY1l;
        G4RunManager::GetRunManager()->ReinitializeGeometry();  
        break; 
        
        case 11:
        this->posicionMuestraRadioactiva = -worldSizeY+2*soporteY1l+41*soporteY2l + 2*suplementobandejaY1l;
        G4RunManager::GetRunManager()->ReinitializeGeometry();  
        break;   
        
        case 12:
        this->posicionMuestraRadioactiva = -worldSizeY+2*soporteY1l+45*soporteY2l + 2*suplementobandejaY1l;
        G4RunManager::GetRunManager()->ReinitializeGeometry();  
        break;       
    }
    
    
}


