//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4NucleiPropertiesTheoreticalTableA.cc 67971 2013-03-13 10:13:24Z gcosmo $
//
//
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
// ------------------------------------------------------------
// Remove "theInstance"  by H.Kurashige (12 Dec. 03)

#include "../include/G4NucleiPropertiesTheoreticalTable.hh"

const double nMassC2=939.56563;
const double electron_mass_c2 = 0.510998910;
const double amu_c2 = 931.494028;
const double eV=1.e-6;
const double keV=1.e-3;

// Determine the table index for a Nuclide with Z protons and A nucleons
int G4NucleiPropertiesTheoreticalTable::GetIndex(int Z, int A)
{

  if(A>339) {
    std::cout << "Nucleon number larger than 339" << std::endl;
    return -1;
  } else if(A<16) {
    std::cout <<  " Nucleon number smaller than 16" << std::endl;
    return -1;
  } else if(Z>136) {
    std::cout << "Proton number larger than 136" << std::endl;
    return -1;
  } else if(Z<8) {
    std::cout << "Proton number smaller than 8" << std::endl;
    return -1;
  } else if(Z>A) {
    std::cout << "Nucleon number smaller than Z" << std::endl;
    return -1;
  }

  int i = shortTable[Z-8];
  while ( i < shortTable[Z-8+1] ) {
    if (indexArray[1][i] != A ) i++;
    else return i;
  }

  return -1;
}



double G4NucleiPropertiesTheoreticalTable::GetMassExcess(int Z, int A)
{
  int i=GetIndex(Z, A);
  if (i >= 0) {
    return AtomicMassExcess[i];
  } else {
    return 0.0;
  }
}

double G4NucleiPropertiesTheoreticalTable::GetBindingEnergy(int Z, int A)
{
  int i=GetIndex(Z, A);
  if (i >= 0){
    const double Mh = 7.289034;  // hydrogen atom mass excess
    const double Mn = 8.071431;  // neutron mass excess
    return double(Z)*Mh + double(A-Z)*Mn - AtomicMassExcess[i];
  } else {
    return 0.0;
  }
}



double  G4NucleiPropertiesTheoreticalTable::GetAtomicMass(int Z, int A)
{
  int i=GetIndex(Z, A);
  if (i >= 0) {
    return AtomicMassExcess[i] + A*amu_c2;
    } else {
      return 0.0;
    }
}



double  G4NucleiPropertiesTheoreticalTable::GetNuclearMass(int Z, int A)
{
  int i=GetIndex(Z, A);
  if (i >= 0) {
    return GetAtomicMass(Z,A) - double(Z)*electron_mass_c2 + ElectronicBindingEnergy(Z);
  } else {
    return 0.0;
  }
}

double G4NucleiPropertiesTheoreticalTable::ElectronicBindingEnergy(int Z) {
  const double ael = 1.433e-5; // electronic-binding constant
  return ael*std::pow(double(Z),2.39);
}

bool G4NucleiPropertiesTheoreticalTable::IsInTable(int Z, int A)
{
  return (Z <= A && A >= 16 && A <= 339 && Z <= 136 && Z >= 8 && GetIndex(Z, A) >= 0);
}








