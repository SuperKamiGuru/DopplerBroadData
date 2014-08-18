// Created by: Wesley Ford
//This program either takes a single cross-section data file, a directory full of cross-section data files
// or a directory contianing the directories Elastic, Inelastic, Fission and Capture which contian cross-section data files,
// it sorts out which files are valid for doppler broadening, it doppler broadens the data based off the given temperature difference
// and then it outputs the doppler broadened data into the given output filename or a given output directory with the same
// internal structure as the input directory

// many of the ideas used in this file and many of the classes used by this file are copied from the GEANT4 source code.
// This was primarily done to make sure the outputed data is consitent with what the user would get when using the on the fly doppler broadening algorythm.
// A some modifications were made to the copied code however to improve the efficiency and portablitiy of this program.

#define Timer 2

using namespace std;

#include "include/Random/RandomEngine.h"
#include "include/Random/RanecuEngine.h"
//#include <CLHEP/Random/Randomize.h>
#include "include/Random/RandGauss.h"
#include "include/Random/RandGaussQ.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include "zlib.h"
#include <fstream>
#include <ctime>
#include "include/G4NucleiPropertiesTableAME03.hh"
#include "include/G4NucleiPropertiesTheoreticalTable.hh"
#include <dirent.h>
#include "include/ElementNames.hh"

// default mass units are m*c^2=MeV, the same as energy
// default momentum units are p*c=MeV, same as energy
// the velocity is kept in terms of v=v/c on purpose for the below calc
// see G4NucleiPropertiesTableAME03 for finding isotope nucleus mass

//use GetDataStream instead
//make the outputfile binary and compressed for speed

//Note the GetDataStream data function used in this program is just a slight modification of
// G4NeutronHPManager::GetDataStream and the doppler broadening algorythm is based off the
// G4NeutronHP[Process]Data::GetCrossSection

#define indexSize 100

const double k_Boltzmann=8.617343*pow(10,-11);
const double nMassC2=939.56563;
const double electron_mass_c2 = 0.510998910;
const double amu_c2 = 931.494028;
const double eV=1.e-6;
const double keV=1.e-3;

#if Timer>=1
    void ConvertDirect(string inDirName, string outDirName, double prevTemp, double newTemp, bool ascii, bool log, std::ofstream* logFile,
    int totalNumFiles, int &index, long totalFileSize2, long &sumFileSize2, double &sumDuration, std::vector<int> &fileSize2List);
    void printProgress(int index, std::vector<int> &fileSize2List, long totalFileSize2, int totalNumFiles, double duration, long sumFileSize2, double sumDuration);
    void GetfileSize2List(string inFileName, std::vector<int> &fileSize2List, int &totalNumFiles, long &totalFileSize2);
    void GetDirectoryFileSize(string inDirName, std::vector<int> &fileSize2List, int &index, long &totalFileSize2);
    bool isApplicable(string fileName);
    int GetFileSize2(string fileName);
    void ExtractZA(string fileName, int &Z, int &A);
#else
    void ConvertDirect(string inDirName, string outDirName, double prevTemp, double newTemp, bool ascii, bool log, std::ofstream* logFile);
#endif
void ConvertFile(string inDirName, string outDirName, string fileName, double prevTemp, double newTemp, bool ascii, bool log, std::ofstream* logFile);
void ConvertFile(string inFileName, string outFileName, double prevTemp, double newTemp, bool ascii, bool log, std::ofstream* logFile);
bool DirectoryExists( const char* pzPath );
void ExtractZA(string fileName, int &Z, int &A, bool log, std::ofstream* logFile);
double GetNuclearMass(double A, double Z, bool log, std::ofstream* logFile);
double NuclearMass(double A, double Z, bool log, std::ofstream* logFile);
double  AtomicMass(double A, double Z);
double  BindingEnergy(double A, double Z);
int DoppBroad(string inFileName, string outFileName, double prevTemp, double newTemp, double isoMassC2, bool ascii, bool log, std::ofstream* logFile);
double findCS(double nKEnerTrans, const double* prevEnVec, const double* prevCSVec, int vecSize);
void GetDataStream( string, std::stringstream&, bool, std::ofstream*);
void SetDataStream( string, std::stringstream&, bool ascii, bool, std::ofstream*);



int main(int argc, char **argv)
{
    #if Timer>=1
        std::clock_t start;
    #endif

    CLHEP::HepRandom *CLHEPRand = new CLHEP::HepRandom();

    CLHEP::RanecuEngine *theEngine = new CLHEP::RanecuEngine(2013092304);

    CLHEPRand->setTheEngine(theEngine);

    CLHEPRand->showEngineStatus();
    CLHEPRand->saveEngineStatus();

    ElementNames elementNames;

    elementNames.SetElementNames();

    //double duration;
    int result=0;
    string inFileName, outFileName;
    string inSubDirName, outSubDirName;
    string outputType="false", createLogFile="false";
    string logFileName;
    double newTemp;
    double prevTemp;
    bool ascii=true, log=false;
    stringstream ss;
    std::ofstream *logFile=NULL;

    if(argc==7)
    {
        ss << argv[1] << ' ' << argv[2] << ' ' << argv[3] << ' ' << argv[4] << ' ' << argv[5] << ' ' << argv[6];
        ss >> inFileName >> outFileName >> prevTemp >> newTemp >> outputType >> createLogFile;
    }
    else if(argc==6)
    {
        ss << argv[1] << ' ' << argv[2] << ' ' << argv[3] << ' ' << argv[4] << ' ' << argv[5];
        ss >> inFileName >> outFileName >> prevTemp >> newTemp >> outputType;
    }
    else
    {
        ss << argv[1] << ' ' << argv[2] << ' ' << argv[3] << ' ' << argv[4];
        ss >> inFileName >> outFileName >> prevTemp >> newTemp;
    }

    ss.str("");
    ss.clear();

    if(outputType == "compressed"||outputType == "compress"||outputType == "Compressed"||outputType == "Compress"||outputType == "Zipped"||outputType == "Zip"||outputType == "zipped"||outputType == "zip" )
        ascii=false;

    if(createLogFile == "true"||createLogFile == "True"||createLogFile == "1"||createLogFile == "On"||createLogFile == "on")
        log=true;

    if(log)
    {
        if (outFileName.back()=='/')
        {
            logFileName=outFileName+"LogDopplerBroadData.txt";
            logFile = new std::ofstream( logFileName.c_str() , std::ios::out | std::ios::app );
            if(!logFile->good())
            {
                cout << "### Error: could not open log file ###" << endl;
                return 1;
            }
        }
        else
        {
            size_t pos = outFileName.find_last_of('/');
            if(pos == std::string::npos)
            {
                logFileName="./LogDopplerBroadData.txt";
                logFile = new std::ofstream( logFileName.c_str() , std::ios::out | std::ios::app );
                if(!logFile->good())
                {
                    cout << "### Error: could not open log file ###" << endl;
                    return 1;
                }
            }
            else
            {
                logFileName=outFileName.substr(0,pos+1)+"LogDopplerBroadData.txt";
                logFile = new std::ofstream( logFileName.c_str() , std::ios::out | std::ios::app );
                if(!logFile->good())
                {
                    cout << "### Error: could not open log file ###" << endl;
                    return 1;
                }
            }
        }

        logFile->fill('#');
        (*logFile) << std::setw(84) << std::left << ""<< endl;
        (*logFile) << "The Input File Name is " << inFileName << endl;
        (*logFile) << "The Output File Name is " << outFileName << endl;
        (*logFile) << "The Temperature Difference is " << (newTemp-prevTemp) << endl;
        (*logFile) << "The Output Format is " << ((ascii)? "ASCII":"Compressed") << endl;

    }

    if (inFileName.back()=='/')
    {
        #if Timer>=1
            std::vector<int> fileSize2List(1400);
            int totalNumFiles=0, index=0;
            long totalFileSize2=0, sumFileSize2=0;
            GetfileSize2List(inFileName, fileSize2List, totalNumFiles, totalFileSize2);
            double sumDuration=0, duration=0;
        #endif

        DIR *dir;
        struct dirent *ent;
        if ((dir = opendir (inFileName.c_str())) != NULL)
        {
          /* print all the files and directories within directory */
          while ((ent = readdir (dir)) != NULL)
          {
            if(string(ent->d_name)=="Elastic")
            {
                inSubDirName = inFileName + "Elastic/CrossSection/";
                outSubDirName = outFileName + "Elastic/CrossSection/";

                #if Timer>=1
                    ConvertDirect(inSubDirName, outSubDirName, prevTemp, newTemp, ascii, log, logFile, totalNumFiles, index, totalFileSize2, sumFileSize2, sumDuration, fileSize2List);
                #else
                    ConvertDirect(inSubDirName, outSubDirName, prevTemp, newTemp, ascii, log, logFile);
                #endif

            }
            else if(string(ent->d_name)=="Inelastic")
            {
                inSubDirName = inFileName + "Inelastic/CrossSection/";
                outSubDirName = outFileName + "Inelastic/CrossSection/";

                #if Timer>=1
                    ConvertDirect(inSubDirName, outSubDirName, prevTemp, newTemp, ascii, log, logFile, totalNumFiles, index, totalFileSize2, sumFileSize2, sumDuration, fileSize2List);
                #else
                    ConvertDirect(inSubDirName, outSubDirName, prevTemp, newTemp, ascii, log, logFile);
                #endif
            }
            else if(string(ent->d_name)=="Fission")
            {
                inSubDirName = inFileName + "Fission/CrossSection/";
                outSubDirName = outFileName + "Fission/CrossSection/";

                #if Timer>=1
                    ConvertDirect(inSubDirName, outSubDirName, prevTemp, newTemp, ascii, log, logFile, totalNumFiles, index, totalFileSize2, sumFileSize2, sumDuration, fileSize2List);
                #else
                    ConvertDirect(inSubDirName, outSubDirName, prevTemp, newTemp, ascii, log, logFile);
                #endif
            }
            else if(string(ent->d_name)=="Capture")
            {
                inSubDirName = inFileName + "Capture/CrossSection/";
                outSubDirName = outFileName + "Capture/CrossSection/";

                #if Timer>=1
                    ConvertDirect(inSubDirName, outSubDirName, prevTemp, newTemp, ascii, log, logFile, totalNumFiles, index, totalFileSize2, sumFileSize2, sumDuration, fileSize2List);
                #else
                    ConvertDirect(inSubDirName, outSubDirName, prevTemp, newTemp, ascii, log, logFile);
                #endif
            }
            else
            {
                #if Timer>=1
                    start = std::clock();
                    ConvertFile(inFileName, outFileName, string(ent->d_name), prevTemp, newTemp, ascii, log, logFile);
                    duration = double(std::clock()-start)/CLOCKS_PER_SEC;
                    printProgress(index, fileSize2List, totalFileSize2, totalNumFiles, duration, sumFileSize2, sumDuration);
                    index++;
                #else
                    ConvertFile(inFileName, outFileName, string(ent->d_name), prevTemp, newTemp, ascii, log, logFile);
                #endif
            }
          }
          closedir(dir);
        }
        else
        {
            cout << "Error: Could not open the given directory" << endl;
            return 1;
        }

    }
    else
    {
        ConvertFile(inFileName, outFileName, prevTemp, newTemp, ascii, log, logFile);
    }


    delete CLHEPRand;
    delete theEngine;
    elementNames.ClearStore();

    return result;

}

#if Timer>=1
void GetfileSize2List(string inFileName, std::vector<int> &fileSize2List, int &totalNumFiles, long &totalFileSize2)
{
    string filename, inSubDirName;
    DIR *dir;
    struct dirent *ent;
    if ((dir = opendir (inFileName.c_str())) != NULL)
    {
      /* print all the files and directories within directory */
      while ((ent = readdir (dir)) != NULL)
      {
        if(string(ent->d_name)=="Elastic")
        {
            inSubDirName = inFileName + "Elastic/CrossSection/";
            GetDirectoryFileSize(inSubDirName, fileSize2List, totalNumFiles, totalFileSize2);
        }
        else if(string(ent->d_name)=="Inelastic")
        {
            inSubDirName = inFileName + "Inelastic/CrossSection/";
            GetDirectoryFileSize(inSubDirName, fileSize2List, totalNumFiles, totalFileSize2);
        }
        else if(string(ent->d_name)=="Fission")
        {
            inSubDirName = inFileName + "Fission/CrossSection/";
            GetDirectoryFileSize(inSubDirName, fileSize2List, totalNumFiles, totalFileSize2);
        }
        else if(string(ent->d_name)=="Capture")
        {
            inSubDirName = inFileName + "Capture/CrossSection/";
            GetDirectoryFileSize(inSubDirName, fileSize2List, totalNumFiles, totalFileSize2);
        }
        else
        {
            if(isApplicable(ent->d_name))
            {
                filename = inFileName + ent->d_name;
                fileSize2List[totalNumFiles]=GetFileSize2(filename);
                totalFileSize2+=fileSize2List[totalNumFiles];
                totalNumFiles++;
            }
        }
      }
      closedir(dir);
    }
    else
    {
        cout << "Error: Could not open the given directory" << endl;
    }
}

void GetDirectoryFileSize(string inDirName, std::vector<int> &fileSize2List, int &index, long &totalFileSize2)
{
    DIR *dir;
    struct dirent *ent;
    string fileName;

    if ((dir = opendir (inDirName.c_str())) != NULL)
    {
      /* print all the files and directories within directory */
      while ((ent = readdir (dir)) != NULL)
      {

        fileName = ent->d_name;
        if(isApplicable(fileName))
        {
            fileName = inDirName + fileName;
            fileSize2List[index]=GetFileSize2(fileName);
            totalFileSize2+=fileSize2List[index];
            index++;
        }

      }
      closedir(dir);
    }
    else
    {
        cout << "### Error: Could not open directory" << inDirName << " ###" << endl;
    }
}

bool isApplicable(string fileName)
{
    int Z=-1, A=-1;
    ExtractZA(fileName, Z, A);
    if((Z==-1)||(A==-1))
        return false;
    else
        return true;
}

int GetFileSize2(string fileName)
{
    std::ifstream in( fileName.c_str() , std::ios::binary | std::ios::ate );
    int file_size2 =0;
    if(in.good())
    {
        file_size2 = (in.tellg());
        file_size2*=file_size2;
        in.seekg( 0 , std::ios::beg );
        in.close();
    }
    else
    {
        cout << "### Error: failed to open file " << fileName << " to determine the filesize ###" << endl;
    }

    return file_size2;

}

void ConvertDirect(string inDirName, string outDirName, double prevTemp, double newTemp, bool ascii, bool log, std::ofstream* logFile,
    int totalNumFiles, int &index, long totalFileSize2, long &sumFileSize2, double &sumDuration, std::vector<int> &fileSize2List)
{
    DIR *dir;
    struct dirent *ent;
    string fileName;
    std::clock_t start;
    double duration;
    int temp=1;

    if ((dir = opendir (inDirName.c_str())) != NULL)
    {
      /* print all the files and directories within directory */
      while ((ent = readdir (dir)) != NULL)
      {
            if(!(DirectoryExists(outDirName.c_str())))
                temp = system( ("mkdir -p "+outDirName).c_str());

            if(temp)
            {
                start = std::clock();
                fileName = ent->d_name;
                ConvertFile(inDirName, outDirName, fileName, prevTemp, newTemp, ascii, log, logFile);
                duration = double(std::clock()-start)/CLOCKS_PER_SEC;
                printProgress(index, fileSize2List, totalFileSize2, totalNumFiles, duration, sumFileSize2, sumDuration);
                index++;
            }
            else
            {
                cout << "### Error: Could not create missing directory" << inDirName << " ###" << endl;
            }
            temp=1;

      }
      closedir(dir);
    }
    else
    {
        cout << "### Error: Could not open directory" << inDirName << " ###" << endl;
    }
}

void printProgress(int index, std::vector<int> &fileSize2List, long totalFileSize2, int totalNumFiles, double duration, long sumFileSize2, double sumDuration)
{
    sumFileSize2+=fileSize2List[index];
    sumDuration+=duration;
    double timeLeft = sumDuration/sumFileSize2*(totalFileSize2-sumFileSize2);
    cout << "Files Converted (" << index+1 << "/" << totalNumFiles << "), Time Remaining " << timeLeft << "s\n"
        << "Progress [";
    for (int i=0; i<int(60*sumFileSize2/totalFileSize2); i++)
    {
        cout << "#";
    }
    for (int i=0; i<int(60-60*sumFileSize2/totalFileSize2); i++)
    {
        cout << "_";
    }
    cout << "]" << endl;
}

#else

void ConvertDirect(string inDirName, string outDirName, double prevTemp, double newTemp, bool ascii, bool log, std::ofstream* logFile)
{
    DIR *dir;
    struct dirent *ent;
    string fileName;
    int temp=1;

    if ((dir = opendir (inDirName.c_str())) != NULL)
    {
      /* print all the files and directories within directory */
      while ((ent = readdir (dir)) != NULL)
      {
            if(!(DirectoryExists(outDirName.c_str())))
                temp = system( ("mkdir -p "+outDirName).c_str());

            if(temp)
            {
                fileName = ent->d_name;
                ConvertFile(inDirName, outDirName, fileName, prevTemp, newTemp, ascii, log, logFile);
            }
            else
            {
                if(log)
                    (*logFile) << "Could not create missing directory" << inDirName << endl;
            }
            temp=1;

      }
      closedir(dir);
    }
    else
    {
        if(log)
            (*logFile) << "Could not open directory" << inDirName << endl;
    }
}

#endif

void ConvertFile(string inDirName, string outDirName, string fileName, double prevTemp, double newTemp, bool ascii, bool log, std::ofstream* logFile)
{
    if(log)
    {
        logFile->fill('-');
        (*logFile) << std::setw(84) << std::left << "Start Broadening "+inDirName+fileName << endl;
    }

    int Z=-1, A=-1, result=0;
    double isoMassC2;
    ExtractZA(fileName, Z, A, log, logFile);
    isoMassC2 = GetNuclearMass(A,Z, log, logFile);
    inDirName += fileName;
    outDirName += fileName;
    if(isoMassC2)
        result = DoppBroad(inDirName, outDirName, prevTemp, newTemp, isoMassC2, ascii, log, logFile);
    else
    {
        if(log)
            (*logFile) << "### Skipped File " << inDirName << " ###" << endl;
        result = 1;
    }
    if(result&&log)
        (*logFile) << "### Error: Broadening File " << inDirName << " ###" << endl;

    if(log)
    {
        logFile->fill('-');
        (*logFile) << std::setw(84) << std::left << "End Broadening "+inDirName << endl;
    }
}

void ConvertFile(string inFileName, string outFileName, double prevTemp, double newTemp, bool ascii, bool log, std::ofstream* logFile)
{

    if(log)
    {
        logFile->fill('-');
        (*logFile) << std::setw(84) << std::left << "Start Broadening "+inFileName << endl;
    }

    int Z=-1, A=-1, result=0;
    double isoMassC2;

    size_t pos = inFileName.find_last_of('/');
    if(pos == std::string::npos)
        pos=0;
    else
        pos++;

    string fileName = inFileName.substr(pos, std::string::npos);
    ExtractZA(fileName, Z, A, log, logFile);
    isoMassC2 = GetNuclearMass(A,Z, log, logFile);
    if(isoMassC2)
        result = DoppBroad(inFileName, outFileName, prevTemp, newTemp, isoMassC2, ascii, log, logFile);
    else
    {
        if(log)
            (*logFile) << "### Skipped File " << inFileName << " ###" << endl;
        result = 1;
    }
    if(result&&log)
        (*logFile) << "### Error: Broadening File " << inFileName << " ###" << endl;

    if(log)
    {
        logFile->fill('-');
        (*logFile) << std::setw(84) << std::left << "End Broadening "+inFileName << endl;
    }
}

bool DirectoryExists( const char* pzPath )
{
    if ( pzPath == NULL) return false;

    DIR *pDir;
    bool bExists = false;

    pDir = opendir (pzPath);

    if (pDir != NULL)
    {
        bExists = true;
        closedir (pDir);
    }

    return bExists;
}
 #if Timer>=1
void ExtractZA(string fileName, int &Z, int &A)
{
        std::size_t startPos=0;
        stringstream ss;
        ElementNames* elementNames;
        while(startPos!=fileName.length() && (fileName[startPos]<'0' || fileName[startPos]>'9'))
            startPos++;

        if(startPos==fileName.length())
        {
            Z=A=-1;
        }
        else
        {
        ////
            std::size_t found1 = fileName.find_first_of('_', startPos);
            if (found1==std::string::npos)
            {
                Z=A=-1;
            }
            else
            {
                std::size_t found2 = fileName.find_first_of('_', found1+1);
                if (found2==std::string::npos)
                {
                    Z=A=-1;
                }
                else
                {

                    ss.str(fileName.substr(startPos, found1));
                    ss >> Z;
                    ss.str("");
                    ss.clear();
                    if(((found2-found1-1) > 2) && (fileName[found2-2] == 'm'))
                        ss.str(fileName.substr(found1+1, found2-found1-3));
                    else
                        ss.str(fileName.substr(found1+1, found2-found1-1));
                    ss >> A;
                    ss.str("");
                    ss.clear();
                    ss.str(fileName.substr(found2+1));
                    if (!(elementNames->CheckName(ss.str(), Z)))
                    {
                        Z=A=-1;
                    }
                    ss.str("");
                    ss.clear();
                }

            }

        }
}
#endif
void ExtractZA(string fileName, int &Z, int &A, bool log, std::ofstream* logFile)
{
        std::size_t startPos=0;
        stringstream ss;
        ElementNames* elementNames;
        while(startPos!=fileName.length() && (fileName[startPos]<'0' || fileName[startPos]>'9'))
            startPos++;

        if(startPos==fileName.length())
        {
            if(log)
                (*logFile) << "### File Name Does Not Contian a Z or an A Value " << fileName << " is Invalid for Broadening ###" << endl;
            Z=A=-1;
        }
        else
        {
        ////
            std::size_t found1 = fileName.find_first_of('_', startPos);
            if (found1==std::string::npos)
            {
                if(log)
                    (*logFile) << "### File Name Does Not Contian a '_', two are needed, one to seperate the Z and A value, \
                and one to seperate the A and the Element name " << fileName << " is Invalid for Broadening ###" << endl;
                Z=A=-1;
            }
            else
            {
                std::size_t found2 = fileName.find_first_of('_', found1+1);
                if (found2==std::string::npos)
                {
                    if(log)
                        (*logFile) << "### File Name Does Not Contian a second '_', two are needed, one to seperate the Z and A value, \
                    and one to seperate the A and the Element name " << fileName << " is Invalid for Broadening ###" << endl;
                    Z=A=-1;
                }
                else
                {

                    ss.str(fileName.substr(startPos, found1));
                    ss >> Z;
                    ss.str("");
                    ss.clear();
                    if(((found2-found1-1) > 2) && (fileName[found2-2] == 'm'))
                        ss.str(fileName.substr(found1+1, found2-found1-3));
                    else
                        ss.str(fileName.substr(found1+1, found2-found1-1));
                    ss >> A;
                    ss.str("");
                    ss.clear();
                    ss.str(fileName.substr(found2+1));
                    if (!(elementNames->CheckName(ss.str(), Z)))
                    {
                        if(log)
                            (*logFile) << "### " << fileName << " does not include the correct element name at the end ###" << endl;
                        Z=A=-1;
                    }
                    ss.str("");
                    ss.clear();
                }

            }

        }
}

double GetNuclearMass(double A, double Z, bool log, std::ofstream* logFile)
{
    double mass;
    if (Z == -1 || A==-1)
        return 0.;

    if (G4NucleiPropertiesTableAME03::IsInTable(Z,A))
    {
      // AME 03 table
      mass = G4NucleiPropertiesTableAME03::GetNuclearMass(Z,A);
    }
    else if (G4NucleiPropertiesTheoreticalTable::IsInTable(Z,A))
    {
      // Theoretical table
      mass = G4NucleiPropertiesTheoreticalTable::GetNuclearMass(Z,A);
    }
    else
    {
      mass = NuclearMass(double(A),double(Z), log, logFile);
    }
    return mass;
}

double NuclearMass(double A, double Z, bool log, std::ofstream* logFile)
{
  if (A < 1 || Z < 0 || Z > A) {

      if(log)
        (*logFile) << "NuclearMass: Wrong values for A = "
	     << A << " and Z = " << Z << endl;
    return 0.0;
  }

  double mass = AtomicMass(A,Z);
  // atomic mass is converted to nuclear mass according formula in  AME03
  mass -= Z*electron_mass_c2;
  mass += ( 14.4381*std::pow ( Z , 2.39 ) + 1.55468*1e-6*std::pow ( Z , 5.35 ) )*eV;

  return mass;
}

double  AtomicMass(double A, double Z)
{
  const double hydrogen_mass_excess = G4NucleiPropertiesTableAME03::GetMassExcess(1,1);
  const double neutron_mass_excess =  G4NucleiPropertiesTableAME03::GetMassExcess(0,1);

  double mass =
      (A-Z)*neutron_mass_excess + Z*hydrogen_mass_excess - BindingEnergy(A,Z) + A*amu_c2;

  return mass;
}

double  BindingEnergy(double A, double Z)
{
  //
  // Weitzsaecker's Mass formula
  //
  int Npairing = int(A-Z)%2;                  // pairing
  int Zpairing = int(Z)%2;
  double binding =
      - 15.67*A                           // nuclear volume
      + 17.23*std::pow(A,2./3.)                // surface energy
      + 93.15*((A/2.-Z)*(A/2.-Z))/A       // asymmetry
      + 0.6984523*Z*Z*std::pow(A,-1./3.);      // coulomb
  if( Npairing == Zpairing ) binding += (Npairing+Zpairing-1) * 12.0 / std::sqrt(A);  // pairing

  return -binding;
}

int DoppBroad(string inFileName, string outFileName, double prevTemp, double newTemp, double isoMassC2, bool ascii, bool log, std::ofstream* logFile)
{

    stringstream ss;
    #if Timer>=2
        std::clock_t loopStart;
        loopStart = std::clock();
        double loopdur;
    #endif
    #if Timer>=3
        std::clock_t start;
        double duration;
    #endif

    GetDataStream(inFileName, ss, log, logFile);

	// Find the temperature difference between the temperature the cross
	// section was evaluated at versus the temperature of the material
	double tempDiff = newTemp - prevTemp;

	if (tempDiff<=0)
	{
        ss.str("");
        if(log)
            (*logFile) << "the new cross-section temperature is the same as the previous cross-section temperature";
        return 1;
    }


    int tableSize=0;
    ss >> tableSize;

    // Physics Vector
    int vType=0;
    ss >>  vType;

    // contents
    int siz=0;
    ss >> siz;
    if (ss.fail())
    { ss.str(""); return 1; }
    if (siz<=0)
    {
        if(log)
            (*logFile) << " Invalid vector size: " << siz << endl;
        ss.str("");
        return 1;
    }

    double tempEn, tempCS;
    double* prevEnVec = new double[siz];
    double* prevCSVec = new double[siz];

    for(int i = 0; i < siz ; i++)
    {
        tempEn = 0.;
        tempCS = 0.;
        ss >> tempEn >> tempCS;

        if (ss.fail())
        { ss.str(""); return 1; }

        prevEnVec[i] = tempEn;
        prevCSVec[i] = tempCS;
    }

     ss.str("");
     ss.clear();

     if(ascii)
     {
        ss.fill(' ');
        ss << std::setw(14) << std::right << tableSize << "\n" << std::setw(14) << vType << "\n" << std::setw(14) << siz << "\n";
        ss.precision(6);
        ss.setf(std::ios::scientific);
     }

     else
        ss << tableSize << "\n" << vType << "\n" << siz << "\n";


    // Declarations for averaging loop
    int initNumIL = int(std::max(100.,tempDiff/6.)), numIL/*, numMaxOL=3000*initNumIL*/;
    double nKEnerTrans;
    double aXsection;             // random xsec value
    double result;            // final averaged cross section
    double buffer;            // temporary value for previous average
    int counter;              // counts the number averaged xsec sets
    double nPx, nPy=0, nPz=0, isoPx, isoPy, isoPz, isoP2, isoEn /*isoGamma*/, a, x, y, z, p2, nEn, nVMag, sigma=std::sqrt(k_Boltzmann*tempDiff*isoMassC2);
    double isoMassC2Sq = isoMassC2*isoMassC2, nMassC2Sq = nMassC2*nMassC2, invIsoMassC2 = 1/isoMassC2, invNMassC2 = 1/nMassC2, invNVMag=0;
    #if Timer>=3
        double avgFindTime=0, /*prevAvgFindTime=0,*/ avgMathTime=0/*, prevAvgMathTime=0*/;
        double avgNumLoops=0/*, tempMathTime, tempFindTime*/;
    #endif

     #if Timer>=2
            cout << "Progress of " << inFileName << ":\n [";
            cout.flush();
            int progressCount=1;
    #endif

    for(int i = 0; i < siz ; i++)
    {
        //prevAvgMathTime=avgMathTime;
        //prevAvgFindTime=avgFindTime;

        #if Timer>=3
            start = std::clock();
        #endif
        // Find the Doppler broadened cross section
        result = 0;
        buffer = 0;
        counter = 0;
        numIL = initNumIL;

        nEn = (prevEnVec[i]/1000000)+nMassC2;
        nPx = sqrt(pow(nEn,2)-nMassC2Sq);
        nVMag = nPx*invNMassC2;
        invNVMag = 1/nVMag;

        while(counter == 0 || (/*counter < numMaxOL && */std::abs(buffer-result/std::max(1,counter)) > 0.003*buffer))
        {
            if(counter)
                buffer = result/counter;

            while(counter < numIL)
            {
                isoPx = (CLHEP::RandGaussQ::shoot())*sigma;
                isoPy = (CLHEP::RandGaussQ::shoot())*sigma;
                isoPz = (CLHEP::RandGaussQ::shoot())*sigma;
                isoP2 = isoPx*isoPx+isoPy*isoPy+isoPz*isoPz;

                isoEn = std::sqrt( isoMassC2Sq+isoP2 );

                // this is only needed for the velocity correction when the isoKE is significant compared to the isoMass wich it is not unless temp >10^6K
                //isoGamma = isoEn/(isoMassC2);

            // Find the neutron in the rest frame of the nucleus
                a = ( (nPx*isoPx)/(isoEn+isoMassC2) - nEn )*invIsoMassC2;
                x = nPx+a*isoPx;
                y = nPy+a*isoPy;
                z = nPz+a*isoPz;
                p2 = x*x+y*y+z*z;
                nKEnerTrans = std::sqrt(nMassC2Sq+p2)-nMassC2;

                // Find the cross-section of the element at the energy of the
                // neutron IN THE REST FRAME OF THE NUCLEUS
                #if Timer>=3
                    avgMathTime+=double(std::clock() - start);
                    start = std::clock();
                #endif

                aXsection = findCS(nKEnerTrans, prevEnVec, prevCSVec, siz);

                #if Timer>=3
                    avgFindTime+=double(std::clock() - start);
                    start = std::clock();
                #endif

                // Velocity correction
                aXsection *= sqrt(pow((nVMag-isoPx*invIsoMassC2),2)+pow((isoPy*invIsoMassC2),2)+pow((isoPz*invIsoMassC2),2))*invNVMag;

                // Add cross section to results
                result += aXsection;

                #if Timer>=3
                    avgMathTime+=double(std::clock() - start);
                    start = std::clock();
                #endif

                counter++;
            }

            numIL +=numIL;

        }

        //tempFindTime = double((avgFindTime-prevAvgFindTime)/counter)/CLOCKS_PER_SEC;
        //tempMathTime = double((avgMathTime-prevAvgMathTime)/counter)/CLOCKS_PER_SEC;
        #if Timer>=3
            avgNumLoops += counter;
        #endif

        if (ascii)
        {
            if(!((i)%3) && i!=0 )
                ss << "\n";
            ss << std::setw(14) << prevEnVec[i] << std::setw(14) << (result/counter);
        }
        else
            ss << prevEnVec[i] << ' ' << (result/counter) << ' ';

        #if Timer>=2
            if(progressCount==int(60*i/(siz-1)))
            {
                cout << "#";
                cout.flush();
                progressCount++;
            }
            if(i==siz-1)
            {
                loopdur = double(std::clock()-loopStart)/CLOCKS_PER_SEC;
                cout << "] Time Taken " << loopdur << "s" << endl;
            }


        #endif
    }

    #if Timer>=3
        avgFindTime /= avgNumLoops;
        avgMathTime /= avgNumLoops;
        avgNumLoops /= siz;

        double findperc, mathperc;
        findperc = avgFindTime*100/(avgFindTime+avgMathTime);
        mathperc = avgMathTime*100/(avgFindTime+avgMathTime);
    #endif

    ss << '\n';

    SetDataStream( outFileName , ss, ascii, log, logFile );

    #if Timer>=3
        if(log)
        {
            duration = (double(std::clock() - timebegin))/CLOCKS_PER_SEC;
            (*logFile) << "\n \n \n ### " << outFileName << " took " << duration <<"s to broaden and it had " << siz << " data points ###" << endl;
            duration = (double(avgMathTime))/CLOCKS_PER_SEC;
            (*logFile) << "### the average time it took to do the math per loop was " << duration << "s, this took " << mathperc << "% of the time ###" << endl;
            duration = (double(avgFindTime))/CLOCKS_PER_SEC;
            (*logFile) << "### the average time it took to find the CS per loop was " << duration << "s, this took " << findperc << "% of the time ###" << endl;
            (*logFile) << "### the average number of loops per data point was " << avgNumLoops << "### \n \n \n" << endl;
        }
    #endif
    ss.str("");
    delete [] prevCSVec;
    delete [] prevEnVec;
    return 0 ;
}

double findCS(double nKEnerTrans, const double* prevEnVec, const double* prevCSVec, int vecSize)
{
    int index=0;
    bool Back, Front;
    int step=int(vecSize/10);

    nKEnerTrans=nKEnerTrans*1000000.;

    if (nKEnerTrans<=prevEnVec[0])
    {
        return  prevCSVec[0];
    }
    else if (nKEnerTrans>=prevEnVec[vecSize-1])
    {
        return  prevCSVec[vecSize-1];
    }

    else
    {
        index = int((vecSize-1)*(nKEnerTrans-prevEnVec[0])
                    /(prevEnVec[vecSize-1]-prevEnVec[0]));

        while(step>1)
        {
            if(index<0)
                index=0;

            if((index+step)>int(vecSize-1))
                index=vecSize-1-step;

            (prevEnVec[index]<=nKEnerTrans)? Back=true : Back=false;
            (prevEnVec[index+step]>nKEnerTrans)? Front=true : Front=false;

            while (!Back)
            {
                index -= step;
                if(index<0)
                    index=0;
                (prevEnVec[index]<=nKEnerTrans)? Back=true : Back=false;
            }

            while (!Front)
            {
                index += step;
                if((index+step)>int(vecSize-1))
                    index=vecSize-1-step;
                (prevEnVec[index+step]>nKEnerTrans)? Front=true : Front=false;
            }

            if(prevEnVec[index]==nKEnerTrans)
                break;

            index = int((step)*(nKEnerTrans-prevEnVec[index])
                    /(prevEnVec[index+step]-prevEnVec[index])+index);

            step/=10;
        }

       while(prevEnVec[index]<nKEnerTrans)
       {
            index++;
       }

       while(prevEnVec[index]>nKEnerTrans)
       {
            index--;
       }

       if (prevEnVec[index]!=nKEnerTrans)
       {
            return ((prevCSVec[index+1]-prevCSVec[index])*(nKEnerTrans-prevEnVec[index])/(prevEnVec[index+1]-prevEnVec[index])+prevCSVec[index]);
       }

       return  prevCSVec[index];

    }

}

void GetDataStream( string filename , std::stringstream& ss, bool log, std::ofstream* logFile)
{
   string* data=NULL;
   std::ifstream* in=NULL;
   //string compfilename(filename);

   if(filename.substr((filename.length()-2),2)==".z")
   {
        in = new std::ifstream ( filename.c_str() , std::ios::binary | std::ios::ate );
   }

   if ( in!=NULL && in->good() )
   {
// Use the compressed file
      uLongf file_size = (uLongf)(in->tellg());
      in->seekg( 0 , std::ios::beg );
      Bytef* compdata = new Bytef[ file_size ];

      while ( *in )
      {
         in->read( (char*)compdata , file_size );
      }

      uLongf complen = (uLongf) ( file_size*4 );
      Bytef* uncompdata = new Bytef[complen];

      while ( Z_OK != uncompress ( uncompdata , &complen , compdata , file_size ) )
      {
         delete[] uncompdata;
         complen *= 2;
         uncompdata = new Bytef[complen];
      }
      delete [] compdata;
      //                                 Now "complen" has uncomplessed size
      data = new string ( (char*)uncompdata , (long)complen );
      delete [] uncompdata;
   }
   else {
// Use regular text file
      std::ifstream thefData( filename.c_str() , std::ios::in | std::ios::ate );
      if ( thefData.good() )
      {
         int file_size = thefData.tellg();
         thefData.seekg( 0 , std::ios::beg );
         char* filedata = new char[ file_size ];
         while ( thefData )
         {
            thefData.read( filedata , file_size );
         }
         thefData.close();
         data = new string ( filedata , file_size );
         delete [] filedata;
      }
      else
      {
// found no data file
//                 set error bit to the stream
         ss.setstate( std::ios::badbit );
         if(log)
            (*logFile) << endl << "### failed to open ascii file " << filename << " ###" << endl;
      }
   }
   if (data != NULL)
   {
        ss.str(*data);
        if(data->back()!='\n')
            ss << "\n";
        ss.seekg( 0 , std::ios::beg );
    }

   if(in!=NULL)
   {
        in->close();
        delete in;
   }

   delete data;
}


void SetDataStream( string filename , std::stringstream& ss, bool ascii, bool log, std::ofstream* logFile )
{
    //bool cond=true;
   if (!ascii)
   {
        string compfilename(filename);

        if(compfilename.back()!='z')
            compfilename += ".z";

       std::ofstream* out = new std::ofstream ( compfilename.c_str() , std::ios::binary | std::ios::trunc);
       if ( ss.good() )
       {
       //
    // Create the compressed file
          ss.seekg( 0 , std::ios::end );
          uLongf file_size = (uLongf)(ss.tellg());
          ss.seekg( 0 , std::ios::beg );
          Bytef* uncompdata = new Bytef[ file_size ];

          while ( ss ) {
              ss.read( (char*)uncompdata , file_size );
          }

          uLongf complen = compressBound(file_size);

          Bytef* compdata = new Bytef[complen];

          if ( Z_OK == compress ( compdata , &complen , uncompdata , file_size ) )
          {
            out->write((char*)compdata, (long)complen);
            if (out->fail())
            {
                if(log)
                    (*logFile) << endl << "writing the compressed data to the output file " << compfilename << " failed" << endl
                    << " may not have permission to delete an older version of the file" << endl;
            }
          }
          else
          {
            if(log)
                (*logFile) << endl << "compressing the data failed" << endl;
          }

          delete [] uncompdata;
          delete [] compdata;
       }
       else
       {
            if(log)
                (*logFile) << endl << "### failed to write to binary file ###" << endl;
       }

       out->close(); delete out;
   }
   else
   {
// Use regular text file
    string compfilename(filename);

    if(compfilename.substr((compfilename.length()-2),2)==".z")
    {
        compfilename.pop_back();
        compfilename.pop_back();
    }

      std::ofstream out( compfilename.c_str() , std::ios::out | std::ios::trunc );
      if ( ss.good() )
      {
         ss.seekg( 0 , std::ios::end );
         int file_size = ss.tellg();
         ss.seekg( 0 , std::ios::beg );
         char* filedata = new char[ file_size ];
         while ( ss ) {
            ss.read( filedata , file_size );
            if(!file_size)
            {
                if(log)
                    (*logFile) << "\n #### Error the size of the stringstream is invalid ###" << endl;
                break;
            }
         }
         out.write(filedata, file_size);
         if (out.fail())
        {
            if(log)
                (*logFile) << endl << "writing the ascii data to the output file " << compfilename << " failed" << endl
                 << " may not have permission to delete an older version of the file" << endl;
        }
         out.close();
         delete [] filedata;
      }
      else
      {
// found no data file
//                 set error bit to the stream
         ss.setstate( std::ios::badbit );

         if(log)
            (*logFile) << endl << "### failed to write to ascii file " << compfilename << " ###" << endl;
      }
   }
   ss.str("");
}
