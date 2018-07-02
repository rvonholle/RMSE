#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <algorithm>

#define SWITCH_ACCURACY_THRESHOLD 10e-4
#define LEGO_ACCURACY_THRESHOLD 10e-3
#define ACCURACY_ERROR_SWITCH 122
#define ACCURACY_ERROR_LEGO 123
#define FILE_ERROR 125

using namespace std;

////////////////////////////////////////////////////////////////////////////////
// Angle = Frequency + Incident Angle + Scatter Angle
////////////////////////////////////////////////////////////////////////////////
struct Angle {
// Data Members
   double freq, theta, phi, stheta, sphi;

// Functions
   bool operator< (Angle const& rhs) const
   {
      if(freq   < rhs.freq  ) return true;
      if(freq   > rhs.freq  ) return false;

      if(phi    < rhs.phi   ) return true;
      if(phi    > rhs.phi   ) return false;

      if(theta  < rhs.theta ) return true;
      if(theta  > rhs.theta ) return false;

      if(sphi   < rhs.sphi  ) return true;
      if(sphi   > rhs.sphi  ) return false;

      if(stheta < rhs.stheta) return true;
      if(stheta > rhs.stheta) return false;

      return false;
   }
   double getFreq()
   {
      return freq;
   }
};

////////////////////////////////////////////////////////////////////////////////
// Field = Phi-Phi and Theta-Theta Phase & Magnitude (in dBsm)
////////////////////////////////////////////////////////////////////////////////
struct Field {
// Data Members
   double pp_dbsm, pp_phase, tt_dbsm, tt_phase;

// Functions
   complex<double> ppField() const
   {
      return field(pp_dbsm, pp_phase);
   }
   complex<double> ttField() const
   {
      return field(tt_dbsm, tt_phase);
   }
   complex<double> field(const double dbsm, const double phase) const
   {
      double pi = acos(-1);

      double RCS = pow(10.0, dbsm/10.0);
      double EMag = sqrt(RCS/(4.0*pi));
      complex<double> Field = EMag * exp(complex<double>(0.0,phase*pi/180.0));
      return Field;
   }
};

////////////////////////////////////////////////////////////////////////////////
// AngleFieldPair = Angle + Associated Field
// Angle = Frequency + Incident Angle + Scatter Angle
// Field = Phi-Phi and Theta-Theta Phase & Magnitude (in dBsm)
////////////////////////////////////////////////////////////////////////////////
class AngleFieldPair {
public:
   Angle angle;
   Field field;
private:
   friend istream& operator>>(istream& in, AngleFieldPair& pair);
};

istream& operator>>(istream &in, AngleFieldPair& pair)
{
   vector<double> floats;
   while(true) {
      if(in.eof()) {
         break;
      }
      if(in.peek() == '\n') {
         in.get();
         break;
      }
      double tmp;
      in >> tmp;
      while(!in.eof() && in.peek() == ' ') {
         // Consume extra whitespace
         in.get();
      }
      if(in.fail()) {
         break;
      }
      floats.push_back(tmp);
   }
   switch (floats.size()) {
      case 7: // Mono-static case
         pair.angle.freq   = floats[0];
         pair.angle.theta  = floats[1];
         pair.angle.phi    = floats[2];
         pair.angle.stheta = floats[1];
         pair.angle.sphi   = floats[2];
         pair.field.pp_dbsm  = floats[3];
         pair.field.pp_phase = floats[4];
         pair.field.tt_dbsm  = floats[5];
         pair.field.tt_phase = floats[6];
         //for (int i = 0; i < 7; i++) {
           // cout << setprecision(10) << floats[i] << endl;
         //}
         break;
      case 9: // Bi-static case
         pair.angle.freq   = floats[0];
         pair.angle.theta  = floats[1];
         pair.angle.phi    = floats[2];
         pair.angle.stheta = floats[3];
         pair.angle.sphi   = floats[4];
         pair.field.pp_dbsm  = floats[5];
         pair.field.pp_phase = floats[6];
         pair.field.tt_dbsm  = floats[7];
         pair.field.tt_phase = floats[8];
         break;
      default:
         cout << "Warning: Wrong number of floats on line. Bad plot file?" << endl;
         break;
   }
}

////////////////////////////////////////////////////////////////////////////////
// ErrorInfoByFreq
////////////////////////////////////////////////////////////////////////////////
struct ErrorInfo {
   int numAnglesCompared, numAnglesNotCompared;
   double ppEFieldError2Sum, ttEFieldError2Sum;
   double ppEFieldRef2Sum, ttEFieldRef2Sum;
   ErrorInfo()
   {
      numAnglesCompared = 0;
      numAnglesNotCompared = 0;
      ppEFieldError2Sum = 0;
      ttEFieldError2Sum = 0;
      ppEFieldRef2Sum = 0;
      ttEFieldRef2Sum = 0;
   }
};

class ErrorInfoByFreq {
public:
   void newEFieldError(double freq, double ppEFieldError, double ttEFieldError,
                                   double ppEFieldRef, double ttEFieldRef)
   {
      freqInfo[freq].numAnglesCompared++;
      freqInfo[freq].ppEFieldError2Sum += pow(ppEFieldError, 2);
      freqInfo[freq].ttEFieldError2Sum += pow(ttEFieldError, 2);
      freqInfo[freq].ppEFieldRef2Sum += pow(ppEFieldRef, 2);
      freqInfo[freq].ttEFieldRef2Sum += pow(ttEFieldRef, 2);
   }
   void angleNotCompared(double freq)
   {
      freqInfo[freq].numAnglesNotCompared++;
   }
   void outputErrors()
   {
      map<double, ErrorInfo>::iterator it;
      for(it = freqInfo.begin(); it != freqInfo.end(); it++) {
         int numCompared = it->second.numAnglesCompared;
         int numNotCompared = it->second.numAnglesNotCompared;
         double ppRMSError = sqrt(it->second.ppEFieldError2Sum / numCompared);
         //cout << it->second.ppEFieldError2Sum << endl;
         double ttRMSError = sqrt(it->second.ttEFieldError2Sum / numCompared);
         double ppRMSRelError = sqrt(it->second.ppEFieldError2Sum / it->second.ppEFieldRef2Sum);
         double ttRMSRelError = sqrt(it->second.ttEFieldError2Sum / it->second.ttEFieldRef2Sum);
         if (numCompared != 0) {
            cout << setw(57) << ""
                 << "+------------------------------" << endl;
            cout << setprecision(5);
            cout << "    Phi-Phi Error (Rel. Error) = " << setw(10) << ppRMSError
                 << " (" << setw(10) << ppRMSRelError << ")"
                 << " | RMS Error at " << it->first << " GHz" << endl;
            cout << "Theta-Theta Error (Rel. Error) = " << setw(10) << ttRMSError
                 << " (" << setw(10) << ttRMSRelError << ")"
                 << " | " << 100.0 * numCompared/(numCompared+numNotCompared)
                 << "% inc/scatter coverage" <<  endl;
            //if (ppRMSRelError > 0.01 || ttRMSRelError > 0.01){
            //   exit(123);
            //}
         } else {
            cout << setw(57) << ""
                 << "+------------------------------" << endl;
            cout << setprecision(5);
            cout << "    Phi-Phi Error (Rel. Error) =        N/A (       N/A)"
                 << " | RMS Error at " << it->first << " GHz" << endl;
            cout << "Theta-Theta Error (Rel. Error) =        N/A (       N/A)"
                 << " | " << 100.0 * numCompared/(numCompared+numNotCompared)
                 << "% inc/scatter coverage" <<  endl;
            //exit(123);
         }
      }
      if(freqInfo.size() != 0)
         cout << setw(57) << ""
              << "+------------------------------" << endl;
   }
   int passFail()
   {
      map<double, ErrorInfo>::iterator it;
      double maxError = 0.0;
      for(it = freqInfo.begin(); it != freqInfo.end(); it++) {
         int numCompared = it->second.numAnglesCompared;
         int numNotCompared = it->second.numAnglesNotCompared;
         double ppRMSError = sqrt(it->second.ppEFieldError2Sum / numCompared);
         double ttRMSError = sqrt(it->second.ttEFieldError2Sum / numCompared);
         double ppRMSRelError = sqrt(it->second.ppEFieldError2Sum / it->second.ppEFieldRef2Sum);
         double ttRMSRelError = sqrt(it->second.ttEFieldError2Sum / it->second.ttEFieldRef2Sum);
         // Check for nothing to compare
         if (numCompared == 0)
            maxError = 2.0*LEGO_ACCURACY_THRESHOLD;
         // This checks for NaN
         if ((ppRMSRelError != ppRMSRelError) || (ttRMSRelError != ttRMSRelError))
            maxError = 2.0*LEGO_ACCURACY_THRESHOLD;
         else
            maxError = max(max(maxError,ttRMSRelError),ppRMSRelError);
      }
      // SWITCH and SWITCH-LEGO pass
      if (maxError < LEGO_ACCURACY_THRESHOLD) {
         return 0;
      }
      // For sure both SWITCH and SWITCH-LEGO failure
      if (maxError > SWITCH_ACCURACY_THRESHOLD) { 
         return ACCURACY_ERROR_SWITCH;
      }
      return ACCURACY_ERROR_LEGO;
   }
private:
   map<double, ErrorInfo> freqInfo;
};

////////////////////////////////////////////////////////////////////////////////
// Main Program
////////////////////////////////////////////////////////////////////////////////
extern "C" {
int calcRMS(int argc, char *argv[])
{
   if(argc != 3) {
      cout << "Usage: " << argv[0] << " [Reference Plot File] [Plot File]" << endl;
      exit(1);
   }
   fstream fileRef(argv[1], ios_base::in);
   if(!fileRef.is_open()) {
      cout << "Error opening reference plot file `" << argv[1] << "'\n";
      exit(FILE_ERROR);
   }
   fstream fileNew(argv[2], ios_base::in);
   if(!fileNew.is_open()) {
      cout << "Error opening plot file `" << argv[2] << "'\n";
      exit(FILE_ERROR);
   }
   string lineToSkip;
   Field valNew, valRef;
   AngleFieldPair keyValuePair;
   ErrorInfoByFreq errorInfo;

   unsigned numAnglesRef = 0;
   unsigned numAnglesNew = 0;

   // Read from reference file (populate map)
   map<Angle, Field> refMap;
   while(true) {
      if(fileRef.peek() == '#') {
         getline(fileRef, lineToSkip);
      } else {
         if(fileRef.eof()) break;
         fileRef >> keyValuePair;
         numAnglesRef++;
         if(fileRef.fail()) {
            cout << "Error while reading reference plot file.\n";
            exit(FILE_ERROR);
         }
         refMap[keyValuePair.angle] = keyValuePair.field;
      }
   }

   // Read from new plot file (calculate RMS)
   while(true) {
      if(fileNew.peek() == '#') {
         getline(fileNew, lineToSkip);
      } else {
         if(fileNew.eof()) break;
         fileNew >> keyValuePair;
         numAnglesNew++;
         if(fileNew.fail()) {
            cerr << "Error while reading new plot file.\n";
            exit(FILE_ERROR);
         }
         map<Angle, Field>::iterator it;
         it = refMap.find(keyValuePair.angle);
         if(it != refMap.end()) {
            valRef = it->second;
            valNew = keyValuePair.field;
            errorInfo.newEFieldError(keyValuePair.angle.getFreq(),
                                     abs(valNew.ppField() - valRef.ppField()),
                                     abs(valNew.ttField() - valRef.ttField()),
                                     abs(valRef.ppField()),
                                     abs(valRef.ttField()));
           // cout << "Item1 pp" << setprecision(10) << valNew.ppField() << endl;
           // cout << "Item2 pp" << setprecision(10) << valRef.ppField() << endl;
           // cout << "Item1 tt" << setprecision(10) << valNew.ttField() << endl;
           // cout << "Item2 tt" << setprecision(10) << valRef.ttField() << endl;
         } else {
            errorInfo.angleNotCompared(keyValuePair.angle.getFreq());
         }
      }
   }
   if(numAnglesRef == 0) {
      cout << "Error: No data found in reference plot file" << endl;
      exit(FILE_ERROR);
   }
   if(numAnglesNew == 0) {
      cout << "Error: No data found in new plot file" << endl;
      exit(FILE_ERROR);
   }
   errorInfo.outputErrors();

   return errorInfo.passFail();
}
}
 
int main(int argc, char *argv[])
{
   calcRMS(argc, argv);
}  
  
  
