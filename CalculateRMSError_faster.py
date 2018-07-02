#!/user/bin/python

from collections import namedtuple
from operator import itemgetter
import cmath
import math

################################################################
# Named Tuples
################################################################
fileText = namedtuple('fileText','freq th ph sth sph ppdbsm ppPh ttdbsm ttPh')

LEGO_ACCURACY_THRESHOLD   = 1e-3
SWITCH_ACCURACY_THRESHOLD = 1e-4
THRESHOLD                 = LEGO_ACCURACY_THRESHOLD

################################################################
# Function Defitions
################################################################
def ParsePlotFile(fName):
   import string

   data = []

   #Read in the file
   fileObj = open(fName,'r')
   with open(fName,'r') as f:
      content = f.read()
   f.close()

#   print content

   #Remove all comment lines
   content = content.split('\n')
   content = [line for line in content if not ('#' in line)]
   content = filter(None,content)

   #break lines into constituent parts
   for line in content:
      tmp = filter(None,line.split(' '))
      if ( (len(tmp)!=7) and (len(tmp)!=9) ):
         print(fName+' has incorrect number of entries on the following line:')
         print(line)
         break
      for index in range(0,len(tmp)):
         tmp[index] = float(tmp[index])
      data.append(tmp)
   
   return data
#
def ProcessData(text):
   import operator
   fileData = {}
   for line in text:
      if ( len(line)==7 ):
         tempAngle = (line[0], line[1], line[2])
         tempField = (cmath.sqrt((10.0**(line[3]/10.0))/(4.0*cmath.pi)) * \
                      cmath.exp(complex(0.0,line[4]*cmath.pi/180)), \
                      cmath.sqrt((10.0**(line[5]/10.0))/(4.0*cmath.pi)) * \
                      cmath.exp(complex(0.0,line[6]*cmath.pi/180)))
      elif ( len(line)==9 ):
         tempAngle = (line[0], line[1], line[2], line[3], line[4])
         tempField = (cmath.sqrt((10.0**(line[5]/10.0))/(4.0*cmath.pi)) * \
                      cmath.exp(complex(0.0,line[6]*cmath.pi/180)), \
                      cmath.sqrt((10.0**(line[7]/10.0))/(4.0*cmath.pi)) * \
                      cmath.exp(complex(0.0,line[8]*cmath.pi/180)))
      else:
         print(fName+' has incorrect number of entries on the following line:')
         print(line)
         break
      fileData.update({tempAngle:tempField})
   
   #Handle incorrect outpt
   return fileData
#
def AlignData(data_a,data_b):
   for key in data_b:
      if not(data_a.has_key(key)):
         del data_a[key]
   return data_a
#
def findFreqs(dataSet):
   freqs = []
   for key in dataSet:
      freqs.append(key[0])
   freqs = set(freqs)
   freqs = list(freqs)
   freqs.sort()
   return freqs
#
def CalcRMSError(data_test,data_ref):
   rmse=complex(0,0)
   if ( (len(data_test)!=0) and (len(data_ref)!=0) ):
      for item1,item2 in zip(data_test,data_ref):
         rmse+=abs(item1-item2)**2
      rmse=cmath.sqrt(rmse/len(data_test))
   return rmse
#
# cf. SWITCH/utils/CalculateRMSError.cpp
def CalcRMSErrorRel(data_test,data_ref):
   rmse_test=complex(0,0)
   rmse_ref =complex(0,0)
   rmse     =complex(0,0)

   if ( (len(data_test)!=0) and (len(data_ref)!=0) ):
      for item1,item2 in zip(data_test,data_ref):
         rmse_test+=(abs(item1-item2))**2
         rmse_ref +=(abs(item2))**2
      rmse=cmath.sqrt(rmse_test/rmse_ref)
   else:
      rmse=complex(float('inf'))
   return rmse
#
def PassFail(RMSE_pp,RMSE_tt,RMSE_pp_rel,RMSE_tt_rel,freqs):
   check = 'Passed'
   i = 0
   #Catch nan errors
   for item in freqs:
      i = i + 1
      if ( (cmath.isnan(RMSE_pp[item])) or (cmath.isnan(RMSE_tt[item])) ):
         print 'NaN Encountered.'
         check = 'Failed'
      #Catch if greater than threshold
      elif ( (abs(RMSE_pp[item]) >= THRESHOLD) or (abs(RMSE_tt[item]) >= THRESHOLD) ):
         print 'At least one of the RMSE values is unexpectedly high'
         print '                                                           +------------------------------'
         print '    Phi-Phi Error (Rel. Error) = {0:.4e}'.format(RMSE_pp[item].real), '( {0:.4e}'.format(RMSE_pp_rel[item].real), ') | RMS Error at', item, 'GHz'
         print 'Theta-Theta Error (Rel. Error) = {0:.4e}'.format(RMSE_tt[item].real), '( {0:.4e}'.format(RMSE_tt_rel[item].real), ') | 100% inc/scatter coverage'
         if ( i == len(freqs) ):
            print '                                                           +------------------------------'
         check = 'Failed'
      #Check values
      elif ( (abs(RMSE_pp[item]) < THRESHOLD) or (abs(RMSE_tt[item]) < THRESHOLD) ):
         print '                                                           +------------------------------'
         print '    Phi-Phi Error (Rel. Error) = {0:.4e}'.format(RMSE_pp[item].real), '( {0:.4e}'.format(RMSE_pp_rel[item].real), ') | RMS Error at', item, 'GHz'
         print 'Theta-Theta Error (Rel. Error) = {0:.4e}'.format(RMSE_tt[item].real), '( {0:.4e}'.format(RMSE_tt_rel[item].real), ') | 100% inc/scatter coverage'
         if ( i == len(freqs) ):
            print '                                                           +------------------------------'
      else:
         print 'Unexpected Error. Review Log Files.'
         check = 'Failed'
   return check
#
########################################
# MAIN FUNCTION
########################################
def CalculateRMSError_RCS(fname_test,fname_ref):
   sigma_pp_R = []
   sigma_tt_R = []
   sigma_pp_T = []
   sigma_tt_T = []
   RMSE_pp = {}
   RMSE_tt = {}
   RMSE_pp_rel = {}
   RMSE_tt_rel = {}
   freqs = []
   
   fileText  = ParsePlotFile(fname_ref)
   data_ref  = ProcessData(fileText)
   fileText  = ParsePlotFile(fname_test)
   data_test = ProcessData(fileText)

#   print 'Finished parsing'

   data_test = AlignData(data_test,data_ref)
   data_ref = AlignData(data_ref,data_test)

   sorted(data_ref.iterkeys())
   sorted(data_test.iterkeys())

   freqs = findFreqs(data_ref)

#   print 'Finished sorting'

   for item in freqs:
      for key in data_ref:
         if ( key[0] == item ):
            sigma_pp_R.append(data_ref[key][0])
            sigma_tt_R.append(data_ref[key][1])
            sigma_pp_T.append(data_test[key][0])
            sigma_tt_T.append(data_test[key][1])
      RMSE_pp.update({item:CalcRMSError(sigma_pp_T,sigma_pp_R)})
      RMSE_tt.update({item:CalcRMSError(sigma_tt_T,sigma_tt_R)})
      RMSE_pp_rel.update({item:CalcRMSErrorRel(sigma_pp_T,sigma_pp_R)})
      RMSE_tt_rel.update({item:CalcRMSErrorRel(sigma_tt_T,sigma_tt_R)})
      sigma_pp_R[:] = []
      sigma_tt_R[:] = []
      sigma_pp_T[:] = []
      sigma_tt_T[:] = []

   return PassFail(RMSE_pp,RMSE_tt,RMSE_pp_rel,RMSE_tt_rel,freqs)
