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
   fileData = []
   for line in text:
      if ( len(line)==7 ):
         tmp = [line[0],line[1],line[2],\
                cmath.sqrt((10.0**line[3]/10.0)/(4.0*cmath.pi)) * \
                cmath.exp(complex(0.0,line[4]*cmath.pi/180)),\
                cmath.sqrt((10.0**line[5]/10.0)/(4.0*cmath.pi)) * \
                cmath.exp(complex(0.0,line[6]*cmath.pi/180))]
      elif ( len(line)==9 ):
         tmp = [line[0],line[1],line[2],\
                line[3],line[4], \
                cmath.sqrt((10.0**line[5]/10.0)/(4.0*cmath.pi)) * \
                cmath.exp(complex(0.0,line[6]*cmath.pi/180)),\
                cmath.sqrt((10.0**line[7]/10.0)/(4.0*cmath.pi)) * \
                cmath.exp(complex(0.0,line[8]*cmath.pi/180))]
      else:
         print(fName+' has incorrect number of entries on the following line:')
         print(line)
         break
      fileData.append(tmp)
   
   #Handle incorrect outpt
   return sorted(fileData,key=operator.itemgetter(0,1,2))
#
def AlignData(long_data,short_data):
   #Iterate through all the data and make sure ref and test have same data points
   new_long_data = []

   #Force nans to be ignored for comparison purposes
   for item in long_data:
      found = False
      for point in short_data:
         if ( len(item)==5 ):
            if ( FieldFlags(item,point,3) ):
               found = True
               new_long_data.append(item)
         elif( len(item)==7 ):
            if ( FieldFlags(item,point,5) ):
               found = True
               new_long_data.append(item)
   #Remove list duplicates
   new_long_data = [list(i) for i in set(tuple(i) for i in new_long_data)]
   return new_long_data
#
def FieldFlags(list1,list2,num):
   flag = True
   for ii in range(0,num):
      flag = (flag and (list1[ii]==list2[ii]))
   return flag
#
def ResortData(listOfLists):
   if ( len(listOfLists[0])==5 ):
      listOfLists = sorted(listOfLists, key=itemgetter(0,2,1))
   elif ( len(listOfLists[0])==7 ):
      listOfLists = sorted(listOfLists, key=itemgetter(0,2,1,4,3))
   return listOfLists
#
def CalcRMSError(data_test,data_ref):
   rmse=complex(0,0)
   if ( (len(data_test)!=0) and (len(data_ref)!=0) ):
      for item1,item2 in zip(data_test,data_ref):
         rmse+=(item1-item2)**2
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
         rmse_test+=(item1-item2)**2
         rmse_ref +=item2**2
      rmse=cmath.sqrt(rmse_test/rmse_ref)
   else:
      rmse=complex(float('inf'))
   return rmse
#
def PassFail(RMSE_pp,RMSE_tt):
   check = 'Failed'
   if ( (cmath.isnan(RMSE_pp)) or (cmath.isnan(RMSE_tt)) ):#Catch nan errors
      print 'NaN Encountered.'
      check = 'Failed'
   elif ( (abs(RMSE_pp) >= THRESHOLD) or (abs(RMSE_tt) >= THRESHOLD) ):# Catch if greather than threshold
      print 'At least one of the RMSE values is incorrect'
      print 'RMSE_pp=%7.5e'%RMSE_pp.real+'+%7.5e'%RMSE_pp.imag+'j'
      print 'RMSE_tt=%7.5e'%RMSE_tt.real+'+%7.5e'%RMSE_tt.imag+'j'
      check = 'Failed'
   elif ( (abs(RMSE_pp) < THRESHOLD) or (abs(RMSE_tt) < THRESHOLD) ):# Check values
      print 'RMSE_pp=%7.5e'%RMSE_pp.real+'+%7.5e'%RMSE_pp.imag+'j'
      print 'RMSE_tt=%7.5e'%RMSE_tt.real+'+%7.5e'%RMSE_tt.imag+'j'
      check = 'Passed'
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
   
   fileText  = ParsePlotFile(fname_ref)
   data_ref  = ProcessData(fileText)
   fileText  = ParsePlotFile(fname_test)
   data_test = ProcessData(fileText)

   print 'Finished parsing'

   if ( len(data_test)>=len(data_ref) ):
      data_test = AlignData(data_test,data_ref)
      data_test = ResortData(data_test)
   else:
      data_ref = AlignData(data_ref,data_test)
      data_ref = ResortData(data_ref)

   data_ref  = ResortData(data_ref)
   data_test = ResortData(data_test)

   print 'Finished sorting'

   for dataPt in data_ref:
      sigma_pp_R.append(dataPt[3])
      sigma_tt_R.append(dataPt[4])

   for dataPt in data_test:
      sigma_pp_T.append(dataPt[3])
      sigma_tt_T.append(dataPt[4])

   RMSE_pp = CalcRMSErrorRel(sigma_pp_T,sigma_pp_R)
   RMSE_tt = CalcRMSErrorRel(sigma_tt_T,sigma_tt_R)

   return PassFail(RMSE_pp,RMSE_tt)
