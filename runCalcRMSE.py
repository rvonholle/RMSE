import CalculateRMSError_faster
import sys
import os

if ( len(sys.argv) != 3 ):
   print "Usage:", sys.argv[0], "[Reference Plot File] [Plot File]"
   exit(1)

success = CalculateRMSError_faster.CalculateRMSError_RCS(sys.argv[1],sys.argv[2])
