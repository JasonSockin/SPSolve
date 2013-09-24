# Functions to check if the matrix a contains any NaN or Inf values
#
# Created by: Jason Sockin
# Date: September 11, 2013

# import numpy package
from numpy import *

# Checks for any NaN in a
def existsNaN(a):
    aTemp = isnan(a)
    aRows, aColumns = aTemp.shape
    result = False
    for i in range(0,aRows):
        for j in range(0,aColumns):
            if aTemp[i,j] == True:
                result = True
                break
        if result == True:
            break
    return result


# Checks for any Inf in a
def existsInf(a):
    aTemp = isinf(a)
    aRows, aColumns = aTemp.shape
    result = False
    for i in range(0,aRows):
        for j in range(0,aColumns):
            if aTemp[i,j] == True:
                result = True
                break
        if result == True:
            break
    return result
    
    
