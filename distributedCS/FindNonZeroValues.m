function NumDomEigVals= FindNonZeroValues(array,threshold)

NumDomEigVals = sum(abs(array > threshold));