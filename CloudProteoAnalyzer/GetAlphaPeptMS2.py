import os
import platform
import sys
from alphapept.pyrawfilereader import RawFileReader

def MSSpectra(fname):
#def MSSpectra(fname):
    Hweight=1.0072
    raw_path = fname
    rawFile = RawFileReader(raw_path)
    length=rawFile.GetNumSpectra()
    msname=fname.replace('.raw','.ms2')
    with open(msname, "w") as f:
        for i in range(1,rawFile.GetNumSpectra()+1):
            scanNumber = i
            if rawFile.GetNumberOfMSOrdersFromScanNum(scanNumber)==0:
                continue
            masses, intens = rawFile.GetCentroidMassListFromScanNum(scanNumber)
            f.write("S\t" + str(scanNumber) + "\t" + str(scanNumber) + "\t" + str(rawFile.GetPrecursorMassForScanNum(scanNumber))+'\n')
            mz=rawFile.GetMS2MonoMzAndChargeFromScanNum(scanNumber)[0]
            z=rawFile.GetMS2MonoMzAndChargeFromScanNum(scanNumber)[1]
            mass=z*mz-(z-1)*Hweight
            f.write("Z\t" + str(z) + "\t" + str(mass)+'\n')
            for index in range(0, len(masses)):
                f.write(str(masses[index]) + " " + str(intens[index])+'\n')
            
            #bar.progress(i/length)
                    
    rawFile.Close()
    return msname

#MSSpectra("OSU_D10_FASP_Elite_03202014_01.raw")
