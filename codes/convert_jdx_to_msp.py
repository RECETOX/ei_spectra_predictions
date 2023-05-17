import os
import pathlib

### get current file path
workPath = os.getcwd() + '/'
os.sys.path.append(workPath)

### should not contain CR/LF at the end of file
### multiple CR/LF are not handled or need to be stripped
metadataFilename = ('input/metadata.txt')
metadataHandle = workPath + metadataFilename
metadataText = pathlib.Path(metadataHandle).read_text()

# read spectrum file from plotms
spectrumFilename = ('accuratemass.jdx')
spectrumHandle = workPath + spectrumFilename
spectrumText = pathlib.Path(spectrumHandle).read_text()

# get number of peaks from accuratemass.jdx (element [5] in line 6)
numPeaks = spectrumText.split('\n')[5].split('=')[1]

# get mz-abd pairs and add all peaks and annotations
peakPairs = spectrumText.split('\n')[7:int(numPeaks)+7]
peakPairsText = '\n'.join(peakPairs)

# open file for writing new MSP (NIST format)
exportFilename = ('result-ms2.msp')
exportHandle = workPath + exportFilename
with open(exportHandle, 'w+') as f:
# write complete metadata block
f.write(metadataText)
# add number of peaks
numPeaksText = 'NumPeaks: ' + str(numPeaks) +'\n'
f.write(numPeaksText)
# write m/z - abundance pairs with annotation
f.write(peakPairsText)
# close the file
f.close()

# Program finished
print('\nResult MSP file converted: result-ms2.msp\n')