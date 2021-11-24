# Disfuse-FORS-Analasis-Functions
Python Module Created for My Masters in Physics Project at Exeter University in Diffuse Fibre Optic Reflectence Spectral Analysis

**ONLY FOR GOOGLE COLAB**

## Dependencies

- numpy
- scipy
- scipy.interpolate
- scipy.misc
- scipy.signal
- matplotlib.pyplot
- glob
- os
- shutil
- math
- google.colab.auth
- gspread
- oauth2client.client.GoogleCredentials

## Google Colab Headers to Import Function and Load Changes to Modual
```
### indicates different cell

!pip install GitPython
###
from google.colab import drive
drive.mount('/content/drive')
###
from git import Repo
import shutil
import sys
import importlib

#shutil.rmtree("/testing/")#uncomment for second run
Repo.clone_from("https://github.com/Haggi1181/Disfuse-FORS-Analasis-Functions.git", "/testing/")

sys.path.insert(1, "/testing/")

import MastersModual as mm
importlib.reload(mm)
```

## Functions
- PlotProcessedSpectra
    - Fuction to plot a already genorated FORS spectra
    - Args DataFilePath: File to plot, MatPlotLibColour: Colour for Mat Plot Lib plotts, HeaderSize: Number of lines to skip for header, title: lable at top of graph

- PlotFolderTxt
    - Fuction to plot all text files in a given dirrectory
    - Args Dir: Directory to plot, MatPlotLibColour: Colour for Mat Plot Lib plotts, HeaderSize: Number of lines to skip for header


- SpectralGenAndSave
    - Function to genorate and save a single FORS spectra
    - Args RawDataDir: Full file path for target data, DarkCountDir: Full file path for dark count, ReflectanceStanderdDir: Full file path for ref standard, SaveDir: Directory to save in, SafeFileName: name of resultant file auto adds .txt, HeaderSize: Number of lines to ohmit at top of reading files

- MatchingDatabaseCreator
    - Function to genorate a spectral database from spectra found in a single folder saving calculated spectra into one folder
    - args: Dir: Directory of data, DarkCountFileName: Full file path for dark count, ReflectenceStandardFileName: Full file path for ref standard, SaveDir: Directory to save in, HeaderSize: Number of lines to ohmit at top of reading files

- DiffuseRefelctencePlot
    - Function to plot FORS Spectra between 400 and 900nm
    - Args LegendName: Name of plot for legend, SpreadSheet: Name of google sheet found in gdrive containg data, DarkCountIndex: page index of dark count data found in spesifide google sheet starting at 0, ReflectenceStandardIndex: page index of reflectence standard data found in spesifide google sheet starting at 0, DataIndex: page index of raw data found in spesifide google sheet starting at 0, MatPlotLibColour: optional colour of plot

- SpectralGen
    - Function to Calculate FORS Spectra between 400 and 900nm
    - Args yRawData: y vals for target pigmnent, yDarkCounts: y vals for dark count, yReflectanceStanderd: y vals for reflectence standard

- DiffuseRefelctencePlotFolder
    - Function to plot FORS Spectra between 400 and 900nm, plots all txt files found in a directory skipping the first HeaderSize lines
    - args: Dir: Directory of data, DarkCountFileName: Full file path for dark count, ReflectenceStandardFileName: Full file path for ref standard, MatPlotLibColour: optional colour of plot, HeaderSize: Number of lines to ohmit at top of reading files

- DiffuseRefelctencePlotTxt
    - Function to plot FORS Spectra between 400 and 900nm, plots all txt files found in a directory skipping the first HeaderSize lines
    - args: RawFileName: Full file path for raw data, DarkCountFileName: Full file path for dark count, ReflectenceStandardFileName: Full file path for ref standard, MatPlotLibColour: optional colour of plot, HeaderSize: Number of lines to ohmit at top of reading files

- GeoMixing
    - Fuction to calculate the geometric mean for 2 spectra
    - args ySpectra1: y vals of one spectra, ySpectra2: y vals of second spectra

- MixingAlgorithm
    - Function to plot 2 FORS Spectra between 400 and 900nm and predict the mixed spectra
    - Args arrays 1by4 in order, first spectra, second spectra, mix spectra, prodicted mixed spectra
    - Args LegendName: Name of plot for legend, SpreadSheet: Name of google sheet found in gdrive containg data, DarkCountIndex: page index of dark count data found in spesifide google sheet starting at 0, ReflectenceStandardIndex: page index of reflectence standard data found in spesifide google sheet starting at 0, DataIndex: page index of raw data found in spesifide google sheet starting at 0, MatPlotLibColour: optional colour of plot, ModelUsed: optional picker for mixing algorithm used (defalt Geo)

- MixingFromSpectraSaveing
    - Function to save 2 FORS Spectra between 400 and 900nm and predict the mixed spectra
    - Args SpectraFilePath1: Full file path of first spectra, SpectraFilePath2: Full file path of second spectra, SaveDir: Save data location, MatPlotLibColour: optional colour of plot, ModelUsed: optional picker for mixing algorithm used (defalt Geo), HeaderSize: Number of lines to skip at top of files

- FullMixingSpectratoTXT
    - Function to take raw txt files from a spectrograph and convert them into a mixing peek matching plot
    - Args Raw1, Raw2, Mixed all indicate the 3 pices of recorded data shown below, SaveDir: Save Directory for Prossesed Data, HeaderSize: Number of skip lines in begginging of text file
    - fp: Raw Data, dc: Dark Count, rf: Reflectence Standard

- MixingLinePlotGenorator
    - Function to genorate a plot of the base pigments, the resultent mixed pigments and a estimation of mixing on one plot dividing up lines for prodicted peeks into quadrants
    - args Base1: File path for processed spectra of base pigment, Base2: File path for prossesed Spectra of other base pigment, Mix: File path for prossesed spectra for mixed pigment, Gen: File path for genorated pigment, Legend: legend names for plotting each line

- MixingPlotFullGen
    - Function to tern 3 raw readings of 2 base pigments and one mixed into a plot of pigments and there respective detected peeks
    - Args Raw1, Raw2, Mixed all indicate the 3 pices of recorded data shown below, SaveDir: Save Directory for Prossesed Data, HeaderSize: Number of skip lines in begginging of text file, Legend: Legends for plots
    - fp: Raw Data, dc: Dark Count, rf: Reflectence Standard

- ProssesedFolderPeekToTXT
    - Function to take prossesed data, genorate the values of the detected peeks using PeekFinderTXT and save the results
    - Args Dir: Target Directory containing prossesed data, SaveDir: Directory to save data in, Headersize: Number of skip lines at top of files

- PeekFinderTXT
    - Peek finder for already processed data
    - Args DataPath: Full file path for raw data, HeaderSize: number of lines to skip at top of file
    - Returns: Array of peek locations

- funPeekFinder
    - Peek finder for FORS reflectence spectra between 400 and 900nm
    - Args arrSpectralX: X values for spectrum, arrSpectralY, Y values for spectrum
    - Returns: Array of peek locations

