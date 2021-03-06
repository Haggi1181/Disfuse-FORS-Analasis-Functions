import matplotlib
import scipy as sp
import numpy as np
from scipy import interpolate
from scipy import misc
import matplotlib.pyplot as plt
import glob
from scipy import misc
import os
import shutil
import math

plt.rcParams['font.size'] = 13
plt.tight_layout()

RunningLocation = os.path.abspath('')

#print("Initialisation Complete")


def PlotProcessedSpectra(DataFilePath, MatPlotLibColour = None, HeaderSize = 0, Label = None, title = None):
    """
    Fuction to plot a already genorated FORS spectra
    Args DataFilePath: File to plot, MatPlotLibColour: Colour for Mat Plot Lib plotts, HeaderSize: Number of lines to skip for header, title: lable at top of graph
    """
    xRaw,yRaw = sp.loadtxt(DataFilePath, unpack = True, skiprows = HeaderSize)
    plt.plot(xRaw, yRaw, c = MatPlotLibColour, label = Label)

    #plotting
    plt.ylim(0,1)
    plt.xlim(400, 900)
    plt.title(title)
    plt.xlabel("Wavelength / nm")
    plt.ylabel("Reflectance")
    plt.legend()

def PlotFolderTxt(Dir, MatPlotLibColour = None, HeaderSize = 0, title = None):
    """
    Fuction to plot all text files in a given dirrectory
    Args Dir: Directory to plot, MatPlotLibColour: Colour for Mat Plot Lib plotts, HeaderSize: Number of lines to skip for header
    """
    arrFilePaths = []
    arrFileName = []

    Dir = Dir + "/*.txt"
    for filepath in (glob.glob(Dir)):
        arrFilePaths.append(filepath)
        arrFileName.append(os.path.basename(filepath))
        #print("File Loaded")
    i = 0
    xRaw = [None]*len(arrFilePaths)
    yRaw = [None]*len(arrFilePaths)
    while i< len(arrFilePaths):
        xRaw[i],yRaw[i] = sp.loadtxt(arrFilePaths[i], unpack = True, skiprows = HeaderSize)
        plt.plot(xRaw[i], yRaw[i], label = arrFileName[i], c = MatPlotLibColour)
        i = i+1
        #print("Data Read")
    i = 0
    #plotting
    plt.ylim(0,1)
    plt.xlim(400, 900)
    plt.title(title)
    plt.xlabel("Wavelength / nm")
    plt.ylabel("Reflectance")
    plt.legend()

def SpectralGenAndSave(RawDataDir, DarkCountDir, ReflectanceStanderdDir, SaveDir, SaveFileName, HeaderSize = 14):
    """
    Function to genorate and save a single FORS spectra
    Args RawDataDir: Full file path for target data, DarkCountDir: Full file path for dark count, ReflectanceStanderdDir: Full file path for ref standard, SaveDir: Directory to save in, SafeFileName: name of resultant file auto adds .txt, HeaderSize: Number of lines to ohmit at top of reading files
    """

    xRaw, yRaw = sp.loadtxt(RawDataDir, unpack = True, skiprows = HeaderSize)
    xDark, yDark = sp.loadtxt(DarkCountDir, unpack = True, skiprows = HeaderSize)
    xRef, yRef = sp.loadtxt(ReflectanceStanderdDir, unpack = True, skiprows = HeaderSize)

    ySpectra = SpectralGen(yRaw, yDark, yRef)

    finalpath = SaveDir + "/" + SaveFileName + ".txt"
    sp.savetxt(finalpath, np.column_stack([xRaw, ySpectra]))

def SpectralDatabaseCreator(Dir, DarkCountFilePathName, ReflectenceStandardFilePathName, SaveDir, HeaderSize = 14):
    """
    Function to genorate a spectral database from spectra found in a single folder saving calculated spectra into one folder
    args: Dir: Directory of data, DarkCountFilePathName: Full file path for dark count, ReflectenceStandardFilePathName: Full file path for ref standard, SaveDir: Directory to save in, HeaderSize: Number of lines to ohmit at top of reading files
    """
    arrFilePaths = []
    arrFileName = []
    Dir = Dir + "/*.txt"
    for filepath in (glob.glob(Dir)):
        arrFilePaths.append(filepath)
        arrFileName.append(os.path.splitext(os.path.basename(filepath))[0])
        #print("File Loaded")


    i=0
    while i != len(arrFilePaths):
        if arrFilePaths[i] == DarkCountFilePathName or arrFilePaths[i] == ReflectenceStandardFilePathName:
            arrFilePaths[i] = "asd"
            arrFileName[i] = "asd"
            #print("FileOmmitted")
        i = i + 1
        #print("File Check Loop")

    temp = filter(lambda c: c != "asd", arrFilePaths)
    arrFilePaths = list(temp)
    temp = filter(lambda c: c != "asd", arrFileName)
    arrFileName = list(temp)
    #print("Files Trimmed")

    i=0
    while i != len(arrFilePaths):
        SpectralGenAndSave(arrFilePaths[i], DarkCountFilePathName, ReflectenceStandardFilePathName, SaveDir, arrFileName[i])
        i=i+1


def SpectralGen(yRawData, yDarkCounts, yReflectanceStanderd):
    """
    Function to Calculate FORS Spectra between 400 and 900nm
    Args yRawData: y vals for target pigmnent, yDarkCounts: y vals for dark count, yReflectanceStanderd: y vals for reflectence standard
    """
    #loops to calculate equation shown in work book
    difference = []
    zip_object = zip(yRawData, yDarkCounts)
    for yRawData_i, yDarkCounts_i in zip_object:
        difference.append(yRawData_i-yDarkCounts_i)
    difference2 = []
    zip_object2 = zip(yReflectanceStanderd, yDarkCounts)
    for yReflectanceStanderd_i, yDarkCounts_i in zip_object2:
        difference2.append(yReflectanceStanderd_i-yDarkCounts_i)
    ySpectra = []
    zip_object3 = zip(difference, difference2)
    for difference_i, difference2_i in zip_object3:
        ySpectra.append((difference_i/difference2_i))
    return ySpectra

def DiffuseRefelctencePlotFolder(Dir, DarkCountFilePathName, ReflectenceStandardFilePathName, MatPlotLibColour=None, HeaderSize = 14):
    """
    Function to plot FORS Spectra between 400 and 900nm, plots all txt files found in a directory skipping the first HeaderSize lines
    args: Dir: Directory of data, DarkCountFilePathName: Full file path for dark count, ReflectenceStandardFilePathName: Full file path for ref standard, MatPlotLibColour: optional colour of plot, HeaderSize: Number of lines to ohmit at top of reading files
    """

    arrFilePaths = []
    arrFileName = []
    Dir = Dir + "/*.txt"
    for filepath in (glob.glob(Dir)):
        arrFilePaths.append(filepath)
        arrFileName.append(os.path.basename(filepath))
        #print("File Loaded")


    i=0
    while i != len(arrFilePaths):
        if arrFilePaths[i] == DarkCountFilePathName or arrFilePaths[i] == ReflectenceStandardFilePathName:
            arrFilePaths[i] = "asd"
            arrFileName[i] = "asd"
            #print("FileOmmitted")
        i = i + 1
        #print("File Check Loop")

    temp = filter(lambda c: c != "asd", arrFilePaths)
    arrFilePaths = list(temp)
    temp = filter(lambda c: c != "asd", arrFileName)
    arrFileName = list(temp)
    #print("Files Trimmed")
    i=0

    xRaw = [None]*len(arrFilePaths)
    yRaw = [None]*len(arrFilePaths)

    while i< len(arrFilePaths):
        xRaw[i],yRaw[i] = sp.loadtxt(arrFilePaths[i], unpack = True, skiprows = HeaderSize)
        i = i+1
        #print("Data Read")
    i = 0

    xDark, yDark = sp.loadtxt(DarkCountFilePathName, unpack = True, skiprows = HeaderSize)

    xRef, yRef = sp.loadtxt(ReflectenceStandardFilePathName, unpack = True, skiprows = HeaderSize)
    #print("Ref Dark Load")
    i=0
    while i != len(arrFilePaths):
        ySpectra = SpectralGen(yRaw[i], yDark, yRef)
        #print("Spectra Genorated")
        plt.plot(xRaw[i], ySpectra, label = arrFileName[i], c = MatPlotLibColour)
        #print(i)
        i=i+1
    #plotting
    plt.ylim(0,1)
    plt.xlim(400, 900)
    plt.title("Final Spectra")
    plt.xlabel("Wavelength / nm")
    plt.ylabel("Reflectance")
    plt.legend()

def DiffuseRefelctencePlotTxt(RawFileName, DarkCountFilePathName, ReflectenceStandardFilePathName, MatPlotLibColour=None, HeaderSize = 14, Legend = ""):
    """
    Function to plot FORS Spectra between 400 and 900nm, plots all txt files found in a directory skipping the first HeaderSize lines
    args: RawFileName: Full file path for raw data, DarkCountFilePathName: Full file path for dark count, ReflectenceStandardFilePathName: Full file path for ref standard, MatPlotLibColour: optional colour of plot, HeaderSize: Number of lines to ohmit at top of reading files
    """
    xRaw, yRaw = sp.loadtxt(RawFileName, unpack = True, skiprows = HeaderSize)

    xDark, yDark = sp.loadtxt(DarkCountFilePathName, unpack = True, skiprows = HeaderSize)

    xRef, yRef = sp.loadtxt(ReflectenceStandardFilePathName, unpack = True, skiprows = HeaderSize)
    
    ySpectra = SpectralGen(yRaw, yDark, yRef)
    plt.plot(xRaw, ySpectra, label = Legend, c = MatPlotLibColour)
    #plotting
    plt.ylim(0,1)
    plt.xlim(400, 900)
    plt.title("Final Spectra")
    plt.xlabel("Wavelength / nm")
    plt.ylabel("Reflectance")
    plt.legend()

def GeoMixing(ySpectra1, ySpectra2):
    """
    Fuction to calculate the geometric mean for 2 spectra
    args ySpectra1: y vals of one spectra, ySpectra2: y vals of second spectra
    """
    i=0
    yGeometricMean = []
    while i < len(ySpectra2):
        yGeometricMean.append((ySpectra1[i]*ySpectra2[i])**(1/2))
        i = i+1
    return yGeometricMean


def MixingFromSpectraSaveing(SpectraFilePath1, SpectraFilePath2, SaveDir, ModelUsed = "Geo", HeaderSize = 0):
    """
    Function to save 2 FORS Spectra between 400 and 900nm and predict the mixed spectra
    Args SpectraFilePath1: Full file path of first spectra, SpectraFilePath2: Full file path of second spectra, SaveDir: Save data location, MatPlotLibColour: optional colour of plot, ModelUsed: optional picker for mixing algorithm used (defalt Geo), HeaderSize: Number of lines to skip at top of files
    """
    

    xSpectra1, ySpectra1 = sp.loadtxt(SpectraFilePath1, unpack = True, skiprows = HeaderSize)
    
    xSpectra2, ySpectra2 = sp.loadtxt(SpectraFilePath2, unpack = True, skiprows = HeaderSize)


    #Runs the selected mixing Algorithm
    if ModelUsed == "Geo":
        yMixed = GeoMixing(ySpectra1, ySpectra2)
    
    finalpath = SaveDir+(os.path.splitext(os.path.basename(SpectraFilePath1))[0])+"&"+(os.path.splitext(os.path.basename(SpectraFilePath2))[0])+".txt"

    sp.savetxt(finalpath, np.column_stack([xSpectra1, yMixed]))

def FullMixingSpectratoTXT(fpRaw1, dcRaw1, rfRaw1, fpRaw2, dcRaw2, rfRaw2, fpRawMixed, dcRawMixed, rfRawMixed, SaveDir, HeaderSize = 14):
    """
    Function to take raw txt files from a spectrograph and convert them into a mixing peek matching plot
    Args Raw1, Raw2, Mixed all indicate the 3 pices of recorded data shown below, SaveDir: Save Directory for Prossesed Data, HeaderSize: Number of skip lines in begginging of text file
    fp: Raw Data, dc: Dark Count, rf: Reflectence Standard
    """

    SpectralGenAndSave(fpRaw1, dcRaw1, rfRaw1, SaveDir, os.path.splitext(os.path.basename(fpRaw1))[0], HeaderSize)
    SpectralGenAndSave(fpRaw2, dcRaw2, rfRaw2, SaveDir, os.path.splitext(os.path.basename(fpRaw2))[0], HeaderSize)
    SpectralGenAndSave(fpRawMixed, dcRawMixed, rfRawMixed, SaveDir, os.path.splitext(os.path.basename(fpRawMixed))[0], HeaderSize)

    SpectraPath1 = SaveDir + os.path.basename(fpRaw1)
    SpectraPath2 = SaveDir + os.path.basename(fpRaw2)

    MixingFromSpectraSaveing(SpectraPath1, SpectraPath2, SaveDir)

def MixingLinePlotGenorator(Base1, Base2, Mix, Gen, Legend):
    """
    Function to genorate a plot of the base pigments, the resultent mixed pigments and a estimation of mixing on one plot dividing up lines for prodicted peeks into quadrants
    args Base1: File path for processed spectra of base pigment, Base2: File path for prossesed Spectra of other base pigment, Mix: File path for prossesed spectra for mixed pigment, Gen: File path for genorated pigment, Legend: legend names for plotting each line
    """
    PlotProcessedSpectra(Base1, MatPlotLibColour = "r", Label = Legend[0])
    plt.vlines(PeekFinderTXT(Base1),ymin = 0, ymax=0.25, colors = "r")

    PlotProcessedSpectra(Base2, MatPlotLibColour = "g", Label = Legend[1])
    plt.vlines(PeekFinderTXT(Base2),ymin = 0.25, ymax=0.5, colors = "g")

    PlotProcessedSpectra(Mix, MatPlotLibColour = "y", Label = Legend[2])
    plt.vlines(PeekFinderTXT(Mix),ymin = 0.5, ymax=0.75, colors = "y")

    PlotProcessedSpectra(Gen, MatPlotLibColour = "b", Label = Legend[3], title = "FPCG")
    plt.vlines(PeekFinderTXT(Gen),ymin = 0.75, ymax=1, colors = "b")

def MixingPlotFullGen(fpRaw1, dcRaw1, rfRaw1, fpRaw2, dcRaw2, rfRaw2, fpRawMixed, dcRawMixed, rfRawMixed, SaveDir, Legend = [None, None, None, None]):
    """
    Function to tern 3 raw readings of 2 base pigments and one mixed into a plot of pigments and there respective detected peeks
    Args Raw1, Raw2, Mixed all indicate the 3 pices of recorded data shown below, SaveDir: Save Directory for Prossesed Data, HeaderSize: Number of skip lines in begginging of text file, Legend: Legends for plots
    fp: Raw Data, dc: Dark Count, rf: Reflectence Standard
    """
    FullMixingSpectratoTXT(fpRaw1, dcRaw1, rfRaw1, fpRaw2, dcRaw2, rfRaw2, fpRawMixed, dcRawMixed, rfRawMixed, SaveDir)
    
    Base1=SaveDir + os.path.basename(fpRaw1)
    Base2=SaveDir + os.path.basename(fpRaw2)
    Mix=SaveDir + os.path.basename(fpRawMixed)
    Gen=SaveDir+(os.path.splitext(os.path.basename(Base1))[0])+"&"+(os.path.splitext(os.path.basename(Base2))[0])+".txt"

    MixingLinePlotGenorator(Base1, Base2, Mix, Gen, Legend)

def ProssesedDataToPeekDatabase(Dir, SaveDir, Headersize = 0, debug = False):
    """
    Function to take prossesed data, genorate the values of the detected peeks using PeekFinderTXT and save the results
    Args Dir: Target Directory containing prossesed data, SaveDir: Directory to save data in, Headersize: Number of skip lines at top of files
    """
    arrFilePaths = []
    arrFileName = []

    Dir = Dir + "/*.txt"
    for filepath in (glob.glob(Dir)):
        arrFilePaths.append(filepath)
        arrFileName.append(os.path.splitext(os.path.basename(filepath))[0])
        #print("File Loaded")
    i = 0
    temp = []    
    while i< len(arrFilePaths):
        temp = PeekFinderTXT(arrFilePaths[i], Headersize, debug)
        finalpath = SaveDir + "/" + arrFileName[i] + "Peeks" + ".txt"
        sp.savetxt(finalpath, temp)
        if debug == True:
            shutil.move(RunningLocation+"/debug.txt", SaveDir + "/" + arrFileName[i] + "Debug.txt")
        i = i+1
        #print("Data Read")
    i = 0

def RawDataToMatchingAlgorithmDatabase(Dir, DarkCountFilePathName, ReflectenceStandardFilePathName, SaveDir, HeaderSize = 14, debug = False):
    """
    Function to create a database from a directory of unprossesed fors data, outputs a series of prossesed spectra and assosiated detected peeks
    args: Dir: Directory of data, DarkCountFilePathName: Full file path for dark count, ReflectenceStandardFilePathName: Full file path for ref standard, SaveDir: Directory to save in, HeaderSize: Number of lines to ohmit at top of reading files, debug: causes additional text files containing debug information in the peek detection step to be saved
    """

    SpectralDatabaseCreator(Dir, DarkCountFilePathName, ReflectenceStandardFilePathName, SaveDir, HeaderSize)
    ProssesedDataToPeekDatabase(SaveDir, SaveDir, 0, debug)

def PeekFinderTXT(DataPath, HeaderSize = 0, debug = False):
    """
    Peek finder for already processed data
    Args DataPath: Full file path for raw data, HeaderSize: number of lines to skip at top of file
    Returns: Array of peek locations
    """

    SpectraX, SpectraY = sp.loadtxt(DataPath, unpack = True, skiprows = HeaderSize)

    a = 0
    while a < len(SpectraX):
        if SpectraX[a] < 400 or SpectraX[a] > 900:
            SpectraX[a] = 0
            SpectraY[a] = 0
        a = a + 1
    tempx = filter(lambda c: c != 0, SpectraX)
    SpectraX = list(tempx)
    tempy = filter(lambda c: c != 0, SpectraY)
    SpectraY = list(tempy)




    arr0 = PeekFinder(SpectraX, SpectraY, debug)
    return arr0

def PeekFinderPlotGenorater(Dir, SaveDir, Headersize = 0, debug = False):
    """
    Function to genorate plots for all prossesed data found in a directory
    Args Dir: Directory containting all data, SaveDir: Location to save plots, Headersize: Number of likes to skip at beggining of reading data
    """
    #loading spectral data base
    arrFilePaths = []


    for filepath in (glob.glob(Dir + "/*.txt")):
        arrFilePaths.append(filepath)

    names = [os.path.basename(x) for x in glob.glob(Dir + "/*.txt")]

    i = 0

    xMixesData = [None]*len(arrFilePaths)
    yMixesData = [None]*len(arrFilePaths)


    while i< len(arrFilePaths):
        xMixesData[i],yMixesData[i] = sp.loadtxt(arrFilePaths[i], unpack = True)
        i = i+1
    i = 0

    while i < len(arrFilePaths):
        a = 0
        while a < len(xMixesData[i]):
            if xMixesData[i][a] < 400 or xMixesData[i][a] > 900:
                xMixesData[i][a] = 0
                yMixesData[i][a] = 0
            a = a + 1
        tempx = filter(lambda c: c != 0, xMixesData[i])
        xMixesData[i] = list(tempx)
        tempy = filter(lambda c: c != 0, yMixesData[i])
        yMixesData[i] = list(tempy)
        plt.title(names[i])
        plt.xlabel("Wavelength/nm")
        plt.ylabel("Reflectance/%")
        plt.plot(xMixesData[i], yMixesData[i])
        test = PeekFinder(xMixesData[i], yMixesData[i], debug)
        if debug == True:
            shutil.move(RunningLocation+"/debug.txt", SaveDir + "/" + names[i] + "Debug.txt")
        plt.text(400,0,test)
        plt.xlim(400,900)
        plt.ylim(0,1)
        plt.vlines(test, ymin = 0, ymax=1, colors = "r")
        plt.savefig(os.path.splitext(names[i])[0] + ".png")
        plt.clf()
        i = i + 1
    i = 0

    for filepath in (glob.glob(RunningLocation+"/*.png")):
        shutil.move(filepath, SaveDir)

def PeekFinder(arrSpectralX, arrSpectralY, debug = False):
    """
    Peek finder for FORS reflectence spectra between 400 and 900nm
    Args arrSpectralX: X values for spectrum, arrSpectralY, Y values for spectrum
    Returns: Array of peek locations
    """
    #trims out nul values, a value that may be genorated by extreamly small floats, fixes issues with splev
    a = 0
    while a != len(arrSpectralX):
        if math.isnan(arrSpectralY[a])==True:
            arrSpectralX[a] = -10000
            arrSpectralY[a] = -10000
        a = a + 1
    tempx = filter(lambda c: c != -10000, arrSpectralX)
    arrSpectralX = list(tempx)
    tempy = filter(lambda c: c != -10000, arrSpectralY)
    arrSpectralY = list(tempy)


    #Smoothes input array
    objSmoothed = sp.interpolate.splrep(arrSpectralX, arrSpectralY, s=0.001)
    arrSmoothedY = sp.interpolate.splev(arrSpectralX,objSmoothed)
    funInterpFunction = sp.interpolate.interp1d(arrSpectralX, arrSmoothedY, fill_value= "extrapolate")

    arrSpectralX = sp.linspace(200, 1000, 400)
    arrSmoothedY = funInterpFunction(arrSpectralX)

    objSmoothed = sp.interpolate.splrep(arrSpectralX, arrSmoothedY, s=0.003)
    arrSmoothedY = sp.interpolate.splev(arrSpectralX,objSmoothed)
    #plt.plot(arrSpectralX, arrSmoothedY, linewidth = 1, color = "k")#Testing plot
    #interpulates the genorated smoothed array to fill in gaps for processing
    funInterpFunction = sp.interpolate.interp1d(arrSpectralX, arrSmoothedY, fill_value= "extrapolate")

    #restricts the array of values down to the range we are intrested with so that the floor test is uneffected by noise outside of this region
    a = 0
    while a < len(arrSmoothedY):
        if arrSpectralX[a] < 400 or arrSpectralX[a] > 900:
            arrSmoothedY[a] = 0
            arrSpectralX[a] = 0
        a = a + 1
    tempx = filter(lambda c: c != 0, arrSpectralX)
    arrSpectralX = list(tempx)
    tempy = filter(lambda c: c != 0, arrSmoothedY)
    arrSmoothedY = list(tempy)

    #initialises veriables for loop below
    arrDerivative = []
    arrDerWave = []
    arr0 = []
    runningpeek = 400.0#typical begginging of ramp up section of spectra

    if debug == True:
        f = open("debug.txt", "w")


    i=400.0#at 400, beggingin of data range
    while i <= 900.0:#Loops for length of spectral data
        #resets inflection test state veriable
        booInflectionTest = False
        if misc.derivative(funInterpFunction,(i-1))>0 and misc.derivative(funInterpFunction,(i+1))<0:#checks for a change in sign of derivative
            booInflectionTest = True#sets state if true
        #performes the following checks, a threshhold check on the value of the derivative that it is close enough to 0, a value to make shore the peek isnt on the "floor" of the spectra, to see if there was a change in sign in the derivative and to see if the values before and after are smaller than the current one
        if abs(misc.derivative(funInterpFunction, i))<0.001 and funInterpFunction(i)>min(arrSmoothedY)*1.1 and booInflectionTest == True and funInterpFunction(i-1)<funInterpFunction(i) and funInterpFunction(i+1)<funInterpFunction(i):
            arr0.append(i)#if all tests are true saves the location of the peek
            #print("Maxima Detected at:", i, " Value:",funInterpFunction(i), " Derivative:", misc.derivative(funInterpFunction, i)) #Debug Print

        #less of a peek and more of a change in "trajectory" of line so function below is specialised
        if funInterpFunction(i)-funInterpFunction(i-1) > 0.005:#looks to see where the change in the lines trajectory is small enough to be considerd a mode change
            runningpeek = i#moves the value of the detected peek up untill the mode changes
            if debug == True:
                f.write("Running Peek Moved \n")
        #Debuging Prints
        if debug == True:
            savetext = "at:" + str(i) + " Value:" + str(funInterpFunction(i)) + " Derivative:" + str(misc.derivative(funInterpFunction, i)) + "\n"
            f.write(savetext)
            if abs(misc.derivative(funInterpFunction, i))<0.001 and funInterpFunction(i)>min(arrSmoothedY)*1.1 and booInflectionTest == True and funInterpFunction(i-1)<funInterpFunction(i) and funInterpFunction(i+1)<funInterpFunction(i):
                f.write("____________________PeekDectected____________________"+ "\n")
            if abs(misc.derivative(funInterpFunction, i))<0.001:
                f.write("   Pass Derivative Threshhold"+ "\n")
            if funInterpFunction(i)>min(arrSmoothedY)*1.1:
                f.write("   Pass Floor"+ "\n")
            if booInflectionTest == True:
                f.write("   Pass Inflection"+ "\n")
            if funInterpFunction(i-1)<funInterpFunction(i):
                f.write("   Pass Turning Low"+ "\n")
            if funInterpFunction(i+1)<funInterpFunction(i):
                f.write("   Pass Turning High"+ "\n")
        i=i+1.0
    #check to see if the found peek is close enough to the max reflectivity of the overall spectra to be considerd the leveling off point
    if (funInterpFunction(runningpeek) <= max(arrSmoothedY)-0.2) or (runningpeek <= 401.0):
        if debug == True:
            f.write("RunningPeek Not Saved \n")
            f.write(str(runningpeek <= 401.0) + "\n")
            f.write(str(funInterpFunction(runningpeek) <= max(arrSmoothedY)-0.2) + "\n")
            f.write("RunningPeek Val:" + str(runningpeek))
        pass
    else:
        if debug == True:
            f.write("RunningPeek Saved \n")
            f.write(str(runningpeek <= 401.0) + "\n")
            f.write(str(funInterpFunction(runningpeek) <= max(arrSmoothedY)-0.2) + "\n")
            f.write("RunningPeek Val:" + str(runningpeek))
        arr0.append(runningpeek)
    return(arr0)

def PeekLoading(Dir, Headersize = 0, debug = False):
    """
    a function to load data from spectral peek files genorated by this modual
    Args Dir: Directory of which to read from containing peek files, Headersize: Number of lines to skip at begginging of document, debug: Adds additional outputs for debugging
    Returns: Array of peeks
    """
    
    arrFilePaths = []
    for filepath in (glob.glob(Dir + "/*Peeks.txt")):
        arrFilePaths.append(filepath)

    if debug == True:
        print("File Path Names: ",arrFilePaths)

    Peeks = [None]*len(arrFilePaths)

    i=0
    while i< len(arrFilePaths):
        Peeks[i] = sp.loadtxt(arrFilePaths[i], skiprows = Headersize)
        i = i+1
    i = 0


    return(Peeks)


def MatchingAlgorithm(PeekDir, PathUnknownSpectra, debug = False):
    """
    Function to match a unknonw spectra to a known one from a database via a peek matching method
    Args PeekDir: Directory containing peek data, PathUnknownSpectra: File path of the pre prossesed spectra of unknonw pigment, debug: flag in function to print additional debugging information
    """
    peeks = PeekLoading(PeekDir)
    UnknownPeeks = PeekFinderTXT(PathUnknownSpectra, debug = debug)
    #peeks = sortRowWise(peeks)
    #UnknownPeeks = UnknownPeeks.sort()

    names = [os.path.basename(x) for x in glob.glob(PeekDir + "/*Peeks.txt")]


    PeeksList = list()
    for row in peeks:
        PeeksList.append(np.atleast_1d(row).tolist())

    #PeeksList = sortRowWise(PeeksList)

    if debug == True:
        print("Peeks",peeks)
        print("UnknownPeeks",UnknownPeeks)
        print("typePeeks",type(peeks))
        print("typePeeks0",type(peeks[0]))
        print("lenPeeks",len(peeks))
        print("peeks",peeks)
        print("PeeksList", PeeksList)
        print("typePeeksList", type(PeeksList))
        print("typePeeksList0", type(PeeksList[0]))
        print("typePeeksList1", type(PeeksList[1]))



    distance = list()
    i=0
    removallist = list()

    while i != len(PeeksList):
        dist = absolute(UnknownPeeks, PeeksList[i], debug=debug)
        if len(dist[0]) == 0:
            removallist.append(i)
        else:
            distance.append(dist)
        i = i + 1
    i = 0

    if debug == True:
        print("Remval List Indexes: ", removallist)
        print("Pre Removal Results: ", PeeksList)
        print("Pre Removal Names: ", names)


    PeeksList = [v for i, v in enumerate(PeeksList) if i not in removallist]
    names = [v for i, v in enumerate(names) if i not in removallist]

    if debug == True:
        print("Post RemovalResults: ", PeeksList)
        print("Post RemovalNames: ", names)

    i = 0
    
    if debug == True:
        print("len of dist: ",len(distance))
        print("len of dist [0]: ",len(distance[0]))
        print("len of dist [0][0]: ",len(distance[0][0]))
        print("dist: ",distance)
        print("dist [0]: ",distance[0])
        print("dist [0][0]: ",distance[0][0])
        print("dist [0][0][0]: ",distance[0][0][0])

    templine = list()
    MinMatrix = list()
    i=0
    while i != len(distance):
        if debug == True:
            print("MatrixLoop i: ", i)
        j=0
        while j != len(distance[i]):
            if debug == True:
                print("MatrixLoop j: ", j)
                print("Dist[i][j]: ", distance[i][j])
                print("Dist Min[i][j]: ",min(distance[i][j]))
            templine.append(min(distance[i][j]))
            if debug == True:
                print("Temp Line: ", templine)
            j=j+1
        MinMatrix.append(templine)
        templine = []
        if debug == True:
            print("MinMatrix: ", MinMatrix)
        i=i+1
    if debug == True:
        print(MinMatrix)

    MinMatrixSum = list()
    i=0
    while i != len(MinMatrix):
        MinMatrixSum.append(sum(MinMatrix[i]))
        i = i + 1
    if debug == True:
        print(MinMatrixSum)

    results = (np.atleast_1d(sp.argsort(MinMatrixSum)).tolist())


    if debug == True:
        print("Results: ", results)
    i=0

    print("Results from most to least likely, for confidence closer to 0 is better:")
    while i != len(results):
        print("Result:",i+1," ",(names[results[i]]).replace("Peeks", ""), " Confidence: ", MinMatrixSum[results[i]])
        i = i + 1


def absolute(m1,m2, debug = False):
    """
    function to create a array of all absolute values between elements inside of 2 given arrays
    args m1: 1d first array, m2: 1d second array, debug: flag to print aditional information for debugging information
    return: 2d array of the differences between all elements of m1 and m2
    """

    m1len = len(m1)

    m2len = len(m2)


    if debug == True:
        print("absolute length of m1: ",m1len)
        print("absolute length of m2: ",m2len)

    val = np.zeros((m1len, m2len))
    vals = val.tolist()
    #vals = [[None]*m1len,[None]*m2len]

    i=0
    j=0
    while i != m1len:
        if debug == True:
            print("Absolute Loop i: ",i, " Element: ", m1[i])
        while j != m2len:
            if debug == True:
                print("Absolute Loop j: ",j, " Element: ", m2[j])
            vals[i][j] = ((m1[i]-m2[j])**2)**(1/2)
            j = j + 1
        j=0
        i = i + 1

    return(vals)


def sortRowWise(m):
     
    # loop for rows of matrix
    temp1 = len(m)
    for i in range(temp1):
        print(type(m))
        print(type(m[i]))
        temp2 = len(m[i])
        # loop for column of matrix
        for j in range(temp2):
             
            # loop for comparison and swapping
            for k in range(temp2 - j - 1):
                 
                if (m[i][k] > m[i][k + 1]):
                     
                    # swapping of elements
                    t = m[i][k]
                    m[i][k] = m[i][k + 1]
                    m[i][k + 1] = t
                     
    return m

def PersentageToFractionalTXT(Dir, SaveDir):
    """
    Function to take persentage y values to fractional (/100) for already prossesed fors data for all txt files in a given folder
    args Dir: Directory containing input files, SaveDir: Directory to save to
    """
    arrFilePaths = []
    arrFileName=[]
    
    for filepath in (glob.glob(Dir + "/*.txt")):
        arrFilePaths.append(filepath)
        arrFileName.append(os.path.basename(filepath))


    xMixesData = [None]*len(arrFilePaths)
    yMixesData = [None]*len(arrFilePaths)

    i=0
    while i< len(arrFilePaths):
        xMixesData,yMixesData = sp.loadtxt(arrFilePaths[i], unpack = True)
        yMixesData = yMixesData/100
        finalpath = SaveDir + "/" + arrFileName[i]
        sp.savetxt(finalpath, np.column_stack([xMixesData, yMixesData]))
        i = i+1


def LineMatchingAlgorithm(SpectraDir, PathUnknownSpectra, debug = False, HeaderSize = 0):
    """
    Function to match a unknonw spectra to a known one from a database via a peek matching method using deviding algorithm
    Args PeekDir: Directory containing peek data, PathUnknownSpectra: File path of the pre prossesed spectra of unknonw pigment, debug: flag in function to print additional debugging information, HeaderSize: number of lines to skip in raw data files
    """
    
    arrFilePaths = []
    arrFileName = []

    SpectraDir = SpectraDir + "/*.txt"
    for filepath in (glob.glob(SpectraDir)):
        arrFilePaths.append(filepath)
        arrFileName.append(os.path.basename(filepath))
        #print("File Loaded")
    i = 0
    xRaw = [None]*len(arrFilePaths)
    yRaw = [None]*len(arrFilePaths)
    while i< len(arrFilePaths):
        xRaw[i],yRaw[i] = sp.loadtxt(arrFilePaths[i], unpack = True, skiprows = HeaderSize)
        i = i+1

    xUnknown,yUnknown = sp.loadtxt(PathUnknownSpectra, unpack = True, skiprows = HeaderSize)

    ydivlines = [None]*len(arrFilePaths)
    i=0
    while i< len(arrFilePaths):
        ydivlines[i] = LineDivider(yRaw[i], yUnknown)
        i = i + 1

    sx = [None]*len(arrFilePaths)
    sy = [None]*len(arrFilePaths)

    i=0
    while i< len(arrFilePaths):
        sx[i], sy[i] = smoother(xRaw[i], ydivlines[i])
        #plt.plot(sx[i], sy[i])
        i = i + 1

    derivatives = [None]*len(arrFilePaths)


    i=0
    while i< len(arrFilePaths):
        derivatives[i] = netabsderivative(sx[i], sy[i])
        #print(derivatives[i])
        i = i + 1
    i=0

    results = sp.argsort(derivatives)
    #print(results)


    names = [os.path.basename(x) for x in glob.glob(SpectraDir)]
    #print(names)

    while i != len(results):
        print("Result:",i+1," ",(names[results[i]]), " Confidence: ", derivatives[results[i]])
        i = i + 1



def netabsderivative(x,y):
    """
    Function to calculate the sum of the absolutes of the gradients of a line
    Args x: list of x vals, y: list of y vals
    """
    # funInterpFunction = sp.interpolate.interp1d(x, y, fill_value= "extrapolate")

    # der = sp.misc.derivative(funInterpFunction)

    dx = x[1] - x[0]
    der = np.gradient(y, dx)

    test = np.sum([abs(a) for a in der])

    plt.plot(x,der)

    return test


def smoother(x,y):
    """
    Function to perform smoothing operations on a FORS spectra
    Args x: list of x vals, y: list of y vals
    """
    #trims out nul values, a value that may be genorated by extreamly small floats, fixes issues with splev
    # print(x)
    # print(y)
    a = 0
    while a != len(x):
        if math.isnan(x[a])==True:
            x[a] = -10000
            y[a] = -10000
        a = a + 1
    tempx = filter(lambda c: c != -10000, x)
    x = list(tempx)
    tempy = filter(lambda c: c != -10000, y)
    y = list(tempy)

    # print(x)
    # print(y)
    a = 0
    while a < len(y):
        if x[a] < 400 or x[a] > 900:
            y[a] = -10000
            x[a] = -10000
        a = a + 1
    tempx = filter(lambda c: c != -10000, x)
    x = list(tempx)
    tempy = filter(lambda c: c != -10000, y)
    y = list(tempy)

    #print(x)
    #print(y)

    funInterpFunction = sp.interpolate.interp1d(x, y, fill_value= "extrapolate")


    x = sp.linspace(400, 900, 500)
    y = funInterpFunction(x)

    #Smoothes input array
    objSmoothed = sp.interpolate.splrep(x, y, s=0.001)
    sy = sp.interpolate.splev(x,objSmoothed)
    
    funInterpFunction = sp.interpolate.interp1d(x, sy, fill_value= "extrapolate")

    fx = sp.linspace(400, 900, 500)
    fy = funInterpFunction(x)

    return fx, fy


def LineDivider(yDatabase, yUnknown):
    """
    Function to find the ratio of y vals for 2 sets of datas y values
    Args yDatabase: Set of y values of known spectra, yUnknown: Set of y values for spectra of unknown pigment
    """
    div = []
    zip_object = zip(yDatabase, yUnknown)
    for yDatabase_i, yUnknown_i in zip_object:
        div.append(yUnknown_i - yDatabase_i)
    
    return div



