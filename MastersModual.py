import scipy as sp
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
import glob
from scipy import misc
from scipy import signal
import os
import shutil
import math


#print("Initialisation Complete")


def PlotProcessedSpectra(DataFilePath, MatPlotLibColour = None, HeaderSize = 0, Label = None, title = None):
    """
    Fuction to plot a already genorated FORS spectra
    Args DataFilePath: File to plot, MatPlotLibColour: Colour for Mat Plot Lib plotts, HeaderSize: Number of lines to skip for header, title: lable at top of graph
    """
    xRaw,yRaw = sp.loadtxt(DataFilePath, unpack = True, skiprows = HeaderSize)
    plt.plot(xRaw, yRaw, Label = Label, c = MatPlotLibColour)

    #plotting
    plt.ylim(0,1)
    plt.xlim(400, 900)
    plt.title(title)
    plt.xlabel("Wavelength / nm")
    plt.ylabel("Reflectence")
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
    plt.ylabel("Reflectence")
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
    plt.title("Final Spectras")
    plt.xlabel("Wavelength / nm")
    plt.ylabel("Reflectence")
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
    plt.title("Final Spectras")
    plt.xlabel("Wavelength / nm")
    plt.ylabel("Reflectence")
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
            yMixesData[i][a] = yMixesData[i][a]
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
        plt.ylabel("Reflectence/%")
        plt.plot(xMixesData[i], yMixesData[i])
        test = PeekFinder(xMixesData[i], yMixesData[i], debug)
        if debug == True:
            shutil.move("/content/debug.txt", SaveDir + "/" + names[i] + "Debug.txt")
        plt.text(400,0,test)
        plt.xlim(400,900)
        plt.ylim(0,1)
        plt.vlines(test, ymin = 0, ymax=1, colors = "r")
        plt.savefig(os.path.splitext(names[i])[0] + ".png")
        plt.clf()
        i = i + 1
    i = 0

    for filepath in (glob.glob("/content/*.png")):
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
    #objSmoothed = sp.interpolate.splrep(arrSpectralX, arrSpectralY, s=0.001)
    #arrSmoothedY = sp.interpolate.splev(arrSpectralX,objSmoothed)
    funInterpFunction = sp.interpolate.interp1d(arrSpectralX, arrSpectralY, fill_value= "extrapolate")

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
    Args: Dir: Directory of which to read from containing peek files, Headersize: Number of lines to skip at begginging of document, Debug: Adds additional outputs for debugging
    """

    arrFilePaths = []
    for filepath in (glob.glob(Dir + "/*Peeks.txt")):
        arrFilePaths.append(filepath)

    Peeks = [None]*len(arrFilePaths)

    i=0
    while i< len(arrFilePaths):
        Peeks[i] = sp.loadtxt(arrFilePaths[i], unpack = True, skiprows = Headersize)
        i = i+1
    i = 0

    if debug == True:
        print(Peeks)


    return(Peeks)

