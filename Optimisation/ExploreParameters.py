# Start of script
print("Script is running...")
# Import all packages
import os
os.environ["MKL_CBWR"] = "SSE2"
os.environ["OMP_NUM_THREADS"] = '1'
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import time
from pyDOE import *
# import odeFunctions
# odeFunctions.use('TkAgg')
from odeFunctions import *
from multiprocessing import Pool
import itertools
import re
import sys
import pandas as pd
from scipy.optimize import least_squares
import contextlib
from scipy.io import loadmat, savemat
import seaborn as sns
# import scikits.odes.ode as sciode
import platform

from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

print("All packages imported")
print("Python version %s.%s.%s" % (sys.version_info[0], sys.version_info[1], sys.version_info[2]))
computerName = platform.node()

###################################################
# Part I - USER DEFINED VARIABLES
###################################################

# Define the model of interest
MODEL = "M019"
TRMNT = "CDDP"

# Define the path to the input folder containing the model output, i.e. the path where the MODEL directory is located
inputPath = "/Users/muriel/Documents/LACDR/Projects/PHH/Models/v2_firstRevision/Output/"
#inputPath = "/Users/murielheldring/TransferredFiles/" # Write down the path to the folder that contains the folders of the
# models and their parameter estimation outputs

# Define the path to the ode system file, i.e. the path where the input files are located
odeModelFolder = "/Users/muriel/Documents/LACDR/Projects/PHH/Models/v2_firstRevision/" # Write down the path to the folder that
# contains the file with the ODE system

# Define the path to the data file, i.e. the path where the input files are located
dataFolder = "/Users/muriel/Documents/LACDR/Projects/PHH/Models/v2_firstRevision/" # Write down the path to the folder that
# contains the file with the data

###################################################
# Part II - Run script
###################################################

# Create the input folder name
inputFolder = inputPath + MODEL + "/" + TRMNT
# Define the name of the output folder
outputFolder = "ParameterIdentifiability/"

# Check output path availability
outPath = inputFolder + "/" + outputFolder
if os.path.isdir(outPath):
    print(outPath + " does already exist. Overwriting/adding output...")
else:
    print(outPath + " does not exist.\nCreating a new directory and writing the output to " + outPath)
    os.mkdir(outPath)


# Define paths to data
Dates = os.listdir(inputFolder)
parmset_paths = []
#setting_paths = []
for Date in Dates:
    path1 = inputFolder + "/" + Date
    if os.path.isdir(path1):
        path2 = path1 + "/"
        # Runs = os.listdir(path1)
        # for Run in Runs:
        #     path2 = path1 + "/" + Run + "/"
        #     if os.path.isdir(path2):
                # Get the files that contain the parameter estimates
        files = os.listdir(path2)
        for file in files:
            if re.search("\d{1}.\d{2}.csv", file):
                parmset_paths.append(path2 + file)
                print(path2+file)
            elif re.search("settings.txt", file):
                #setting_paths.append(path2 + file)
                settingPath = path2 + file
                print(settingPath)
            elif re.search("metadata.txt", file):
                metadataPath = path2 + file
                print(metadataPath)

# Define Nstar and Nso from the settings file
# Read in the metadata from the metadata.txt file
with open(metadataPath) as f:
    lines = f.readlines()
    for line in lines:
        print(line)
        exec(line)

# Read in the parameters from a file
try:
    fh = open(odeModelFolder + odeFile, 'r')
    print("Path to text file with ODEs and parms exists")
    with open(odeModelFolder + odeFile) as f:
        lines = f.readlines()
    startReading = False
    stopReading = False

    print("The system parameters:")
    for i, line in enumerate(lines):
        if "####START_PARAMETERS####\n" == lines[i - 1]:
            startReading = True
        elif "####END_PARAMETERS####\n" == line:
            stopReading = True

        if startReading and not stopReading:
            print(line)
            exec(line)
except FileNotFoundError:
    sys.exit("Path to text file with ODEs and parameters does not exist.")

# Read in the settings from the settings.txt file
with open(settingPath) as g:
    lines = g.readlines()
    for line in lines:
        print(line)
        exec(line)

# Define paths to data
allDataFiles = os.listdir(inputFolder)
# Select the parameter files
parFiles = []
for file in allDataFiles:
    # if str.split(file,".")[-1] == "csv":
    if re.search(".\d{2}.csv", file): #str.split(file, ".")[-1] == "csv":
        parFiles = parFiles + [file]
print("These files will be read: \n%s" % (',\n'.join([str(i) for i in parFiles])))

# Paths to output folder
pathToOutput = inputFolder

# Read in the ODEs from a file
try:
    fh = open(odeModelFolder + odeFile, 'r')
    print("Path to text file with ODEs exists")
    with open(odeModelFolder + odeFile) as f:
        lines = f.readlines()
    startReading = False
    stopReading = False
    # Make a dictionary that will contain the ODEs as specified in the file
    odeDict = {}
    obsDict = {}
    print("The system of ODES:")
    for i, line in enumerate(lines):
        if "####START_ODES####\n" == lines[i - 1]:
            startReading = True
        elif "####END_ODES####\n" == line:
            stopReading = True

        if startReading and not stopReading:
            print(line)
            exec(line)
            # Make a dictionary that contains the ODEs
            # It is a match if the first chars of line are f0, f1, ..., fn
            match = re.search('^f[0-9]+', line)
            if match:
                parts = line.split(" = ")
                odeDict[parts[0]] = parts[1][:-1]

            # Make a dictionary that contains the observables
            # It is a match if the first chars of line are g0, g1, ..., gn
            match = re.search('^g[0-9]+', line)
            if match:
                parts = line.split(" = ")
                obsDict[parts[0]] = parts[1][:-1]
except FileNotFoundError:
    sys.exit("Path to text file with ODEs does not exist.")

# Read in the data
data = pd.read_csv(dataFolder + dataFile)

# Make a separate dataTemp data frame per cell line, because every cell line can have a different number
# of replicates and save those data frames in a list:
dfListInterpol = []
dfListReal = []
# Loop over the reporters
for cell_linei, celllinet in enumerate(obsStateVarList):
    # Create NrOfDoses matrices with NrOfReplicates rows and NrOfTimepoints columns
    # For example: 22 time points, 3 replicates, 4 reporter cell lines, 5 doses
    dataTempInterpol = np.zeros((len(doseList), len(plateID_list[cell_linei]), Nt))
    dataTempReal = np.zeros((len(doseList), len(plateID_list[cell_linei]), Nt))

    # Loop over the doses
    for dosei, stresslevel in enumerate(doseList):

        # Loop over the replicates
        for replicatei, replicateID in enumerate(plateID_list[cell_linei]):

            # Loop over the time points
            for tpi in np.arange(0, Nt, 1):
                timeIDinData = tpi + 1
                datasample = data[(data.StateVar == celllinet) & (data.replID == replicateID) &
                                  (data.dose_uMadj == stresslevel) & (data.timeID == timeIDinData)]

                dataTempReal[dosei, replicatei, tpi] = datasample[realDataColName]
                dataTempInterpol[dosei, replicatei, tpi] = datasample[interpolDataColName]
    dfListInterpol.append(dataTempInterpol)
    dfListReal.append(dataTempReal)

dfListMeanReal = []
dfListStdReal = []
dfListMean = []
dfListStd = []
# Loop over the dataframes, i.e. the cell lines
for cell_linei, celllinet in enumerate(obsStateVarList):
    df = dfListInterpol[cell_linei]
    dfReal = dfListReal[cell_linei]

    dfTemp = np.array([]).reshape(len(plateID_list[cell_linei]), 0)
    dfTempReal = np.array([]).reshape(len(plateID_list[cell_linei]), 0)
    for dosei in np.arange(0, len(doseList), 1):
        # Bind the doses horizontally to obtain a NrOfReplicates x (NrOfDoses x NrOfTimepoints) matrix
        dfTemp = np.hstack((dfTemp, df[dosei]))
        dfTempReal = np.hstack((dfTempReal, dfReal[dosei]))

    # Calculate the mean and std per column, i.e. the mean and std of the replicates
    dfTempM = np.mean(dfTemp, axis=0)
    dfTempS = np.std(dfTemp, axis=0)
    dfTempMReal = np.mean(dfTempReal, axis=0)
    dfTempSReal = np.std(dfTempReal, axis=0)

    # Make a list of the means, where the length of the list is NrOfCelllines
    dfListMean.append(dfTempM)
    dfListStd.append(dfTempS)
    dfListMeanReal.append(dfTempMReal)
    dfListStdReal.append(dfTempSReal)

# Bind the mean lists as rows to form a matrix with dimensions NrOfCellines x (NrOfDoses x NrOfTimepoints)
dataM = np.vstack(dfListMean)
dataS = np.vstack(dfListStd)
dataMReal = np.vstack(dfListMeanReal)
dataSReal = np.vstack(dfListStdReal)

# Transpose the matrices to form a matrix with dimensions (NrOfDoses x NrOfTimepoints) x NrOfCelllines
dataMT = np.transpose(dataM)
dataST = np.transpose(dataS)

# Append all the rows to form a matrix with dimensions 1 x (NrOfCelllines x NrOfDoses x NrOfTimepoints)
dataMTR = dataMT.reshape(Nm * Nt * len(doseList), )
dataSTR = dataST.reshape(Nm * Nt * len(doseList), )

dataNoise = [dataMTR, dataSTR]

# Make an array of the time points as they are dictated by the data
dataTimePoints = data[(data.StateVar == obsStateVarList[0]) & (data.replID == plateID_list[0][0]) &
                      (data.dose_uMadj == doseList[0])][interpolTimeColName][0:Nt]

#########

idx = 1
if len(parmset_paths) == 1:
    path0 = parmset_paths[0]
    df = pd.read_csv(path0, index_col=0)
    r1 = df.iloc[:, [1]].T
    cost = re.findall("cost_\d+.\d{2}",path0)[0].split("_")[1]
    r1.rename(index={"est_value": cost}, inplace=True)
elif len(parmset_paths) > 1:
    path0 = parmset_paths[0]
    df = pd.read_csv(path0, index_col=0)
    r1 = df.iloc[:, [1]].T
    r1.rename(index={"est_value": 0}, inplace=True)
    cost = float(re.findall("cost_\d+.\d{2}",path0)[0].split("_")[1])
    c0 = pd.DataFrame({'Cost':[cost]})
    r1 = r1.join(c0)
    #print(r1)
    for path in parmset_paths[1:]:
        print(path)
        df = pd.read_csv(path, index_col=0)
        r2 = df.iloc[:, [1]].T
        cost = float(re.findall("cost_\d+.\d{2}", path)[0].split("_")[1])
        r2.rename(index={"est_value": 0}, inplace=True)
        c0 = pd.DataFrame({'Cost': [cost]})
        r2 = r2.join(c0)
        #print(r2)
        r1 = r1.append(r2)
    r1.sort_values(by = ["Cost"], inplace = True)

else:
    print("No files found.")

r1.set_index("Cost", inplace=True)

print(r1)

# Save the cost overview
r1.to_csv(outPath + "ParameterOverview.csv")

# Convert the row names to floats
r1.index = r1.index.astype("float")

# Make a plot with all the costs
plt.scatter(range(1,r1.shape[0]+1), r1.index, color = 'k')
plt.plot(range(1,r1.shape[0]+1), r1.index, color = 'k')
plt.xlabel('Index')
plt.ylabel('Cost')
plt.tight_layout()
plt.savefig(outPath + "CostFigure.pdf")
plt.close()

# Get unique costs
r1.index = r1.index.astype("int")
rownames = r1.index.unique().tolist()

# Make a subset of the r1 matrix for the boxplot
r1_withSO = r1#abs(r1) # r1.iloc[:,np.concatenate([np.arange(0,Np-Nso),np.arange(Np,Np+Nstar)])]
r1_withoutSO = r1.iloc[:,np.concatenate([np.arange(0,Np-Nso),np.arange(Np,Np+Nstar)])]
r1_SO = r1.iloc[:,Np-Nso:Np]
#r1_noStarSO = r1.iloc[:,:-(Nstar+Nso)]

# Make boxplot for the parameter estimates of the first three lowest costs
dfList = []
for i in range(0,len(rownames)):
    sText1 = "dfSub" + str(i) + " = r1_SO.loc[rownames[i],:]"
    exec(sText1)
    sText2 = "dfList.append(dfSub" + str(i) + ")"
    exec(sText2)
# dfSub1 = r1_noSO.loc[rownames[i],:]
# dfSub2 = r1_noSO.loc[rownames[1],:]
# dfSub3 = r1_noSO.loc[rownames[2],:]
# dfList = [dfSub1, dfSub2, dfSub3]

order = ["1_lowest","2_lower","3_low"] + [str(j) for j in range(4,len(dfList)+1)]
for i,dfSub in enumerate(dfList):
    if len(dfSub.shape) > 1:
        sns.set_style("whitegrid")
        ax = sns.boxplot(data = dfSub, orient = "v", width = 0.3, palette = "viridis",fliersize=2)
        ax.xaxis.grid(True, color = "lightgrey")
        ax.yaxis.grid(False)
        plt.xlabel('Parameter')
        plt.xticks(rotation='vertical')
        plt.ylabel('Estimated value')
        plt.title(str(dfSub.shape[0]) + ' parameter sets')
        plt.tight_layout()
        plt.savefig(outPath + "Figure" + "_" + order[i] + ".pdf")
        plt.close()

# Make boxplot for the parameter estimates of the first three lowest costs
dfList = []
for i in range(0,len(rownames)):
    sText1 = "dfSub" + str(i) + " = r1_withoutSO.loc[rownames[i],:]"
    exec(sText1)
    sText2 = "dfList.append(dfSub" + str(i) + ")"
    exec(sText2)

order = ["1_lowest","2_lower","3_low"] + [str(j) for j in range(4,len(dfList)+1)]
for i,dfSub in enumerate(dfList):
    if len(dfSub.shape) > 1:
        sns.set_style("whitegrid")
        ax = sns.boxplot(data = dfSub, orient = "v", palette = "viridis", width = 0.5, fliersize=2)
        ax.xaxis.grid(True, color = "lightgrey")
        ax.yaxis.grid(False)
        plt.yscale("log")
        plt.xlabel('Parameter')
        plt.xticks(rotation='vertical')
        plt.ylabel('Estimated value')
        plt.title(str(dfSub.shape[0]) + ' parameter sets')
        plt.tight_layout()
        plt.savefig(outPath + "Figure_Log" + "_" + order[i] + ".pdf")
        plt.close()

# Make a color palette
#palette = sns.color_palette(None, Nm)

def replaceStateVars(stateList, d, string=""):
    dictExec = d
    # Replace the state variables by z[0], z[1], ..., z[n]
    for i, stateVar in enumerate(stateList):
        result = []
        for j, key in enumerate(dictExec.keys()):
            splitres = re.split(r"([^a-zA-Z]" + stateVar + "[^a-zA-Z])", dictExec[key])
            subres = []
            for res in splitres:
                if re.search(r"((^|[^a-zA-Z])" + stateVar + "([^a-zA-Z]|$))", res):
                    subres = subres + [re.sub(stateVar, string + str(i) + "]", res)]
                else:
                    subres = subres + [res]
            joinedres = "".join(subres)
            dictExec[key] = joinedres
    return dictExec


def modelSimulation(z, t, args):
    [odeDictExec, Stressautowrap, p] = args

    # Stress level
    sText = "Stressautowrap("
    for i, var in enumerate(p):
        if i < len(p) - 1:
            sText = sText + str(var) + ","
        else:
            sText = sText + str(var) + ")"
    S = eval(sText)
    # print("Executed: S = %s" % sText)

    # Evaluate the ODEs
    output = []
    for i, key in enumerate(sorted(odeDictExec.keys())):
        output.append(eval(odeDictExec[key]))

    # Output
    dzdt = output
    return dzdt

# Define t
t = np.array(tspan)

# Make simulations with the parameter estimates of dfSub4
allSimus = []
dfSub4 = r1.loc[rownames[0],:]

try:
    for row in range(0,dfSub4.shape[0]):
        pAllNew = np.array(dfSub4.iloc[row,:])

        # Make the parameters global variables
        for i,par in enumerate(np.array(dfSub4.columns)):
            sText = par + " = " + str(pAllNew[i])
            # print(sText)
            exec(sText)

        # initial condition
        # Make a global variable z that contains the initial conditions, for example:
        sText = "z = ["
        # allIniStates = stateOList+stateOKnownNameList
        allIniStates = allIniStatesList
        for stateVari in np.arange(0,len(allIniStates),1):
            if stateVari < len(allIniStates) - 1:
                sText = sText + allIniStates[stateVari] + ","
            else:
                sText = sText + allIniStates[stateVari] + "]"
        exec(sText)
        # print("Executed: %s" % sText)

        # Replace the state variables in odeDict by z[0], z[1], ..., z[n]
        odeDictExec = replaceStateVars(stateList, odeDict, string="z[")
        # Replace the state variables in odeDict by z[0], z[1], ..., z[n]
        obsDictExec = replaceStateVars(stateList, obsDict, string="simu[:,")

        # Simulate the model with the estimated parameters (simulations) and transform the output
        # by using the mapping g (mappedSimu)
        simulations = np.zeros((len(doseList), len(t), len(stateList)))
        mappedSimu = np.zeros((len(doseList), len(t), len(obsList)))
        for i, dose in enumerate(doseList):  # enumerate(outputParms[Ns:Ns+Nd]):
            # utargs = outputParms[Ns:Ns+Nd].tolist() + ["t"] + [doseList[i]] + outputParms[Ns+Nd:Ns+Nd+Npk].tolist()
            # states = pAllNew[0:Ns].tolist()
            if Nstar == 0:
                p = pAllNew[Nse:].tolist() + ["t"] + [dose]
            else:
                p = pAllNew[-Nstar:].tolist() + pAllNew[Nse:-Nstar].tolist() + ["t"] + [dose]
            # p = ["t"] + [doseList[i]]
            args = [odeDictExec, Stressautowrap, p]
            simu = odeint(modelSimulation, z, t, args=(args,))
            simulations[i, :, :] = simu
            if len(obsList) > 1:
                sText = "mappedSimu[i,:,:] = np.transpose(np.array(["
                for j, key in enumerate(sorted(obsDictExec.keys())):
                    if j < len(obsList) - 1:
                        sText = sText + obsDictExec[key] + ","
                    else:
                        sText = sText + obsDictExec[key] + "]))"
                exec(sText)
                # print("Executed: %s" % sText)
            else:
                sText = "mappedSimu[i,:,0] = " + list(obsDictExec.values())[0]
                exec(sText)
                # print("Executed: %s" % sText)
        # print(mappedSimu)
        allSimus.append(mappedSimu)
except:
    sys.exit("There is only one best fit. Exiting script.")
    #pass

meanSimu = np.mean(allSimus, axis=0)

# Make a plot that plots the best fits and the mean of those fits
# Plot the simulations of the state variables together with the data mean
fig, ax = plt.subplots(Nm, len(doseList), sharex='col', sharey='row')
if len(doseList) > 1:
    plot_titles = [(str(i) + " uM") for i in doseList]
else:
    plot_titles = str(doseList[0]) + " uM"

begin = 0
if len(doseList) > 1:
    for i, dose in enumerate(doseList):
        for j, stateVar in enumerate(obsStateVarList):
            # Plot the data
            ax[j, i].plot(dataTimePoints, dataMReal[j, begin:begin + Nt], 'o', label=stateVar,
                          color=palette[j])
            # Plot the errorbars around the data
            ax[j, i].errorbar(dataTimePoints, dataMReal[j, begin:begin + Nt], yerr=dataSReal[j, begin:begin + Nt],
                              capsize=4, fmt="none",
                              color=palette[j])
            # Plot the model simulation
            for mappedsimu in allSimus:
                ax[j, i].plot(dataTimePoints, mappedsimu[i, :, j], label=stateVar, color="#D3D3D3", linewidth=3)
            ax[j, i].plot(dataTimePoints, meanSimu[i, :, j], '--',label=stateVar, color="#000000", linewidth=3)
            ax[j, Nd - 1].set_ylabel(stateVar, rotation=-90, fontsize=24)
            ax[j, Nd - 1].yaxis.set_label_coords(1.03, 0.5)
            ax[j, i].tick_params(labelsize=16)
        begin += Nt
        ax[0, i].set_title(plot_titles[i],
                           fontdict={'fontsize': 24})
else:
    for j, stateVar in enumerate(obsStateVarList):
        # Plot the data
        ax[j].plot(dataTimePoints, dataMReal[j, begin:begin + Nt], 'o', label=stateVar,
                   color=palette[j])
        # Plot the errorbars around the data
        ax[j].errorbar(dataTimePoints, dataMReal[j, begin:begin + Nt], yerr=dataSReal[j, begin:begin + Nt],
                       capsize=4, fmt="none",
                       color=palette[j])
        # Plot the model simulation
        ax[j].plot(dataTimePoints, mappedsimu[0, :, j], label=stateVar, color="#000000", linewidth=3)
        ax[j].set_ylabel(stateVar, rotation=-90, fontsize=24)
        ax[j].yaxis.set_label_coords(1.03, 0.5)
        ax[j].tick_params(labelsize=16)
    begin += 43
    ax[j].set_title(plot_titles,
                    fontdict={'fontsize': 24})

# add a big axes, hide frame
fig.add_subplot(111, frameon=False)
fig.subplots_adjust(hspace=0.2)
fig.set_figheight(12)
fig.set_figwidth(10)
fig.suptitle('Model vs. Data', fontsize=28)
fig.text(0.02, 0.5, 'Intensity (a.u.)', va='center', rotation='vertical', fontsize=28)

# hide tick and tick label of the big axes
plt.grid(False)
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.xlabel("Time (h)", fontsize=28)
plt.savefig(outPath + MODEL + "_" + TRMNT +"_PlotOf" + str(dfSub4.shape[0]) + "BestModels.pdf")
plt.close()

print("Plotting done! Script ended.")

# # Check which treatments are in the inputFolder
# treatments = os.listdir(inputPath + MODEL)
# for trmnt in treatments:
#     treatmentInputFolder = inputPath + MODEL + trmnt
#
# # Read in the metadata from the metadata.txt file
# with open(inputFolder + "metadata.txt") as f:
#     lines = f.readlines()
#     for line in lines:
#         print(line)
#         exec(line)
#
# # Read in the parameters from a file
# try:
#     fh = open(odeModelFolder + odeFile, 'r')
#     print("Path to text file with ODEs and parms exists")
#     with open(odeModelFolder + odeFile) as f:
#         lines = f.readlines()
#     startReading = False
#     stopReading = False
#
#     print("The system parameters:")
#     for i, line in enumerate(lines):
#         if "####START_PARAMETERS####\n" == lines[i - 1]:
#             startReading = True
#         elif "####END_PARAMETERS####\n" == line:
#             stopReading = True
#
#         if startReading and not stopReading:
#             print(line)
#             exec(line)
# except FileNotFoundError:
#     sys.exit("Path to text file with ODEs and parameters does not exist.")
#
# # Read the settings from a file
# with open(path2 + "settings.txt") as g:
#     lines = g.readlines()
#     for line in lines:
#         print(line)
#         exec(line)

# from sklearn.feature_extraction.text import TfidfVectorizer
#
# documents = [open(f) for f in text_files]
# tfidf = TfidfVectorizer().fit_transform(documents)
# # no need to normalize, since Vectorizer will return normalized tf-idf
# pairwise_similarity = tfidf * tfidf.T