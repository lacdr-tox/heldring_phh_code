import matplotlib
import numpy as np
from sklearn.linear_model import LinearRegression
import pandas as pd
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib.axes as axes
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
from sympy import symbols
from sympy import Piecewise
from sympy import exp
from sympy.utilities.autowrap import autowrap
import re
import time
import multiprocessing as mp
import sys
from itertools import chain
from scipy.optimize import fsolve

###################################################
# Part I - USER DEFINED VARIABLES
###################################################

#inputFolder = "/data/muriel/Projects/PHH/Modeling/Output/M033/CDDP/20220326_110349/"
inputFolder = "/data/muriel/Projects/PHH/Modeling/Output/M033/CDDP/20220326_110349/" #  # Write down the path to the folder that
# contains the parameter estimates of all the finished runs
#inputFolder = "/Users/muriel/Documents/LACDR/Projects/DDP/Modeling/02_Output/withRNA/M029/CDDP/20210329/Run1/" #

# Define the path to the ode system file
odeModelFolder = "/home/muriel/Modeling/Input/withRNA/" # Write down the path to the folder that
# contains the file with the ODE system
#odeModelFolder = "/Users/muriel/Documents/LACDR/Projects/DDP/Modeling/01_Input/withRNA/" #

# Define output folder
outputFolder = "/data/muriel/Projects/PHH/Modeling/Output/VaryParmsData/"
#outputFolder = "/Users/muriel/Documents/LACDR/Projects/DDP/Modeling/02_Output/VirtualData/" #

# Choose one file that contains the parameter estimates that need to be altered
#fileName = "20220326_110349_MH_M033Model_parameterEstimates_parmset13_cost_31.37.csv"
fileName = "20220326_110349_MH_M033Model_parameterEstimates_parmset13_cost_31.37.csv" # "20220314_111222_MH_M019Model_parameterEstimates_parmset3_cost_80.80.csv" 
#"20210406_MH_M031Model_parameterEstimates_parmset43_cost_81.59.csv" #"20210329_MH_M029Model_parameterEstimates_parmset1_cost_81.98.csv" # "20201114_MH_M022Model_parameterEstimates_parmset4_cost_81.59.csv"

# Overwrite date
CURDATE = "20220331"
OUTPUTNAME = "MDM2synth"
SDVAR = [0.2] # [0.001,0.01] #,0.1] # [0.2]
RUNS = 1000
SAMPLES = 50
CORES = 50
withRNA = True
CONC = 3.3
Nutlin = True
rlist = [0.1,10]
fixed_states = [[0, ["DD_init", 0, 0]],
                [1, ["P53rna_init", 1, 1]], # [position in Ns, [Name, init_value, est_value]]
                [2, ["P53_init", 1, 1]],
                [3, ["P53P_init", 0, 0]],
                [4, ["MDM2rna_init", 1, 1]],
                [5, ["MDM2_init", 1, 1]],
                [6, ["MDM2S395_init", 0, 0]],
                [7, ["C_init", 0, 0]],
                [8, ["P21rna_init", 1, 1]],
                [9, ["P21_init", 1, 1]],
                [10, ["BTG2rna_init", 1, 1]],
                [11, ["BTG2_init", 1, 1]],
                [15, ["kd_p53rna", 1, 1]],
                [16, ["kd_btg2", 1, 1]],
                [17, ["kd_p21", 1, 1]],
                [18, ["kd_mdm2", 1, 1]],
                [19, ["kd_p53", 1, 1]]]
fixed_parms = ["kd_p53rna", "kd_btg2", "kd_p21", "kd_mdm2","kd_p53"]

# fixed_states = [[1, ["P53rna_init", 1, 1]], # [position in Ns, [Name, init_value, est_value]]
#                 [2, ["P53_init", 1, 1]],
#                 [4, ["MDM2rna_init", 1, 1]],
#                 [5, ["MDM2_init", 1, 1]],
#                 [6, ["P21rna_init", 1, 1]],
#                 [7, ["P21_init", 1, 1]],
#                 [8, ["BTG2rna_init", 1, 1]],
#                 [9, ["BTG2_init", 1, 1]],
#                 [13, ["kd_p53rna", 1, 1]],
#                 [14, ["kd_btg2", 1, 1]],
#                 [15, ["kd_p21", 1, 1]],
#                 [16, ["kd_mdm2", 1, 1]]]
# fixed_parms = ["kd_p53rna", "kd_btg2", "kd_p21", "kd_mdm2"]

# Select the time point
tp8hr = 8/1.5
tp24hr = 24/1.5

# Make colnames
ROWNAMES_8hr = np.array(["DD_8hr","P53rna_8hr","P53_8hr","P53P_8hr","MDM2rna_8hr","MDM2_8hr","MDM2S395_8hr","C_8hr","P21rna_8hr","P21_8hr","BTG2rna_8hr","BTG2_8hr"])
ROWNAMES_24hr = np.array(["DD_24hr","P53rna_24hr","P53_24hr","P53P_24hr","MDM2rna_24hr","MDM2_24hr","MDM2S395_24hr","C_24hr","P21rna_24hr","P21_24hr","BTG2rna_24hr","BTG2_24hr"])
#ROWNAMES_8hr = np.array(["DD_8hr","P53rna_8hr","P53_8hr","P53P_8hr","MDM2rna_8hr","MDM2_8hr","P21rna_8hr","P21_8hr","BTG2rna_8hr","BTG2_8hr"])
#ROWNAMES_24hr = np.array(["DD_24hr","P53rna_24hr","P53_24hr","P53P_24hr","MDM2rna_24hr","MDM2_24hr","P21rna_24hr","P21_24hr","BTG2rna_24hr","BTG2_24hr"])


# Output name
outputFactors1 = list(fileName.split("_")[i] for i in [0,1,3])
outputFactors2 = [OUTPUTNAME] #[inputFolder.split("/")[-2]]
outputFactors3 = [fileName.split("_")[5]]
outputName = "varyR_" + "_".join(outputFactors2 + outputFactors1 + outputFactors3) # + outputFactors3)
print("This is the output file name: %s" % outputName)

###################################################
# Part II - Define Functions
###################################################

def interpolConc(doseList, doseNames, concentration):
    # Make an array of the applied concentrations
    x = np.array(doseList).reshape((-1, 1))

    # Make a list with the doses
    fittedDose = []
    for dn in doseNames:
        fittedDose.append(eval(dn))
    y = np.array(fittedDose)

    # Interpolate the stress levels
    reg = LinearRegression().fit(x, y)
    ECint = float(reg.predict(np.array([concentration]).reshape((-1, 1))))

    # Append the interpolated stress level
    ynew = np.append(y, ECint)
    ynew = list(np.sort(ynew))

    return ynew


def makeNewRows(df, parsToVary, parsToKeep, samples, std):
    # new_rows = np.empty([0,samples*len(sd)*runs])
    np.random.seed()
    new_rows = []
    for index, row in df.iterrows():
        new_values = []
        # for run in range(1,runs+1):
        for proportionalSD in std:
            if row.ParName in parsToVary:
                sd = row.est_value * proportionalSD
                new_values = new_values + list(abs(row.est_value + np.random.normal(0, sd, samples)))
            elif row.ParName in parsToKeep:
                new_values = new_values + [row.est_value for i in range(samples)]
            else:
                new_values = new_values + [row.est_value for i in range(samples)]
        new_rows.append(new_values)
    return (new_rows)


# Stress level decrease function instead of Stressautowrap
def sFunction(s_init, r, t):
    out = s_init * np.exp(-r * t)
    return out


def simulateModel(stressLevel, t, parmsNamesList, parmsList, fromTp1, odeDict, colname, idc):
    # Set the initial states and parameters
    if fromTp1:
        # parms = parmsList[Ns+Npk+Nd:]
        z = parmsList[0:Ns]
    else:
        # parms = parmsList[Ns+Npk+Nd:-Ns]
        z = parmsList[-Ns::]

    # Set the pkparms
    tau1 = parmsList[Ns + Nd]

    # solve ODE
    dzdt = odeint(modelSimulation, z, t, args=([sFunction, tau1, stressLevel, parmsList, odeDict],))
    
    if (idc+1)%500 == 0:
        # Print when done
        print("Done with column %s (parmset %d out of %d)" % (colname,(idc+1), RUNS*SAMPLES*len(SDVAR)))

    return (dzdt[-1])


def replaceStateVars(stateList, d, string=""):
    dictExec = d
    # Replace the state variables by z[0], z[1], ..., z[n]
    for i, stateVar in enumerate(stateList):
        result = []
        for j, key in enumerate(dictExec.keys()):
            splitres = re.split(r"([^a-zA-Z]" + stateVar + "[^a-zA-Z])", dictExec[key])
            # print(splitres)
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
    [sFunction, r, s_init, p, odeDict] = args
    # [odeDictExec, sFunction, p] = args

    S = sFunction(s_init, r, t)

    # Get the correct order of f equations:
    ordr = []
    for i in range(0, len(stateList)):
        ordr.append('f' + str(i))

    # Evaluate the ODEs
    output = []
    for i, key in enumerate(ordr):
        output.append(eval(odeDict[key]))

    # Output
    dzdt = output
    return dzdt


def replaceParameters(parList, d, string=""):
    dictExec = d
    # Replace the parameters by p[0], p[1], ..., p[n]
    for i, par in enumerate(parList):
        result = []
        for j, key in enumerate(dictExec.keys()):
            splitres = re.split(r"((?:^|[^a-zA-Z_])" + par + "[^a-zA-Z_])", dictExec[key])
            subres = []
            for res in splitres:
                if re.search(r"((^|[^a-zA-Z&^_])" + par + "([^a-zA-Z&^_]|$))", res):
                    subres = subres + [re.sub(par, string + str(i) + "]", res)]
                else:
                    subres = subres + [res]
            joinedres = "".join(subres)
            dictExec[key] = joinedres
    return dictExec

def calcEquil(stressLevel, parmsNamesList, parmsList, fromTp1, odeDict):

    # Set the initial states and parameters
    if fromTp1:
        z = parmsList[0:Ns]
    else:
        z = parmsList[-Ns::]
                                                 
    # Set the pkparms
    tau1 = parmsList[Ns+Nd]
                                                                   
    # calc equilibrium
    eq = fsolve(modelSimulation, z, args=(0, [sFunction, tau1, stressLevel, parmsList, odeDict]))
                                                                                           
    return(eq)


# Function to insert row in the dataframe
def Insert_row(row_number, df, row_value):
    # Starting value of upper half
    start_upper = 0

    # End value of upper half
    end_upper = row_number

    # Start value of lower half
    start_lower = row_number

    # End value of lower half
    end_lower = df.shape[0]

    # Create a list of upper_half index
    upper_half = [*range(start_upper, end_upper, 1)]

    # Create a list of lower_half index
    lower_half = [*range(start_lower, end_lower, 1)]

    # Increment the value of lower half by 1
    lower_half = [x.__add__(1) for x in lower_half]

    # Combine the two lists
    index_ = upper_half + lower_half

    # Update the index of the dataframe
    df.index = index_

    # Insert a row at the end
    df.loc[row_number] = row_value

    # Sort the index labels
    df = df.sort_index()

    # return the dataframe
    return df

###################################################
# Part III - Read data from files
###################################################

parmsDFlist = []
expressionDFlist = []

for r in rlist:

    # Get the date, treatment and model name from the path
    date = str.split(inputFolder,"/")[-2]
    trmnt = str.split(inputFolder,"/")[-3]
    model = str.split(inputFolder,"/")[-4]

    # Read in the metadata from the metadata.txt file
    with open(inputFolder + "metadata.txt") as f:
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
    with open(inputFolder + "settings.txt") as g:
        lines = g.readlines()
        for line in lines:
            print(line)
            exec(line)

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

    ###################################################
    # Part IV - Vary parms
    ###################################################

    # Read in the file with the parameter estimates
    file = inputFolder+fileName

    # Define the parameters
    df_parmsOrig = pd.read_csv(file)

    # Manually add fixed parameters in their appropriate positions
    # df_parmsComplete = Insert_row(1, df_parmsOrig, ["P53rna_init", 1, 1])
    # df_parmsComplete = Insert_row(2, df_parmsComplete, ["P53_init", 1, 1])
    if fixed_states:
        for fs in fixed_states:
            df_parmsOrig = Insert_row(fs[0], df_parmsOrig, fs[1])

    # Vary the star parameters and snow parameters without the scaling and offset
    # Keep kd_dd and ks_dd the same
    parsToVary = paraStarList + fixed_parms + paraSnowList[0:Nsnow - Nso]
    parsToKeep = doseParmsList + PKparms + paraSnowList[Nsnow-Nso:]

    # Rename parameterEst data frame
    df_parmsOrig = df_parmsOrig.rename(columns={"Unnamed: 0": "ParName"})
    #print(', '.join(map(str, list(df_parmsOrig.loc[:,"ParName"]))))

    if Nutlin:
        # Change the Mdm2-dependent degradation with factor r
        #df_parmsOrig.loc[df_parmsOrig["ParName"] == "kd_p53_mdm2","est_value"] = r * df_parmsOrig.loc[df_parmsOrig["ParName"] == "kd_p53_mdm2","est_value"]
        #df_parmsOrig.loc[df_parmsOrig["ParName"] == "kd_p53p_mdm2","est_value"] = r * df_parmsOrig.loc[df_parmsOrig["ParName"] == "kd_p53p_mdm2","est_value"]
        #df_parmsOrig.loc[df_parmsOrig["ParName"] == "k_dp","est_value"] = r * df_parmsOrig.loc[df_parmsOrig["ParName"] == "k_dp","est_value"]
        df_parmsOrig.loc[df_parmsOrig["ParName"] == "ks_mdm2_p53p","est_value"] = r * df_parmsOrig.loc[df_parmsOrig["ParName"] == "ks_mdm2_p53p","est_value"]

    # Get the stress levels
    # Make the concentrations global variables
    doseDF = df_parmsOrig[Ns:Ns+Nd]
    for i in doseDF.iterrows():
        sText = i[1].ParName + " = " + str(i[1].est_value)
        print(sText)
        exec(sText)

    # Check whether 0 is in the dose list. If not, add 0
    if 0 in doseList:
        doseListNew = doseList
        doseNamesNew = doseNames
        intDoses = interpolConc(doseList,doseNames,CONC)
    else:
        EC0 = 0
        doseListNew = [0] + doseList
        doseNamesNew = ["EC0"] + doseNames
        intDoses = interpolConc(doseListNew,doseNamesNew,CONC)

    print("These are the stress levels: %s" % (', '.join(map(str, intDoses))))

    doseDFnew = pd.DataFrame({"Applied": list(np.sort(np.append(doseListNew,CONC))),
                             "Effective": intDoses})

    # Get the concentration of interest
    concOfInterest = float(doseDFnew.loc[doseDFnew["Applied"] == CONC,"Effective"])
    print("This is the stress level of interest: %f" % concOfInterest)

    # Make new columns with new parameter values
    numberOfRuns = RUNS
    numberOfSamples = SAMPLES
    proportionalSD = SDVAR
    newColNames = [["r" + str(k+1) + "_newValue" + str(i+1) + "_v" + str(j) + "_rfactor" + str(r) for j in SDVAR for i in range(numberOfSamples) ] for k in range(numberOfRuns)  ]

    ### Do parameter perturbation ###
    print("Perturbing parameters...")

    # perturbations = [makeNewRows(df_parmsOrig,parsToVary, parsToKeep, numberOfSamples, SDVAR)]
    # perturbations.append(makeNewRows(df_parmsOrig,parsToVary, parsToKeep, numberOfSamples, SDVAR))

    pool = mp.Pool(CORES)
    perturbations = [pool.apply(makeNewRows, args = (df_parmsOrig,parsToVary, parsToKeep, numberOfSamples, SDVAR)) for run in range(1,numberOfRuns+1)]
    pool.close()

    print("Done perturbing parameters...")

    # Bind the perturbations
    df_parmsPHH = df_parmsOrig.reset_index(drop=True)
    for run in range(numberOfRuns):
        df_tmp = pd.DataFrame.from_records(perturbations[run], columns = newColNames[run])
        df_parmsPHH = pd.concat([df_parmsPHH.reset_index(drop=True), df_tmp], axis=1)

    ###################################################
    # Part V - Simulate until steady state
    ###################################################

    # Replace the state variables in odeDict by z[0], z[1], ..., z[n]
    # This is necessary because they are changing during simulation
    replaceStateVars(stateList, odeDict, string="z[");
    # Replace the state variables in odeDict by simu[0], simu[1], ..., simu[n]
    replaceStateVars(stateList, obsDict, string="simu[:,");

    # Replace the parameters in odeDict by p[i], p[i+1], ..., p[n]
    replaceParameters(list(df_parmsPHH.loc[:,"ParName"]), odeDict, string="p[");

    ### Do simulation until steady state ###
    print("Simulating until steady state...")
    time_start = time.time()

    # time points
    fromTp1 = True

    # sols = [calcEquil(0, list(df_parmsPHH.loc[:,"ParName"]), list(df_parmsPHH.loc[:,list(chain.from_iterable(newColNames))[0]]), fromTp1, odeDict)]
    # sols.append(calcEquil(0, list(df_parmsPHH.loc[:,"ParName"]), list(df_parmsPHH.loc[:,list(chain.from_iterable(newColNames))[1]]), fromTp1, odeDict))
    # sols.append(calcEquil(0, list(df_parmsPHH.loc[:,"ParName"]), list(df_parmsPHH.loc[:,list(chain.from_iterable(newColNames))[2]]), fromTp1, odeDict))
    # sols.append(calcEquil(0, list(df_parmsPHH.loc[:,"ParName"]), list(df_parmsPHH.loc[:,list(chain.from_iterable(newColNames))[3]]), fromTp1, odeDict))
    # sols.append(calcEquil(0, list(df_parmsPHH.loc[:,"ParName"]), list(df_parmsPHH.loc[:,list(chain.from_iterable(newColNames))[4]]), fromTp1, odeDict))
    # sols.append(calcEquil(0, list(df_parmsPHH.loc[:,"ParName"]), list(df_parmsPHH.loc[:,list(chain.from_iterable(newColNames))[5]]), fromTp1, odeDict))
    # sols.append(calcEquil(0, list(df_parmsPHH.loc[:,"ParName"]), list(df_parmsPHH.loc[:,list(chain.from_iterable(newColNames))[6]]), fromTp1, odeDict))
    # sols.append(calcEquil(0, list(df_parmsPHH.loc[:,"ParName"]), list(df_parmsPHH.loc[:,list(chain.from_iterable(newColNames))[7]]), fromTp1, odeDict))

    # Simulate ...
    pool = mp.Pool(CORES)
    sols = [pool.apply(calcEquil, args = (0, list(df_parmsPHH.loc[:,"ParName"]), list(df_parmsPHH.loc[:,colName]), fromTp1, odeDict)) for colName in list(chain.from_iterable(newColNames))]
    pool.close()

    time_end = time.time()
    print("Spent %f seconds on simulating %i parameter sets" % ((time_end - time_start), (RUNS*SAMPLES*len(SDVAR))))

    ###################################################
    # Part VI - Simulate with stress for first 8 hours
    ###################################################

    print("Replacing state variable values with steady states...")

    # Replace old steady state values with the new ones
    df_parmsPHH_med = df_parmsPHH.copy()
    for i, colName in enumerate(list(chain.from_iterable(newColNames))):
        for j, sol in enumerate(sols[i]):
            df_parmsPHH_med.iloc[j, i + 3] = sol

    parmsDFlist.append(df_parmsPHH_med)

    dfSubset = df_parmsPHH_med.iloc[0:Ns,[0]+list(range(3, df_parmsPHH_med.shape[1]))]
    dfSubset.set_index('ParName', inplace=True)
    dfSubsetMedium = dfSubset.T

    print("Simulating with stress for first 8 hours...")

    # Get the steady states after simulation
    # time points
    t1 = np.linspace(1, tp8hr, int(tp8hr*2))
    t2 = np.linspace(tp8hr, tp24hr, int((tp24hr-tp8hr)*2))
    #print(t1)
    #print(t2)

    time_start = time.time()
    fromTp1 = True

    # sols_8hr = [simulateModel(concOfInterest, t1, list(df_parmsPHH_med.loc[:,"ParName"]), list(df_parmsPHH_med.loc[:,list(chain.from_iterable(newColNames))[0]]),fromTp1,odeDict)]
    # sols_8hr.append(simulateModel(concOfInterest, t1, list(df_parmsPHH_med.loc[:,"ParName"]), list(df_parmsPHH_med.loc[:,list(chain.from_iterable(newColNames))[1]]),fromTp1,odeDict))
    # sols_8hr.append(simulateModel(concOfInterest, t1, list(df_parmsPHH_med.loc[:,"ParName"]), list(df_parmsPHH_med.loc[:,list(chain.from_iterable(newColNames))[2]]),fromTp1,odeDict))
    # sols_8hr.append(simulateModel(concOfInterest, t1, list(df_parmsPHH_med.loc[:,"ParName"]), list(df_parmsPHH_med.loc[:,list(chain.from_iterable(newColNames))[3]]),fromTp1,odeDict))
    # sols_8hr.append(simulateModel(concOfInterest, t1, list(df_parmsPHH_med.loc[:,"ParName"]), list(df_parmsPHH_med.loc[:,list(chain.from_iterable(newColNames))[4]]),fromTp1,odeDict))
    # sols_8hr.append(simulateModel(concOfInterest, t1, list(df_parmsPHH_med.loc[:,"ParName"]), list(df_parmsPHH_med.loc[:,list(chain.from_iterable(newColNames))[5]]),fromTp1,odeDict))
    # sols_8hr.append(simulateModel(concOfInterest, t1, list(df_parmsPHH_med.loc[:,"ParName"]), list(df_parmsPHH_med.loc[:,list(chain.from_iterable(newColNames))[6]]),fromTp1,odeDict))
    # sols_8hr.append(simulateModel(concOfInterest, t1, list(df_parmsPHH_med.loc[:,"ParName"]), list(df_parmsPHH_med.loc[:,list(chain.from_iterable(newColNames))[7]]),fromTp1,odeDict))

    pool = mp.Pool(CORES)
    sols_8hr = [pool.apply(simulateModel, args = (concOfInterest, t1, list(df_parmsPHH_med.loc[:,"ParName"]), list(df_parmsPHH_med.loc[:,colName]),fromTp1,odeDict,colName,idc)) for idc,colName in enumerate(list(chain.from_iterable(newColNames)))]
    pool.close()

    time_end = time.time()
    print("Spent %f seconds on simulating %i parameter sets" % ((time_end - time_start), (RUNS*SAMPLES*len(SDVAR))))

    # Make columns for the new df
    df_parmsPHH_cpt5 = df_parmsPHH_med
    # 8hr rows
    newNames = ROWNAMES_8hr

    newCols_8hr = [newNames, np.array([np.NaN]*len(newNames)), np.array([np.NaN]*len(newNames))] + sols_8hr

    # Make a data frame of the columns
    df_cpt5_8hr = pd.DataFrame([])
    for i, cName in enumerate(list(df_parmsPHH_cpt5.columns)):
        df_cpt5_8hr[cName] = newCols_8hr[i]

    # Append new rows with the values at 8 hr
    df_parmsPHH_cpt5 = df_parmsPHH_cpt5.append(df_cpt5_8hr)
    df_parmsPHH_cpt5.iloc[list(range(0,Ns))+list(range(len(df_parmsPHH_cpt5.index)-Ns,len(df_parmsPHH_cpt5.index)))]

    ###################################################
    # Part VII - Simulate with stress for 24 hours
    ###################################################
    print("Simulating with stress for 24 hours...")
    # Make 24 hr rows
    fromTp1 = False

    # sols_24hr = [simulateModel(concOfInterest, t2, list(df_parmsPHH_cpt5.loc[:,"ParName"]), list(df_parmsPHH_cpt5.loc[:,list(chain.from_iterable(newColNames))[0]]),fromTp1,odeDict)]
    # sols_24hr.append(simulateModel(concOfInterest, t2, list(df_parmsPHH_cpt5.loc[:,"ParName"]), list(df_parmsPHH_cpt5.loc[:,list(chain.from_iterable(newColNames))[1]]),fromTp1,odeDict))
    # sols_24hr.append(simulateModel(concOfInterest, t2, list(df_parmsPHH_cpt5.loc[:,"ParName"]), list(df_parmsPHH_cpt5.loc[:,list(chain.from_iterable(newColNames))[2]]),fromTp1,odeDict))
    # sols_24hr.append(simulateModel(concOfInterest, t2, list(df_parmsPHH_cpt5.loc[:,"ParName"]), list(df_parmsPHH_cpt5.loc[:,list(chain.from_iterable(newColNames))[3]]),fromTp1,odeDict))
    # sols_24hr.append(simulateModel(concOfInterest, t2, list(df_parmsPHH_cpt5.loc[:,"ParName"]), list(df_parmsPHH_cpt5.loc[:,list(chain.from_iterable(newColNames))[4]]),fromTp1,odeDict))
    # sols_24hr.append(simulateModel(concOfInterest, t2, list(df_parmsPHH_cpt5.loc[:,"ParName"]), list(df_parmsPHH_cpt5.loc[:,list(chain.from_iterable(newColNames))[5]]),fromTp1,odeDict))
    # sols_24hr.append(simulateModel(concOfInterest, t2, list(df_parmsPHH_cpt5.loc[:,"ParName"]), list(df_parmsPHH_cpt5.loc[:,list(chain.from_iterable(newColNames))[6]]),fromTp1,odeDict))
    # sols_24hr.append(simulateModel(concOfInterest, t2, list(df_parmsPHH_cpt5.loc[:,"ParName"]), list(df_parmsPHH_cpt5.loc[:,list(chain.from_iterable(newColNames))[7]]),fromTp1,odeDict))

    pool = mp.Pool(CORES)
    sols_24hr = [pool.apply(simulateModel, args = (concOfInterest, t2, list(df_parmsPHH_cpt5.loc[:,"ParName"]), list(df_parmsPHH_cpt5.loc[:,colName]),fromTp1,odeDict,colName,idc)) for idc,colName in enumerate(list(chain.from_iterable(newColNames)))]
    pool.close()

    # 24 hr rows
    newNames = ROWNAMES_24hr

    newCols_24hr = [newNames, np.array([np.NaN]*len(newNames)), np.array([np.NaN]*len(newNames))] + sols_24hr

    # Make a data frame of the columns
    df_cpt5_24hr = pd.DataFrame([])
    for i, cName in enumerate(list(df_parmsPHH_cpt5.columns)):
        df_cpt5_24hr[cName] = newCols_24hr[i]

    # Append new rows with the values at 8 hr
    df_parmsPHH_cpt5 = df_parmsPHH_cpt5.append(df_cpt5_24hr)
    df_parmsPHH_cpt5.iloc[list(range(0,Ns))+list(range(len(df_parmsPHH_cpt5.index)-Ns,len(df_parmsPHH_cpt5.index)))]

    dfSubset = df_parmsPHH_cpt5.iloc[list(range(0,Ns))+list(range(len(df_parmsPHH_cpt5.index)-(Ns*2),len(df_parmsPHH_cpt5.index))),[0]+list(range(3, df_parmsPHH_cpt5.shape[1]))]
    dfSubset.set_index('ParName', inplace=True)
    dfSubsetCPT = dfSubset.T

    expressionDFlist.append(dfSubsetCPT)

# Make one data frame from list of data frames:
parmsDF = parmsDFlist[0]
for df in parmsDFlist[1:len(parmsDFlist)]:
    parmsDF = pd.merge(parmsDF, df.drop(columns = ["est_value"]), on=['ParName', 'init_value'])
expressionDF = pd.concat(expressionDFlist)


if withRNA:
    fileName = CURDATE + "_" + outputName + "_" + str(numberOfRuns) + "times" + str(
        numberOfSamples) + "SampleSimulations_parameterPerturbations_withRNA_" + str(CONC) + "uM.csv"
else:
    fileName = CURDATE + "_" + outputName + "_" + str(numberOfRuns) + "times" + str(
        numberOfSamples) + "SampleSimulations_parameterPerturbations_" + str(CONC) + "uM.csv"

# Save the output
print("Will save df here: %s" % (outputFolder + fileName))
parmsDF.to_csv(outputFolder + fileName)

if withRNA:
    fileName = CURDATE + "_" + outputName + "_" + str(numberOfRuns) + "times" + str(numberOfSamples) + "SampleSimulations_withRNA_" + str(CONC) + "uM.csv"
else:
    fileName = CURDATE + "_" + outputName + "_" + str(numberOfRuns) + "times" + str(numberOfSamples) + "SampleSimulations_" + str(CONC) + "uM.csv"

# Save the output
print("Will save df here: %s" % (outputFolder + fileName))
expressionDF.to_csv(outputFolder + fileName)

# Done!
print("Files saved and script ended!")
