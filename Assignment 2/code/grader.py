import numpy as np
import pdb
import argparse
import subprocess # For executing c++ executable
import pandas as pd
from timeit import default_timer as timer

###############################################################
################### Util Functions Below ######################

def convertPIs(aString):
    """ Input: A comma seperated string like "pi/2,pi/4,pi/2,pi/4,pi/2,"
            or 1.2,4,5.3 etc
    Output: string replacing the pis 1.57079,...,...
    """
    if aString[-1] == ",":  # Remove training comma if there is one
        aString = aString[:-1]
    aString = aString.replace("pi", "3.141592") # Replace pi with 3.14... if needed
    vecOfStrings = aString.split(",")
    ans = []
    for anExpression in vecOfStrings:
        ans.append(str(eval(anExpression))) # Evaluate expressions if needed
    return ans

###############################################################
################### Main Functions Below ######################


def graderMain(executablePath, gradingCSV):
    problems = [#["./map1.txt", "pi/2,pi/4,pi/2,0.785398,1.570796","0.392699,2.356194,pi,2.8274328,pi*3/2"],
            #["./map2.txt", "0.392699,2.356194,3.141592","1.570796,0.785398,1.570796"],
            ["./map2.txt", "1.46787,1.08484,5.01259,1.87066,4.84654","1.21053,1.85669,4.59725,0.84908,3.76929"], 
            ["./map2.txt", "1.32385,5.83845,1.85629,2.00166,5.38377","1.26862,1.11962,4.27879,1.71452,3.54118"],
            ["./map2.txt", "1.36943,2.39311,0.55128,1.74171,3.07647","1.67695,0.58518,3.91363,0.15771,4.33940"],
            ["./map2.txt", "0.54123,1.60231,4.30597,0.25755,1.27960","0.26064,1.53907,3.46564,1.72787,2.61059"],
            ["./map2.txt", "1.61502,1.93441,3.11048,1.43682,4.74141","1.17505,0.96505,3.62719,5.86555,0.00437"],                   ]
    scores = []
    for z in range(4):
        for aPlanner in [2]:
            for i, data in enumerate(problems):
                inputMap, startPos, goalPos = [*data]
                numDOFs = len(startPos.split(","))
                outputSolutionFile = "tmp.txt"
                startPosString = ",".join(convertPIs(startPos))
                goalPosString = ",".join(convertPIs(goalPos))
                commandPlan = "{} {} {} {} {} {} {}".format(
                    executablePath,
                    inputMap, numDOFs, startPosString, goalPosString,
                    aPlanner, outputSolutionFile)
                commandVerify = "./verifier.exe {} {} {} {} {}".format(
                    inputMap, numDOFs, startPosString, goalPosString,
                    outputSolutionFile)
                try:
                    start = timer()
                    subprocess.run(commandPlan.split(" "), check=True) # True if want to see failure errors
                    timespent = timer() - start
                    returncode = subprocess.run(commandVerify.split(" "), check=False).returncode
                    if returncode != 0:
                        print("Returned an invalid solution")
                    
                    ### Calculate the cost from their solution
                    with open(outputSolutionFile) as f:
                        line = f.readline().rstrip()  # filepath of the map
                        solution = []
                        for line in f:
                            solution.append(line.split(",")[:-1]) # :-1 to drop trailing comma
                        solution = np.asarray(solution).astype(float)
                        numSteps = solution.shape[0]

                        ## Cost is sum of all joint angle movements
                        difsPos = np.abs(solution[1:,]-solution[:-1,])
                        cost = np.minimum(difsPos, np.abs(2*np.pi - difsPos)).sum()

                        success = returncode == 0
                        scores.append([aPlanner, inputMap, i, numSteps, cost, timespent, success])
                
                    ### Visualize their results
                    # commandViz = "python visualizer.py tmp.txt --gifFilepath=grades_{}.gif".format(i)
                    # commandViz += " --incPrev=1"
                    # subprocess.run(commandViz.split(" "), check=True) # True if want to see failure errors
                except Exception as exc:
                    print("Failed: {} !!".format(exc))
                    scores.append([aPlanner, inputMap, i, -1, -1, timespent, False])

    ### Save all the scores into a csv to compute overall grades
    df = pd.DataFrame(scores, columns=["planner", "mapName", "problemIndex", "numSteps", "cost", "timespent", "success"])
    df.to_csv(gradingCSV, index=False)
            

if __name__ == "__main__":
    graderMain("./planner.exe", "test.csv")