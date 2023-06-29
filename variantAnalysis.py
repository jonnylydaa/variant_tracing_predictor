import sys
import argparse
import pandas as pd
import csv
import matplotlib.pyplot as plt
from datetime import datetime
import numpy as np

"""Finds the average variant emergence time based on the information present in a disease's 
metadata file. Can be called with: python variantAnalysis.py -i (metadata file name)"""

def argParser():
    """Command line argument parser, takes in input metadata file."""
    parser = argparse.ArgumentParser(description="P.")
    parser.add_argument("-i", "--input", required=True, help='metadata tsv file')
    args = parser.parse_args()
    return args

def validate(date_text):
    """
    Checks if datetime is valid. Takes in a date as input. 
    Code from: https://stackoverflow.com/questions/16870663/how-do-i-validate-a-date-string-format-in-python
    """
    try:
        if date_text != datetime.strptime(date_text, "%Y-%m-%d").strftime('%Y-%m-%d'):
            raise ValueError
        return True
    except ValueError:
        return False
    
def casesPerStrain(spd, oldSD):
    """Finds how many postive cases per each strain. Takes spd (samples per day) dict and the oldSD (old sample dictionary)."""
    cpd = []
    strainCount = 1
    for day, cases in spd:
        if day in oldSD:
            strainCount += 1
        casesPerStrain = cases/strainCount
        cpd.append((day, casesPerStrain))
    return cpd
        
        

def parseMetadata(inFile):
    """Parses the inputted metadata file and returns dictionaries that have information
    on each sample, recorded date of the sample, and the variant the sample belongs to. Also
    includes how many positive cases per day every day."""
    samples = {}
    strainDict = {}
    samplesPerDay = {}   # virus wide number
    with open(inFile) as f:
        data = csv.reader(f, delimiter="\t")
        count = 0
        oldestDate = None
        for line in data:
            if count != 0:
                # print(line)
                if validate(line[2]):
                    currDate = datetime.strptime(line[2], "%Y-%m-%d")
                    if oldestDate == None:
                        oldestDate = currDate
                    elif currDate < oldestDate:
                        oldestDate = currDate
                    
                    if currDate not in samplesPerDay:
                        samplesPerDay[currDate] = 1
                    else:
                        samplesPerDay[currDate] += 1
                    
                    samples[line[0]] = (line[2], line[9])
                    if line[9] not in strainDict:
                        strainDict[line[9]] = []
                    strainDict[line[9]].append(line[0])
                
            count += 1
    return samples, strainDict, oldestDate, samplesPerDay


def findEmergenceRate(strainDict, sampleDict, oldestSample):
    """Finds the variant emergence rate of a disease. a dictionary with strain information (strainDict), 
    a dictionary with sample information (sampleDict), and the oldest sample (oldestSample) (Wuhan/Hu-1) are inputted.
    Returns the average variant emergence rate, a strain dictionary that gives the oldest sample for each strain, and
    a list of tuples that include a strain name and the oldest recorded date of a positve case of that strain (strainList)."""
    strainEmergence = {}
    oldestStrainDict = {}
    for strain, nodes in strainDict.items():
        oldestDate = None
        for node in nodes:
            if oldestDate is None:
                oldestDate = datetime.strptime(sampleDict[node][0], "%Y-%m-%d")
            else: 
                newdate = datetime.strptime(sampleDict[node][0], "%Y-%m-%d")
                # if sampleDict[node][0]   # date
                if newdate < oldestDate:
                    oldestDate = newdate
        # print("oldest date, strain: ", oldestDate, strain, len(nodes))
        oldestStrainDict[oldestDate] = strain
        dateDifference = oldestDate - oldestSample
        strainEmergence[strain] = dateDifference.days
    # print("strain emergence times: ", strainEmergence)
    newStrainEmergence = []
    newStrainDict = {}
    strainList = []
    for strain, days in strainEmergence.items():
        newStrainEmergence.append((strain, days))
    newStrainEmergence.sort(key=lambda a: a[1])
    for x in range(0, len(newStrainEmergence)):
        totalDays = 0
        count = 1
        currentStrainDays = newStrainEmergence[x][1]
        totalDays += currentStrainDays
        for y in range(x, 0, -1):
            if y != x:
                totalDays += currentStrainDays - newStrainEmergence[y][1]
                count += 1
        avgEmergence = totalDays/count
        newStrainDict[newStrainEmergence[x][0]] = avgEmergence
        strainList.append((newStrainEmergence[x][0], avgEmergence))
    # print(strainList)
    strainAmount = 0
    totalDays = 0
    for k, v in newStrainDict.items():
        totalDays += v
        strainAmount += 1
    emergenceRate = totalDays/strainAmount
    return emergenceRate, oldestStrainDict, strainList

def printPlots(spd):
    """Prints the number of cases per day. Inputted samples per day."""
    days = []
    infected = []
    for k, v in spd.items():
        days.append(k)
        infected.append(v)
    dataframe = pd.DataFrame({'date': np.array(days), 'cases': infected})
    plt.plot_date(dataframe.date, dataframe.cases, lw=.02, markersize=1.8)
    plt.title('Number of Covid-19 Cases Per Day')
    plt.xticks(rotation=30, ha='right')
    plt.xlabel('Date')
    plt.ylabel('Cases')
    plt.show()
    return 

def printEmergenceRatePlot(erList, avgeRate):
    """Prints the variant emergence rate plot. Includes the average emergence rate to compare against the average emergence
    rate for each variant."""
    strain, rates = zip(*erList)
    plt.plot(strain, rates, marker='o', markersize=3.2, label="average days until variant emerged")
    plt.title('Number of Days for Variants to Emerge')
    plt.axhline(y=avgeRate, color='r', linestyle='--', label="predicted average variant emergence rate")
    plt.xticks(rotation=70, ha='right')
    plt.xlabel('Variant')
    plt.ylabel('Days')
    plt.legend()
    plt.show()
    return 

def main():
    """Finds the variant emergence rate based on a disease's metadata file."""
    # parser = argParser()
    # args = parser.parse_args()
    args = argParser()
    sampleDict, strainDict, oldestSample, samplesPD = parseMetadata(args.input)
    avgEmergenceRate, oldStrainDict, erList = findEmergenceRate(strainDict, sampleDict, oldestSample)
    print("avg emergence rate: ", avgEmergenceRate)
    printEmergenceRatePlot(erList, avgEmergenceRate)
    
if __name__ == "__main__":
    """Calls main."""
    main()