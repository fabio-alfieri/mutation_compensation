'''
@author: fabio alfieri
'''

import pandas as pd
import os

wd = os.getcwd() + "/tangoResults/"

os.system("ls " + wd + " | grep aggregation | grep HUMAN.txt >" + wd + "files.txt")
listFiles = open(wd + "files.txt", "r")

for line in listFiles:
    print(line)
    file = open(wd + line.split('\n')[0])
    data = pd.read_csv(file, header = 0, delimiter = "\t")
    wtAggregation = data.iloc[0,1]
    FC = data.iloc[:,1]/wtAggregation
    FC_df = pd.DataFrame(FC)
    fileContent = pd.concat([data,FC_df], axis = 1)
    fileContent.to_csv(wd + line.split('\n')[0].split('.')[0] + "_result.txt", index = True)
    

os.system("cd " + wd)
os.system("ls " + wd + " | grep '_result.txt' > " + os.getcwd( + "/result_files.txt")
listFiles = open(os.getcwd( + "/result_files.txt", "r")
results = pd.DataFrame()

for line in listFiles:
    print(line)
    file = pd.read_csv(wd + line.split('\n')[0])
    results = pd.concat([results, file], axis = 0)
    

results[['UniProtKB','Species','Mutation', ""]] = results.Sequence.str.split("_", expand = True)
results = results.drop([""], axis = 1)

results.to_csv(os.getcwd() + "/results.csv", index = True)

