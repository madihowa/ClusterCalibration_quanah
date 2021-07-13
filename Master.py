import glob, os, argparse, codecs, sys
import pandas as pd
import matplotlib.pyplot as plt
import shutil, os
from CSV_read_write import *
from RNN import FitRNetwork, NetworkRPredict

#Set the right input files
#MIH: reads the argument from the command line and saves it as variables in the code
parser = argparse.ArgumentParser(description='Creates a NN and trains it to calibrate topo-cluster energies')
parser.add_argument("-r", type=int, help="Pick a run number to call previous data else will make new network.", default=-1)
parser.add_argument("-d", type=str, help="Pick a date which the run is from, default today", default="today")
args = parser.parse_args()
run = args.r
day = args.d


csv_dir = "/lustre/work/madihowa/CERN/ClusterCalibration/pi_plus_minus_zero_data"
#MIH: determines whether the run is new or not
if run == -1:
    new_run = True
else:
    new_run = False

directory = get_directory(day, run)
file = open("NetworkHistory.txt", "a")
if run == -1:
    run = 0
    for root, dirs, files in os.walk(os.getcwd()):
        for dir in dirs:
            if day in dir:
                run = run + 1
file.write('\n{:}\t|{:1}\t|'.format(day, run))

print("\nRun is {} and Date is {}\n\n".format(run,day))

#MIH: copies RNN.py into the newly made Run directory
shutil.copy("RNN.py", directory)

#MIH: if the run is a new run, it trains the network
if new_run:
    FitRNetwork(directory, csv_dir)  #training
else:
    file.write('\t\t')

#MIH: if not new or after trained, it now predicts using the network
NetworkRPredict(directory, csv_dir) #testing
