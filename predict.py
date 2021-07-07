import glob, os, argparse, codecs, sys
import pandas as pd
import matplotlib.pyplot as plt
import shutil, os
from CSV_read_write import *
from RNN import FitRNetwork, NetworkRPredict

directory = sys.argv[1]

#MIH: if not new or after trained, it now predicts using the network
NetworkRPredict(directory) #testing
