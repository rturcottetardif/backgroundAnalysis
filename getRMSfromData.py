#!/usr/bin/env python

"""
========================================================================================= #
Author : Roxanne Turcotte

Description : Takes the rms of the background traces and the time and saves
                it. Will be use to see the galaxie.
========================================================================================= #
"""
import numpy as np
# import pandas as pd

from icecube.icetray import I3Units
from I3Tray import *
from icecube import icetray, dataclasses, radcube, dataio
from icecube.icetray.i3logging import log_info, log_warn, log_fatal
from datetime import datetime

bandLimits = [70 * I3Units.megahertz, 350 * I3Units.megahertz]

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('inputName', type=str, nargs="+", default=[], help='i3Files soft triggers')
parser.add_argument('outputName', type=str, default=[], help='output file')

args = parser.parse_args()

class GalacticBackground(icetray.I3Module):
    def __init__(self, ctx):
        icetray.I3Module.__init__(self, ctx)
        self.AddParameter('InputName', 'InputName', "InputName" )
        self.AddParameter('TimeOutput', 'TimeOutput', "ResultsArrayTime")
        self.AddParameter('RmsOutput', 'RmsOutput', "ResultsArrayRms")

    def Configure(self):
        self.inputName = self.GetParameter('InputName')
        self.timeOutput = self.GetParameter('TimeOutput')
        self.rmsOutput = self.GetParameter('RmsOutput')
        print("... I am starting")

    def DAQ(self, frame):
        try:
            time = frame["TaxiTime"]
            time_np = np.datetime64(time.date_time)
            time_new = np.datetime64(time_np).astype(datetime)

            rms = []
            antennaDataMap = frame[self.inputName]
            for iant, antkey in enumerate(antennaDataMap.keys()):
                channelMap = antennaDataMap[antkey]
                for ichan, chkey in enumerate(channelMap.keys()):
                    fft = channelMap[ichan].GetFFTData()
                    timeSeries = fft.GetTimeSeries()
                    rms.append(radcube.GetRMS(timeSeries))
            self.timeOutput.append(time_new)
            self.rmsOutput.append(rms)
            self.PushFrame(frame)
        except:
            print("this frame was skipped...")

    def Finish(self):
        print("I'm done doing stuff")

""" Run for one file """
timeOutput, rmsOutput = [], []
tray = I3Tray()
tray.Add('I3Reader', "FileReader",
         FilenameList=args.inputName
         )

tray.AddModule("BandpassFilter", "BoxFilter",
               # InputName="TAXIRadioWaveform",
               InputName="TAXIRadioWaveform",
               OutputName="FilteredMap",
               FilterType=radcube.eBox,
               FilterLimits=bandLimits,
               # ButterworthOrder=radcube.defaults.butterworthOrder,
               ApplyInDAQ=True
               )

tray.AddModule(GalacticBackground, "TheGalaxyObserver",
               InputName="FilteredMap",
               TimeOutput=timeOutput,
               RmsOutput=rmsOutput,
               )

tray.Execute()
timeOutput = np.asarray(timeOutput)
rmsOutput = np.asarray(rmsOutput)
datadir= "/data/user/rturcotte/analysis/background/data/"
print(timeOutput.shape)
print(rmsOutput.shape)
np.savez(datadir + args.outputName, time=timeOutput, rms10=rmsOutput[:, 0], rms11=rmsOutput[:, 1], rms20=rmsOutput[:, 2], rms21=rmsOutput[:, 3], rms30=rmsOutput[:, 4], rms31=rmsOutput[:, 5])

# data = np.load("/home/rturcotte/work/scripts/crazyIdeas/data/BGTest.npz", allow_pickle=True)
# print(data.files)
# print(data["rms"])
# print(data["rms"][0][1])
