#!/usr/bin/env python

"""
========================================================================================= #
Author : Roxanne Turcotte

Description : Takes the rms of the background traces and the time and saves
                it. Will be use to see the galaxie.
========================================================================================= #
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
# import pandas as pd

from icecube.icetray import I3Units
from I3Tray import *
from icecube import icetray, dataclasses, radcube, dataio
from icecube.icetray.i3logging import log_info, log_warn, log_fatal
from datetime import datetime

bandLimits = [70 * I3Units.megahertz, 350 * I3Units.megahertz]
ant = 1
pol = 0
plotdir = "/data/user/rturcotte/analysis/background/plot/"

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('inputName', type=str, nargs="+", default=[], help='i3Files scint triggers')
args = parser.parse_args()


SMALL_SIZE = 14
MEDIUM_SIZE = 16
BIGGER_SIZE = 18

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


class PlotRMS(icetray.I3Module):
    def __init__(self, ctx):
        icetray.I3Module.__init__(self, ctx)
        self.AddParameter('InputName', 'InputName', "InputName" )

    def Configure(self):
        self.inputName = self.GetParameter('InputName')
        self.timeOutput = []
        self.baselineOutput = []

        self.fig = plt.figure(figsize=[16, 9])
        self.spec = GridSpec(3, 2)

    def DAQ(self, frame):
        #left, bottom, width, height = [0.1, 0.5, 0.1, 0.1]
        print("Imma plot the RMS of that file ...")
        time = frame["TaxiTime"]
        print(time)
        time_np = np.datetime64(time.date_time)
        time_new = np.datetime64(time_np).astype(datetime)

        baseline = []
        antennaDataMap = frame[self.inputName]
        for iant, antkey in enumerate(antennaDataMap.keys()):
            channelMap = antennaDataMap[antkey]
            for ichan, chkey in enumerate(channelMap.keys()):
                fft = channelMap[pol].GetFFTData()
                timeSeries = fft.GetTimeSeries()
                time, timeSeries_np = radcube.RadTraceToPythonList(timeSeries)
                #rms = radcube.GetPower(timeSeries, 1)
                #rms = radcube.GetIntegralBackgroundSubtracted(timeSeries, 1)
                rms = radcube.GetRMS(timeSeries)
                baseline.append(radcube.GetRMS(timeSeries))
                ax = self.fig.add_subplot(self.spec[iant, ichan])
                ax.plot_date(time_new, rms, marker=".", ls="", color="k", alpha=0.5)
                ax.set_ylabel("ADC counts")
                ax.set_xlabel("Date")

        self.timeOutput.append(time_new)
        self.baselineOutput.append(baseline)
        self.PushFrame(frame)

    def Finish(self):
        print("I'm done doing stuff")
        self.fig.suptitle("{0}".format(args.inputName[0].split(".")[0]))
        self.fig.autofmt_xdate()
        plt.tight_layout()
        self.fig.savefig(plotdir + "RMSvsDate_{0}.pdf".format(args.inputName[0].split(".")[0]))


class PlotWaveforms(icetray.I3Module):
    def __init__(self, ctx):
        icetray.I3Module.__init__(self, ctx)
        self.AddParameter('InputName', 'InputName', "InputName" )

    def Configure(self):
        self.inputName = self.GetParameter('InputName')
        self.timeOutput = []
        self.baselineOutput = []
        self.rois_low = []
        self.rois_high = []
        self.spectrums_low = []
        self.spectrums_high = []

        self.fig = plt.figure(figsize=[16, 9])
        self.fig2 = plt.figure(figsize=[16, 9])
        widths = [3, 1, 0.5]
        heights = [3, 3]
        self.spec = GridSpec(2, 3, width_ratios=widths, height_ratios=heights)
        self.spec2 = GridSpec(1, 1)

    def DAQ(self, frame):
        #left, bottom, width, height = [0.1, 0.5, 0.1, 0.1]
        print("Imma plot some ugly waveforms...")
        time = frame["TaxiTime"]
        print(time)
        time_np = np.datetime64(time.date_time)
        time_new = np.datetime64(time_np).astype(datetime)

        baseline = []
        antennaDataMap = frame[self.inputName]
        antkey = antennaDataMap.keys()[ant]
        channelMap = antennaDataMap[antkey]
        #for ichan, chkey in enumerate(channelMap.keys()):
        fft = channelMap[pol].GetFFTData()
        timeSeries = fft.GetTimeSeries()
        time, timeSeries_np = radcube.RadTraceToPythonList(timeSeries)
        rms = radcube.GetRMS(timeSeries)
        baseline.append(radcube.GetRMS(timeSeries))
        # Low RMS
        if rms < 100:
            self.rois_low.append(frame["AntennaROI"][ant])
            self.spectrums_low.append(fft.GetSpectrum())
            ax = self.fig.add_subplot(self.spec[0, 0])
            ax.plot(time, timeSeries_np, ls="-", color="k",
                    lw=0.5, alpha=0.5)
            ax.set_xlim(500, 1500)
            ax.set_ylim(-500, 500)
            ax.set_ylabel("ADC counts")
            ax.set_xlabel("Time [ns]")

            ax2 = self.fig.add_axes([0.16, 0.57, 0.2, 0.1])
            ax2.plot(time[500:600], timeSeries_np[500:600], ls="-", color="k",
                     lw=0.5, alpha=0.5)
            ax2.set_ylim(-300, 300)

            ax = self.fig.add_subplot(self.spec[0, 1])
            ax.plot(time, timeSeries_np, ls="-", color="k",
                    lw=0.5, alpha=0.5)
            ax.axvline(1024, lw=0.5, color="red", ls=":")
            ax.set_xlim(980, 1060)
            ax.set_ylim(-500, 500)
            ax.set_ylabel("")
            ax.set_xlabel("Time [ns]")

        # High RMS
        elif rms > 100:
            self.rois_high.append(frame["AntennaROI"][ant])
            ax = self.fig.add_subplot(self.spec[1, 0])
            ax.plot(time, timeSeries_np, ls="-", color="purple",
                    lw=0.2, alpha=0.5)
            ax.set_xlim(500, 1500)
            ax.set_ylabel("ADC counts")
            ax.set_xlabel("Time [ns]")

            ax2 = self.fig.add_axes([0.16, 0.15, 0.2, 0.1])
            ax2.plot(time[500:600], timeSeries_np[500:600], ls="-", color="purple",
                     lw=0.5, alpha=0.5)
            ax2.set_ylim(-300, 300)

            ax = self.fig.add_subplot(self.spec[1, 1])
            ax.plot(time, timeSeries_np, ls="-", color="purple",
                    lw=0.2, alpha=0.5)
            ax.axvline(1024, lw=0.5, color="red", ls=":")
            ax.set_xlim(980, 1060)
            ax.set_ylabel("")
            ax.set_xlabel("Time [ns]")
        self.timeOutput.append(time_new)
        self.baselineOutput.append(baseline)
        self.PushFrame(frame)

    def Finish(self):
        print("I'm done doing stuff")
        bins=np.arange(-0.5, 1024.5, 1)
        ax = self.fig.add_subplot(self.spec[0, 2])
        ax.hist(self.rois_low, histtype="step", bins=bins, color="k")
        ax.set_xlabel("ROI")

        ax = self.fig.add_subplot(self.spec[1, 2])
        ax.hist(self.rois_high, histtype="step", bins=bins, color="purple")
        ax.set_xlabel("ROI")

        self.fig.suptitle("Antenna {0} Polarisation {1}".format(ant+1, pol))
        self.fig.autofmt_xdate()
        #plt.tight_layout()
        self.fig.savefig(plotdir + "WfHighRMS_NF_{0}.pdf".format(args.inputName[0].split(".")[0]))


        #'GetDbmHzFromFourierAmplitude'



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
               ApplyInDAQ=True
               )

# tray.AddModule(PlotRMS, "TheRMSPlotter",
#                #InputName="FilteredMap",
#                InputName="TAXIRadioWaveform"
#                )

tray.AddModule(PlotWaveforms, "TheWaveformsPlotter",
               InputName="FilteredMap",
               #InputName="TAXIRadioWaveform"
               )

tray.Execute()


