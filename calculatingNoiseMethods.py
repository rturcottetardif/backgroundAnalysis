#!/bin/usr/env python3

#!/bin/usr/env python3
import numpy as np
import os
import matplotlib.pyplot as plt
from datetime import datetime

import sys
sys.path.insert(1, '/home/rturcotte/work/modules')
from TriggerPicker import TriggerPicker

from icecube import dataio, dataclasses, radcube
from icecube.icetray.i3logging import log_warn, log_info, log_fatal
from icecube.icetray import I3Frame, I3Units
#import utils.extractI3Variables as i3var
from I3Tray import *

SMALL_SIZE = 12
MEDIUM_SIZE = 14
BIGGER_SIZE = 16

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


class NoiseCalculation(icetray.I3Module):
    def __init__(self, ctx):
        icetray.I3Module.__init__(self, ctx)
        self.AddParameter('InputName', 'InputName', 'InputName')
        self.AddParameter('NoiseWindow', 'NoiseWindow', 32)
        self.AddParameter('ApplyInDAQ', 'ApplyInDAQ', True)
        self.AddParameter('OutFig', 'name of the output figure',
                          "calculatingNoiseMethods.png")
        self.AddParameter('OutFile', 'name of the .npz file', "calculatingNoiseMethods.npz")

    def Configure(self):
        self.inputName = self.GetParameter('InputName')
        self.noiseWindow = self.GetParameter('NoiseWindow')
        self.applyInDAQ = self.GetParameter('ApplyInDAQ')
        self.outfig = self.GetParameter('OutFig')
        self.outfile = self.GetParameter('OutFile')
        self.counts = 0
        self.colors = ["c", "b", "m", "r", "y", "g"]
        self.plotting = False

        self.noise_rms_step = []
        self.noise_rms_window = []
        self.time = []

        log_info("calculating the noise with a {0} window".format(self.noiseWindow))

    def GetNoise(self, frame):
        if frame.Has("TaxiTime"):
            time = frame["TaxiTime"]
        elif frame.Has("RadioTaxiTime"):
            time = frame["RadioTaxiTime"]
        else:
            print("no time found...")
        # time_np = np.datetime64(time.date_time)
        # time_new = np.datetime64(time_np).astype(datetime)
        self.time.append(np.datetime64(time.date_time).astype(datetime))
        self.counts += 1
        antennaDataMap = frame[self.inputName]
        for iant, antkey in enumerate(antennaDataMap.keys()):
            channelMap = antennaDataMap[antkey]
            for ichan, chkey in enumerate(channelMap.keys()):
                fft = channelMap[ichan].GetFFTData()
                timeSeries = fft.GetTimeSeries()

                # Only taking a window at the beginning
                noise_timeSeries = timeSeries.GetSubset(10, 410)
                # Noise RMS is the one that should be used for SNR
                self.noise_rms_window.append(radcube.GetRMS(noise_timeSeries))
                # noise_power = radcube.GetPower(noise_timeSeries, 50*I3Units.ohm)

                # Takes sections, calculate power (or RMS), remove the maximums, take the median
                # take the whole WF, pulse will be removed with the removal of maximums
                noise_rms_step, noise_power_step = [], []
                steps = np.arange(0, len(timeSeries), self.noiseWindow)
                numberOfSubTraces = len(timeSeries)//self.noiseWindow
                for i in range(numberOfSubTraces - 1):
                    noise_chopped = timeSeries.GetSubset(
                        int(steps[i]), int(steps[i + 1])
                        )
                    noise_rms_step.append(
                        radcube.GetRMS(
                            noise_chopped
                            )
                        )
                    # noise_power_step.append(
                    #     radcube.GetPower(
                    #         noise_chopped, 50*I3Units.ohm
                    #         )
                    #     )

                # We keep only 10 minimum value and average them.
                # The rational is that RFI screw up the noise, but we can hardly have lower noise
                noise_rms_step.sort()
                self.noise_rms_step.append(
                    np.median(noise_rms_step[:int(numberOfSubTraces/2)]))
                # noise_power_step.sort()
                # self.noise_power_step.append(np.median(noise_power_step[:10]))

    def DAQ(self, frame):
        if self.applyInDAQ:
            self.GetNoise(frame)
            self.PushFrame(frame)

    def Physics(self, frame):
        if not self.applyInDAQ:
            self.GetNoise(frame)
            self.PushFrame(frame)

    def Finish(self):
        noise_std = np.asarray(self.noise_rms_window).reshape(self.counts, 3, 2)
        noise_step = np.asarray(self.noise_rms_step).reshape(self.counts, 3, 2)
        time = np.asarray(self.time)
        print("counts", self.counts)
        if self.plotting:
            fig, ax = plt.subplots(
                figsize=[16, 5], nrows=1, ncols=1, tight_layout=True
                # gridspec_kw={'width_ratios': [2, 1]}
                )
            for iant in range(3):
                for ich in range(2):
                    ax.plot_date(
                        time, noise_std[:, iant, ich], marker=".", alpha=0.5,
                        color=self.colors[2*iant + ich]
                        )
                    ax.plot_date(
                        time, noise_step[:, iant, ich], marker="x", alpha=0.5,
                        color=self.colors[2*iant + ich]
                        )
            ax.set_ylabel("amplitude / ADC")
            fig.savefig(self.outfig)

        # Saving the file :
        np.savez(self.outfile, time=time, noiseStd=noise_std, noiseStep=noise_step)

if __name__ == '__main__':
    print("I've been called !")
    noiseWindow = 64
    path = "/data/user/rturcotte/analysis/background/noiseMethod/"
    filename = "/data/exp/IceCube/2022/unbiased/surface/radio/V5/eventData_1645525547_2022-02-22_10-25-47.i3.gz"
    filename = "/data/user/rturcotte/poleData/CleanedBackground_2020-12.i3.gz"

    tray = I3Tray()
    tray.Add(
        'I3Reader', "dataReader",
        FilenameList=[filename]
        )
    # tray.AddModule(
    #     TriggerPicker, "TheSoftTriggerPicker",
    #     Trigger="Soft",
    #     ApplyInDAQ=True
    #     )
    # tray.Add(
    #     radcube.modules.RemoveTAXIArtifacts, "ArtifactRemover",
    #     InputName="TAXIRadioWaveform",
    #     OutputName="CleanedWaveform",
    #     BaselineValue=0,
    #     medianOverCascades=False,
    #     RemoveBinSpikes=True,
    #     BinSpikeDeviance=int(2**12),
    #     RemoveNegativeBins=True
    #     )
    # Filtering just to put the baseline to zero
    tray.AddModule(
        "BandpassFilter", "BoxFilter",
        InputName="CleanedWaveform",
        OutputName="FilteredMap",
        FilterType=radcube.eBox,
        FilterLimits=[1 * I3Units.megahertz, 500 * I3Units.megahertz],
        ApplyInDAQ=True
        )

    tray.AddModule(
        NoiseCalculation, "TheNoiseCalculator",
        InputName="FilteredMap",
        noiseWindow=noiseWindow,
        OutFig=path + "noiseCalculation_{0}.png".format(noiseWindow),
        OutFile=path + "noiseCalculation_CleanedBackground_2020-12_{0}.npz".format(noiseWindow)
        )
    tray.Execute()
