#!/bin/env python3

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import numpy as np

# Custom modules
import sys
sys.path.insert(1, '/home/rturcotte/work/modules')
from MedianNonCascaded import MedianNonCascaded
from TriggerPicker import TriggerPicker

sys.path.insert(1, '/home/rturcotte/work/modules')

from icecube.dataclasses import I3AntennaGeo
from icecube.taxi_reader import taxi_tools
from icecube import radcube, dataclasses, icetray
from icecube.radcube import defaults
from icecube.icetray import I3Units
from I3Tray import *


SMALLER_SIZE = 12
SMALL_SIZE = 18
MEDIUM_SIZE = 20
BIGGER_SIZE = 22

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALLER_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

# You can request waveforms with various binning and length
WaveformLengths = [1024]
AverageOver = 1
Binning = [1 * I3Units.ns]
antType = I3AntennaGeo.AntennaType.SKALA2
electronicName = "electronicResponseName"
antennaName = "antennaResponseName"
year = "2022"

# Custom module to make waveforms with zero amplitude to boostrap the
# BringTheNoise module
class MakeEmptyWaveforms(icetray.I3Module):
    def __init__(self, ctx):
        icetray.I3Module.__init__(self, ctx)

    def Configure(self):
        self.NEntries = 0

    def Process(self):
        print("I made an empty waveform")
        # Making many frames with empty waveforms and only one antenna
        for k in range(AverageOver):
            frame = icetray.I3Frame(icetray.I3Frame.DAQ)
            antMap = dataclasses.I3AntennaDataMap()
            geometry = dataclasses.I3Geometry()

            if self.NEntries < len(WaveformLengths):

                antkey = radcube.GetGenericAntennaKey(self.NEntries)

                antGeo = dataclasses.I3AntennaGeo()
                antGeo.cableLength = 50 * I3Units.meter
                geometry.antennageo[antkey] = antGeo

                antChMap = dataclasses.I3AntennaChannelMap()

                fft = dataclasses.FFTData()
                timeSeries = fft.GetTimeSeries()
                timeSeries.binning = 1. * I3Units.ns
                timeSeries.offset = 0.

                for ibin in range(WaveformLengths[self.NEntries]):
                    timeSeries.PushBack(0.)

                antCh = dataclasses.I3AntennaChannel(fft)
                antChMap[0] = antCh
                antCh2 = dataclasses.I3AntennaChannel(fft)
                antChMap[1] = antCh2
                antMap[antkey] = antChMap

            frame["EmptyWaveform"] = antMap
            frame[radcube.GetDefaultGeometryName()] = geometry

            self.NEntries += 1
            if self.NEntries == len(WaveformLengths):
                self.RequestSuspension()
        self.PushFrame(frame)


class Expectation(icetray.I3Module):
    def __init__(self, ctx):
        icetray.I3Module.__init__(self, ctx)
        self.AddParameter('OutFig', 'OutFig', "spectrum.png")
        self.AddParameter('OutFile', 'OutFile', "spectrum.npy")
        self.AddParameter('Plotting', 'Plotting', False)
        self.AddParameter('Saving', 'Saving', False)
        self.NStartBins = 10
        self.NTimeBins = 1015
        self.NFreqBins = int((self.NTimeBins-self.NStartBins) / 2 + 1)
        self.NEntries = 0
        self.df = 1*I3Units.ns
        # self.NFreqBins_empty = int(WaveformLengths[0] / 2 + 1)
        self.plotdBm = np.zeros(self.NFreqBins)
        self.plotFreqs = np.zeros(self.NFreqBins)
        # self.AverageSpectrum = np.zeros((self.NAntennas, self.NArms, self.NFreqBins))

        self.fig, self.ax = plt.subplots(figsize=[8, 5], tight_layout=True,
                                         nrows=1, ncols=1)

    def Configure(self):
        self.outfig = self.GetParameter('OutFig')
        self.outfile = self.GetParameter('OutFile')
        self.plotting = self.GetParameter('Plotting')
        self.saving = self. GetParameter('Saving')

    def MakeExpectation(self, ax, name, frame):
        antDataMap = frame[name]
        antkey = antDataMap.keys()[0]

        chDataMap = antDataMap[antkey]
        chkey = chDataMap.keys()[0]

        fft = chDataMap[chkey].GetFFTData()
        spectrum = fft.GetFrequencySpectrum()
        freqs, amps = radcube.RadTraceToPythonList(spectrum)

        amps = [radcube.GetDbmHzFromFourierAmplitude(
            abs(thisAmp), spectrum.binning, 50 * I3Units.ohm) for thisAmp in amps]

        plotFreqs = []
        plotdBm = []

        for (freq, amp) in zip(freqs, amps):
            if amp > -200:
                plotFreqs.append(freq)
                plotdBm.append(amp)
                # print("{0}, {1:0.3f},".format(freq/I3Units.megahertz, amp))

        if self.plotting:
            print("I'm plotting the spectrum expectation .... ", self.outfig)
            ax.plot(np.array(plotFreqs) / I3Units.megahertz, plotdBm, c="gray", label="Galactic Noise")
            ax.set_xlabel("Frequency / MHz")
            ax.set_ylabel("Spectral power / dBm Hz$^-1$")
            ax.set_xlim(0, max(np.array(plotFreqs) / I3Units.megahertz) * 1.1)
        if self.saving:
            print("I'm saving the file ... ", self.outfile)
            spectrum_array = np.append(np.array(plotFreqs) / I3Units.megahertz, np.array(plotdBm))
            spectrum_array = spectrum_array.reshape(2, -1)
            np.save(self.outfile, spectrum_array)

    def DAQ(self, frame):
        self.MakeExpectation(self.ax, "FoldedWaveforms", frame)
        self.PushFrame(frame)

    def Finish(self):
        print("i'm finishing")
        self.fig.savefig(self.outfig, transparent=True)


saveplotpath = "/data/user/rturcotte/analysis/background/comparisonDAQs/"


tray = I3Tray()
# tray.Add('I3Reader', "SimReader",
#          Filename="/data/exp/IceCube/2022/unbiased/surface/radio/V5/eventData_1645525547_2022-02-22_10-25-47.i3.gz"
#          )
electronicName = "electronicResponse"
# import os
# file = os.getenv('I3_BUILD') + \
#     "/radcube/resources/data/test-data/electronic-response/Dummy1DResponse.txt"
if year == "2022":
    print("Taking TAXI 3.2 ...")
    tray.AddService("I3ElectronicsResponseFactory", electronicName,
                    AntennaType=antType,
                    IncludeLNA=True,
                    IncludeCables=True,
                    CableTemperature=radcube.defaults.cableTemp,
                    IncludeRadioBoard=False,
                    IncludeTaxi=False,
                    InstallServiceAs=electronicName,
                    CustomResponseFiles=["/cvmfs/icecube.opensciencegrid.org/users/acoleman/radcube-datasets/electronic-response/TAXIv3.2_Radioboardv2_2021.04.27.dat"],
                    )


else:
    tray.AddService("I3ElectronicsResponseFactory", electronicName,
                    IncludeLNA=True,
                    IncludeCables=True,
                    IncludeRadioBoard=True,
                    AdditionalGain=0,
                    IncludeTaxi=True,
                    CustomResponseFiles=[],
                    InstallServiceAs=electronicName
                    )


antennaName = radcube.CreateDefaultAntennaResponse(tray)
# antType = I3AntennaGeo.AntennaType.SKALA2

tray.AddService("I3GSLRandomServiceFactory", "gslRandom",
                Seed=666,
                InstallServiceAs="gslRandom")


tray.AddModule(MakeEmptyWaveforms, "Empty")

tray.AddModule("BringTheNoise", "NoiseGenerator",
               AntennaResponseName=antennaName,
               UseThermalNoise=True,
               ThermalNoiseTemp=40 * I3Units.kelvin,
               UseCaneNoise=True,
               RandomServiceName="gslRandom",
               InputName="EmptyWaveform",
               OutputName="NoisyWaveform",
               # ApplyInDAQ=False
               )


tray.AddModule("ElectronicResponseAdder", "AddElectronics",
               InputName="NoisyWaveform",
               OutputName="FoldedWaveforms",
               ElectronicsResponse=electronicName,
               # ApplyInDAQ=False
               )

tray.AddModule(Expectation, "CaneExpectation",
               Plotting=True,
               Saving=True,
               OutFig=saveplotpath + "spectrum_test.png",
               OutFile=saveplotpath + "spectrum_expectation_{0}.npy".format(year))


tray.Execute(1)

