#!/bin/env python3

# This script makes a plot of the expected amplitude according to the Cane model
# plus thermal noise. Outputs a plot called CaneExpectation.pdf

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import numpy as np

from icecube.dataclasses import I3AntennaGeo
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
Binning = [1 * I3Units.ns]
antType = I3AntennaGeo.AntennaType.SKALA2
electronicName = "electronicResponseName"
antennaName = "antennaResponseName"


# Custom module to make waveforms with zero amplitude to boostrap the
# BringTheNoise module
class MakeEmptyWaveforms(icetray.I3Module):

    def __init__(self, ctx):
        icetray.I3Module.__init__(self, ctx)

    def Configure(self):
        self.NEntries = 0

    def Process(self):
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
        self.PushFrame(frame)
        if self.NEntries == len(WaveformLengths):
            self.RequestSuspension()


class PlotExpectation(icetray.I3Module):

    def __init__(self, ctx):
        icetray.I3Module.__init__(self, ctx)
        self.fig = plt.figure(figsize=(8, 5))
        self.gs = gridspec.GridSpec(
            len(WaveformLengths), 1, wspace=0.1, hspace=0.1)
        self.gsCount = 0

    def MakePlot(self, ax, name, frame):
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

        ax.plot(np.array(plotFreqs) / I3Units.megahertz, plotdBm, c="gray", label="Galactic Noise")
        ax.set_xlabel("Frequency [MHz]")
        ax.set_ylabel("Spectral power [dBm/Hz]")
        ax.set_xlim(0, max(np.array(plotFreqs) / I3Units.megahertz) * 1.1)

    def DAQ(self, frame):

        ax = self.fig.add_subplot(self.gs[self.gsCount])
        self.gsCount += 1

        self.MakePlot(ax, "FoldedWaveforms", frame)
        print("Making file CaneExpectation.pdf")


tray = I3Tray()
electronicName = "Tester"
import os
file = os.getenv('I3_BUILD') + \
    "/radcube/resources/data/test-data/electronic-response/Dummy1DResponse.txt"
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
               OutputName="NoisyWaveform"
               )


tray.AddModule("ElectronicResponseAdder", "AddElectronics",
               InputName="NoisyWaveform",
               OutputName="FoldedWaveforms",
               ElectronicsResponse=electronicName
               )

tray.AddModule(PlotExpectation, "Plotter")

tray.Execute(1)

spectrumFilesDir = "/home/acoleman/analyses/taxi-studies/data/median-spectra/"
filename = "2020-04-17.npy"
spectrum = np.load(spectrumFilesDir + filename)
freqs = np.array(spectrum[0])
medianSpectrum = np.array(spectrum[1:, :])
lw = 0.7

#for i, spec in enumerate(medianSpectrum):
plt.plot(freqs[2:]/I3Units.megahertz, medianSpectrum[0][2:], lw=lw, c="m", label="antenna 1, polarisation 0")
plt.plot(freqs[2:]/I3Units.megahertz, medianSpectrum[1][2:], ls="--", lw=lw, c="m", label="antenna 1, polarisation 1")
plt.plot(freqs[2:]/I3Units.megahertz, medianSpectrum[2][2:], lw=lw, c="y", label="antenna 2, polarisation 0")
plt.plot(freqs[2:]/I3Units.megahertz, medianSpectrum[3][2:], ls="--", lw=lw, c="y", label="antenna 2, polarisation 1")
plt.plot(freqs[2:]/I3Units.megahertz, medianSpectrum[4][2:], lw=lw, c="c", label="antenna 3, polarisation 0")
plt.plot(freqs[2:]/I3Units.megahertz, medianSpectrum[5][2:], ls="--", lw=lw, c="c", label="antenna 3, polarisation 1")
plt.legend()
plt.savefig("SpectrumAndCane_withTempNoise_ForMarie.pdf", bbox_inches='tight', transparent=True)
plt.close()
