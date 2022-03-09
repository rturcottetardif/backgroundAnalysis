#!/usr/bin/env python3

"""
========================================================================================= #
Author : Roxanne Turcotte

Description : Takes the rms of the background traces and the time and saves
                it. Will be use to see the galaxie.
========================================================================================= #
"""
# /data/exp/IceCube/2022/unbiased/surface/radio/V5/eventData_1643191520_2022-01-26_10-05-20.i3.gz test2022.npz
# /data/exp/IceCube/2021/unbiased/surface/radio/V5/*2021-01-26* test21.npz

import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import sys
sys.path.insert(1, '/home/rturcotte/work/modules')
from MedianNonCascaded import MedianNonCascaded

from icecube.icetray import I3Units
# from icecube.taxi_reader import taxi_tools
from icecube import icetray, dataclasses, radcube, dataio
from icecube.icetray.i3logging import log_info, log_warn, log_fatal
from icecube.dataclasses import I3AntennaGeo
from I3Tray import *

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('inputName', type=str, nargs="+", default=[], help='i3Files soft triggers')
# parser.add_argument('outputName', type=str, default=[], help='output file')
args = parser.parse_args()

headDirOverride = ""
antType = I3AntennaGeo.AntennaType.SKALA2
electronic30 = "electronicResponse30"
electronic32 = "electronicResponse32"
antennaName = "antennaResponseName"
GCDFile = "/data/user/rturcotte/gcd-files/GCD-AntennaSurvey_2020.06.15.i3.gz"
datadir = "/data/user/rturcotte/analysis/background/comparisonDAQs/"


name = "2022"             # for the plots
color = "red"
softTrig = False


def cutTraces(radTrace, lengthSubTraces=16, mode="rms"):
    steps = np.arange(0, len(radTrace), lengthSubTraces)
    nbSubTraces = len(radTrace) / lengthSubTraces
    temp = []
    for i in range(int(nbSubTraces)-1):
        chopped = radTrace.GetSubset(int(steps[i]), int(steps[i + 1]))
        if mode == "rms":
            temp.append(radcube.GetRMS(chopped))
        elif mode == "power":
            temp.append(radcube.GetPower(chopped, 50*I3Units.ohm))
        temp.sort()
    return temp


class GalacticBackground(icetray.I3Module):
    def __init__(self, ctx):
        icetray.I3Module.__init__(self, ctx)
        self.AddParameter('InputName', 'InputName', "InputName")
        self.AddParameter('Output', 'Output', "Output")
        self.AddParameter('SoftTrigger', 'SoftTrigger', 'SoftTrigger')
        self.AddParameter("ApplyInDAQ", "ApplyInDAQ", False)

    def Configure(self):
        self.inputName = self.GetParameter('InputName')
        self.output = self.GetParameter('Output')
        self.softTrigger = self.GetParameter('SoftTrigger')
        self.applyinDAQ = self.GetParameter("ApplyInDAQ")
        self.timeOutput = []
        self.baselineRms = []
        self.baselinePower = []
        print("... I am starting")

    def RunForOneFrame(self, frame):
        try:
        # if 1:
            if self.softTrigger:
                if frame.Has("SurfaceFilters"):
                    soft_trig = frame["SurfaceFilters"].get("soft_flag").condition_passed
                    if soft_trig:
                        print("Got a soft trig...")
                else:
                    soft_trig = True
                    print("you are not using V5 processing....")
                    print("I just assume you took the right file and I'll take every trigger")
            else:
                soft_trig = True
                print("I will be taking all trigger types...")
            if soft_trig:
                if frame.Has("TaxiTime"):
                    time = frame["TaxiTime"]
                elif frame.Has("RadioTaxiTime"):
                    time = frame["RadioTaxiTime"]
                else:
                    print("no time found...")
                time_np = np.datetime64(time.date_time)
                time_new = np.datetime64(time_np).astype(datetime)

                rmsTraces, powerTraces = [], []
                antennaDataMap = frame[self.inputName]
                for iant, antkey in enumerate(antennaDataMap.keys()):
                    channelMap = antennaDataMap[antkey]
                    for ichan, chkey in enumerate(channelMap.keys()):
                        fft = channelMap[ichan].GetFFTData()
                        timeSeries = fft.GetTimeSeries()
                        # timeSeries_chopped = timeSeries.GetSubset(10, 1001)
                        noises = cutTraces(timeSeries, lengthSubTraces=64, mode="rms")
                        rmsTraces.append(np.mean(noises[:10]))
                        noises = cutTraces(timeSeries, lengthSubTraces=64, mode="power")
                        powerTraces.append(np.mean(noises[:10]))


                self.timeOutput.append(time_new)
                self.baselineRms.append(rmsTraces)
                self.baselinePower.append(powerTraces)
        except:
            log_warn("this frame was skipped...")

    def DAQ(self, frame):
        if self.applyinDAQ:
            self.RunForOneFrame(frame)
            self.PushFrame(frame)

    def Physics(self, frame):
        if not self.applyinDAQ:
            self.RunForOneFrame(frame)
            self.PushFrame(frame)

    def Finish(self):
        print("Saving the data ...", self.output)
        timeOutput = np.asarray(self.timeOutput)
        baselineRms = np.asarray(self.baselineRms)
        baselinePower = np.asarray(self.baselinePower)
        print(timeOutput.shape)
        print(baselineRms.shape)
        print(baselinePower.shape)

        np.savez(datadir + self.output,
                 time=timeOutput,
                 rms10=baselineRms[:, 0],
                 rms11=baselineRms[:, 1],
                 rms20=baselineRms[:, 2],
                 rms21=baselineRms[:, 3],
                 rms30=baselineRms[:, 4],
                 rms31=baselineRms[:, 5],
                 power10=baselinePower[:, 0],
                 power11=baselinePower[:, 1],
                 power20=baselinePower[:, 2],
                 power21=baselinePower[:, 3],
                 power30=baselinePower[:, 4],
                 power31=baselinePower[:, 5],
                 )


class WaveformPlotter(icetray.I3Module):
    def __init__(self, ctx):
        icetray.I3Module.__init__(self, ctx)
        self.AddParameter('InputName', 'InputName', "InputName")
        self.AddParameter('OutFig', 'OutFig', "OutFig")

    def Configure(self):
        self.inputname = self.GetParameter('InputName')
        self.outfig = self.GetParameter('OutFig')
        self.fig, self.axs = plt.subplots(
            figsize=[16, 9], tight_layout=True,
            nrows=3, ncols=2
            )

    def Physics(self, frame):
        antennaDataMap = frame[self.inputname]
        for iant, antkey in enumerate(antennaDataMap.keys()):
            channelMap = antennaDataMap[antkey]
            for ichan, chkey in enumerate(channelMap.keys()):
                fft = channelMap[ichan].GetFFTData()
                timeSeries = fft.GetTimeSeries()
                hilbertEnvelope = dataclasses.fft.GetHilbertEnvelope(
                    timeSeries)
                hilbert_np = np.array(
                        radcube.RadTraceToPythonList(hilbertEnvelope))
                timeSeries_np = np.array(
                        radcube.RadTraceToPythonList(timeSeries))

                self.axs[iant, ichan].plot(hilbert_np[1], lw=0.5, c="gray")
                self.axs[iant, ichan].plot(timeSeries_np[1], lw=0.5)
                self.axs[iant, ichan].set_xlabel("Time /ns")
                self.axs[iant, ichan].set_ylabel("Deconvolved waveform")
                self.axs[iant, ichan].set_xlim(10, 1000)

        self.PushFrame(frame)

    def Finish(self):
        print('You did it Chief !!')
        outfig = datadir + self.outfig
        self.fig.savefig(outfig)


class SpectrumPlotter(icetray.I3Module):
    def __init__(self, ctx):
        icetray.I3Module.__init__(self,ctx)
        self.AddParameter('InputName', 'InputName', "InputName")
        self.AddParameter('OutFig', 'OutFig', "OutFig")

    def Configure(self):
        self.inputname = self.GetParameter('InputName')
        self.outfig = self.GetParameter('OutFig')

        self.NAntennas = 3
        self.NArms = 2
        self.NStartBins = 10
        self.NTimeBins = 1015
        self.NFreqBins = int((self.NTimeBins-self.NStartBins) / 2 + 1)

        self.AverageSpectrum = np.zeros((self.NAntennas, self.NArms, self.NFreqBins))
        # self.AverageTimeSeries = np.zeros((self.NAntennas, self.NArms, self.NTimeBins))

        self.NEntries = 0
        self.df = -1

    def Physics(self, frame):
        antennaDataMap = frame[self.inputname]
        for iant, antkey in enumerate(antennaDataMap.keys()):
            channelMap = antennaDataMap[antkey]
            for ichannel, chkey in enumerate(channelMap.keys()):
                timeSeries_chopped = []
                fft = channelMap[ichannel].GetFFTData()
                timeSeries = fft.GetTimeSeries()
                truncatedTimeSeries = timeSeries.GetSubset(self.NStartBins, self.NTimeBins)
                fft.LoadTimeSeries(truncatedTimeSeries)
                spectrum = fft.GetFrequencySpectrum()
                self.df = spectrum.binning
                truncatedFreqSpec = spectrum.GetSubset(0, self.NFreqBins - 1)
                freqs, pythonSpectrum = radcube.RadTraceToPythonList(truncatedFreqSpec)
                pythonSpectrum = np.abs(pythonSpectrum)
                self.AverageSpectrum[iant][ichannel] += pythonSpectrum
                self.NEntries += 1.
        self.PushFrame(frame)

    def Finish(self):
        log_info("We have grabbed information from all the events!")
        self.AverageSpectrum /= self.NEntries
        # self.AverageSpectrum = np.median(self.AverageSpectrum, axis=0)
        fig, axs = plt.subplots(
            figsize=[16, 9], tight_layout=True,
            nrows=3, ncols=2
            )

        for iant in range(len(self.AverageSpectrum)):
            for ichannel in range(len(self.AverageSpectrum[iant])):
                print("For antenna {0} and channel {1}".format(iant + 1, ichannel))
                dBmHz = []
                for i in range(len(self.AverageSpectrum[iant][ichannel])):
                    fourierAmp = self.AverageSpectrum[iant][ichannel][i]
                    dBmHz.append(radcube.GetDbmHzFromFourierAmplitude(fourierAmp, self.df, 50 * I3Units.ohm))
                axs[iant, ichannel].plot(dBmHz, lw=1.5, c=color)
                axs[iant, ichannel].set_xlabel("Frequency / MHz")
                axs[iant, ichannel].set_ylabel("Amplitude / dBm/Hz")
                axs[iant, ichannel].set_ylim(-200, -140)
                # np.savez(datadir + "spectrum_{0}-{1}_{2}.npz".format(iant, ichannel, name),
                #          spectrum=dBmHz)
        outfig = datadir + self.outfig#"spectrum_DeconvolvedCut_{0}.png".format(name)
        fig.savefig(outfig, transparent=True)

""" Run for one file """
tray = I3Tray()

if "2022" in name:
    electronicname = electronic32
    print("Taking TAXI 3.2 ...")
    tray.AddService("I3ElectronicsResponseFactory", electronic32,
                    AntennaType=antType,
                    IncludeLNA=True,
                    IncludeCables=True,
                    CableTemperature=radcube.defaults.cableTemp,
                    IncludeRadioBoard=False,
                    IncludeTaxi=False,
                    InstallServiceAs=electronic32,
                    # CustomResponseFiles=["/cvmfs/icecube.opensciencegrid.org/users/acoleman/radcube-datasets/electronic-response/TAXI_Plus_Radioboard_200603.dat"]
                    CustomResponseFiles=["/cvmfs/icecube.opensciencegrid.org/users/acoleman/radcube-datasets/electronic-response/TAXIv3.2_Radioboardv2_2021.04.27.dat"],
                   )

else:
    print("Taking TAXI 3.0 ...")
    electronicname = electronic30
    tray.AddService("I3ElectronicsResponseFactory", electronic30,
                    AntennaType=antType,
                    IncludeLNA=True,
                    IncludeCables=True,
                    CableTemperature=radcube.defaults.cableTemp,
                    IncludeRadioBoard=True,
                    IncludeTaxi=True,
                    InstallServiceAs=electronic30,
                    OverrideHeadDir=headDirOverride
                    )


tray.Add('I3Reader', "FileReader",
         FilenameList=[GCDFile, *args.inputName]
         )


tray.AddModule("I3NullSplitter", "splitter",
               SubEventStreamName="RadioEvent"
               )


tray.AddModule('PedestalRemover', "PedestalRemover",
               InputName="TAXIRadioWaveform",
               OutputName="PedestalRemovedMap",
               ConvertToVoltage=True,
               ElectronicsResponse=electronicname
               )

# if "2022" in name:
#     tray.AddModule(MedianNonCascaded, "TheMediantor",
#                    InputName="PedestalRemovedMap",
#                    # InputName="TAXIRadioWaveform",
#                    OutputName="MedianMap",
#                    ApplyInDAQ=False)

#     frameName = "MedianMap"

# else:
frameName = "PedestalRemovedMap"

# tray.AddModule(WaveformPlotter, "TheMedianWFPlotter",
#                InputName="MedianMap",
#                OutFig="waveform_CleanedCut_{0}.png".format(name))

tray.AddModule("BandpassFilter", "BoxFilter",
               InputName=frameName,
               OutputName="FilteredMap",
               FilterType=radcube.eBox,
               FilterLimits=[100*I3Units.megahertz, 300*I3Units.megahertz],
               ApplyInDAQ=False
               )

tray.AddModule("ElectronicResponseRemover", "ElectronicResponseRemover",
               # InputName="PedestalRemovedMap",
               InputName="FilteredMap",
               OutputName="ElectronicResponseRemoved",
               ElectronicsResponse=electronicname)

# tray.AddModule(WaveformPlotter, "TheRawWFPlotter",
#                InputName="TAXIRadioWaveform",
#                OutFig="waveform_RawCut_{0}.png".format(name))

# tray.AddModule(WaveformPlotter, "ThePedestalWFPlotter",
#                InputName="PedestalRemovedMap",
#                OutFig="waveform_PedestalCut_{0}.png".format(name))

# tray.AddModule(WaveformPlotter, "TheDeconvolvedWFPlotter",
#                InputName="ElectronicResponseRemoved",
#                OutFig="waveform_CleanedCut_{0}.png".format(name))


tray.AddModule(GalacticBackground, "TheGalaxyObserverDeconvolved",
               InputName="ElectronicResponseRemoved",
               # InputName="PedestalRemovedMap",
               # InputName="FilteredMap",
               SoftTrigger=softTrig,
               Output="GalOscillation_Deconvolved_{0}{1}.npz".format(name, softTrig*"_soft"),
               ApplyInDAQ=False
               )


# tray.AddModule(SpectrumPlotter, "TheSpectrumPlotter",
#                InputName="ElectronicResponseRemoved",
#                OutFig="spectrum_CleanedCut_{0}.png".format(name)
#                # OutFig="spectrum_DeconvolvedCut_{0}.png".format(name)
#                )

# tray.AddModule(SpectrumPlotter, "TheSpectrumFilteredPlotter",
#                InputName="MedianMap",
#                OutFig="spectrum_FilteredCut_{0}.png".format(name)
#                # OutFig="spectrum_DeconvolvedCut_{0}.png".format(name)
#                )

tray.Execute()

