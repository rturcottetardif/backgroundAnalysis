import numpy as np
from icecube.icetray import I3Units
from icecube.taxi_reader import taxi_tools
from I3Tray import *
from icecube import icetray, dataio, dataclasses, taxi_reader, radcube
from icecube.icetray.i3logging import *

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('input', type=str, default=[], help='Input data files. (assumed to be I3 unless --istaxifile is set')
parser.add_argument('--istaxifile', action='store_true', help='Use this flag if you are inputting the .bin files directly from taxi')
args = parser.parse_args()

filename = args.input.split("/")[-1].split("_")[1]
path = "/data/user/rturcotte/work/i3files/means/"

def convertFile(filename):
    raw_data, time, roi = get_waveform_array(filename)
    np.savez(filename[:-4]+'.npz', data=raw_data.astype(np.int16), time=time.astype(np.int64), roi=roi.astype(np.int16))

def convertToDb(amps):
    amps = 20 * np.log10(abs(np.array(amps)))
    amps[0] = 0
    return amps
    # thisFreqs = []
    # thisAmps = []
    # for iamp, amp in enumerate(amps):
    #     if amp > 0:
    #         thisFreqs.append(freqs[iamp])
    #         thisAmps.append(20 * np.log10(amp))

def baselineZero(wf):
    return wf - np.mean(wf)



class Averager(icetray.I3Module):
    def __init__(self, ctx):
        icetray.I3Module.__init__(self, ctx)
        self.timeUnit = I3Units.nanosecond
        self.freqUnit = I3Units.hertz * 10**6
    
    def Configure(self):
        # This function run before looking at the frames
        self.spectrums = []
        self.timeSeries = []
        self.freq = []
        print("begin ...")

    def DAQ(self, frame):

        if "TAXIRadioWaveform" in frame:
            antennaMap = frame["TAXIRadioWaveform"]  #This is a container that holds all the antenna info
        else:
            fatal_log("No waveform !!")


        for iant, antkey in enumerate(antennaMap.keys()):
            channelMap = antennaMap[antkey]
            for ichan, chankey in enumerate(channelMap.keys()):

                fft = channelMap[ichan].GetFFTData() #This container holds time series and spectrum
                
                # stored as 'I3RadVector' classes
                timeSerie = fft.GetTimeSeries()
                spectrum = fft.GetFrequencySpectrum()

                # Hilbert of the mean of the waveform
                df = spectrum.binning
                
                #Baseline to zero
                # for i in range(4):
                #     subset = timeSerie.GetSubset(i, 1024*i)
                #     subset = subset - radcube.GetMean(subset)
                #     wf = wf.append(subset)

                #Cascading
                #Convert to dB
                #Convert to mV -> Radcube function

                #Convert to python list
                pyTimeSerie = [timeSerie[ibin] for ibin in range(len(timeSerie))]
                pySpectrum = [spectrum[ibin] for ibin in range(len(spectrum))]
                #pythonSpectrum = ConvertToDb(pythonSpectrum)
                
                self.spectrums.append(pySpectrum)
                self.timeSeries.append(pyTimeSerie)
        self.freq = df*np.arange(0, len(spectrum))
        #print("spectrums", len(self.spectrums))
        #print("timeSeries", len(self.timeSeries))


    def Finish(self):
        log_info("Finished")
        print("finish")
        meanSpectrum = np.mean(self.spectrums, axis=0)
        meanTimeSeries = np.mean(self.timeSeries, axis=0)
        freq = self.freq
        plt.plot(freq, meanSpectrum)
        np.savez(path + 'average_' + filename + '.npz', spectrum=meanSpectrum, timeSerie=meanTimeSeries.astype(np.int16), freq=freq.astype(np.int16))



                
tray = I3Tray()

if args.istaxifile:
      tray.AddModule(taxi_reader.modules.TaxiToI3Converter.TaxiToI3Converter, "Framer",
                IsSteering = True,
                DataFile = args.input
               )  
else:
    tray.Add('I3Reader', "MyReader",
            FilenameList = [args.input]
            )


tray.AddModule(Averager, "TheAverager")

tray.Execute(10)

