# ========================================================================================= #
# Author : Roxanne Turcotte
# Description : Test reading to see if the new variable is there 
# ========================================================================================= #

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

# Define the parameters here! 
# ============================
path = "/data/user/rturcotte/work/"
i3Filename = args.input.split("/")[-1].split("_")[0]

class Reader(icetray.I3Module):
    def __init__(self, ctx):
        icetray.I3Module.__init__(self, ctx)
        self.AddParameter("AddedParameter", "AddedParameter", "AddedParameter") # test parameter
    
    def Configure(self):
        self.exampleVariable= self.GetParameter('AddedParameter')

    def DAQ(self, frame):
        frame[self.exampleVariable] = dataclasses.I3VectorBool([])
 
        antennaMap = frame["Traces"] 
        log_info("....Reading new output map")
        
        for iant, antkey in enumerate(antennaMap.keys()):
            channelMap = antennaMap[antkey]
            for ichan, chankey in enumerate(channelMap.keys()):
                fft = channelMap[ichan].GetFFTData() #This container holds time series and spectrum
                timeSeries = fft.GetTimeSeries()
                print(radcube.RadTraceToPythonList(timeSeries)[1][:10])
                print(len(timeSeries))

                spec = fft.GetFrequencySpectrum()
                # print(spec[:10])
                print(len(spec))

    def Finish(self):
        log_info("Finished")

tray = I3Tray()

tray.Add('I3Reader', "TheReader",
        FilenameList = [args.input]
        )

tray.AddModule(Reader, "MyReader")

tray.Execute(5)
