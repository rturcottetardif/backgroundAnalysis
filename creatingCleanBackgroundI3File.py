#!/bin/usr/env python3
import sys
sys.path.insert(1, '/home/rturcotte/work/modules')
from TriggerPicker import TriggerPicker
from icecube import radcube
from I3Tray import *
import argparse
parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('input', nargs='+',
                    help='input files')
args = parser.parse_args()
filename = args.input
# filename = "/data/exp/IceCube/2022/unbiased/surface/radio/V5/eventData_1645525547_2022-02-22_10-25-47.i3.gz"

tray = I3Tray()
tray.Add(
    'I3Reader', "dataReader",
    FilenameList=filename
    )
tray.AddModule(
    TriggerPicker, "TheSoftTriggerPicker",
    Trigger="Soft",
    ApplyInDAQ=True
    )
tray.Add(
    radcube.modules.RemoveTAXIArtifacts, "ArtifactRemover",
    InputName="TAXIRadioWaveform",
    OutputName="CleanedWaveform",
    BaselineValue=0,
    medianOverCascades=False,
    RemoveBinSpikes=True,
    BinSpikeDeviance=int(2**12),
    RemoveNegativeBins=True
    )
tray.AddModule(
    "I3Writer", "I3Write",
    Filename="/data/user/rturcotte/poleData/CleanedBackground_2020-12.i3.gz")
tray.Execute()
