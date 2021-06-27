# Call this to have the right numbers 

def getWaveformLength():
    return int(2**10)

def getWaveformLengthCascaded():
    return int(getChannelPerPolarisation()*getWaveformLength())

def getAntennas():
    return int(3)

def getBuffers():
    return int(8)

def getChannels():
    return int(2)

def getBuffersPerChannel():
    return int(getBuffers()/getChannels())
