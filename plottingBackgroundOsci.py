import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.dates
from datetime import datetime
import pandas as pd
from scipy.optimize import leastsq
from icecube.icetray import I3Units

from astropy.visualization import astropy_mpl_style, quantity_support
plt.style.use(astropy_mpl_style)
quantity_support()

import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, Angle

import matplotlib.pyplot as plt

SMALL_SIZE = 18
MEDIUM_SIZE = 22
BIGGER_SIZE = 24

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
#plt.rcParams["font.family"] = "cursive"

year="2021"
dataDir = "/data/user/rturcotte/analysis/background/comparisonDAQs/"
plotDir = "/data/user/rturcotte/analysis/background/comparisonDAQs/"

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('inputName', type=str, default="", help='npz with rms of traces')
args = parser.parse_args()

#################################################
#
#           Let's follow SgrA
#
#################################################
# zenGain_filename = dataDir + "antennaZenithGain.npy"
zenGain_filename = "/data/user/rturcotte/analysis/background/data/antennaZenithGain.npy"


def readAntennaZenithResponse(filename):
    return np.load(filename)


def getGalacticCenter(startTime='2022-01-26', endTime='2022-01-30', nbPoints=500):
    southPole = EarthLocation("89d59m24s", "-63d27m11s", 2800 * u.m)
    sgrA = SkyCoord.from_name('sgrA')
    # time = Time(['2020-4-17 23:00:00', '2020-4-18 4:00:00'])
    # sgrAaltaz = sgrA.transform_to(AltAz(obstime=time, location=southPole))
    start = pd.Timestamp(startTime)
    end = pd.Timestamp(endTime)
    t = np.linspace(start.value, end.value, nbPoints)
    t = pd.to_datetime(t)
    april = Time(t)
    frame_timeWindow = AltAz(obstime=april, location=southPole)
    galacticCenterData = sgrA.transform_to(frame_timeWindow)
    #sgrAairmasss_July13night = sgrAaltazs_July13night.secz
    return t, galacticCenterData


def convoluteAntennaAndZenith(antennaResponse, galacticCenterData):
    GainPerZen = antennaResponse[galacticCenterData.zen.astype(int)]
    return GainPerZen/np.mean(GainPerZen)*np.sin(galacticCenterData.az)


def plotGalacticCenter(ax, t, galacticCenterToPlot, color="slateblue", label="Sgr A*"):
    ax.plot_date(t, galacticCenterToPlot, marker="", ls="-", lw=4, c=color, label=label)


class BackgroundOscillation():
    def __init__(self, filename):
        self.dataDir = dataDir
        self.plotDir = plotDir
        self.ant = 3
        self.pol = 2

        self.time = []
        self.rms = []
        # self.data = []
        self.window = 300
        self.rmsMax = 30
        self.winType = None
        self.df = pd.DataFrame()
        # self.startTIme = self.df["time"][0]
        # self.startTIme = self.df["time"][-1]

        self.readNpz(self.dataDir + filename)
        print("reading the data ... from ", self.dataDir + filename)

    def readNpz(self, filename):
        data = np.load(filename, allow_pickle=True)
        self.time = data["time"]
        self.rms = np.array([
            data["rms10"],
            data["rms11"],
            data["rms20"],
            data["rms21"],
            data["rms30"],
            data["rms31"]])
        self.power = np.array([
            data["power10"],
            data["power11"],
            data["power20"],
            data["power21"],
            data["power30"],
            data["power31"]])
        # self.data = np.array([data["time"], data["rms10"], data["rms11"], data["rms20"], data["rms21"], data["rms30"], data["rms31"]])

    def processData(self):
        self.toPandas()
        self.movingAverage()
        print("processing the data...")

    def setMaxRms(self, value):
        self.rmsMax = value

    def setWindow(self, value):
        self.window = value

    def setWinType(self, value):
        self.winType = value

    def cleanData(self, n=5):
        for iant in range(self.ant):
            for ich in range(self.pol):
                idx = "rms{0}{1}".format(1+iant, ich)
                self.df.loc[bg.df[idx] > np.min(bg.df[idx]) + n] = float("NaN")

    def setTimeWindow(self, startTime, endTime):
        self.startTime = startTime
        self.endTime = endTime
        self.df = self.df.sort_values("time", axis=0, ascending=True)
        self.df = self.df.loc[bg.df['time'] > startTime]
        self.df = self.df.loc[bg.df['time'] < endTime]

    def toPandas(self):
        self.df["time"] = self.time
        for iant in range(self.ant):
            for ich in range(self.pol):
                idx = "rms{0}{1}".format(1+iant, ich)
                self.df[idx] = self.rms[2*iant + ich]
                idx = "power{0}{1}".format(1+iant, ich)
                self.df[idx] = self.rms[2*iant + ich]

    def movingAverage(self):
        for iant in range(self.ant):
            for ich in range(self.pol):
                idx = "{0}{1}".format(1+iant, ich)
                self.df["average"+idx] = self.df["rms"+idx].rolling(self.window, win_type=self.winType, min_periods=1).mean()

    def fitSinus(self, time, rms, phase=20):
        guess_mean = np.mean(rms)
        guess_std = np.std(rms)  # 3*np.std(rms)/(2**0.5)/(2**0.5)
        guess_phase = 1
        guess_freq = 0.00002314814 * 6
        guess_amp = guess_std / np.sqrt(2)
        t = [val.timestamp() for val in time]
        t = np.asarray(t)

        # we'll use this to plot our first estimate. This might already be good enough for you
        data_first_guess = guess_std * \
            np.sin(t * guess_freq + guess_phase) + guess_mean

        print(data_first_guess.shape)
        print(rms.shape)
        xx = np.arange(0, len(rms))
        # # Define the function to optimize, in this case, we want to minimize the difference
        # # between the actual data and our "guessed" parameters
        def optimize_func(x): return x[0] * np.sin(x[1] * t + x[2]) + x[3] - rms
        est_amp, est_freq, est_phase, est_mean = leastsq(
            optimize_func, [guess_amp, guess_freq, guess_phase, guess_mean])[0]

        # # recreate the fitted curve using the optimized parameters
        data_fit = est_amp * np.sin(est_freq * t + est_phase) + est_mean
        print("estimated amp {0}, freq {1}, phase {2}".format(
            est_amp, est_freq, est_phase))

        # # recreate the fitted curve using the optimized parameters

        # fine_t = np.arange(0,max(t),0.1)
        # data_fit=est_amp*np.sin(est_freq*fine_t+est_phase)+est_mean
        # plt.figure(figsize=[20,12])
        #plt.plot(t, rms, '.')
        #plt.plot(t, data_first_guess, label='first guess', lw=0.5)
        return data_fit

    def plotSinusFit(self, ax, time, data_fit, c="k"):
        #plt.plot_date(time, rms, '.', color=color, alpha=0.1, label="moving average")
        #plt.plot_date(time, data_first_guess, ls="--", lw=1.5, marker="", label='first guess', color=color)
        ax.plot_date(time, data_fit, ls="-", lw=3, marker="", label='sinus fitting', color=c)

    def plotOnePol(self, ax, iant, ich):
        color = ["purple", "teal"]#"goldenrod",]
        idx = "{0}{1}".format(1+iant, ich)
        averageRmsNormalized = self.df["average"+idx] - np.mean(self.df["average"+idx])
        ax.plot_date(self.df["time"], self.df["rms"+idx]/I3Units.mV,# - np.mean(self.df["average"+idx]),
                     marker=".", c=color[ich], alpha=0.3, label="Ant. {0}, Pol. {1}".format(1+iant, ich))
        # ax.plot_date(self.df["time"], averageRmsNormalized,
        #              marker=".", c=color[ich], label="Moving average")

    def plotAll(self):
        fig = plt.figure(figsize=[20, 15])#, constrained_layout=True)
        spec = gridspec.GridSpec(ncols=1, nrows=3)
        color = ["goldenrod", "purple"]
        for iant in range(self.ant):
            for ich in range(self.pol):
                ax = fig.add_subplot(spec[iant])
                self.plotOnePol(ax, iant, ich)
                ax.tick_params(axis='y', labelcolor="k")

            ax.set_ylabel("RMS / mV")
            # ax1.set_xlim([datetime(2020, 4, 17), datetime(2020, 4, 21)])
            # ax.set_ylim(0, 0.3)
            #ax.legend(["test 1", "test 2"])
            plt.legend(loc="upper right")
            plt.setp(ax.get_xticklabels(), visible=False)
        plt.setp(ax.get_xticklabels(), visible=True)
        plt.xticks(rotation=45)
        plt.tight_layout()

        # plt.savefig(plotDir + "AllAntennas.png")
        # plt.close()

    def plotRmsAndGalaxy(self):
        zenAntennaResponse = readAntennaZenithResponse(zenGain_filename)
        t, GC = getGalacticCenter(startTime=self.startTime, endTime=self.endTime, nbPoints=len(self.df['time']))
        GC_convoluted = convoluteAntennaAndZenith(zenAntennaResponse, GC)


        fig = plt.figure(figsize=[20, 15])
        spec = gridspec.GridSpec(ncols=1, nrows=3)
        color = ["goldenrod", "purple"]
        for iant in range(self.ant):
            mean = (np.mean(self.df["average{0}0".format(1+iant)]) +  np.mean(self.df["average{0}1".format(1+iant)]))/2
            for ich in range(self.pol):
                ax1 = fig.add_subplot(spec[iant])
                self.plotOnePol(ax1, iant, ich)
                ax1.tick_params(axis='y', labelcolor="k")

            ax2 = ax1.twinx()
            phase_delay = -35*u.deg
            print(phase_delay)
            plotGalacticCenter(ax2, t, np.sin(2*GC.az), color="gray", label="Sgr A* azimuth")
            # plotGalacticCenter(ax2 , t, (GC_convoluted), color="lavender", label="Sgr A* corrected")
            ax2.tick_params(axis='y', labelcolor='goldenrod', color="goldenrod")
            ax2.set_ylabel("$\sin (\phi)$")

            ax1.set_ylabel("amplitude RMS /mV")
            ax1.set_xlim([datetime(2022, 1, 26), datetime(2022, 1, 30)])
            # ax1.set_ylim(-2.5, 2.5)
            ax2.set_ylim(-2, 2)
            # ax1.get_legend().remove()
            # ax1.legend()
            # ax2.legend()
            plt.setp(ax1.get_xticklabels(), visible=False)
        plt.setp(ax1.get_xticklabels(), visible=True)
        plt.xticks(rotation=45)
        plt.tight_layout()

    def correlationRmsGalaxy(self):
        color = ["m", "y", "c"]
        zenAntennaResponse = readAntennaZenithResponse(zenGain_filename)
        t, GC = getGalacticCenter(startTime=self.startTime, endTime=self.endTime, nbPoints=len(self.df['time']))
        GC_convoluted = convoluteAntennaAndZenith(zenAntennaResponse, GC)
        fig = plt.figure(figsize=[20, 15])
        spec = gridspec.GridSpec(ncols=1, nrows=4)
        for iant in range(self.ant):
        # for ich in range(self.pol):
            idx = "{0}{1}".format(1+iant, 0)
            ax = fig.add_subplot(spec[0])
            averageRmsNormalized = self.df["average"+idx] - np.mean(self.df["average"+idx])
            ax.scatter(abs(GC_convoluted), averageRmsNormalized, marker=",", alpha=0.2)#, c=color[iant], label="ant. {0}, pol. {1}".format(1+iant, ich))
            ax = fig.add_subplot(spec[1])
            ax.scatter(abs(np.cos(GC.az)), averageRmsNormalized, marker=",", alpha=0.2)#, c=color[iant], label="ant. {0}, pol. {1}".format(1+iant, ich))
            ax = fig.add_subplot(spec[2])
            ax.scatter(abs(np.sin(2*GC.az)), averageRmsNormalized, marker=",", alpha=0.2)#, c=color[iant], label="ant. {0}, pol. {1}".format(1+iant, ich))
        ax = fig.add_subplot(spec[3])
        theta=np.linspace(0, 360, len(averageRmsNormalized))
        ax.scatter(np.sin(theta), np.cos(theta))
        ax.scatter(np.sin(theta), np.sin(theta))
        # ax.scatter(t, abs(GC.az) - averageRmsNormalized , marker=",", alpha=0.2)#, c=color[iant], label="ant. {0}, pol. {1}".format(1+iant, ich))



    #     plt.xticks(rotation=45)
    #     plt.tight_layout()

# works only for bg_april for now
bg = BackgroundOscillation(args.inputName)
startTime = '{0}-01-26'.format(year)
endTime = '{0}-01-30'.format(year)
bg.processData()

bg.plotAll()
plt.savefig(plotDir + "movingAverage_all_{0}.png".format(year))
plt.close()


bg.setTimeWindow(startTime=startTime, endTime=endTime)
bg.plotAll()
plt.savefig(plotDir + "movingAverage_all_timeWindow_{0}.png".format(year))
plt.close()

# bg.cleanData(n=2)
# bg.movingAverage()
# bg.setTimeWindow(startTime=startTime, endTime=endTime)
# bg.plotAll()
# plt.savefig(plotDir + "movingAverage_all_cleaned.png")
# plt.close()

# bg.plotRmsAndGalaxy()
# plt.savefig(plotDir + "rmsAndGalaxy_2021.png")
# plt.close()

# bg.correlationRmsGalaxy()
# plt.savefig(plotDir + "correlationRmsAndGalaxy.png")
# plt.close()










