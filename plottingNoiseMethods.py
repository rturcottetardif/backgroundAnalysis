import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

noiseWindow = 64
what = "medianHalf"
path = "/data/user/rturcotte/analysis/background/noiseMethod/"
filename = "noiseCalculation_{1}_{0}.npz".format(noiseWindow, what)
# filename = "noiseCalculation_CleanedBackground_2020-12_64.npz"


def readOutData(file):
    data = np.load(file, allow_pickle=True)
    return data["time"], data["noiseStd"], data["noiseStep"]


time, noise_std, noise_step = readOutData(path+filename)
colors = ["c", "b", "m", "r", "y", "g"]


# Plotting the Noise Level over time
fig, axs = plt.subplots(
    figsize=[16, 8], nrows=2, ncols=1, tight_layout=True,
    sharex=True, sharey=True
    # gridspec_kw={'width_ratios': [2, 1]}
    )

for iant in range(3):
    for ich in range(2):
        ax = axs[0]
        ax.plot_date(
            time, noise_std[:, iant, ich], marker=".", alpha=0.2,
            color=colors[2*iant + ich]
            )
        ax = axs[1]
        ax.plot_date(
            time, noise_step[:, iant, ich], marker=".", alpha=0.2,
            color=colors[2*iant + ich]
            )
for ax in axs:
    ax.set_ylabel("amplitude / ADC")
    # ax.set_xlim(datetime(2022, 2, 22, 12, 00, 00), datetime(2022, 2, 22, 13, 00, 00))
    # ax.set_xlim(datetime(2020, 12, 21, 15, 00, 00), datetime(2020, 12, 21, 17, 1, 00))
    ax.set_ylim(top=300)
axs[0].legend(["window method"])
axs[1].legend(["subtrace method"])
fig.savefig(path+"noiseOverTime_{1}_{0}.png".format(noiseWindow, what))


# Plotting the Histograms of the noise level
fig, ax = plt.subplots(
    figsize=[10, 5], nrows=1, ncols=1, tight_layout=True
    # gridspec_kw={'width_ratios': [2, 1]}
    )
for iant in range(3):
    for ich in range(2):
        ax.hist(
            noise_step[:, iant, ich], histtype="stepfilled", alpha=0.45,
            color=colors[2*iant + ich], bins="fd",
            label="subtraces - ant.{0}, pol.{1}".format(iant+1, ich),
            density=True
            )
        ax.hist(
            noise_std[:, iant, ich], histtype="step", alpha=1,
            color=colors[2*iant + ich], bins="fd",
            label="400ns window - ant.{0}, pol.{1}".format(iant+1, ich),
            density=True
            )
ax.set_xlabel("amplitude / ADC")
ax.set_ylabel("normilized counts")
ax.set_xlim(right=300)
ax.legend(ncol=2)
fig.savefig(path+"histNoise_{1}_{0}.png".format(noiseWindow, what))

# Comparing window sizes
windowSizes = [1, 8, 16, 32, 64, 128, 256, 512]

fig, axs = plt.subplots(
    figsize=[10, 5], nrows=1, ncols=2, tight_layout=True,
    gridspec_kw={'width_ratios': [1, 1]}
    )

for w, win in enumerate(windowSizes):
    filename = "noiseCalculation_{1}_{0}.npz".format(win, what)
    time, noise_std, noise_step = readOutData(path+filename)
    ax = axs[0]
    for iant in range(3):
        for ich in range(2):
            ax.scatter(
                win, np.std(noise_step[:, iant, ich]), marker="s",
                color=colors[2*iant + ich],
                )
            if w == 0:
                ax.axhline(
                    np.std(noise_std[:, iant, ich]),
                    color=colors[2*iant + ich], ls=":"
                    )

    # ax.legend(["ant.1 pol.0", "ant.1 pol.1", "ant.2 pol.0",
    #           "ant.2 pol.1", "ant.3 pol.0", "ant.3 pol.1"], ncol=3)
    ax.set_ylabel("standard deviation / ADC")
    ax.set_xlabel("subtraces length size / ns")

    ax = axs[1]
    for iant in range(3):
        for ich in range(2):
            ax.errorbar(
                win, np.mean(noise_step[:, iant, ich]), np.std(noise_step[:, iant, ich])/2, marker="s",
                color=colors[2*iant + ich],
                )
            if w == 0:
                ax.axhline(
                    np.mean(noise_std[:, iant, ich]),
                    color=colors[2*iant + ich], ls=":"
                    )
    # ax.legend(["ant.1 pol.0", "ant.1 pol.1", "ant.2 pol.0",
              # "ant.2 pol.1", "ant.3 pol.0", "ant.3 pol.1"], ncol=3)
    ax.set_ylabel("average / ADC")
    ax.set_xlabel("subtraces length size / ns")
fig.savefig(path + "comparingWindows_{0}.png".format(what))
