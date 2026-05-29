import pandas as pd
import matplotlib.pyplot as plt
import os
import ast

target_folder = "foldername" # specify folder name in "param_estimation_result"

path_base = os.path.join("param_search", "param_estimation_result", target_folder)
# Get a SEED_num
filename = os.path.join(path_base, "hyperparams.csv")
df = pd.read_csv(filename)
SEED_num = ast.literal_eval(df["GA"].iloc[0])["SEED_num"]

# plot
path_history = os.path.join(path_base, "GA_loss_history")
fig1 = plt.figure()
ax1 = fig1.add_subplot(1,1,1)
cmap = plt.get_cmap('jet')
for seed in range(SEED_num):
    filename = os.path.join(path_history, f"seed_{seed:03d}.csv")
    df = pd.read_csv(filename)
    color = cmap(seed / (SEED_num-1))
    ax1.plot(df["gen"], df["min"], linewidth=0.2, color=color, alpha=0.5, rasterized=True)

sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=SEED_num-1))
plt.colorbar(sm, ax=ax1)

filename = os.path.join(path_base, "history1.png")
plt.savefig(filename, dpi=1200, bbox_inches="tight")

fig1 = plt.figure()
ax1.set_ylim(22.,25.)
cmap = plt.get_cmap('jet')
for seed in range(SEED_num):
    filename = os.path.join(path_history, f"seed_{seed:03d}.csv")
    df = pd.read_csv(filename)
    color = cmap(seed / (SEED_num-1))
    ax1.plot(df["gen"], df["min"], linewidth=1., color=color, alpha=0.5, rasterized=True)

sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=SEED_num-1))
plt.colorbar(sm, ax=ax1)

filename = os.path.join(path_base, "history2.png")
plt.savefig(filename, dpi=1200, bbox_inches="tight")
plt.show()