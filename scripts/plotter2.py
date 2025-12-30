import pandas as pd
import matplotlib.pyplot as plt

# Path to CSV file 
csv_file = "/home/kjw2kor/Backups/Raw_Script/csv/Filtered.csv"  
# Read CSV
df = pd.read_csv(csv_file)

# Use dataset points instead of time
x = range(len(df))

# Set up 2x3 subplots
fig, axes = plt.subplots(2, 3, figsize=(15, 8), sharex=True)

# Define signals
signals = ["Ax", "Ay", "Az", "Gx", "Gy", "Gz"]

# Plot each signal
for i, signal in enumerate(signals):
    row = i // 3
    col = i % 3
    axes[row, col].plot(x, df[signal], label=signal)
    axes[row, col].set_title(signal)
    axes[row, col].grid(True)
    axes[row, col].legend()

# Common labels
fig.suptitle(csv_file, fontsize=16)
fig.supxlabel("Dataset Points (Index)")
fig.supylabel("Sensor Values")

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.show()
