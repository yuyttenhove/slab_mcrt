import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


def main():
    fname = sys.argv[1]
    basename = fname.split(".")[0]
    output_dir = Path(__file__).parent / "output"

    data = np.loadtxt(str(output_dir / fname))
    mu = data[:, 0]
    intensity = data[:, 1]

    fig, ax = plt.subplots()
    ax.semilogy(mu, intensity)
    plt.savefig(f"{output_dir / basename}.png")
    plt.show()
    print("finished")


if __name__ == "__main__":
    main()
