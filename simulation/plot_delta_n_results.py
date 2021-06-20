import os
import sys
from bisect import bisect_left
from glob import glob

import numpy as np
from matplotlib import pyplot as plt


def main() -> None:
    res = []
    step = 0
    for path in glob("delta_results/delta_results_*.csv"):
        n = int(os.path.splitext(os.path.basename(path))[0].split("_")[-1])
        data = np.genfromtxt(
            path,
            delimiter=",",
            names=["x", "y"],
            skip_header=1,
        )
        i = bisect_left(data["y"], 0.95)
        res.append((n, data["x"][i]))
        step = data["x"][1] - data["x"][0]
    res.sort(key=lambda x: x[0])
    x, y = zip(*res)
    plt.plot(x, y)

    plt.title(f"Delta simulation n results, step={round(step, 3)}")
    plt.xlabel("n_seq")
    plt.ylabel("q_point")

    plt.savefig(sys.argv[1])
    plt.close()


if __name__ == "__main__":
    main()
