import sys

import numpy as np
from matplotlib import pyplot as plt


def main() -> None:
    data = np.genfromtxt(
        sys.argv[1],
        delimiter=",",
        names=["x", "y"],
        skip_header=1,
    )

    plt.plot(data["x"], data["y"])

    plt.title("Delta simulation results")
    plt.xlabel("p_j_cond_i")
    plt.ylabel("acc")
    plt.savefig(sys.argv[2])
    plt.close()


if __name__ == "__main__":
    main()
