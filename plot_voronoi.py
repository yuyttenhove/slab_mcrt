from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
import seaborn as sns
from parse import parse
import pandas as pd
import numpy as np
from phd_python_scripts.utils.plotting_utils import subplots


def clean_lines(lines):
    return [l for l in lines.split("\n") if len(l) > 0]


def read_points_from_lines(lines):
    df = pd.DataFrame(index=range(len(lines)), columns=["x", "y"], dtype=float)
    for line in lines:
        res = parse("{}\t({}, {})", line)
        i = int(res[0])
        df.loc[i, "x"] = float(res[1])
        df.loc[i, "y"] = float(res[2])
    return df


def plot_tesselation(vertices, cells, centroids):
    _, ax = subplots(figsize=(10, 10))
    # ax = sns.scatterplot(data=centroids, x="x", y="y", s=16, color="red", ax=ax, marker="P")
    # ax = sns.scatterplot(data=vertices, x="x", y="y", s=8, color="blue", ax=ax)
    cells_xy_list = [vertices.values[cell_idx, :] for cell_idx in cells]
    patches = [plt.Polygon(xy, closed=True, edgecolor="blue", facecolor="none", linewidth=1) for xy in cells_xy_list]
    patches.append(
        plt.Polygon([[0, 0], [1, 0], [1, 1], [0, 1]], closed=True, edgecolor="red", facecolor="none", ls=":")
    )
    ax.add_collection(PatchCollection(patches, match_original=True))
    ax.set_xlim([-.1, 1.1])
    ax.set_ylim([-.1, 1.1])
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    ax.set_aspect("equal")
    plt.tight_layout()
    plt.savefig("voronoi.pdf")


def read_file(fname: Path):
    with open(fname, "r") as file:
        lines = file.read()
    lines = lines.split("# Vertices #")[-1]
    vertex_lines, cell_lines = lines.split("# Cells #")
    cell_lines, centroid_lines = cell_lines.split("# Centroids #")
    vertex_lines = clean_lines(vertex_lines)
    cell_lines = clean_lines(cell_lines)
    centroid_lines = clean_lines(centroid_lines)

    vertex_df = read_points_from_lines(vertex_lines)
    centroid_df = read_points_from_lines(centroid_lines)

    cell_list = []
    for triangle_line in cell_lines:
        res = parse("{}\t({})", triangle_line)
        vertices = np.array([int(i) for i in res[1].split(", ")])
        cell_list.append(vertices)

    return vertex_df, cell_list, centroid_df


def main():
    base_path = Path(__file__).parent
    fname = base_path / "vor_test.txt"
    vertices, cells, centroid_df = read_file(fname)

    plot_tesselation(vertices, cells, centroid_df)


if __name__ == "__main__":
    main()
