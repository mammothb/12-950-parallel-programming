import contextlib
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from PIL import Image

sns.set_theme()


def main():
    output_dir = Path(__file__).resolve().parents[1] / "images"
    output_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv("result.csv")
    nt = df.shape[0]
    nx = df.shape[1] - 1

    x = np.linspace(0.0, 1.0, nx)

    for i in range(nt):
        fig = plt.figure()
        plt.plot(x, df.iloc[i, 1:])
        fig.savefig(output_dir / f"image_{i:06d}.png")
        plt.close()

    with contextlib.ExitStack() as stack:
        images = (
            stack.enter_context(Image.open(file))
            for file in sorted(output_dir.glob("*.png"))
        )
        image = next(images)
        image.save(
            output_dir / "image.gif",
            format="GIF",
            append_images=images,
            save_all=True,
            duration=100,
            loop=0,
        )


if __name__ == "__main__":
    main()
