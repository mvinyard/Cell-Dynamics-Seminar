
import numpy as np
from tqdm.notebook import tqdm
import matplotlib.pyplot as plt
import vinplots


from ._transformers import Transformers

transform = Transformers()

class TrajectoryPlot:
    def __init__(self, adata, X_pred: np.array):

        """
        X_pred
            type: np.ndarray
            shape:(n_steps x n_cells x n_dim)
        """

        self._X_umap = adata.obsm["X_umap"]
        X_pred_ = []
        for i in tqdm(range(X_pred.shape[1])):
            print("UMAP: ({} / {})".format(i + 1, X_pred.shape[1]), end="\r")
            X_pred_.append(transform.umap.transform(X_pred[:, i, :]))

        self._X_pred = np.array(X_pred_)

    def init_plot(self):

        self.fig, self.axes = vinplots.quick_plot(figsize=1.5, rm_ticks=True)
        self.fig.modify_spines(ax="all", spines_to_delete=["bottom", "left"])

    def manifold(self):
        self.axes[0].scatter(
            self._X_umap[:, 0],
            self._X_umap[:, 1],
            c="lightgrey",
            zorder=0,
            s=10,
            alpha=0.8,
        )

    def t0(self, X_pred):

        self.axes[0].scatter(X_pred[0, 0], X_pred[0, 1], c="w", zorder=1, s=80)
        self.axes[0].scatter(
            X_pred[0, 0], X_pred[0, 1], c="red", s=50, zorder=3, label="$t_0$"
        )

    def path(self, cell=0, color="k"):

        xu = self._X_pred[cell]
        self.t0(xu)
        self.axes[0].plot(xu[:, 0], xu[:, 1], c=color, zorder=1, alpha=0.5)
        self.axes[0].scatter(xu[:, 0], xu[:, 1], c=color, zorder=2, s=10, alpha=0.5)
        self.tf(xu)

    def tf(self, X_pred):
        self.axes[0].scatter(X_pred[-1, 0], X_pred[-1, 1], c="w", zorder=1, s=80)
        self.axes[0].scatter(
            X_pred[-1, 0], X_pred[-1, 1], c="blue", s=50, zorder=3, label="$t_f$"
        )

    def __call__(self, cell=0, color="k"):

        self.init_plot()
        self.manifold()

        if isinstance(cell, int):
            cell = [cell]

        for i, c in enumerate(cell):
            self.path(cell=c, color=color)
            if i == 0:
                self.axes[0].legend(edgecolor="w")