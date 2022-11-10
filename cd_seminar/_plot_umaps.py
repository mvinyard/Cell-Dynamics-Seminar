import vinplots
import matplotlib.pyplot as plt


def set_z(group, n):
    if group == "Undifferentiated":
        return 0
    else:
        return n


def to_list(item):
    if isinstance(item, list):
        return item
    if isinstance(item, str):
        return [item]


def plot_groupby(
    adata, ax=None, groupby="Cell type annotation", color_dict={}, legend_loc=None
):
    if not ax:
        ax = plt
    grouped = adata.obs.copy().groupby(groupby)
    for n, (group, group_df) in enumerate(grouped):
        xu = adata[group_df.index].obsm["X_umap"].toarray()
        ax.scatter(
            xu[:, 0],
            xu[:, 1],
            c=color_dict[group],
            s=5,
            alpha=0.5,
            zorder=set_z(group, n),
            label=group,
        )

    ax.legend(edgecolor="w", markerscale=2, loc=legend_loc)


def overview_plot(adata, groupby: list([str, str]), ncols: int = 4):

    groupby = to_list(groupby)
    nplots = len(groupby)

    if nplots < ncols:
        ncols = nplots

    if ncols == 2:
        legend_locs = [1, (1, 0.5)]

    fig, axes = vinplots.quick_plot(
        nplots=nplots,
        ncols=ncols,
        figsize=1.2,
        wspace=0.1,
        spines_to_delete=["top", "bottom", "left", "right"],
        rm_ticks=True,
    )

    group_colors = adata.uns["cmaps"]

    for n, (group, cmap) in enumerate(group_colors.items()):
        plot_groupby(
            adata, ax=axes[n], groupby=group, color_dict=cmap, legend_loc=legend_locs[n]
        )
        axes[n].set_title(group, fontsize=14)