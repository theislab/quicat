from typing import Dict, List, Optional, Tuple, Union

import anndata
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def stacked_barplot(
    adata: anndata.AnnData,
    groupby: str,
    obs_key: str,
    palette: Optional[Union[str, List[str]]] = None,
    cmap: Optional[Dict[str, str]] = None,
    figsize: Tuple[float, float] = (10, 6),
    xlabel: Optional[str] = None,
    ylabel: str = "Percentage (%)",
    title: Optional[str] = None,
    legend_title: Optional[str] = None,
    show: bool = True,
    save: Optional[str] = None,
    dpi: int = 150,
    **kwargs,
) -> Optional[plt.Figure]:
    """
    Creates a stacked bar plot showing the percentage distribution of one categorical variable over another.

    Parameters:
        adata (anndata.AnnData): The input AnnData object containing observations.
        groupby (str): The categorical variable in `adata.obs` used for grouping on the x-axis.
        obs_key (str): The categorical variable in `adata.obs` to stack.
        palette (Union[str, List[str]], optional): Seaborn palette or list of colors. Defaults to None.
        cmap (Optional[Dict[str, str]], optional): Dictionary mapping categories in `obs` to specific colors.
            Overrides `palette` if provided. Defaults to None.
        figsize (Tuple[float, float], optional): Figure size. Defaults to (10, 6).
        xlabel (Optional[str], optional): Label for the x-axis. Defaults to `groupby`.
        ylabel (str, optional): Label for the y-axis. Defaults to 'Percentage (%)'.
        title (Optional[str], optional): Title of the plot. Defaults to None.
        legend_title (Optional[str], optional): Title for the legend. Defaults to the `obs` variable.
        show (bool, optional): If True, shows the plot. If False, returns the figure. Defaults to True.
        save (Optional[str], optional): Path to save the figure. Defaults to None.
        dpi (int, optional): DPI for the saved figure. Defaults to 150.
        **kwargs: Additional arguments passed to `plt.bar`.

    Returns:
        Optional[plt.Figure]: The matplotlib Figure object if `show` is False, otherwise None.

    Raises:
        ValueError: If `groupby` or `obs` are not found in `adata.obs` or if both `palette` and `cmap` are provided.
    """
    if groupby not in adata.obs.columns:
        e = f"'{groupby}' not found in adata.obs columns."
        raise ValueError(e)
    if obs_key not in adata.obs.columns:
        e = f"'{obs_key}' not found in adata.obs columns."
        raise ValueError(e)
    if cmap is not None and palette is not None:
        e = "Specify either 'cmap' or 'palette', not both."
        raise ValueError(e)

    data = adata.obs[[groupby, obs_key]]
    pivot_table = data.pivot_table(index=groupby, columns=obs_key, aggfunc="size", fill_value=0)
    relative_abundances = pivot_table.div(pivot_table.sum(axis=1), axis=0) * 100

    fig, ax = plt.subplots(figsize=figsize)

    if palette is None and cmap is None:
        uns_palette = adata.uns.get(f"{obs_key}_colors", None)
        if uns_palette is not None and isinstance(uns_palette, list):
            used_palette = uns_palette
        else:
            used_palette = sns.color_palette("Set2", n_colors=len(relative_abundances.columns))

    elif cmap is not None:
        categories = relative_abundances.columns
        missing_categories = set(categories) - set(cmap.keys())
        if missing_categories:
            e = f"The following categories are missing in cmap: {missing_categories}"
            raise ValueError(e)
        used_palette = cmap
    else:
        used_palette = list(sns.color_palette(palette, n_colors=len(relative_abundances.columns)))

    bottom = pd.Series([0] * relative_abundances.shape[0], index=relative_abundances.index)
    for i, col in enumerate(relative_abundances.columns):
        color = used_palette[i] if isinstance(used_palette, list) else used_palette[col]
        ax.bar(relative_abundances.index, relative_abundances[col], bottom=bottom, color=color, label=col, **kwargs)
        bottom += relative_abundances[col]

    ax.set_xlabel(xlabel if xlabel else groupby)
    ax.set_ylabel(ylabel)

    if title:
        ax.set_title(title)

    _, labels = ax.get_legend_handles_labels()
    if isinstance(used_palette, list):
        custom_handles = [
            mlines.Line2D([], [], color=color, marker="o", linestyle="None", markersize=10) for color in used_palette
        ]
    elif isinstance(used_palette, dict):
        custom_handles = [
            mlines.Line2D([], [], color=color, marker="o", linestyle="None", markersize=10)
            for color in used_palette.values()
        ]
    ax.legend(
        handles=custom_handles,
        labels=labels,
        title=legend_title if legend_title else obs_key,
        bbox_to_anchor=(1.05, 1),
        loc="upper left",
        alignment="left",
        frameon=False,
    )

    plt.tight_layout()
    if save:
        fig.savefig(save, bbox_inches="tight", dpi=dpi)
    if show:
        plt.show()
        return None
    else:
        return fig
