from typing import Dict, List, Optional, Tuple, Union

import anndata
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.lines import Line2D


def boxplot(
    adata: anndata.AnnData,
    groupby: str,
    obs_key: str,
    hue: Optional[str] = None,
    palette: Optional[Union[str, List[str], Dict[str, str]]] = None,
    cmap: Optional[Dict[str, str]] = None,
    figsize: Tuple[float, float] = (10, 6),
    xlabel: Optional[str] = None,
    ylabel: Optional[str] = None,
    title: Optional[str] = None,
    legend_title: Optional[str] = None,
    show: bool = True,
    save: Optional[str] = None,
    dpi: int = 150,
    **kwargs,
) -> Optional[plt.Figure]:
    """
    Creates a box plot showing the distribution of a variable across categories.

    Parameters:
        adata (anndata.AnnData): The input AnnData object containing observations.
        groupby (str): The categorical variable in `adata.obs` used for grouping and color encoding on the x-axis.
        obs_key (str): The variable in `adata.obs` to plot on the y-axis.
        hue (str, optional): The categorical variable in `adata.obs` used for color encoding. Defaults to None.
        palette (Union[str, List[str], Dict[str, str]], optional): Seaborn palette name, list of colors, or a dictionary mapping categories to colors.
            Defaults to None.
        cmap (Optional[Dict[str, str]], optional): Dictionary mapping categories in `hue` to specific colors.
            Overrides `palette` if provided. Defaults to None.
        figsize (Tuple[float, float], optional): Figure size. Defaults to (10, 6).
        xlabel (Optional[str], optional): Label for the x-axis. Defaults to `groupby`.
        ylabel (Optional[str], optional): Label for the y-axis. Defaults to `obs_key`.
        title (Optional[str], optional): Title of the plot. Defaults to None.
        legend_title (Optional[str], optional): Title for the legend. Defaults to the `hue` variable.
        show (bool, optional): If True, shows the plot. If False, returns the figure. Defaults to True.
        save (Optional[str], optional): Path to save the figure. Defaults to None.
        dpi (int, optional): DPI for the saved figure. Defaults to 150.
        **kwargs: Additional arguments passed to `sns.boxplot`.

    Returns:
        Optional[plt.Figure]: The matplotlib Figure object if `show` is False, otherwise None.

    Raises:
        ValueError: If `groupby`, `obs_key`, or `hue` are not found in `adata.obs` or if both `palette` and `cmap` are provided.
    """
    if groupby not in adata.obs.columns:
        e = f"'{groupby}' not found in adata.obs columns."
        raise ValueError(e)
    if obs_key not in adata.obs.columns:
        e = f"'{obs_key}' not found in adata.obs columns."
        raise ValueError(e)
    if hue is not None and hue not in adata.obs.columns:
        e = f"'{hue}' not found in adata.obs columns."
        raise ValueError(e)
    if cmap is not None and palette is not None:
        e = "Specify either 'cmap' or 'palette', not both."
        raise ValueError(e)

    columns = [groupby, obs_key] + ([hue] if hue else [])
    data = adata.obs[columns].copy()

    used_palette = None
    if palette is None and cmap is None:
        if hue is not None:
            uns_palette = adata.uns.get(f"{hue}_colors", None)
            categories = data[hue].unique()
            if uns_palette is not None and isinstance(uns_palette, list):
                if len(uns_palette) >= len(categories):
                    used_palette = dict(zip(categories, uns_palette))
                else:
                    used_palette_colors = sns.color_palette("Set2", n_colors=len(categories))
                    used_palette = dict(zip(categories, used_palette_colors))
            else:
                used_palette_colors = sns.color_palette("Set2", n_colors=len(categories))
                used_palette = dict(zip(categories, used_palette_colors))
        else:
            used_palette = sns.color_palette("Set2")
    elif cmap is not None:
        used_palette = cmap
        if hue is not None:
            categories = data[hue].unique()
            missing_categories = set(categories) - set(used_palette.keys())
            if missing_categories:
                e = f"The following categories are missing in cmap: {missing_categories}"
                raise ValueError(e)
    else:
        # When palette is specified
        if hue is not None:
            categories = data[hue].unique()
            if isinstance(palette, str):
                used_palette_colors = sns.color_palette(palette, n_colors=len(categories))
                used_palette = dict(zip(categories, used_palette_colors))
            elif isinstance(palette, list):
                if len(palette) < len(categories):
                    e = "Not enough colors in the palette for the number of categories."
                    raise ValueError(e)
                used_palette = dict(zip(categories, palette))
            elif isinstance(palette, dict):
                used_palette = palette
            else:
                e = "Invalid type for 'palette'. Must be a string, list, or dictionary."
                raise ValueError(e)
        else:
            # No hue specified, use the palette directly
            if isinstance(palette, str):
                used_palette = sns.color_palette(palette)
            elif isinstance(palette, list):
                used_palette = palette
            else:
                e = "Invalid type for 'palette'. Must be a string or list when 'hue' is None."
                raise ValueError(e)

    if hue is not None and "dodge" not in kwargs:
        kwargs["dodge"] = False

    fig, ax = plt.subplots(figsize=figsize)
    sns.boxplot(
        data=data,
        x=groupby,
        y=obs_key,
        hue=hue,
        palette=used_palette if hue else None,
        ax=ax,
        showfliers=False,
        **kwargs,
    )

    ax.set_xlabel(xlabel if xlabel else groupby)
    ax.set_ylabel(ylabel if ylabel else obs_key)
    if title:
        ax.set_title(title)

    if hue is not None:
        handles, labels = ax.get_legend_handles_labels()
        # Remove duplicate labels and handles
        unique = dict(zip(labels, handles))
        labels = list(unique.keys())
        handles = list(unique.values())
        custom_handles = [
            Line2D(
                [0],
                [0],
                color=used_palette[label],
                marker="o",
                linestyle="None",
                markersize=10,
            )
            for label in labels
        ]
        if legend_title is not None:
            ax.legend(custom_handles, labels, title=legend_title, frameon=False)
        else:
            ax.legend(custom_handles, labels, title=hue, frameon=False)
    else:
        if ax.legend_ is not None:
            ax.legend_.remove()

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
