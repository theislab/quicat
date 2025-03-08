from typing import List, Optional, Tuple, Union

import anndata
import matplotlib.pyplot as plt
import seaborn as sns


def barplot(
    adata: anndata.AnnData,
    groupby: str,
    order: Optional[List[str]] = None,
    color: Union[str, Tuple[float, float, float]] = "blue",
    figsize: Tuple[float, float] = (10, 6),
    xlabel: Optional[str] = None,
    ylabel: str = "Percentage (%)",
    title: Optional[str] = None,
    show: bool = True,
    save: Optional[str] = None,
    dpi: int = 150,
    **kwargs,
) -> Optional[plt.Figure]:
    """
    Creates a bar plot showing the percentage of each category in a specified variable.

    Parameters:
        adata (anndata.AnnData): The input AnnData object.
        groupby (str): Name of the variable (column in `adata.obs`) to plot.
        order (Optional[List[str]], optional): Specific order of the categories. If None,
            categories are ordered from highest to lowest percentage. Defaults to None.
        color (Union[str, Tuple[float, float, float]], optional): Color for the bars. Can be a color name, an RGB tuple, or a hex code.
            Defaults to 'blue'.
        figsize (Tuple[float, float], optional): Size of the figure. Defaults to (10, 6).
        xlabel (Optional[str], optional): Label for the x-axis. Defaults to None.
        ylabel (str, optional): Label for the y-axis. Defaults to 'Percentage (%)'.
        title (Optional[str], optional): Title of the plot. Defaults to None.
        show (bool, optional): If True, displays the plot. If False, returns the figure object. Defaults to True.
        save (Optional[str], optional): Path to save the figure. If None, the figure is not saved.
            Defaults to None.
        dpi (int, optional): Resolution of the saved figure. Defaults to 150.
        **kwargs: Additional keyword arguments to pass to `sns.barplot`.

    Returns:
        Optional[plt.Figure]: The matplotlib Figure object if `show` is False, otherwise None.
    """
    data = adata.obs[groupby].value_counts().reset_index()
    data.columns = [groupby, "count"]

    data["percentage"] = 100 * data["count"] / data["count"].sum()

    if order is not None:
        data = data.set_index(groupby).loc[order].reset_index()
    else:
        data = data.sort_values(by="percentage", ascending=False)

    fig, ax = plt.subplots(figsize=figsize)
    sns.barplot(data=data, x=groupby, y="percentage", color=color, ax=ax, **kwargs)

    ax.set_xlabel(xlabel if xlabel else groupby)
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)

    plt.tight_layout()

    if save:
        fig.savefig(save, bbox_inches="tight", dpi=dpi)

    if show:
        plt.show()
        return None
    else:
        return fig
