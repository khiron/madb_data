from plotly.graph_objects import Figure

species_colors = {
    "Chimpanzee": "#CC0000",
    "Gorilla":    "#E69F00",
    "Macaque":    "#0072B2"
}

def update_layout(
    fig: Figure,
    in_panel: bool = False,
    x_title: str = None,
    y_title: str = None
) -> Figure:
    """
    Applies consistent formatting to a Plotly figure for thesis visuals.

    Args:
        fig (Figure): The Plotly figure to update.
        in_panel (bool): Use smaller dimensions for panel layout.
        x_title (str | None): Optional override for x-axis title.
        y_title (str | None): Optional override for y-axis title.

    Returns:
        Figure: The updated Plotly figure.
    """

    # Apply consistent marker formatting
    for trace in fig.data:
        if hasattr(trace, "marker") and trace.marker is not None:
            if not hasattr(trace.marker, "size") or trace.marker.size is None:
                trace.marker.size = 10

            name = getattr(trace, "name", None)
            if name:
                for species, color in species_colors.items():
                    if species in name:
                        trace.marker.color = color
                        break

                if not hasattr(trace.marker, "symbol") or trace.marker.symbol is None:
                    trace.marker.symbol = "x" if "cycle" in name.lower() else "circle"


    # Set figure size
    width = 400 if in_panel else 800
    height = 400 if in_panel else 400

    # Base layout
    layout_updates = {
        "font": dict(size=14),
        "width": width,
        "height": height,
        "margin": dict(l=40, r=20, t=20, b=40),
        "xaxis": dict(
            showgrid=False,
            zeroline=False,
            ticks="outside",
            ticklen=5,
            tickwidth=1,
            title_standoff=10,
        ),
        "yaxis": dict(
            showgrid=False,
            zeroline=False,
            ticks="outside",
            ticklen=5,
            tickwidth=1,
            title_standoff=10,
        ),
        "showlegend": False,
    }

    # Override axis titles if provided
    if x_title is not None:
        layout_updates["xaxis"]["title"] = x_title
    if y_title is not None:
        layout_updates["yaxis"]["title"] = y_title

    fig.update_layout(**layout_updates)
    return fig
