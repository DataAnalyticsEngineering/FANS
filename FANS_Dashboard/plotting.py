import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

def plot_subplots(data1, data2, labels_x=None, labels_y=None, subplot_titles=None, title="Subplot Grid", nrows=None, ncols=None,
                  linewidth=1, markersize=4, linecolor="blue", markercolor="red", aspect_ratio=None, fontsize=12):
    """
    Plot a grid of subplots using Plotly, handling both single-component (scalar vs scalar) and multi-component data.

    Parameters:
    - data1: numpy array, first set of data to plot (e.g., strain, time) with shape (n_datapoints, n_plots)
    - data2: numpy array, second set of data to plot (e.g., stress) with shape (n_datapoints, n_plots)
    - labels_x: list of strings, labels for the x axes of each subplot (optional, default=None)
    - labels_y: list of strings, labels for the y axes of each subplot (optional, default=None)
    - subplot_titles: list of strings, titles for each subplot (optional, default=None)
    - title: string, title of the overall plot
    - nrows: int, number of rows in the subplot grid (optional)
    - ncols: int, number of columns in the subplot grid (optional)
    - linewidth: int, line width for the plots (optional, default=1)
    - markersize: int, size of the markers (optional, default=4)
    - linecolor: string, color of the lines (optional, default="blue")
    - markercolor: string, color of the markers (optional, default="red")
    - aspect_ratio: tuple, aspect ratio of the subplots as (x, y) (optional, default=None)
    - fontsize: int, font size for axis labels, subplot titles, and tick labels (optional, default=12)
    """
    # Validate data shapes
    if not isinstance(data1, np.ndarray) or not isinstance(data2, np.ndarray):
        raise ValueError("data1 and data2 must be numpy arrays.")

    if data1.shape[0] != data2.shape[0]:
        raise ValueError("data1 and data2 must have the same number of data points (rows).")

    if data1.shape[1] != data2.shape[1]:
        raise ValueError("data1 and data2 must have the same number of components (columns).")

    # Set the number of components based on data shape
    n_components = data1.shape[1]
    
    # If nrows or ncols is not specified, determine an optimal grid layout
    if nrows is None or ncols is None:
        nrows = int(np.ceil(np.sqrt(n_components)))
        ncols = int(np.ceil(n_components / nrows))
    
    # Handle subplot titles
    if subplot_titles is None:
        subplot_titles = [f'Component {i+1}' for i in range(n_components)]
    elif len(subplot_titles) != n_components:
        raise ValueError(f"The length of subplot_titles must match the number of components ({n_components}).")
    
    # Handle labels_x and labels_y
    if labels_x is None:
        labels_x = [""] * n_components
    elif len(labels_x) != n_components:
        raise ValueError(f"The length of labels_x must match the number of components ({n_components}).")

    if labels_y is None:
        labels_y = [""] * n_components
    elif len(labels_y) != n_components:
        raise ValueError(f"The length of labels_y must match the number of components ({n_components}).")
    
    # Create the subplot figure
    fig = make_subplots(rows=nrows, cols=ncols, subplot_titles=subplot_titles)

    # Add traces for each component
    for i in range(n_components):
        row = i // ncols + 1
        col = i % ncols + 1
        fig.add_trace(go.Scatter(
            x=data1[:, i],
            y=data2[:, i],
            mode='lines+markers',
            marker=dict(symbol='x', size=markersize, color=markercolor),
            line=dict(width=linewidth, color=linecolor),
            name=f'Component {i+1}'
        ), row=row, col=col)

        # Update axes with text labels
        fig.update_xaxes(title_text=labels_x[i], row=row, col=col, showgrid=True, mirror=True, 
                         ticks="inside", tickwidth=2, ticklen=6, title_font=dict(size=fontsize),
                         tickfont=dict(size=fontsize))
        fig.update_yaxes(title_text=labels_y[i], row=row, col=col, showgrid=True, mirror=True, 
                         ticks="inside", tickwidth=2, ticklen=6, title_font=dict(size=fontsize),
                         tickfont=dict(size=fontsize))

    # Update layout with the overall plot title and styling
    fig.update_layout(
        height=600,
        width=900,
        title_text=title,
        title_font=dict(size=fontsize),
        showlegend=False,
        template="plotly_white",
        margin=dict(l=20, r=20, t=50, b=20),
        title_x=0.5,
        autosize=False,
    )

    # Add a box outline around all subplots
    for i in range(1, nrows * ncols + 1):
        fig.update_xaxes(showline=True, linewidth=2, linecolor='black', row=(i-1)//ncols + 1, col=(i-1)%ncols + 1)
        fig.update_yaxes(showline=True, linewidth=2, linecolor='black', row=(i-1)//ncols + 1, col=(i-1)%ncols + 1)

    # Update subplot titles with the specified fontsize
    for annotation in fig['layout']['annotations']:
        annotation['font'] = dict(size=fontsize)

    # Set aspect ratio if provided
    if aspect_ratio:
        fig.update_layout(
            # xaxis=dict(scaleanchor="y", scaleratio=aspect_ratio[0]/aspect_ratio[1]),
            # yaxis=dict(scaleanchor="x", scaleratio=aspect_ratio[1]/aspect_ratio[0])
        )



    # # Increase the DPI/quality of the plots by setting higher resolution dimensions
    # fig.update_layout(
    #     height=1200,
    #     width=1800,
    #     autosize=False,
    # )

    # Show the plot
    fig.show()


    



