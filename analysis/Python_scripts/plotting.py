import textwrap
import pandas as pd
from matplotlib import pyplot as plt
from typing import List
import seaborn as sns
import plotly.graph_objs as go
import numpy as np


def make_boxplot(grouped_df: pd.DataFrame, colname: str, legend: str) -> None:
    """Create a boxplot for each group in a DataFrame.

    Args:
        grouped_df (pd.DataFrame): The DataFrame grouped by some criteria.
        colname (str): The name of the column to create boxplots for.
        legend (str): The legend for the plot.
    """
    # Create the plot with a width of 10 inches
    fig, ax = plt.subplots(figsize=(17, 5))

    # Create a list of DataFrames, one for each group
    grouped_dfs = [grouped_df.get_group(group) for group in grouped_df.groups]

    # Create a boxplot for each DataFrame
    ax.boxplot([df[colname] for df in grouped_dfs])

    labels = [f"{group} ({len(x)})" for group, x in grouped_df]
    # Set the x-axis tick labels and increase the distance between ticks to 0.5 inches
    ax.set_xticklabels(labels, rotation=90)
    ax.tick_params(axis='x', which='major', pad=0.5)

    # Set the x-axis label, y-axis label, and plot title
    # ax.set_ylim([0, 1])
    ax.set_xlabel('Class', fontsize=16)
    ax.set_ylabel(colname, fontsize=16)
    ax.legend([legend], fontsize=14)

    # Show the plot
    plt.show()


def plot_histograms_sidebyside(
    same_query_ref1: pd.DataFrame,
    same_query_ref2: pd.DataFrame,
    same_query_ref3: pd.DataFrame,
    same_query_ref4: pd.DataFrame,
    same_query_ref5: pd.DataFrame,
    column_name: str,
    xaxis_title: str = '',
    title: str = ''
) -> None:
    """Create and display a histogram for each of five input DataFrames.

    Args:
        same_query_ref1 (pd.DataFrame): The first input DataFrame.
        same_query_ref2 (pd.DataFrame): The second input DataFrame.
        same_query_ref3 (pd.DataFrame): The third input DataFrame.
        same_query_ref4 (pd.DataFrame): The fourth input DataFrame.
        same_query_ref5 (pd.DataFrame): The fifth input DataFrame.
        column_name (str): The name of the column to create histograms for.
        xaxis_title (str): The title for the x-axis.
        title (str): The title for the plot.
    """
    # Define number of bins
    n_bins = 10

    # Create data arrays
    x1 = same_query_ref1[column_name]
    x2 = same_query_ref2[column_name]
    x3 = same_query_ref3[column_name]
    x4 = same_query_ref4[column_name]
    x5 = same_query_ref5[column_name]

    # Create histogram traces
    trace1 = go.Histogram(x=x1, nbinsx=n_bins, name='All peaks ({})'.format(
        len(x1)), xbins=dict(start=0, end=1))
    trace2 = go.Histogram(x=x2, nbinsx=n_bins, name='Top 3 peaks ({})'.format(
        len(x2)), xbins=dict(start=0, end=1))
    trace3 = go.Histogram(x=x3, nbinsx=n_bins, name='Top 5 peaks ({})'.format(
        len(x3)), xbins=dict(start=0, end=1))
    trace4 = go.Histogram(x=x4, nbinsx=n_bins, name='Top 10 peaks ({})'.format(
        len(x4)), xbins=dict(start=0, end=1))
    trace5 = go.Histogram(x=x5, nbinsx=n_bins, name='Top 20 peaks ({})'.format(
        len(x5)), xbins=dict(start=0, end=1))

    # Create layout
    layout = go.Layout(title=title,
                       xaxis=dict(title=xaxis_title, range=[
                                  0, 1], dtick=0.1, tickfont=dict(size=15)),
                       yaxis=dict(title='Frequency', range=[
                                  0, 140], tickfont=dict(size=15)),
                       legend=dict(x=0.81, y=1.0),
                       font=dict(size=17))

    # Create figure
    fig = go.Figure(data=[trace1, trace2, trace3,
                    trace4, trace5], layout=layout)

    # Display the plot
    fig.show()


def plot_histogram(x: List[float], xaxis_title: str = '', title: str = '') -> None:
    """Create and display a histogram for the input data.

    Args:
        x (List[float]): The input data.
        xaxis_title (str): The title for the x-axis.
        title (str): The title for the plot.
    """
    # Define number of bins
    n_bins = 20

    # Create histogram traces
    trace1 = go.Histogram(x=x, nbinsx=n_bins, name='All peaks ({})'.format(
        len(x)), xbins=dict(start=0, end=1))

    # Create layout
    layout = go.Layout(title=title,
                       xaxis=dict(title=xaxis_title, range=[
                                  0, 1], dtick=0.1, tickfont=dict(size=15)),
                       yaxis=dict(title='Frequency', range=[
                                  0, 140], tickfont=dict(size=15)),
                       legend=dict(x=0.81, y=1.0),
                       font=dict(size=17))

    # Create figure
    fig = go.Figure(data=[trace1], layout=layout)

    # Display the plot
    fig.show()


def create_plot(df: pd.DataFrame,
                grouping_column: str,
                xlabel: str = None,
                showlegend: bool = True,
                normalized_matches: bool = True,
                nist_scale: bool = True,
                hide_labels: bool = False,
                order: List[str] = None,
                colors : List[str] = ['yellow', 'deepskyblue'],
                median_colors : List[str] = ['darkgreen', 'blue']) -> plt.Figure:
    """ Create a boxplot with two y-axes for the input DataFrame.

    Args:
        df (pd.DataFrame): The input DataFrame.
        grouping_column (str): The name of the column to group by.
        xlabel (str): The label for the x-axis.
        showlegend (bool): Whether to show the legend.
        normalized_matches (bool): Whether to normalize the matches.
        nist_scale (bool): Whether to use the NIST scale.
        hide_labels (bool): Whether to hide the labels.

    Returns:
        fig (plt.Figure): The plot.
    """
    matches_col, scores_col, label_fontsize, tick_fontsize, text_width, flierprops = init()

    # Set the style
    sns.set_palette(sns.color_palette(colors))
    sns.set_style(style='white')

    # Calculate the number of bars and figure size
    plot_width = calc_plot_width(df, grouping_column, 2)
    fig = plt.figure(figsize=(plot_width, 5))

    ax = sns.boxplot(x=grouping_column,
                     y="value",
                     hue="Number",
                     data=df,
                     hue_order=[matches_col, np.nan],
                     medianprops={'color': median_colors[0], 'linewidth': 4.0},
                     flierprops={'marker': 'o', 'markersize': 10, 'markerfacecolor': 'none'},
                     order=order)
    ax.legend_.remove()

    rescale_matches_axes(df, normalized_matches, matches_col, label_fontsize, tick_fontsize, ax)

    ax2 = ax.twinx()

    sns.boxplot(ax=ax2, x=grouping_column, y='value', hue='Number',
                data=df, hue_order=[np.nan, scores_col],
                medianprops={'color': median_colors[1], 'linewidth': 4.0},
                flierprops={'marker': 'o', 'markersize': 10, 'markerfacecolor': 'none'})

    ax2.legend_.remove()

    rescale_scores_axes(nist_scale, label_fontsize, tick_fontsize, ax2)

    if xlabel:
        # Set x-axis label and font size
        ax.set_xlabel(xlabel, fontsize=label_fontsize)
    else:
        ax.set_xlabel("", fontsize=0)
        ax2.set_xlabel("", fontsize=0)

    # Create a count for each x-axis label
    count_data = df[grouping_column].value_counts()

    # # Add the count labels to the x-axis
    make_xticklabels(text_width, ax, count_data)

    # Change the legend labels
    if showlegend:
        handles, labels = ax.get_legend_handles_labels()
        labels[0] = ax.get_ylabel()
        labels[1] = ax2.get_ylabel()
        ax.legend(handles, labels, loc='upper right', fontsize=label_fontsize)

    if hide_labels:
        ax.set_ylabel("", fontsize=0)
        ax2.set_ylabel("", fontsize=0)

    return fig

def init():
    matches_col = 'matches'
    scores_col = 'scores'

    label_fontsize = 22
    tick_fontsize = 22
    text_width = 21
    flierprops={
        'marker': 'o',
        'markersize': 10,
        'markerfacecolor': 'none'
    }
    return matches_col,scores_col,label_fontsize,tick_fontsize,text_width, flierprops


def rescale_scores_axes(nist_scale, label_fontsize, tick_fontsize, ax2):
    if nist_scale:
        ax2.set_ylim([-100, 1100])
    else:
        ax2.set_ylim([-0.1, 1.1])  # Set y-axis limits

    ax2.set_ylabel('scores', fontsize=label_fontsize)  # Set y-axis label
    # Set font size of y-axis tick labels
    ax2.tick_params(axis='y', labelsize=tick_fontsize)
    ax2.tick_params(axis='x', labelsize=tick_fontsize)


def rescale_matches_axes(df, normalized_matches, matches_col, label_fontsize, tick_fontsize, ax):
    if normalized_matches:
        ax.set_ylim([-10, 110])  # Set y-axis limits
        ax.set_ylabel('ions matching reference', fontsize=label_fontsize)  # Set y-axis label
        ax.set_yticklabels([str(f"{int(x)}%") for x in ax.get_yticks()])
    else:
        ax.set_ylabel('ion matches', fontsize=label_fontsize)
        top = df.loc[df['Number'] == matches_col, 'value'].max()
        if top <= 6:
            ax.set_ylim(-0.5, 5.5)
        else:
            ax.set_ylim(0 - 0.1 * top, top=1.1 * top)  # Set y-axis limits

    # Set font size of x-axis tick labels
    ax.tick_params(axis='x', labelsize=tick_fontsize)
    # Set font size of y-axis tick labels
    ax.tick_params(axis='y', labelsize=tick_fontsize)


def make_xticklabels(text_width, ax, count_data, ha='right', rotation=45):
    xlabels = [label.get_text() for label in ax.get_xticklabels()]
    xlabels = ['\n'.join(textwrap.wrap(
        label + ' (' + str((count_data.loc[label] // 2)) + ')', width=text_width)) for label in xlabels]
    ax.set_xticklabels(xlabels, rotation=rotation, ha=ha)


def scatterplot_matplotlib(df: pd.DataFrame) -> plt.Figure:
    """ Create a scatterplot with the input DataFrame.

    Args:
        df (pd.DataFrame): The input DataFrame.

    Returns:
        fig (plt.Figure): The plot.
    """
    _, _, label_fontsize, tick_fontsize, _, _ = init()

    fig = plt.figure(figsize=(24, 8))
    scatter = plt.scatter(
        df['scores'],
        df['matches'],
        # Adjust the size scaling factor as needed
        s=df['FractionQuery'] * 200,
        c=df['FractionReference'] * 100,
        cmap='viridis',  # change the colorscale as needed
        alpha=0.5,
        vmin=0,
        vmax=100
    )
    for ax in fig.get_axes():
        ax.tick_params(axis='x', labelsize=tick_fontsize)
        ax.tick_params(axis='y', labelsize=tick_fontsize)

    cbar = plt.colorbar(scatter)
    cbar.set_label('ions matching reference (%)', size=label_fontsize)
    cbar.ax.tick_params(labelsize=tick_fontsize)
    plt.xlabel('scores', fontsize = label_fontsize)
    plt.ylabel('ion matches', fontsize = label_fontsize)

    # Add a legend for the size
    sizes = [1, 50, 100]
    for size in sizes:
        plt.scatter([], [], c='c', alpha=0.5, s=size * 2, label=str(size))
    plt.legend(scatterpoints=1, title='ions matching query (%)', title_fontsize=label_fontsize,
               labelspacing=1, loc='upper left', fontsize = tick_fontsize, ncols=3)
    return fig

def create_dual_plot(df1: pd.DataFrame, df2: pd.DataFrame,
                     grouping_column: str,
                     order1: List[str] = None,
                     order2: List[str] = None
                     ) -> plt.Figure:
    """ Create a boxplot with two y-axes for the input DataFrame.

    Args:
        df1, df2 (pd.DataFrame): The input DataFrames.
        grouping_column (str): The name of the column to group by.
        order1, order2 (List[str]): The order of the bars.
    Returns:
        fig (plt.Figure): The plot.
    """
    medianprops={
        'linewidth': 2.0
    }

    # Set the style
    sns.set_style(style='white')
    sns.set_palette(sns.color_palette(['salmon', 'crimson']))

    # Calculate the number of bars and figure size
    plot_width = calc_plot_width(df1, grouping_column, 1.5)

    # Create the figure and axes
    fig, ax1 = plt.subplots(figsize=(plot_width, 5))

    # Create the boxplot for the first DataFrame
    ax1, ax1_2 = create_boxplot_subplot(
        df1,
        grouping_column,
        order1,
        medianprops,
        ax1,
        'right',
    )

    # Create the second axes for the second DataFrame
    ax2 = ax1.twiny()
    sns.set_palette(sns.color_palette(['skyblue', 'deepskyblue']))

    # Create the boxplot for the second DataFrame
    ax2, ax2_2 = create_boxplot_subplot(
        df2,
        grouping_column,
        order2,
        medianprops,
        ax2,
        'left',
    )

    lines, labels = ax1.get_legend_handles_labels()
    labels[0] = ax1.get_ylabel() + ' (with N)'
    labels[1] = ax1_2.get_ylabel() + ' (with N)'
    lines2, labels2 = ax2.get_legend_handles_labels()
    labels2[0] = ax2.get_ylabel()
    labels2[1] = ax2_2.get_ylabel()
    ax1.legend(lines + lines2, labels + labels2, loc='upper right', ncols=2)

    ax1_2.get_legend().remove()
    ax2_2.get_legend().remove()
    ax2.get_legend().remove()
    

    return fig

def calc_plot_width(df, grouping_column, bar_width):
    n_bars = df[grouping_column].nunique()
    plot_width = n_bars * bar_width
    return plot_width

def create_boxplot_subplot(df,
                           grouping_column,
                           order,
                           medianprops,
                           ax,
                           ha):
    matches_col, scores_col, label_fontsize, tick_fontsize, text_width, flierprops = init()

    sns.boxplot(
        ax=ax,
        data=df,
        order=order,
        x=grouping_column,
        y='value',
        hue="Number",
        hue_order=[matches_col, np.nan],
        medianprops=medianprops,
        flierprops=flierprops,
        boxprops=dict(alpha=0.5)
    )
   
    rescale_matches_axes(
        df,
        True,
        matches_col,
        label_fontsize,
        tick_fontsize,
        ax
    )

    ax1_2 = ax.twinx()
    sns.boxplot(
        ax=ax1_2,
        data=df,
        order=order,
        x=grouping_column,
        y='value',
        hue='Number',
        hue_order=[np.nan, scores_col],
        medianprops=medianprops,
        flierprops=flierprops,
        boxprops=dict(alpha=0.5)
    )
    rescale_scores_axes(
        True,
        label_fontsize,
        tick_fontsize,
        ax1_2
    )

    count_data1 = df[grouping_column].value_counts()
    make_xticklabels(text_width, ax, count_data1, ha=ha)

    ax.set_xlabel("", fontsize=0)
    ax1_2.set_xlabel("", fontsize=0)

    return ax, ax1_2


def boxplot_comparison(df, order, df_N, order_N, col, colors = ['crimson', 'deepskyblue']):
    matches_col, scores_col, label_fontsize, tick_fontsize, text_width, flierprops = init()
    medianprops = {
        'linewidth': 2.0
    }

    grouping_column = 'composition'	
    plot_width = calc_plot_width(df, grouping_column, 1)
    sns.set_style(style='white')
    sns.set_palette(sns.color_palette(colors))

    scores_df = df[df['Number'] == col]
    scores_df_N = df_N[df_N['Number'] == col]


    fig, ax1 = plt.subplots(figsize=(plot_width, 5))
    sns.boxplot(
        ax=ax1,
        data=scores_df_N,
        order=order_N,
        x=grouping_column,
        y='value',
        hue="Number",
        hue_order=[col, np.nan],
        flierprops=flierprops,
        boxprops=dict(alpha=0.5),
        medianprops= {
            'linewidth': 2.0,
            'color': colors[0]
        }
    )
    ax1.set_xlabel("", fontsize=0)

    if col == 'scores':
        rescale_scores_axes(
            True,
            label_fontsize,
            tick_fontsize,
            ax1
        )
    else:
        rescale_matches_axes(
            scores_df_N,
            True,
            matches_col,
            label_fontsize,
            tick_fontsize,
            ax1
        )

    count_data1 = df_N[grouping_column].value_counts()
    make_xticklabels(text_width, ax1, count_data1, ha='right')
    
    ax2 = ax1.twiny()
    sns.boxplot(
        ax=ax2,
        data=scores_df,
        order=order,
        x=grouping_column,
        y='value',
        hue="Number",
        hue_order=[np.nan, col],
        flierprops=flierprops,
        boxprops=dict(alpha=0.5),
        medianprops= {
            'linewidth': 2.0,
            'color': colors[1]
        }    )
    ax2.set_xlabel("", fontsize=0)

    if col == 'scores':
        rescale_scores_axes(
            True,
            label_fontsize,
            tick_fontsize,
            ax2
        )
    else:
        rescale_matches_axes(
            scores_df,
            True,
            matches_col,
            label_fontsize,
            tick_fontsize,
            ax2
        )

    count_data2 = df[grouping_column].value_counts()
    make_xticklabels(text_width, ax2, count_data2, ha='left')

    handles, labels = ax1.get_legend_handles_labels()
    labels[0] = 'nitrogen'
    labels[1] = 'baseline'
    ax1.legend(handles, labels, loc='upper right', fontsize=14)
    ax2.get_legend().remove()
    return fig

def make_simple_boxplot(df: pd.DataFrame, x: str, y: str, order: List[str] = None, text_width = 21):
    data = df.groupby(x).filter(lambda x: len(x) > 2)
    _, _, label_fontsize, tick_fontsize, _, _ = init()

    plot_width = calc_plot_width(data, x, 1.5)
    fig = plt.figure(figsize=(plot_width, 6))
    ax = sns.boxplot(
    x=x,
    y=y,
    data=data,
    order = order,
    color="deepskyblue",
    medianprops={'color': "blue", 'linewidth': 4.0},
    flierprops={'marker': 'o', 'markersize': 10, 'markerfacecolor': 'none'})

    count_data = data[x].value_counts() * 2

    # # Add the count labels to the x-axis
    make_xticklabels(text_width, ax, count_data, ha="center", rotation=90)
    plt.ylim([0, 1000])
    # plt.ylabel('scores', fontsize = label_fontsize)
    plt.ylabel("", fontsize=0)

    plt.tick_params(axis='y', labelsize=tick_fontsize)
    plt.tick_params(axis='x', labelsize=tick_fontsize)
    plt.xlabel("", fontsize=0)

    return fig