import textwrap
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import plotly.graph_objs as go
import numpy as np


def make_boxplot(grouped_df: pd.DataFrame, colname: str, legend: str):
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


def plot_histograms_sidebyside(same_query_ref1, same_query_ref2, same_query_ref3, same_query_ref4, same_query_ref5, column_name, xaxis_title='', title=''):
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


def plot_histogram(x, xaxis_title='', title=''):
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
                hide_labels: bool = False):
    matches_col = 'CosineHungarian_0.01_0.0_1.0_matches'
    scores_col = 'CosineHungarian_0.01_0.0_1.0_scores'

    label_fontsize = 20
    tick_fontsize = 13
    text_width = 22

    # Set the style
    colors = ['yellow', 'deepskyblue']
    sns.set_palette(sns.color_palette(colors))
    sns.set_style(style='white')

    # Calculate the number of bars and figure size
    n_bars = df[grouping_column].nunique()
    bar_width = 2
    plot_width = n_bars * bar_width
    fig = plt.figure(figsize=(plot_width, 5))


    ax = sns.boxplot(x=grouping_column, y="value", hue="Number",
                     data=df, hue_order=[matches_col, np.nan],
                     medianprops={'color': 'darkgreen', 'linewidth': 4.0},
                     flierprops={'marker': 'o', 'markersize': 10, 'markerfacecolor': 'none'})  # RUN PLOT
    ax.legend_.remove()

    if normalized_matches:
        ax.set_ylim([0, 100])  # Set y-axis limits
        ax.set_ylabel('ions matching reference', fontsize=label_fontsize)  # Set y-axis label
        ax.set_yticklabels([str(f"{x}%") for x in range(0, 101, 20)])
    else:
        ax.set_ylabel('ion matches', fontsize=label_fontsize)
        ax.set_ylim(0, top=df.loc[df['Number'] == matches_col, 'value'].max())  # Set y-axis limits


    # Set font size of x-axis tick labels
    ax.tick_params(axis='x', labelsize=tick_fontsize)
    # Set font size of y-axis tick labels
    ax.tick_params(axis='y', labelsize=tick_fontsize)

    ax2 = ax.twinx()

    sns.boxplot(ax=ax2, x=grouping_column, y='value', hue='Number',
                data=df, hue_order=[np.nan, scores_col],
                medianprops={'color': 'b', 'linewidth': 4.0},
                flierprops={'marker': 'o', 'markersize': 10, 'markerfacecolor': 'none'})

    ax2.legend_.remove()

    if nist_scale:
        ax2.set_ylim([0, 1000])
    else:
        ax2.set_ylim([0, 1])  # Set y-axis limits

    ax2.set_ylabel('scores', fontsize=label_fontsize)  # Set y-axis label
    # Set font size of y-axis tick labels
    ax2.tick_params(axis='y', labelsize=tick_fontsize)

    if xlabel:
        # Set x-axis label and font size
        ax.set_xlabel(xlabel, fontsize=label_fontsize)
    else:
        ax.set_xlabel("", fontsize=0)
        ax2.set_xlabel("", fontsize=0)
    

    # Create a count for each x-axis label
    count_data = df[grouping_column].value_counts()

    # # Add the count labels to the x-axis
    xlabels = [label.get_text() for label in ax.get_xticklabels()]
    xlabels = ['\n'.join(textwrap.wrap(
        label + ' (' + str((count_data.loc[label] // 2)) + ')', width=text_width)) for label in xlabels]
    ax.set_xticklabels(xlabels, rotation=45, ha='right')

    # Change the legend labels
    if showlegend:
        handles, labels = ax.get_legend_handles_labels()
        labels[0] = ax.get_ylabel()
        labels[1] = ax2.get_ylabel()
        ax.legend(handles, labels, loc='upper right', fontsize=14)
    
    if hide_labels:
        ax.set_ylabel("", fontsize=0)
        ax2.set_ylabel("", fontsize=0)

    return fig


def scatterplot_matplotlib(df):
    fig = plt.figure(figsize=(18, 6))
    scatter = plt.scatter(
        df['CosineHungarian_0.01_0.0_1.0_scores'],
        df['CosineHungarian_0.01_0.0_1.0_matches'],
        # Adjust the size scaling factor as needed
        s=df['FractionQuery'] * 200,
        c=df['FractionReference'] * 100,
        cmap='viridis',  # change the colorscale as needed
        alpha=0.5,
        vmin=0,
        vmax=100
    )
    plt.colorbar(scatter).set_label('Reference Matched %')
    plt.xlabel('scores')
    plt.ylabel('ion matches')

    # Add a legend for the size
    sizes = [1, 50, 100]
    for size in sizes:
        plt.scatter([], [], c='c', alpha=0.5, s=size * 2, label=str(size))
    plt.legend(scatterpoints=1, title='Query Matched %',
               labelspacing=1, loc='upper left')
    return fig
