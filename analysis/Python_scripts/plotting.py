import pandas as pd
from matplotlib import pyplot as plt
import plotly.graph_objs as go

def make_boxplot(grouped_df: pd.DataFrame, colname: str, legend: str):
    #Create the plot with a width of 10 inches
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
    #ax.set_ylim([0, 1])
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
    trace1 = go.Histogram(x=x1, nbinsx=n_bins, name='All peaks ({})'.format(len(x1)), xbins=dict(start=0, end=1))
    trace2 = go.Histogram(x=x2, nbinsx=n_bins, name='Top 3 peaks ({})'.format(len(x2)), xbins=dict(start=0, end=1))
    trace3 = go.Histogram(x=x3, nbinsx=n_bins, name='Top 5 peaks ({})'.format(len(x3)), xbins=dict(start=0, end=1))
    trace4 = go.Histogram(x=x4, nbinsx=n_bins, name='Top 10 peaks ({})'.format(len(x4)), xbins=dict(start=0, end=1))
    trace5 = go.Histogram(x=x5, nbinsx=n_bins, name='Top 20 peaks ({})'.format(len(x5)), xbins=dict(start=0, end=1))

    # Create layout
    layout = go.Layout(title=title,
                       xaxis=dict(title=xaxis_title, range=[0, 1], dtick=0.1, tickfont=dict(size=15)),
                       yaxis=dict(title='Frequency', range=[0, 140], tickfont=dict(size=15)),
                       legend=dict(x=0.81, y=1.0),
                       font=dict(size=17))

    # Create figure
    fig = go.Figure(data=[trace1, trace2, trace3, trace4, trace5], layout=layout)

    # Display the plot
    fig.show()

def plot_histogram(x, xaxis_title='', title=''):
    # Define number of bins
    n_bins = 20
  

    # Create histogram traces
    trace1 = go.Histogram(x=x, nbinsx=n_bins, name='All peaks ({})'.format(len(x)), xbins=dict(start=0, end=1))
  

    # Create layout
    layout = go.Layout(title=title,
                       xaxis=dict(title=xaxis_title, range=[0, 1], dtick=0.1, tickfont=dict(size=15)),
                       yaxis=dict(title='Frequency', range=[0, 140], tickfont=dict(size=15)),
                       legend=dict(x=0.81, y=1.0),
                       font=dict(size=17))

    # Create figure
    fig = go.Figure(data=[trace1], layout=layout)

    # Display the plot
    fig.show()