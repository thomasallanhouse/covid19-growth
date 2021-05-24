# Copyright (c) 2021 Helena Stage

# See LICENCE for licensing information

#

# Analysis of literature citations for growth rate and R estimates

# Figure 4 of main text



#!/usr/bin/env python3

# -*- coding: utf-8 -*-

"""

Created on Wed Feb 10 23:30:34 2021



@author: helena stage

"""

from matplotlib import pyplot as plt

import pandas

import seaborn as sns

import matplotlib.dates as mdates

from matplotlib import markers

import matplotlib.lines as mlines



summary = pandas.read_excel('Literature_data.xlsx')

summary_uom_sep = pandas.read_excel('UoM_estimates.xlsx')

size=summary['Citations']

size_uom_sep=summary_uom_sep['Citations']

alpha = 0.8



review = summary['Preprint']



peer = summary.loc[review=='No']

pre = summary.loc[(review=='Yes') & (summary['Data from Asia']=='Yes')]

pre1 = summary.loc[(review=='Yes') & (summary['Data from Asia']=='No') & (summary['Estimate']=='Yes')]

pre2 = summary.loc[(review=='Yes') & (summary['Data from Asia']=='No') & (summary['Estimate']=='No')]



fig_dims = (6, 4)

fig, ax = plt.subplots(figsize=fig_dims)

g = sns.scatterplot(data=peer, x='Days', y='Double', hue='Data from Asia', palette="colorblind", style='Estimate', size=peer['Citations'], sizes=(20, 1000), alpha=alpha, legend='brief')

ax.set_ylabel("Doubling time (days)", fontsize=8)

ax.xaxis.set_major_locator(mdates.MonthLocator(interval=1))

ax.xaxis.set_major_formatter(mdates.DateFormatter('%d-%m'))

plt.gcf().autofmt_xdate()

ax.tick_params(labelsize=8)

plt.ylim(0, 14)

sns.despine()

g.set(xlabel=None)

g1 = sns.scatterplot(data=summary_uom_sep, x='Days', y='Double', hue='Source', palette=[sns.color_palette("colorblind")[2]], style='Source', markers='s', size=size_uom_sep, sizes=(20, 1000), alpha=alpha, legend='brief', ec=sns.color_palette("colorblind")[2], fc="none", linewidth=1)

g1.set(xlabel=None)

plt.ylim(0, 10)

g2 = sns.scatterplot(data=pre, x='Days', y='Double', hue='Data from Asia', palette=[sns.color_palette("colorblind")[0]], style='Estimate', markers='o', size=pre['Citations'], sizes=(12, 58), alpha=alpha, legend='brief', ec=sns.color_palette("colorblind")[0], fc="none", linewidth=1)

g2.set(xlabel=None)

g3 = sns.scatterplot(data=pre1, x='Days', y='Double', hue='Data from Asia', palette=[sns.color_palette("colorblind")[1]], style='Estimate', markers='o', size=pre1['Citations'], sizes=(12, 22), alpha=alpha, legend='brief', ec=sns.color_palette("colorblind")[1], fc="none", linewidth=1)

g3.set(xlabel=None)

g4 = sns.scatterplot(data=pre2, x='Days', y='Double', hue='Data from Asia', palette=[sns.color_palette("colorblind")[1]], style='Estimate', markers='X', size=pre2['Citations'], sizes=(207, 207), alpha=alpha, legend='brief', ec=sns.color_palette("colorblind")[1], fc="none", linewidth=1)

g4.set(xlabel=None)

h,l = ax.get_legend_handles_labels()





a = mlines.Line2D([], [], color=sns.color_palette("colorblind")[0], marker='o', linestyle='None', markersize=5, label='Asian estimates (peer-reviewed)')

b = mlines.Line2D([], [], color=sns.color_palette("colorblind")[0], marker='o', linestyle='None', markersize=5, label='Asian estimates (not peer-reviewed)', mec=sns.color_palette("colorblind")[0], mfc="none")

c = mlines.Line2D([], [], color=sns.color_palette("colorblind")[1], marker='o', linestyle='None', markersize=5, label='Other estimates (peer-reviewed)')

d = mlines.Line2D([], [], color=sns.color_palette("colorblind")[1], marker='o', linestyle='None', markersize=5, label='Other estimates (not peer-reviewed)', mec=sns.color_palette("colorblind")[1], mfc="none")

e = mlines.Line2D([], [], color=sns.color_palette("colorblind")[1], marker='X', linestyle='None', markersize=5, label='Other assumptions (peer-reviewed)')

f = mlines.Line2D([], [], color=sns.color_palette("colorblind")[1], marker='X', linestyle='None', markersize=5, label='Other assumptions (not peer-reviewed)', mec=sns.color_palette("colorblind")[1], mfc="none")

g = mlines.Line2D([], [], color=sns.color_palette("colorblind")[2], marker='s', linestyle='None', markersize=5, label='Estimates from this work (Figure 2)', mec=sns.color_palette("colorblind")[2], mfc="none")



col_lgd = plt.legend(handles = [a, b, c, d, e, f, g], loc='upper right', 

                     bbox_to_anchor=(1.05, 1.0), fancybox=True, shadow=False, ncol=1, fontsize=7)

plt.savefig("Fig_4a.pdf")



fig_dims = (6, 4)

fig, ax = plt.subplots(figsize=fig_dims)

g = sns.scatterplot(data=peer, x='Days', y='R', hue='Data from Asia', palette="colorblind", style='Estimate', size=peer['Citations'], sizes=(20, 1000), alpha=alpha, legend=False)

ax.set_ylabel("Reproduction number", fontsize=8)

ax.xaxis.set_major_locator(mdates.MonthLocator(interval=1))

ax.xaxis.set_major_formatter(mdates.DateFormatter('%d-%m'))

plt.gcf().autofmt_xdate()

ax.tick_params(labelsize=8)

plt.ylim(0, 7)

sns.despine()

g.set(xlabel=None)

g2 = sns.scatterplot(data=pre, x='Days', y='R', hue='Data from Asia', palette=[sns.color_palette("colorblind")[0]], style='Estimate', markers='o', size=pre['Citations'], sizes=(12, 58), alpha=alpha, legend=False, ec=sns.color_palette("colorblind")[0], fc="none", linewidth=1)

g2.set(xlabel=None)

g3 = sns.scatterplot(data=pre1, x='Days', y='R', hue='Data from Asia', palette=[sns.color_palette("colorblind")[1]], style='Estimate', markers='o', size=pre1['Citations'], sizes=(12, 22), alpha=alpha, legend=False, ec=sns.color_palette("colorblind")[1], fc="none", linewidth=1)

g3.set(xlabel=None)

g4 = sns.scatterplot(data=pre2, x='Days', y='R', hue='Data from Asia', palette=[sns.color_palette("colorblind")[1]], style='Estimate', markers='X', size=pre2['Citations'], sizes=(207, 207), alpha=alpha, legend=False, ec=sns.color_palette("colorblind")[1], fc="none", linewidth=1)

g4.set(xlabel=None)



plt.savefig("Fig_4b.pdf")

