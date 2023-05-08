#!/bin/env python

#######################################################################
# Copyright (C) 2022 Julian Dosch
#
# This file is part of SpICE.
#
#  SpICE is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  SpICE is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with expNet.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px
import pandas as pd
import math
import numpy


def calc_inter_rmsd(calc):
    rmsds = []
    for x in range(len(calc[0])):
        if x < len(calc[0]) - 1:
            for y in range(x + 1, len(calc[0])):
                rmsds.append(0.0)
                for i in calc:
                    rmsds[-1] += (i[x] - i[y]) ** 2
                rmsds[-1] = math.sqrt(rmsds[-1] / len(calc))
    rmsd = (numpy.mean(rmsds), numpy.mean(rmsds) + numpy.std(rmsds), max(rmsds))
    return rmsd


def prepare_FAS_graph(FAS_data, isoforms, expression, c1, c2, directional, toggle_zero, f_space):
    fas_graph = []
    done = []
    for seed in isoforms:
        c1_mean = expression[(expression['transcriptid'] == seed)
                             & (expression['Condition'] == c1)]['expression'].mean()
        c2_mean = expression[(expression['transcriptid'] == seed)
                             & (expression['Condition'] == c2)]['expression'].mean()
        color = 'black'
        if seed in f_space:
            label = seed + '*'
        else:
            label = seed
        if c1_mean >= c2_mean + 1.0:
            color = 'red'
        if c1_mean <= c2_mean - 1.0:
            color = 'blue'
        fas_graph.append({'data': {
            'id': seed,
            'label': label,
            'color': color,
            'position': {'x': 0, 'y': 0},
        }})
        for query in isoforms:
            if not seed == query and directional:
                if not toggle_zero or (toggle_zero and FAS_data[seed][query] > 0):
                    fas_graph.append({'data': {
                        'source': seed,
                        'target': query,
                        'weight': FAS_data[seed][query]*5+1,
                        'label': f'{float(FAS_data[seed][query]):.2}'
                    }})
            elif not seed == query and not directional:
                if (query, seed) not in done:
                    if not toggle_zero or (toggle_zero and ((FAS_data[seed][query] + FAS_data[query][seed]) / 2) > 0):
                        fas_graph.append({'data': {
                            'source': seed,
                            'target': query,
                            'weight': ((FAS_data[seed][query] + FAS_data[query][seed]) / 2) * 5 + 1,
                            'label': f'{(float(FAS_data[seed][query]) + float(FAS_data[query][seed])) / 2:.2}'
                        }})
                        done.append((query, seed))
    return fas_graph


def mov_figure_polygon(mov_data, gene_id, c1, c2):
    fig = None
    r1 = mov_data[c1]['mean']
    r2 = mov_data[c2]['mean']
    labels = mov_data[c1]['prot_ids']
    i = 0
    while len(r1) < 3:
        r1.append(0.5)
        r2.append(0.5)
        labels.append('D' + str(i))
        i += 1
    r1.append(r1[0])
    r2.append(r2[0])
    labels.append(labels[0])
    fig = make_subplots(rows=1, cols=2, specs=[[{'type': 'polar'}] * 2])
    fig.add_trace(go.Scatterpolar(r=r1,
                                  theta=labels,
                                  fill='toself'),
                  row=1, col=1)
    fig.add_trace(go.Scatterpolar(r=r2,
                                  theta=labels,
                                  fill='toself'),
                  row=1, col=2)
    fig.update_layout(
        polar=dict(
            radialaxis=dict(
                visible=True,
                range=[0, 1]
            )),
        polar2=dict(
            radialaxis=dict(
                visible=True,
                range=[0, 1]
            )),
        showlegend=False,
        title_text=c1 + ' | ' + c2
    )
    return fig


def add_log_fold(r_table, exp_data):
    new_r_table = []
    for entry in r_table:
        new = entry
        new['logFoldChange'] = exp_data[entry['geneid']]['logFoldChange']
        new['minExp'] = exp_data[entry['geneid']]['minExp']
        new_r_table.append(new)
    return new_r_table


def mov_figure_ring(mov_data, gene_id, c1, c2):
    fig = None
    r = {c1: {}, c2: {}}
    for x in ['min', 'l_std', 'mean', 'u_std', 'max']:
        for i in [c1, c2]:
            r[i][x] = mov_data[i][x]
    labels = mov_data[c1]['prot_ids']
    y = 0
    while len(labels) < 3:
        for i in [c1, c2]:
            r[i]['mean'].append(0.5)
            r[i]['min'].append(0.25)
            r[i]['max'].append(0.75)
            r[i]['l_std'].append(0.325)
            r[i]['u_std'].append(0.625)
        labels.append('D' + str(y))
        y += 1
    for x in ['min', 'l_std', 'mean', 'u_std', 'max']:
        for i in [c1, c2]:
            r[i][x].append(r[i][x][0])
    labels.append(labels[0])
    fig = make_subplots(rows=1, cols=2, specs=[[{'type': 'polar'}] * 2])
    tmp = {c1: ['red', 1], c2: ['blue', 2]}
    for i in tmp:
        fig.add_trace(go.Scatterpolar(
            r=r[i]['min'],
            theta=labels,
            marker=dict(color=tmp[i][0], opacity=0.0),
            line=dict(color=None, width=0, shape='spline'),
            fillcolor=tmp[i][0],
            opacity=0.2,
            name='Min ' + i),
            row=1, col=tmp[i][1]
        )
        fig.add_trace(go.Scatterpolar(
            r=r[i]['l_std'],
            theta=labels,
            marker=dict(color=tmp[i][0], opacity=0.0),
            line=dict(color=None, width=0, shape='spline'),
            fillcolor=tmp[i][0],
            fill='tonext',
            opacity=0.4,
            name='Low STD ' + i),
            row=1, col=tmp[i][1]
        )
        fig.add_trace(go.Scatterpolar(
            r=r[i]['mean'],
            theta=labels,
            marker=dict(color=tmp[i][0], opacity=0.0),
            line=dict(color=tmp[i][0], width=3, shape='spline'),
            fillcolor=tmp[i][0],
            fill='tonext',
            opacity=0.4,
            name='Mean ' + i),
            row=1, col=tmp[i][1]
        )
        fig.add_trace(go.Scatterpolar(
            r=r[i]['u_std'],
            theta=labels,
            marker=dict(color=tmp[i][0], opacity=0.0),
            line=dict(color=None, width=0, shape='spline'),
            fillcolor=tmp[i][0],
            fill='tonext',
            opacity=0.2,
            name='High STD ' + i),
            row=1, col=tmp[i][1]
        )
        fig.add_trace(go.Scatterpolar(
            r=r[i]['max'],
            theta=labels,
            marker=dict(color=tmp[i][0], opacity=0.0),
            line=dict(color=None, width=0, shape='spline'),
            fillcolor=tmp[i][0],
            fill='tonext',
            opacity=0.2,
            name='Max ' + i),
            row=1, col=tmp[i][1]
        )
        fig.update_layout(
            polar=dict(
                hole=0.25,
                radialaxis=dict(
                    visible=True,
                    range=[-0.1, 1]
                )),
            polar2=dict(
                hole=0.25,
                radialaxis=dict(
                    visible=True,
                    range=[-0.1, 1]
                )),
            showlegend=False,
            title_text=c1 + ' | ' + c2
        )
    return fig


def organize_fa_data(fa_data, path):
    length, features = fa_data['length'], fa_data['fmap']
    final_fa_data = {}
    if path is None:
        path = list(features.keys())
    for i in path:
        instance = features[i]
        tool = instance[0].split('_')[0]
        if tool not in final_fa_data:
            final_fa_data[tool] = {}
        if instance[0] not in final_fa_data[tool]:
            final_fa_data[tool][instance[0]] = {'x': [], 'y': []}
        final_fa_data[tool][instance[0]]['x'].append(instance[1])
        final_fa_data[tool][instance[0]]['x'].append(instance[2])
        final_fa_data[tool][instance[0]]['x'].append(None)
    lane = 1
    for tool in final_fa_data:
        for feature in final_fa_data[tool]:
            last = [-1]
            entry = 0
            while entry < len(final_fa_data[tool][feature]['x']):
                if final_fa_data[tool][feature]['x'][entry] is None:
                    final_fa_data[tool][feature]['y'].append(None)
                    entry += 1
                else:
                    tmp = True
                    x = 0
                    while tmp:
                        if x >= len(last):
                            final_fa_data[tool][feature]['y'].extend([lane + x, lane + x])
                            last.append(final_fa_data[tool][feature]['x'][entry+1])
                            tmp = False
                        elif final_fa_data[tool][feature]['x'][entry] > last[x]:
                            final_fa_data[tool][feature]['y'].extend([lane + x, lane + x])
                            last[x] = final_fa_data[tool][feature]['x'][entry+1]
                            tmp = False
                        else:
                            x += 1
                    entry += 2
            lane = lane + len(last)
    return final_fa_data, length, lane


def create_fa_plot_input(fa_data, length, isoforms):
    stepsize = 1
    maxlen = 1
    for le in length:
        if not le == None:
            if le > maxlen:
                maxlen = le
    if maxlen >= 100:
        stepsize = int(maxlen / 100)

    tools = []
    for i in range(len(isoforms)):
        if fa_data[isoforms[i]] == None:
            fa_data[isoforms[i]] = {}
        tools.extend(list(fa_data[isoforms[i]].keys()))
    tools = set(tools)
    x, y, labels = [[], []], [[], []], [[], []]
    for tool in tools:
        features = []
        for i in range(len(isoforms)):
            if tool not in fa_data[isoforms[i]]:
                fa_data[isoforms[i]][tool] = {}
            features.extend(list(fa_data[isoforms[i]][tool].keys()))
        features = set(features)
        for feature in features:
            for i in range(len(isoforms)):
                if feature in fa_data[isoforms[i]][tool]:
                    entry = 0
                    while entry < len(fa_data[isoforms[i]][tool][feature]['x']):
                        start = fa_data[isoforms[i]][tool][feature]['x'][entry]
                        stop = fa_data[isoforms[i]][tool][feature]['x'][entry + 1]
                        x[i].append(start)
                        y[i].append(fa_data[isoforms[i]][tool][feature]['y'][entry])
                        labels[i].append(feature)
                        step = start + stepsize
                        while step < stop:
                            x[i].append(step)
                            y[i].append(fa_data[isoforms[i]][tool][feature]['y'][entry])
                            labels[i].append(feature)
                            step += stepsize
                        x[i].extend([stop, None])
                        y[i].extend([fa_data[isoforms[i]][tool][feature]['y'][entry]] * 2)
                        labels[i].extend([feature] * 2)
                        entry += 3
                else:
                    x[i].append(None)
                    y[i].append(None)
                    labels[i].append(feature)
    isoform_labels = []
    for i in range(len(length)):
        if not length[i] == None:
            x[i].extend([0, length[0]])
            y[i].extend([0, 0])
            labels[i].extend(['Protein Length', 'Protein Length'])
            isoform_labels.append(isoforms[i])
        else:
            x[i].extend([0, maxlen, 0, 0, maxlen, 0])
            y[i].extend([-0.25, 1.25, 1.25, -0.25, -0.25, 1.25])
            labels[i].extend(['Non-coding', 'Non-coding', 'Non-coding', 'Non-coding', 'Non-coding', 'Non-coding'])
            isoform_labels.append(isoforms[i] + ' (Non-Coding)')
    return x, y, labels, maxlen, isoform_labels


def create_fa_plot(x, y, labels, maxlen, isoforms, line_size, lane):
    fig = make_subplots(rows=2, cols=1, specs=[[{'type': 'scatter'}], [{'type': 'scatter'}]],
                        subplot_titles=isoforms)
    df1 = pd.DataFrame({'x': x[0], 'y': y[0], 'labels': labels[0], 'fids': labels[0]})
    df2 = pd.DataFrame({'x': x[1], 'y': y[1], 'labels': labels[1], 'fids': labels[1]})
    tmpfig1 = px.line(df1, x='x', y='y', color='labels', custom_data=("fids",))
    tmpfig2 = px.line(df2, x='x', y='y', color='labels', custom_data=("fids",))
    tmp = []
    for trace in range(len(tmpfig1["data"])):
        tmp.append(tmpfig1["data"][trace]['legendgroup'])
        new_trace = tmpfig1["data"][trace]
        new_trace['line']['width'] = line_size
        if new_trace['name'] == 'Protein Length':
            new_trace['line']['color'] = 'black'
            new_trace['line']['width'] = 2
        elif new_trace['name'] == 'Non-coding':
            new_trace['showlegend'] = False
            new_trace['line']['color'] = 'black'
            new_trace['line']['width'] = 2
        fig.add_trace(new_trace, row=1, col=1)
    for trace in range(len(tmpfig2["data"])):
        new_trace = tmpfig2["data"][trace]
        new_trace['line']['width'] = line_size
        if new_trace['legendgroup'] in tmp:
            new_trace['showlegend'] = False
        if new_trace['name'] == 'Protein Length' or new_trace['name'] == 'Non-coding':
            new_trace['showlegend'] = False
            new_trace['line']['color'] = 'black'
            new_trace['line']['width'] = 2
        fig.add_trace(new_trace, row=2, col=1)
    fig.update_traces(connectgaps=False)
    fig.update_yaxes(showticklabels=False)
    fig.update_layout(yaxis_title='', xaxis=dict(range=[0, maxlen]), xaxis2=dict(range=[0, maxlen]),
                      yaxis=dict(range=[-0.5, lane[0]+0.5]), yaxis2=dict(range=[-0.5, lane[1]+0.5]))
    return fig


def prepare_pca_plot_data(pca_data, pc1, pc2, pc3):
    pc1_id = pc1.split(' ')[0]
    pc2_id = pc2.split(' ')[0]
    pc3_id = pc3.split(' ')[0]
    plot_data = {pc1: [], pc2: [], pc3: [], 'Condition': [], 'Replicate': []}
    for condition in pca_data:
        for replicate in pca_data[condition]:
            plot_data[pc1].append(pca_data[condition][replicate][pc1_id])
            plot_data[pc2].append(pca_data[condition][replicate][pc2_id])
            plot_data[pc3].append(pca_data[condition][replicate][pc3_id])
            plot_data['Condition'].append(condition)
            plot_data['Replicate'].append(replicate)
    return plot_data


def create_pca_plot(pc, markersize, pc1, pc2, pc3):
    df = pd.DataFrame(pc)
    fig2 = make_subplots(rows=2, cols=3,
                         specs=[[{'rowspan': 1, 'colspan': 1, 'type': 'scatter'}, None,
                                 {'rowspan': 1, 'colspan': 1, 'type': 'scatter'}],
                                [None, {'rowspan': 1, 'colspan': 1, 'type': 'scatter'}, None],
                                ])

    fig1 = px.scatter_3d(df, x=pc1, y=pc2, z=pc3, color='Condition', hover_data=['Condition', 'Replicate'])
    fig1.update_traces(marker=dict(size=markersize))
    tmp2d1 = px.scatter(df, x=pc1, y=pc2, color='Condition', hover_data=['Condition', 'Replicate'])
    tmp2d2 = px.scatter(df, x=pc1, y=pc3, color='Condition', hover_data=['Condition', 'Replicate'])
    tmp2d3 = px.scatter(df, x=pc2, y=pc3, color='Condition', hover_data=['Condition', 'Replicate'])

    for trace in range(len(tmp2d1["data"])):
        new_trace = tmp2d1["data"][trace]
        new_trace['marker_size'] = markersize*3
        new_trace['showlegend'] = False
        fig2.add_trace(new_trace, row=1, col=1)
    for trace in range(len(tmp2d2["data"])):
        new_trace = tmp2d2["data"][trace]
        new_trace['marker_size'] = markersize*3
        new_trace['showlegend'] = False
        fig2.add_trace(new_trace, row=1, col=3)
    for trace in range(len(tmp2d3["data"])):
        new_trace = tmp2d3["data"][trace]
        new_trace['marker_size'] = markersize*3
        new_trace['showlegend'] = False
        fig2.add_trace(new_trace, row=2, col=2)
    fig2.update_layout(xaxis_title=pc1, yaxis_title=pc2, xaxis2_title=pc1, yaxis2_title=pc3,
                       xaxis3_title=pc2, yaxis3_title=pc3)
    return fig1, fig2
