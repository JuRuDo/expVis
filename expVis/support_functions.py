

import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px
import pandas as pd


def prepare_movement(movement, geneid):
    mov_table = []
    figure_data = {'isoforms': []}
    for isoform in movement[geneid]:
        figure_data['isoforms'].append(isoform)
        for sample in movement[geneid][isoform]['scaled']:
            mov_table.append({
                'transcriptid': isoform,
                'Condition': sample,
                'Scaled': movement[geneid][isoform]['scaled'][sample],
                'Unscaled': movement[geneid][isoform]['unscaled'][sample]
            })
            if sample not in figure_data:
                figure_data[sample] = {'scaled': [], 'unscaled': []}
            figure_data[sample]['scaled'].append(movement[geneid][isoform]['scaled'][sample])
            figure_data[sample]['unscaled'].append(movement[geneid][isoform]['unscaled'][sample])
    for sample in figure_data:
        if not sample == 'isoforms':
            figure_data[sample]['scaled'].append(figure_data[sample]['scaled'][0])
            figure_data[sample]['unscaled'].append(figure_data[sample]['unscaled'][0])
        else:
            figure_data[sample].append(figure_data[sample][0])
    return mov_table, figure_data


def prepare_FAS_graph(FAS_data, isoforms, expression, c1, c2, directional, toggle_zero):
    fas_graph = []
    done = []
    for seed in isoforms:
        c1_mean = expression[(expression['transcriptid'] == seed)
                             & (expression['Condition'] == c1)]['expression'].mean()
        c2_mean = expression[(expression['transcriptid'] == seed)
                             & (expression['Condition'] == c2)]['expression'].mean()
        color = 'black'
        if c1_mean >= c2_mean + 0.1:
            color = 'red'
        if c1_mean <= c2_mean - 0.1:
            color = 'blue'
        fas_graph.append({'data': {
            'id': seed,
            'label': seed,
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
    lane = 0
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


def create_fa_plot_input(fa_data, length, isoforms, lanes):
    if len(isoforms) == 2:
        tools = set(list(fa_data[isoforms[0]].keys()) + list(fa_data[isoforms[1]].keys()))
    else:
        tools = list(fa_data[isoforms[0]].keys())
    x, y, labels = [[], []], [[], []], [[], []]
    for tool in tools:
        if len(isoforms) == 2:
            if tool not in fa_data[isoforms[0]]:
                fa_data[isoforms[0]][tool] = {}
            if tool not in fa_data[isoforms[1]]:
                fa_data[isoforms[1]][tool] = {}
            features = set(list(fa_data[isoforms[0]][tool].keys()) + list(fa_data[isoforms[1]][tool].keys()))
        else:
            features = set(list(fa_data[isoforms[0]][tool].keys()))
        for feature in features:
            for i in range(len(isoforms)):
                if feature in fa_data[isoforms[i]][tool]:
                    x[i].extend(fa_data[isoforms[i]][tool][feature]['x'])
                    y[i].extend(fa_data[isoforms[i]][tool][feature]['y'])
                    labels[i].extend([feature]*len(fa_data[isoforms[i]][tool][feature]['x']))
                else:
                    x[i].append(None)
                    y[i].append(None)
                    labels[i].append(feature)
    if len(length) == 2:
        if length[0] > length[1]:
            x[1].extend([length[1]]*2)
            y[1].extend([0, lanes[1]])
            labels[1].extend(['END', 'END'])
        elif length[1] > length[0]:
            x[0].extend([length[0]]*2)
            y[0].extend([0, lanes[0]])
            labels[0].extend(['END', 'END'])
    return x, y, labels


def create_fa_plot(x, y, labels, lengths, isoforms, line_size, lane):
    fig = make_subplots(rows=2, cols=1, specs=[[{'type': 'scatter'}], [{'type': 'scatter'}]],
                        subplot_titles=isoforms)
    df1 = pd.DataFrame({'x': x[0], 'y': y[0], 'labels': labels[0]})
    df2 = pd.DataFrame({'x': x[1], 'y': y[1], 'labels': labels[1]})
    tmpfig1 = px.line(df1, x='x', y='y', color='labels')
    tmpfig2 = px.line(df2, x='x', y='y', color='labels')
    tmp = []
    for trace in range(len(tmpfig1["data"])):
        tmp.append(tmpfig1["data"][trace]['legendgroup'])
        new_trace = tmpfig1["data"][trace]
        new_trace['line']['width'] = line_size
        if new_trace['name'] == 'END':
            new_trace['showlegend'] = False
            new_trace['line']['dash'] = 'dot'
            new_trace['line']['color'] = 'black'
            new_trace['line']['width'] = 2
        fig.add_trace(new_trace, row=1, col=1)
    for trace in range(len(tmpfig2["data"])):
        new_trace = tmpfig2["data"][trace]
        new_trace['line']['width'] = line_size
        if new_trace['legendgroup'] in tmp:
            new_trace['showlegend'] = False
        if new_trace['name'] == 'END':
            new_trace['showlegend'] = False
            new_trace['line']['dash'] = 'dot'
            new_trace['line']['color'] = 'black'
            new_trace['line']['width'] = 2
        fig.add_trace(new_trace, row=2, col=1)
    fig.update_traces(connectgaps=False)
    fig.update_yaxes(showticklabels=False)
    fig.update_layout(yaxis_title='', xaxis=dict(range=[0, max(lengths)]), xaxis2=dict(range=[0, max(lengths)]),
                      yaxis=dict(range=[-0.5, lane[0]+0.5]), yaxis2=dict(range=[-0.5, lane[1]+0.5]))
    return fig
