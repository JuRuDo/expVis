





import pandas as pd
import json


def read_main_results(path):
    result_data = []
    start = True
    samples = None
    movement = {}
    genes = []
    isoform_dict = {}
    with open(path, 'r') as infile:
        for line in infile.readlines():
            cells = line.rstrip('\n').split('\t')
            if not cells[0] == 'gene_id':
                if start:
                    samples = cells[1].split(';')
                    start = False
                genes.append(cells[0])
                isoforms = cells[2].split(';')
                isoform_dict[cells[0]] = isoforms
                unscaled0, unscaled1 = read_main_results_sub(cells[3].split(';'))
                scaled0, scaled1 = read_main_results_sub(cells[4].split(';'))
                movement[cells[0]] = {}
                for i in range(len(isoforms)):
                    movement[cells[0]][isoforms[i]] = {'unscaled':
                                                           {samples[0]: unscaled0[i], samples[1]: unscaled1[i]},
                                                       'scaled':
                                                           {samples[0]: scaled0[i], samples[1]: scaled1[i]}
                                                       }
                tmp = {'geneid': cells[0], 'isoforms': cells[2], '#isoforms': len(isoforms),
                       'unscaled_rmsd': round(float(cells[5]), 4), 'scaled_rmsd': round(float(cells[6]), 4),
                       'max_tsl': int(cells[7])}
                result_data.append(tmp)
    return result_data, samples, movement, genes, isoform_dict


def read_main_results_sub(movement):
    set0 = movement[0].split(':')
    set1 = movement[1].split(':')
    for i in range(len(set0)):
        set0[i] = round(float(set0[i]), 4)
        set1[i] = round(float(set1[i]), 4)
    return set0, set1


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


def read_json(path):
    with open(path, 'r') as infile:
        in_dict = json.loads(infile.read())
    return in_dict


def read_exp_input(path, sample, exp_data):
    with open(path, 'r') as infile:
        for line in infile.readlines():
            cells = line.rstrip('\n').split('\t')
            if not cells[0] == 'geneID':
                if cells[0] not in exp_data:
                    exp_data[cells[0]] = {'table': []}
                exp_data[cells[0]][sample] = float(cells[2])
                isoforms = cells[1].split(';')
                for isoform in isoforms:
                    values = isoform.split(':')
                    exp_data[cells[0]]['table'].append({
                        'transcriptid': values[0],
                        'Condition': sample,
                        'expression': round(float(values[1]), 4)
                    })
    return exp_data


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
                        'label': f'{FAS_data[seed][query]:.2}'
                    }})
            elif not seed == query and not directional:
                if (query, seed) not in done:
                    if not toggle_zero or (toggle_zero and ((FAS_data[seed][query] + FAS_data[query][seed]) / 2) > 0):
                        fas_graph.append({'data': {
                            'source': seed,
                            'target': query,
                            'weight': ((FAS_data[seed][query] + FAS_data[query][seed]) / 2) * 5 + 1,
                            'label': f'{(FAS_data[seed][query] + FAS_data[query][seed]) / 2:.2}'
                        }})
                        done.append((query, seed))
    return fas_graph
