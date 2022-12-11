





import pandas as pd
import json


def read_config_file(path):
    configDict = read_json(path)
    conditions = list(configDict['conditions'].keys())
    species = configDict['conditions'][conditions[0]]['species']
    release = configDict['conditions'][conditions[0]]['release']
    fas_modes = configDict['conditions'][conditions[0]]['FAS_modes']
    replicates = {}
    for condition in conditions:
        replicates[condition] = configDict['conditions'][condition]['replicates']
    return conditions, species, release, fas_modes, replicates


def read_results_main(path):
    result_data = []
    genes = []
    isoform_dict = {}
    with open(path, 'r') as infile:
        for line in infile.readlines():
            cells = line.rstrip('\n').split('\t')
            if not (cells[0] == 'gene_id' or cells[0][0] == '!'):
                gene = cells[0]
                isoforms = cells[1].split(';')
                rmsd = cells[11]
                tsl = cells[12]
                std_check = 'Yes'
                if int(cells[8]) or int(cells[10]):
                    std_check = 'No'
                max_check = 'Yes'
                if int(cells[7]) or int(cells[9]):
                    max_check = 'No'
                genes.append(gene)
                isoform_dict[gene] = isoforms
                result_data.append({'geneid': gene, '#isoforms': len(isoforms), 'rmsd': float(rmsd),
                                    'max_tsl': int(tsl), 'std_check': std_check, 'max_check': max_check})
    return result_data, genes, isoform_dict


def read_results_exp(path, c1, c2):
    exp1 = read_json(path + '/expression/expression_' + c1 + '.json')
    exp2 = read_json(path + '/expression/expression_' + c2 + '.json')
    exp_data = {}
    for gene in exp1['expression']:
        exp_data[gene] = {'table': [], c1: {}, c2: {}}
        exp_data = read_results_exp_sub(exp1, exp_data, gene, c1)
        exp_data = read_results_exp_sub(exp2, exp_data, gene, c2)
    return exp_data


def read_results_exp_sub(exp, exp_data, gene, condition):
    for replicate in exp['expression'][gene]:
        for prot in exp['expression'][gene][replicate]:
            if prot == 'total':
                exp_data[gene][condition][replicate] = exp['expression'][gene][replicate]['total']
            else:
                if (float(exp['expression'][gene][replicate]['total']) > 0
                        and exp['expression'][gene][replicate][prot] > 0):
                    expression = round(float(exp['expression'][gene][replicate][prot])
                                       / float(exp['expression'][gene][replicate]['total']), 4)
                    exp_data[gene]['table'].append({
                        'transcriptid': prot,
                        'Condition': condition,
                        'Replicate': replicate,
                        'expression': expression
                    })
    return exp_data


def read_results_mov(path, c1, c2):
    mov1 = read_json(path + 'movement/movement_' + c1 + '.json')
    mov2 = read_json(path + 'movement/movement_' + c2 + '.json')
    mov_data = {}
    for gene in mov1['movement']:
        mov_data[gene] = {'table': [],
                          c1: {'intersample_rmsd': mov1['movement'][gene]['intersample_rmsd_mean']},
                          c2: {'intersample_rmsd': mov2['movement'][gene]['intersample_rmsd_mean']}}
        mov_data = read_results_mov_sub(mov1, mov_data, gene, c1)
        mov_data = read_results_mov_sub(mov2, mov_data, gene, c2)
    return mov_data


def read_results_mov_sub(mov, mov_data, gene, condition):
    for i in range(len(mov['movement'][gene]['protein_ids'])):
        mov_data[gene]['table'].append(
            {'transcript': mov['movement'][gene]['protein_ids'][i],
             'condition': condition,
             'min': mov['movement'][gene]['min_mov'][i],
             'mean': mov['movement'][gene]['mean_mov'][i],
             'max': mov['movement'][gene]['max_mov'][i]
             })
    mov_data[gene][condition]['prot_ids'] = mov['movement'][gene]['prot_ids']
    mov_data[gene][condition]['min'] = mov['movement'][gene]['min_mov']
    mov_data[gene][condition]['l_std'] = mov['movement'][gene]['minus_std_mov']
    mov_data[gene][condition]['mean'] = mov['movement'][gene]['mean_mov']
    mov_data[gene][condition]['u_std'] = mov['movement'][gene]['plus_std_mov']
    mov_data[gene][condition]['max'] = mov['movement'][gene]['max_mov']


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
