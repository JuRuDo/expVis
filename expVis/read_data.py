

import json
import numpy


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


def read_library_config(path):
    l_config = {}
    with open(path, 'r') as infile:
        for line in infile.readlines():
            cells = line.rstrip('\n').split('\t')
            if len(cells) > 1:
                l_config[cells[0]] = cells[1]
    return l_config


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


def read_results_exp(path, c1, c2, isoform_data):
    exp1 = read_json(path + '/expression/expression_' + c1 + '.json')
    exp2 = read_json(path + '/expression/expression_' + c2 + '.json')
    exp_data = {}
    for gene in exp1['expression']:
        exp_data[gene] = {'table': [], c1: {}, c2: {}}
        g = []
        if gene in isoform_data:
            g = isoform_data[gene]
        exp_data, c1_exp = read_results_exp_sub(exp1, exp_data, gene, c1, g)
        exp_data, c2_exp = read_results_exp_sub(exp2, exp_data, gene, c2, g)
        log_fold_change = numpy.log2((numpy.mean(c2_exp) + 0.01) / (numpy.mean(c1_exp) + 0.01))
        exp_data[gene]['logFoldChange'] = round(log_fold_change, 4)
        exp_data[gene]['minExp'] = round(min(c1_exp + c2_exp))
        exp_data[gene][c1]['mean'] = round(numpy.mean(c1_exp), 4)
        exp_data[gene][c2]['mean'] = round(numpy.mean(c2_exp), 4)
    return exp_data


def read_results_exp_sub(exp, exp_data, gene, condition, isoform_data):
    exp_l = []
    for replicate in exp['expression'][gene]:
        for prot in exp['expression'][gene][replicate]:
            if prot == 'total':
                exp_data[gene][condition][replicate] = exp['expression'][gene][replicate]['total']
                exp_l.append(exp['expression'][gene][replicate]['total'])
            else:
                if prot in isoform_data:
                    exp_data[gene]['table'].append({
                        'transcriptid': prot,
                        'Condition': condition,
                        'Replicate': replicate,
                        'expression': round(float(exp['expression'][gene][replicate][prot]), 4)
                    })
    return exp_data, exp_l


def read_results_mov(path, c1, c2, fmode, isoform_data):
    mov1 = read_json(path + '/movement/movement_' + c1 + '_' + fmode + '.json')
    mov2 = read_json(path + '/movement/movement_' + c2 + '_' + fmode + '.json')
    mov_data = {}
    for gene in mov1['movement']:
        mov_data[gene] = {'table': [],
                          c1: {'intersample_rmsd': mov1['movement'][gene]['intersample_rmsd_mean']},
                          c2: {'intersample_rmsd': mov2['movement'][gene]['intersample_rmsd_mean']}}
        mov_data = read_results_mov_sub(mov1, mov_data, gene, c1, isoform_data[gene])
        mov_data = read_results_mov_sub(mov2, mov_data, gene, c2, isoform_data[gene])
    return mov_data


def read_results_mov_sub(mov, mov_data, gene, condition, isoform_data):
    keep = []
    for i in range(len(mov['movement'][gene]['prot_ids'])):
        if mov['movement'][gene]['prot_ids'][i] in isoform_data:
            mov_data[gene]['table'].append(
                {'Transcript': mov['movement'][gene]['prot_ids'][i],
                 'Condition': condition,
                 'Min': round(mov['movement'][gene]['min_mov'][i], 4),
                 'Mean': round(mov['movement'][gene]['mean_mov'][i], 4),
                 'Max': round(mov['movement'][gene]['max_mov'][i], 4)
                 })
            keep.append(i)
    mov_data[gene][condition]['prot_ids'] = []
    mov_data[gene][condition]['min'] = []
    mov_data[gene][condition]['l_std'] = []
    mov_data[gene][condition]['mean'] = []
    mov_data[gene][condition]['u_std'] = []
    mov_data[gene][condition]['max'] = []
    for i in keep:
        mov_data[gene][condition]['prot_ids'].append(mov['movement'][gene]['prot_ids'][i])
        mov_data[gene][condition]['min'].append(mov['movement'][gene]['min_mov'][i])
        mov_data[gene][condition]['l_std'].append(mov['movement'][gene]['minus_std_mov'][i])
        mov_data[gene][condition]['mean'].append(mov['movement'][gene]['mean_mov'][i])
        mov_data[gene][condition]['u_std'].append(mov['movement'][gene]['plus_std_mov'][i])
        mov_data[gene][condition]['max'].append(mov['movement'][gene]['max_mov'][i])
    return mov_data

def read_json(path):
    with open(path, 'r') as infile:
        in_dict = json.loads(infile.read())
    return in_dict
