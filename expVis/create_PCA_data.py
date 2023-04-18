

from read_data import read_json
import pandas

def read_mov(path, mov_data):
    mov = read_json(path)
    mov_data[mov['name']] = {}
    for gene in mov['data']:
        for i in range(len(mov['data'][gene]['ids'])):
            mov_data[mov['name']][mov['data'][gene]['ids'][i]] = mov['data'][gene]['ewfd_rel_expr'][i]
    return mov_data


def mov_main(path, replicates):
    mov_data = {}
    for replicate in replicates:
        mov_data = read_mov(path + replicate + '.json', mov_data)
    return pandas.DataFrame(mov_data).T


def read_exp(path, t_exp_data, g_exp_data, rel_exp_data):
    exp = read_json(path)
    t_exp_data[exp['name']] = {}
    g_exp_data[exp['name']] = {}
    rel_exp_data[exp['name']] = {}
    for gene in exp['data']:
        for i in range(len(exp['data'][gene]['ids'])):
            t_exp_data[exp['name']][exp['data'][gene]['ids'][i]] = exp['data'][gene]['expression'][i]
            rel_exp_data[exp['name']][exp['data'][gene]['ids'][i]] = exp['data'][gene]['expression_rel'][i]
        g_exp_data[gene] = sum(exp['data'][gene]['expression'])
    return t_exp_data, g_exp_data


def exp_main(path, replicates):
    t_exp_data = {}
    g_exp_data = {}
    rel_exp_data = {}
    for replicate in replicates:
        exp_data = read_mov(path + replicate + '.json', t_exp_data, g_exp_data, rel_exp_data)
    return pandas.DataFrame(t_exp_data).T, pandas.DataFrame(g_exp_data).T, pandas.DataFrame(rel_exp_data).T


def calc_pca_data(data):
    
    pass
