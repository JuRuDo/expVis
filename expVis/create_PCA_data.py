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


import json
import pandas as pd
import argparse
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from pkg_resources import get_distribution
from expVis.read_data import read_json



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
    return pd.DataFrame(mov_data).T


def read_exp(path, t_exp_data, g_exp_data, rel_exp_data):
    exp = read_json(path)
    t_exp_data[exp['name']] = {}
    g_exp_data[exp['name']] = {}
    rel_exp_data[exp['name']] = {}
    for gene in exp['data']:
        for i in range(len(exp['data'][gene]['ids'])):
            t_exp_data[exp['name']][exp['data'][gene]['ids'][i]] = exp['data'][gene]['expression'][i]
            rel_exp_data[exp['name']][exp['data'][gene]['ids'][i]] = exp['data'][gene]['expression_rel'][i]
        g_exp_data[exp['name']][gene] = sum(exp['data'][gene]['expression'])
    return t_exp_data, g_exp_data, rel_exp_data


def exp_main(path, replicates):
    t_exp_data = {}
    g_exp_data = {}
    rel_exp_data = {}
    for replicate in replicates:
        t_exp_data, g_exp_data, rel_exp_data = read_exp(path + replicate + '.json', t_exp_data, g_exp_data,
                                                        rel_exp_data)
    return pd.DataFrame(t_exp_data).T, pd.DataFrame(g_exp_data).T, pd.DataFrame(rel_exp_data).T


def do_pca(pca_DF, replicates):
    """
    This function applies PCA on the given DataFrame containing expression
    data. The overall information content vector and principal components
    for all replicates under all conditions will be returned.

    DataFrame: the DataFrame containing expression value.
    replicates: the list of replicates for each condition.

    """
    # Applying the PCA approach to generate the PC dataframe.
    # features are FPKM/relEXP/EWFD of all genes/transcripts
    features = list(pca_DF.columns)
    x = pca_DF.loc[:, features].values
    x = StandardScaler().fit_transform(x)
    n_repl = len(replicates)
    pca = PCA(n_components=n_repl)
    column = []
    for i in range(n_repl):
        column.append("PC" + str(i + 1))
    principalComponents = pca.fit_transform(x)
    principalDf = pd.DataFrame(data=principalComponents
                               , columns=column, index=replicates)
    # Generate the dictionary for general information content of PCs.
    info_content = pca.explained_variance_ratio_
    info_contentDf = pd.DataFrame(data=info_content
                                  , columns=["information_content"]
                                  , index=column).T
    result_info = info_contentDf.to_dict("index")["information_content"]
    # Add replicate PC coordinates
    repl = principalDf.iloc[0:len(replicates)].to_dict("index")
    return result_info, repl


def read_cmap(path):
    conditions = []
    with open(path, 'r') as infile:
        for line in infile:
            conditions.append(line.rstrip('\n'))
    return conditions


def main(path, infop, cpath, outpath):
    conditions = read_cmap(cpath)
    outdict = {
        'Gene Expression': {'conditions': {}},
        'Transcript Expression': {'conditions': {}},
        'Relative Transcript Expression': {'conditions': {}},
        'Transcript EWFD': {'conditions': {}}
    }
    tmp_data = {}
    info = read_json(infop)
    replicates = list(info['expression_imports']['replicates'].keys())
    t_exp_data, g_exp_data, rel_exp_data = exp_main(path + 'expression/replicates/expression_', replicates)
    outdict['Transcript Expression']['information_content'], tmp_data['Transcript Expression'] = do_pca(t_exp_data,
                                                                                                        replicates)
    outdict['Gene Expression']['information_content'], tmp_data['Gene Expression'] = do_pca(g_exp_data, replicates)
    outdict['Relative Transcript Expression']['information_content'], tmp_data['Relative Transcript Expression'] = \
        do_pca(rel_exp_data, replicates)
    ewfd_data = mov_main(path + 'ewfd/replicates/ewfd_', replicates)
    outdict['Transcript EWFD']['information_content'], tmp_data['Transcript EWFD'] = do_pca(ewfd_data, replicates)
    tmp = ['Transcript Expression', 'Gene Expression', 'Relative Transcript Expression', 'Transcript EWFD']
    for condition in conditions:
        for x in tmp:
            outdict[x]['conditions'][condition] = {}
            for replicate in info['expression_imports']['conditions'][condition]['replicates']:
                outdict[x]['conditions'][condition][replicate] = tmp_data[x][replicate]
    jsonOut = json.dumps(outdict, indent='    ')
    f = open(outpath + '/principle_components.json', 'w')
    f.write(jsonOut)
    f.close()



def get_options():
    version = get_distribution('SpICE').version
    parser = argparse.ArgumentParser(description='You are running SpICE version ' + str(version) + '.',
                                     epilog="This script does the PCA from the results of the SpICE pipeline.")
    required = parser.add_argument_group('required arguments')
    parser.add_argument('--version', action='version', version=str(version))
    required.add_argument("-i", "--resultpath", default=None, type=str, required=True,
                          help="path to the results folder")
    required.add_argument("-c", "--conditions", default=None, type=str, required=True,
                          help="A file containing a list of conditions to be considered for the PCA.")
    required.add_argument("-o", "--outpath", default=None, type=str, required=True,
                          help="path to output directory")
    args = parser.parse_args()
    main(args.resultpath + '/', args.resultpath + '/info.json', args.conditions, args.outpath + '/')


if __name__ == '__main__':
    get_options()
