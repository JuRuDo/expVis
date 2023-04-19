

import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from read_data import read_json

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
        g_exp_data[gene] = sum(exp['data'][gene]['expression'])
    return t_exp_data, g_exp_data


def exp_main(path, replicates):
    t_exp_data = {}
    g_exp_data = {}
    rel_exp_data = {}
    for replicate in replicates:
        exp_data = read_mov(path + replicate + '.json', t_exp_data, g_exp_data, rel_exp_data)
    return pd.DataFrame(t_exp_data).T, pd.DataFrame(g_exp_data).T, pd.DataFrame(rel_exp_data).T


def PCA_tool_expr(pca_DF, replicates, conditions):
    """
    This function applies PCA on the given DataFrame containing expression
    data. The overall information content vector and principal components
    for all replicates under all conditions will be returned.

    DataFrame: the DataFrame containing expression value.
    replicates: the list of replicates for each condition.
    conditions: the list of involved conditions.

    """
    # Applying the PCA approach to generate the PC dataframe.
    tdf = pca_DF.T
    ids = tdf.columns
    x = tdf.loc[:, ids].values
    y = pd.DataFrame(tdf.index.values, columns=["repl"])
    x = StandardScaler().fit_transform(x)
    # Number of all replicates = number of DF columns
    # column is the ids.
    n_repl = len(pd.columns)
    pca = PCA(n_components=n_repl)
    column = []
    for i in range(n_repl):
        column.append("PC" + str(i + 1))
    principalComponents = pca.fit_transform(x)
    # Columns are the PCs and index replicates.
    principalDf = pd.DataFrame(data=principalComponents
                               , columns=column, index=y["repl"])

    # Iterativly add condition and its corresponded replicates with CP to
    # a dictionary.
    temp_pos = 0
    cond_dict = {}
    for i in range(len(conditions)):
        repl_dict = {}
        repl = principalDf.iloc[temp_pos:temp_pos + len(replicates[i])]. \
            to_dict("index")
        repl_dict["replicates"] = repl
        cond_dict[conditions[i]] = repl_dict
        temp_pos += len(replicates[i])

    # Generate the dictionary for general information content.
    info_content = pca.explained_variance_ratio_
    info_contentDf = pd.DataFrame(data=info_content
                                  , columns=["information_content"]
                                  , index=column).T
    result_info = info_contentDf.to_dict("index")["information_content"]

    return result_info, cond_dict


def calc_pca_data(data):
    
    pass
