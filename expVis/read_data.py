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
import numpy
import yaml
import support_functions
import math
from yaml.loader import BaseLoader


def read_config_file(path):
    configDict = read_json(path)
    conditions = list(configDict['expression_imports']['conditions'].keys())
    species = configDict['species']
    release = configDict['release']
    replicates = {}
    for condition in conditions:
        replicates[condition] = configDict['expression_imports']['conditions'][condition]['replicates']
    return conditions, species, release, replicates, configDict


def read_library_yaml(path):
    with open(path, 'r') as infile:
        l_config = yaml.load(infile, Loader=BaseLoader)
        return l_config


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


def read_results_exp(path, c1, c2, min_exp, tags, max_tsl):
    # Read in exp data for both conditions
    exp1 = read_json(f'{path}/expression/conditions/expression_{c1}.json')
    exp2 = read_json(f'{path}/expression/conditions/expression_{c2}.json')

    # initialize exp_data and rel_isoforms
    exp_data = {}
    rel_isoforms = {}
    transcript_tags = {}

    # loop through all genes
    for gene in exp1['data']:
        # initialize dictionary for gene
        exp_data[gene] = {'table': [], c1: {}, c2: {}}

        # all the read_results_exp_sub function to process expression data for the current gene in both conditions
        exp_data, c1_exp, isof1, transcript_tags_a, tmp1 = read_results_exp_sub(exp1, exp_data, gene, c1, min_exp,
                                                                                tags, max_tsl)
        exp_data, c2_exp, isof2, transcript_tags_b, tmp2 = read_results_exp_sub(exp2, exp_data, gene, c2, min_exp,
                                                                                tags, max_tsl)
        # add transcript tags
        transcript_tags[gene] = {}
        for t in transcript_tags_a:
            if t not in transcript_tags[gene]:
                transcript_tags[gene][t] = transcript_tags_a[t]
        for t in transcript_tags_b:
            if t not in transcript_tags[gene]:
                transcript_tags[gene][t] = transcript_tags_b[t]

        # For each isoform present in c1 but not in c2, add with expression 0.0 for each replicate, and vice versa
        for i in isof1[1]:
            if i not in isof2[1]:
                tmp2[i] = 0.0
                for r in exp2['replicates']:
                    exp_data[gene]['table'].append({
                        'transcriptid': i,
                        'Condition': c2,
                        'Replicate': r,
                        'expression': 0.0
                    })
        for i in isof2[1]:
            if i not in isof1[1]:
                tmp1[i] = 0.0
                for r in exp1['replicates']:
                    exp_data[gene]['table'].append({
                        'transcriptid': i,
                        'Condition': c1,
                        'Replicate': r,
                        'expression': 0.0
                    })

        # Add the data of expressed and filtered by tags isoforms and calculate relative expression change
        rel_isoforms[gene] = (list(set(isof1[0] + isof2[0])), list(set(isof1[1] + isof2[1])))
        rel_change = 0.0
        for isoform in rel_isoforms[gene][1]:
            rel_change += abs(tmp1[isoform] - tmp2[isoform])
        rel_change = rel_change / 2.0

        #  Calculate the log fold change between the mean expression values of the two conditions
        log_fold_change = numpy.log2((numpy.mean(c2_exp) + 0.01) / (numpy.mean(c1_exp) + 0.01))

        # Add additional data to exp_data
        exp_data[gene]['logFoldChange'] = round(log_fold_change, 4)
        exp_data[gene]['minExp'] = round(min(c1_exp + c2_exp))
        exp_data[gene]['relExpChange'] = round(rel_change, 4)
        exp_data[gene][c1]['mean'] = round(numpy.mean(c1_exp), 4)
        exp_data[gene][c2]['mean'] = round(numpy.mean(c2_exp), 4)

    return exp_data, rel_isoforms, list(exp1['data'].keys()), transcript_tags


def read_results_exp_sub(exp, exp_data, gene, condition, min_exp, tags, max_tsl):
    # Initialize variables and lists
    exp_l = []
    replicates = exp['replicates']
    transcripts = exp['data'][gene]['ids']
    rel_isoforms = [[], []]  # Store isoforms with and without filtering tags
    transcript_tags = {}
    tmp = {}

    # Loop through transcripts
    for x in range(len(transcripts)):
        # Check if maximum expression value is greater than 0
        if max(exp['data'][gene]['expression_all'][x]) > 0:
            rel_isoforms[1].append(transcripts[x])
            tmp[transcripts[x]] = exp['data'][gene]['expression_rel_avg'][x]
            # Loop through replicates and add expression values to dictionary
            for i in range(len(replicates)):
                if not replicates[i] in exp_data[gene][condition]:
                    exp_data[gene][condition][replicates[i]] = 0.0
                exp_data[gene][condition][replicates[i]] += exp['data'][gene]['expression_all'][x][i]
                exp_data[gene]['table'].append({
                    'transcriptid': transcripts[x],
                    'Condition': condition,
                    'Replicate': replicates[i],
                    'expression': round(float(exp['data'][gene]['expression_all'][x][i]), 4)
                })

            # Add data on transcript tags
            transcript_tags[transcripts[x]] = {
                'biotype': exp['data'][gene]['biotypes'][x], 'tags': []}
            for tag in exp['data'][gene]['tags'][x]:
                transcript_tags[transcripts[x]]['tags'].append(tag)
            # Check if maximum expression value is greater than minimum expression cutoff
            if max(exp['data'][gene]['expression_all'][x]) > min_exp:
                # Check if transcript passes filtering tags and transcript support level
                ptags = exp['data'][gene]['tags'][x]
                ptags.append(exp['data'][gene]['biotypes'][x])
                pluscheck = False
                minuscheck = True
                if not tags['+']:
                    pluscheck = True
                else:
                    for tag in tags['+']:
                        if tag in ptags:
                            pluscheck = True
                            break
                if tags['-']:
                    for tag in tags['-']:
                        if tag in ptags:
                            minuscheck = False
                            break
                if pluscheck and minuscheck and exp['data'][gene]['transcript_support_levels'][x] <= max_tsl:
                    rel_isoforms[0].append(transcripts[x])

    # Append expression values to list
    if len(exp_data[gene][condition]) > 0:
        for i in exp_data[gene][condition]:
            exp_l.append(exp_data[gene][condition][i])
    else:
        exp_l.append(0.0)
    return exp_data, exp_l, rel_isoforms, transcript_tags, tmp


def read_results_mov(path, c1, c2, rel_isoforms, exp_data):
    # Read the two ewfd files for the given conditions
    mov1 = read_json(f"{path}/ewfd/conditions/ewfd_{c1}.json")
    mov2 = read_json(f"{path}/ewfd/conditions/ewfd_{c2}.json")

    # initialize mov_data and results
    mov_data = {}
    results = []

    # Loop over all genes in mov1 data (should be synonymous with mov2)
    for gene in mov1["data"]:
        # Create a new dictionary entry to store the data for the current gene
        mov_data[gene] = {
            "table": [],
            c1: {},
            c2: {},
        }

        # Read the results for both conditions and update mov_data
        mov_data, tmp1, max_tsl = read_results_mov_sub(mov1, mov_data, gene, c1, rel_isoforms[gene])
        mov_data, tmp2, max_tsl = read_results_mov_sub(mov2, mov_data, gene, c2, rel_isoforms[gene])

        # Calculate the root-mean-square deviation between the two conditions
        rmsd = 0.0
        for prot in tmp1:
            rmsd += (tmp1[prot] - tmp2[prot])**2
        if len(tmp1) > 0:
            rmsd = round(math.sqrt(rmsd / len(tmp1)), 4)

        # Check if the RMSD between conditions is higher than the std&max of intersample_rmsds for both conditions
        std_check = "No"
        max_check = "No"
        if mov_data[gene][c1]["intersample_rmsd"][1] < rmsd > mov_data[gene][c2]["intersample_rmsd"][1]:
            std_check = "Yes"
            if mov_data[gene][c1]["intersample_rmsd"][2] < rmsd > mov_data[gene][c2]["intersample_rmsd"][2]:
                max_check = "Yes"

        # Add the results for the current gene to the results list
        results.append({
            "geneid": gene,
            "#isoforms": len(rel_isoforms[gene][1]),
            "rmsd": rmsd,
            "std_check": std_check,
            "max_check": max_check,
            "logFoldChange": exp_data[gene]["logFoldChange"],
            "minExp": exp_data[gene]["minExp"],
            "max_tsl": max_tsl,
            "relExpChange": exp_data[gene]["relExpChange"]
        })
    return mov_data, results


def read_results_mov_sub(mov, mov_data, gene, condition, rel_isoforms):
    calc = []
    tmp = {}
    tsl = []

    # Initialize empty lists in mov_data dictionary for the given gene and condition
    mov_data[gene][condition]['prot_ids'] = []
    mov_data[gene][condition]['min'] = []
    mov_data[gene][condition]['l_std'] = []
    mov_data[gene][condition]['mean'] = []
    mov_data[gene][condition]['u_std'] = []
    mov_data[gene][condition]['max'] = []

    # Loop through transcript IDs for the given gene
    for i in range(len(mov['data'][gene]['ids'])):
        # Check if transcript ID is expressed in at least one replicate
        if mov['data'][gene]['ids'][i] in rel_isoforms[1]:
            # Append data to the mov_data table
            mov_data[gene]['table'].append(
                {'Transcript': mov['data'][gene]['ids'][i],
                 'Condition': condition,
                 'Min': round(mov['data'][gene]['ewfd_min'][i], 4),
                 'Mean': round(mov['data'][gene]['avg_ewfd'][i], 4),
                 'Max': round(mov['data'][gene]['ewfd_max'][i], 4)
                 })
            mov_data[gene][condition]['prot_ids'].append(mov['data'][gene]['ids'][i])
            mov_data[gene][condition]['min'].append(mov['data'][gene]['ewfd_min'][i])
            mov_data[gene][condition]['l_std'].append(mov['data'][gene]['avg_ewfd-std'][i])
            mov_data[gene][condition]['mean'].append(mov['data'][gene]['avg_ewfd'][i])
            mov_data[gene][condition]['u_std'].append(mov['data'][gene]['avg_ewfd-std'][i])
            mov_data[gene][condition]['max'].append(mov['data'][gene]['ewfd_max'][i])
            tsl.append(mov['data'][gene]['transcript_support_levels'][i])
            # If transcript is in the list of isoforms for rmsd calculation, append its ewfd_all value to the calc list
            if mov['data'][gene]['ids'][i] in rel_isoforms[0]:
                calc.append(mov['data'][gene]['ewfd_all'][i])
                tmp[mov['data'][gene]['ids'][i]] = mov['data'][gene]['avg_ewfd'][i]

    # Calculate intersample_rmsd value for the calc list
    if len(calc) > 0:
        mov_data[gene][condition]['intersample_rmsd'] = support_functions.calc_inter_rmsd(calc)
    else:
        mov_data[gene][condition]['intersample_rmsd'] = [0.0, 0.0, 0.0]

    if len(tsl) > 0:
        max_tsl = max(tsl)
    else:
        max_tsl = 6

    return mov_data, tmp, max_tsl


def read_json(path):
    with open(path, 'r') as infile:
        in_dict = json.loads(infile.read())
    return in_dict
