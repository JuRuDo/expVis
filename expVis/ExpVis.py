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

import io
import os
import dash
import dash_daq as daq
import dash_bootstrap_components as dbc
import pandas as pd
import plotly.express as px
import webbrowser
import dash_cytoscape as cyto
from dash import html, dcc, dash_table
from dash.dependencies import Input, Output, State
from expVis import read_data
from expVis import support_functions
from base64 import b64encode



##


fas_style = [
            {
                'selector': 'node',
                'style': {
                    'label': 'data(label)',
                    'width': 'data(size)',
                    'height': 'data(size)',
                    'background-color': 'data(color)',
                    'background-blacken': 'data(blacken)'
                }
            },
            {
                'selector': 'edge',
                'style': {
                    'label': 'data(label)',
                    'width': 'data(weight)',
                    'curve-style': 'bezier',
                    'arrow-scale': 1,
                    "target-arrow-shape": "triangle",
                }
            }
        ]

cyto.load_extra_layouts()
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.CERULEAN])


########################################

styles = {
    'output': {
        'overflow-y': 'scroll',
        'overflow-wrap': 'break-word',
        'height': 'calc(100% - 25px)',
        'border': 'thin lightgrey solid'
    },
    'tab': {'height': 'calc(98vh - 115px)'}
}

########### Blocks #############

### Library Card

library_card = dbc.Card([
    dbc.CardHeader('Library Information', className="bg-primary text-white"),
    html.Div(dbc.Label('Path')),
    dbc.Row([
        dbc.Col([
            library_input := dbc.Input(
                type='text',
                value='',
                placeholder='Enter Path to library folder here',
            ),
        ], width=10),
        dbc.Col([
            library_load_button := dbc.Button('Load', color="primary", className="me-1"),
        ], width=2),
    ]),
    dcc.Loading(type='default', children=[
        dbc.Row([
            dbc.Label('Species', className='bg-secondary'),
            library_Species := html.P(''),
            dbc.Label('Release', className='bg-secondary'),
            library_Release := html.P(''),
        ]),
    ]),
])

### Result Card

result_card = dbc.Card([
    dbc.CardHeader('Result Information', className="bg-primary text-white"),
    dbc.Label('Path', className='bg-secondary'),
    dbc.Row([
        dbc.Col([
            result_input := dbc.Input(
                type='text',
                value='',
                placeholder='Enter Path to result folder here'
            ),
        ], width=10),
        dbc.Col([
            result_load_button := dbc.Button('Load', color="primary", className="me-1"),
        ], width=2),
    ]),
    dcc.Loading(type='default', children=[
        dbc.Row([
            dbc.Label('Conditions', className='bg-secondary'),
            condition_table := dash_table.DataTable(
                columns=[
                    {'name': 'Condition', 'id': 'id', 'type': 'text'},
                    {'name': 'Replicates', 'id': 'replicates', 'type': 'numeric'},
                ],
                style_data={'textAlign': 'center'},
                data=[],
                page_size=10,
            ),
        ]),
    ]),
])


### PCA Card

pca_card = dbc.Card([
    dbc.CardHeader('PCA Results', className="bg-primary text-white"),
    pca_dropdown := dcc.Dropdown(
        value=None,
        clearable=False,
        options=[]
    ),
    dbc.CardHeader('Components'),
    pc_x_dropdown := dcc.Dropdown(
        value=None,
        clearable=False,
        options=[]
    ),
    pc_y_dropdown := dcc.Dropdown(
        value=None,
        clearable=False,
        options=[]
    ),
    pc_z_dropdown := dcc.Dropdown(
        value=None,
        clearable=False,
        options=[]
    ),
    dbc.CardHeader('Marker Size'),
    pca_point_size := dcc.Input(type="number", min=1, max=30, step=1, value=2),
])


#### Result Loader

main_page = dcc.Tab(label='Main', children=[
    dbc.Row([
        dbc.Col([
            dbc.Row([
                library_card,
                result_card,
            ]),
        ], width=3),
        dbc.Col([
            pca_plot_container := html.Div([
                dbc.Row([
                    pca_plot_3d := dcc.Graph(),
                    ]),
                dbc.Row([
                    pca_plot_2d := dcc.Graph(),
                ]),
            ], hidden=True),
        ], width=6),
        dbc.Col([
            dbc.Row([
                pca_card
            ]),
        ], width=3),
    ])
])


#####
### Gene Selector

filter_options = html.Div([
    dbc.Row([
        dbc.Col([
            dbc.Row([
                dbc.Col([
                    html.Div(dbc.Label("Order"), style={'textAlign': 'center'}),
                    sort2_drop := dcc.Dropdown(['Ascending', 'Descending'], value='Descending')
                ], width=2),
                dbc.Col([
                    html.Div(dbc.Label("Functional Disturbance RMSD"), style={'textAlign': 'center'}),
                    rmsd_slider_0 := dcc.RangeSlider(-0.1, 1, 0.1,
                                                     value=[0, 1],
                                                     allowCross=False,
                                                     tooltip={"placement": "bottom"})
                ], width=4),
                dbc.Col([
                    html.Div(dbc.Label("Replicate Coherence"), style={'textAlign': 'center'}),
                    coherence_drop := dcc.Dropdown(['Coherence [std]', 'Coherence [max]'], value=None)
                ], width=4),
            ]),
            dbc.Row([
                dbc.Col([
                    html.Div(dbc.Label("Fold Change"), style={'textAlign': 'center'}),
                    fold_slider := dcc.RangeSlider(
                        min=0,
                        max=8,
                        step=None,
                        marks={
                            0: '-',
                            1: '-2',
                            2: '-1',
                            3: '-0.5',
                            4: '0',
                            5: '0.5',
                            6: '1',
                            7: '2',
                            8: '+'
                        },
                        value=[0, 8],
                        allowCross=False)
                ], width=4),
                dbc.Col([
                    html.Div(dbc.Label("Min FPKM"), style={'textAlign': 'center'}),
                    min_fpkm := dcc.Input(type="number", min=0, max=9000, step=1, value=1),
                ], width=3),
            ]),
        ], width=9),
        dbc.Col([
            html.Div([
                dbc.Label("Filter by Feature"),
                ao_switch := dcc.RadioItems(
                    ['or', 'and'],
                    'or',
                ),
            ]),
            feature_input := dcc.Dropdown([], multi=True),
        ], width=3),
    ], justify="between", className='mt-3 mb-4')
])


selector_options = html.Div([
    dbc.Row([
        dbc.Col([
            dbc.Label("Show number of rows"),
        ], width=8),
        dbc.Col(html.Div(dbc.Label("Select Gene"), style={'textAlign': 'center'}, className="bg-primary text-white"),
                width=3),
        dbc.Col([
            gene_url := html.A('Ensembl', href='https://www.ensembl.org/', target="_blank")
        ], width=1),
    ]),
    dbc.Row([
        dbc.Col([
            row_drop := dcc.Dropdown(value=10, clearable=False, style={'width': '35%'},
                                     options=[10, 25, 50, 100]),
        ]),
        dbc.Col([
                gene_input := dbc.Input(
                    type='text',
                    list='list-genes',
                    value='',
                ),
        ], width=3),
        dbc.Col([
            gene_select := dbc.Button('Load', color="primary", className='me-1'),
        ], width=1),
    ])
]),


gene_selector = dcc.Tab(label='Gene Selector', children=[
    dbc.Row([
        dbc.Col([
            html.H2("Controller", style={'textAlign': 'center'}),
            html.Div(children=[
                html.Div(dbc.Label('Condition 1'), style={'textAlign': 'center'}, className="bg-primary text-white"),
                c1_drop := dcc.Dropdown(['C1', 'C2', 'C3'], clearable=False, value='C1'),
                html.Div(dbc.Label('Condition 2'), style={'textAlign': 'center'}, className="bg-primary text-white"),
                c2_drop := dcc.Dropdown(['C1', 'C2', 'C3'], clearable=False, value='C2'),
                html.Div(dbc.Label('Filter [RMSD Calculation]'), style={'textAlign': 'center'},
                         className="bg-primary text-white"),
                html.Div(dbc.Label('[+] Tags'), style={'textAlign': 'center'}),
                pos_tags_input := dcc.Dropdown(['protein_coding', 'nonsense_mediated_decay', 'CCDS', 'basic',
                                                'complete', "incomplete", "cds_start_NF", "mRNA_start_NF",
                                                "start_incomplete", "cds_end_NF", "mRNA_end_NF", "end_incomplete",
                                                "NOVEL"],
                                               value=['protein_coding'],
                                               multi=True),
                html.Div(dbc.Label('[-] Tags'), style={'textAlign': 'center'}),
                neg_tags_input := dcc.Dropdown(['protein_coding', 'nonsense_mediated_decay', 'CCDS', 'basic',
                                                'complete', "incomplete", "cds_start_NF", "mRNA_start_NF",
                                                "start_incomplete", "cds_end_NF", "mRNA_end_NF", "end_incomplete",
                                                "NOVEL"],
                                               multi=True),
                html.Div(dbc.Label('Max TSL'), style={'textAlign': 'center'}),
                max_tsl_transcript := dcc.Input(type="number", min=1, max=6, step=1, value=6),
                html.Div(dbc.Label('Min FPKM'), style={'textAlign': 'center'}),
                min_fpkm_transcript := dcc.Input(type="number", min=1, max=9000, step=1, value=1),
                html.Div(dbc.Label('Incomplete Filter'), style={'textAlign': 'center'},
                         className="bg-primary text-white"),
                toggle_incomplete := daq.BooleanSwitch(on=True, color="blue", label='Toggle Incomplete',
                                                       labelPosition="right"),
                incomplete_cutoff := dcc.Input(type="number", min=0.0, max=1.0, step=0.05, value=0.5),
                controller_load := dbc.Button('Load Data', color="primary", className="me-1"),
            ]),
        ], width=2),
        dbc.Col([
            dcc.Loading(type="default", children=[
                dbc.Row(filter_options),
                dbc.Row(selector_options),
                gene_table := dash_table.DataTable(
                    columns=[
                        {'name': 'GeneID', 'id': 'geneid', 'type': 'text'},
                        {'name': 'EWFD RMSD', 'id': 'rmsd', 'type': 'numeric'},
                        {'name': 'Log Fold Change', 'id': 'logFoldChange', 'type': 'numeric'},
                        {'name': 'p-value [FDR-BH]', 'id': 'p-value [FDR-BH]', 'type': 'numeric'},
                        {'name': 'Rel Exp Change', 'id': 'relExpChange', 'type': 'numeric'},
                        {'name': 'Expressed Isoforms', 'id': '#isoforms', 'type': 'numeric'},
                        {'name': 'Replicate Coherence [max]', 'id': 'max_check', 'type': 'text'},
                        {'name': 'Replicate Coherence [std]', 'id': 'std_check', 'type': 'text'},
                    ],
                    data=[],
                    filter_action='native',
                    page_size=10,
                    page_current=0,
                    style_data={
                        'width': '150px', 'minWidth': '150px', 'maxWidth': '150px',
                        'overflow': 'hidden',
                        'textOverflow': 'ellipsis',
                    },
                    export_format="csv",
                ),
                result_data := dcc.Store(data=[{'geneid': 'None', '#isoforms': 0, 'rmsd': 0.0, 'logFoldChange': 0.0,
                                                'std_check': 'No', 'max_check': 'No', 'minExp': 0,
                                                '-log10(p)': 1.0}], id='result_data'),
            ]),
        ], width=10),
    ]),
])

### Volcano Plot
volcano_card = dbc.Card([
    dbc.CardHeader('Options', className="bg-primary text-white"),
    volcano_switch := daq.BooleanSwitch(on=True, color="blue", label='3D', labelPosition="right"),
    dbc.CardHeader('Marker Thresholds', className="bg-primary text-white"),
    dbc.CardHeader('Fold Change'),
    volc_foldC := dcc.Input(type="number", min=0.0, max=5.0, value=2.0),
    dbc.CardHeader('P value'),
    volc_pValue := dcc.Input(type="number", min=0.0, max=1.0, value=0.005),
    dbc.CardHeader('RMSD'),
    volc_rmsd := dcc.Input(type="number", min=0.0, max=1.0,  value=0.5),
    dbc.CardHeader('Marker Size'),
    volcano_point_size := dcc.Input(type="number", min=1, max=30, step=1, value=2),
    volc_html := html.A(
        dbc.Button("Download HTML", color="primary", className="me-1"),
        id="volc_html",
        href="",
        download="Volcano_plot.html"
    ),
    volc_svg := html.A(
        dbc.Button("Download SVG", color="primary", className="me-1"),
        id="volc_svg",
        href="",
        download="Volcano_plot.svg"
    ),
])


volcano_tab = dcc.Tab(label='Volcano Plot', children=[
    dbc.Row([
        dbc.Col([
            volcano_card
        ], width=2),
        dbc.Col([
            volcano_plot := dcc.Graph(style={'height': '85vh'})
        ]),
    ]),
], disabled=True)


### Expression Statistics

expression_stats = dcc.Tab(label="Expression Statistics", children=[
    dbc.Row([
        html.H2(children=['',
                          html.Div(id='exp_header',
                                   style={'display': 'inline', 'textAlign': 'center'})],
                style={'textAlign': 'center'})
    ]),
    dbc.Row([
        dbc.Row([
            dcc.Loading([
                dbc.Row([
                    exp_graph := dcc.Graph(),
                ]),
                dbc.Row([
                    dbc.Col([
                        exp_html := html.A(
                            dbc.Button("Download HTML", color="primary", className="me-1"),
                            id="exp_html",
                            href="",
                            download="Expression_graph.html"
                        ),
                        exp_svg := html.A(
                            dbc.Button("Download SVG", color="primary", className="me-1"),
                            id="exp_svg",
                            href="",
                            download="Expression_graph.svg"
                        ),
                    ], width=6),
                ]),
            ]),
        ]),
        dbc.Row([
            dbc.Row([
                dbc.Col([
                    html.Div(dbc.Label("Sort By"), style={'textAlign': 'center'}),
                    sort_exp_drop := dcc.Dropdown(['Transcript ID', 'Condition', 'Replicate', 'Expression'],
                                                  value='Transcript ID', )
                ], width=2),
                dbc.Col([
                    html.Div(dbc.Label("Order"), style={'textAlign': 'center'}),
                    order_exp_drop := dcc.Dropdown(['Ascending', 'Descending'],
                                                   value='Descending')
                ], width=2),
            ]),
            dbc.Row([
                exp_table := dash_table.DataTable(
                    columns=[
                        {'name': 'TranscriptID', 'id': 'transcriptid', 'type': 'text'},
                        {'name': 'Condition', 'id': 'Condition', 'type': 'text'},
                        {'name': 'Replicate', 'id': 'Replicate', 'type': 'text'},
                        {'name': 'Expression [FPKM]', 'id': 'expression', 'type': 'numeric'}
                    ],
                    data=[],
                    filter_action='native',
                    page_size=10,
                    export_format="csv",
                    style_data={
                        'width': '150px', 'minWidth': '150px', 'maxWidth': '150px',
                        'overflow': 'hidden',
                        'textOverflow': 'ellipsis',
                    }
                ),
            ]),
        ]),
    ]),
], disabled=True)

### Movement Visualisation

mov_vis = dcc.Tab(label="Functional Disturbance", children=[
    dbc.Row([
        html.H2(children=['',
                          html.Div(id='mov_header',
                                   style={'display': 'inline', 'textAlign': 'center'})],
                style={'textAlign': 'center'})
    ]),
    dbc.Row([
        dbc.Col([
            mov_dropdown := dcc.Dropdown(
                value='Mean',
                clearable=False,
                options=['Mean', 'Min/Max']
            ),
            mov_html := html.A(
                dbc.Button("Download HTML", color="primary", className="me-1"),
                id="mov_html",
                href="",
                download="Functional_Disturbance_graph.html"
            ),
            mov_svg := html.A(
                dbc.Button("Download SVG", color="primary", className="me-1"),
                id="mov_svg",
                href="",
                download="Functional_Disturbance_graph.svg"
            ),
        ], width=2),
        dbc.Col([
            mov_graph := dcc.Graph(),
        ]),
    ]),
    dbc.Row([
        dbc.Col([
            dbc.Row([
                dbc.Col([
                    html.Div(dbc.Label("Sort By"), style={'textAlign': 'center'}),
                    sort_mov_drop := dcc.Dropdown(['Transcript ID',
                                                   'Condition',
                                                   'FD (Min)',
                                                   'FD (Mean)',
                                                   'FD (Max)'],
                                                  value='Transcript ID', )
                ]),
                dbc.Col([
                    html.Div(dbc.Label("Order"), style={'textAlign': 'center'}),
                    order_mov_drop := dcc.Dropdown(['Ascending', 'Descending'],
                                                   value='Descending')
                ]),
            ]),
            dbc.Row([
                mov_table := dash_table.DataTable(
                    columns=[
                        {'name': 'Transcript', 'id': 'Transcript', 'type': 'text'},
                        {'name': 'Condition', 'id': 'Condition', 'type': 'text'},
                        {'name': 'FD (Min)', 'id': 'Min', 'type': 'numeric'},
                        {'name': 'FD (Mean)', 'id': 'Mean', 'type': 'numeric'},
                        {'name': 'FD (Max)', 'id': 'Max', 'type': 'numeric'},
                    ],
                    data=[],
                    filter_action='native',
                    page_size=10,
                    export_format="csv",
                    style_data={
                        'width': '150px', 'minWidth': '150px', 'maxWidth': '150px',
                        'overflow': 'hidden',
                        'textOverflow': 'ellipsis',
                    }
                ),
            ]),
        ], width=9),
    ]),
], disabled=True)


### Isoform FAS

##

tap_node = dbc.Card([
    tap_node_header := dbc.CardHeader("Isoform: "),
    html.Div(
        [
            tap_node_url := html.A("Ensembl", href='https://www.ensembl.org/', target="_blank"),
            tap_node_biotype := html.P('Biotype: '),
            tap_node_tags := html.P('Tags: '),
            tap_node_fpkm := html.P('FPKM: '),
            tap_node_ewfd := html.P('EWFD: '),
        ], className="p-4"
    )
], className="m-4")

##

tap_edge = dbc.Card([
    tap_edge_header := dbc.CardHeader("Edge"),
    html.Div(
        [
            tap_edge_seed := html.P('Seed ID: '),
            tap_edge_target := html.P('Target ID: '),
            tap_edge_fas := html.P('Score: ')
        ], className="p-4"
    )
], className="m-4")

###

fa_options = dbc.Card([
    dbc.CardHeader('Score'),
    fa_directional := daq.BooleanSwitch(on=True, color="blue", label='Directional scores',
                                        labelPosition="right"),
    fa_toggle_zero := daq.BooleanSwitch(on=True, color="blue", label='Toggle score 0 edges',
                                        labelPosition="right"),
    dbc.CardHeader('Isoforms'),
    transcript_dropdown := dcc.Dropdown(
        value=['t1', 't2', 't3'],
        multi=True,
        options=['t1', 't2', 't3'],
    ),
    dbc.CardHeader('Label Options'),
    fa_node_labels := daq.BooleanSwitch(on=True, color="blue", label='Node Labels', labelPosition="right"),
    fa_edge_labels := daq.BooleanSwitch(on=True, color="blue", label='Edge Labels', labelPosition="right"),
    dbc.CardHeader('Label Size'),
    label_size := dcc.Input(type="number", min=10, max=50, step=5, value=15),
    fas_svg := dbc.Button("Download", color="primary", className="me-1"),
], className="m-4")

###
iso_fas = dcc.Tab(label="Isoform FAS Graph", children=[
    dbc.Row([
        html.H2(children=['',
                          html.Div(id='FAS_header',
                                   style={'display': 'inline', 'textAlign': 'center'})],
                style={'textAlign': 'center'})
    ]),
    dbc.Row([
        dbc.Col([
            fa_options
        ], width={'size': 3}),
        dbc.Col([
            dbc.Row([
                dbc.Col([
                    fas_figure := cyto.Cytoscape(
                        id='fas_figure',
                        layout={'name': 'circle', 'directed': True},
                        style={'width': '100%', 'height': '40vh'},
                        stylesheet=fas_style,
                        elements={}
                    )
                ], width={'size': 6}),
                dbc.Col([
                    tap_node,
                    tap_edge,
                ], width={'size': 6}),
            ]),
        ], width={'size': 9}),
    ]),
], disabled=True)

#### Feature Architecture

feature_architecture_options = dbc.Card([
    dbc.CardHeader('Isoforms'),
    fa_i1_dropdown := dcc.Dropdown(
        value='T1',
        clearable=False,
        options=['T1', 'T2', 'T3'],
    ),
    fa_i2_dropdown := dcc.Dropdown(
        value=None,
        options=[],
    ),
    fa_linearized := daq.BooleanSwitch(on=True, color='blue', label='Multilayered Architecture',
                                       labelPosition='right', disabled=True),
    dbc.CardHeader('Line Width'),
    line_width := dcc.Input(type="number", min=1, max=30, step=1, value=2),
    fa_html := html.A(
        dbc.Button("Download HTML", color="primary", className="me-1"),
        id="fa_html",
        href="",
        download="Feature_Architecture_graph.html"
    ),
    fa_svg := html.A(
        dbc.Button("Download SVG", color="primary", className="me-1"),
        id="fa_svg",
        href="",
        download="Feature_Architecture_graph.svg"
    ),
], className="m-4")


feature_architecture = dcc.Tab(label="Feature Architecture", children=[
    dbc.Row([
        dbc.Col([
            feature_architecture_options
        ], width=3),
        dbc.Col([
            dcc.Loading(
                fa_plot := dcc.Graph(
                    figure=px.line()
                ),
            ),
        ]),
    ]),
], disabled=True)


#### Exp Analysis

exp_an = dcc.Tab(id='exp_analysis', label='Expression Analysis', children=[
    dbc.Row([
        dcc.Tabs([
            gene_selector,
            volcano_tab,
            expression_stats,
            mov_vis,
            iso_fas,
            feature_architecture,
        ]),
    ]),
], disabled=True)


####### Main #######

app.layout = html.Div([
    dbc.Row(
        dcc.Tabs([
            main_page,
            exp_an,
        ]),
    ),
    exp_store := dcc.Store(data=[], id='exp_store'),
    exp_store2 := dcc.Store(data={}, id='exp_store2'),
    result_details := dcc.Store(data={'path': 'Not selected', 'conditions': [], 'species': 'None',
                                      'version': 'None', 'replicates': [], 'configData': {}},
                                id='result_detail'),
    exp_data := dcc.Store(data={}, id='exp_data'),
    isoform_data := dcc.Store(data={}, id='isoform_data'),
    c_data := dcc.Store(data=['c1', 'c2'], id='c_data'),
    mov_data := dcc.Store(data={'gene1': {'table': [],
                                'c1': {'intersample_rmsd': [0.0, 0.0, 0.0], 'prot_ids': [], 'min': [], 'l_std': [],
                                       'mean': [], 'u_std': [], 'max': []},
                                'c2': {'intersample_rmsd': [0.0, 0.0, 0.0], 'prot_ids': [], 'min': [], 'l_std': [],
                                       'mean': [], 'u_std': [], 'max':[]}}},
                          id='mov_data'),
    mov_gene_data := dcc.Store(data={'table': []}, id='mov_table_data'),
    genes := dcc.Store(data=[], id='genes'),
    fas_data := dcc.Store(data={}, id='fas_data'),
    i_features_data := dcc.Store(data={}, id='i_features_data'),
    library_details := dcc.Store(data={}, id='library_details'),
    fa_map_data := dcc.Store(data={}, id='fa_map_data'),
    fa_ids_data := dcc.Store(data={}, id='fa_ids_data'),
    fa_paths_data := dcc.Store(data={}, id='fa_paths_data'),
    fa_index := dcc.Store(data={}, id='fa_index'),
    library_path := dcc.Store(data='', id='library_path'),
    html.Datalist(id='list-genes',
                  children=[html.Option(value=word) for word in genes]),
    mocks := dcc.Store(data=None, id='mocks'),
    pca_data := dcc.Store(data={}, id='pca_data'),
    transcript_tags := dcc.Store(data={}, id='transcript_tags'),
    current_gene := dcc.Store(data='', id='current_gene'),
    fas_index := dcc.Store(data={}, id='fas_index'),
])




########### Callbacks

### main page
@app.callback(
    Output(result_details, 'data'),
    Output(result_input, 'valid'),
    Output(result_input, 'invalid'),
    Output('exp_analysis', 'disabled'),
    Output(pca_data, 'data'),
    State(result_input, 'value'),
    Input(result_load_button, 'n_clicks'),
)
def load_result_data(path, button):
    config_path = path + '/info.json'
    pca_path = path + '/principle_components.json'
    pca_data = {}
    ctx = dash.callback_context
    if ctx.triggered:
        if os.path.exists(config_path):
            conditions, species, release, replicates, configData = read_data.read_config_file(config_path)
            if os.path.exists(pca_path):
                pca_data = read_data.read_json(pca_path)
            return ({'path': path, 'conditions': conditions, 'species': species, 'version': release,
                    'replicates': replicates, 'configData': configData},
                    True, False, False, pca_data)
        else:
            return ({'path': 'Not selected', 'conditions': [], 'species': 'None', 'version': 'None',
                    'replicates': [], 'configData': {}},
                    False, True, True, pca_data)
    else:
        return ({'path': 'Not selected', 'conditions': [], 'species': 'None', 'version': 'None',
                'replicates': [], 'configData': {}},
                False, True, True, pca_data)


@app.callback(
    Output(library_details, 'data'),
    Output(library_Species, 'children'),
    Output(library_Release, 'children'),
    Output(library_input, 'valid'),
    Output(library_input, 'invalid'),
    Output(fa_ids_data, 'data'),
    Output(fa_index, 'data'),
    Output(fas_index, 'data'),
    Output(library_path, 'data'),
    State(library_input, 'value'),
    Input(library_load_button, 'n_clicks'),
)
def load_library_data(path, button):
    config_path = path + '/info.yaml'
    ctx = dash.callback_context
    if ctx.triggered:
        if os.path.exists(config_path):
            l_config = read_data.read_library_yaml(config_path)
            l_config['path'] = path
            fa_ids = read_data.read_json(path + '/fas_data/annotations.json')['inteprotID']
            fas_index = read_data.read_json(path + '/fas_data/fas_index.json')
            fa_index = read_data.read_json(path + '/fas_data/architectures/index.json')['genes']
            return (l_config, l_config['info']['species'] + ' | ' + l_config['info']['taxon_id'],
                    l_config['info']['release'], True, False, fa_ids, fa_index, fas_index, path)
        else:
            return {}, '', '', False, True, {}, {}, {}, ''
    else:
        return {}, '', '', False, True, {}, {}, {}, ''


@app.callback(
    Output(pca_dropdown, 'options'),
    Output(pca_dropdown, 'value'),
    Input(pca_data, 'data'),
)
def pca_data_action(pca_data):
    r_options = []
    r_value = None
    if pca_data:
        r_options = list(pca_data.keys())
        r_value = r_options[0]
    return r_options, r_value


@app.callback(
    Output(pc_x_dropdown, 'options'),
    Output(pc_x_dropdown, 'value'),
    Input(pca_dropdown, 'value'),
    State(pca_data, 'data')
)
def pca_dropdown_action(pca, pca_data):
    pcs = []
    pc_value = None
    if pca and pca_data:
        for pc in pca_data[pca]['information_content']:
            pc_data = float(pca_data[pca]['information_content'][pc])
            pcs.append(pc + f' {pc_data:.4f}')
        pc_value = pcs[0]
    return pcs, pc_value


@app.callback(
    Output(pc_y_dropdown, 'options'),
    Output(pc_y_dropdown, 'value'),
    Input(pc_x_dropdown, 'value'),
    State(pc_x_dropdown, 'options')
)
def pc_x_dropown_action(x_value, x_data):
    pc_data = []
    pc_value = None
    if x_value and x_data:
        if len(x_data) >= 2:
            x_data.remove(x_value)
            pc_data = x_data
            pc_value = x_data[0]
    return pc_data, pc_value


@app.callback(
    Output(pc_z_dropdown, 'options'),
    Output(pc_z_dropdown, 'value'),
    Input(pc_y_dropdown, 'value'),
    State(pc_y_dropdown, 'options')
)
def pc_y_dropown_action(y_value, y_data):
    pc_data = []
    pc_value = None
    if y_value and y_data:
        if len(y_data) >= 2:
            y_data.remove(y_value)
            pc_data = y_data
            pc_value = y_data[0]
    return pc_data, pc_value


@app.callback(
    Output(pca_plot_container, 'hidden'),
    Output(pca_plot_3d, 'figure'),
    Output(pca_plot_2d, 'figure'),
    Input(pc_z_dropdown, 'value'),
    Input(pca_point_size, 'value'),
    State(pc_x_dropdown, 'value'),
    State(pc_y_dropdown, 'value'),
    State(pca_data, 'data'),
    State(pca_dropdown, 'value')
)
def get_pca_figure(pc_z, point_size, pc_x, pc_y, pca_data, pca_type):
    hidden = True
    if pc_z:
        plot_data = support_functions.prepare_pca_plot_data(pca_data[pca_type]['conditions'], pc_x, pc_y, pc_z)
        hidden = False
    else:
        plot_data = {'pc1': [], 'pc2': [], 'pc3': [], 'Condition': [], 'Replicate': []}
    fig1, fig2 = support_functions.create_pca_plot(plot_data, point_size, pc_x, pc_y, pc_z)
    return hidden, fig1, fig2


### Result Explorer
@app.callback(
    Output(condition_table, 'data'),
    Output(c1_drop, 'options'),
    Output(c1_drop, 'value'),
    Input(result_details, 'data'),
)
def update_result_data(data):
    datatable = []
    c1 = None
    for condition in data['conditions']:
        datatable.append({'id': condition, 'replicates': len(data['replicates'][condition])})
    if len(data['conditions']) > 0:
        c1 = data['conditions'][0]
    return datatable, data['conditions'], c1


@app.callback(
    Output(result_data, 'data'),
    Output(exp_data, 'data'),
    Output(isoform_data, 'data'),
    Output(c_data, 'data'),
    Output(genes, 'data'),
    Output(mov_data, 'data'),
    Output(i_features_data, 'data'),
    Output(feature_input, 'options'),
    Output(transcript_tags, 'data'),
    Input(controller_load, 'n_clicks'),
    State(c1_drop, 'value'),
    State(c2_drop, 'value'),
    State(pos_tags_input, 'value'),
    State(neg_tags_input, 'value'),
    State(max_tsl_transcript, 'value'),
    State(min_fpkm_transcript, 'value'),
    State(result_details, 'data'),
    State(result_data, 'data'),
    State(exp_data, 'data'),
    State(isoform_data, 'data'),
    State(c_data, 'data'),
    State(genes, 'data'),
    State(mov_data, 'data'),
    State(i_features_data, 'data'),
    State(library_details, 'data'),
    State(transcript_tags, 'data'),
    State(toggle_incomplete, 'on'),
    State(incomplete_cutoff, 'value'),
)
def load_result_table(nclicks, c1, c2, pos_tags, neg_tags, max_tsl, min_fpkm, r_details, old_r_data, old_exp_data,
                      old_iso_data, old_c_data, old_gene_data, old_mov_data, old_i_f_data, library_data,
                      old_transcript_tags, toggle_incomplete, incomplete_cutoff):
    # button
    ctx = dash.callback_context
    tags = {'+': pos_tags, '-': neg_tags}

    # initialize with old data
    result_data = old_r_data
    exp_data = old_exp_data
    rel_isoforms = old_iso_data
    c_data = old_c_data
    genes = old_gene_data
    mov_data = old_mov_data
    i_f_data = old_i_f_data
    transcript_tags = old_transcript_tags

    # If button is pressed
    if ctx.triggered:
        # check if both conditions have been selected
        if c1 and c2:
            c_data = [c1, c2]
            # get result path, since directory can be named either c1@c2 or c2@c1
            path = r_details['path']
            fpath = None
            if os.path.isdir(f'{path}/main_comparison/{c1}@{c2}'):
                fpath = f'{path}/main_comparison/{c1}@{c2}/result_{c1}@{c2}_important_features.json'
            elif os.path.isdir(f'{path}/main_comparison/{c2}@{c1}'):
                fpath = f'{path}/main_comparison/{c2}@{c1}/result_{c2}@{c1}_important_features.json'

            # read in exp_data
            exp_data, rel_isoforms, genes, transcript_tags = read_data.read_results_exp(path, c1, c2, min_fpkm, tags,
                                                                                        max_tsl, incomplete_cutoff)
            # read in mov_data and calculate result table
            mov_data, result_data = read_data.read_results_mov(path, c1, c2, rel_isoforms, exp_data, toggle_incomplete)

            # read in important features either from results or library
            if fpath:
                i_f_data = read_data.read_json(fpath)
            elif library_data:
                i_f_data = read_data.read_json(library_data['path'] + '/fas_data/important_features.json')

    return (result_data, exp_data, rel_isoforms, c_data, genes, mov_data, i_f_data, list(i_f_data.keys()),
            transcript_tags)


# GeneSelector Callback
@app.callback(
    Output(c2_drop, 'options'),
    Input(c1_drop, 'value'),
    State(c1_drop, 'options'),
)
def update_con2(value, data):
    if value:
        data.remove(value)
    return data


@app.callback(
    Output(gene_table, 'data'),
    Output(gene_table, 'page_size'),
    Output(volcano_tab, 'disabled'),
    Output(volcano_plot, 'figure'),
    Output(volc_html, 'href'),
    Output(volc_svg, 'href'),
    Input(row_drop, 'value'),
    Input(rmsd_slider_0, 'value'),
    Input(coherence_drop, 'value'),
    Input(fold_slider, 'value'),
    Input(feature_input, 'value'),
    Input(sort2_drop, 'value'),
    Input(ao_switch, 'value'),
    Input(result_data, 'data'),
    Input(min_fpkm, 'value'),
    Input(volcano_switch, 'on'),
    Input(volc_foldC, 'value'),
    Input(volc_pValue, 'value'),
    Input(volc_rmsd, 'value'),
    Input(volcano_point_size, 'value'),
    State(i_features_data, 'data')
)
def update_table_options(row_v, rmsd_v, coherence, fold_v, feature_v, sort2_v, ao_switch, result_data, min_fpkm,
                         volcano_switch, volc_foldC, volc_pValue, volc_rmsd, volcano_point_size, f_dict):
    dff = pd.DataFrame(result_data)
    fig = None
    volc_tab = True
    href_html = ''
    img_b64 = ''
    if result_data:
        dff = dff[(dff['rmsd'] >= rmsd_v[0]) & (dff['rmsd'] <= rmsd_v[1])]
        dff = dff[(dff['minExp'] >= min_fpkm)]
        if coherence == 'Coherence [std]':
            dff = dff[(dff['std_check'] == 'Yes')]
        elif coherence == 'Coherence [max]':
            dff = dff[(dff['max_check'] == 'Yes')]
        fslider_d = {1: -2, 2: -1, 3: -0.5, 4: 0, 5: 0.5, 6: 1, 7: 2}
        if not fold_v[0] == 0:
            dff = dff[dff['logFoldChange'] >= fslider_d[fold_v[0]]]
        if not fold_v[1] == 8:
            dff = dff[dff['logFoldChange'] <= fslider_d[fold_v[1]]]
        if feature_v:
            features = []
            if ao_switch == 'or':
                for feature in feature_v:
                    if feature in f_dict:
                        features.extend(f_dict[feature])
                features = set(features)
            elif ao_switch == 'and':
                for feature in feature_v:
                    if feature in f_dict:
                        features.append(f_dict[feature])
                features = set.intersection(*map(set, features))
            dff = dff[dff['geneid'].isin(features)]

        direct = (sort2_v == 'Ascending')
        dff = dff.sort_values('rmsd', ascending=direct)
        fig = support_functions.volcano_plot(dff, volcano_switch, volc_foldC, volc_pValue, volc_rmsd,
                                             volcano_point_size)
        volc_tab = False

        buffer = io.StringIO()
        fig.write_html(buffer)
        html_bytes = buffer.getvalue().encode()
        encoded_html = b64encode(html_bytes).decode()
        href_html = "data:text/html;base64," + encoded_html
        img_bytes = fig.to_image(format="svg", width=945)
        encoding = b64encode(img_bytes).decode()
        img_b64 = "data:image/svg;base64," + encoding
    return dff.to_dict('records'), row_v, volc_tab, fig, href_html, img_b64


@app.callback(
    Output(gene_input, 'value'),
    Output(gene_input, 'valid', allow_duplicate=True),
    Output(gene_input, 'invalid', allow_duplicate=True),
    Input(gene_table, 'active_cell'),
    State(gene_table, 'derived_virtual_data'),
    State(gene_table, 'page_current'),
    State(gene_table, 'page_size'),
    State(gene_input, 'value'),
    State(gene_input, 'valid'),
    State(gene_input, 'invalid'),
    prevent_initial_call=True
)
def select_gene_table(active_cell, gene_table, current_page, page_size, old_val, old_valid, old_invalid):
    if active_cell:
        if active_cell['column_id'] == 'geneid':
            return gene_table[active_cell['row'] + (current_page*page_size)][active_cell['column_id']], False, False
        else:
            return old_val, old_valid, old_invalid
    else:
        return old_val, old_valid, old_invalid


@app.callback(
    Output(gene_url, 'href'),
    Input(gene_input, 'value'),
    State(genes, 'data'),
)
def create_gene_url(geneid, genes):
    if geneid in genes:
        return 'http://www.ensembl.org/id/' + geneid
    else:
        return 'http://www.ensembl.org/'


@app.callback(
    Output(current_gene, 'data'),
    Output(mov_gene_data, 'data'),
    Output(exp_store, 'data'),
    Output(exp_store2, 'data'),
    Output(transcript_dropdown, 'options'),
    Output(transcript_dropdown, 'value'),
    Output(fa_i1_dropdown, 'options'),
    Output(fa_i1_dropdown, 'value'),
    Output(fa_map_data, 'data'),
    Output(fa_paths_data, 'data'),
    Output(gene_input, 'valid', allow_duplicate=True),
    Output(gene_input, 'invalid', allow_duplicate=True),
    Output(mov_vis, 'disabled'),
    Output(iso_fas, 'disabled'),
    Output(feature_architecture, 'disabled'),
    Output(fas_data, 'data'),
    Input(gene_select, 'n_clicks'),
    State(gene_input, 'value'),
    State(exp_data, 'data'),
    State(isoform_data, 'data'),
    State(c_data, 'data'),
    State(genes, 'data'),
    State(mov_data, 'data'),
    State(fas_index, 'data'),
    State(fa_index, 'data'),
    State(library_path, 'data'),
    prevent_initial_call=True
)
def load_button(n_clicks, geneid, exp_data, isoform_dict, samples, genes, mov_data, fas_index, fa_index, library_path):
    ctx = dash.callback_context
    mock_data = ('Example', {'table': []}, [], {samples[0]: 0.0, samples[1]: 0.0}, ['t1', 't2', 't3'],
                 ['t1', 't2', 't3'], ['t1', 't2', 't3'], 't1', {}, {}, False, True, True, True, True, {})
    if ctx.triggered:
        fa_map, fa_paths, fas_data = {}, {}, {}
        if geneid in genes:
            iso_fas_dis = True
            if geneid in fas_index:
                iso_fas_dis = False
                fas_data = read_data.read_json(library_path + '/fas_data/fas_scores/' + fas_index[geneid])
            fa_i1_start = None
            if geneid in fa_index:
                fileid = str(fa_index[geneid]).rjust(9, '0')
                fa_map = read_data.read_json(library_path + '/fas_data/architectures/' + fileid + '.json')
                fa_paths = read_data.read_json(library_path + '/fas_data/architectures/' + fileid + '_paths.json')
            if isoform_dict[geneid][1]:
                fa_i1_start = isoform_dict[geneid][1][0]
            return (geneid, mov_data[geneid], exp_data[geneid]['table'],
                    {samples[0]: exp_data[geneid][samples[0]], samples[1]: exp_data[geneid][samples[1]]},
                    isoform_dict[geneid][1], isoform_dict[geneid][1], isoform_dict[geneid][1],
                    fa_i1_start, fa_map, fa_paths, True, False, False, iso_fas_dis, iso_fas_dis, fas_data)
        else:
            return mock_data
    else:
        return mock_data

@app.callback(
    Output('FAS_header', 'children'),
    Output('mov_header', 'children'),
    Output('exp_header', 'children'),
    Input(current_gene, 'data'),
)
def update_gene_id(gene_id):
    return gene_id, gene_id, gene_id


# Expression callbacks
@app.callback(
    Output(exp_graph, 'figure'),
    Output(exp_html, 'href'),
    Output(exp_svg, 'href'),
    Output(expression_stats, 'disabled'),
    Input(exp_store, 'data'),
    State(c_data, 'data'),
    State(exp_store2, 'data')
)
def generate_chart(exp_data, conditions, exp_data2):
    fig = None
    href_html = ''
    img_b64 = ''
    disabled = True
    if exp_data:
        dff = pd.DataFrame(exp_data)
        title = conditions[0] + ': ' + str(exp_data2[conditions[0]]['mean']) + ' | ' + conditions[1] + ': '\
                + str(exp_data2[conditions[1]]['mean'])
        fig = px.box(dff, x='transcriptid', y='expression', color='Condition')
        fig.update_layout(title_text=title, yaxis_title='Expression [FPKM]', xaxis_title='Transcript ID',
                          font=dict(size=16))
        buffer = io.StringIO()
        fig.write_html(buffer)
        html_bytes = buffer.getvalue().encode()
        encoded_html = b64encode(html_bytes).decode()
        href_html = "data:text/html;base64," + encoded_html
        img_bytes = fig.to_image(format="svg", width=945)
        encoding = b64encode(img_bytes).decode()
        img_b64 = "data:image/svg;base64," + encoding
        disabled = False
    return fig, href_html, img_b64, disabled


@app.callback(
    Output(exp_table, 'data'),
    Input(sort_exp_drop, 'value'),
    Input(order_exp_drop, 'value'),
    Input(exp_store, 'data')
)
def generate_exp_table(sort_exp, order_exp, exp_data):
    if exp_data:
        sdrop_d = {'Transcript ID': 'transcriptid', 'Condition': 'Condition', 'Expression': 'expression'}
        direct = (order_exp == 'Ascending')
        exp_data = pd.DataFrame(exp_data)
        if sort_exp:
            exp_data = exp_data.sort_values(sdrop_d[sort_exp], ascending=direct)
        exp_data = exp_data.to_dict('records')
    return exp_data


# Movement callbacks

@app.callback(
    Output(mov_graph, 'figure'),
    Output(mov_html, 'href'),
    Output(mov_svg, 'href'),
    Input(mov_dropdown, 'value'),
    Input(mov_gene_data, 'data'),
    State(c_data, 'data'),
)
def generate_mov_figure(figure_type, mov_data, conditions):
    fig = None
    if conditions[0] in mov_data:
        if figure_type == 'Mean':
            fig = support_functions.mov_figure_polygon(mov_data, conditions[0], conditions[1])
        else:
            fig = support_functions.mov_figure_ring(mov_data, conditions[0], conditions[1])
    else:
        fig = support_functions.mov_figure_polygon({conditions[0]: {'mean': [], 'prot_ids': []},
                                                    conditions[1]: {'mean': [], 'prot_ids': []}},
                                                   conditions[0], conditions[1])
    buffer = io.StringIO()
    fig.write_html(buffer)
    html_bytes = buffer.getvalue().encode()
    encoded_html = b64encode(html_bytes).decode()
    href_html = "data:text/html;base64," + encoded_html
    img_bytes = fig.to_image(format="svg", width=945)
    encoding = b64encode(img_bytes).decode()
    img_b64 = "data:image/svg;base64," + encoding
    return fig, href_html, img_b64


@app.callback(
    Output(mov_table, 'data'),
    Input(mov_gene_data, 'data'),
    Input(sort_mov_drop, 'value'),
    Input(order_mov_drop, 'value'),
)
def generate_mov_table(mov_gene_data, sort_mov, order_mov):
    mov_table = []
    if mov_gene_data['table']:
        mov_table = pd.DataFrame(mov_gene_data['table'])
        sdrop_d = {'Transcript ID': 'Transcript', 'Condition': 'Condition', 'EWFD (Min)': 'Min',
                   'EWFD (Mean)': 'Mean', 'EWFD (Max)': 'Max'}
        direct = (order_mov == 'Ascending')
        if sort_mov:
            mov_table = mov_table.sort_values(sdrop_d[sort_mov], ascending=direct)
        mov_table = mov_table.to_dict('records')
    return mov_table


# FAS Page
@app.callback(
    Output(fas_figure, 'elements'),
    Output(fas_figure, 'stylesheet'),
    Input(transcript_dropdown, 'value'),
    Input(fa_node_labels, 'on'),
    Input(fa_edge_labels, 'on'),
    Input(fa_directional, 'on'),
    Input(label_size, 'value'),
    Input(fa_toggle_zero, 'on'),
    State(current_gene, 'data'),
    State(fas_data, 'data'),
    State(c_data, 'data'),
    State(genes, 'data'),
    State(exp_store, 'data'),
    State(isoform_data, 'data'),
)
def fas_cytoscape_figure(transcripts, node_labels, edge_labels, directional, label_size, toggle_zero, geneid, fas_dict,
                         conditions, genes, exp_data, isoform_data):
    if geneid in genes and geneid in fas_dict:
        fas_graph = support_functions.prepare_FAS_graph(fas_dict[geneid], transcripts, pd.DataFrame(exp_data),
                                                        conditions[0], conditions[1], directional, toggle_zero,
                                                        isoform_data[geneid][0])
    else:
        fas_graph = {}
    nodes = {
            'selector': 'node',
            'style': {
                'width': 'data(size)',
                'height': 'data(size)',
                'background-color': 'data(color)',
                'background-blacken': 'data(blacken)',
                'font-size': label_size,
            }
        }
    edges = {
            'selector': 'edge',
            'style': {
                'width': 'data(weight)',
                'font-size': label_size,
            }
        }
    if node_labels:
        nodes['style']['label'] = 'data(label)'
    if edge_labels:
        edges['style']['label'] = 'data(label)'
    if directional:
        edges['style']['curve-style'] = 'bezier'
        edges['style']['arrow-scale'] = 1
        edges['style']['target-arrow-shape'] = 'triangle'
    return fas_graph, [nodes, edges]


@app.callback(
    Output(tap_node_header, 'children'),
    Output(tap_node_url, 'href'),
    Output(tap_node_fpkm, 'children'),
    Output(tap_node_ewfd, 'children'),
    Output(tap_node_biotype, 'children'),
    Output(tap_node_tags, 'children'),
    Input('fas_figure', 'tapNodeData'),
    State(exp_store2, 'data'),
    State(exp_store, 'data'),
    State(c_data, 'data'),
    State(mov_gene_data, 'data'),
    State(transcript_tags, 'data'),
    State(current_gene, 'data'),
)
def display_tap_node_data(data, exp_gene, exp_isoforms, conditions, mov_gene_data, transcript_tags, gene_id):
    if data and gene_id:
        exp_isoforms = pd.DataFrame(exp_isoforms)
        exp_isoforms_c1 = exp_isoforms[(exp_isoforms['transcriptid'] == data['id'])
                                       & (exp_isoforms['Condition'] == conditions[0])]
        exp_isoforms_c2 = exp_isoforms[(exp_isoforms['transcriptid'] == data['id'])
                                       & (exp_isoforms['Condition'] == conditions[1])]
        isomean = (exp_isoforms_c1['expression'].mean(), exp_isoforms_c2['expression'].mean())
        mov_table = pd.DataFrame(mov_gene_data['table'])
        mov_mean_c1 = mov_table[(mov_table['Transcript'] == data['id'])
                                & (mov_table['Condition'] == conditions[0])]['Mean'].mean()
        mov_mean_c2 = mov_table[(mov_table['Transcript'] == data['id'])
                                & (mov_table['Condition'] == conditions[1])]['Mean'].mean()
        biotype = 'NA'
        tags = 'NA'
        if gene_id in transcript_tags:
            if data['id'] in transcript_tags[gene_id]:
                biotype = transcript_tags[gene_id][data["id"]]["biotype"]
                tags = ", ".join(transcript_tags[gene_id][data["id"]]["tags"])
        return ("Isoform: " + data['label'], 'http://www.ensembl.org/id/' + data['id'],
                f'FPKM: {isomean[0]:.4f} (of {exp_gene[conditions[0]]["mean"]:.4f})  /  {isomean[1]:.4f}'
                f' (of {exp_gene[conditions[1]]["mean"]:.4f})', f'EWFD: {mov_mean_c1:.4f} / {mov_mean_c2:.4f}',
                f'Biotype: {biotype}',
                f'Tags: {tags}')
    else:
        return 'Isoform: NA', 'http://www.ensembl.org/', 'FPKM: NA', 'EWFD: NA', 'Biotype: NA', 'Tags: NA'


@app.callback(
    Output(tap_edge_seed, 'children'),
    Output(tap_edge_target, 'children'),
    Output(tap_edge_fas, 'children'),
    Input('fas_figure', 'tapEdgeData'),
)
def display_tap_edge_data(data):
    if data:
        return 'Seed: ' + data['source'], 'Target: ' + data['target'], 'Score: ' + data['label']
    else:
        return 'Seed: NA', 'Target: NA', 'Score: NA'


@app.callback(
    Output(fas_figure, "generateImage"),
    Input(fas_svg, "n_clicks"),
)
def get_image(get_svg_clicks):
    ctx = dash.callback_context
    action = 'store'
    ftype = 'png'
    if ctx.triggered:
        action = "download"
        ftype = 'svg'
    return {
        'type': ftype,
        'action': action
        }



@app.callback(
    Output(fa_i2_dropdown, 'options'),
    Output(fa_i2_dropdown, 'value'),
    State(fa_i1_dropdown, 'options'),
    Input(fa_i1_dropdown, 'value'),
)
def fill_fa_i2_dropdown(isoforms, i1_value):
    if i1_value:
        isoforms.remove(i1_value)
        i2_value = None
        if isoforms:
            i2_value = isoforms[0]
    return isoforms, i2_value


@app.callback(
    Output(fa_linearized, 'disabled'),
    Output(fa_linearized, 'on'),
    Input(fa_i1_dropdown, 'value'),
    Input(fa_i2_dropdown, 'value'),
    Input(fa_paths_data, 'data'),
    State(fa_linearized, 'on')
)
def enable_linearized(i1, i2, paths, on):
    if i1 and i2 and paths:
        return False, on
    else:
        return True, on


@app.callback(
    Output(fa_plot, 'figure'),
    Output(fa_html, 'href'),
    Output(fa_svg, 'href'),
    Input(fa_i1_dropdown, 'value'),
    Input(fa_i2_dropdown, 'value'),
    Input(fa_linearized, 'on'),
    Input(line_width, 'value'),
    State(fa_map_data, 'data'),
    State(gene_input, 'value'),
    State(fa_paths_data, 'data')
)
def generate_fa_plot(isoform1, isoform2, linearized, line_width, fa_data, gid, paths):
    href = ''
    if gid in fa_data and isoform1:
        lin = [None, None]
        if (not linearized) and paths:
            lin = paths[gid][isoform1 + '@' + isoform2]
        data = {}
        length = []
        lanes = [1, 1]
        if isoform1 in fa_data[gid]:
            data[isoform1], tmp, tmp2 = support_functions.organize_fa_data(fa_data[gid][isoform1], lin[0])
            length.append(tmp)
            lanes[0] = tmp2
        else:
            data[isoform1] = None
            length.append(None)
        isoforms = [isoform1]

        if isoform2:
            if isoform2 in fa_data[gid]:
                data[isoform2], tmp, tmp2 = support_functions.organize_fa_data(fa_data[gid][isoform2], lin[1])
                length.append(tmp)
                lanes[1] = tmp2
            else:
                data[isoform2] = None
                length.append(None)
            isoforms.append(isoform2)
        x, y, labels, maxlen, isoform_labels = support_functions.create_fa_plot_input(data, length, isoforms)
        fig = support_functions.create_fa_plot(x, y, labels, maxlen, isoform_labels, line_width, lanes)
        buffer = io.StringIO()
        fig.write_html(buffer)
        html_bytes = buffer.getvalue().encode()
        encoded_html = b64encode(html_bytes).decode()
        href_html = "data:text/html;base64," + encoded_html
        img_bytes = fig.to_image(format="svg", width=945)
        encoding = b64encode(img_bytes).decode()
        img_b64 = "data:image/svg;base64," + encoding
        return fig, href_html, img_b64
    else:
        return px.line(x=[1], y=[1]), href, href


@app.callback(
    Output(mocks, 'data'),
    Input(fa_plot, 'clickData'),
    State(fa_ids_data, 'data'),
)
def fa_open_url(clickData, fa_ids):
    if clickData:
        fname = clickData['points'][0]['customdata'][0]
        if fname in fa_ids:
            url = 'https://www.ebi.ac.uk/interpro/entry/' + fname.split('_')[0] + '/' \
                  + '.'.join(fa_ids[fname].split('.')[0:-1])
            webbrowser.open_new_tab(url)
    return None


def main():
    app.run_server(debug=True)


if __name__ == '__main__':
    main()

