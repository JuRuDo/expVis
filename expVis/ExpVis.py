import os
import dash
import dash_cytoscape as cyto
from dash import html
from dash import dcc
from dash import dash_table
from dash.dependencies import Input, Output, State
import dash_daq as daq
import dash_bootstrap_components as dbc
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import read_data


######################

genedata, samples, movement, genes, isoform_dict = read_data.read_main_results('/home/jd/Data/expVis_data/polygonFAS_all_H1xall_embryo_ED_sorted.tsv')
f_dict = read_data.read_json('/home/jd/Data/expVis_data/isoforms_important_features.json')
exp_data = read_data.read_exp_input('/home/jd/Data/expVis_data/mov_exp_data/polygonFAS_H1.tsv', samples[0], {})
exp_data = read_data.read_exp_input('/home/jd/Data/expVis_data/mov_exp_data/polygonFAS_embryo_ED.tsv', samples[1], exp_data)
fas_dict = read_data.read_json('/home/jd/Data/expVis_data/FAS_scores/distance_master.json')

## Mock Data ##
condition_mock = [{'id': 'c1', 'replicates': 3}, {'id': 'c2', 'replicates': 2}]
fas_mode_mock = [{'id': 'all'}, {'id': 'TMHMM+SignalP'}]
exp_mock = [
    {'transcriptid': 'example1', 'Condition': samples[0], 'expression': 0.5},
    {'transcriptid': 'example2', 'Condition': samples[0], 'expression': 0.5},
    {'transcriptid': 'example3', 'Condition': samples[0], 'expression': 0.5},
    {'transcriptid': 'example1', 'Condition': samples[1], 'expression': 0.2},
    {'transcriptid': 'example2', 'Condition': samples[1], 'expression': 0.7},
    {'transcriptid': 'example3', 'Condition': samples[1], 'expression': 0.5},
]
m_mock0 = [
    {'transcriptid': 'example1', 'Condition': samples[0], 'Unscaled': 0.5, 'Scaled': 0.5},
    {'transcriptid': 'example2', 'Condition': samples[0], 'Unscaled': 0.5, 'Scaled': 0.5},
    {'transcriptid': 'example3', 'Condition': samples[0], 'Unscaled': 0.5, 'Scaled': 0.5},
    {'transcriptid': 'example1', 'Condition': samples[1], 'Unscaled': 0.5, 'Scaled': 0.5},
    {'transcriptid': 'example2', 'Condition': samples[1], 'Unscaled': 0.5, 'Scaled': 0.5},
    {'transcriptid': 'example3', 'Condition': samples[1], 'Unscaled': 0.5, 'Scaled': 0.5},
]
m_mock1 = {
    'isoforms': ['example1', 'example2', 'example3', 'example1'],
    samples[0]: {'scaled': [0.5, 0.5, 0.5, 0.5], 'unscaled': [0.5, 0.5, 0.5, 0.5]},
    samples[1]: {'scaled': [0.5, 0.5, 0.5, 0.5], 'unscaled': [0.5, 0.5, 0.5, 0.5]},
}
fas_mock = [
    {'data': {'id': 't1',
              'label': 't1',
              'position': {'x': 0, 'y': 0}
              }
     },
    {'data': {'id': 't2',
              'label': 't2',
              'position': {'x': 20, 'y': 0}
              }
     },
    {'data': {'id': 't3',
              'label': 't3',
              'position': {'x': 0, 'y': 20}
              }
     },
    {'data': {'source': 't1',
              'target': 't2',
              'weight': 0.0*5+1,
              'label': '0.0',
              }
     },
    {'data': {'source': 't2',
              'target': 't1',
              'weight': 0.1*5+1,
              'label': '0.1',
              }
     },
    {'data': {'source': 't2',
              'target': 't3',
              'weight': 0.1*5+1,
              'label': '0.1',
              }
     },
    {'data': {'source': 't3',
              'target': 't2',
              'weight': 0.0*5+1,
              'label': '0.0',
              }
     },
    {'data': {'source': 't1',
              'target': 't3',
              'weight': 0.7*5+1,
              'label': '0.7',
              }
     },
    {'data': {'source': 't3',
              'target': 't1',
              'weight': 0.8*5+1,
              'label': '0.8',
              }
     },
]


fa_mock = {
    'isoform1': dict(
        y=[1, 1, None, 1.5, 1.5, None, 1, 1, None, 1, 1, 2, 2, None, 2, 2],
        position=[3, 5, None, 4, 7, None, 9, 12, None, 14, 17, 5, 13, None, 15, 20],
        color=['Feature1', 'Feature1', 'Feature1', 'Feature1', 'Feature1', 'Feature1', 'Feature1', 'Feature1',
               'Feature1', 'Feature1', 'Feature1', 'Feature2', 'Feature2', 'Feature2', 'Feature2', 'Feature2'],
        ),
    'isoform2': dict(
        y=[1, 1, None, 1.5, 1.5, None, 1, 1, None, 1, 1, 2, 2],
        position=[3, 5, None, 4, 7, None, 9, 12, None, 14, 17, 5, 13],
        color=['Feature1', 'Feature1', 'Feature1', 'Feature1', 'Feature1', 'Feature1', 'Feature1', 'Feature1',
               'Feature1', 'Feature1', 'Feature1', 'Feature2', 'Feature2'],
        )
}

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

features = list(f_dict.keys())

conditions = samples + ['both']



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
            dbc.Button('Load', color="primary", className="me-1"),
        ], width=2),
    ]),
    dbc.Label('Species', className='bg-secondary'),
    dbc.Label('Release', className='bg-secondary'),
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
    dbc.Label('Conditions', className='bg-secondary'),
    condition_table := dash_table.DataTable(
        columns=[
            {'name': 'Condition', 'id': 'id', 'type': 'text'},
            {'name': 'Replicates', 'id': 'replicates', 'type': 'numeric'},
        ],
        style_data={'textAlign': 'center'},
        data=condition_mock,
        page_size=10,
    ),
    dbc.Label('FAS Modes', className='bg-secondary'),
    fas_modes_table := dash_table.DataTable(
        columns=[
            {'name': 'Mode', 'id': 'id', 'type': 'text'},
        ],
        style_data={'textAlign': 'center'},
        data=fas_mode_mock,
        page_size=10,
    ),
])

#### Result Loader

main_page = dcc.Tab(label='Main', children=[
    dbc.Col([
        dbc.Row([
            library_card,
            result_card,
        ]),
    ], width=3),
])


#####
#### Library Explorer
lib_expl = dcc.Tab(label='Library Explorer', children=[
    dcc.Tabs([
        dcc.Tab(label='Statistics', children=[]),
        dcc.Tab(label='FAS graph', children=[
            dbc.Row([
                dbc.Col([
                    fas_library_dropdown := dcc.Dropdown(
                        value=conditions[-1],
                        clearable=False,
                        options=[
                            {'label': name, 'value': name}
                            for name in ['Full FAS', 'tmhmm+signalp', 'lcr']
                        ]
                    ),
                    fas_library_png := dbc.Button("Download .png", color="primary", className="me-1"),
                    fas_library_svg := dbc.Button("Download .svg", color="primary", className="me-1"),
                ], width={'size': 2}),
                dbc.Col([
                    fas_library_figure := cyto.Cytoscape(
                        id='fas_library_figure',
                        layout={'name': 'circle', 'directed': True},
                        style={'width': '100%', 'height': '40vh'},
                        stylesheet=fas_style,
                        elements=fas_mock
                    )
                ], width={'size': 6}),
            ]),
        ])
    ])
], disabled=True)

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
                    html.Div(dbc.Label("Movement RMSD"), style={'textAlign': 'center'}),
                    rmsd_slider_0 := dcc.RangeSlider(0, 1, 0.1,
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
                    html.Div(dbc.Label("Transcript Support Level"), style={'textAlign': 'center'}),
                    tsl_slider := dcc.RangeSlider(
                        1, 5, 1,
                        value=[1, 5],
                        allowCross=False,
                        tooltip={"placement": "bottom"})
                ], width=4),
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
            feature_input := dcc.Dropdown(features, multi=True),
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
                html.Div(dbc.Label('Feature Space'), style={'textAlign': 'center'}, className="bg-primary text-white"),
                fas_mov_drop := dcc.Dropdown(['FAS default', 'TMHMM&SignalP', 'LCR', 'Disorder'],
                                             clearable=False, value='FAS default'),
                controller_load := dbc.Button('Load Table', color="primary", className="me-1"),
            ]),
        ], width=2),
        dbc.Col([
            dbc.Row(filter_options),
            dbc.Row(selector_options),
            dcc.Loading(id='gene_table_loading', type="default", children=[
                gene_table := dash_table.DataTable(
                    columns=[
                        {'name': 'GeneID', 'id': 'geneid', 'type': 'text'},
                        {'name': 'Movement RMSD', 'id': 'rmsd', 'type': 'numeric'},
                        {'name': 'Expressed Isoforms', 'id': '#isoforms', 'type': 'numeric'},
                        {'name': 'Replicate Coherency [max]', 'id': 'max_check', 'type': 'text'},
                        {'name': 'Replicate Coherency [std]', 'id': 'std_check', 'type': 'text'},
                        {'name': 'Max TSL', 'id': 'max_tsl', 'type': 'numeric'}
                    ],
                    data=[],
                    filter_action='native',
                    page_size=10,

                    style_data={
                        'width': '150px', 'minWidth': '150px', 'maxWidth': '150px',
                        'overflow': 'hidden',
                        'textOverflow': 'ellipsis',
                    }
                ),
            ]),

        ], width=10),
    ]),
])

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
            dbc.Row([
                exp_graph := dcc.Graph(),
            ]),
            dbc.Row([
                dbc.Col([
                    exp_png := dbc.Button("Download .png", color="primary", className="me-1"),
                ], width=6),
                dbc.Col([
                    exp_svg := dbc.Button("Download .svg", color="primary", className="me-1"),
                ], width=6),
            ]),
        ]),
        dbc.Row([
            dbc.Row([
                dbc.Col([
                    html.Div(dbc.Label("Sort By"), style={'textAlign': 'center'}),
                    sort_exp_drop := dcc.Dropdown(['Transcript ID', 'Condition', 'Replicate'
                                                   'Expression'],
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
                        {'name': 'Expression [%]', 'id': 'expression', 'type': 'numeric'}
                    ],
                    data=exp_mock,
                    filter_action='native',
                    page_size=10,

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

mov_vis = dcc.Tab(label="Movement Visualization", children=[
    dbc.Row([
        html.H2(children=['',
                          html.Div(id='mov_header',
                                   style={'display': 'inline', 'textAlign': 'center'})],
                style={'textAlign': 'center'})
    ]),
    dbc.Row([
        dbc.Col([
            mov_dropdown := dcc.Dropdown(
                value='Unscaled',
                clearable=False,
                options=['Unscaled', 'Scaled', 'Minmax_example']
            ),
            mov_png := dbc.Button("Download .png", color="primary", className="me-1"),
            mov_svg := dbc.Button("Download .svg", color="primary", className="me-1"),
        ], width=2),
        dbc.Col([
            mov_graph := dcc.Graph(),
        ]),
    ]),
    dbc.Row([
        dbc.Col([
            mov_graph2 := dcc.Graph(
                figure=px.box(pd.DataFrame(m_mock0),
                              x='transcriptid', y='Unscaled',
                              color='Condition', range_y=[0, 1])
            ),
        ], width=5),
        dbc.Col([
            dbc.Row([
                dbc.Col([
                    html.Div(dbc.Label("Sort By"), style={'textAlign': 'center'}),
                    sort_mov_drop := dcc.Dropdown(['Transcript ID', 'Condition',
                                                   'Movement [Unscaled]',
                                                   'Movement [Scaled]'],
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
                        {'name': 'TranscriptID', 'id': 'transcriptid', 'type': 'text'},
                        {'name': 'Condition', 'id': 'Condition', 'type': 'text'},
                        {'name': 'Scaled Movement', 'id': 'Scaled', 'type': 'numeric'},
                        {'name': 'Unscaled Movement', 'id': 'Unscaled', 'type': 'numeric'}
                    ],
                    data=m_mock0,
                    filter_action='native',
                    page_size=10,

                    style_data={
                        'width': '150px', 'minWidth': '150px', 'maxWidth': '150px',
                        'overflow': 'hidden',
                        'textOverflow': 'ellipsis',
                    }
                ),
            ]),
        ], width=4),
    ]),
], disabled=True)


### Isoform FAS

##

tap_node = dbc.Card([
    tap_node_header := dbc.CardHeader("Isoform: "),
    html.Div(
        [
            tap_node_url := html.A("Ensembl", href='https://www.ensembl.org/', target="_blank"),
#            tap_node_t_id := html.P("Transcript ID: "),
#            tap_node_tsl := html.P("TSL: "),
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
    fas_dropdown := dcc.Dropdown(
        value=conditions[-1],
        clearable=False,
        options=[
            {'label': name, 'value': name}
            for name in ['Full FAS', 'tmhmm+signalp', 'lcr']
        ]
    ),
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
    fas_png := dbc.Button("Download .png", color="primary", className="me-1"),
    fas_svg := dbc.Button("Download .svg", color="primary", className="me-1"),
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
                        elements=fas_mock
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

feature_architecture = dcc.Tab(label="Feature Architecture", children=[
    dbc.Row([
        fa_plot := dcc.Graph(
            figure=px.line()
        ),
    ]),
], disabled=True)

#### Exp Analysis

exp_an = dcc.Tab(id='exp_analysis', label='Expression Analysis', children=[
    dbc.Row([
        dcc.Tabs([
            gene_selector,
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
            lib_expl,
            exp_an,
        ]),
    ),
    html.Div(id="hidden-data-value", style=dict(display="none"), **{
          "data-value-1": "hello",
          "data-value-2": "false"
    }),
    mov_store := dcc.Store(data=m_mock0, id='mov_store'),
    mov_store2 := dcc.Store(data=m_mock1, id='mov_store2'),
    exp_store := dcc.Store(data=exp_mock, id='exp_store'),
    exp_store2 := dcc.Store(data={samples[0]: 0.0, samples[1]: 0.0}, id='exp_store2'),
    fa_store := dcc.Store(data=fa_mock, id='fa_store'),
    result_details := dcc.Store(data={'path': 'Not selected', 'conditions': [], 'species': 'None',
                                      'version': 'None', 'FAS modes': 'None', 'replicates': []},
                                id='result_detail'),
    result_data := dcc.Store(data=[{'geneid': 'None', '#isoforms': 0, 'rmsd': 0.0, 'max_tsl': 0,
                                    'std_check': 'No', 'max_check': 'No'}], id='result_data'),
    exp_data := dcc.Store(data={}, id='exp_data'),
    isoform_data := dcc.Store(data={}, id='isoform_data'),
    c_data := dcc.Store(data=['c1', 'c2'], id='c_data'),
    genes := dcc.Store(data=[], id='genes'),
    html.Datalist(id='list-features',
                  children=[html.Option(value=word) for word in features]),
    html.Datalist(id='list-genes',
                  children=[html.Option(value=word) for word in genes])
])




########### Callbacks

### main page
@app.callback(
    Output(result_details, 'data'),
    Output(result_input, 'valid'),
    Output(result_input, 'invalid'),
    Output('exp_analysis', 'disabled'),
    State(result_input, 'value'),
    Input(result_load_button, 'n_clicks'),
)
def load_result_data(path, button):
    config_path = path + '/result_config.json'
    ctx = dash.callback_context
    if ctx.triggered:
        if os.path.exists(config_path):
            conditions, species, release, fas_modes, replicates = read_data.read_config_file(config_path)
            return ({'path': path, 'conditions': conditions, 'species': species, 'version': release,
                    'FAS modes': fas_modes, 'replicates': replicates},
                    True, False, False)
        else:
            return ({'path': 'Not selected', 'conditions': [], 'species': 'None', 'version': 'None', 'FAS modes': [],
                    'replicates': []},
                    False, True, True)
    else:
        return ({'path': 'Not selected', 'conditions': [], 'species': 'None', 'version': 'None', 'FAS modes': [],
                'replicates': []},
                False, True, True)


@app.callback(
    Output(condition_table, 'data'),
    Output(fas_modes_table, 'data'),
    Output(c1_drop, 'options'),
    Output(c1_drop, 'value'),
    Output(fas_mov_drop, 'options'),
    Output(fas_mov_drop, 'value'),
    Input(result_details, 'data'),
)
def update_result_data(data):
    datatable = []
    datatable2 = []
    c1 = None
    fmode = None
    for condition in data['conditions']:
        datatable.append({'id': condition, 'replicates': len(data['replicates'][condition])})
    for mode in data['FAS modes']:
        datatable2.append({'id': mode})
    if len(data['conditions']) > 0:
        c1 = data['conditions'][0]
    if len(data['FAS modes']) > 0:
        fmode = data['FAS modes'][0]
    return datatable, datatable2, data['conditions'], c1, data['FAS modes'], fmode


@app.callback(
    Output(result_data, 'data'),
    Output(exp_data, 'data'),
    Output(isoform_data, 'data'),
    Output(c_data, 'data'),
    Output(genes, 'data'),
    Input(controller_load, 'n_clicks'),
    State(c1_drop, 'value'),
    State(c2_drop, 'value'),
    State(fas_mov_drop, 'value'),
    State(result_details, 'data'),
    State(result_data, 'data'),
    State(exp_data, 'data'),
    State(isoform_data, 'data'),
    State(c_data, 'data'),
    State(genes, 'data')
)
def load_result_table(nclicks, c1, c2, fmode, r_details, old_r_data, old_exp_data, old_iso_data, old_c_data,
                      old_gene_data):
    ctx = dash.callback_context
    result_data = old_r_data
    exp_data = old_exp_data
    isoform_data = old_iso_data
    c_data = old_c_data
    genes = old_gene_data
    if ctx.triggered:
        if c1 and c2 and fmode:
            path = r_details['path']
            filepath = None
            if os.path.isdir(path + '/main_comparison/' + c1 + '@' + c2):
                filepath = path + '/main_comparison/' + c1 + '@' + c2 + '/result_' + c1 + '@' + c2 + '_' + fmode + '.tsv'
            elif os.path.isdir(path + '/main_comparison/' + c2 + '@' + c1):
                filepath = path + '/main_comparison/' + c2 + '@' + c1 + '/result_' + c2 + '@' + c1 + '_' + fmode + '.tsv'
            if filepath:
                result_data, genes, isoform_data = read_data.read_results_main(filepath)
                exp_data = read_data.read_results_exp(path, c1, c2)
                c_data = [c1, c2]
    return result_data, exp_data, isoform_data, c_data, genes


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
    Input(row_drop, 'value'),
    Input(rmsd_slider_0, 'value'),
    Input(coherence_drop, 'value'),
    Input(fold_slider, 'value'),
    Input(feature_input, 'value'),
    Input(sort2_drop, 'value'),
    Input(tsl_slider, 'value'),
    Input(ao_switch, 'value'),
    Input(result_data, 'data'),
)
def update_table_options(row_v, rmsd_v, coherence, fold_v, feature_v, sort2_v, tsl_v, ao_switch,
                         result_data):
    dff = pd.DataFrame(result_data)
    dff = dff[(dff['rmsd'] >= rmsd_v[0]) & (dff['rmsd'] <= rmsd_v[1])]
    dff = dff[(dff['max_tsl'] >= tsl_v[0]) & (dff['max_tsl'] <= tsl_v[1])]
    if coherence == 'Coherence [std]':
        dff = dff[(dff['std_check'] == 'Yes')]
    elif coherence == 'Coherence [max]':
        dff = dff[(dff['max_check'] == 'Yes')]

#    fslider_d = {1: -2, 2: -1, 3: -0.5, 4: 0, 5: 0.5, 6: 1, 7: 2}
#    if not fold_v[0] == 0:
#        dff = dff[dff['foldchange'] >= fslider_d[fold_v[0]]]
#    if not fold_v[1] == 8:
#        dff = dff[dff['foldchange'] <= fslider_d[fold_v[1]]]

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

    return dff.to_dict('records'), row_v


@app.callback(
    Output(gene_input, 'value'),
    Input(gene_table, 'active_cell'),
    State(gene_table, 'derived_virtual_data'),
)
def select_gene_table(active_cell, gene_table):
    if active_cell:
        if active_cell['column_id'] == 'geneid':
            return gene_table[active_cell['row']][active_cell['column_id']]


@app.callback(
    Output(gene_url, 'href'),
    Input(gene_input, 'value'),
)
def create_gene_url(geneid):
    if geneid in genes:
        return 'http://www.ensembl.org/id/' + geneid
    else:
        return 'http://www.ensembl.org/'


@app.callback(
    Output('FAS_header', 'children'),
    Output('mov_header', 'children'),
    Output('exp_header', 'children'),
    Output(mov_store, 'data'),
    Output(mov_store2, 'data'),
    Output(exp_store, 'data'),
    Output(exp_store2, 'data'),
    Output(transcript_dropdown, 'options'),
    Output(transcript_dropdown, 'value'),
    Output(gene_input, 'valid'),
    Output(gene_input, 'invalid'),
    Output(expression_stats, 'disabled'),
    Input(gene_select, 'n_clicks'),
    State(gene_input, 'value'),
    State(exp_data, 'data'),
    State(isoform_data, 'data'),
    State(c_data, 'data'),
    State(genes, 'data')
)
def load_button(n_clicks, geneid, exp_data, isoform_dict, samples, genes):
    ctx = dash.callback_context
    if ctx.triggered:
        if geneid in genes:
            # mov_table, mov_fig = read_data.prepare_movement(movement, geneid)
            return (geneid, geneid, geneid, m_mock0, m_mock1, exp_data[geneid]['table'],
                    {samples[0]: exp_data[geneid][samples[0]], samples[1]: exp_data[geneid][samples[1]]},
                    isoform_dict[geneid], isoform_dict[geneid], True, False, False)
        else:
            return ('Example', 'Example', 'Example', m_mock0, m_mock1, exp_mock,
                    {samples[0]: 0.0, samples[1]: 0.0}, ['t1', 't2', 't3'], ['t1', 't2', 't3'], False, True, True)
    else:
        return ('Example', 'Example', 'Example', m_mock0, m_mock1, exp_mock,
                {samples[0]: 0.0, samples[1]: 0.0}, ['t1', 't2', 't3'], ['t1', 't2', 't3'], False, False, True)


# Expression callbacks
@app.callback(
    Output(exp_graph, 'figure'),
    Input(exp_store, 'data'),
    Input(exp_store2, 'data'),
    State(gene_input, 'value')
)
def generate_chart(exp_data, exp_gene, gid):
    dff = pd.DataFrame(exp_data)
    fig = px.box(dff, x='transcriptid', y='expression', color='Condition', range_y=[0, 1])
    fig.update_layout(title_text=gid)
    return fig


@app.callback(
    Output(exp_table, 'data'),
    Input(sort_exp_drop, 'value'),
    Input(order_exp_drop, 'value'),
    Input(exp_store, 'data')
)
def generate_exp_table(sort_exp, order_exp, exp_data):
    sdrop_d = {'Transcript ID': 'transcriptid', 'Condition': 'Condition', 'Expression': 'expression'}
    direct = (order_exp == 'Ascending')
    exp_data = pd.DataFrame(exp_data)
    if sort_exp:
        exp_data = exp_data.sort_values(sdrop_d[sort_exp], ascending=direct)
    return exp_data.to_dict('records')


# Movement callbacks
@app.callback(
    Output(mov_graph, 'figure'),
    Input(mov_dropdown, 'value'),
    Input(mov_store2, 'data'),
    State(gene_input, 'value')
)
def generate_polygon_figure(plottype, mov_fig, geneid):
    fig = None
    r1 = []
    r2 = []
    gid = 'Example'
    if geneid:
        gid = geneid
    if mov_fig:
        if plottype == 'Scaled':
            r1 = mov_fig[samples[0]]['scaled']
            r2 = mov_fig[samples[1]]['scaled']
        elif plottype == 'Unscaled':
            r1 = mov_fig[samples[0]]['unscaled']
            r2 = mov_fig[samples[1]]['unscaled']
        fig = make_subplots(rows=1, cols=2, specs=[[{'type': 'polar'}] * 2], subplot_titles=samples)
        fig.add_trace(go.Scatterpolar(r=r1,
                                      theta=mov_fig['isoforms'],
                                      fill='toself'),
                      row=1, col=1)
        fig.add_trace(go.Scatterpolar(r=r2,
                                      theta=mov_fig['isoforms'],
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
            title_text=gid
        )

        if plottype == 'Minmax_example':
            categories = ['t1', 't2', 't3', 't1']

            fig = make_subplots(rows=1, cols=2, specs=[[{'type': 'polar'}] * 2], subplot_titles=['c1', 'c2'])

            fig.add_trace(go.Scatterpolar(
                r=[0, 0.4, 0.3, 0],
                theta=categories,
                marker=dict(color="red", opacity=0.0),
                line=dict(color=None, width=0, shape='spline'),
                fillcolor='red',
                opacity=0.2,
                name='Min A'),
                row=1, col=1
            )
            fig.add_trace(go.Scatterpolar(
                r=[0.15, 0.45, 0.365, 0.15],
                theta=categories,
                marker=dict(color="red", opacity=0.0),
                line=dict(color=None, width=0, shape='spline'),
                fillcolor='red',
                fill='tonext',
                opacity=0.4,
                name='Low STD A'),
                row=1, col=1
            )
            fig.add_trace(go.Scatterpolar(
                r=[0.3, 0.5, 0.45, 0.3],
                theta=categories,
                marker=dict(color="red", opacity=0.0),
                line=dict(color='red', width=3, shape='spline'),
                fillcolor='red',
                fill='tonext',
                opacity=0.4,
                name='Mean A'),
                row=1, col=1
            )
            fig.add_trace(go.Scatterpolar(
                r=[0.45, 0.55, 0.525, 0.45],
                theta=categories,
                marker=dict(color="red", opacity=0.0),
                line=dict(color=None, width=0, shape='spline'),
                fillcolor='red',
                fill='tonext',
                opacity=0.2,
                name='High STD A'),
                row=1, col=1
            )
            fig.add_trace(go.Scatterpolar(
                r=[0.6, 0.6, 0.6, 0.6],
                theta=categories,
                marker=dict(color="red", opacity=0.0),
                line=dict(color=None, width=0, shape='spline'),
                fillcolor='red',
                fill='tonext',
                opacity=0.2,
                name='Max A'),
                row=1, col=1
            )

            fig.add_trace(go.Scatterpolar(
                r=[0.2, 0, 0.5, 0.2],
                theta=categories,
                marker=dict(color="blue", opacity=0.0),
                line=dict(color=None, width=0, shape='spline'),
                fillcolor='blue',
                opacity=0.2,
                name='Min B'),
                row=1, col=2
            )
            fig.add_trace(go.Scatterpolar(
                r=[0.35, 0, 0.565, 0.35],
                theta=categories,
                marker=dict(color="blue", opacity=0.0),
                line=dict(color=None, width=0, shape='spline'),
                fillcolor='blue',
                fill='tonext',
                opacity=0.4,
                name='Low STD B'),
                row=1, col=2
            )
            fig.add_trace(go.Scatterpolar(
                r=[0.5, 0, 0.65, 0.5],
                theta=categories,
                marker=dict(color="blue", opacity=0.0),
                line=dict(color='blue', width=3, shape='spline'),
                fillcolor='blue',
                fill='tonext',
                opacity=0.4,
                name='Mean B'),
                row=1, col=2
            )
            fig.add_trace(go.Scatterpolar(
                r=[0.65, 0, 0.725, 0.65],
                theta=categories,
                marker=dict(color="blue", opacity=0.0),
                line=dict(color=None, width=0, shape='spline'),
                fillcolor='blue',
                fill='tonext',
                opacity=0.2,
                name='High STD B'),
                row=1, col=2
            )
            fig.add_trace(go.Scatterpolar(
                r=[0.8, 0, 0.8, 0.8],
                theta=categories,
                marker=dict(color="blue", opacity=0.0),
                line=dict(color=None, width=0, shape='spline'),
                fillcolor='blue',
                fill='tonext',
                opacity=0.2,
                name='Max B'),
                row=1, col=2
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
                title_text=gid
            )
    return fig


@app.callback(
    Output(mov_graph2, 'figure'),
    Output(mov_table, 'data'),
    Input(mov_dropdown, 'value'),
    Input(mov_store, 'data'),
    Input(sort_mov_drop, 'value'),
    Input(order_mov_drop, 'value'),
)
def generate_boxplot(plottype, mov_table, sort_mov, order_mov):
    figure = px.box(pd.DataFrame(mov_table),
                    x='transcriptid', y='Unscaled',
                    color='Condition', range_y=[0, 1])
    if not plottype == 'Minmax_example':
        figure = px.box(pd.DataFrame(mov_table),
                        x='transcriptid', y=plottype,
                        color='Condition', range_y=[0, 1])
    figure.update_layout(
        xaxis_title='Isoform',
        yaxis_title='Movement'
    )
    mov_table = pd.DataFrame(mov_table)
    sdrop_d = {'Transcript ID': 'transcriptid', 'Condition': 'Condition', 'Movement [Unscaled]': 'Unscaled',
               'Movement [Scaled]': 'Scaled'}
    direct = (order_mov == 'Ascending')
    if sort_mov:
        mov_table = mov_table.sort_values(sdrop_d[sort_mov], ascending=direct)
    return figure, mov_table.to_dict('records')



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
    State(gene_input, 'value'),
)
def fas_cytoscape_figure(transcripts, node_labels, edge_labels, directional, label_size, toggle_zero, geneid):
    if geneid in genes:
        fas_graph = read_data.prepare_FAS_graph(fas_dict[geneid], transcripts, pd.DataFrame(exp_data[geneid]['table']),
                                                'all_H1', 'all_embryo_ED', directional, toggle_zero)
    else:
        fas_graph = fas_mock
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
    Input('fas_figure', 'tapNodeData')
)
def display_tap_node_data(data):
    if data:
        return "Isoform: " + data['label'], 'http://www.ensembl.org/id/' + data['label']
    else:
        return 'Isoform: NA', 'http://www.ensembl.org/'


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
    Output(fa_plot, 'figure'),
    Input(fas_dropdown, 'value'),
    Input(fa_store, 'data'),
    State(gene_input, 'value')
)
def generate_fa_plot(scoring, fa_data, gid):
    df1 = pd.DataFrame(fa_data['isoform1'])
    df2 = pd.DataFrame(fa_data['isoform2'])
    fig = make_subplots(rows=2, cols=1, specs=[[{'type': 'scatter'}], [{'type': 'scatter'}]],
                        subplot_titles=['isoform1', 'isoform2'])
    tmpfig1 = px.line(df1, x='position', y='y', color='color')
    tmpfig2 = px.line(df2, x='position', y='y', color='color')
    for trace in range(len(tmpfig1["data"])):
        fig.add_trace(tmpfig1["data"][trace], row=1, col=1)
    for trace in range(len(tmpfig2["data"])):
        fig.add_trace(tmpfig2["data"][trace], row=2, col=1)
    fig.update_traces(connectgaps=False, line=dict(width=20))
    fig.update_yaxes(showticklabels=False)
    fig.update_layout(yaxis_title='')
    return fig


def main():
    app.run_server(debug=True)


if __name__ == '__main__':
    main()

