import dash
import dash_cytoscape as cyto
from dash import html
from dash import dcc
from dash import dash_table
from dash.dependencies import Input, Output, State
import dash_bootstrap_components as dbc
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from sys import argv
from expNet.help_functions import get_colorscale
from expNet import read_data


cyto.load_extra_layouts()

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.CERULEAN])

################################## Mock data
genedata = [
    {'geneid': '001', 'rmsd': 0.5, 'expchange': 0.8, 'foldchange': 0.1, 'isoforms': 3},
    {'geneid': '002', 'rmsd': 0.05, 'expchange': 0.1, 'foldchange': 1.0, 'isoforms': 3},
    {'geneid': '003', 'rmsd': 0.2, 'expchange': 0.7, 'foldchange': -1.1, 'isoforms': 3},
    {'geneid': '004', 'rmsd': 0.8, 'expchange': 0.9, 'foldchange': -0.3, 'isoforms': 3},
    {'geneid': '005', 'rmsd': 0.9, 'expchange': 0.9, 'foldchange': 0.4, 'isoforms': 3},
    {'geneid': '006', 'rmsd': 0.3, 'expchange': 0.8, 'foldchange': 1.1, 'isoforms': 3},
    {'geneid': '007', 'rmsd': 0.4, 'expchange': 0.8, 'foldchange': -0.1, 'isoforms': 3},
    {'geneid': '008', 'rmsd': 0.1, 'expchange': 0.8, 'foldchange': 0.1, 'isoforms': 3},
    {'geneid': '009', 'rmsd': 0.8, 'expchange': 0.8, 'foldchange': 2.1, 'isoforms': 3},
    {'geneid': '010', 'rmsd': 0.2, 'expchange': 0.3, 'foldchange': 0.1, 'isoforms': 3},
    {'geneid': '011', 'rmsd': 0.5, 'expchange': 0.8, 'foldchange': 0.1, 'isoforms': 3},
    {'geneid': '012', 'rmsd': 0.05, 'expchange': 0.1, 'foldchange': 1.0, 'isoforms': 3},
    {'geneid': '013', 'rmsd': 0.2, 'expchange': 0.7, 'foldchange': -1.1, 'isoforms': 3},
    {'geneid': '014', 'rmsd': 0.8, 'expchange': 0.9, 'foldchange': -0.3, 'isoforms': 3},
    {'geneid': '015', 'rmsd': 0.9, 'expchange': 0.9, 'foldchange': 0.4, 'isoforms': 3},
    {'geneid': '016', 'rmsd': 0.3, 'expchange': 0.8, 'foldchange': 1.1, 'isoforms': 3},
    {'geneid': '017', 'rmsd': 0.4, 'expchange': 0.8, 'foldchange': -0.1, 'isoforms': 3},
    {'geneid': '018', 'rmsd': 0.1, 'expchange': 0.8, 'foldchange': 0.1, 'isoforms': 3},
    {'geneid': '019', 'rmsd': 0.8, 'expchange': 0.8, 'foldchange': 2.1, 'isoforms': 3},
    {'geneid': '020', 'rmsd': 0.2, 'expchange': 0.3, 'foldchange': 0.1, 'isoforms': 3}
]
df = pd.DataFrame(genedata)

f_dict = {
          'pfam_EF_hand': ['001', '002'],
          'smart_EF_hand': ['001', '003'],
          'tmhmm_transmembrane': ['015', '016', '020', '005'],
          'pfam_domain_2': ['010', '011']
}

features = list(f_dict.keys())
genes = ['001', '002', '003', '004', '005', '006', '007', '008', '009', '010']

conditions = ['c1', 'c2', 'both']

exp_data = [
    {'transcriptid': 'T1', 'expression': 0.5, 'expressionp': 50, 'movement': 0, 'condition': 'c1', 'sample': 'c1_1'},
    {'transcriptid': 'T2', 'expression': 0.4, 'expressionp': 40, 'movement': 0.4, 'condition': 'c1', 'sample': 'c1_1'},
    {'transcriptid': 'T3', 'expression': 0.1, 'expressionp': 10, 'movement': 0.3, 'condition': 'c1', 'sample': 'c1_1'},
    {'transcriptid': 'T1', 'expression': 0.6, 'expressionp': 60, 'movement': 0.3, 'condition': 'c1', 'sample': 'c1_2'},
    {'transcriptid': 'T2', 'expression': 0.25, 'expressionp': 25, 'movement': 0.5, 'condition': 'c1', 'sample': 'c1_2'},
    {'transcriptid': 'T3', 'expression': 0.15, 'expressionp': 15, 'movement': 2.25/5, 'condition': 'c1', 'sample': 'c1_2'},
    {'transcriptid': 'T1', 'expression': 0.55, 'expressionp': 55, 'movement': 0.6, 'condition': 'c1', 'sample': 'c1_3'},
    {'transcriptid': 'T2', 'expression': 0.45, 'expressionp': 45, 'movement': 0.6, 'condition': 'c1', 'sample': 'c1_3'},
    {'transcriptid': 'T3', 'expression': 0.0, 'expressionp': 0, 'movement': 0.6, 'condition': 'c1', 'sample': 'c1_3'},

    {'transcriptid': 'T1', 'expression': 0.1, 'expressionp': 10, 'movement': 1/5, 'condition': 'c2', 'sample': 'c2_1'},
    {'transcriptid': 'T2', 'expression': 0.1, 'expressionp': 10, 'movement': 0, 'condition': 'c2', 'sample': 'c2_1'},
    {'transcriptid': 'T3', 'expression': 0.8, 'expressionp': 80, 'movement': 0.5, 'condition': 'c2', 'sample': 'c2_1'},
    {'transcriptid': 'T1', 'expression': 0.15, 'expressionp': 15, 'movement': 0.5, 'condition': 'c2', 'sample': 'c2_2'},
    {'transcriptid': 'T2', 'expression': 0.0, 'expressionp': 0, 'movement': 0, 'condition': 'c2', 'sample': 'c2_2'},
    {'transcriptid': 'T3', 'expression': 0.85, 'expressionp': 85, 'movement': 3.25/5, 'condition': 'c2', 'sample': 'c2_2'},
    {'transcriptid': 'T1', 'expression': 0.0, 'expressionp': 0, 'movement': 0.8, 'condition': 'c2', 'sample': 'c2_3'},
    {'transcriptid': 'T2', 'expression': 0.1, 'expressionp': 10, 'movement': 0, 'condition': 'c2', 'sample': 'c2_3'},
    {'transcriptid': 'T3', 'expression': 0.9, 'expressionp': 90, 'movement': 0.8, 'condition': 'c2', 'sample': 'c2_3'},
]


fas_graph = [
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


app.layout = html.Div([
    dbc.Row(
        html.H1("ExpVis", style={'textAlign': 'center'})
    ),
    dbc.Row(
        dcc.Tabs([
            dcc.Tab(label='Gene Selector', children=[
                dbc.Row([
                    dbc.Col([
                        html.Div(dbc.Label("Sort By"), style={'textAlign': 'center'}),
                        sort_drop := dcc.Dropdown(['RMSD', 'Net Exp Change', 'Fold Change'], value='RMSD', )
                    ]),
                    dbc.Col([
                        html.Div(dbc.Label("Order"), style={'textAlign': 'center'}),
                        sort2_drop := dcc.Dropdown(['Ascending', 'Descending'], value='Descending')
                    ]),
                    dbc.Col([
                        html.Div(dbc.Label("Movement RMSD"), style={'textAlign': 'center'}),
                        rmsd_slider := dcc.RangeSlider(0, 1, 0.1,
                                                       value=[0, 1],
                                                       allowCross=False,
                                                       tooltip={"placement": "bottom"})
                    ], width=3),
                    dbc.Col([
                        html.Div(dbc.Label("Net Exp Change"), style={'textAlign': 'center'}),
                        exp_slider := dcc.RangeSlider(
                            0, 1, 0.1,
                            value=[0, 1],
                            allowCross=False,
                            tooltip={"placement": "bottom"})
                    ], width=3),
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
                    ], width=3),
                ], justify="between", className='mt-3 mb-4'),
                dbc.Row([
                    dbc.Col([
                        dbc.Label("Show number of rows"),
                        row_drop := dcc.Dropdown(value=10, clearable=False, style={'width': '35%'},
                                                 options=[10, 25, 50, 100]),
                    ]),
                ]),

                gene_table := dash_table.DataTable(
                    columns=[
                        {'name': 'GeneID', 'id': 'geneid', 'type': 'text'},
                        {'name': 'Movement RMSD', 'id': 'rmsd', 'type': 'numeric'},
                        {'name': 'Net Exp Change', 'id': 'expchange', 'type': 'numeric'},
                        {'name': 'Fold Change', 'id': 'foldchange', 'type': 'numeric'},
                        {'name': 'Expressed Isoforms', 'id': 'isoforms', 'type': 'numeric'}
                    ],
                    data=genedata,
                    filter_action='native',
                    page_size=10,

                    style_data={
                        'width': '150px', 'minWidth': '150px', 'maxWidth': '150px',
                        'overflow': 'hidden',
                        'textOverflow': 'ellipsis',
                    }
                ),
                dbc.Row([
                    dbc.Col([
                        html.Div(dbc.Label("Filter by Feature")),
                        feature_input := dcc.Input(
                            type='text',
                            list='list-features',
                            value=''
                        ),
                    ], width=3),
                    dbc.Col([
                        html.Div(dbc.Label("Select Gene")),
                        dbc.Col([
                            gene_input := dcc.Input(
                                type='text',
                                list='list-genes',
                                value=''
                            ),
                            gene_select := html.Button("Load"),
                        ]),
                        gene_select_container := html.Div(children='None selected')
                    ], width=4),
                ]),
            ]),
            dcc.Tab(label="Expression Statistics", children=[
                dbc.Row([
                    dbc.Col([
                        exp_dropdown := dcc.Dropdown(
                            value=conditions[-1],
                            clearable=False,
                            options=[
                                {'label': name.capitalize(), 'value': name}
                                for name in conditions
                            ]
                        ),
                        exp_png := html.Button("Download .png"),
                        exp_svg := html.Button("Download .svg"),
                    ], width=2),
                    dbc.Col([
                        exp_graph := dcc.Graph(),
                    ], width=5),
                    dbc.Col([
                        exp_table := dash_table.DataTable(
                            columns=[
                                {'name': 'TranscriptID', 'id': 'transcriptid', 'type': 'text'},
                                {'name': 'Condition', 'id': 'condition', 'type': 'text'},
                                {'name': 'Sample', 'id': 'sample', 'type': 'text'},
                                {'name': 'Expression [%]', 'id': 'expressionp', 'type': 'numeric'}
                            ],
                            data=exp_data,
                            filter_action='native',
                            page_size=10,

                            style_data={
                                'width': '150px', 'minWidth': '150px', 'maxWidth': '150px',
                                'overflow': 'hidden',
                                'textOverflow': 'ellipsis',
                            }
                        ),
                    ], width=4),
                ]),
                dbc.Row([
                    dbc.Col([

                    ], width=2),
                    dbc.Col([

                    ], width=6),

                ]),
            ]),
            dcc.Tab(label="Movement Visualization", children=[
                dbc.Row([
                    dbc.Col([
                        mov_dropdown := dcc.Dropdown(
                            value='Radarplot',
                            clearable=False,
                            options=['Radarplot', 'Boxplot']
                        ),
                        mov_png := html.Button("Download .png"),
                        mov_svg := html.Button("Download .svg"),
                    ], width=2),
                    dbc.Col([
                        mov_graph := dcc.Graph(),
                    ]),
                ]),
                dbc.Row([
                    dbc.Col([

                    ], width=2),
                    dbc.Col([
                        mov_graph2 := dcc.Graph(
                            figure = px.box(pd.DataFrame(exp_data), x='transcriptid', y='movement', color='condition', range_y=[0, 1])
                        ),
                    ], width=5),
                    dbc.Col([
                        mov_table := dash_table.DataTable(
                            columns=[
                                {'name': 'TranscriptID', 'id': 'transcriptid', 'type': 'text'},
                                {'name': 'Condition', 'id': 'condition', 'type': 'text'},
                                {'name': 'Sample', 'id': 'sample', 'type': 'text'},
                                {'name': 'Movement', 'id': 'movement', 'type': 'numeric'}
                            ],
                            data=exp_data,
                            filter_action='native',
                            page_size=10,

                            style_data={
                                'width': '150px', 'minWidth': '150px', 'maxWidth': '150px',
                                'overflow': 'hidden',
                                'textOverflow': 'ellipsis',
                            }
                        ),
                    ], width=4),
                ]),
            ]),
            dcc.Tab(label="Isoform FAS", children=[
                dbc.Row([
                    dbc.Col([
                        fas_dropdown := dcc.Dropdown(
                            value=conditions[-1],
                            clearable=False,
                            options=[
                                {'label': name, 'value': name}
                                for name in ['Full FAS', 'tmhmm+signalp', 'lcr']
                            ]
                        ),
                        fas_png := html.Button("Download .png"),
                        fas_svg := html.Button("Download .svg"),
                    ], width={'size': 2}),
                    dbc.Col([
                        cyto.Cytoscape(
                            id='cytoscape-C',
                            layout={'name': 'circle', 'directed': True},
                            style={'width': '100%', 'height': '40vh'},
                            stylesheet=[
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
                            ],
                            elements=fas_graph
                        )
                    ], width={'size': 6}),
                ]),
            ]),
        ])
    ),
    html.Div(id="hidden-data-value", style=dict(display="none"), **{
          "data-value-1": "hello",
          "data-value-2": "false"
    }),
    dcc.Store(id='gene_data', data=genedata),
    html.Datalist(id='list-features',
                  children=[html.Option(value=word) for word in features]),
    html.Datalist(id='list-genes',
                  children=[html.Option(value=word) for word in genes])
])


# GeneSelector Callback
@app.callback(
    Output(gene_table, 'data'),
    Output(gene_table, 'page_size'),
    Input(row_drop, 'value'),
    Input(rmsd_slider, 'value'),
    Input(exp_slider, 'value'),
    Input(fold_slider, 'value'),
    Input(feature_input, 'value'),
    Input(sort_drop, 'value'),
    Input(sort2_drop, 'value')
)
def update_table_options(row_v, rmsd_v, exp_v, fold_v, feature_v, sort_v, sort2_v):
    dff = df.copy()

#    if cont_v:
#        dff = dff[dff.continent==cont_v]
#    if country_v:
#        dff = dff[dff.country.isin(country_v)]

    dff = dff[(dff['rmsd'] >= rmsd_v[0]) & (dff['rmsd'] <= rmsd_v[1])]
    dff = dff[(dff['expchange'] >= exp_v[0]) & (dff['expchange'] <= exp_v[1])]

    fslider_d = {1: -2, 2: -1, 3: -0.5, 4: 0, 5: 0.5, 6: 1, 7: 2}
    if not fold_v[0] == 0:
        dff = dff[dff['foldchange'] >= fslider_d[fold_v[0]]]
    if not fold_v[1] == 8:
        dff = dff[dff['foldchange'] <= fslider_d[fold_v[1]]]

    if feature_v:
        if feature_v in f_dict:
            dff = dff[dff['geneid'].isin(f_dict[feature_v])]

#    dff = dff[(dff['lifeExp'] >= life_v) & (dff['lifeExp'] < 100)]
    sdrop_d = {'RMSD': 'rmsd', 'Net Exp Change': 'expchange', 'Fold Change': 'foldchange'}
    direct = (sort2_v == 'Ascending')
    if sort_v:
        dff = dff.sort_values(sdrop_d[sort_v], ascending=direct)

    return dff.to_dict('records'), row_v


@app.callback(
    Output(gene_select, 'hidden'),
    Input(gene_input, 'value')
)
def show_load_button(g_input_v):
    if g_input_v in genes:
        return False
    else:
        return True


@app.callback(
    Output(gene_select_container, 'children'),
    Input(gene_select, 'n_clicks'),
    State(gene_input, 'value')
)
def load_gene_button(g_select_v, g_input_v):
    ctx = dash.callback_context
    if ctx.triggered:
        return 'Selected {}'.format(g_input_v)


# Expression callbacks
@app.callback(
    Output(exp_graph, "figure"),
    Input(exp_dropdown, "value"),
)
def generate_chart(condition):
    dff = pd.DataFrame(exp_data)
    if not condition == 'both':
        dff = dff = dff[dff['condition'] == condition]
    fig = px.box(dff, x='transcriptid', y='expression', color='condition', range_y=[0, 1])
    return fig


# Movement callbacks
@app.callback(
    Output(mov_graph, "figure"),
    Input(mov_dropdown, "value"),
)
def generate_chart(plottype):
    dff = pd.DataFrame(exp_data)
    fig = None
    if plottype == 'Boxplot':
        fig = px.box(dff, x='transcriptid', y='movement', color='condition', range_y=[0, 1])
    elif plottype == 'Radarplot':
        categories = ['t1', 't2', 't3', 't1']

        fig = make_subplots(rows=1, cols=2, specs=[[{'type': 'polar'}] * 2], subplot_titles=['c1', 'c2'])
        # fig = go.Figure()

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
            showlegend=False
        )
    return fig


def main():
    app.run_server(debug=True)


if __name__ == '__main__':
    main()

