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
import read_data


from sys import argv



genedata, samples, movement, genes, isoform_dict = read_data.read_main_results('/home/julian/Data/expVis_data/polygonFAS_all_H1xall_embryo_ED_sorted.tsv')
f_dict = read_data.read_json('/home/julian/Data/expVis_data/isoforms_important_features.json')
exp_data = read_data.read_exp_input('/home/julian/Data/expVis_data/mov_exp_data/polygonFAS_H1.tsv', samples[0], {})
exp_data = read_data.read_exp_input('/home/julian/Data/expVis_data/mov_exp_data/polygonFAS_embryo_ED.tsv', samples[1], exp_data)
fas_dict = read_data.read_json('/home/julian/Data/expVis_data/FAS_scores/distance_master.json')

## Mock Data ##
exp_mock = [
    {'transcriptid': 'example1', 'Condition': samples[0], 'expression': 0.5},
    {'transcriptid': 'example2', 'Condition': samples[0], 'expression': 0.5},
    {'transcriptid': 'example3', 'Condition': samples[0], 'expression': 0.5},
    {'transcriptid': 'example1', 'Condition': samples[1], 'expression': 0.5},
    {'transcriptid': 'example2', 'Condition': samples[1], 'expression': 0.5},
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
##

cyto.load_extra_layouts()

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.CERULEAN])


df = pd.DataFrame(genedata)

#f_dict = {
#          'pfam_EF_hand': ['001', '002'],
#          'smart_EF_hand': ['001', '003'],
#          'tmhmm_transmembrane': ['015', '016', '020', '005'],
#          'pfam_domain_2': ['010', '011']
#}

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
                        sort_drop := dcc.Dropdown(['RMSD [unscaled]', 'RMSD [scaled]'], value='RMSD [unscaled]', )
                    ]),
                    dbc.Col([
                        html.Div(dbc.Label("Order"), style={'textAlign': 'center'}),
                        sort2_drop := dcc.Dropdown(['Ascending', 'Descending'], value='Descending')
                    ]),
                    dbc.Col([
                        html.Div(dbc.Label("Movement RMSD [unscaled]"), style={'textAlign': 'center'}),
                        rmsd_slider_0 := dcc.RangeSlider(0, 1, 0.1,
                                                         value=[0, 1],
                                                         allowCross=False,
                                                         tooltip={"placement": "bottom"})
                    ], width=3),
                    dbc.Col([
                        html.Div(dbc.Label("Movement RMSD [scaled]"), style={'textAlign': 'center'}),
                        rmsd_slider_1 := dcc.RangeSlider(
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
                    dbc.Col([
                        html.Div(dbc.Label("Transcript Support Level"), style={'textAlign': 'center'}),
                        tsl_slider := dcc.RangeSlider(
                            1, 5, 1,
                            value=[1, 5],
                            allowCross=False,
                            tooltip={"placement": "bottom"})
                    ], width=3),
                ]),

                gene_table := dash_table.DataTable(
                    columns=[
                        {'name': 'GeneID', 'id': 'geneid', 'type': 'text'},
                        {'name': 'Movement RMSD [unscaled]', 'id': 'unscaled_rmsd', 'type': 'numeric'},
                        {'name': 'Movement RMSD [scaled]', 'id': 'scaled_rmsd', 'type': 'numeric'},
                        {'name': 'Expressed Isoforms', 'id': '#isoforms', 'type': 'numeric'},
                        {'name': 'Max TSL', 'id': 'max_tsl', 'type': 'numeric'}
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
                                {'label': name, 'value': name}
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
                                {'name': 'Condition', 'id': 'Condition', 'type': 'text'},
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
                            value='Unscaled',
                            clearable=False,
                            options=['Unscaled', 'Scaled', 'Minmax_example']
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
                            figure=px.box(pd.DataFrame(m_mock0),
                                          x='transcriptid', y='Unscaled',
                                          color='Condition', range_y=[0, 1])
                        ),
                    ], width=5),
                    dbc.Col([
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
                        fas_figure := cyto.Cytoscape(
                            id='fas_figure',
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
                            elements=fas_mock
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
    mov_store := dcc.Store(data=m_mock0, id='mov_store'),
    mov_store2 := dcc.Store(data=m_mock1, id='mov_store2'),
    exp_store := dcc.Store(data=exp_mock, id='exp_store'),
    exp_store2 := dcc.Store(data={samples[0]: 0.0, samples[1]: 0.0}, id='exp_store2'),
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
    Input(rmsd_slider_0, 'value'),
    Input(rmsd_slider_1, 'value'),
    Input(fold_slider, 'value'),
    Input(feature_input, 'value'),
    Input(sort_drop, 'value'),
    Input(sort2_drop, 'value'),
    Input(tsl_slider, 'value')
)
def update_table_options(row_v, rmsd_unscaled_v, rmsd_scaled_v, fold_v, feature_v, sort_v, sort2_v, tsl_v):
    dff = df.copy()

    dff = dff[(dff['unscaled_rmsd'] >= rmsd_unscaled_v[0]) & (dff['unscaled_rmsd'] <= rmsd_unscaled_v[1])]
    dff = dff[(dff['scaled_rmsd'] >= rmsd_scaled_v[0]) & (dff['scaled_rmsd'] <= rmsd_scaled_v[1])]
    dff = dff[(dff['max_tsl'] >= tsl_v[0]) & (dff['max_tsl'] <= tsl_v[1])]

    fslider_d = {1: -2, 2: -1, 3: -0.5, 4: 0, 5: 0.5, 6: 1, 7: 2}
#    if not fold_v[0] == 0:
#        dff = dff[dff['foldchange'] >= fslider_d[fold_v[0]]]
#    if not fold_v[1] == 8:
#        dff = dff[dff['foldchange'] <= fslider_d[fold_v[1]]]

    if feature_v:
        if feature_v in f_dict:
            dff = dff[dff['geneid'].isin(f_dict[feature_v])]

    sdrop_d = {'RMSD [unscaled]': 'unscaled_rmsd', 'RMSD [scaled]': 'scaled_rmsd'}
    direct = (sort2_v == 'Ascending')
    if sort_v:
        dff = dff.sort_values(sdrop_d[sort_v], ascending=direct)

    return dff.to_dict('records'), row_v


@app.callback(
    Output(gene_select, 'hidden'),
    Input(gene_input, 'value')
)
def show_load_button(g_input_v):
    if g_input_v:
        if g_input_v in genes:
            return False
        else:
            return True
    else:
        return True


@app.callback(
    Output(gene_select_container, 'children'),
    Output(mov_store, 'data'),
    Output(mov_store2, 'data'),
    Output(exp_store, 'data'),
    Output(exp_store2, 'data'),
    Output(fas_figure, 'elements'),
    Input(gene_select, 'n_clicks'),
    State(gene_input, 'value'),
)
def load_button(n_clicks, geneid):
    ctx = dash.callback_context
    if ctx.triggered:
        fas_graph = read_data.prepare_FAS_graph(fas_dict[geneid], isoform_dict[geneid])
        mov_table, mov_fig = read_data.prepare_movement(movement, geneid)
        return ('Selected {}'.format(geneid), mov_table, mov_fig, exp_data[geneid]['table'],
                {samples[0]: exp_data[geneid][samples[0]], samples[1]: exp_data[geneid][samples[1]]}, fas_graph)
    else:
        return None, m_mock0, m_mock1, exp_mock, {samples[0]: 0.0, samples[1]: 0.0}, fas_mock


# Expression callbacks
@app.callback(
    Output(exp_graph, 'figure'),
    Input(exp_dropdown, 'value'),
    Input(exp_store, 'data'),
    Input(exp_store2, 'data'),
    State(gene_input, 'value')
)
def generate_chart(condition, exp_data, exp_gene, gid):
    dff = pd.DataFrame(exp_data)
    if not condition == 'both':
        dff = dff = dff[dff['Condition'] == condition]
    fig = px.box(dff, x='transcriptid', y='expression', color='Condition', range_y=[0, 1])
    fig.update_layout(title_text=gid)
    return fig


@app.callback(
    Output(exp_table, 'data'),
    Input(exp_store, 'data')
)
def generate_chart(exp_data):
    return exp_data


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
                showlegend=False,
                title_text=gid
            )
    return fig


@app.callback(
    Output(mov_graph2, 'figure'),
    Output(mov_table, 'data'),
    Input(mov_dropdown, 'value'),
    Input(mov_store, 'data'),
)
def generate_boxplot(plottype, mov_table):
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
    return figure, mov_table


def main():

    app.run_server(debug=True)


if __name__ == '__main__':
    main()

