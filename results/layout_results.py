import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc

from organtox_backend import *
import callbacks

################################################### Results Layout ###################################################

layout = html.Div([

#### OrganTox Profile ####
    dbc.Container(
    html.Center(children=[
        html.Br(),
        html.Br(),
        html.Hr(className="my-2", style={'width': '80vw'}),
        html.Br(),

# OrganTox Options #
    html.Center(children=[
        html.Div(children=[
            html.H3('OrganTox Profile')
            ], style={'display': 'block'}),
        html.Div(children=[
            html.Label('Data Source', style={'font-weight':'bold', 'display': 'inline-block'}),
            dbc.FormGroup([
                dbc.Checklist(
                    options=[
                        {'label': 'DrugMatrix', 'value': 'DrugMatrix'},
                        {'label': 'TG-GATEs', 'value': 'TG-GATEs'}
                    ], value=['DrugMatrix', 'TG-GATEs'], id='data-source-profile', inline=True)
                ],  style={'verticalAlign': 'middle', 'textAlign': 'left',
                           'padding': '5px 0px 0px 0px', 'display': 'block'})
        ], style={'verticalAlign': 'top', 'textAlign': 'left',
                  'margin-bottom': '1vw', 'display': 'inline-block'}),


        html.Div(children=[
            html.Label('Descriptor Type', style={'font-weight':'bold'}),
            dcc.Dropdown(
                id='descriptor-type-profile',
                options=[{'label': 'Morgan', 'value': 'morgan'},
                         {'label': 'MACCS', 'value': 'maccs'},
                         {'label': 'Avalon', 'value': 'avalon'},
                         {'label': 'Atom Pair', 'value': 'atom_pair'}],
                value='morgan')
            ], style={'verticalAlign': 'top', 'textAlign': 'left',
                          'width': '286px', 'padding': '0px 10px 0px 0px', 'margin-bottom': '1vw', 'display': 'inline-block'})
        ]),

# OrganTox graphs #

    html.Center([
        html.Div([
        dbc.Row([
            dbc.Col(
                    html.Div([
                        html.Center([
                            html.Div(children=[
                                html.Div(children=[
                                    dcc.Graph(id='graph-profile', config={'responsive': True})
                                ], style={'border':'1px solid', 'border-radius': 4, 'border-color': '#F4F6F6',
                                          'display': 'block',
                                          'height': 700}),
                                html.Div(children=[
                                    html.H6('Time (days):', style={'textAlign':'right'}),
                                    dcc.Slider(
                                        id='slider-profile',
                                        min=2,
                                        max=8,
                                        step=1,
                                        value=5,
                                        included=False,
                                        marks={2: {'label': '2'},
                                               3: {'label': '3'},
                                               4: {'label': '4'},
                                               5: {'label': '5'},
                                               6: {'label': '6'},
                                               7: {'label': '7'},
                                               8: {'label': '8'}})
                                ], style={'display':'grid', 'padding': '0px 0px 0px 0px', 'grid-template-columns': '9% 91%'})
                            ])
                        ])
                    ]), lg=10)
            ], justify='center')
        ])
    ]),

    html.Br(),
    html.Br(),


#### Tissue Profile ####

    html.Hr(className="my-2", style={'width': '80vw'}),
    html.Br(),

    # Tissue Options #
    html.Center([
        html.Div(children=[
            html.H3('Tissue Profile')
            ], style={'display': 'block'}),
        html.Div(children=[
            html.Label('Data Source',  style={'font-weight':'bold', 'display': 'inline-block'}),
            dbc.FormGroup([
                dbc.Checklist(
                    options=[
                        {'label': 'DrugMatrix', 'value': 'DrugMatrix'},
                        {'label': 'TG-GATEs', 'value': 'TG-GATEs'}
                    ], value=['DrugMatrix', 'TG-GATEs'], id='data-source-tissue', inline=True)
                ],  style={'verticalAlign': 'middle', 'textAlign': 'left',
                           'padding': '5px 0px 0px 0px', 'display': 'block'})
            ], style={'verticalAlign': 'top', 'textAlign': 'left',
                      'margin-bottom': '1vw', 'display': 'inline-block'}),

        html.Div(children=[
              html.Label('Descriptor Type', style={'font-weight':'bold'}),
              dcc.Dropdown(
                  id='descriptor-type-tissue',
                  options=[{'label': 'Morgan', 'value': 'morgan'},
                           {'label': 'MACCS', 'value': 'maccs'},
                           {'label': 'Avalon', 'value': 'avalon'},
                           {'label': 'Atom Pair', 'value': 'atom_pair'}],
                  value='morgan')
              ], style={'verticalAlign': 'top', 'textAlign': 'left',
                        'width': '286px', 'padding': '0px 10px 0px 0px', 'margin-bottom': '1vw', 'display': 'inline-block'}),

          html.Div(children=[
              html.Label('Tissues', style={'font-weight':'bold'}),
              dcc.Dropdown(
                  id='tissue-type',
                  options=[{'label': i, 'value': i} for i in sorted(weight_list)],
                  value='Body weight')
              ], style={'verticalAlign': 'top', 'textAlign': 'left',
                        'width': '286px', 'padding': '0px 0px 0px 0px', 'margin-bottom': '1vw', 'display': 'inline-block'})

    ]),

    # Tissue Profile Graphs #

        html.Center([
            dbc.Row([
            # Weight Effect #
                dbc.Col(
                    html.Div([
                        html.Center([
                            html.Div(children=[
                                html.Div(children=[
                                    html.H5('Chemical Similarity Map of Weight Change')
                                ], style={'display': 'block'}),
                                html.Div(children=[
                                    dcc.Graph(id='graph-weight')
                                ], style={'border':'1px solid', 'border-radius': 4, 'border-color': '#F4F6F6',
                                          'display': 'block'}),
                                html.Div(children=[
                                    html.H6('Similarity: ', style={'textAlign':'right'}),
                                    dcc.RangeSlider(
                                        id='slider-similarity-tissue',
                                        min=0,
                                        max=1,
                                        step=0.1,
                                        value=[0,1],
                                        pushable=0.1,
                                        allowCross=False,
                                        marks={0: {'label': 'Low'},
                                               0.1: {'label': '0.1'},
                                               0.2: {'label': '0.2'},
                                               0.3: {'label': '0.3'},
                                               0.4: {'label': '0.4'},
                                               0.5: {'label': '0.5'},
                                               0.6: {'label': '0.6'},
                                               0.7: {'label': '0.7'},
                                               0.8: {'label': '0.8'},
                                               0.9: {'label': '0.9'},
                                               1: {'label': 'High'}})
                                ], style={'display':'grid', 'grid-template-columns': '10% 90%'})
                            ])
                        ])
                    ]), lg=5),

                # MuDRA #
                dbc.Col(
                    html.Div([
                        html.Center([
                            html.Div(children=[
                                html.Div(children=[
                                    html.H5('Multi-Descriptor Read Across')
                                ], style={'display': 'block'}),
                                html.Div(children=[
                                    dcc.Graph(id='graph-mudra-tissue')
                                ], style={'border':'1px solid', 'border-radius': 4, 'border-color': '#F4F6F6',
                                          'display': 'block'}),
                                html.Div(children=[
                                    html.H6('Neighbors: ', style={'textAlign':'left'}),
                                    dcc.Slider(
                                        id='slider-mudra-tissue',
                                        min=0,
                                        max=100,
                                        step=1,
                                        value=10,
                                        marks={0: {'label': 'Min'},
                                               100: {'label': 'Max'}}
                                    )
                                ], style={'display':'grid', 'grid-template-columns': '12% 88%'})
                            ])
                        ])
                    ]), lg=5),

                ], justify='center')
        ]),

    html.Br(),
    html.Br(),


    #### Clinical Chemistry Profile ####
    html.Hr(className="my-2", style={'width': '80vw'}),
    html.Br(),

    # Options #
    html.Center([
        html.Div(children=[
            html.H3('Clinical Chemistry Profile')
            ], style={'display': 'block'}),
        html.Div(children=[
            html.Label('Data Source',  style={'font-weight':'bold', 'display': 'inline-block'}),
            dbc.FormGroup([
                dbc.Checklist(
                    options=[
                        {'label': 'DrugMatrix', 'value': 'DrugMatrix'},
                        {'label': 'TG-GATEs', 'value': 'TG-GATEs'}
                    ], value=['DrugMatrix', 'TG-GATEs'], id='data-source-clin-chem', inline=True)
                ],  style={'verticalAlign': 'middle', 'textAlign': 'left',
                           'padding': '5px 0px 0px 0px', 'display': 'block'})
            ], style={'verticalAlign': 'top', 'textAlign': 'left',
                      'margin-bottom': '1vw', 'display': 'inline-block'}),

        html.Div(children=[
              html.Label('Descriptor Type', style={'font-weight':'bold'}),
              dcc.Dropdown(
                  id='descriptor-type-clin-chem',
                  options=[{'label': 'Morgan', 'value': 'morgan'},
                           {'label': 'MACCS', 'value': 'maccs'},
                           {'label': 'Avalon', 'value': 'avalon'},
                           {'label': 'Atom Pair', 'value': 'atom_pair'}],
                  value='morgan')
              ], style={'verticalAlign': 'top', 'textAlign': 'left',
                        'width': '286px', 'padding': '0px 10px 0px 0px', 'margin-bottom': '1vw', 'display': 'inline-block'}),

          html.Div(children=[
              html.Label('Clinical Chemistry', style={'font-weight':'bold'}),
              dcc.Dropdown(
                  id='clin-chem-type',
                  options=[{'label': i, 'value': i} for i in sorted(clin_chem_list)],
                  value='Creatinine')
              ], style={'verticalAlign': 'top', 'textAlign': 'left',
                        'width': '286px', 'padding': '0px 0px 0px 0px', 'margin-bottom': '1vw', 'display': 'inline-block'})
    ]),

        # Tissue Profile Graphs #
        html.Center([
            dbc.Row([
            # Clinical Chemistry Effect #
                dbc.Col(
                    html.Div([
                        html.Center([
                            html.Div(children=[
                                html.H5('Chemical Similarity Map of Clinical Chemistry')
                            ], style={'display': 'block'}),
                            html.Div(children=[
                                dcc.Graph(id='graph-clin-chem')
                            ], style={'border':'1px solid', 'border-radius': 4, 'border-color': '#F4F6F6',
                                          'display': 'block'}),
                            html.Div(children=[
                                html.H6('Similarity: ', style={'textAlign':'right'}),
                                    dcc.RangeSlider(
                                        id='slider-similarity-clin-chem',
                                        min=0,
                                        max=1,
                                        step=0.1,
                                        value=[0,1],
                                        pushable=0.1,
                                        allowCross=False,
                                        marks={0: {'label': 'Low'},
                                               0.1: {'label': '0.1'},
                                               0.2: {'label': '0.2'},
                                               0.3: {'label': '0.3'},
                                               0.4: {'label': '0.4'},
                                               0.5: {'label': '0.5'},
                                               0.6: {'label': '0.6'},
                                               0.7: {'label': '0.7'},
                                               0.8: {'label': '0.8'},
                                               0.9: {'label': '0.9'},
                                               1: {'label': 'High'}})
                            ], style={'display':'grid', 'grid-template-columns': '10% 90%'})
                        ])
                    ]), lg=5),


            # MuDRA #
                dbc.Col(
                    html.Div([
                        html.Center([
                            html.Div(children=[
                                html.Div(children=[
                                    html.H5('Multi-Descriptor Read Across')
                                ], style={'display': 'block'}),
                                html.Div(children=[
                                    dcc.Graph(id='graph-mudra-clin-chem')
                                ], style={'border':'1px solid', 'border-radius': 4, 'border-color': '#F4F6F6',
                                          'display': 'block'}),
                                html.Div(children=[
                                    html.H6('Neighbors: ', style={'textAlign':'left'}),
                                    dcc.Slider(
                                        id='slider-mudra-clin-chem',
                                        min=0,
                                        max=100,
                                        step=1,
                                        value=10,
                                        marks={0: {'label': 'Min'},
                                               100: {'label': 'Max'}}
                                    )
                                ], style={'display':'grid', 'grid-template-columns': '12% 88%'})
                            ])
                        ])
                    ]), lg=5),
                    ], justify='center')
        ]),

    html.Br(),
    html.Br(),


        #### Hematology Profile ####
    html.Hr(className="my-2", style={'width': '80vw'}),

    html.Br(),

    # Options #
    html.Center([
        html.Div(children=[
            html.H3('Hematology Profile')
            ], style={'display': 'block'}),
        html.Div(children=[
            html.Label('Data Source',  style={'font-weight':'bold', 'display': 'inline-block'}),
            dbc.FormGroup([
                dbc.Checklist(
                    options=[
                        {'label': 'DrugMatrix', 'value': 'DrugMatrix'},
                        {'label': 'TG-GATEs', 'value': 'TG-GATEs'}
                    ], value=['DrugMatrix', 'TG-GATEs'], id='data-source-hemato', inline=True)
                ],  style={'verticalAlign': 'middle', 'textAlign': 'left',
                           'padding': '5px 0px 0px 0px', 'display': 'block'})
            ], style={'verticalAlign': 'top', 'textAlign': 'left',
                      'margin-bottom': '1vw', 'display': 'inline-block'}),

        html.Div(children=[
              html.Label('Descriptor Type', style={'font-weight':'bold'}),
              dcc.Dropdown(
                  id='descriptor-type-hemato',
                  options=[{'label': 'Morgan', 'value': 'morgan'},
                           {'label': 'MACCS', 'value': 'maccs'},
                           {'label': 'Avalon', 'value': 'avalon'},
                           {'label': 'Atom Pair', 'value': 'atom_pair'}],
                  value='morgan')
              ], style={'verticalAlign': 'top', 'textAlign': 'left',
                        'width': '286px', 'padding': '0px 10px 0px 0px', 'margin-bottom': '1vw', 'display': 'inline-block'}),

          html.Div(children=[
              html.Label('Hematology', style={'font-weight':'bold'}),
              dcc.Dropdown(
                  id='hemato-type',
                  options=[{'label': i, 'value': i} for i in sorted(hemato_list)],
                  value='Hemoglobin')
              ], style={'verticalAlign': 'top', 'textAlign': 'left',
                        'width': '286px', 'padding': '0px 0px 0px 0px', 'margin-bottom': '1vw', 'display': 'inline-block'})
    ]),

        # Hematology Graphs #
        html.Center([
            dbc.Row([
            # Hematology Effect #
                dbc.Col(
                    html.Div([
                        html.Center([
                            html.Div(children=[
                                html.H5('Chemical Similarity Map of Hematology')
                            ], style={'display': 'block'}),
                            html.Div(children=[
                                dcc.Graph(id='graph-hemato')
                            ], style={'border':'1px solid', 'border-radius': 4, 'border-color': '#F4F6F6',
                                          'display': 'block'}),
                            html.Div(children=[
                                html.H6('Similarity: ', style={'textAlign':'right'}),
                                    dcc.RangeSlider(
                                        id='slider-similarity-hemato',
                                        min=0,
                                        max=1,
                                        step=0.1,
                                        value=[0,1],
                                        pushable=0.1,
                                        allowCross=False,
                                        marks={0: {'label': 'Low'},
                                               0.1: {'label': '0.1'},
                                               0.2: {'label': '0.2'},
                                               0.3: {'label': '0.3'},
                                               0.4: {'label': '0.4'},
                                               0.5: {'label': '0.5'},
                                               0.6: {'label': '0.6'},
                                               0.7: {'label': '0.7'},
                                               0.8: {'label': '0.8'},
                                               0.9: {'label': '0.9'},
                                               1: {'label': 'High'}})
                            ], style={'display':'grid', 'grid-template-columns': '10% 90%'})
                        ])
                    ]), lg=5),


            # MuDRA #
                dbc.Col(
                    html.Div([
                        html.Center([
                            html.Div(children=[
                                html.Div(children=[
                                    html.H5('Multi-Descriptor Read Across')
                                ], style={'display': 'block'}),
                                html.Div(children=[
                                    dcc.Graph(id='graph-mudra-hemato')
                                ], style={'border':'1px solid', 'border-radius': 4, 'border-color': '#F4F6F6',
                                          'display': 'block'}),
                                html.Div(children=[
                                    html.H6('Neighbors: ', style={'textAlign':'left'}),
                                    dcc.Slider(
                                        id='slider-mudra-hemato',
                                        min=0,
                                        max=100,
                                        step=1,
                                        value=10,
                                        marks={0: {'label': 'Min'},
                                               100: {'label': 'Max'}}
                                    )
                                ], style={'display':'grid', 'grid-template-columns': '12% 88%'})
                            ])
                        ])
                    ]), lg=5),
                    ], justify='center')
        ]),
        html.Br(),
        html.Br(),

    ])
    , fluid=True)
])