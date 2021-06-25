#import dash
#import dash_jsme

import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
import numpy as np
import pandas as pd
#import dash_daq as daq
#
import plotly.express as px
import plotly.graph_objects as go
#
#import time
from app import app
#from dash import no_update
from dash.dependencies import Input, Output, State
#from jupyter_dash import JupyterDash

#from index import similarity_df
import visdcc
#import layout_results
from organtox_backend import *
#

############################################################ Results Layout ############################################################

layout = html.Div([

    html.Div(id='similarity_df'),
   
#### OrganTox Profile #### 
    
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
        

        html.Div(children=[
            visdcc.Run_js('javascript'),
            dbc.Button('EVALUATE NEW COMPOUND', id='refresh-page', className='mr-1', color='info', style={'width': '580px'})
        ], style={'margin-bottom':'1vw'}),
          
        html.Br(),
    ])
])

############################################################ Callbacks ############################################################

#### OrganTox Profile ####
@app.callback(
    Output('graph-profile', 'figure'),
    [Input('data-source-profile', 'value'), Input('descriptor-type-profile', 'value'), Input('slider-profile', 'value'), Input('similarity_df', 'children')])
def update_figure(source_profile, descriptor_profile, slider_value_profile, json_df):
    try:
        df = pd.read_json(json_df, orient='split')
        df_0 = df[df.source.isin(source_profile)]
        df_1 = df_0[df_0.descriptor == descriptor_profile]
        df_2 = df_1[df_1.time == slider_value_profile]
        
        df_3 = df_2.groupby(['compound_name', 'parameter'], as_index=False).agg({'outcome': 'max', 'similarity': 'median'})
        df_4 = df_3.sort_values(by='similarity', ascending=False).head(1000)
         
        df_pivot = df_4.groupby(['compound_name', 'parameter'], sort=False, as_index=True)['outcome'].median().unstack()
        sort_list = [k for k in param_list if k in df_pivot.columns.values]
        df_pivot = df_pivot[sort_list]
        
        df_to_plotly = {'x': df_pivot.columns.tolist(), 'y': df_pivot.index.tolist(), 'z': df_pivot.values}
    
        fig = go.Figure(data=go.Heatmap(df_to_plotly,
                                        colorscale = 'Cividis_r',
                                        hovertemplate='Parameter: %{x}<br>Compound: %{y}<br>Effect: %{z}<extra></extra>'))
        
        fig.update_layout(height=700,
                          showlegend=False,
                          transition_duration=500,
                          xaxis=dict(
                              tickangle=30,
                              tickmode='array',
                              tickvals=[x for x in sorted(df_pivot.columns.tolist())]),
                          margin = dict(t=30, b=30),
                          #coloraxis_colorbar = dict(
                          #    thicknessmode="pixels", thickness=50,
                          #    yanchor="top",
                          #    ticks="outside"),
                          yaxis=dict(
                              title='Nearest Neighbors',
                              tickmode='array',
                              tickvals=[y for y in sorted(df_pivot.index.tolist())]),
                          yaxis_nticks=len(df_pivot.index.tolist())
                         )
        
        fig.update_traces(colorbar_title='Effect',
                          colorbar_tickmode='array',
                          colorbar_tickvals=[1, 0],
                          colorbar_ticktext=['Observed', 'Non-observed'],
                          selector=dict(type='heatmap'))
        
        return fig
    
    except:
        fig = go.Figure().add_annotation(text='No data to display', font=dict(color='#EAEAEA', size=25),
                                         showarrow=False)
        
        fig.update_layout(height=700,
                         xaxis={'showgrid': False, 'zeroline': False, 'showticklabels':False},
                         yaxis={'showgrid': False, 'zeroline': False, 'showticklabels':False},
                         margin=dict(t=30,b=30),
                         plot_bgcolor='white')
        return fig


    #### Tissue Profile ####
@app.callback(
    Output('graph-weight', 'figure'),
    [Input('data-source-tissue', 'value'), Input('descriptor-type-tissue', 'value'), Input('tissue-type', 'value'), Input('slider-similarity-tissue', 'value'), Input('similarity_df', 'children')])
def update_figure(source_tissue, descriptor_tissue, tissue, slider_value_tissue, json_df):
    try:
        df = pd.read_json(json_df, orient='split')
        df_0 = df[df.parameter_type == 'Weight']
        df_1 = df_0[df_0.source.isin(source_tissue)]
        df_2 = df_1[df_1.descriptor == descriptor_tissue]
        df_3 = df_2[df_2.parameter == tissue]
        df_4 = df_3[df_3.similarity.between(slider_value_tissue[0], slider_value_tissue[1])]
        df_5 = df_4.round(2)
    
        fig = px.scatter(df_5,
                         x='time',
                         y='value',
                         color='similarity',
                         color_continuous_scale='Inferno_r',
                         symbol=f'compound_name',
                         range_color=[slider_value_tissue[0], slider_value_tissue[1]],
                         hover_data=['value', 'value_unit', 'dose'], #control_val
                         labels={'compound_name': 'Compound Name',
                                 'time': 'Days',
                                 'dose': 'Dose (mg/kg/day)',
                                 'value': f'Experiment Value',
                                 'value_unit': 'Value Unit',
                                 'similarity': 'Tanimoto'}
                        )
    
        fig.update_layout(showlegend=False,
                          transition_duration=500,
                          xaxis = dict(
                              tickmode = 'array',
                              tickvals = [x for x in sorted(df_4['time'].unique())]),
                          margin=dict(t=30,b=30)
                         )
        
        fig.update_coloraxes(colorbar_tickmode='array',
                             colorbar_tickvals=np.arange(0,1.1,0.1),
                             colorbar_ticktext=['0.0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1.0']
                            )
    
        return fig
    except:
        fig = go.Figure().add_annotation(x=2.5, y=1.5, text='No data to display', font=dict(color='#EAEAEA', size=25),
                                         showarrow=False)
        
        fig.update_layout(xaxis={'showgrid': False, 'zeroline': False, 'showticklabels':False},
                          yaxis={'showgrid': False, 'zeroline': False, 'showticklabels':False},
                          margin=dict(t=30,b=30),
                          plot_bgcolor='white')
        return fig

    #### Tissue MuDRA plot ####
@app.callback(
    Output('graph-mudra-tissue', 'figure'),
    [Input('data-source-tissue', 'value'), Input('tissue-type', 'value'), Input('slider-mudra-tissue', 'value'), Input('similarity_df', 'children')])
def update_mudra(source, tissue, slider_value, json_df):
    try:
        df = pd.read_json(json_df, orient='split')
        df0 = df[df.parameter_type == 'Weight']
        df1 = df0[df0.source.isin(source)]
        df2 = df1[df1.parameter == tissue]
        
        df_morgan = df2[df2.descriptor == 'morgan'].groupby('InChIKey').first().sort_values(by='similarity', ascending=False).head(slider_value)#.sample(frac=1)
        df_maccs = df2[df2.descriptor == 'maccs'].groupby('InChIKey').first().sort_values(by='similarity', ascending=False).head(slider_value)#.sample(frac=1)
        df_atom_pair = df2[df2.descriptor == 'atom_pair'].groupby('InChIKey').first().sort_values(by='similarity', ascending=False).head(slider_value)#.sample(frac=1)
        df_avalon = df2[df2.descriptor == 'avalon'].groupby('InChIKey').first().sort_values(by='similarity', ascending=False).head(slider_value)#.sample(frac=1)

        fig = go.Figure()
        fig.add_trace(go.Scatterpolargl(r = df_morgan.similarity,
                                        theta = np.arange(92, 179, 87/df_morgan.shape[0]),
                                        name = 'Morgan',
                                        marker=dict(size=10, color='mediumseagreen'),
                                        customdata=df_morgan['compound_name']
                                       )
                     )
        
        fig.add_trace(go.Scatterpolargl(r= df_maccs.similarity,
                                        theta = np.arange(2, 89, 87/df_maccs.shape[0]),
                                        name = 'MACCS',
                                        marker=dict(size=10, color='darkorange'),
                                        customdata=df_maccs['compound_name']
                                       )
                     )

        fig.add_trace(go.Scatterpolargl(r = df_atom_pair.similarity,
                                        theta = np.arange(182, 269, 87/df_avalon.shape[0]),
                                        name = 'AtomPair',
                                        marker=dict(size=10, color='mediumpurple'),
                                        customdata=df_avalon['compound_name']
                                       )
                     )
         
        fig.add_trace(go.Scatterpolargl(r = df_avalon.similarity,
                                        theta = np.arange(272, 359, 87/df_atom_pair.shape[0]),
                                        name = 'Avalon',
                                        marker=dict(size=10, color = 'magenta'),
                                        customdata=df_atom_pair['compound_name']
                                       )
                     )
        
        fig.update_traces(mode='markers',
                          marker=dict(line_color='white', opacity=0.7),
                          hovertemplate='<br>'.join(['Compound Name: %{customdata}',
                                                     'Similarity: %{r}',
                                                    ])
                         )
        
        fig.update_layout(
            transition_duration=500,
            font_size = 16,
            showlegend = True,
            polar = dict(angularaxis = dict(showticklabels=True,
                                            showgrid=True,
                                            ticks='',
                                            direction = 'counterclockwise',
                                            tickmode='array',
                                            tickvals=[0, 90, 180, 270],
                                            ticktext=['','', '', ''],
                                            linewidth=1,
                                            showline=True,
                                            linecolor='black'
                                           ),
                         radialaxis = dict(tickfont_size=12,
                                           tickangle=30,
                                           ticks='outside',
                                           showline=True,
                                           range=[0, 1]
                                          )                         
                        ), margin=dict(t=30,b=30),
            
        )

    
        return fig
    except:
        fig = go.Figure().add_annotation(x=2.5, y=1.5, text='No data to display', font=dict(color='#EAEAEA', size=25),
                                         showarrow=False)
        
        fig.update_layout(xaxis={'showgrid': False, 'zeroline': False, 'showticklabels':False},
                          yaxis={'showgrid': False, 'zeroline': False, 'showticklabels':False},
                          margin=dict(t=30,b=30),
                          plot_bgcolor='white')
        return fig

    
    #### Clinical Hemistry Profile ####
@app.callback(
    Output('graph-clin-chem', 'figure'),
    [Input('data-source-clin-chem', 'value'), Input('descriptor-type-clin-chem', 'value'),
     Input('clin-chem-type', 'value'), Input('slider-similarity-clin-chem', 'value'),
     Input('similarity_df', 'children')])
def update_figure(source_tissue, descriptor_tissue, tissue, slider_value_tissue, json_df):
    try:
        df = pd.read_json(json_df, orient='split')
        df_0 = df[df.parameter_type == 'Clinical chemistry']
        df_1 = df_0[df_0.source.isin(source_tissue)]
        df_2 = df_1[df_1.descriptor == descriptor_tissue]
        df_3 = df_2[df_2.parameter == tissue]
        df_4 = df_3[df_3.similarity.between(slider_value_tissue[0], slider_value_tissue[1])]
    
        fig = px.scatter(df_4,
                         x='time',
                         y='value',
                         color='similarity',
                         color_continuous_scale='Inferno_r',
                         symbol=f'compound_name',
                         range_color=[slider_value_tissue[0], slider_value_tissue[1]],
                         hover_data=['value', 'value_unit', 'dose'],
                         labels={'compound_name': 'Compound Name',
                                 'time': 'Days',
                                 'dose': 'Dose (mg/kg/day)',
                                 'value': 'Experiment Value',
                                 'value_unit': 'Value Unit',
                                 'similarity': 'Tanimoto'}
                        )
    
        fig.update_layout(showlegend=False,
                          transition_duration=500,
                          xaxis = dict(
                              tickmode = 'array',
                              tickvals = [x for x in sorted(df_4['time'].unique())]),
                          margin=dict(t=30,b=30)
                         )
        
        fig.update_coloraxes(colorbar_tickmode='array',
                             colorbar_tickvals=np.arange(0,1.1,0.1),
                             colorbar_ticktext=['0.0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1.0']
                            )
    
        return fig
    except:
        fig = go.Figure().add_annotation(x=2.5, y=1.5, text='No data to display', font=dict(color='#EAEAEA', size=25),
                                         showarrow=False)
        
        fig.update_layout(xaxis={'showgrid': False, 'zeroline': False, 'showticklabels':False},
                          yaxis={'showgrid': False, 'zeroline': False, 'showticklabels':False},
                          margin=dict(t=30,b=30),
                          plot_bgcolor='white')
        return fig
    

    #### Clinical Chemistry MuDRA plot ####
@app.callback(
    Output('graph-mudra-clin-chem', 'figure'),
    [Input('data-source-clin-chem', 'value'), Input('clin-chem-type', 'value'),
     Input('slider-mudra-clin-chem', 'value'), Input('similarity_df', 'children')])
def update_mudra(source, clinchem, slider_value, json_df):
    try:
        df = pd.read_json(json_df, orient='split')
        df0 = df[df.parameter_type == 'Clinical chemistry']
        df1 = df0[df0.source.isin(source)]
        df2 = df1[df1.parameter == clinchem]
        
        df_morgan = df2[df2.descriptor == 'morgan'].groupby('InChIKey').first().sort_values(by='similarity', ascending=False).head(slider_value)#.sample(frac=1)
        df_maccs = df2[df2.descriptor == 'maccs'].groupby('InChIKey').first().sort_values(by='similarity', ascending=False).head(slider_value)#.sample(frac=1)
        df_atom_pair = df2[df2.descriptor == 'atom_pair'].groupby('InChIKey').first().sort_values(by='similarity', ascending=False).head(slider_value)#.sample(frac=1)
        df_avalon = df2[df2.descriptor == 'avalon'].groupby('InChIKey').first().sort_values(by='similarity', ascending=False).head(slider_value)#.sample(frac=1)

        fig = go.Figure()
        
        fig.add_trace(go.Scatterpolargl(r = df_morgan.similarity,
                                        theta = np.arange(92, 179, 87/df_morgan.shape[0]),
                                        name = 'Morgan',
                                        marker=dict(size=10, color='mediumseagreen'),
                                        customdata=df_morgan['compound_name']
                                       )
                     )
        
        fig.add_trace(go.Scatterpolargl(r= df_maccs.similarity,
                                        theta = np.arange(2, 89, 87/df_maccs.shape[0]),
                                        name = 'MACCS',
                                        marker=dict(size=10, color='darkorange'),
                                        customdata=df_maccs['compound_name']
                                       )
                     )

        fig.add_trace(go.Scatterpolargl(r = df_atom_pair.similarity,
                                        theta = np.arange(182, 269, 87/df_avalon.shape[0]),
                                        name = 'AtomPair',
                                        marker=dict(size=10, color='mediumpurple'),
                                        customdata=df_avalon['compound_name']
                                       )
                     )
         
        fig.add_trace(go.Scatterpolargl(r = df_avalon.similarity,
                                        theta = np.arange(272, 359, 87/df_atom_pair.shape[0]),
                                        name = 'Avalon',
                                        marker=dict(size=10, color = 'magenta'),
                                        customdata=df_atom_pair['compound_name']
                                       )
                     )
        
        fig.update_traces(mode='markers',
                          marker=dict(line_color='white', opacity=0.7),
                          hovertemplate='<br>'.join(['Compound Name: %{customdata}',
                                                     'Similarity: %{r}',
                                                    ])
                         )
        
        fig.update_layout(
            transition_duration=500,
            font_size = 16,
            showlegend = True,
            polar = dict(angularaxis = dict(showticklabels=True,
                                            showgrid=True,
                                            ticks='',
                                            direction = 'counterclockwise',
                                            tickmode='array',
                                            tickvals=[0, 90, 180, 270],
                                            ticktext=['','', '', ''],
                                            linewidth=1,
                                            showline=True,
                                            linecolor='black'
                                           ),
                         radialaxis = dict(tickfont_size=12,
                                           tickangle=30,
                                           ticks='outside',
                                           showline=True,
                                           range=[0, 1]
                                          )                         
                        ), margin=dict(t=30,b=30),
            
        )

    
        return fig
    except:
        fig = go.Figure().add_annotation(x=2.5, y=1.5, text='No data to display', font=dict(color='#EAEAEA', size=25),
                                         showarrow=False)
        
        fig.update_layout(xaxis={'showgrid': False, 'zeroline': False, 'showticklabels':False},
                          yaxis={'showgrid': False, 'zeroline': False, 'showticklabels':False},
                          margin=dict(t=30,b=30),
                          plot_bgcolor='white')
        return fig
    
    
        #### Hematology Profile ####
@app.callback(
    Output('graph-hemato', 'figure'),
    [Input('data-source-hemato', 'value'), Input('descriptor-type-hemato', 'value'), Input('hemato-type', 'value'),
     Input('slider-similarity-hemato', 'value'), Input('similarity_df', 'children')])
def update_figure(source_tissue, descriptor_tissue, tissue, slider_value_tissue, json_df):
    try: 
        df = pd.read_json(json_df, orient='split')
        df_0 = df[df.parameter_type == 'Hematology']
        df_1 = df_0[df_0.source.isin(source_tissue)]
        df_2 = df_1[df_1.descriptor == descriptor_tissue]
        df_3 = df_2[df_2.parameter == tissue]
        df_4 = df_3[df_3.similarity.between(slider_value_tissue[0], slider_value_tissue[1])]
    
        fig = px.scatter(df_4,
                         x='time',
                         y='value',
                         color='similarity',
                         color_continuous_scale='Inferno_r',
                         symbol=f'compound_name',
                         range_color=[slider_value_tissue[0], slider_value_tissue[1]],
                         hover_data=['value', 'value_unit', 'dose'],
                         labels={'compound_name': 'Compound Name',
                                 'time': 'Days',
                                 'dose': 'Dose (mg/kg/day)',
                                 'value': 'Experiment Value',
                                 'value_unit': 'Value Unit',
                                 'similarity': 'Tanimoto'}
                        )
    
        fig.update_layout(showlegend=False,
                          transition_duration=500,
                          xaxis = dict(
                              tickmode = 'array',
                              tickvals = [x for x in sorted(df_4['time'].unique())]),
                          margin=dict(t=30,b=30)
                         )
        
        fig.update_coloraxes(colorbar_tickmode='array',
                             colorbar_tickvals=np.arange(0,1.1,0.1),
                             colorbar_ticktext=['0.0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1.0']
                            )
    
        return fig
    except:
        fig = go.Figure().add_annotation(x=2.5, y=1.5, text='No data to display', font=dict(color='#EAEAEA', size=25),
                                         showarrow=False)
        
        fig.update_layout(xaxis={'showgrid': False, 'zeroline': False, 'showticklabels':False},
                          yaxis={'showgrid': False, 'zeroline': False, 'showticklabels':False},
                          margin=dict(t=30,b=30),
                          plot_bgcolor='white')
        return fig
    
    
    
        #### Hematology MuDRA plot ####
@app.callback(
    Output('graph-mudra-hemato', 'figure'),
    [Input('data-source-hemato', 'value'), Input('hemato-type', 'value'), Input('slider-mudra-hemato', 'value'), Input('similarity_df', 'children')])
def update_mudra(source, hemato, slider_value, json_df):
    try:
        df = pd.read_json(json_df, orient='split')
        df0 = df[df.parameter_type == 'Hematology']
        df1 = df0[df0.source.isin(source)]
        df2 = df1[df1.parameter == hemato]
        
        df_morgan = df2[df2.descriptor == 'morgan'].groupby('InChIKey').first().sort_values(by='similarity', ascending=False).head(slider_value)#.sample(frac=1)
        df_maccs = df2[df2.descriptor == 'maccs'].groupby('InChIKey').first().sort_values(by='similarity', ascending=False).head(slider_value)#.sample(frac=1)
        df_atom_pair = df2[df2.descriptor == 'atom_pair'].groupby('InChIKey').first().sort_values(by='similarity', ascending=False).head(slider_value)#.sample(frac=1)
        df_avalon = df2[df2.descriptor == 'avalon'].groupby('InChIKey').first().sort_values(by='similarity', ascending=False).head(slider_value)#.sample(frac=1)

        fig = go.Figure()
        fig.add_trace(go.Scatterpolargl(r = df_morgan.similarity,
                                        theta = np.arange(92, 179, 87/df_morgan.shape[0]),
                                        name = 'Morgan',
                                        marker=dict(size=10, color='mediumseagreen'),
                                        customdata=df_morgan['compound_name']
                                       )
                     )
        
        fig.add_trace(go.Scatterpolargl(r= df_maccs.similarity,
                                        theta = np.arange(2, 89, 87/df_maccs.shape[0]),
                                        name = 'MACCS',
                                        marker=dict(size=10, color='darkorange'),
                                        customdata=df_maccs['compound_name']
                                       )
                     )

        fig.add_trace(go.Scatterpolargl(r = df_atom_pair.similarity,
                                        theta = np.arange(182, 269, 87/df_avalon.shape[0]),
                                        name = 'AtomPair',
                                        marker=dict(size=10, color='mediumpurple'),
                                        customdata=df_avalon['compound_name']
                                       )
                     )
         
        fig.add_trace(go.Scatterpolargl(r = df_avalon.similarity,
                                        theta = np.arange(272, 359, 87/df_atom_pair.shape[0]),
                                        name = 'Avalon',
                                        marker=dict(size=10, color = 'magenta'),
                                        customdata=df_atom_pair['compound_name']
                                       )
                     )
        
        fig.update_traces(mode='markers',
                          marker=dict(line_color='white', opacity=0.7),
                          hovertemplate='<br>'.join(['Compound Name: %{customdata}',
                                                     'Similarity: %{r}',
                                                    ])
                         )
        
        fig.update_layout(
            transition_duration=500,
            font_size = 16,
            showlegend = True,
            polar = dict(angularaxis = dict(showticklabels=True,
                                            showgrid=True,
                                            ticks='',
                                            direction = 'counterclockwise',
                                            tickmode='array',
                                            tickvals=[0, 90, 180, 270],
                                            ticktext=['','', '', ''],
                                            linewidth=1,
                                            showline=True,
                                            linecolor='black'
                                           ),
                         radialaxis = dict(tickfont_size=12,
                                           tickangle=30,
                                           ticks='outside',
                                           showline=True,
                                           range=[0, 1]
                                          )                         
                        ), margin=dict(t=30,b=30),
            
        )
    
        return fig
    except:
        fig = go.Figure().add_annotation(x=2.5, y=1.5, text='No data to display', font=dict(color='#EAEAEA', size=25),
                                         showarrow=False)
        
        fig.update_layout(xaxis={'showgrid': False, 'zeroline': False, 'showticklabels':False},
                          yaxis={'showgrid': False, 'zeroline': False, 'showticklabels':False},
                          margin=dict(t=30,b=30),
                          plot_bgcolor='white')
        return fig

#### Refresh Page ####
@app.callback(
    Output('javascript', 'run'),
    [Input('refresh-page', 'n_clicks')])
def refresh_button_clicked(data):
    if data is not None:
        #fetch_new_data_function()
        return "location.reload();"
    return ''