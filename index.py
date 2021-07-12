import dash_jsme
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc

from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate

from app import app
from results import layout_results
from organtox_backend import *
import callbacks

############################################################# Index Layout #############################################################

layout_index = html.Div([
    html.Center(dbc.Jumbotron([html.H2('OrganTox', className='display-4'),
                               html.P('Organ Toxicity Web Application v. 0.1', className="lead"),
                               html.Hr(className="my-2"),
                               html.P('The OrganTox app is an user-friendly tool for profiling '
                                      'the subchronic toxicity of similar chemicals.')
                              ], style={'horizontalAlign':'middle','height': '260px', 'width': '90%','margin-top': '2px'})),

#### Molecule Editor ####
# Instructions Buttons #

    html.Center(children=[
        html.Div(children=[
                html.Div(children=[
                    dbc.Button('INSTRUCTIONS', id='open-instructions', className='mr-1', color='light', style={'width': '286px'}),
                    dbc.Button('APP CHARACTERISTICS', id='open-app-characheristics', color='light', style={'width': '286px'}),
                    dbc.Modal(children=[
                        dbc.ModalHeader('Instructions'),
                        dbc.ModalBody([dcc.Markdown('''
                                       ##### Insert SMILES

                                       Directly paste the [SMILES](https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system)
                                       representation of the desired chemical structure.

                                       ##### ...or Draw Chemical
                                       Draw the chemical structure using the [Molecular Editor](https://jsme-editor.github.io/)

                                       ##### Evaluate
                                       Hit the SUBMIT button for profiling similar chemicals.

                                       ##### Important
                                       Current version of this app can only handle small organic molecules.
                                       ''')                                      ]),
                        dbc.ModalFooter(dbc.Button('Close', id='close-instructions', className='ml-auto'))
                    ], id='modal-instructions'),
                    dbc.Modal(children=[
                        dbc.ModalHeader('APP CHARACTERISTICS'),
                        dbc.ModalBody([dcc.Markdown('''
                                       OrganTox is an open-source web application to facilitate the visualization and read-across assessment of
                                       affected organs in sub-acute and sub-chronic toxicity tests.
                                       The web app leverages curated public data on compounds tested orally in
                                       short-term in vivo assays (1-8 days) that contain body and organ weight, clinical chemistry, and hematology data.
                                       The data were collected from the [DrugMatrix](https://www.sciencedirect.com/science/article/abs/pii/S0168165605001513)
                                       and [TG-GATEs](https://academic.oup.com/nar/article/43/D1/D921/2439524) database.
                                       Curated data comprised 63,814 and 14,971 data points for
                                       1-8 days and 514 and 126 chemicals, respectively. The integration of clinical chemistry,
                                       hematology, organ and body weights provided a better indication of appropriate dose selection
                                       (i.e., identifying a maximum tolerated dose). The tool integrates all the available data implementing the
                                       [Multi-Descriptor Read-Across (MuDRA)](https://pubs.acs.org/doi/10.1021/acs.jcim.8b00124)
                                       approach for similarity analysis. Doing this allows the user to perform a structure-based search with their target
                                       chemical and then identify source chemicals to inform the potential target organ effects.
                                       The curated data included within the OrganTox platform is freely available for download.
                                       The proposed application is unique due to the composition and rigor in the underlying data collection and curation and
                                       the application of best practices in predictive toxicity model development.
                                       ''')
                                      ], style={'textAlign': 'justify'}),
                        dbc.ModalFooter(dbc.Button('Close', id='close-app-characteristics', className='ml-auto'))
                    ], id='modal-app-characteristics', scrollable=True)
                ])
        ]),

        html.Br(),


# JSME #

            html.Div(children=[
                html.Div(children=[
                    dbc.Input(id='smiles-input', placeholder='Insert SMILES', type='text',
                             style={'height':'110%', 'width': '580px', 'font-size':'110%', 'margin-bottom': '1vw'}),
                    html.Div(id='smiles-output', style={'display':'none'})
                ]),

                html.Div(children=[
                    dash_jsme.DashJsme(id='mol-input', label='my-label', height='400px', width='580px'), html.Div(id='mol-output', style={'display':'none'})
                    ], style={'textAlign': 'left', 'display': 'inline-block', 'margin-bottom': 'auto'}),

                html.Div(id='submit-output', style={'textAlign': 'center', 'display': 'inherit', 'margin-bottom': '1vw'}),

                html.Div(children=[
                    dbc.Button('SUBMIT', id='submit-button', className='mr-1', color='info', style={'width': '580px'})
                ], style={'margin-bottom':'1vw'}),

                html.Div(children=[
                    dbc.Spinner(html.Div(id='loading-output'),
                                color='secondary', size='sm', type='grow')
                ]),

                # Store the generated df with similarity values
                dcc.Store(id='similarity_df', storage_type='memory')

                ])
            ])
    ])

## App layout and validation
app.title = 'OrganTox'
app.layout = layout_index

############################################################# Callbacks #############################################################

#### Instructions & App Characteristics ####
def toggle_modal(n1, n2, is_open):
    if n1 or n2:
        return not is_open
    return is_open

app.callback(Output('modal-instructions', 'is_open'),
              [Input('open-instructions', 'n_clicks'), Input('close-instructions', 'n_clicks')],
              [State('modal-instructions', 'is_open')])(toggle_modal)

app.callback(Output('modal-app-characteristics', 'is_open'),
              [Input('open-app-characheristics', 'n_clicks'), Input('close-app-characteristics', 'n_clicks')],
              [State('modal-app-characteristics', 'is_open')])(toggle_modal)


#### JSME ####
@app.callback(Output('submit-output', 'children'),
              [Input('smiles-input', 'value'), Input('mol-input', 'value')])
def display_output(smiles, mol):
    if smiles:
        #chem = smiles
        return 'SMILES inserted.'
    elif mol:
        #chem = mol
        return 'Chemical structure drawn.'
    else:
        return 'Paste SMILES or Draw a Chemical.'

#### Submit Button ####
@app.callback([Output('loading-output', 'children'), Output('similarity_df', 'data')],
              [Input('submit-button', 'n_clicks')], state=[State('smiles-input', 'value'), State('mol-input', 'value')])
def load_output(click, smiles, mol):
    if click:
        if smiles:
            smol = get_mol(smiles)
            df = get_similarity_df(smol)
            return layout_results.layout, df.to_json(date_format='iso', orient='split')
        elif mol:
            smol = get_mol(mol)
            df = get_similarity_df(smol)
            return layout_results.layout, df.to_json(date_format='iso', orient='split')
    #return print('No structures submitted yet.'), None
    else:
        return []

if __name__ == '__main__':
    app.run_server(debug=True, dev_tools_hot_reload=False)
