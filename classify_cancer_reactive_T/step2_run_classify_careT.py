import papermill as pm

# Define the notebook path
notebook_path = 'step2_classify_careT.ipynb'

studies = ["Chen_2024", "Chow_2023", "Liu_2022", "Liu_2025", "Zheng_2021"]

# Define different sets of parameters
parameter_sets = [
    {'cell_type': "CD4", 'learn_rate': 0.001, 'batch_size': 32},
    {'cell_type': "CD8", 'learn_rate': 0.001, 'batch_size': 32}
]

for study in studies:
    for i, params in enumerate(parameter_sets):
        if study == "Chow_2023" and params['cell_type'] == "CD4":
            continue

        params['test_study'] = study
        output_path = f'{params["cell_type"]}_{study}_lr_{str(params["learn_rate"])}_bs_{str(params["batch_size"])}.ipynb'
        pm.execute_notebook(
            notebook_path,
            output_path,
            parameters=params
        )
        print(f'Notebook executed with parameters {params} and saved to {output_path}')
