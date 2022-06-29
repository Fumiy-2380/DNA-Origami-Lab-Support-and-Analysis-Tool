import json; import test; from DOLSAModule import idtapi
import numpy as np; import pandas as pd
from .dolsat import DataCleaningModule
from .dolsat import AnalysisModule; from DOLSAModule import dolsat


class Main:
    @staticmethod
    def mainscript():
        Main.output_thermal_ramp_csv()
        strand_mapping = AnalysisModule.analysis_matrix(Main.staple_strands)
        try:
            with open('Outputs/Strand_mapping.txt', 'w') as strand_map:
                for key, value in strand_mapping.items():
                    strand_map.write(f"{key}: {value} \n")
        except FileNotFoundError:
            with open('Outputs/Strand_mapping.txt', 'x') as strand_map:
                for key, value in strand_mapping.items():
                    strand_map.write(f"{key}: {value} \n")                

    @staticmethod
    def output_thermal_ramp_analysis(datamodule: DataCleaningModule, DNA_content:tuple, idtapi_contents:tuple, conditions:tuple) -> list:
        # Unpack input parameters
        ordered_stap_coords, scaffold_dataframe = DNA_content
        idtapi_obj, access_key = idtapi_contents
        Na_conc, Mg_conc, oligo_conc, dNTP, strand_type = conditions
        # Define tuple of condition to calculate melting temperature 
        strand_condition = conditions
        # Get staple strand instances
        staple_strands = datamodule.accumulate_staples(ordered_stap_coords, scaffold_dataframe, idtapi_obj, access_key, conditions = strand_condition)
        Main.staple_strands = staple_strands
        # list up melting temperatures of the staple strands
        list_of_melt_temps = [staple_obj.melting_temperature for staple_obj in staple_strands]
        # Create the thermal ramp
        thermal_ramp = AnalysisModule.make_thermal_ramp(list_of_melt_temps)
        thermal_ramp_axes = list(zip(thermal_ramp[0], thermal_ramp[1]))
        # Make the temporary dictionary to be passed into the dataframe input
        df_input = []
        for time, temp in thermal_ramp_axes:
            melt_temp_dict = {}
            melt_temp_dict["Na"] = Na_conc
            melt_temp_dict['Mg'] = Mg_conc
            melt_temp_dict['Strand'] = "DNA"
            melt_temp_dict['Timestamp'] = time
            melt_temp_dict['Temperature'] = temp
            df_input.append(melt_temp_dict)
        
        return df_input

    @staticmethod
    def output_thermal_ramp_csv() -> None:
        with open("DOLSAModule/idt_api_credentials.json", "r") as credentials_file:
            cfd = json.loads(credentials_file.read())
        
        credentials_tup = (cfd["client_id"], cfd['client_secret'], cfd['idt_username'], cfd['idt_password'])
        idtapi_obj = idtapi.IDTApiDataExtractor(baseurl = "https://www.idtdna.com/Restapi/v1/OligoAnalyzer/Analyze/", *credentials_tup)
        access_key = idtapi_obj.find_access_token_cache()
        # Pack into tuple to make data input easier
        idt_contents_tup = (idtapi_obj, access_key)

        # Initialize data cleaning module
        datamodule = dolsat.DataCleaningModule()
        # Clean scaffold data
        scaffold_dataframe = datamodule.create_scaffold_dataframe('p8064')
        # Clean staple data
        staple_fragment_pool = datamodule.clean_and_fragment('stap')
        ordered_staple_coords = datamodule.get_all_staple_coordinates(staple_fragment_pool)
        # Pack into tuple to make data input easier
        dna_content_tup = (ordered_staple_coords, scaffold_dataframe)

        df_inputs = []
        for Na_conc in range(30, 90, 20):
            for Mg_conc in range(8, 28, 4):
                condition_for_melt_temp = (Na_conc, Mg_conc, 0.45, 0, "DNA")
                ramps_for_cond = Main.output_thermal_ramp_analysis(datamodule, dna_content_tup, idt_contents_tup, condition_for_melt_temp) # list of dictionary
                df_inputs.append(ramps_for_cond)
        
        df_inputs = sum(df_inputs, [])
        thermal_ramp_df = pd.DataFrame(df_inputs)
        thermal_ramp_df.to_csv('Outputs/ThermalRampOutput.csv')

        return None
    
if __name__ == "__main__":
    main = Main()
    main.mainscript()