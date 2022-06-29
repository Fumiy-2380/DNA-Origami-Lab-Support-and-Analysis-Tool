'''
DNA Origami Lab Support and Analysis Tool 
@Author Fumi inaba
'''
import json; import copy; import colorama; import math
import tkinter as tk; from tkinter import filedialog as fd
import numpy as np; import pandas as pd; from dataclasses import dataclass

def progress_bar(progress:int, total:int, string:str, color = colorama.Fore.YELLOW):
    percent = 100 * (progress / total)
    bar = '#' * int(percent) + '_' * (100 - int(percent))
    print(color + f'\r|{bar}| {percent:.2f}%', end = '\r')
    if progress == total:
        print(colorama.Fore.GREEN + f'\r|{bar}| {percent:.2f}%')



class OligoOperator:
    '''
    A class describing an object that will do a variety of operations on an oligonucleotide sequence such as getting the 
    complementary sequence or getting the melting temperature of an oligonucleotide sequence using REST API.
    '''
    def __init__(self, sequence: str):
        self.sequence = sequence
    
    @staticmethod
    def get_complement(oligo_seq:str) -> str:
        # Clean input sequence a little bit
        oligo_seq = oligo_seq.strip().upper()
        # The resulting complementary strand sequence
        complement = ""
        # Iterate through the 'nucleotides' in the string sequence 
        for nt in oligo_seq:
            if nt == "A":
                complement += "T"
            elif nt == "T":
                complement += "A"
            elif nt == "C":
                complement += "G"
            elif nt == "G":
                complement += "C"
            elif nt == "?":
                complement += "?"
        # Reverse so that the sequence is in the 5' -> 3' directionality
        complement = complement[::-1]
        return complement
    
    @staticmethod
    def get_melting_temperature(oligo_seq:str, idtAPI_obj, api_token:str, condition:tuple):
        Na, Mg, oligo_conc, dNTP_conc, nt_type = condition
        # Check the json cache file to see if the sequence has been queried to the online IDT API under conditions
        with open('DOLSAModule/melting_temp_cache.json', 'r') as melt_cache:
            melt_cache_contents = melt_cache.read()
        # If the melt cache is empty create an empty dictionary that will be written into the cache file
        if melt_cache_contents == "":
            melt_cache_dict = {}
        # If there are cache contents, read it in as a Python dictionary
        else:
            melt_cache_dict = json.loads(melt_cache_contents)
        
        # Check if the current request exists in the json cache file
        already_checked = f'{oligo_seq}-{Na}{Mg}{oligo_conc}{dNTP_conc}' in melt_cache_dict.keys()
        # Return the melting temperature if it's already in the cache
        if already_checked:
            melting_temperature = melt_cache_dict[f"{oligo_seq}-{Na}{Mg}{oligo_conc}{dNTP_conc}"]['Melting_temp']
            return melting_temperature
        # Otherwise use the IDT REST api object to get the melting temperature
        else:
            # Make request using the IDT API object method
            output = idtAPI_obj.get_melting_temp(oligo_seq, *condition, api_token)
            # Make the dictionary key of the request
            melt_cache_dict[f'{oligo_seq}-{Na}{Mg}{oligo_conc}{dNTP_conc}'] = {
                "Sequence": output['Sequence'],
                "Melting_temp": round(output['MeltTemp'], 1), 'Na_conc': output['NaConc'],
                'Mg_conc': output['MgConc'], 'Oligo_conc': output['OligoConc']
            }
            # Write results into the json cache file
            with open('DOLSAModule/melting_temp_cache.json', 'w') as cache_write:
                cache_write.write(json.dumps(melt_cache_dict, indent = 4))
            # Return just the melting temperature
            return round(output['MeltTemp'], 1)

@dataclass
class StapleStrand:
    oligo_sequence:str
    coordinates:list
    melting_temperature:float
    length: int

@dataclass
class Scaffold:
    dataframe: pd.DataFrame

@dataclass
class DNAnostructure:
    scaffold: Scaffold
    stapleList: list

    def __repr__(self):
        return "This is a DNA nanostructure"

class DataCleaningModule:
    def __init__(self):
        self.vstrands, self.oligo_seq_list = self.select_files()

    def select_files(self) -> tuple:
        # Create a root tkinter widget to select the json and csv files exported from caDNAno using tkinter GUI
        root = tk.Tk()
        root.withdraw() # Hide the root widget
        # First select the json file
        json_path = fd.askopenfilename(title = "Select JSON file exported from caDNAno: ")
        # Then select the csv file
        csv_path = fd.askopenfilename(title = 'Select CSV file exported from caDNAno: ')
        # Open the files and read in their contents
        with open(json_path, 'r') as jsonf:
            jsonf_contents = json.loads(jsonf.read())
        with open(csv_path, 'r') as csvf:
            csvf_contents = [line.split(',')[2] for line in csvf.readlines()[1:]]
        # Return the file contents
        return jsonf_contents['vstrands'], csvf_contents

    def clean_and_fragment(self, strand_type:str = 'scaf') -> list:
        # Container for nucleotide (nt) coordinates contained in the json file - removes empty coordinates ([-1, -1, -1, -1])
        cleaned_helix_pool = []
        # Make a copy of self.vstrands so there are no side effects of list operations on self.vstrands
        virtual_helices = copy.deepcopy(self.vstrands)
        # Iterate through each virtual helix 
        for vHelix in virtual_helices:
            # Change orientation to 5' -> 3' directionality to make subsequent cleaning steps easier
            virtual_helix = vHelix[strand_type]
            vHelix_num = vHelix['num']
            even = True if vHelix_num % 2 == 0 else False
            if strand_type == 'scaf':
                if not even:
                    virtual_helix.reverse()
            elif strand_type == 'stap':
                if even:
                    virtual_helix.reverse()
            # Container for 'fragments' which are regions of DNA between crossovers or between the end of a strand
            strand_fragment = []
            for nt in virtual_helix:
                # If the nt coordinate stored in the list is not empty (empty coordinates stored as [-1, -1, -1, -1])
                if nt != [-1, -1, -1, -1]:
                    # Coordinates of 'nt' are in the form of [a, b, c, d] where a, b denote the base on the 5' side of the nt
                    # and c, d denote the base on the 3' side of the nt
                    helix_5_side, len_5_side, helix_3_side, len_3_side = nt
                    # A parameter to account for the offset of some base coordinates when a crossover happens
                    offset = 1 if strand_type == 'scaf' else -1 # direction of offset differs between scaffold (scaf) and staple (stap) strand_type
                    # When helix_5_side != helix_3_side, the bases on the left and right of the nt are on different virtual helices
                    # If there is a cross over:
                    if not even: 
                        offset *= -1
                    if helix_5_side != helix_3_side:
                        if helix_5_side != -1 and helix_3_side != -1: 
                            if helix_5_side == vHelix_num:
                                sorting_tag = (vHelix_num, len_5_side + offset, (helix_3_side, (len_5_side + offset)))
                                strand_fragment.append((nt, sorting_tag))
                                cleaned_helix_pool.append(strand_fragment); strand_fragment = []
                            elif helix_3_side == vHelix_num:
                                sorting_tag = (vHelix_num, len_3_side - offset, (helix_5_side, (len_3_side - offset)))
                                strand_fragment.append((nt, sorting_tag))
                        # If the nt base coordinate is a 5' end of a strand
                        elif helix_5_side == -1:
                            sorting_tag = (vHelix_num, len_3_side - offset, (-1, (len_3_side - offset)))
                            strand_fragment.append((nt, sorting_tag))
                        # If the nt base coordinate is a 3' end of a strand
                        elif helix_3_side == -1:
                            sorting_tag = (vHelix_num, len_5_side + offset, (-1, (len_5_side + offset)))
                            strand_fragment.append((nt, sorting_tag))
                            cleaned_helix_pool.append(strand_fragment); strand_fragment = []
                    # If the base is a part of a continuous portion of the strand
                    else:
                        sorting_tag = (vHelix_num, len_5_side + offset, (None, None))
                        strand_fragment.append((nt, sorting_tag))
        
        return cleaned_helix_pool

    def search_continous_nt(self, nt_fragments:list, nt_to_find:list) -> int:
        # unpack the coordinates fo nt_to_find
        helix_5, len_5, helix_3, len_3 = nt_to_find
        # Linear search through the list of fragments (nt_fragments) to find the next fragment to connect to
        for fragment in nt_fragments:
            head = fragment[0][0]
            #print(f"The head is {head}")
            #print(f"The nt to find is {nt_to_find}")
            # If the inputted base coordinate matches the head of a fragment, the fragment is the next fragment to connect to
            if nt_to_find == head:
                index = nt_fragments.index(fragment)
                return index
            # Search for a certain pattern of base coordinate that indicates it is the next fragment to connect to
            elif head == [helix_5, len_5 + 1, helix_3, len_3 - 1] or head == [helix_5, len_5 - 1, helix_3, len_3 + 1]:
                index = nt_fragments.index(fragment)
                return index

        return None

    def connect_fragments(self, fragment_pool:list, linearize:bool = False) -> list:
        # Initialize a list to store the reordered fragments from 5' -> 3'
        reordered_frags = []
        # Lienar search to get the first starting 5' fragment 
        #print("starting loop")
        for fragment in fragment_pool: # fragment is a list of tuples
            # Condition for the fragment to be the starting 5' end
            #print(fragment)
            if fragment[0][0][0] == -1:
                #print("condition found")
                # Rename as the first 5' fragment
                initial_5_fragment = fragment
                reordered_frags.append(initial_5_fragment)
                #print(initial_5_fragment)
                # Select the last base coordinate and store in tail_nt
                tail_nt = initial_5_fragment[-1][0]
                # Delete the selected fragment from the list to make subsequent linear searches faster
                fragment_pool.pop(fragment_pool.index(initial_5_fragment))
                # Check if the fragment is in itself a complete strand - if so return the strand
                if tail_nt[-1] == -1:
                    return initial_5_fragment
                # Break the for loop and move on to next step
                break
        
        # Continue while the tail_nt is not -1, or not the 3' end 
        while tail_nt[-1] != -1:
            # Call search_continous_nt() method to get the index of the next fragment to append to reordered_frags
            next_index = self.search_continous_nt(fragment_pool, tail_nt)
            # Update the tail_nt
            tail_nt = fragment_pool[next_index][-1][0]
            # Append the resulting next fragment to reordered_frags list
            reordered_frags.append(fragment_pool[next_index])
            # Remove selected elements to make subsequent linear searches faster
            fragment_pool.pop(next_index)
            # If linearize = False, conserve the crossover shapes and return a 3D list
            if tail_nt[-1] == -1 and linearize == False:
                return reordered_frags
            # If linearize is true, return a 2D list, a list of all the base coordinates which are lists themselves
            elif tail_nt[-1] == -1 and linearize == True:
                linearized = []
                for frag in reordered_frags:
                    for nt in frag:
                        linearized.append(nt)
                return linearized
    
    #Scaffold processing methods
    def longest_scaffold(self, unordered_fragments:list) -> list:
        '''
        Sometimes there are multiple conntinous scaffold strands in a design JSON file, this selects the longest of them.
        '''
        # list to store linearized scaffolds
        scaffold_coordinates = []
        # Continue while the length of the input list is not empty
        while len(unordered_fragments) > 0:
            # Get coordinates for a single scaffold strand
            scaffold = self.connect_fragments(unordered_fragments, linearize = True)
            scaffold_coordinates.append(scaffold)
        # Select the longest scaffold    
        longest = max(scaffold_coordinates, key = len)
        return longest
    
    def assign_scaffold_sequence(self, longest_scaffold:list, scaf_id:str) -> list:
        # Open the file containing scaffold sequences
        with open("DOLSAModule/scaffolds.json", 'r') as scaf:
            scaf_seq = json.loads(scaf.read())[scaf_id] # Select one sequence of choice
        # Make the complementary sequence using the static method from OligoOperator class
        c_scaf_seq = OligoOperator.get_complement(scaf_seq)[::-1]
        # The zipped version with the coordinate, nt base and the complementary base for each base coordinate
        scaf_complete = list(zip(longest_scaffold, scaf_seq, c_scaf_seq))
        # Sort each base by the sorting_tag made earlier steps
        sorted_scaf_complete = sorted(scaf_complete, key = lambda nt: nt[0][1])
        return sorted_scaf_complete
    
    def create_scaffold_dataframe(self, scaf_id:str) -> pd.DataFrame:
        # The list of dictionaries that will be the input for the Dataframe
        df_input = []
        # Call self methods to prepare the data 
        # Make the pool of fragments
        fragment_pool = self.clean_and_fragment()
        # Select the longest scaffold
        longest_scaffold = self.longest_scaffold(fragment_pool)
        # Assign a nucleotide and complementary nucleotide base 
        scaf_complete = self.assign_scaffold_sequence(longest_scaffold, scaf_id)
        # For each nculeotide in the scaffold
        for nt in scaf_complete:
            # Dictionary to store data about the nucleotide base 
            tempDict = {}
            tempDict['Nucleotide'] = nt[1] # The nucleotide base sequence
            tempDict['Complement'] = nt[2] # The complementary base sequence
            tempDict['Sorting Tag - Helix'] = nt[0][1][0] # The vHelix number of the sorting tag
            tempDict['Sorting Tag - Length'] = nt[0][1][1] # The legnth position of the sorting tabg
            tempDict['Sorting Tag - Xover'] = nt[0][1][2] # The crossover tag
            tempDict['Sorting Tag'] = (nt[0][1][0], nt[0][1][1]) # The sorting tag as a tuple
            tempDict['Coordinate'] = nt[0][0] # The coordinate of the nucleotide base
            df_input.append(tempDict) # Append the dictionary to the df_input list
        # Create dataframe
        df = pd.DataFrame(df_input)
        return df            

    # Staple analyzing methods
    def get_all_staple_coordinates(self, unordered_staple_fragments:list) -> list: 
        # To accumulate all the staple coordinates that have been reordered
        ordered_staple_coords = []
        # While there are still fragments left in the list passed in as a parameter
        while len(unordered_staple_fragments) > 0:
            # Reorganize the staple coordinates to make them a continuous strand
            staple = self.connect_fragments(unordered_staple_fragments, linearize = True)
            ordered_staple_coords.append(staple)
        return ordered_staple_coords
    
    def search_complement_sequence(self, staple_coords:list, df:pd.DataFrame) -> str:
        '''
        Binary search algorithm for complementary base
        '''
        # The nucleotide sequence of the staple coordinates
        stap_strand_seq = ""

        for base in staple_coords:
            # Search the dataframe with scaffold information using the sorting tag
            search_tag = (base[1][0], base[1][1])
            try:
                query_index = df.index[df['Sorting Tag'] == search_tag]
                stap_strand_seq = stap_strand_seq + df.at[query_index[0], 'Complement']
            # If the staple strand is single stranded, and not paired to a scaffold strand
            except IndexError:
                stap_strand_seq += "?"

        return stap_strand_seq
    
    def accumulate_staples(self, ordered_staple_fragments:list, scaffold_df:pd.DataFrame, idtapiobj, access_key:str, conditions:tuple = None):
        # A list that will contain all of the staple strand dataclass objects
        staples = []
        # Initialize progress bar
        progress_bar(0, len(ordered_staple_fragments),'')
        # For each staple strand coordinate in the passed in list
        for staple_coord in ordered_staple_fragments:
            # Applya sequence to it
            sequence = self.search_complement_sequence(staple_coord, scaffold_df)
            # Conditions is for the melting temperature. 
            if conditions != None:
                # Request for the melting temperature using REST API
                melt_temp = OligoOperator.get_melting_temperature(sequence, idtapiobj, access_key, conditions)
            # If conditions are not defined, use the defaults outlined in idtAPI_obj.get_melting_temperature() method
            else:
                melt_temp = OligoOperator.get_melting_temperature(sequence, idtapiobj, access_key)
            # Initialize staple strand object
            staple_obj = StapleStrand(sequence, staple_coord, melt_temp, len(sequence))
            staples.append(staple_obj)
            # To keep track of progress, print progress to the terminal
            progress_bar(len(staples), len(ordered_staple_fragments), staple_obj.oligo_sequence)
        print(colorama.Fore.RESET)
        return staples


class AnalysisModule:
    @staticmethod
    def strand_similarity(staple_object1:StapleStrand, staple_object2:StapleStrand, check_frame:int = 2, complementary:bool = False) -> int:
        # Define the shorter and longer strands
        shorter = staple_object1 if staple_object1.length <= staple_object2.length else staple_object2
        longer = staple_object2 if staple_object2.length > staple_object1.length else staple_object1

        # If complementary == True, the compare the complementary strand of one staple object to the sequence of the other
        if complementary == True:
            # Get the complementary sequence
            compl_seq = OligoOperator.get_complement(shorter.oligo_sequence)
            # Iterate through each position in the shorter string and take a slice of length check_frame
            for i in range(len(compl_seq)-check_frame):
                # If that slice is found in the longer string, recursively check a longer check_frame
                if compl_seq[i:i + check_frame] in longer.oligo_sequence:
                    return AnalysisModule.strand_similarity(staple_object1, staple_object2, check_frame + 1, complementary)
                # If not, return the largest value of check_frame when a match was found
                else:
                    return check_frame - 1
            return 0
        # Same as above, but compare the similarity between the raw sequences of the two staple strand objects
        else:
            for i in range(shorter.length - check_frame):
                if shorter.oligo_sequence[i:i + check_frame] in longer.oligo_sequence:
                    return AnalysisModule.strand_similarity(staple_object1, staple_object2, check_frame + 1, complementary)
                else:
                    return check_frame - 1
            return 0
    
    @staticmethod
    def detect_sandwich_strands(staple_obj:StapleStrand):
        fragments_len = []
        for stap_frag in staple_obj.coordinates:
            fragments_len.append(len(stap_frag))
        forward_iter = 0
        backward_iter = -1 
        # Sandwich strands occur when a short fragment of a staple is flanked by two longer fragments
        while fragments_len[forward_iter] != fragments_len[backward_iter]:
            # Check length of fragments from the front
            if fragments_len[forward_iter] > fragments_len[forward_iter + 1]:
                return True
            # Check length of fragments from the back
            if fragments_len[backward_iter] > fragments_len[backward_iter - 1]:
                return True
            # Update iteration trackers
            forward_iter += 1
            backward_iter -= 1
            # If a sandwich strand is not detected and forward and backward iteration reaches the same fragment, it is not a sandwich strand
            if fragments_len[forward_iter] is fragments_len[backward_iter]:
                return False
    
    @staticmethod
    def find_melting_temp_binds(melt_temps:list):
        # turn input parameter list into a NumPy array
        arrayMT = np.asarray(melt_temps)
        # Get lowest and highest melting temperatures
        lowest_MT = math.trunc(arrayMT.min())
        highest_MT = math.ceil(arrayMT.max())
        # Define the bins of the melting temperatures with increments of 1 deg Celsius
        bins = range(lowest_MT, highest_MT + 1)
        # Make the histogram
        hist = np.histogram(arrayMT, bins)

        return hist
    
    @staticmethod
    def make_thermal_ramp(melt_temps:list):
        anneal_time_per_strand = 1 # Minutes per strand
        # Make a histogram with bins of strand melting temperatures
        hist = AnalysisModule.find_melting_temp_binds(melt_temps)

        # Determine the range of the x axis. The x axis is in minutes. 
        x_axis_increments = np.sum(hist[0]*anneal_time_per_strand)
        # Multiplying all counts of strand melting temps and multiplying by anneal_time_per_strand gives the total time it will take
        x_axis = np.linspace(0, x_axis_increments, x_axis_increments)

        y_axis = np.array([])
        # For length of bins in the histogram, make y points in between the bin markers
        for i in range(len(hist[1]) - 1):
            # Make the y points in between the x_axis bins (annealing time per strand * the frequency of strands in each bin)
            y_point = np.linspace(hist[1][i], hist[1][i + 1], hist[0][i] * anneal_time_per_strand, endpoint = True)
            # This makes it so that the more strands there are with melting temperatures between bins, the longer it will take to decrease system temp by 1C. 
            y_axis = np.append(y_axis, y_point)
        y_axis = np.flip(y_axis)
        
        return x_axis.tolist(), y_axis.tolist()

    @staticmethod
    def analysis_matrix(staple_objs:list, staple_objs2:list = None) -> dict:
        if staple_objs2 == None:
            staple_objs2 = staple_objs
        strand_id = {}
        df_input = []
        for staple_obj in staple_objs:
            indexA = staple_objs.index(staple_obj)
            if strand_id.get(indexA, "") == "":
                strand_id[indexA] = staple_obj.oligo_sequence
            for staple_obj2 in staple_objs2:
                indexB = staple_objs2.index(staple_obj2)
                # Dictionary for sequence ID mapping to sequence
                if strand_id.get(indexB, "") == "":
                    strand_id[indexB] = staple_obj2.oligo_sequence
                # Dictionary for dataframe input
                dict_comparison = {}
                num_similar = AnalysisModule.strand_similarity(staple_obj, staple_obj2)
                num_complement = AnalysisModule.strand_similarity(staple_obj, staple_obj2, complementary = True)
                dict_comparison['SequenceA'] = indexA
                dict_comparison['SequenceB'] = indexB
                dict_comparison['Similarity'] = num_similar
                dict_comparison['Complementarity'] = num_complement
                df_input.append(dict_comparison)

        returnMatrix = pd.DataFrame(df_input)
        returnMatrix.to_csv("Outputs/ComparisonMatrix.csv")

        return strand_id