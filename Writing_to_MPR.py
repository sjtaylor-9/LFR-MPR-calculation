"""
Writing_to_MPR.py

This code is responsible for outputting the nuclide information and reaction cross-sections to a .txt file in the specific format required for the MPR file in ORION.
This is done by loading the ZAIDs of the parent/daughter nuclides from a .csv, the cross-sections of each reaction from a .csv, and the ORION IDs of the all the nuclides from a .xlsx. These are all stored as numpy arrays or pandas dataframes.
The code first sorts the pandas dataframe containing the nuclide information into ascending ORION IDs and then iterates through each row of the dataframe and writes all the information into the .txt file.

A description of how to use the scripts is given in the README of this repo: https://github.com/sjtaylor-9/LFR-MPR-calculation.git
Author: Sam Taylor (sam.taylor@newcleo.com)
Last Edited: 13/01/2025
"""
# ---------  Import Libraries  --------- #
import numpy as np
import pandas as pd
import argparse
import os
# -------------  Functions  ------------ #
def parse_arguments():
    """
    Parses the arguments needed along the code. Arguments:

    --reactor               Used to specify which reactor out of the LFR200 and LFR30 is being examined.
                            The argument must be one of: [LFR30, LFR200].
    --enrichment_zone       If the LFR200 is being examined, then the enrichment zone must be specified. The data for the LFR30 neutron spectrum only describes one enrichment zone, the fuel assembly.
                            The argument must be one of: [inner, middle, outer, FA].
    --path                  Used to specify the directory in which the one-group cross-section data should be written. It is not required,
                            in the case it is not specified, the default path is the current working directory.
    --newcleo_input         Used to specify the directory in which internal newcleo inpput data should be found. This includes data such as the neutron spectra.
    --cross_section_data    Used to specify the directory in which the public data regarding the ENDF and JEFF reaction cross-sections can be found.
    --reaction_data
    --burnup            
    
    Returns the parsed arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--reactor",
        type = str,
        choices = ['LFR30', 'LFR200'],
        required = True,
        help = 'flag to set whether the LFR30 or LFR200 is to be examined'
    )
    parser.add_argument(
        "--enrichment_zone",
        type = str,
        choices = ['inner', 'middle', 'outer', 'FA'],
        required = True,
        help = 'flag to set which enrichment zone in the LFR200 is to be examined'
    )
    parser.add_argument(
        "--path",
        type=dir_path,
        required=False,
        default=os.getcwd(),
        help="flag to set the path where the output files should be written to"
    )
    parser.add_argument(
        "--cross_section_data",
        type=dir_path,
        required=True,
        help="flag to set the path where the input cross-section data should be found"
    )
    parser.add_argument(
        "--reaction_data",
        type=dir_path,
        required=True,
        help="flag to set the path where the input reaction data should be found"
    )
    parser.add_argument(
        "--input",
        type=dir_path,
        required=True,
        help="flag to set the path where the input ZAID data should be found"
    )
    return parser.parse_args()

def dir_path(string):
    '''
    Checks if a given string is the path to a directory.
    If affirmative, returns the string. If negative, gives an error.
    '''
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)

def calculate_parent_orion(parent_zaid, orion_data):
    """
    Determines the ORION ID of the parent nuclide.
    This is done by reading and comparing the parent nuclide name to the list of nuclide names in the ORION database, and finding the ORION ID as the index of the row in the ORION array.
    The Excel spreadsheet containing the nuclides in the ORION database is loaded into a pandas dataframe and the ORION ID is determined by the cell row number.

    Args:
        parent_zaid (int): The ZAID of the parent nuclide.
        orion_data (pandas dataframe): The pandas dataframe that contains all of the ORION nuclides and the ORION IDs are the respective indexes.

    Returns:
        parent_orion (int): The ORION ID of the parent nuclide.
    """
    parent_name = build_parent_daughter_name(parent_zaid)
    # In the ORION database the nuclide names are in the form 12-MG33. 
    # Removes the N- from each cell so that the nuclide names are in the same format as in ZAID_results.csv file 
    df_nuclide_without_prefix = [element.split('-', 1)[1] for element in orion_data]

    # Iterates through each row in the array containing the ORION IDs
    for row in df_nuclide_without_prefix:
        # Compares the current row of the ORION ID array to the name of the parent, if these match then the for loop is exited
        if row == parent_name[0]:
            parent_orion = df_nuclide_without_prefix.index(row) # The ORION ID of the parent is found as the index of the current row in the array
            break
    return parent_orion

def calculate_daughter_orion(daughter_zaid, orion_data):
    """
    Determines the ORION ID of the daughter nuclide.
    This is done by reading and comparing the daughter nuclide name to the list of nuclide names in the ORION database, and finding the ORION ID as the index of the row in the ORION array.
    The Excel spreadsheet containing the nuclides in the ORION database is loaded into a pandas dataframe and the ORION ID is determined by the cell row number.


    Args:
        daughter_zaid (int): The ZAID of the daughter nuclide.
        orion_data (pandas dataframe): The pandas dataframe that contains all of the ORION nuclides and the ORION IDs are the respective indexes.

    Returns:
        daughter_orion: int
    """
    daughter_name = build_parent_daughter_name(daughter_zaid)
    df_nuclide_without_prefix = [element.split('-', 1)[1] for element in orion_data]
    
    for row in df_nuclide_without_prefix:
        if row == daughter_name[0]:
            daughter_orion = df_nuclide_without_prefix.index(row)
            break

    return daughter_orion

def build_parent_daughter_name(id_number):
    """
    Constructs the name of the parent/daughter nuclide in the form PU243 from it's ZAID.
    The ZAID is in the form ZZAAAI, so the function first removes I(=0), and then sets Z to the first 2 digits and A as the rest, whilst removing all leading/trailing zeroes.

    Args:
        id_number (float): The ZAID number of the parent/daughter nuclide, in the form 10010.0.

    Returns:
        nuclear_notation (str): The name of the parent/daughter nuclide in the form PU243.
        atomic_number (int): The atomic number of the parent/daughter nuclide.
        mass_number (int): The mass number of the parent/daughter nuclide.
    """
    # Format the ZAIDs correctly, in both str and int types, such that .0 is removed and the trailing (I=0) is removed as well
    id_number = int(id_number)
    unconverted_id_number = id_number
    id_number = id_number // 10
    id_string = str(id_number)
    
    # Set Z to the first digit if the second digit is '0' and the nuclide ID number is less than 9999
    if id_number < 10000:
        atomic_number = int(id_string[0])
    else:
        # Otherwise, set Z to the first two digits
        atomic_number = int(id_string[:2])
   
    # Mass number is set to figure left after removing first 2 digits
    mass_number = id_string[2:]

    # Remove leading zeros from A; if A is empty after stripping, set it to 0
    mass_number = int(float(mass_number.lstrip('0'))) if mass_number.lstrip('0') else 0    
    
    # Nuclide name in form PU243
    nuclear_notation = element_symbols[atomic_number-1].upper() + str(mass_number)

    # If the parent or daughter is metastable, include the excited state in the nuclide name
    if (unconverted_id_number % 10 == 1):
        nuclear_notation = nuclear_notation + "(1)"
    if (unconverted_id_number % 10 == 2):
        nuclear_notation = nuclear_notation + "(2)"

    return nuclear_notation, atomic_number, mass_number

def get_cross_section(zaid_in_MPR, reaction, burnup):
    """
    Finds the one-group cross-section of the reaction from the list of cross-sections in the .csv file.

    Args:
        zaid_in_MPR (int): The ZAID of the nuclide.
        reaction (int): The MT value of the reaction (16, 17, 18, 102, 103).

    Returns:
        xs (float): The one-group cross-section of the specific reaction.
    """
    # Opens the csv file containing the cross-section data for each nuclide and stores it in an array
    xs_file = f'{args.cross_section_data}/{burnup}/{args.reactor}_{args.enrichment_zone}_{burnup}_all_reactiondata.csv'
    xs_data = np.genfromtxt(xs_file, comments = '%', delimiter = ',')
    # Iterates through the csv file until the required reaction is found and the cross-section is outputted
    for row in xs_data:
        if row[0] == zaid_in_MPR and row[1] == reaction:
            xs = row[3]
            break
        
    return xs

def save_to_MPR(ORION, 
                nuclide_ID, 
                name, 
                NumberReactions, 
                daughter_16,
                daughter_17,
                daughter_18,
                daughter_102,
                daughter_103,
                NumberParents,
                parent_16,
                parent_17,
                parent_102,
                parent_103,
                parent_16_1m,
                parent_17_1m,
                parent_102_1m,
                parent_103_1m,
                parent_16_2m,
                parent_17_2m,
                parent_102_2m,
                parent_103_2m,
                orion_dataframe):
    """
    Writes all of the nuclide information and reaction data to the .txt file in the format of the MPR file.

    Args:
        ORION (int): The ORION ID of the nuclide.
        nuclide_ID (int): The ZAID of the nuclide.
        name (string): The name of the nuclide in the form PU243.
        NumberReactions (int): The total number of reactions in the database that the nuclide can undergo.
        daughter_16 (float): The ZAID of the daughter nuclide of reaction 16.
        daughter_17 (float): The ZAID of the daughter nuclide of reaction 17.
        daughter_18 (float): The ZAID of the daughter nuclide of reaction 18.
        daughter_102 (float): The ZAID of the daughter nuclide of reaction 102.
        daughter_103 (float): The ZAID of the daughter nuclide of reaction 103.
        NumberParents (int): The total number of parents that the nuclide has in the database.
        parent_16 (float): The ZAID of the parent nuclide of reaction 16.
        parent_17 (float): The ZAID of the parent nuclide of reaction 17.
        parent_102 (float): The ZAID of the parent nuclide of reaction 102.
        parent_103 (float): The ZAID of the parent nuclide of reaction 103.
    """
    # The ORION and ZAIDs are right aligned with total widths of 5 and 8, respectively, and the nuclide name is left aligned with a total width of 11
    width_orion = 5
    aligned_orion = "{:>{width}}".format(ORION, width = width_orion)
    width_nuclide = 8
    aligned_nuclide = "{:>{width}}".format(nuclide_ID, width = width_nuclide)
    width_name = 11
    aligned_name = "{:<{width}}".format(name, width = width_name)
    
    # NumberReactions and NumberParents should always be < 10 but have right aligned as validation
    width_num_reactions = 5
    aligned_reactions = "{:>{width}}".format(NumberReactions, width = width_num_reactions)
    width_num_parents = 4
    aligned_parents = "{:>{width}}".format(NumberParents, width = width_num_parents)
    file.write(f"{aligned_orion}{aligned_nuclide}  {aligned_name}NumberReactions{aligned_reactions} NumberParents{aligned_parents}\n")
    
    # Adds the parent IDs in the required format
    # Order of reactions in increasing ORION IDs is 102, 16, 17, 103, no ORION ID for reaction 18
    # Text for both IDs must be right aligned both with total width of 12
    width_parent_orion = 12
    width_parent_zaid = 12
    # If the parent ZAID is not 'NA' then the parent information is written into the MPR file
    if not np.isnan(parent_102):
        parent_orion_102 = calculate_parent_orion(parent_102, orion_dataframe)
        # Right aligns the text with the stated given width
        aligned_orion = "{:>{width}}".format(parent_orion_102, width = width_parent_orion)
        aligned_zaid = "{:>{width}}".format(int(parent_102), width = width_parent_zaid) # The parent ZAID is changed to type int to remove the .0
        file.write(f"{aligned_orion}{aligned_zaid}\n")
    if not np.isnan(parent_102_1m):
        parent_orion_102_1m = calculate_parent_orion(parent_102_1m, orion_dataframe)
        aligned_orion = "{:>{width}}".format(parent_orion_102_1m, width = width_parent_orion)
        aligned_zaid = "{:>{width}}".format(int(parent_102_1m), width = width_parent_zaid) # The parent ZAID is changed to type int to remove the .0
        file.write(f"{aligned_orion}{aligned_zaid}\n")
    if not np.isnan(parent_102_2m):
        parent_orion_102_2m = calculate_parent_orion(parent_102_2m, orion_dataframe)
        aligned_orion = "{:>{width}}".format(parent_orion_102_2m, width = width_parent_orion)
        aligned_zaid = "{:>{width}}".format(int(parent_102_2m), width = width_parent_zaid) # The parent ZAID is changed to type int to remove the .0
        file.write(f"{aligned_orion}{aligned_zaid}\n")
    if not np.isnan(parent_16):
        parent_orion_16 = calculate_parent_orion(parent_16, orion_dataframe)
        aligned_orion = "{:>{width}}".format(parent_orion_16, width = width_parent_orion)
        aligned_zaid = "{:>{width}}".format(int(parent_16), width = width_parent_zaid)
        file.write(f"{aligned_orion}{aligned_zaid}\n")
    if not np.isnan(parent_16_1m):
        parent_orion_16_1m = calculate_parent_orion(parent_16_1m, orion_dataframe)
        aligned_orion = "{:>{width}}".format(parent_orion_16_1m, width = width_parent_orion)
        aligned_zaid = "{:>{width}}".format(int(parent_16_1m), width = width_parent_zaid)
        file.write(f"{aligned_orion}{aligned_zaid}\n")
    if not np.isnan(parent_16_2m):
        parent_orion_16_2m = calculate_parent_orion(parent_16_2m, orion_dataframe)
        aligned_orion = "{:>{width}}".format(parent_orion_16_2m, width = width_parent_orion)
        aligned_zaid = "{:>{width}}".format(int(parent_16_2m), width = width_parent_zaid)
        file.write(f"{aligned_orion}{aligned_zaid}\n")
    if not np.isnan(parent_17):
        parent_orion_17 = calculate_parent_orion(parent_17, orion_dataframe)
        aligned_orion = "{:>{width}}".format(parent_orion_17, width = width_parent_orion)
        aligned_zaid = "{:>{width}}".format(int(parent_17), width = width_parent_zaid)
        file.write(f"{aligned_orion}{aligned_zaid}\n")
    if not np.isnan(parent_17_1m):
        parent_orion_17_1m = calculate_parent_orion(parent_17_1m, orion_dataframe)
        aligned_orion = "{:>{width}}".format(parent_orion_17_1m, width = width_parent_orion)
        aligned_zaid = "{:>{width}}".format(int(parent_17_1m), width = width_parent_zaid)
        file.write(f"{aligned_orion}{aligned_zaid}\n")
    if not np.isnan(parent_17_2m):
        parent_orion_17_2m = calculate_parent_orion(parent_17_2m, orion_dataframe)
        aligned_orion = "{:>{width}}".format(parent_orion_17_2m, width = width_parent_orion)
        aligned_zaid = "{:>{width}}".format(int(parent_17_2m), width = width_parent_zaid)
        file.write(f"{aligned_orion}{aligned_zaid}\n")
    if not np.isnan(parent_103):
        parent_orion_103 = calculate_parent_orion(parent_103, orion_dataframe)
        aligned_orion = "{:>{width}}".format(parent_orion_103, width = width_parent_orion)
        aligned_zaid = "{:>{width}}".format(int(parent_103), width = width_parent_zaid)
        file.write(f"{aligned_orion}{aligned_zaid}\n")
    if not np.isnan(parent_103_1m):
        parent_orion_103_1m = calculate_parent_orion(parent_103_1m, orion_dataframe)
        aligned_orion = "{:>{width}}".format(parent_orion_103_1m, width = width_parent_orion)
        aligned_zaid = "{:>{width}}".format(int(parent_103_1m), width = width_parent_zaid)
        file.write(f"{aligned_orion}{aligned_zaid}\n")
    if not np.isnan(parent_103_2m):
        parent_orion_103_2m = calculate_parent_orion(parent_103_2m, orion_dataframe)
        aligned_orion = "{:>{width}}".format(parent_orion_103_2m, width = width_parent_orion)
        aligned_zaid = "{:>{width}}".format(int(parent_103_2m), width = width_parent_zaid)
        file.write(f"{aligned_orion}{aligned_zaid}\n")
    # Adds the relevant data for the daughter nuclides in the required format
    # The ORION ID of the daughter nuclide must be right aligned with a total width of 7 (6 preceeding spaces for 1 digit) and ZAID must be right aligned with a total width of 8
    width_orion = 7
    width_zaid = 8
    # If the daughter ZAID is not 'NA' then the daughter reaction information is written into the MPR file
    if not np.isnan(daughter_16):
        daughter_orion_16 = calculate_daughter_orion(daughter_16, orion_dataframe)
        aligned_orion = "{:>{width}}".format(daughter_orion_16, width = width_orion)
        aligned_zaid = "{:>{width}}".format(int(daughter_16), width = width_zaid)
        file.write(f"{aligned_orion}{aligned_zaid} Reaction  16 (n,2n)\n")
        # This will need to be generalised in future if multiple burnups to be used
        if args.reactor == 'LFR30':
            file.write(f"{get_cross_section(nuclide_ID, reaction = 16, burnup = 'single')}\n")
        elif args.reactor == 'LFR200':
            for burnup_step in ['BoC', 'EoC']:
                file.write(f"{get_cross_section(nuclide_ID, reaction = 16, burnup = burnup_step)}\n")

    if not np.isnan(daughter_17):
        daughter_orion_17 = calculate_daughter_orion(daughter_17, orion_dataframe)
        aligned_orion = "{:>{width}}".format(daughter_orion_17, width = width_orion)
        aligned_zaid = "{:>{width}}".format(int(daughter_17), width = width_zaid)
        file.write(f"{aligned_orion}{aligned_zaid} Reaction  17 (n,3n)\n")
        if args.reactor == 'LFR30':
            file.write(f"{get_cross_section(nuclide_ID, reaction = 17, burnup = 'single')}\n")
        elif args.reactor == 'LFR200':
            for burnup_step in ['BoC', 'EoC']:
                file.write(f"{get_cross_section(nuclide_ID, reaction = 17, burnup = burnup_step)}\n")
    
    if not np.isnan(daughter_18):
        daughter_orion_18 = -1 # ORION ID for fission product is -1
        aligned_orion = "{:>{width}}".format(daughter_orion_18, width = width_orion)
        aligned_zaid = "{:>{width}}".format(int(daughter_18), width = width_zaid)
        file.write(f"{aligned_orion}{aligned_zaid} Reaction  18 (n,fission)\n")
        if args.reactor == 'LFR30':
            file.write(f"{get_cross_section(nuclide_ID, reaction = 18, burnup = 'single')}\n")
        elif args.reactor == 'LFR200':
            for burnup_step in ['BoC', 'EoC']:
                file.write(f"{get_cross_section(nuclide_ID, reaction = 18, burnup = burnup_step)}\n")

    if not np.isnan(daughter_102):
        daughter_orion_102 = calculate_daughter_orion(daughter_102, orion_dataframe)
        aligned_orion = "{:>{width}}".format(daughter_orion_102, width = width_orion)
        aligned_zaid = "{:>{width}}".format(int(daughter_102), width = width_zaid)
        file.write(f"{aligned_orion}{aligned_zaid} Reaction 102 (n,gamma)\n")
        if args.reactor == 'LFR30':
            file.write(f"{get_cross_section(nuclide_ID, reaction = 102, burnup = 'single')}\n")
        elif args.reactor == 'LFR200':
            for burnup_step in ['BoC', 'EoC']:
                file.write(f"{get_cross_section(nuclide_ID, reaction = 102, burnup = burnup_step)}\n")

    if not np.isnan(daughter_103):
        daughter_orion_103 = calculate_daughter_orion(daughter_103, orion_dataframe)
        aligned_orion = "{:>{width}}".format(daughter_orion_103, width = width_orion)
        aligned_zaid = "{:>{width}}".format(int(daughter_103), width = width_zaid)
        file.write(f"{aligned_orion}{aligned_zaid} Reaction 103 (n,proton)\n")
        if args.reactor == 'LFR30':
            file.write(f"{get_cross_section(nuclide_ID, reaction = 103, burnup = 'single')}\n")
        elif args.reactor == 'LFR200':
            for burnup_step in ['BoC', 'EoC']:
                file.write(f"{get_cross_section(nuclide_ID, reaction = 103, burnup = burnup_step)}\n")

    return
# -------------  Main Code  ------------ #
args = parse_arguments()
# Define the Element Symbols. Eisteinium included to prevent list index iss ues when calculating the parent nuclide name of Cf.
element_symbols = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es']

# Reads in the csv file containing the data to be saved to the MPR file and assigns it to an array
file_dir = f'{args.input}/ZAID_results.csv'
df = pd.read_csv(file_dir, header = 0)

# Assigns each column of csv file to a variable
nuclide = df[df.columns[0]]
zaid = df[df.columns[1]]
orion_id = df[df.columns[2]]
NumberReactions = df[df.columns[3]]
daughter_16 = df[df.columns[4]]
daughter_17 = df[df.columns[5]]
daughter_18 = df[df.columns[6]]
daughter_102 = df[df.columns[7]]
daughter_103 = df[df.columns[8]]
NumberParents = df[df.columns[9]]
parent_16 = df[df.columns[10]]
parent_17 = df[df.columns[11]]
parent_102 = df[df.columns[12]]
parent_103 = df[df.columns[13]]
parent_16_1m = df[df.columns[14]]
parent_16_2m = df[df.columns[15]]
parent_17_1m = df[df.columns[16]]
parent_17_2m = df[df.columns[17]]
parent_102_1m = df[df.columns[18]]
parent_102_2m = df[df.columns[19]]
parent_103_1m = df[df.columns[20]]
parent_103_2m = df[df.columns[21]]

# Sort the dataframe so that it is in the order of ascending ORION IDs
df = df.sort_values(by = df.columns[2])

file_path = f'{args.path}/{args.reactor}_{args.enrichment_zone}_MPR.orion'

# Opens the excel file containing the ORION IDs for each nuclide and loads it in a pandas dataframe
ORION_ID_dir = f'{args.reaction_data}/orion_nuclides_list.xlsx'
df_ORION_ID = pd.read_excel(ORION_ID_dir, header = None)
# Columns in excel file do not contain headers so create them here
df_ORION_ID.columns = ['Nuclide Name', 'Buffer Mass']
nuclide_in_df = df_ORION_ID[df_ORION_ID.columns[0]]

# This needs to be modified if more burnup steps are to be used in future analyses
if args.reactor == 'LFR30':
    num_burnup_steps = 1
elif args.reactor == 'LFR200':
    num_burnup_steps = 2

width_nnuclides = 7
aligned_nuclide = "{:>{width}}".format(len(nuclide), width = width_nnuclides)
# Write the preliminary information into the MPR file
with open(file_path, "w") as file:
    file.write(f'BurnupSteps     {num_burnup_steps}\n')
    file.write('0\n')  # 0 MWd/tHM (fresh fuel, BoC)
    # This should be generalised in future for multiple bunrup steps
    if args.reactor == 'LFR200':
        file.write('70000\n') # EoC burnup
    file.write(f'NNuclides{aligned_nuclide}\n')
    file.write(f'CrossSections\n')
    file.write(f'MaxReactions    {max(NumberReactions)}\n') # Maximum number of reactions for a single nuclide, set to 5 as only interested in MT = 16, 17, 18, 102, 103

    # Iterates through each row of the pandas dataframe and writes the data into the .txt file in the MPR file format
    for i, row in df.iterrows():
        save_to_MPR(orion_id[i], 
                    zaid[i], 
                    nuclide[i], 
                    NumberReactions[i], 
                    daughter_16[i], 
                    daughter_17[i], 
                    daughter_18[i], 
                    daughter_102[i], 
                    daughter_103[i], 
                    NumberParents[i], 
                    parent_16[i], 
                    parent_17[i],
                    parent_102[i], 
                    parent_103[i],
                    parent_16_1m[i],
                    parent_17_1m[i],
                    parent_102_1m[i],
                    parent_103_1m[i],
                    parent_16_2m[i],
                    parent_17_2m[i],
                    parent_102_2m[i],
                    parent_103_2m[i],
                    nuclide_in_df)
        

