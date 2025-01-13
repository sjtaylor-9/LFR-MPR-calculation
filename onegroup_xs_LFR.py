"""
onegroup_xs_LFR.py

This script is responsible for calculating the one-group cross-sections for reactions with MT IDs 16, 17, 18, 102, and 103 in the ENDF, JEFF (3.3), and TENDL (2019) databases.
For a parent nuclide x and daughter nuclide(s) y, reaction 16 is x(n,2n)y, reaction 17 is x(n,3n)y, reaction 18 is x(n,fission)y, reaction 102 is x(n,gamma)y, and reaction 103 is x(n,p)y.
The code takes in parsers for the reactor, enrichment zone and burnup so that the relevant neutron spectra can be read in. Also directories for the public ENDF/JEFF/TENDL data, newcleo neutron spectra data, and the list of examined reactions are required. 

The one-group cross-section is calculated from the integral of the energy dependent cross-section * neutron flux at given energy w.r.t energy divided by the integrated neutron flux spectrum.
Since, the neutron spectra are grouped in energy bins, this integral is calculation is discretised as an approximation.
The code reads in the respetive list of ENDF, JEFF, or TENDL reactions and assigns them to a numpy array. It loops through each row of the array and for each reaction calculates the one-group cross-section and appends it to a new list along with the parent ZAID, reaction ID, and daughter ZAID.
Firstly, the energy-dependent cross-section is grouped into the same energy binning scheme as used for the neutron flux spectra. If a bin contains no data, then the script interpolates the cross-sections using data either side.
Then, the one-group cross-section for the bin is calculated using the defined equation.
The daughter nuclides are found by reading in the list of reactions with the NNL example MPR and extracting the daughter ZAID for the relevant reaction. This is done as the ENDF/JEFF/TENDL libraries do not specify if the daughter nuclide is metastable or not; hence, the only way to deduce this is to check to see if the daughter is metastable in the NNL list.
If the reaction does not exist in the NNL database, then the daughter nuclide is assumed to be a ground state nuclide and calculated manually.

The one-group cross-sections for all reactions are stored within a .csv file.

Authors: Michael Weekes (michael.weekes@newcleo.com) and Sam Taylor (sam.taylor@newcleo.com)
Last Edited: 12/01/2025
"""
# ---------  Import Libraries  --------- #
import numpy as np
from numpy import genfromtxt
import endf
import pandas as pd
import os
import argparse
# --------- Global Variables ----------- #
element_symbols = ['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh',	'Pd','Ag','Cd','In','Sn','Sb','Te',	'I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No']
# Set up units (working in CGS system)
barn = 1.e-24 #cm^2
# MT IDs of reactions under consideration
reaction_ids = [16, 17, 18, 102, 103]
reaction_strs = ["(n,2n)", "(n,3n)", "(n,fission)", "(n,$\gamma$)", "(n,p)"]
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
    --cross_section_data    Used to specify the directory in which the public data regarding the ENDF, JEFF and TENDL reaction cross-sections can be found.
    --reaction_data         Used to specify the directory in which the lists of ENDF, JEFF and TENDL reactions to be examined can be found.
    --burnup                This specifies which burnup step is being examined. This is necessary so that the correct neutron spectrum can be found and the burnup of the output cross-sections is specified.
                            For the LFR200 the argument must be one of BoC, EoC, BoL, EoL, whereby BoC is 'Beginning of Cycle', BoL is 'Beginning of Life', and E represents the respective end burnup step.
                            For the LFR30, currently only a single burnup step is available.
    
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
        "--newcleo_input",
        type=dir_path,
        required=True,
        help="flag to set the path where internal newcleo input data should be found"
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
        "--burnup",
        type = str,
        choices = ['BoL', 'EoL', 'single', 'BoC', 'EoC'],
        required = True,
        help = 'flag to set which burnup step is to be examined'
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

def get_neutron_spectrum():
    """
    This function finds the relative neutron spectrum data based on the arguments passed in argparse.
    The neutron spectrum data contains the edges of the energy bins (eV for LFR30 and MeV for LFR200) and the neutron fluxes (n/cm^2/s) within these energy bins.

    Returns:
        neutronspec_bin_edges (numpy array): A list containing the edges of the neutron energy (eV) binning scheme used within the neutron flux spectara. The values are the lower edges of the bins apart from the final array element which is the final bin edge (20 MeV).
        num_energy_bins (integer): The number of energy bins used in the binning scheme.
        central_value_neutron_flux (numpy array): A list containing the neutron fluxes within each bin.
    """
    if args.reactor == 'LFR30':
        # csv file with the edges of the energy groups
        neutronspec_bin_edges = genfromtxt(f'{args.newcleo_input}/LFR30_neutronspec_groupedges.csv', delimiter = ',')
        # group cross-section data according to the neutron spec group structure
        central_value_neutron_flux = genfromtxt(f'{args.newcleo_input}/LFR30_neutronspec_values.csv', delimiter = ',')
    elif args.reactor == 'LFR200':
        if args.enrichment_zone == 'inner':
            neutronspec_bin_edges = 1e6 * pd.read_excel(f'{args.newcleo_input}/Flux_LFR_AS_200_PHASE_2_04_187_211_246.xlsx',
                          sheet_name = 'INNER FUEL',
                          usecols = ['Emin'],
                          header = 2,
                          nrows = 252)
            neutronspec_bin_edges = neutronspec_bin_edges.to_numpy() # Converts pd dataframe to numpy array
            neutronspec_bin_edges = neutronspec_bin_edges.reshape(-1) # Ensures the numpy array is 1D
            # Find the last bin edge in the Emax column of the spreadsheet. This is needed as pandas groupby() requires the final bin edge
            e_max_values = 1e6 * pd.read_excel(f'{args.newcleo_input}/Flux_LFR_AS_200_PHASE_2_04_187_211_246.xlsx',
                          sheet_name = 'OUTER FUEL',
                          usecols = ['Emax'],
                          header = 2)
            maximum_e_max = e_max_values.max()
            neutronspec_bin_edges = np.append(neutronspec_bin_edges, maximum_e_max)
            central_value_neutron_flux = pd.read_excel(f'{args.newcleo_input}/Flux_LFR_AS_200_PHASE_2_04_187_211_246.xlsx',
                          sheet_name = 'INNER FUEL',
                          usecols = [args.burnup],
                          header = 2,
                          nrows = 252)
            central_value_neutron_flux = central_value_neutron_flux.to_numpy()
        elif args.enrichment_zone == 'middle':
            neutronspec_bin_edges = 1e6* pd.read_excel(f'{args.newcleo_input}/Flux_LFR_AS_200_PHASE_2_04_187_211_246.xlsx',
                          sheet_name = 'MIDDLE FUEL',
                          usecols = ['Emin'],
                          header = 2,
                          nrows = 252)
            neutronspec_bin_edges = neutronspec_bin_edges.to_numpy()
            neutronspec_bin_edges = neutronspec_bin_edges.reshape(-1)
            e_max_values = 1e6 * pd.read_excel(f'{args.newcleo_input}/Flux_LFR_AS_200_PHASE_2_04_187_211_246.xlsx',
                          sheet_name = 'MIDDLE FUEL',
                          usecols = ['Emax'],
                          header = 2)
            maximum_e_max = e_max_values.max()
            neutronspec_bin_edges = np.append(neutronspec_bin_edges, maximum_e_max)
            central_value_neutron_flux = 1e6 * pd.read_excel(f'{args.newcleo_input}/Flux_LFR_AS_200_PHASE_2_04_187_211_246.xlsx',
                          sheet_name = 'MIDDLE FUEL',
                          usecols = [args.burnup],
                          header = 2,
                          nrows = 252)
            central_value_neutron_flux = central_value_neutron_flux.to_numpy()
        elif args.enrichment_zone == 'outer':
            neutronspec_bin_edges = 1e6 * pd.read_excel(f'{args.newcleo_input}/Flux_LFR_AS_200_PHASE_2_04_187_211_246.xlsx',
                          sheet_name = 'OUTER FUEL',
                          usecols = ['Emin'],
                          header = 2,
                          nrows = 252)
            neutronspec_bin_edges = neutronspec_bin_edges.to_numpy()
            neutronspec_bin_edges = neutronspec_bin_edges.reshape(-1)
            e_max_values = 1e6 * pd.read_excel(f'{args.newcleo_input}/Flux_LFR_AS_200_PHASE_2_04_187_211_246.xlsx',
                          sheet_name = 'OUTER FUEL',
                          usecols = ['Emax'],
                          header = 2)
            maximum_e_max = e_max_values.max()
            neutronspec_bin_edges = np.append(neutronspec_bin_edges, maximum_e_max)
            central_value_neutron_flux = pd.read_excel(f'{args.newcleo_input}/Flux_LFR_AS_200_PHASE_2_04_187_211_246.xlsx',
                          sheet_name = 'OUTER FUEL',
                          usecols = [args.burnup],
                          header = 2,
                          nrows = 252)
            central_value_neutron_flux = central_value_neutron_flux.to_numpy()

    num_energy_bins = len(neutronspec_bin_edges) - 1

    return neutronspec_bin_edges, num_energy_bins, central_value_neutron_flux

def get_jeff_reaction_xs(bin_edges, num_bins, neutron_flux, nnl_parents, nnl_reactions, nnl_daughters):
    """
    This function finds the JEFF 3.3 reactions to be used within the MPR file and outputs the one-group cross-section for each of these reactions.
    The energy dependent cross-sections are stored within the .ACE files and so are read in from the user's local storage using the cross_section_data parser.
    The LFR200 core will be approximately 540 C so the 900 K .ACE files are a good approximation.

    Args:
        bin_edges (numpy array): A list containing the edges of the neutron energy (eV) binning scheme used within the neutron flux spectara. The values are the lower edges of the bins apart from the final array element which is the final bin edge (20 MeV).
        num_bins (integer): The number of energy bins used in the binning scheme.
        neutron_flux (numpy array): A list containing the neutron fluxes within each bin.

    Returns:
        jeff_cross_sections (array): A 2D array containing the parent ZAID, reaction ID, daughter ZAID, and one-group cross-section.
    """
    # Select element, isotope and reaction (by MT ID)
    lib_dir = f'{args.cross_section_data}/ace_900/ace_900/'
    jeff_reactions = genfromtxt(f'{args.reaction_data}/JEFF_reactions.csv', delimiter = ',')
    n_reactions = len(jeff_reactions[:, 0])
    jeff_cross_sections = []

    # Get Z from element symbol
    for i in range(n_reactions):

        parent_ZAI = jeff_reactions[i, 0]

        test_reaction_id = jeff_reactions[i, 1]

        meta_flag = np.mod(parent_ZAI, 10)
        #print(meta_flag)
        #print("%i, %s, %i" % (parent_ZAI,ZAI_str,meta_flag))
        parent_ZA = parent_ZAI - meta_flag

        parent_Z = int(parent_ZA / 10000)
        parent_A = int(((parent_ZA) - (10000 * parent_Z)) / 10)

        if(meta_flag == 1): meta_str = "m"
        else: meta_str = "g"

        element = element_symbols[parent_Z - 1]
        filename = "%i-%s-%i%s-900.ace" % (parent_Z, element, parent_A, meta_str)
        infile_str = lib_dir + filename
        table = endf.ace.get_table(infile_str)
        test_reaction = (table.interpret().reactions[test_reaction_id].xs)['900K']

        d1 = {'Energy': test_reaction.x, 'Sigma': test_reaction.y}
        test_df = pd.DataFrame(data = d1)

        # group cross-section data according to the neutron spec group structure
        grouped_xs_df = test_df.groupby(pd.cut(test_df["Energy"], bin_edges), observed = False).mean()
        grouped_xs_df = grouped_xs_df.Sigma.reset_index()
        grouped_xs = grouped_xs_df["Sigma"].copy()

        # if all NaN (i.e. no overlap between xs data and neutronspec data) skip reaction
        if(sum(np.isnan(grouped_xs)) == len(grouped_xs)): continue

        # if last value is NaN - replace with zeros, count back and replace with zero until a non-NaN value is found
        if(np.isnan(grouped_xs[len(grouped_xs) - 1]) == 1):
            stop_flag = 0
            index = len(grouped_xs) - 1
            print(index)
            while(stop_flag == 0):
                grouped_xs[index] = 0
                if(np.isnan(grouped_xs[index - 1]) == 1):
                    index -= 1
                else:
                    stop_flag = 1

        # if first value is NaN - replace with zeros, count forward and replace with zero until a non-NaN value is found
        if(np.isnan(grouped_xs[0]) == 1):
            stop_flag = 0
            index = 0
            while(stop_flag == 0):
                grouped_xs[index] = 0
                if(np.isnan(grouped_xs[index + 1]) == 1):
                    index += 1
                else:
                    stop_flag = 1


        # check if any NaN values remain
        any_nan = 0
        if(np.sum(np.isnan(grouped_xs)) > 0): any_nan = 1

        while(any_nan == 1):
            # count through NaN values, replace each with interpolation from the two closest values
            for i in range(1, len(grouped_xs) - 1):
                if np.isnan(grouped_xs[i]):
                    lb = 1
                    while((np.isnan(grouped_xs[i - lb]) == 1) & ((i - lb) > 1)):
                        lb += 1
                    ub = 1
                    while((np.isnan(grouped_xs[i + ub]) == 1) & ((i + ub)<len(grouped_xs) - 1)):
                        ub += 1
                    grouped_xs[i] = 0.5 * (grouped_xs[i - lb] + grouped_xs[i + ub])
            if(np.sum(np.isnan(grouped_xs)) == 0): any_nan = 0

        # Get daughter ZAIs from NNL MPR file
        nnl_found = 0
        for i in range(len(nnl_parents)):
            if (nnl_found == 1): break
            if ((nnl_parents[i] == parent_ZAI) & (nnl_reactions[i] == test_reaction_id)):
                nnl_found = 1
                daughter_ZAI = nnl_daughters[i]
        
        prod_sum = 0
        flux_sum = 0
        for i in range(len(grouped_xs)):
            if(np.isnan(grouped_xs[i])): grouped_xs[i] = 0
            prod_sum += grouped_xs[i] * neutron_flux[i]
            flux_sum += neutron_flux[i]

        # If the reaction doesn't exist in the NNL database then the daughter ZAID is calculated manually 
        if(nnl_found == 0):
            if test_reaction_id == 16:
                daughter_Z = parent_Z
                daughter_A = parent_A - 1
            if test_reaction_id == 17:
                daughter_Z = parent_Z
                daughter_A = parent_A - 2
            if test_reaction_id == 102:
                daughter_Z = parent_Z
                daughter_A = parent_A + 1
            if test_reaction_id == 103:
                daughter_Z = parent_Z - 1
                daughter_A = parent_A
            if test_reaction_id == 18:
                daughter_Z = 0
                daughter_A = 0
            daughter_ZAI = (daughter_Z * 10000) + (daughter_A * 10)
        #print("%i,%i,%i,%.5e" % (parent_ZAI, test_reaction_id, daughter_ZAI, (prod_sum/flux_sum)))
        jeff_cross_sections.append([parent_ZAI, test_reaction_id, daughter_ZAI, (prod_sum / flux_sum)])

    return jeff_cross_sections

def get_endf_reaction_xs(bin_edges, num_bins, neutron_flux, nnl_parents, nnl_reactions, nnl_daughters):
    """
    This function finds the ENDF reactions that are not in JEFF 3.3 to be used within the MPR file and outputs the one-group cross-section for each of these reactions.
    The energy dependent cross-sections are stored within the .ACE files and so are read in from the user's local storage using the cross_section_data parser.
    The LFR200 core will be approximately 540 C so the 900 K .ACE files are a good approximation.

    Args:
        bin_edges (numpy array): A list containing the edges of the neutron energy (eV) binning scheme used within the neutron flux spectara. The values are the lower edges of the bins apart from the final array element which is the final bin edge (20 MeV).
        num_bins (integer): The number of energy bins used in the binning scheme.
        neutron_flux (numpy): A list containing the neutron fluxes within each bin.

    Returns:
        endf_cross_sections (array): A 2D array containing the parent ZAID, reaction ID, daughter ZAID, and one-group cross-section.
    """
    # Select element, isotope and reaction (by MT ID)
    lib_dir = f'{args.cross_section_data}/Lib80x/Lib80x/'
    missing_reactions = genfromtxt(f'{args.reaction_data}/reactions_in_endf_not_jeff.csv', delimiter = ',')
    n_reactions = len(missing_reactions[:, 0])
    endf_cross_sections = []
    
    # Get Z from element symbol
    for i in range(n_reactions):
        parent_ZAI = missing_reactions[i, 0]

        test_reaction_id = missing_reactions[i, 1]

        meta_flag = np.mod(parent_ZAI, 10)
        #print("%i, %s, %i" % (parent_ZAI,ZAI_str,meta_flag))
        file_ZA = int(parent_ZAI / 10)
        parent_Z = int(parent_ZAI / 10000)
        parent_A = int((parent_ZAI - (10000 * parent_Z)) / 10)

        if(meta_flag == 1): meta_str = "m1_"
        else: meta_str = ""

        element = element_symbols[parent_Z - 1]

        filename = "%s/%s%i.802nc" % (element, meta_str, file_ZA)
        infile_str = lib_dir + filename
        table = endf.ace.get_table(infile_str)
        test_reaction = (table.interpret().reactions[test_reaction_id].xs)['900K']

        d1 = {'Energy': test_reaction.x, 'Sigma': test_reaction.y}
        test_df = pd.DataFrame(data = d1)

        # group cross-section data according to the neutron spec group structure
        grouped_xs_df = test_df.groupby(pd.cut(test_df["Energy"], bin_edges), observed = False).mean()
        grouped_xs_df = grouped_xs_df.Sigma.reset_index()
        grouped_xs = grouped_xs_df["Sigma"].copy()

         # if all NaN (i.e. no overlap between xs data and neutronspec data) skip reaction
        if(sum(np.isnan(grouped_xs)) == len(grouped_xs)): continue

        # if last value is NaN - replace with zeros, count back and replace with zero until a non-NaN value is found
        if(np.isnan(grouped_xs[len(grouped_xs) - 1]) == 1):
            stop_flag = 0
            index = len(grouped_xs) - 1
            print(index)
            while(stop_flag == 0):
                grouped_xs[index] = 0
                if(np.isnan(grouped_xs[index - 1]) == 1):
                    index -= 1
                else:
                    stop_flag = 1

        # if first value is NaN - replace with zeros, count forward and replace with zero until a non-NaN value is found
        if(np.isnan(grouped_xs[0]) == 1):
            stop_flag = 0
            index = 0
            while(stop_flag == 0):
                grouped_xs[index] = 0
                if(np.isnan(grouped_xs[index + 1]) == 1):
                    index += 1
                else:
                    stop_flag = 1


        # check if any NaN values remain
        any_nan = 0
        if(np.sum(np.isnan(grouped_xs)) > 0): any_nan = 1

        while(any_nan == 1):
            # count through NaN values, replace each with interpolation from the two closest values
            for i in range(1, len(grouped_xs) - 1):
                if np.isnan(grouped_xs[i]):
                    lb = 1
                    while((np.isnan(grouped_xs[i - lb]) == 1) & ((i - lb) > 1)):
                        lb += 1
                    ub = 1
                    while((np.isnan(grouped_xs[i + ub]) == 1) & ((i + ub) < len(grouped_xs) - 1)):
                        ub += 1
                    grouped_xs[i] = 0.5 * (grouped_xs[i - lb] + grouped_xs[i + ub])
            if(np.sum(np.isnan(grouped_xs)) == 0): any_nan = 0   

        # Get daughter ZAIs from NNL MPR file
        nnl_found = 0
        for i in range(len(nnl_parents)):
            if (nnl_found == 1): break
            if ((nnl_parents[i] == parent_ZAI) & (nnl_reactions[i] == test_reaction_id)):
                nnl_found = 1
                daughter_ZAI = nnl_daughters[i]
        

        prod_sum = 0
        flux_sum = 0
        for i in range(len(grouped_xs)):
            if(np.isnan(grouped_xs[i])): grouped_xs[i] = 0
            prod_sum += grouped_xs[i] * neutron_flux[i]
            flux_sum += neutron_flux[i]

        # If the reaction doesn't exist in the NNL database then the daughter ZAID is calculated manually 
        if(nnl_found == 0):
            if test_reaction_id == 16:
                daughter_Z = parent_Z
                daughter_A = parent_A - 1
            if test_reaction_id == 17:
                daughter_Z = parent_Z
                daughter_A = parent_A - 2
            if test_reaction_id == 102:
                daughter_Z = parent_Z
                daughter_A = parent_A + 1
            if test_reaction_id == 103:
                daughter_Z = parent_Z - 1
                daughter_A = parent_A
            if test_reaction_id == 18:
                daughter_Z = 0
                daughter_A = 0
            daughter_ZAI = (daughter_Z * 10000) + (daughter_A * 10)
        #print("%i,%i,%i,%.5e" % (parent_ZAI, test_reaction_id, daughter_ZAI, (prod_sum / flux_sum)))
        endf_cross_sections.append([parent_ZAI, test_reaction_id, daughter_ZAI, (prod_sum / flux_sum)])

    return endf_cross_sections

def get_tendl_reaction_xs(bin_edges, num_bins, neutron_flux, nnl_parents, nnl_reactions, nnl_daughters):
    """
    This function finds the TENDL (2019) reactions that are not in JEFF 3.3 or ENDF to be used within the MPR file and outputs the one-group cross-section for each of these reactions.
    The energy dependent cross-sections are stored within the .ACE files and so are read in from the user's local storage using the cross_section_data parser.

    Args:
        bin_edges (numpy array): A list containing the edges of the neutron energy (eV) binning scheme used within the neutron flux spectara. The values are the lower edges of the bins apart from the final array element which is the final bin edge (20 MeV).
        num_bins (integer): The number of energy bins used in the binning scheme.
        neutron_flux (numpy array): A list containing the neutron fluxes within each bin.

    Returns:
        tendl_cross_sections (array): A 2D array containing the parent ZAID, reaction ID, daughter ZAID, and one-group cross-section.
    """
    # Select element, isotope and reaction (by MT ID)
    lib_dir = f'{args.cross_section_data}/tendl19c/tendl19c/'
    tendl_reactions = genfromtxt(f'{args.reaction_data}/TENDL_reactions.csv', delimiter = ',')
    n_reactions = len(tendl_reactions[:, 0])
    tendl_cross_sections = []
    
    # Get Z from element symbol
    for i in range(n_reactions):
        parent_ZAI = tendl_reactions[i, 0]

        test_reaction_id = tendl_reactions[i, 1]

        meta_flag = np.mod(parent_ZAI, 10)
        #print("%i, %s, %i" % (parent_ZAI,ZAI_str,meta_flag))
        file_ZA = int(parent_ZAI / 10)
        parent_Z = int(parent_ZAI / 10000)
        parent_A = int((parent_ZAI - (10000 * parent_Z)) / 10)

        if(meta_flag == 0): meta_str = ""
        else: meta_str = "m"

        if (parent_A >= 100): A_str="%i" % parent_A
        elif parent_A >= 10: A_str="0%i" % parent_A
        else: A_str = "00%i" % parent_A

        element = element_symbols[parent_Z - 1]

        filename="%s%s%s" % (element, A_str, meta_str)
        infile_str = lib_dir + filename
        table = endf.ace.get_table(infile_str)
        test_reaction = (table.interpret().reactions[test_reaction_id].xs)['294K']

        d1 = {'Energy': test_reaction.x, 'Sigma': test_reaction.y}
        test_df = pd.DataFrame(data = d1)

        # group cross-section data according to the neutron spec group structure
        grouped_xs_df = test_df.groupby(pd.cut(test_df["Energy"], bin_edges), observed = False).mean()
        grouped_xs_df = grouped_xs_df.Sigma.reset_index()
        grouped_xs = grouped_xs_df["Sigma"].copy()

         # if all NaN (i.e. no overlap between xs data and neutronspec data) skip reaction
        if(sum(np.isnan(grouped_xs)) == len(grouped_xs)): continue

        # if last value is NaN - replace with zeros, count back and replace with zero until a non-NaN value is found
        if(np.isnan(grouped_xs[len(grouped_xs) - 1]) == 1):
            stop_flag = 0
            index = len(grouped_xs) - 1
            print(index)
            while(stop_flag == 0):
                grouped_xs[index] = 0
                if(np.isnan(grouped_xs[index - 1]) == 1):
                    index -= 1
                else:
                    stop_flag = 1

        # if first value is NaN - replace with zeros, count forward and replace with zero until a non-NaN value is found
        if(np.isnan(grouped_xs[0]) == 1):
            stop_flag = 0
            index = 0
            while(stop_flag == 0):
                grouped_xs[index] = 0
                if(np.isnan(grouped_xs[index + 1]) == 1):
                    index += 1
                else:
                    stop_flag = 1


        # check if any NaN values remain
        any_nan = 0
        if(np.sum(np.isnan(grouped_xs)) > 0): any_nan = 1

        while(any_nan == 1):
            # count through NaN values, replace each with interpolation from the two closest values
            for i in range(1, len(grouped_xs) - 1):
                if np.isnan(grouped_xs[i]):
                    lb = 1
                    while((np.isnan(grouped_xs[i - lb]) == 1) & ((i - lb) > 1)):
                        lb += 1
                    ub = 1
                    while((np.isnan(grouped_xs[i + ub]) == 1) & ((i + ub) < len(grouped_xs) - 1)):
                        ub += 1
                    grouped_xs[i] = 0.5 * (grouped_xs[i - lb] + grouped_xs[i + ub])
            if(np.sum(np.isnan(grouped_xs)) == 0): any_nan = 0   

        # Get daughter ZAIs from NNL MPR file
        nnl_found = 0
        for i in range(len(nnl_parents)):
            if (nnl_found == 1): break
            if ((nnl_parents[i] == parent_ZAI) & (nnl_reactions[i] == test_reaction_id)):
                nnl_found = 1
                daughter_ZAI = nnl_daughters[i]

        prod_sum = 0
        flux_sum = 0
        for i in range(len(grouped_xs)):
            if(np.isnan(grouped_xs[i])): grouped_xs[i] = 0
            prod_sum += grouped_xs[i] * neutron_flux[i]
            flux_sum += neutron_flux[i]

        # If the reaction doesn't exist in the NNL database then the daughter ZAID is calculated manually 
        if(nnl_found == 0):
            if test_reaction_id == 16:
                daughter_Z = parent_Z
                daughter_A = parent_A - 1
            if test_reaction_id == 17:
                daughter_Z = parent_Z
                daughter_A = parent_A - 2
            if test_reaction_id == 102:
                daughter_Z = parent_Z
                daughter_A = parent_A + 1
            if test_reaction_id == 103:
                daughter_Z = parent_Z - 1
                daughter_A = parent_A
            if test_reaction_id == 18:
                daughter_Z = 0
                daughter_A = 0
            daughter_ZAI = (daughter_Z * 10000) + (daughter_A * 10)
        #print("%i,%i,%i,%.5e" % (parent_ZAI, test_reaction_id, daughter_ZAI, (prod_sum / flux_sum)))
        tendl_cross_sections.append([parent_ZAI, test_reaction_id, daughter_ZAI, (prod_sum / flux_sum)])

    return tendl_cross_sections

def save_results(jeff_xs, endf_xs, tendl_xs):
    """
    This function saves the output arrays storing the respective one-group cross-sections to .csv files. A singular .csv file containing all reactions is also created with the list in order of parent ZAID.

    Args:
        jeff_xs (array): A 2D array containing the parent ZAID, reaction ID, daughter ZAID, and one-group cross-section for the JEFF reactions.
        endf_xs (array): A 2D array containing the parent ZAID, reaction ID, daughter ZAID, and one-group cross-section for the ENDF reactions.
        tendl_xs (array): A 2D array containing the parent ZAID, reaction ID, daughter ZAID, and one-group cross-section for the TENDL reactions.
    """
    # Convert the cross-sections to standard form with 2 d.p.
    for i in range(len(endf_xs)):
        xs = endf_xs[i][3]
        formatted_xs = np.format_float_scientific(xs, precision = 2)
        endf_xs[i][3] = formatted_xs
    for i in range(len(jeff_xs)):
        xs = jeff_xs[i][3]
        formatted_xs = np.format_float_scientific(xs, precision = 2)
        jeff_xs[i][3] = formatted_xs
    for i in range(len(tendl_xs)):
        xs = tendl_xs[i][3]
        formatted_xs = np.format_float_scientific(xs, precision = 2)
        tendl_xs[i][3] = formatted_xs

    # Convert numpy arrays to pandas dataframes
    jeff_data_df = pd.DataFrame(jeff_xs)
    endf_data_df = pd.DataFrame(endf_xs)
    tendl_data_df = pd.DataFrame(tendl_xs)

    # Format the dataframes such that first 2 columns don't have .0 decimals
    jeff_data_df[0] = jeff_data_df[0].apply(lambda x: f'{x: .0f}')
    jeff_data_df[1] = jeff_data_df[1].apply(lambda x: f'{x: .0f}')
    endf_data_df[0] = endf_data_df[0].apply(lambda x: f'{x: .0f}')
    endf_data_df[1] = endf_data_df[1].apply(lambda x: f'{x: .0f}')
    tendl_data_df[0] = tendl_data_df[0].apply(lambda x: f'{x: .0f}')
    tendl_data_df[1] = tendl_data_df[1].apply(lambda x: f'{x: .0f}')
    
    # Create directories for output files
    jeff_file = f'{args.path}/{args.reactor}_{args.enrichment_zone}_{args.burnup}_JEFF_reactiondata.csv'
    endf_file = f'{args.path}/{args.reactor}_{args.enrichment_zone}_{args.burnup}_ENDF_reactiondata.csv'
    tendl_file = f'{args.path}/{args.reactor}_{args.enrichment_zone}_{args.burnup}_TENDL_reactiondata.csv'

    # Outputs the JEFF and ENDF reaction cross-sections to .csv files
    jeff_data_df.to_csv(jeff_file, index = False, header = False)
    endf_data_df.to_csv(endf_file, index = False, header = False)
    tendl_data_df.to_csv(tendl_file, index = False, header = False)

    # Creates a single numpy array with both the JEFF and ENDF reaction cross-sections and converts this to a pandas dataframe
    all_reactions = np.vstack((jeff_xs, endf_xs, tendl_xs))
    all_reactions_df = pd.DataFrame(all_reactions)

    # Sort the dataframe so that it is in the order of ascending ZAIDs
    all_reactions_df = all_reactions_df.sort_values(by = all_reactions_df.columns[0])
    all_reactions_file = f'{args.path}/{args.reactor}_{args.enrichment_zone}_{args.burnup}_all_reactiondata.csv'
    all_reactions_df.to_csv(all_reactions_file, index = False, header = False)

    return
# -------------  Main Code  ------------ #
args = parse_arguments()

NNL_file_reactions = genfromtxt(f'{args.reaction_data}/NNLMPR_allreactions.csv', delimiter =',')
NNL_parents = NNL_file_reactions[1:, 0]
NNL_reactionids = NNL_file_reactions[1:, 1]
NNL_daughters = NNL_file_reactions[1:, 2]

neutronspec_edges, n_energy_bins, neutronspec_central = get_neutron_spectrum()
print(f"Calculating the cross-sectional data for the JEFF reactions")
jeff_data = get_jeff_reaction_xs(neutronspec_edges, n_energy_bins, neutronspec_central, NNL_parents, NNL_reactionids, NNL_daughters)
print(f"Calculating the cross-sectional data for the ENDF reactions")
endf_data = get_endf_reaction_xs(neutronspec_edges, n_energy_bins, neutronspec_central, NNL_parents, NNL_reactionids, NNL_daughters)
print(f"Calculating the cross-sectional data for the TENDL reactions")
tendl_data = get_tendl_reaction_xs(neutronspec_edges, n_energy_bins, neutronspec_central, NNL_parents, NNL_reactionids, NNL_daughters)

save_results(jeff_data, endf_data, tendl_data)
print(f"Calculated the cross-sectional data for the {args.reactor} in the enrichment zone: {args.enrichment_zone} at the {args.burnup} burnup step")
