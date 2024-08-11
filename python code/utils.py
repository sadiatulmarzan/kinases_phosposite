import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import logomaker
from collections import Counter
import pandas as pd


all_sequences = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
def plotting_pssm(psp, johnson, common_family, title1 = 'PSP', title2 = 'Johnson'):

    for family in common_family:
        pssm1 = psp[family]
        pssm2 = johnson[family]

        # Add missing columns to the PSSMs
        pssm1 = add_missing_columns(pssm1, all_sequences)
        pssm2 = add_missing_columns(pssm2, all_sequences)

        pssm1['position'] = list(range(-7, 8))
        pssm1.set_index('position', inplace=True)

        pssm2['position'] = list(range(-7, 8))
        pssm2.set_index('position', inplace=True)

        pssm1 = pssm1.T.sort_index()
        pssm2 = pssm2.T.sort_index()  


        if pssm1.shape[1] != pssm2.shape[1]:
            print(family, 'PSP and Johson data have different lengths')
            print(f"Psp is {pssm1.shape[1]} and Johnson is {pssm2.shape[1]}")
            print(f"psp columns: {pssm1.columns} and johnson columns: {pssm2.columns}")
            print("----------------------------------------")
            continue

        # make the columns match
        try:
            # euclidean distance between the two PSSMs
            distances = np.sqrt(np.sum((pssm1 - pssm2) ** 2, axis=1))
            # Plotting the PSSMs and distances
            fig, axes = plt.subplots(1, 3, figsize=(15, 5))
            # set figure title 
            fig.suptitle(family)

            sns.heatmap(pssm1, ax=axes[0], cmap='viridis', cbar=True)
            axes[0].set_title(title1)

            sns.heatmap(pssm2, ax=axes[1], cmap='viridis', cbar=True)
            axes[1].set_title(title2)

            axes[2].plot(distances, marker='o')
            axes[2].set_title('Euclidean Distance per Position')
            axes[2].set_xlabel('Amino Acid')
            axes[2].set_ylabel('Euclidean Distance')

            # Set x-ticks to the row indices of pssm1 or pssm2
            row_indices = pssm1.index
            axes[2].set_xticks(range(len(row_indices)))
            axes[2].set_xticklabels(row_indices, rotation=45)

            plt.tight_layout()
            plt.show()
        except Exception as e:
            print(family, e)
    

def add_missing_columns(pssm1, all_sequences):
    for col in all_sequences:
        if col not in pssm1.columns:
            pssm1[col] = -1
    if '_' in pssm1.columns:
        pssm1.drop('_', axis=1, inplace=True)
    return pssm1



def generate_kinase_logos(df, need_family= None):
    """
    Generate sequence logos for each kinase in the provided dataframe.
    
    Parameters:
    df (pandas.DataFrame): DataFrame containing kinase data with columns 'ranked_1' and 'SITE_+/-7_AA'.
    need_family (set): Set of kinase families that need to be processed.
    
    Returns:
    dict: Dictionary where keys are kinase names and values are logomaker.Logo objects.
    """
    # Group by 'ranked_1' and combine all sequences into a list
    kinase_sequences = df.groupby('ranked_1')['SITE_+/-7_AA'].apply(list)

    # Initialize an empty dictionary to store the Logos for each kinase
    kinase_logos = {}

    # For each kinase, create a position frequency matrix and then a sequence logo
    for kinase, sequences in kinase_sequences.items():
        if need_family is not None and kinase not in need_family:
            continue

        # Combine sequences into a single string for each position
        aligned_sequences = [''.join(seq) for seq in zip(*sequences)]
        
        # Create a DataFrame where each row corresponds to one position
        position_df = pd.DataFrame([Counter(pos) for pos in aligned_sequences]).fillna(0)
        position_df['position'] = range(-7, 8)
        position_df.set_index('position', inplace=True)
        
        # Normalize the counts to get frequencies
        position_freq_matrix = position_df.div(position_df.sum(axis=0), axis=1)

        # Create the sequence logo
        logo = logomaker.Logo(position_freq_matrix, color_scheme='skylign_protein', stack_order='small_on_top')
        
        # Set title and axis labels
        logo.ax.set_title(f'Sequence Logo for {kinase} from Jonson')
        logo.ax.set_xlabel('Position')
        logo.ax.set_ylabel('Frequency')
        
        # Store the logo in the dictionary
        kinase_logos[kinase] = logo

    # Display the logos
    plt.show()
    
    return kinase_logos

def convert_to_normalize_matrix(df, columns_to_keep, suffix = ''): 
    # Flatten the DataFrame
    flattened_data = {}
    for family_name, df in psp.items():
        df = add_missing_columns(df)
        df['position'] = range(-7, 8)
        df.set_index('position', inplace=True)
        # sort the columns 
        df = df[sorted(df.columns)]
        flattened_df = df.unstack().reset_index()
        flattened_df['feature'] = flattened_df['position'].astype(str) + flattened_df['level_0']
        flattened_df = flattened_df[['feature', 0]].set_index('feature').T
        flattened_data[family_name] = flattened_df

    # Concatenate all families into a single DataFrame
    final_df = pd.concat(flattened_data.values(), keys=flattened_data.keys())

    # Reset index to have family names as a column
    final_df.reset_index(level=1, drop=True, inplace=True)
    final_df.reset_index(inplace=True)
    final_df.rename(columns={'index': 'Family'}, inplace=True)
    # add suffix to the Family column
    final_df['Family'] = final_df['Family'] + suffix

    final_df.set_index('Family', inplace=True)
    final_df = final_df[columns_to_keep]

    # drop the feature column 
    # final_df.reset_index( inplace=True)
    # final_df.set_index('Family', inplace=True)
    
    return final_df

def plotting_pssm_only_heatmap(psp, johnson, scope3p, common_family, title1 = 'PSP', title2 = 'Johnson', title3 = 'Scope3P'):   

    for family in common_family:
        pssm1 = psp[family]
        pssm2 = johnson[family]
        pssm3 = scope3p[family]

        # Add missing columns to the PSSMs
        pssm1 = add_missing_columns(pssm1, all_sequences)
        pssm2 = add_missing_columns(pssm2, all_sequences)
        pssm3 = add_missing_columns(pssm3, all_sequences)

        pssm1['position'] = list(range(-7, 8))
        pssm1.set_index('position', inplace=True)

        pssm2['position'] = list(range(-7, 8))
        pssm2.set_index('position', inplace=True)

        pssm3['position'] = list(range(-7, 8))
        pssm3.set_index('position', inplace=True)

        pssm1 = pssm1.T.sort_index()
        pssm2 = pssm2.T.sort_index()  
        pssm3 = pssm3.T.sort_index()



        # make the columns match
        try:
            # Plotting the PSSMs and distances
            fig, axes = plt.subplots(1, 3, figsize=(15, 5))
            # set figure title 
            fig.suptitle(family)

            sns.heatmap(pssm1, ax=axes[0], cmap='viridis', cbar=True)
            axes[0].set_title(title1)

            sns.heatmap(pssm2, ax=axes[1], cmap='viridis', cbar=True)
            axes[1].set_title(title2)

            sns.heatmap(pssm3, ax=axes[2], cmap='viridis', cbar=True)
            axes[2].set_title(title3)
            
            plt.tight_layout()
            plt.show()
        except Exception as e:
            print(family, e)