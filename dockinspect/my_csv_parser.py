"""
Parses prediction and residue data from CSV files for use in pocket and interaction analysis.

Functions:
    get_df: Reads a CSV file and returns a cleaned DataFrame.
    parse_predictions: Converts a prediction CSV to a pocket-to-residues dictionary.
    parse_residues: Converts a residue CSV to a residue-location-to-amino-acid dictionary.
"""

import pandas as pd

def get_df(predictions_file: str) -> pd.DataFrame:
    """
    Parses a CSV file containing PDB predictions to generate an according dataframe.

    :param predictions_file: Path to the CSV file containing PDB predictions.
    :return: Dataframe from the given CSV file.
    """
    df = pd.read_csv(predictions_file)
    df.columns = df.columns.str.strip()
    return df

def parse_predictions(predictions_file: str) -> dict:
    """
    Parses a CSV file containing predictions to generate a dictionary mapping pocket IDs 
    to lists of residue IDs.

    :param predictions_file: Path to the CSV file containing PDB predictions.
    :return: A dictionary where keys are pocket IDs and values are lists of residue IDs.
    """
    
    predictions_df = get_df(predictions_file)
    pocket_locations_dict = {}

    for index, row in predictions_df.iterrows():
        pocket_name = row["name"].strip() 
        residue_ids = row["residue_ids"].split()  
        pocket_locations_dict[pocket_name] = residue_ids

    return pocket_locations_dict

def parse_residues(residues_file: str) -> dict:
    """
    Parses a CSV file containing PDB residue information and generates a dictionary mapping  each chain and residue ID 
    to the corresponding amino acid three-letter abbreviation.

    :param residues_file: Path to the CSV file containing PDB residue data.
    :return: A dictionary where keys are residue IDs and values are the three-letter amino acid abbreviations at those locations.
    """

    residues_df = get_df(residues_file)
    location_aa_dict = {}

    for _, row in residues_df.iterrows():
        chain = row["chain"].strip()  
        residue_id = row["residue_label"]  
        amino_acid = row["residue_name"].strip()  

        location_aa_dict[f"{chain}_{residue_id}"] = amino_acid

    return location_aa_dict