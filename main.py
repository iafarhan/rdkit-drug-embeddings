# i've set drugs list by default to ["sertraline", "chloroquine"] but it can work with any compound with SMILES rep

import numpy as np
import pubchempy as pcp
from rdkit import Chem
from descriptastorus.descriptors import MakeGenerator
import argparse

def canonicalize_smiles(smiles):
    if smiles:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return Chem.MolToSmiles(
            mol,
            isomericSmiles=True,
            canonical=True,
            doRandom=False,
            allBondsExplicit=False,
            allHsExplicit=False,
            kekuleSmiles=False
        )
    else:
        return None

def embed_smile(smile):
    generator = MakeGenerator(["RDKit2D"])
    result = generator.process(smile)
    if result is None or result[0] is False:
        print(f"Failed to process SMILES: {smile}")
        return None
    return result[1:]

def get_smiles_from_name(drug_name):
    try:
        compounds = pcp.get_compounds(drug_name, 'name')
        if not compounds:
            raise ValueError(f"no SMILES found for {drug_name}!!")
            return None
        smiles = compounds[0].canonical_smiles
        canon_smiles = canonicalize_smiles(smiles)
        if canon_smiles:
            return canon_smiles
        else:
            print(f"Invalid SMILES for {drug_name}: {smiles}")
            return None
    except Exception as e:
        print(f"Error fetching SMILES for {drug_name}: {e}")
        return None

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate RDKit2D embeddings from drug names")
    parser.add_argument(
        "--drugs",
        nargs="+",
        default=["sertraline", "chloroquine"],
        help="List of drug names (default: sertraline, chloroquine)"
    )
    args = parser.parse_args()

    drug_names = args.drugs
    for name in drug_names:
        smiles = get_smiles_from_name(name)
        if smiles:
            embedding = embed_smile(smiles)
            if embedding is not None:
                embedding = np.array(embedding)
                print(f"{name} -> SMILES: {smiles}")
                print(f"{name} -> embedding shape: {embedding.shape}")
                print(f"{name} -> embeddings: {embedding[0:10]}")
                print("\n\n")
