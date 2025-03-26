### `canonicalize_smiles`

- Standardizes SMILES strings by parsing them into RDKit molecules and re-encoding them with `Chem.MolToSmiles` to make sure chemically identical molecules (e.g., sertraline variants) yield the same embeddings.

---

### `embed_smile(smile)`

- Computes a 200-d embeddings for a single SMILES string that maps molecular features relevant to drug activity by using RDKit2D descriptors. This embedding is used in downstream analysis in scFATE.

- For a **combination of drugs with dosages**, we compute the embedding as a weighted sum:

  **Embedding_combo = Σ (dosageᵢ × Embedding_drugᵢ)**

  where `dosageᵢ` is the normalized dosage in μM of drug *i*, and `Embedding_drugᵢ` is the 200-d embedding of that drug.

---

### `get_smiles_from_name(drug_name)`

- Automates retrieval of a drug's SMILES string from PubChem. The resulting SMILES is then canonicalized using `canonicalize_smiles`.

---

### Post-processing

- Not detailed here, but we perform a series of post-processing steps, which typically include:
  - Dropping descriptors with standard deviation ≤ 0.01
  - Normalizing descriptors using z-score normalization
