from rdkit import Chem

# Stereo-aware SMILES for Artemether
smiles = "CO[C@H]1C[C@@H]2[C@@H]3CCC(C)[C@H](O2)[C@@]4(OOC(C)C)[C@H]3CC[C@@H]14"

# Convert SMILES to molecule
mol = Chem.MolFromSmiles(smiles)

if mol is None:
    print("Invalid SMILES string!")
    exit()

# Add hydrogens
mol = Chem.AddHs(mol)

# Assign stereochemistry
Chem.AssignStereochemistry(mol, force=True, cleanIt=True)

# Find chiral centers (only assigned ones)
chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False)

# Print formatted output
print("Chiral Centers (Atom Index, Configuration):")

for atom_index, config in chiral_centers:
    print(f"Atom {atom_index}: {config}")

print("\nTotal number of chiral centers:", len(chiral_centers))