import streamlit as st
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, Descriptors
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# Set page config
st.set_page_config(
    page_title="Artemether Chiral Analyzer",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 Artemether Chiral Center Analysis")
st.markdown("---")

# SMILES for Artemether (β-artemether)
ARTEMETHER_SMILES = "CC(C)C1CC2=CC(=O)C(=C3CCOC(C)(C)C3C2C1)C"

# Create molecule object
mol = Chem.MolFromSmiles(ARTEMETHER_SMILES)

if mol is None:
    st.error("Could not parse Artemether structure")
else:
    # Sidebar information
    with st.sidebar:
        st.header("📋 Artemether Info")
        st.write("""
        **Common Name:** Artemether
        
        **Type:** Semi-synthetic derivative of artemisinin
        
        **Source:** Artemisia annua (Sweet wormwood)
        
        **Use:** Antimalarial drug
        
        **Form:** β-artemether (active isomer)
        """)
    
    # Main layout
    col1, col2 = st.columns(2)
    
    # Column 1: 2D Structure
    with col1:
        st.subheader("2D Molecular Structure")
        
        # Highlight chiral centers
        chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        
        img = Draw.MolToImage(mol, size=(400, 400))
        st.image(img, use_column_width=True)
    
    # Column 2: Molecular Properties
    with col2:
        st.subheader("🔬 Molecular Properties")
        
        properties = {
            "Molecular Weight": f"{Descriptors.MolWt(mol):.2f} g/mol",
            "Molecular Formula": Chem.rdMolDescriptors.CalcMolFormula(mol),
            "LogP": f"{Descriptors.MolLogP(mol):.2f}",
            "H-Bond Donors": f"{Descriptors.NumHDonors(mol)}",
            "H-Bond Acceptors": f"{Descriptors.NumHAcceptors(mol)}",
            "Rotatable Bonds": f"{Descriptors.NumRotatableBonds(mol)}",
            "Aromatic Rings": f"{Descriptors.NumAromaticRings(mol)}"
        }
        
        for prop, value in properties.items():
            st.metric(prop, value)
    
    st.markdown("---")
    
    # Chiral Center Analysis
    st.subheader("⚛️ Chiral Center Analysis")
    
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    
    if chiral_centers:
        col1, col2, col3 = st.columns(3)
        
        with col1:
            st.metric("Number of Chiral Centers", len(chiral_centers))
        
        with col2:
            st.metric("Possible Stereoisomers", f"2^{len(chiral_centers)} = {2**len(chiral_centers)}")
        
        with col3:
            st.metric("Natural Isomer", "β-artemether")
        
        st.write("**Chiral Center Positions:**")
        chiral_data = {
            "Atom Index": [cc[0] for cc in chiral_centers],
            "Chirality": [cc[1] for cc in chiral_centers]
        }
        st.dataframe(pd.DataFrame(chiral_data), use_container_width=True)
    else:
        st.info("No chiral centers detected (achiral molecule)")
    
    st.markdown("---")
    
    # Stereoisomer Information
    st.subheader("🔄 Stereoisomer Information")
    
    st.write(f"""
    **Artemether has {len(chiral_centers)} stereogenic centers.**
    
    This means there are **{2**len(chiral_centers)} possible stereoisomers** (including enantiomers and diastereomers).
    
    The naturally active form is **β-artemether**, which is the specific stereoisomer used as an antimalarial drug.
    
    The stereochemistry is critical for:
    - ✅ Antimalarial activity
    - ✅ Pharmacokinetics (how the body processes it)
    - ✅ Safety and efficacy
    """)
    
    # 3D Structure Section
    st.subheader("🌐 3D Structure Visualization")
    
    # Generate 3D coordinates
    mol_3d = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol_3d, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol_3d)
    
    # Extract coordinates
    conf = mol_3d.GetConformer()
    positions = conf.GetPositions()
    
    # Create 3D scatter plot
    fig = go.Figure()
    
    # Add atoms
    for i, atom in enumerate(mol_3d.GetAtoms()):
        pos = positions[i]
        atomic_num = atom.GetAtomicNum()
        
        # Color by element
        if atomic_num == 6:  # Carbon
            color = 'gray'
        elif atomic_num == 1:  # Hydrogen
            color = 'white'
        elif atomic_num == 8:  # Oxygen
            color = 'red'
        elif atomic_num == 7:  # Nitrogen
            color = 'blue'
        else:
            color = 'green'
        
        fig.add_trace(go.Scatter3d(
            x=[pos[0]], y=[pos[1]], z=[pos[2]],
            mode='markers',
            marker=dict(size=8, color=color),
            text=f"{atom.GetSymbol()} (Index: {i})",
            hoverinfo='text'
        ))
    
    # Add bonds
    for bond in mol_3d.GetBonds():
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        begin_pos = positions[begin_idx]
        end_pos = positions[end_idx]
        
        fig.add_trace(go.Scatter3d(
            x=[begin_pos[0], end_pos[0]],
            y=[begin_pos[1], end_pos[1]],
            z=[begin_pos[2], end_pos[2]],
            mode='lines',
            line=dict(color='lightgray', width=2),
            hoverinfo='skip',
            showlegend=False
        ))
    
    fig.update_layout(
        title="3D Artemether Structure",
        scene=dict(
            xaxis_title="X",
            yaxis_title="Y",
            zaxis_title="Z"
        ),
        width=800,
        height=600,
        showlegend=False
    )
    
    st.plotly_chart(fig, use_container_width=True)
    
    st.markdown("---")
    st.info("💡 **Tip:** Rotate the 3D structure with your mouse to explore the molecule from different angles!")
