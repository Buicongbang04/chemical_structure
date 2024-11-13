from rdkit import Chem
from rdkit.Chem import AllChem
import pyvista as pv
import numpy as np

def create_pyvista_molecule(mol):
    # Tạo plotter (thằng này dùng để tạo ra không gian 3D chứa các nguyên tử và liên kết)
    plotter = pv.Plotter()
    
    # Lấy tọa độ của các nguyên tử từ RDKit
    conf = mol.GetConformer()
    atom_positions = np.array([conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())])
    
    # Thêm các nguyên tử vào mô hình
    for i, pos in enumerate(atom_positions):
        atom = mol.GetAtomWithIdx(i) # Lấy nguyên tử từ RDKit
        if atom.GetSymbol() == "O":
            color = "red"  # Màu đỏ cho oxy
            radius = 0.35
        elif atom.GetSymbol() == "C":
            color = "black"  # Màu đen cho carbon
            radius = 0.3
        elif atom.GetSymbol() == "H":
            color = "white"  # Màu trắng cho hydro
            radius = 0.2
        else:
            color = "gray"  # Màu xám cho các nguyên tử khác
            radius = 0.25
        sphere = pv.Sphere(radius=radius, center=pos) # tạo nguyên tử
        plotter.add_mesh(sphere, color=color)
    
    # Thêm các liên kết giữa các nguyên tử
    for bond in mol.GetBonds():
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        start = atom_positions[begin_idx]
        end = atom_positions[end_idx]
        cylinder = pv.Cylinder(center=(start + end) / 2, direction=(end - start),
                               radius=0.1, height=np.linalg.norm(end - start))
        plotter.add_mesh(cylinder, color="gray")  # Màu xám cho liên kết
    
    return plotter

# α-D-Glucose mạch vòng, b-D-Glucose mạch vòng và D-Glucose mạch hở
alpha_glucose_smiles = "OC1C(CO)OC(O)C(O)C1O"
beta_glucose_smiles = "OC1C(CO)OC(CO)C1O"
open_chain_glucose_smiles = "O=C(C(O)C(O)C(O)C(O)C(O))"

# Tạo phần tử glucose alpha
alpha_glucose = Chem.MolFromSmiles(alpha_glucose_smiles)
alpha_glucose = Chem.AddHs(alpha_glucose)  # Thêm các nguyên tử hydro
AllChem.EmbedMolecule(alpha_glucose)
AllChem.UFFOptimizeMolecule(alpha_glucose)

# Tạo phần tử glucose beta
beta_glucose = Chem.MolFromSmiles(beta_glucose_smiles)
beta_glucose = Chem.AddHs(beta_glucose)  # Thêm các nguyên tử hydro
AllChem.EmbedMolecule(beta_glucose)
AllChem.UFFOptimizeMolecule(beta_glucose)

# Tạo phần tử glucose mạch hở
open_chain_glucose = Chem.MolFromSmiles(open_chain_glucose_smiles)
open_chain_glucose = Chem.AddHs(open_chain_glucose)  # Thêm các nguyên tử hydro
AllChem.EmbedMolecule(open_chain_glucose)
AllChem.UFFOptimizeMolecule(open_chain_glucose)

# Trực quan hóa bằng PyVista
print("Hiển thị α-D-Glucose (mạch vòng - alpha):")
plotter_alpha = create_pyvista_molecule(alpha_glucose)
plotter_alpha.show()

print("Hiển thị b-D-Glucose (mạch vòng - beta):")
plotter_open_chain = create_pyvista_molecule(beta_glucose)
plotter_open_chain.show()

print("Hiển thị D-Glucose (mạch hở):")
plotter_open_chain = create_pyvista_molecule(open_chain_glucose)
plotter_open_chain.show()