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

#  Maltose
maltose_smiles = "OC1C(O)C(O)C(CO)OC1OC1C(CO)OC(O)C(O)C1(O)"

# Tạo cấu trúc Maltose
maltose = Chem.MolFromSmiles(maltose_smiles)
maltose = Chem.AddHs(maltose)  # Thêm các nguyên tử hydro
AllChem.EmbedMolecule(maltose)
AllChem.UFFOptimizeMolecule(maltose)

print("Hiển thị Maltose:")
plotter_alpha = create_pyvista_molecule(maltose)
plotter_alpha.show()