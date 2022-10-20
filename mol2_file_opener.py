import os
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
#from rdkit.Chem.Draw import IpythonConsole
from rdkit.Chem import Descriptors

def Mol2MolSupplier (file=None, sanitize=False):
    mols=[]
    with open(file, 'r') as mol2_file:
        mol2_lines = mol2_file.readline()
        while not mol2_file.tell() == os.fstat(mol2_file.fileno()).st_size:
            if mol2_lines.startswith("@<TRIPOS>MOLECULE"):
                mol = []
                mol.append(mol2_lines)
                mol2_lines = mol2_file.readline()
                while not mol2_lines.startswith("@<TRIPOS>MOLECULE"):
                    mol.append(mol2_lines)
                    mol2_lines = mol2_file.readline()
                    if mol2_file.tell() == os.fstat(mol2_file.fileno()).st_size:
                        mol.append(mol2_lines)
                        break
                mol[-1] = mol[-1].rstrip()
                block = ",".join(mol).replace(',', '')
                m = Chem.MolFromMol2Block(block, sanitize=sanitize)
            mols.append(m)
    return(mols)

from rdkit.Chem.MolStandardize import rdMolStandardize
salt_remover = rdMolStandardize.FragmentRemover()

if __name__ == '__main__':
    mols = Mol2MolSupplier(file='aligned_hets.mol2')
    for mol in mols:

        smiles = Chem.MolToSmiles(mol, kekuleSmiles=False)
        smiles_frame = Chem.MolToSmiles(Chem.RemoveHs(Chem.MolFromSmiles(smiles, sanitize=False), sanitize=False))
        print(smiles_frame)
        new_mol = Chem.MolFromSmiles(smiles_frame.replace("+", ""))
        AllChem.ComputeGasteigerCharges(new_mol)
        print([new_mol.GetAtomWithIdx(atom.GetIdx()).GetDoubleProp('_GasteigerCharge') for atom in new_mol.GetAtoms()])









    #m = Chem.MolFromSmiles('c1ccccc1C(=O)O')
#AllChem.ComputeGasteigerCharges(m)
#m.GetAtomWithIdx(0).GetDoubleProp('_GasteigerCharge')
#print(m)