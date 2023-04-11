from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdDepictor
from IPython.display import SVG

rdDepictor.SetPreferCoordGen(True)

def _get_avg_bond_length(mol):
    lengths = []
    conf = mol.GetConformer()
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        p1 = conf.GetAtomPosition(i)
        p2 = conf.GetAtomPosition(j)

        length = ((p1.x - p2.x) ** 2 + (p1.y - p2.y) ** 2) ** 0.5
        lengths.append(length)

    return sum(lengths) / len(lengths)

def draw_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol.GetNumConformers() == 0:
        rdDepictor.Compute2DCoords(mol, nSample=100, nFlipsPerSample=100)

    draw = rdMolDraw2D.MolDraw2DSVG(-1, -1)
    draw.drawOptions().padding = 0.1
    draw.drawOptions().atomLabelDeuteriumTritium = True
    draw.drawOptions().additionalAtomLabelPadding = 0.1
    draw.drawOptions().scaleBondWidth = False
    draw.drawOptions().fixedFontSize = 16
    draw.drawOptions().fixedBondLength = 42
    draw.drawOptions().drawMolsSameScale = False
    draw.drawOptions().scalingFactor = 20
    draw.DrawMolecule(mol)
    draw.FinishDrawing()
    svg = SVG(draw.GetDrawingText())

    return svg
