import requests
import re
import xml.etree.ElementTree as ET
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors
from rdkit.Chem import rdchem


def get_smiles_and_formula_from_name(compound_name, anion, cation):
    """
    Функция для получения SMILES-строки и молекулярной формулы соединения
    по его названию из базы данных PubChem.

    Аргументы:
    compound_name (str): Название химического соединения.

    Возвращает:
    tuple: Кортеж, содержащий SMILES-строку и молекулярную формулу соединения.
           Если соединение не найдено, возвращает (None, None).
    """
    # URL для поиска соединения по названию
    search_url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}/property/CanonicalSMILES,MolecularFormula/JSON'.format(
        compound_name)

    # Отправка запроса к PubChem для поиска соединения
    search_response = requests.get(search_url)

    if search_response.status_code == 200:
        search_response.raise_for_status()

        # Получение списка идентификаторов соединений из ответа
        search_data = search_response.json()
        compound_ids = [ident['CID'] for ident in search_data['PropertyTable']['Properties']]

        if compound_ids:
            # Получение SMILES-строки и молекулярной формулы для первого найденного соединения
            cid = compound_ids[0]
            smiles_formula_url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/CanonicalSMILES,MolecularFormula/JSON'
            smiles_formula_response = requests.get(smiles_formula_url)
            smiles_formula_response.raise_for_status()

            smiles_formula_data = smiles_formula_response.json()
            smiles = smiles_formula_data['PropertyTable']['Properties'][0]['CanonicalSMILES']
            formula = smiles_formula_data['PropertyTable']['Properties'][0]['MolecularFormula']
            return smiles, formula
        else:
            print(f"Соединение '{compound_name}' не найдено в базе данных PubChem.")
            return None, None

    else:
        search_url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}/property/CanonicalSMILES,MolecularFormula/JSON'.format(
            anion)
        search_response = requests.get(search_url)
        if search_response.status_code == 200:

            if cation.find("mim"):

                tempreq = "1-butyl-3-methylimidazolium" + " " + anion
                search_url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}/property/CanonicalSMILES,MolecularFormula/JSON'.format(
                    tempreq)
                search_response = requests.get(search_url)
                search_response.raise_for_status()
                # Получение списка идентификаторов соединений из ответа
                search_data = search_response.json()
                compound_ids = [ident['CID'] for ident in search_data['PropertyTable']['Properties']]

                if compound_ids:
                    # Получение SMILES-строки и молекулярной формулы для первого найденного соединения
                    cid = compound_ids[0]
                    smiles_formula_url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/CanonicalSMILES,MolecularFormula/JSON'
                    smiles_formula_response = requests.get(smiles_formula_url)
                    smiles_formula_response.raise_for_status()

                    smiles_formula_data = smiles_formula_response.json()
                    smiles = smiles_formula_data['PropertyTable']['Properties'][0]['CanonicalSMILES']

                    # ТУТ НУЖНО СПЛИТНУТЬ ТК НЕ РАБОТАЕТ (".")
                    number_str = re.search(r'\d+', cation).group()
                    numberalkyl = int(number_str)
                    print(numberalkyl)

                    pattern = r'^CCCC'
                    smiles = re.sub(pattern, 'C' * numberalkyl, smiles)

                    mol = Chem.MolFromSmiles(smiles)

                    formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
                    #

                    #formula = smiles_formula_data['PropertyTable']['Properties'][0]['MolecularFormula']
                    return smiles, formula
                else:
                    print(f"Соединение '{compound_name}' не найдено в базе данных PubChem.")
                    return None, None

            if True:
                search_response.raise_for_status()
                # Получение списка идентификаторов соединений из ответа
                search_data = search_response.json()
                compound_ids = [ident['CID'] for ident in search_data['PropertyTable']['Properties']]

                if compound_ids:
                    # Получение SMILES-строки и молекулярной формулы для первого найденного соединения
                    cid = compound_ids[0]
                    smiles_formula_url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/CanonicalSMILES,MolecularFormula/JSON'
                    smiles_formula_response = requests.get(smiles_formula_url)
                    smiles_formula_response.raise_for_status()

                    smiles_formula_data = smiles_formula_response.json()
                    smiles = smiles_formula_data['PropertyTable']['Properties'][0]['CanonicalSMILES']
                    formula = smiles_formula_data['PropertyTable']['Properties'][0]['MolecularFormula']
                    return smiles, formula
                else:
                    print(f"Соединение '{compound_name}' не найдено в базе данных PubChem.")
                    return None, None

        else:
            print("Anion doesn't exist")
            return None, None









def molar_mass(compound: str, decimal_places=None) -> float:
    MM_of_Elements = {'H': 1.00794, 'He': 4.002602, 'Li': 6.941, 'Be': 9.012182, 'B': 10.811, 'C': 12.0107,
                      'N': 14.0067,
                      'O': 15.9994, 'F': 18.9984032, 'Ne': 20.1797, 'Na': 22.98976928, 'Mg': 24.305, 'Al': 26.9815386,
                      'Si': 28.0855, 'P': 30.973762, 'S': 32.065, 'Cl': 35.453, 'Ar': 39.948, 'K': 39.0983,
                      'Ca': 40.078,
                      'Sc': 44.955912, 'Ti': 47.867, 'V': 50.9415, 'Cr': 51.9961, 'Mn': 54.938045,
                      'Fe': 55.845, 'Co': 58.933195, 'Ni': 58.6934, 'Cu': 63.546, 'Zn': 65.409, 'Ga': 69.723,
                      'Ge': 72.64,
                      'As': 74.9216, 'Se': 78.96, 'Br': 79.904, 'Kr': 83.798, 'Rb': 85.4678, 'Sr': 87.62, 'Y': 88.90585,
                      'Zr': 91.224, 'Nb': 92.90638, 'Mo': 95.94, 'Tc': 98.9063, 'Ru': 101.07, 'Rh': 102.9055,
                      'Pd': 106.42,
                      'Ag': 107.8682, 'Cd': 112.411, 'In': 114.818, 'Sn': 118.71, 'Sb': 121.760, 'Te': 127.6,
                      'I': 126.90447, 'Xe': 131.293, 'Cs': 132.9054519, 'Ba': 137.327, 'La': 138.90547, 'Ce': 140.116,
                      'Pr': 140.90465, 'Nd': 144.242, 'Pm': 146.9151, 'Sm': 150.36, 'Eu': 151.964, 'Gd': 157.25,
                      'Tb': 158.92535, 'Dy': 162.5, 'Ho': 164.93032, 'Er': 167.259, 'Tm': 168.93421, 'Yb': 173.04,
                      'Lu': 174.967, 'Hf': 178.49, 'Ta': 180.9479, 'W': 183.84, 'Re': 186.207, 'Os': 190.23,
                      'Ir': 192.217,
                      'Pt': 195.084, 'Au': 196.966569, 'Hg': 200.59, 'Tl': 204.3833, 'Pb': 207.2, 'Bi': 208.9804,
                      'Po': 208.9824, 'At': 209.9871, 'Rn': 222.0176, 'Fr': 223.0197, 'Ra': 226.0254, 'Ac': 227.0278,
                      'Th': 232.03806, 'Pa': 231.03588, 'U': 238.02891, 'Np': 237.0482, 'Pu': 244.0642, 'Am': 243.0614,
                      'Cm': 247.0703, 'Bk': 247.0703, 'Cf': 251.0796, 'Es': 252.0829, 'Fm': 257.0951, 'Md': 258.0951,
                      'No': 259.1009, 'Lr': 262, 'Rf': 267, 'Db': 268, 'Sg': 271, 'Bh': 270, 'Hs': 269, 'Mt': 278,
                      'Ds': 281, 'Rg': 281, 'Cn': 285, 'Nh': 284, 'Fl': 289, 'Mc': 289, 'Lv': 292, 'Ts': 294, 'Og': 294,
                      '': 0}
    is_polyatomic = end = multiply = False
    polyatomic_mass, m_m, multiplier = 0, 0, 1
    element = ''

    for e in compound:
        if is_polyatomic:
            if end:
                is_polyatomic = False
                m_m += int(e) * polyatomic_mass if e.isdigit() else polyatomic_mass + MM_of_Elements[e]
            elif e.isdigit():
                multiplier = int(str(multiplier) + e) if multiply else int(e)
                multiply = True
            elif e.islower():
                element += e
            elif e.isupper():
                polyatomic_mass += multiplier * MM_of_Elements[element]
                element, multiplier, multiply = e, 1, False
            elif e == ')':
                polyatomic_mass += multiplier * MM_of_Elements[element]
                element, multiplier = '', 1
                end, multiply = True, False
        elif e == '(':
            m_m += multiplier * MM_of_Elements[element]
            element, multiplier = '', 1
            is_polyatomic, multiply = True, False
        elif e.isdigit():
            multiplier = int(str(multiplier) + e) if multiply else int(e)
            multiply = True
        elif e.islower():
            element += e
        elif e.isupper():
            m_m += multiplier * MM_of_Elements[element]
            element, multiplier, multiply = e, 1, False
    m_m += multiplier * MM_of_Elements[element]
    if decimal_places is not None:
        return round(m_m, decimal_places)
    return m_m


def count_bonds_mc(smiles):
    """
    Calculates the number of bonds in a molecule based on the McGowan method.

    Args:
        smiles (str): SMILES string representing the molecule.

    Returns:
        int: The number of bonds in the molecule.
    """

    # Создаем молекулярный объект из SMILES-строки
    mol = Chem.MolFromSmiles(smiles)

    # Получаем количество атомов
    num_atoms = mol.GetNumAtoms()

    # Получаем количество циклов
    num_rings = mol.GetRingInfo().NumRings()

    # Рассчитываем количество связей по формуле B = N - 1 + Rg
    num_bonds = num_atoms - 1 + num_rings

    return num_bonds


def allcount_bonds(smiles):
    # Создаем молекулярный объект из SMILES-строки
    mol = Chem.MolFromSmiles(smiles)

    # Добавляем водородные атомы
    mol = Chem.AddHs(mol)

    # Инициализируем счетчик связей
    bond_count = 0

    # Проходимся по всем связям в молекуле
    for bond in mol.GetBonds():
        # Получаем тип связи
        bond_type = bond.GetBondType()

        # Учитываем ароматические связи
        if bond_type == Chem.BondType.AROMATIC:
            bond_count += 1.5
        # Учитываем одинарные связи
        elif bond_type == Chem.BondType.SINGLE:
            bond_count += 1
        # Учитываем двойные связи
        elif bond_type == Chem.BondType.DOUBLE:
            bond_count += 2
        # Учитываем тройные связи
        elif bond_type == Chem.BondType.TRIPLE:
            bond_count += 3
        # Учитываем другие типы связей
        else:
            bond_order = rdchem.BondDict.lookup(bond.GetBondType())
            bond_count += bond_order

    return int(bond_count)


def calculate_mcgowan_volume(smiles, molar_mass, num_bonds):
    MCGOWAN_VOLUMES = {
        'H': 8.71, 'He': 6.75, 'Li': 22.23, 'Be': 20.27, 'B': 18.31, 'C': 16.35,
        'N': 14.39, 'O': 12.43, 'F': 10.47, 'Ne': 8.51, 'Na': 32.71, 'Mg': 30.75,
        'Al': 28.79, 'Si': 26.83, 'P': 24.87, 'S': 22.91, 'Cl': 20.95, 'Ar': 18.99,
        'K': 51.89, 'Ca': 50.28, 'Sc': 48.68, 'Ti': 47.07, 'V': 45.47, 'Cr': 43.86,
        'Mn': 42.26, 'Fe': 40.65, 'Co': 39.05, 'Ni': 37.44, 'Cu': 35.84, 'Zn': 34.23,
        'Ga': 32.63, 'Ge': 31.02, 'As': 29.42, 'Se': 27.81, 'Br': 26.21, 'Kr': 24.60,
        'Rb': 60.22, 'Sr': 58.61, 'Y': 57.01, 'Zr': 55.40, 'Nb': 53.80, 'Mo': 52.19,
        'Tc': 50.59, 'Ru': 48.98, 'Rh': 47.38, 'Pd': 45.77, 'Ag': 44.17, 'Cd': 42.56,
        'In': 40.96, 'Sn': 39.35, 'Sb': 37.75, 'Te': 36.14, 'I': 34.54, 'Xe': 32.93,
        'Cs': 77.25, 'Ba': 76.00, 'La': 74.75, 'Hf': 55.97, 'Ta': 54.71, 'W': 53.46,
        'Re': 52.21, 'Os': 50.96, 'Ir': 49.71, 'Pt': 48.45, 'Au': 47.20, 'Hg': 45.95,
        'Tl': 44.70, 'Pb': 43.45, 'Bi': 42.19, 'Po': 40.94, 'At': 39.69, 'Rn': 38.44,
        'Fr': 75.59, 'Ra': 74.34, 'Ac': 73.09, 'Ce': 73.49, 'Pr': 72.24, 'Nd': 70.99,
        'Pm': 69.74, 'Sm': 68.49, 'Eu': 67.23, 'Gd': 65.98, 'Tb': 64.73, 'Dy': 63.48,
        'Ho': 62.23, 'Er': 60.97, 'Tm': 59.72, 'Yb': 58.47, 'Lu': 57.22, 'Th': 71.83,
        'Pa': 70.58, 'U': 69.33, 'Np': 68.08, 'Pu': 66.83, 'Am': 65.57, 'Cm': 64.32,
        'Bk': 63.07, 'Cf': 61.82, 'Es': 60.57, 'Fm': 59.31, 'Md': 58.06, 'No': 56.81,
        'Lr': 55.56
    }

    mol = Chem.MolFromSmiles(smiles)

    # Подсчитываем общий объем МакГована
    total_volume = 0
    for atom in mol.GetAtoms():
        element = atom.GetSymbol()
        total_volume += MCGOWAN_VOLUMES[element]

    # Вычитаем 6.56 за каждую связь
    # num_bonds = mol.GetNumBonds()

    total_volume -= 6.56 * num_bonds

    # Добавляем вклад ионов (если есть)
    # total_charge = calculate_total_charge(smiles, molar_mass)
    # if total_charge:
    # total_volume += calculate_ion_volume(total_charge)

    return total_volume


compound_name = "C32mim acetate"

anion_name = compound_name[compound_name.find(" ") + 1:]
cation_name = compound_name[:compound_name.find(" ")]

smiles, formula = get_smiles_and_formula_from_name(compound_name, anion_name, cation_name)

if smiles and formula:
    molarweight_compound = molar_mass(formula)

    num_bonds_mc = count_bonds_mc(smiles)
    mcvolume = calculate_mcgowan_volume(smiles, molarweight_compound, num_bonds_mc)

    print(f"SMILES-строка для {compound_name}:   {smiles}")
    print(f"Молекулярная формула для {compound_name}:   {formula}")
    print(f"Молярная масса для {compound_name}:   {molarweight_compound}")
    print(f"Обьем Макгована {compound_name}:   {mcvolume}")

else:
    print(f"{compound_name} не найдена в базе данных.")

mol = Chem.MolFromSmiles(smiles)
img = Draw.MolToImage(mol)
img.show()

print(f"All bonds : {allcount_bonds(smiles)}")
print(f"All bonds mc : {count_bonds_mc(smiles)}")
