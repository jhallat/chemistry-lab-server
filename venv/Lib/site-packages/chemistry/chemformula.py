from collections import namedtuple
from decimal import Decimal

from chemistry.formula.parser import parse_formula, parse_ion_equation, FormulaNodeType
from chemistry.measurement import Measurement, grams, validate_measurement, GRAMS, MOLES
from chemistry.mole import molar_mass, moles_from_grams
from chemistry.periodictable import PeriodicTable, Ions, Ion, CommonCompounds

Component = namedtuple("Component", "count symbol mass mass_percent")

class Composition:

    def __init__(self):
        self._components = []

    def __len__(self):
        return len(self._components)

    def __getitem__(self, item):
        if str(item).isnumeric():
            return self._components[item]
        else:
            return next(component for component in self._components if component.symbol == item)

    def __mul__(self, other):
        composition = Composition()
        for count, symbol, mass, mass_percent in self._components:
            mult_component = Component(count * other, symbol, mass * other, mass_percent)
            composition.append(mult_component)
        return composition

    def append(self, component):
        self._components.append(component)

    def __repr__(self):
        return str(self._components)

    def _formula(self):
        formula = ''
        for count, symbol, _, _ in self._components:
            formula += symbol
            formula += str(int(count)) if count != 1 else ''
        return formula

    def __format__(self, format_spec):
        if format_spec == 'f':
            return self._formula()
        else:
            #TODO return a more readable representation
            return str(self)

# def _recursive_composition(_composition, _molar_mass, mass:Measurement = grams('1.000') ) -> [Component]:
#     composition = Composition()
#     for element in _composition:
#         if isinstance(element, list):
#             composition.append(_recursive_composition(element, _molar_mass, mass))
#         else:
#             count, symbol = element
#             _mass_percent = molar_mass(symbol)/_molar_mass
#             _mass_percent = _mass_percent.value.decimal() * int(count)
#             _mass = mass * _mass_percent
#             _mass_percent = round(_mass_percent * 100, 2)
#             composition.append(Component(count, symbol, _mass, _mass_percent))
#
#     return composition

def composition(formula: str, mass:Measurement = grams('1.000') ) -> [Component]:

    composition = Composition()
    elements = parse_formula(formula).flatten()
    _molar_mass = molar_mass(elements)
    for element in elements:
       count, symbol = element
       _mass_percent = molar_mass(symbol)/_molar_mass
       _mass_percent = _mass_percent.value.decimal() * int(count)
       _mass = mass * _mass_percent
       _mass_percent = round(_mass_percent * 100, 2)
       composition.append(Component(count, symbol, _mass, _mass_percent))

    return composition


def formula_from_percent(elements: [(str, Decimal)], molar_mass = None) -> str:
    """Creates the simplest formula based on elements and percentages"""
    return formula_from_mass(elements)


def formula_from_mass(elements: [(str, Measurement)], molar_mass = None) -> str:
#    _masses = [validate_measurement(mass, GRAMS) for (_, mass) in elements]
#    _elements = [element for (element, _) in elements]

    _elements = [element for (element, _) in elements]
    _masses = [moles_from_grams(grams(percent), element) for (element, percent) in elements]

    min_mass = min(_masses)

    _subscripts = []
    for mass in _masses:
        subscript = mass.value.decimal() / min_mass.value.decimal()
        subscript = round(subscript, 4)
        _subscripts.append(subscript)

    formula = ''
    _subscripts = _simplify(_subscripts)
    for element, subscript in zip(_elements, _subscripts):
        formula += element + (str(subscript) if subscript > 1 else '')

    if molar_mass:
        return mol_formula_from_simple_formula(formula, molar_mass)
    return formula

def mol_formula_from_simple_formula(simple: str, ex_molar_mass: Measurement) -> str:
    _ex_molar_mass = validate_measurement(ex_molar_mass, GRAMS / MOLES)
    _sm_molar_mass = molar_mass(simple)
    multiplier = int(_ex_molar_mass.value / _sm_molar_mass.value)
    simple_composition = composition(simple)
    return f'{(simple_composition * multiplier):f}'

def _simplify(numbers: [Decimal]) -> Decimal:

    _numbers = numbers[0:]
    with_decimal = [number for number in _numbers if number - int(number) >= 0.1]
    prevent_infinite = 0
    while (len(with_decimal) > 0 and prevent_infinite < 100):
        prevent_infinite += 1
        decimal_part = with_decimal[0] - int(with_decimal[0])
        multiplier = round(1 / decimal_part, 0)
        _numbers = [number * multiplier for number in _numbers]
        with_decimal = [number for number in _numbers if number - int(number) >= 0.1]

    return [round(number,0) for number in _numbers]

def is_polyatomic(ion: str):
    return len([char for char in ion if char.isupper()]) > 1

def _binary_molecular_formula(products, ions):
    prefixes = [('2','di'),
                ('3','tri'),
                ('4','tetra'),
                ('4','tetr'),
                ('5','penta'),
                ('5','pent'),
                ('6','hexa'),
                ('6','hex'),
                ('7','hepta'),
                ('7','hept'),
                ('8','octa'),
                ('8','oct')]
    subscript_one = ''
    subscript_two = ''
    _products = products.split()
    for count, prefix in prefixes:
        if _products[0].startswith(prefix):
            subscript_one = count
        if _products[1].startswith(prefix):
            subscript_two = count
        if len(subscript_one + subscript_two) == 2:
            break
    return ions[0].symbol + subscript_one + ions[1].symbol + subscript_two


def predict_formula(products):
    """Predicts the formula of combining two ions.

    Ions can be added as a string:
    >>> predict_formula('Li + O')
    'Li2O'

    or as a list:
    >>> predict_formula(['Ni', 'S'])
    'NiS'

    if the ions are valid, but a formula cannot be predicted,
    the original ions are returned
    >>> predict_formula(['Zn','Ag'])
    'Zn + Ag'

    formulas can also be predicted based on the name of the
    compound
    >>> predict_formula('dinitrogen pentoxide')
    'N2O5'

    common names for some compounds may also be used
    >>> predict_formula('water')
    'H2O'

    """
    common_compounds = CommonCompounds()
    if isinstance(products, list):
        ions = parse_ion_equation(' + '.join(products))
    else:
        if products in common_compounds:
            return common_compounds[products]
        ions = parse_ion_equation(products)
    if len(ions) == 1:
        return ions[1].symbol
    if ions[0].charge > 0 and ions[1].charge > 0:
        return f'{ions[0].symbol} + {ions[1].symbol}'
    if ions[0].charge < 0 and ions[1].charge < 0:
        return _binary_molecular_formula(products, ions)

    positive, negative = (ions[0], ions[1]) if ions[0].charge > ions[1].charge else (ions[1], ions[0])
    product_charge = abs(positive.charge * negative.charge)
    positive_subscript = product_charge // positive.charge
    negative_subscript = abs(product_charge // negative.charge)
    if positive_subscript == negative_subscript:
        positive_subscript, negative_subscript = 1,1
    p_symbol = positive.symbol
    n_symbol = negative.symbol
    if is_polyatomic(p_symbol) and positive_subscript > 1:
        p_symbol = '(' + p_symbol + ')'
    if is_polyatomic(n_symbol) and negative_subscript > 1:
        n_symbol = '(' + n_symbol + ')'
    formula = p_symbol
    formula += str(positive_subscript) if positive_subscript > 1 else ''
    formula += n_symbol
    formula += str(negative_subscript) if negative_subscript > 1 else ''
    return formula

def compound_name(formula):
    periodic_table = PeriodicTable()
    ions = Ions()
    common_compounds = CommonCompounds()
    if formula in common_compounds:
        return common_compounds[formula]
    root = parse_formula(formula)
    components = root[0].children
    if len(components) == 2:
        cation = _find_ion(components[0], periodic_table, ions)
        anion = _find_ion(components[1], periodic_table, ions)
        first = cation.name
        if isinstance(cation.charge, list):
            total_charge = abs(int(components[1].count) * anion.charge)
            first += ('(' + ('I' * total_charge) + ')')
        second = anion.name
        if '-' in first:
            first = first.replace('-','')
        if '-' in second:
            second = second.split('-')[0] + 'ide'
        #TODO there is probably a better way to do this
        if anion.charge != -1 and cation.charge != 1:
            first = prefix(components[0].count, cation.symbol[0]) + first
            second = prefix(components[1].count, anion.symbol[0]) + second
        return  first + ' ' + second
    else:
        return 'unknown'

def prefix(count, first_letter):
    vowels = {'a', 'e', 'i', 'o', 'u'}
    _count = int(count)
    if _count < 2 or _count > 8:
        return ''
    prefix = ['di','tri','tetra','penta','hexa','hepta','octa'][_count - 2]
    if first_letter.lower() in vowels and prefix[-1] == 'a':
        prefix = prefix[:-1]
    return prefix

def _find_ion(component, ptable, ions):
    if component.type == FormulaNodeType.ATOM:
        if component.symbol in ions:
            return ions[component.symbol]
        elif component.symbol in ptable:
            atom = ptable[component.symbol]
            #TODO at some point periodic table should store charges, but 0
            #is ok at this point
            return Ion(atom.symbol, atom.name, 0)
    if component.type == FormulaNodeType.POLYATOMIC_ION:
        if component.symbol in ions:
            return ions[component.symbol]
    return Ion('','unknown',0)
