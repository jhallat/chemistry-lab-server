from chemistry.formula.parser import parse_formula
from chemistry.measurement import Measurement, ATOM, MOLES, validate_measurement, GRAMS, LITERS
from chemistry.periodictable import PeriodicTable
from chemistry.scinotation import Scinot

AVAGADROS_NUMBER = Scinot('6.022x10^23')

def mass(value: str) -> Measurement:
    _molar_mass = molar_mass(value)
    return _molar_mass / Measurement(AVAGADROS_NUMBER, ATOM/MOLES)

def molecules_in_mass(formula: str, grams: Measurement) -> Measurement:
    _molar_mass = molar_mass(formula)
    _grams = validate_measurement(grams, GRAMS)
    _mass = _grams * Measurement(AVAGADROS_NUMBER, ATOM/MOLES) / _molar_mass
    if _mass.value == 1:
        return Measurement(_mass.value, 'atom')
    else:
        return Measurement(_mass.value, 'atoms')

def molar_mass(formula: str) -> Measurement:
    """Calculates the molar mass of a compound"""
    if isinstance(formula, list):
        atoms = formula
    else:
        atoms = parse_formula(formula).flatten()
    total_mass = Measurement("0", GRAMS / MOLES)
    table = PeriodicTable()
    for index, atom in enumerate(atoms):
        count, symbol = atom
        mass = Scinot(table[symbol].atomic_mass)
        product = Measurement(mass * count, GRAMS / MOLES)
        total_mass += product
    return total_mass


def moles_from_grams(grams, formula):
    _mmass = molar_mass(formula)
    _grams = validate_measurement(grams, GRAMS)
    value = _grams / _mmass
    return Measurement(value, MOLES)


def molarity(moles, liters=None):
    _moles = validate_measurement(moles, MOLES)
    if liters == None:
        return Measurement(_moles.value, MOLES / LITERS)

    _liters = validate_measurement(liters, LITERS)
    molarity = _moles / _liters
    return molarity

def grams_in_solution(molarity: Measurement, molar_mass: Measurement, liters: Measurement) -> Measurement:
    """ Calculates the mass in grams of a compound in a solution  of a specified molarity
    """
    _molarity = validate_measurement(molarity, MOLES / LITERS)
    _molar_mass = validate_measurement(molar_mass, GRAMS / MOLES)
    _liters = validate_measurement(liters, LITERS)

    return _liters * _molar_mass * _molarity


def molarity_of_solution(grams: Measurement, molar_mass: Measurement, liters:Measurement) -> Measurement:
    """ Calculates the molarity of a solution
    """

    _grams = validate_measurement(grams, GRAMS)
    _molar_mass = validate_measurement(molar_mass, GRAMS / MOLES)
    _liters = validate_measurement(liters, LITERS)

    _moles = _grams / _molar_mass

    return _moles / _liters


def moles_in_solution(molarity, liters):
    """Finds the number of moles in a solution of a specified molarity
    """

    _molarity = validate_measurement(molarity, MOLES / LITERS)
    _liters = validate_measurement(liters, LITERS)

    return _molarity * _liters


def volume_from_solution(moles: Measurement, molarity: Measurement) -> Measurement:
    """Finds the volume of a solution that will contain the specified number of moles when taken
       from a solution of a specified molarity
    """

    _moles = validate_measurement(moles, MOLES)
    _molarity = validate_measurement(molarity, MOLES / LITERS)

    return _moles / _molarity



