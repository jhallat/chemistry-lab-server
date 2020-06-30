from decimal import Decimal
from enum import Enum

from chemistry.formula.tokenizer import FormulaTokenType, FormulaToken, tokenize
from chemistry.periodictable import Ions, PeriodicTable, Ion
from chemistry.scinotation import Count


class FormulaParseError(Exception):
    pass

class FormulaNodeType(Enum):
    COMPOUND = 0
    ATOM = 1
    MONATOMIC_ION = 2
    POLYATOMIC_ION = 3

class ParseState(Enum):
    START = 1
    COEFFICIENT = 2
    SYMBOL = 3
    SUBSCRIPT = 4

class CompoundState(Enum):
    UNDEFINED = 0
    GAS = 1
    LIQUID = 2
    SOLID = 3
    AQUEOUS = 4

class FormulaRoot:

    def __init__(self, symbol, children):
        self.symbol = symbol
        self.children = children

    def __eq__(self, other):
        if isinstance(other, FormulaRoot):
            return self.symbol == other.symbol and self.children == other.children
        else:
            return False

    def __repr__(self):
        return f"root(symbol='{self.symbol},children={self.children}'"

    def __getitem__(self, item):
        return self.children[item]

    def flatten(self):
        atoms = []
        for compound in self.children:
            atoms += self._flatten_compound(compound.count, compound.children)

        atom_map = {}
        for atom in atoms:
            atom_sum = atom_map.get(atom[1], Count(0))
            atom_sum += atom[0]
            atom_map[atom[1]] = atom_sum

        return [(value, key) for key, value in atom_map.items()]

    def _flatten_compound(self, count, compound):

        atoms = []
        for node in compound:
            if node.type == FormulaNodeType.ATOM:
                atoms.append((node.count * count, node.symbol))
            else:
                atoms += self._flatten_compound(node.count, node.children)
        return atoms

#TODO need to convert to a class and implement __eq__
class FormulaNode:

    def __init__(self, count, symbol, type, children,):
        self.count = count
        self.symbol = symbol
        self.type = type
        self.children = children

    def __eq__(self, other):
        if isinstance(other, FormulaNode):
            return self.count == other.count \
                   and self.symbol == other.symbol \
                   and self.type == other.type \
                   and self.children == other.children
        else:
            return False

    def __repr__(self):
        return f"node(count={self.count},symbol={self.symbol},type={self.type},children={self.children}"

class CompoundNode(FormulaNode):

    def __init__(self, count, symbol, children, state=CompoundState.UNDEFINED):
        super().__init__(count, symbol, FormulaNodeType.COMPOUND, children)
        self.state = state

    def __eq__(self, other):
        if isinstance(other, CompoundNode):
            return super().__eq__(other) and self.state == other.state
        if isinstance(other, FormulaNode):
            return super().__eq__(other) and self.state == CompoundState.UNDEFINED
        return False

    def __repr__(self):
        return f"compound(count={self.count},symbol={self.symbol},state={self.state},children={self.children})"

class AtomNode(FormulaNode):

    def __init__(self, count, symbol):
        super().__init__(count, symbol, FormulaNodeType.ATOM, [])

    def __repr__(self):
        return f"atom(count={self.count},symbol={self.symbol})"

class PolyatomicIonNode(FormulaNode):

    def __init__(self, count, symbol, children):
        super().__init__(count, symbol, FormulaNodeType.POLYATOMIC_ION, [])

    def __repr__(self):
        return f"polyatomic_ion(count={self.count},symbol={self.symbol},children={self.children})"

def _reduce(elements: [(Decimal, str)]) -> [(Decimal, str)]:
    element_map = {}
    reduced_elements = []
    for element in elements:
        amount = element_map.get(element[1], Decimal(0.000))
        amount += element[0]
        element_map[element[1]] = amount

    for element in elements:
        if element[1] in element_map:
            reduced_elements.append((element_map.pop(element[1]), element[1]))

    return reduced_elements

def state_start(index, token, coefficient, symbol, polyatomic, compound_state):
    if token.type == FormulaTokenType.COEFFICIENT:
        coefficient = Decimal(token.value)
        state = ParseState.COEFFICIENT
    elif token.type == FormulaTokenType.SYMBOL:
        symbol = token.value
        state = ParseState.SYMBOL
    elif token.type == FormulaTokenType.SUBSCRIPT and polyatomic:
        coefficient = Decimal(token.value)
        state = ParseState.COEFFICIENT
    elif token.type == FormulaTokenType.STATE:
        if token.value == 's':
            compound_state = CompoundState.SOLID
        if token.value == 'l':
            compound_state = CompoundState.LIQUID
        if token.value == 'g':
            compound_state = CompoundState.GAS
        if token.value == 'aq':
            compound_state = CompoundState.AQUEOUS
        state = ParseState.START
    else:
        raise FormulaParseError(f"unexpected token '{token.value}' at index {index}, expected coefficient or symbol")
    return coefficient, symbol, state, compound_state

def state_coefficient(index, token):
    if token.type == FormulaTokenType.SYMBOL:
        symbol = token.value
        state = ParseState.SYMBOL
    else:
        raise FormulaParseError(f"unexpected token '{token.value}' at index {index}, expected symbol")
    return symbol, state

def state_symbol(index, token, coefficient, symbol, atoms, compound):
    if token.type == FormulaTokenType.SYMBOL:
        atoms.append(FormulaNode(Count(1), symbol, FormulaNodeType.ATOM, []))
        compound += symbol
        symbol = token.value
        state = ParseState.SYMBOL
    elif token.type == FormulaTokenType.SUBSCRIPT:
        atoms.append(FormulaNode(Count(token.value), symbol, FormulaNodeType.ATOM, []))
        compound += symbol + token.value
        symbol = ''
        state = ParseState.SUBSCRIPT
    else:
        raise FormulaParseError(f"unexpected token '{token.value}' at index {index}, expected symbol, subscript or '*'")
    return coefficient, symbol, atoms, compound, state

def state_subscript(index, token):
    if token.type == FormulaTokenType.SYMBOL:
        symbol = token.value
        state = ParseState.SYMBOL
    else:
        raise FormulaParseError(f"unexpected token '{token.value}' at index {index}, expected symbol or '*'")
    return symbol, state

def _parse_pass_two(tokens: [FormulaToken], polyatomic = False):

    atoms = []
    state = ParseState.START
    coefficient = Decimal('1.000')
    symbol = ''
    compound = ''
    compound_state = CompoundState.UNDEFINED

    for index, token in enumerate(tokens):
        if state == ParseState.START:
            if isinstance(token, list):
                atoms.append(_parse_pass_two(token, True))
            else:
                coefficient, symbol, state, compound_state = \
                    state_start(index, token, coefficient, symbol, polyatomic, compound_state)

        elif state == ParseState.COEFFICIENT:
            if isinstance(token, list):
                atoms.append(_parse_pass_two(token, True))
            else:
                symbol, state = state_coefficient(index, token)
        
        elif state == ParseState.SYMBOL:
            if isinstance(token, list):
                atoms.append(FormulaNode(Count(1), symbol, FormulaNodeType.ATOM, []))
                ion = _parse_pass_two(token, True)
                atoms.append(ion)
                compound += symbol
                symbol = ''
                if ion.count > 1:
                    compound += f'({ion.symbol}){str(ion.count)}'
                else:
                    compound += ion.symbol
            else:
                coefficient, symbol, atoms, compound, state =\
                    state_symbol(index, token, coefficient, symbol, atoms, compound)
        
        elif state == ParseState.SUBSCRIPT:
            if isinstance(token, list):
                ion = _parse_pass_two(token, True)
                atoms.append(ion)
                if ion.count > 1:
                    compound += f'({ion.symbol}){str(ion.count)}'
                else:
                    compound += ion.symbol
            else:
                symbol, state = state_subscript(index, token)
        else:
            raise FormulaParseError(f"unexpected state '{state}")    

    if len(symbol) > 0 and state == ParseState.SYMBOL:
        compound += symbol
        atoms.append(FormulaNode(Count(1), symbol, FormulaNodeType.ATOM, []))

    if polyatomic:
        return FormulaNode(Count(coefficient), compound, FormulaNodeType.POLYATOMIC_ION, atoms)
    else:
        return CompoundNode(Count(coefficient), compound, atoms, compound_state)

def _parse_pass_one(tokens):

    COMPOUND = 0
    POLYATOMIC = 1
    POLYATOMIC_END = 2

    root = []
    compound = []
    polyatomic = []
    state = COMPOUND
    for token in tokens:
        if token.type == FormulaTokenType.COMBINDED_WITH:
            root.append(compound)
            compound = []
        if token.type == FormulaTokenType.POLYATOMIC_START:
            state = POLYATOMIC
        if token.type == FormulaTokenType.POLYATOMIC_END:
            state = POLYATOMIC_END
        if token.type == FormulaTokenType.COEFFICIENT:
            if state == COMPOUND:
                compound.append(token)
            if state == POLYATOMIC:
                polyatomic.append(token)
            if state == POLYATOMIC_END:
                compound.append(polyatomic)
                compound.append(token)
                polyatomic = []
                state = COMPOUND
        if token.type == FormulaTokenType.SUBSCRIPT:
            if state == COMPOUND:
                compound.append(token)
            if state == POLYATOMIC:
                polyatomic.append(token)
            if state == POLYATOMIC_END:
                polyatomic.insert(0, token)
                compound.append(polyatomic)
                polyatomic = []
                state = COMPOUND
        if token.type == FormulaTokenType.SYMBOL:
            if state == COMPOUND:
                compound.append(token)
            if state == POLYATOMIC:
                polyatomic.append(token)
            if state == POLYATOMIC_END:
                compound.append(polyatomic)
                polyatomic = []
                compound.append(token)
                state = COMPOUND
        if token.type == FormulaTokenType.STATE:
            if state == COMPOUND:
                compound = [token] + compound
    if polyatomic:
        compound.append(polyatomic)
    root.append(compound)
    return root


#TODO add a precision parameter for multiply coefficient with subscript. Default = 3
def parse_formula(formula: str):

    tokens = tokenize(formula)
    print(f'tokens = {tokens}')
    root = _parse_pass_one(tokens)
    print(f'pass one = {root}')
    children = [_parse_pass_two(compound) for compound in root]
    return FormulaRoot(formula, children)


def _find_ion(symbol):
    prefixes = ['di', 'tri', 'tetra', 'tetr', 'penta', 'pent', 'hexa', 'hex', 'hepta', 'hept', 'octa', 'oct']
    ions = Ions()
    periodic_table = PeriodicTable()
    if symbol in ions:
        return ions[symbol]
    if symbol in periodic_table:
        atom = periodic_table[symbol]
        return Ion(atom.symbol, atom.name, 0)
    for prefix in prefixes:
        if symbol.startswith(prefix):
            return(_find_ion(symbol[len(prefix):]))

    raise FormulaParseError(f'unknown symbol: {symbol}')

def parse_ion_equation(equation: str):

    parts = equation.split('->')
    if '+' in parts[0]:
        reactants = [_find_ion(reactant.strip()) for reactant in parts[0].split('+')]
    else:
        reactants = [_find_ion(reactant.strip()) for reactant in parts[0].split()]
    return reactants


