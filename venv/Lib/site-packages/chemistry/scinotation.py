# -*- coding: utf-8 -*-
"""
Created on Sun Apr 12 14:06:21 2020

@author: John Hallat
"""
from math import floor
from re import match, search
from fractions import Fraction
from decimal import Decimal

class Count:
    """Class for numbers that represent an exact count. These values do not affect
    the count of significant digits in calculations with Scinot (scientific notation) values"""

    def __init__(self, value):
        self.value = Fraction(value)

    def __repr__(self):
        return str(self.value)

    def __eq__(self, other):
        if isinstance(other, Count):
            return self.value == other.value
        else:
            return self.value == Fraction(other)

    def __gt__(self, other):
        if isinstance(other, Count):
            return self.value > other.value
        else:
            return self.value > Fraction(other)

    def __add__(self, other):
        if isinstance(other, Count):
            return Count(self.value + other.value)
        else:
            return Count(self.value + Fraction(other))

    def __mul__(self, other):
        if isinstance(other, Count):
            return Count(self.value * other.value)
        else:
            return Count(self.value * Fraction(other))

    def __int__(self):
        return floor(self.value)

## TODO Should be immutable
class Scinot:
    """Class for representing numbers in scientific notation"""
    
    def __init__(self, value):
        if isinstance(value, Scinot):
            self._integral = value._integral
            self._decimal = value._decimal
            self._exponent = value._exponent
        else:
            self._integral, self._decimal, self._exponent = self.parse(value)

    def parse(self, value):

        """Parses a string in the form of scientific notation"""
        _value = str(value)
        result = match('-{0,}[0-9]((\.[0-9]{1,}){0,1}x10\^|E-{0,1}[0-9]{1,}){0,1}', _value)
        if not result:
            message = f"Invalid format for scientific notation '{_value}'"
            raise Exception(message)

        integral_match = search('-{0,1}[0-9]{1,}', _value)
        integral = integral_match[0]

        decimal_match = search('(?<=\.)[0-9]{1,}', _value)
        decimal = decimal_match[0] if decimal_match else ''
    
        exponent_match = search('(?<=[\^,E])-{0,1}[0-9]{1,}', _value)
        exponent = exponent_match[0] if exponent_match else '0'

        while len(integral) > 1:
            trailer = integral[-1]
            integral = integral[:-1]
            decimal = trailer + decimal
            exponent = str(int(exponent) + 1)

        while integral == '0' and len(decimal) > 0:
            integral = decimal[0]
            decimal = decimal[1:]
            exponent = str(int(exponent) - 1)

        return (integral, decimal, exponent)  

    
    def sig_digits(self) -> int:
        "Determine the amount of significant digits from scientific notation"
        if self._decimal:
            return len(self._integral) + len(self._decimal)
        else:
            return 1


    def decimal(self) -> Decimal:
        """Returns a decimal representation of the scientific notation"""

        if self._decimal:
            _decver = Decimal(self._integral + '.' + self._decimal)
        else:
            _decver = Decimal(self._integral)
        return _decver * Decimal(10 ** int(self._exponent))

    # TODO Try to simplify
    def __round__(self, digits):
        _svalue = self._integral + self._decimal
        _exponent = self._exponent
        _svalue = _svalue[0:digits + 1]
        _len = len(_svalue)
        if _len < digits:
            _svalue += '0' * (digits - _len)
        _svalue = _svalue[0] + '.' + _svalue[1:] if len(_svalue) > 1 else _svalue[0]
        if _len <= digits:
            return Scinot(_svalue + 'x10^' + _exponent)
        return self._round_last(Scinot(_svalue + 'x10^' + _exponent))

    # TODO Try to simplify, may not be needed
    def _round_last(self, value):
        _value = Scinot(value)
        _number = _value._integral + _value._decimal
        _exponent = int(_value._exponent)
        _last = _number[-1]
        _number = _number[:-1]

        _trail = ''

        if int(_last) >= 5:
            if _number[-1] == '9':
                while _number[-1] == '9' and len(_number) > 1:
                    _trail = '0' + _trail
                    _number = _number[:-1]
            _next_last = int(_number[-1]) + 1
            if _next_last == 10:
                _trail = _trail[:-1] if len(_trail) > 0 else ''
                _exponent += 1
            _last = _number[-1]
            _number = _number[:-1] + str(_next_last) + _trail

        return Scinot(_number[0] + '.' + _number[1:] + 'x10^' + str(_exponent))


    def __int__(self):
        dec_value = self.decimal()
        return int(dec_value)

    def __repr__(self):
        if self._decimal:
            return self._integral + '.' + self._decimal + 'x10^' + self._exponent
        return self._integral + 'x10^' + self._exponent

    def __mul__(self, other):
        multiplicand = self.decimal()
        if isinstance(other, Scinot):
            multiplier = other.decimal()
            digits = min(self.sig_digits(), other.sig_digits())
        elif isinstance(other, Count):
            value = other.value.limit_denominator()
            multiplier = Decimal(value.numerator / value.denominator)
            digits = self.sig_digits()
        else:
            multiplier = Decimal(other)
            digits = self.sig_digits()

        product = multiplicand * multiplier

        return round(Scinot(product), digits)

    def __truediv__(self, other):
        dividend = self.decimal()
        if isinstance(other, Scinot):
            divisor = other.decimal()
            digits = min(self.sig_digits(), other.sig_digits())
        else:
            divisor = Decimal(other)
            digits = self.sig_digits()

        quotient = dividend / divisor
        rounded = round(Scinot(quotient), digits)

        return rounded

    def __add__(self, other):
        total = Scinot(self.decimal() + other.decimal())
        is_self_zero = self == 0
        is_other_zero = other == 0
        non_zero = [value.sig_digits() for value in [self, other] if value != 0]
        digits = min(non_zero)
        return round(total, digits)

    def __sub__(self, other):
        total = Scinot(self.decimal() - other.decimal())
        digits = min(self.sig_digits(), other.sig_digits())
        return round(total, digits)

    def __eq__(self, other):
        if isinstance(other, int):
            return self._eq_integer(other)
        if not isinstance(other, Scinot):
            return False
        return self._integral == other._integral and self._decimal == other._decimal and self._exponent == other._exponent

    def __lt__(self, other):
        if not isinstance(other, Scinot):
            return self.decimal() < other
        return self.decimal() < other.decimal()

    def _eq_integer(self, value: int) -> bool:
        """Return true if notation equals 1"""

        if not self._integral == str(value) or not self._exponent == '0':
            return False;

        if len(self._decimal) == 0 or len(self._decimal) == self._decimal.count('0'):
            return True
        else:
            return False



