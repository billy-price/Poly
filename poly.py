import itertools
from tabulate import tabulate

class Poly():
    '''
    Polynomials are represented as n-tuples
    where (a,b,c,d,e,f) represents a + bx + cx^2 + dx^3 + ex^4 + fx^6
    We also impose 
    '''


    def __init__(self, *coeffs, char=None, deg_mult=0, ):
        self.char = char
        self.coeffs = (0,)*deg_mult + coeffs
        self.normalise()

    def normalise(self):

        # apply mod char operation to scalars
        if self.char!=None:
            coeffs = list(x % self.char for x in self.coeffs)

        # dump excess zeroes
        while len(coeffs) and coeffs[-1]==0: coeffs.pop()

        self.coeffs = tuple(coeffs)
        
    def __iter__(self):
        return iter(self.coeffs)

    def __getitem__(self, key):
        return self.coeffs.__getitem__(key)

    def __add__(self, other):
        assert(self.char == other.char)
        return Poly(*(x+y for x,y in itertools.zip_longest(self, other, fillvalue=0)), char=self.char)
    
    def __sub__(self, other):
        assert(self.char == other.char)
        return Poly(*(x-y for x,y in itertools.zip_longest(self, other, fillvalue=0)), char=self.char)

    def __mul__(self, other):
        assert(self.char == other.char)
        result = [0] * (len(self) + len(other))
        for i, x in enumerate(self):
            for j, y in enumerate(other):
                result[i+j] += x*y
        
        return Poly(*result, char=self.char)

    def __rmul__(self, scalar):
        return Poly(*(coeff*scalar for coeff in self), char=self.char)

    def degree(self):
        return len(self)-1

    def __len__(self):
        return len(self.coeffs)

    def __repr__(self):
        if (len(self) == 0):
            return "0"
        else:
            return " + ".join((f"{coeff}" if (coeff!=1 or i==0) else "") +  # coefficient part
                              ((f"x^{i}" if i!=1 else "x") if i!=0 else "") # x part
                                    for i,coeff in enumerate(self) if coeff!=0)
    
    def __eq__(self, other):
        return self.coeffs == other.coeffs

    def get_rule_reduced(self, deg_to_replace, equal_poly):
        result = Poly(*self.coeffs, char=self.char)
        while (deg_to_replace <= result.degree()):
            # e.g if deg_to_replace=3 and equal_poly=x+1
            # this would replace the highest term a*x^n in self with a*x^(n-3)*(x+1), and then tack on the rest of it
            result = Poly(self[-1], deg_mult=(result.degree()-deg_to_replace), char=self.char) * equal_poly \
                        + Poly(*result.coeffs[:-1], char=self.char)
        
        return result

    @classmethod
    def iter_polynomials(cls, char, max_deg):
        
        for choice in itertools.product(range(char), repeat=max_deg+1):
            yield cls(*(reversed(choice)), char=char)

    @classmethod
    def plus_table(cls, char, max_deg):
        table = [[f]+[f + g for g in Poly.iter_polynomials(char, max_deg)] for f in Poly.iter_polynomials(char, 2)]

        table_headers = itertools.chain(("+",), Poly.iter_polynomials(char, max_deg))
        
        return tabulate(table, headers=table_headers, tablefmt='grid', stralign='center')
    
    @classmethod
    def times_table(cls, char, max_deg, deg_to_replace, equal_poly):
        table = [[f]+[(f * g).get_rule_reduced(deg_to_replace, equal_poly) for g in Poly.iter_polynomials(char, max_deg)] for f in Poly.iter_polynomials(char, max_deg)]

        table_headers = itertools.chain(("*",), Poly.iter_polynomials(char, max_deg))
        
        return tabulate(table, headers=table_headers, tablefmt='grid')

class FieldTimesTable():
    
    def __init__(self, char, max_deg, deg_to_replace, equal_poly):
        self.char = char
        self.max_deg = max_deg
        self.degree_to_replace = deg_to_replace
        self.equal_poly = equal_poly
        self.table = [[(f * g).get_rule_reduced(deg_to_replace, equal_poly) for g in Poly.iter_polynomials(char, max_deg)] for f in Poly.iter_polynomials(char, max_deg)]

    def __repr__(self):
        print(f"Times Table. characteristic={self.char}. max degree={self.max_deg}. rule : x^{self.degree_to_replace} = {self.equal_poly} ")
        return tabulate(self.table, tablefmt='grid')


class FieldPlusTable():

    def __init__(self, char, max_deg):
        self.char = char
        self.max_deg = max_deg
        self.table = [[f + g for g in Poly.iter_polynomials(char, max_deg)] for f in Poly.iter_polynomials(char, max_deg)]

    def __repr__(self):
        print(f"Plus Table. characteristic={self.char}. max degree={self.max_deg}.")
        return tabulate(self.table, tablefmt='grid')

if __name__ == "__main__":

    # This represents the rule x^3 = 1 + x
    equal_poly = Poly(1,1,char=2)
    deg_to_replace = 3
    
    CHAR = 2

    tt = FieldTimesTable(char=CHAR, max_deg=2, deg_to_replace=deg_to_replace, equal_poly=equal_poly)
    print(tt)
    pt = FieldPlusTable(char=CHAR, max_deg=2)
    print(pt)

    

    


    
