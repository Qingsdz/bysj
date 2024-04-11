"""
此脚本实现对类似于变换多项式进行化简，合并，乘积等计算


"""
from itertools import chain
default_cal = 7

neg_degree = -7

class CoffTerm:
    """
        exp
    word, der
    """
    def __init__(self, word, exponent, derivatives):
        self.exponent = exponent
        self.word = word
        self.derivatives = derivatives

    def __repr__(self):
        return f"Trem({self.word}, {self.exponent}, {self.derivatives})"
    
    def __eq__(self, other):
        if isinstance(other, CoffTerm):
            return self.word == other.word and self.exponent == other.exponent and self.derivatives == other.derivatives
        return False
    
    def __str__(self):
        deriva = ""
        if self.derivatives != 0:
            deriva = ","
            deriva += f"x({self.derivatives})"
        if self.word == "":
            return str("")
        else:
            if self.exponent == 0:
                return str("")
            elif self.exponent == 1:
                return f"{self.word}{deriva}"
            else:
                return f"{self.word}{deriva}^{self.exponent}"

class Term:
    """
    front_cofficient   cofficient    exponent
           a           [CT, CT, ..]   exp
          exp    
    a CT p
    """
    def __init__(self, front_cofficient, coefficient, exponent): #f_c float  coeff [c_t. c_t, ...]  exponent int
        self.front_coefficient = front_cofficient
        
        coeff = []
        for coeff_itr in coefficient:
            flag = 1
            for itr in coeff:
                if coeff_itr.word == itr.word and coeff_itr.derivatives == itr.derivatives:
                    itr.exponent += coeff_itr.exponent
                    flag = 0
            if flag == 1:
                coeff.append(coeff_itr)
        self.coefficient = sorted(coeff, key = lambda coeff:coeff.word) #[c_t. c_t, ...] 可能具有重复项，下面消除重复项为指数
        self.exponent = exponent

    def __repr__(self):
        return f"Term({self.front_coefficient}, {self.coefficient}, {self.exponent})"

    def __str__(self):
        if self.exponent == 0:
            words = [str(coeff) for coeff in self.coefficient]
            wordss = ""
            for w in words:
                wordss += (w+"*")
            return f"{self.front_coefficient}{wordss}"
        elif self.exponent == 1:
            words = [str(coeff) for coeff in self.coefficient]
            wordss = ""
            for w in words:
                wordss += (w+"*")
            return f"{self.front_coefficient}{wordss}p"
        else:
            words = [str(coeff) for coeff in self.coefficient]
            wordss = ""
            for w in words:
                wordss += (w+"*")
            return f"{self.front_coefficient}{wordss}p**{self.exponent}"
        
    def is_combine(self, other):
        if isinstance(other, Term):
            return self.coefficient == other.coefficient and self.exponent == other.exponent
        return False

class Polynomial:
    """
    T1 + T2 + T3 + ...
        terms    
    [T1, T2, T3, ...]
    """
    
    def __init__(self, terms):
        self.terms = terms #[term11, term12, ..., term21, term22, ..., ...]
        self.combine_like_terms()

    def combine_like_terms(self):
        combined_terms = []
        for term in self.terms:
            flag = 1
            for ite in combined_terms:
                if term.is_combine(ite):
                    ite.front_coefficient+=term.front_coefficient
                    flag = 0
                    break
            if flag == 1:
                combined_terms.append(term)

        self.terms = [term for term in combined_terms if term.front_coefficient != 0]
        if self.terms:
            self.terms.sort(key=lambda term: term.exponent, reverse=True)
    def __repr__(self):
        return f"Polynomial({self.terms})"

    def __str__(self):
        terms_str = [str(term) for term in self.terms]
        return " + ".join(terms_str)
    
    def __add__(self, other):
        combine_term = []
        combine_term.extend(self.terms)
        combine_term.extend(other.terms)
        return Polynomial(combine_term)

    def __sub__(self, other):
        combine_term = []
        combine_term.extend([Term(term.front_coefficient, [CoffTerm(coeff.word, coeff.exponent, coeff.derivatives) for coeff in term.coefficient], term.exponent) for term in self.terms  if term.exponent > neg_degree])
        combine_term.extend([Term(-term.front_coefficient, [CoffTerm(coeff.word, coeff.exponent, coeff.derivatives) for coeff in term.coefficient], term.exponent) for term in other.terms  if term.exponent > neg_degree])
        return Polynomial(combine_term)
    
    def __mul__(self, other):
        mul_term = []
        if isinstance(other, float) or isinstance(other, int):
            mul_term.extend([Term(other * term.front_coefficient, [CoffTerm(coeff.word, coeff.exponent, coeff.derivatives) for coeff in term.coefficient], term.exponent) for term in self.terms  if term.exponent > neg_degree])
            return Polynomial(mul_term)
        elif isinstance(other, Polynomial) and other.is_single():
            mul_term.extend([Term(other.terms[0].front_coefficient * term.front_coefficient, [CoffTerm(coeff.word, coeff.exponent, coeff.derivatives) for coeff in term.coefficient], term.exponent + other.terms[0].exponent) for term in self.terms  if term.exponent > neg_degree])
            return Polynomial(mul_term)
        elif isinstance(other, Polynomial):
            for sterm in self.terms:
                if sterm.exponent == 0:
                    mul_term.extend([Term(sterm.front_coefficient * term.front_coefficient, 
                                          [CoffTerm(coeff.word, coeff.exponent, coeff.derivatives) for coeff in chain(term.coefficient, sterm.coefficient)], 
                                          term.exponent) 
                                          for term in other.terms  if term.exponent > neg_degree])
                elif sterm.exponent > 0:
                    new_poly = postive_leibnitz_func(other, sterm.exponent)
                    mul_term.extend([
                        Term(sterm.front_coefficient * term.front_coefficient, 
                             [CoffTerm(coeff.word, coeff.exponent, coeff.derivatives) for coeff in chain(term.coefficient, sterm.coefficient)],
                             term.exponent)
                        for term in new_poly.terms if term.exponent > neg_degree])
                elif sterm.exponent < 0:
                    new_poly = negtive_leibnitz_func_order(other, sterm.exponent, default_cal)
                    mul_term.extend([
                        Term(sterm.front_coefficient * term.front_coefficient, 
                             [CoffTerm(coeff.word, coeff.exponent, coeff.derivatives) for coeff in chain(term.coefficient, sterm.coefficient)],
                             term.exponent)
                        for term in new_poly.terms if term.exponent > neg_degree])
        return Polynomial(mul_term)


    def is_single(self):
        if len(self.terms) == 1 and len(self.terms[0].coefficient) == 1 and self.terms[0].coefficient[0].word == "":
            return True
        return False
        
    def derivative(self):
        der_terms = []
        for term in self.terms:
            for coff in term.coefficient:
                if coff.word != "":
                    der_coff = []
                    der_coff = [CoffTerm(coff1.word, coff1.exponent, coff1.derivatives) for coff1 in term.coefficient if coff1 != coff]
                    if coff.exponent-1 != 0:
                        der_coff.append(CoffTerm(coff.word, coff.exponent-1, coff.derivatives))
                    der_coff.append(CoffTerm(coff.word, 1, coff.derivatives+1))
                    der_terms.append(Term(term.front_coefficient * coff.exponent, der_coff, term.exponent))

        return Polynomial(der_terms)

    def derivatives(self, n):
        if n == 0:
            return Polynomial([Term(term.front_coefficient, [CoffTerm(coeff.word, coeff.exponent, coeff.derivatives) for coeff in term.coefficient], term.exponent) for term in (self.terms)]) 
        elif n == 1:
            return self.derivative()
        return self.derivative().derivatives(n-1)
    
    def dell(self, word = "u_0"):
        new_term = []
        for term in self.terms:
            flag = 0
            for coeff in term.coefficient:
                if coeff.word == word:
                    flag = 1
                    break
            if flag == 0:
                new_term.extend([Term(term.front_coefficient, [CoffTerm(coeff.word, coeff.exponent, coeff.derivatives) for coeff in term.coefficient],term.exponent)])
        return Polynomial(new_term)
    
    def remain_pos(self, pos = 0):
        return Polynomial([Term(term.front_coefficient, [CoffTerm(coeff.word, coeff.exponent, coeff.derivatives) for coeff in term.coefficient], term.exponent) for term in (self.terms) if term.exponent >= pos]) 



def times(i, j):
    result = 1
    for t in range(i, j+1):
        result *= t
    return result



def postive_leibnitz_func(poly, n):
    if n == 0:
        result = Polynomial([Term(term.front_coefficient, term.coefficient, term.exponent) for term in (poly.terms) if term.exponent > neg_degree]) 
    elif n > 0:
        result = Polynomial([Term(term.front_coefficient, term.coefficient, term.exponent + n) for term in (poly.terms) if term.exponent > neg_degree]) 
        for i in range(1, n+1):
            terms = [Term(term.front_coefficient, term.coefficient, term.exponent + n-i) for term in poly.derivatives(i).terms]
            der_polys = Polynomial(terms)
            result += der_polys * (times(n-i+1, n) / times(1, i))
    return result

def none_n(n):
    if n % 2 == 0:
        return 1
    else:
        return -1

def negtive_leibnitz_func(poly, items_num):
    result = Polynomial([Term(0, [CoffTerm("", 0, 0)], 0)])
    for i in range(1, items_num+1):
        der_poly = poly.derivatives(i-1)
        result += der_poly * Polynomial([Term(1, [CoffTerm("", 0, 0)], -i)]) * (-1 * none_n(i))
    return result


def negtive_leibnitz_func_order(poly, n, items_num = default_cal):
    if n == 0:
        return Polynomial([Term(term.front_coefficient, term.coefficient, term.exponent) for term in (poly.terms) if term.exponent > neg_degree]) 
    elif n == -1:
        return negtive_leibnitz_func(poly, items_num)
    return negtive_leibnitz_func_order(negtive_leibnitz_func(poly, items_num), n+1, items_num)


terms_L = [Term(1, [CoffTerm("", 0, 0)], 1)]
terms_L.extend([Term(1, [CoffTerm(f"u_{i}", 1, 0)], -i) for i in range(0, 5)])
Poly_L = Polynomial(terms_L)

terms_1 = [Term(1, [CoffTerm("u_1", 1, 0), CoffTerm("u_2", 1, 0)], 1)]
Poly_L = Polynomial(terms_1)
print(Poly_L)
print(Poly_L * Poly_L)


terms_L = [Term(1, [CoffTerm("u_-1", 1, 0)], 1)]
terms_L.extend([Term(1, [CoffTerm("u_-2", 1, 0)], 2)])
terms_L.extend([Term(1, [CoffTerm(f"u_{i}", 1, 0)], -i) for i in range(0, 3)])
Poly_L = Polynomial(terms_L)

print(Poly_L)
print(Poly_L * Poly_L)


#print(poly_L_1 * Poly_L)
#print(poly_L_0 * Poly_L)
#print(poly_L_m1 * Poly_L)

#P_P = Poly_L * Poly_L

#B_1 = Poly_L.remain_pos(pos=1)
#B_2 = P_P.dell().remain_pos(pos=1)
#B_3 = (P_P * Poly_L).dell().remain_pos(pos=1)
#B_4 = (P_P*P_P).dell().remain_pos(pos=1)

#B = [B_1, B_2, B_3, B_4]

#print(B_1)
#print(B_2)
#print(B_3)
#print(B_4)

#Lax_L = [Term(1, [CoffTerm("", 0, 0)], 1)]
#Lax_L.extend([Term(1, [CoffTerm(f"u_{i}", 1, 0)], -i) for i in range(1, 5)])
#Poly_Lax_L = Polynomial(Lax_L)

#print("(p_1 L:)\n", (B_1*Poly_Lax_L - Poly_Lax_L*B_1).dell())
#print("(p_2 L:)\n", (B_2*Poly_Lax_L - Poly_Lax_L*B_2).dell())
#print("(p_3 L:)\n", (B_3*Poly_Lax_L - Poly_Lax_L*B_3).dell())

#phow = [Term(1, [CoffTerm("", 0, 0)], 0)]
#phow.extend([Term(1, [CoffTerm(f"w_{i}", 1, 0)], -i) for i in range(1, 5)])
#Poly_Lax_phow = Polynomial(phow)

#phov = [Term(1, [CoffTerm("", 0, 0)], 0)]
#phov.extend([Term(1, [CoffTerm(f"v_{i}", 1, 0)], -i) for i in range(1, 5)])
#Poly_Lax_phov = Polynomial(phov)

#print("(pho * pho^-1:)\n", Poly_Lax_phow * Poly_Lax_phov)
#print("(L o pho:)\n", Poly_Lax_L * Poly_Lax_phow)