
'''
为了保持一致性 本代码不对已经生成的多项式进行修改，而是生成新的多项式 所以全部为浅拷贝 节省内存 尝试
'''
from itertools import chain
import copy
default_cal = -5

class CoffTerm:
    '''
    实现算子的系数项的单项
       exponent
    word                  
       derivatives

    Atteributes:
        word:单项的字母
        exponent:单项的指数
        derivatives:单项的对变量导数的次数
    '''
    def __init__(self, word = "", exponent = 0, derivatives = 0):
        self.exponent = exponent
        self.word = word
        self.derivatives = derivatives

    def __repr__(self) -> str:
        return f"Term({self.word}, {self.exponent}, {self.derivatives})"

    def __eq__(self, other):
        if isinstance(other, CoffTerm):
            return self.word == other.word and self.exponent == other.exponent and self.derivatives == other.derivatives
        return False
    
    def __str__(self):
        '''
        返回单项的字符串表示
        在生成系数项的时候会对单项中为空的单项进行去除，所以不需要考虑单项为空时的输出情况   
        '''
        deriva = ""
        if self.derivatives != 0:
            deriva = ","
            deriva += f"x({self.derivatives})"
        if self.exponent == 1:
            return f"{self.word}{deriva}"
        else:
            return f"{self.word}{deriva}^{self.exponent}"

    def is_null(self):
        '''
        判断单项是否为空 定义为空的时候单项为1 即为乘法单位元
        '''
        return self.word == "" or self.exponent == 0
    
class Term:
    '''
    实现算子的系数项 只有前系数为0的时候才会为空 系数单项列表为空时候为1 
    front_cofficient   cofficient    exponent
        a             [CT, CT, ..]     exp
    Attributes:
        front_coefficient:系数项前的系数
        coefficient:系数项的单项 包含系数项单项的列表 
        exponent:系数项的指数
    '''
    def __init__(self, front_coefficient = 0, coefficient=[], exponent = 0):
        self.front_coefficient = float(front_coefficient)
        self.exponent = exponent
        temp_coefficient = self.combine_same_coefficient(coefficient)
        temp_coefficient = self.clear_null_coefficient(temp_coefficient)
        temp_coefficient = self.sort_coefficient(temp_coefficient)
        self.coefficient = temp_coefficient

    @staticmethod
    def clear_null_coefficient(coefficient):
        '''
        清除系数项中的空单项
        '''
        return [copy.deepcopy(term) for term in coefficient if not term.is_null()]
    
    @staticmethod
    def combine_same_coefficient(coefficient):
        '''
        合并系数项中相同的单项 注意此时要将相同单项的系数相加
        '''
        coeff = []
        for coeff_itr in coefficient:
            flag = 1
            for itr in coeff:
                if coeff_itr.word == itr.word and coeff_itr.derivatives == itr.derivatives:
                    itr.exponent += coeff_itr.exponent
                    flag = 0
                    break
            if flag == 1:
                coeff.append(coeff_itr)
        return coeff
    
    @staticmethod
    def sort_coefficient(coefficient):
        return sorted(coefficient, key = lambda coeff:coeff.word)
    
    def __repr__(self):
        return f"Term({self.front_coefficient}, {self.coefficient}, {self.exponent})"
    
    def __str__(self):
        words = ""
        if self.front_coefficient == 0:
            return words
        words += f"{self.front_coefficient}*"

        for coeff in self.coefficient:
            words += f"{coeff}*"

        #此处必定会留存一个*， 所以只需要补充一个*即可
        if self.exponent == 1:
            words += f"p"
        elif self.exponent == 0:
            words = words[:-1]
        else:
            words += f"p^{self.exponent}"                       
        return words

    def is_combine(self, other):
        '''
        判断两个项是否可以合并
        '''
        if isinstance(other, Term):
            return self.coefficient == other.coefficient and self.exponent == other.exponent
        
    def is_null(self):
        return self.front_coefficient == 0
    
    def is_single(self):
        return self.coefficient == []
    
class Polynomial:
    '''
    实现多项式项，多项式项是由系数项组成的
    T1 + T2 + T3 + ...
    Attributes:
    terms:系数项的列表 [T1, T2, T3, ...]
    '''
    def __init__(self, terms = []):
        temp_terms = self.combine_same_term(terms)
        temp_terms = self.clear_null_term(temp_terms)
        temp_terms = self.sort_terms(temp_terms)
        self.terms = temp_terms

    @staticmethod
    def clear_null_term(terms):
        '''
        清除多项式中的空项
        '''
        return [term for term in terms if not term.is_null()]
    @staticmethod
    def combine_same_term(terms):
        '''
        合并相同的项，与合并系数项单项的方法类似，但是这里是前系数相加
        '''    
        combined_terms = []
        for term in terms:
            flag = 1
            for itr in combined_terms:
                if term.is_combine(itr):
                    itr.front_coefficient += term.front_coefficient
                    flag = 0
                    break
            if flag == 1:
                combined_terms.append(term)

        return combined_terms
    @staticmethod
    def sort_terms(terms):
        return sorted(terms, key=lambda term: term.exponent, reverse=True)
    
    def __repr__(self):
        return f"Polynomial({self.terms})"
    
    def __str__(self):
        words = ""
        for term in self.terms:
            words += f"{term} + "
        return words[:-3]
    
    def __add__(self, other):
        combine_terms = []
        combine_terms.extend(copy.deepcopy(self.terms))
        combine_terms.extend(other.terms)
        return Polynomial(combine_terms)

    def __mul__(self, other):
        mul_result = Polynomial([])
        for term1 in self.terms:
            for term2 in other.terms:
                mul_result += tmult(term1, term2)
        return mul_result   

    def remain_order(self, target_order):
        return Polynomial([term for term in self.terms if term.exponent >= target_order]) 

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
            return self
        elif n == 1:
            return self.derivative()
        return self.derivative().derivatives(n-1)
    
    
def tmult(term1, term2):
    order1 = term1.exponent
    order2 = term2.exponent

    if order1 == 0:
        front_coeffi = term1.front_coefficient * term2.front_coefficient
        coeffi = []
        coeffi.extend(copy.deepcopy(term1.coefficient))    
        coeffi.extend(term2.coefficient)
        return Polynomial([Term(front_coeffi, coeffi, order2)])
    
    elif term2.is_single():
        return Polynomial([Term(term1.front_coefficient * term2.front_coefficient, term1.coefficient, order1 + order2)])

    elif order1 > 0:
        expland_p_term2 = get_postive_leibniz_plus_term(order1, term2)
        return Polynomial([Term(term1.front_coefficient, term1.coefficient, 0)]) * expland_p_term2
    
    elif order1 < 0:
        expland_p_term2 = get_negative_leibniz_times_poly(order1, Polynomial([term2]))
        return Polynomial([Term(term1.front_coefficient, term1.coefficient, 0)]) * expland_p_term2
        

def times(i, j):
    result = 1
    for t in range(i, j+1):
        result *= t
    return result


def get_postive_leibniz_plus_term(order, term):
    result_Poly = Polynomial([])
    for i in range(0, order+1):
        C = times(order-i+1, order) / times(1, i) 
        result_Poly +=  Polynomial([Term(C, [], 0)]) * Polynomial([term]).derivatives(i) * Polynomial([Term(1, [], order - i)])
    return result_Poly


def none_n(n):
    if n % 2 == 0:
        return 1
    else:
        return -1


def negtive_leibnitz_func(poly, items_num):
    result = Polynomial([])
    for i in range(1, items_num+1):
        der_poly = poly.derivatives(i-1)
        result += Polynomial([Term((-1 * none_n(i)), [], 0)]) * der_poly * Polynomial([Term(1, [], -i)])
    return result


def get_negative_leibniz_times_poly(order, poly, items_num = -default_cal):
    if order == 0:
        return poly
    elif order == -1:
        return negtive_leibnitz_func(poly, items_num)
    return get_negative_leibniz_times_poly(order+1, negtive_leibnitz_func(poly, items_num), items_num)       

def get_target_order(poly1, poly2, target_order):
    pass


def specified_mul(poly1, poly2, target_order = -5):
    order1, order2 = get_target_order(poly1, poly2, target_order)
    return (poly1.remain_order(order1)) * (poly2.remain_order(order2))


terms_L = Term(1, [CoffTerm("a", 1, 0), CoffTerm("b", 1, 0)], -1)
Poly_L = Polynomial([terms_L])

#print(Poly_L)
#print(Poly_L * Poly_L)

terms_m1 = Term(1, [], -3)
terms_f = Term(1, [CoffTerm("f", 1, 0)], 0)
print(Polynomial([terms_m1]) * Polynomial([terms_f]))


        



