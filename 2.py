'''
为了保持一致性 本不对已经生成的多项式进行修改，而是生成新的多项式 所以全部为浅拷贝 只有设计到列表拓展的时候对其进行深拷贝
'''
from itertools import chain
import copy
default_cal = -6
debug = False


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
        

    def print_lex(self):
        words = ""
        deriva = ""
        exponent_text = ""
        if self.derivatives != 0:
            deriva += ","
            deriva += f"{'x' * self.derivatives}"
        
        if self.exponent == 1:
            exponent_text += ""
        else:
            exponent_text += f"^{{{self.exponent}}}"
        
        count_ = self.word.count('_')
        if count_ >= 1:
            first_word, _, second_word = self.word.rpartition("_")
            words += first_word
            words += exponent_text
            words += "_{"
            words += second_word
            words += deriva
            words += "}"

        else:
            words += self.word
            words += exponent_text
            words += "_{"
            words += deriva
            words += "}"

        return words
    

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
        word_sorted = sorted(coefficient, key = lambda coeff:coeff.word)
        derivatives_sorted = sorted(word_sorted, key = lambda coeff:coeff.derivatives)
        return derivatives_sorted
    
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
    
    def print_format_term(self):
        words = ""
        if self.front_coefficient == 0:
            return words
        words += f"{self.front_coefficient}*"

        for coeff in self.coefficient:
            words += f"{coeff}*"

        return words

    def print_lex(self):
        '''
        实现lex格式的打印 此模式为修正系数项符号的打印 所以系数项符号在Term中实现打印
        '''
        words = ""
        if self.front_coefficient == 0:
            return words     

        if self.front_coefficient == 1.0:
            words+= f"+ "
        elif self.front_coefficient == -1.0:
            words += f"- "
        elif self.front_coefficient > 0:
            words+=f"+ {int(self.front_coefficient)}"
        elif self.front_coefficient < 0:
            words+=f"- {-int(self.front_coefficient)}"
        
        for coeff in self.coefficient:
            words += f"{coeff.print_lex()}"

        if self.exponent == 0:
            pass
        elif self.exponent == 1:
            words+= f"\partial"
        else:
            words+= f"\partial^{{{self.exponent}}}"
        
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

        if len(self.terms) == 0:
            self.left = 0
            self.right = 0
        else:
            self.left = self.terms[0].exponent
            self.right = self.terms[-1].exponent


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
    
    def print_format_ploy(self):
        words = ""
        if len(self.terms) == 0:
            return words
        i = self.terms[0].exponent
        words+=f"p^{i:} + {self.terms[0].print_format_term()}"
        for term in self.terms[1:]:
            if term.exponent == i:
                words += f" + {term.print_format_term()}"
            else:
                i = term.exponent
                words += f" + \n\np^{i}: {term.print_format_term()}"
        return words
    
    def print_lex(self):
        words = ""
        for term in self.terms:
            words += f" {term.print_lex()} "
        if words != '' and words[1] == '+':
            words = words[2:]
        return words


    def __add__(self, other):
        combine_terms = []
        combine_terms.extend(copy.deepcopy(self.terms))
        combine_terms.extend(other.terms)
        return Polynomial(combine_terms)

    def __mul__(self, other):
        '''
        实现多项式的乘法，不对其精度做限制
        使用分配率进行依次计算
        '''
        mul_result = Polynomial([])
        for term1 in self.terms:
            for term2 in other.terms:
                mul_result += tmult(term1, term2)
        return mul_result  
 
    def __sub__(self, other):
        return self + Polynomial([Term(-1, [], 0)]) * other

    def remain_order(self, target_order):
        return Polynomial([term for term in self.terms if term.exponent >= target_order]) 

    def derivative(self):
        '''
        对多项式进行求导 返回多项式的导数 返回类型为多项式
        '''
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
        '''
        递归调用多项式求导 求取高阶导数
        '''
        if n == 0:
            return self
        elif n == 1:
            return self.derivative()
        return self.derivative().derivatives(n-1)
    
    def dell(self, word = "u_0"):
        '''
        去掉含有word中的多项式项
        未验证正确与否 
        '''
        new_term = []
        for term in self.terms:
            flag = 0
            for coeff in term.coefficient:
                if coeff.word == word:
                    flag = 1
                    break
            if flag == 0:
                new_term.extend(copy.deepcopy(term))
        return Polynomial(new_term)




def target_mul(poly1, poly2, target):
    '''
    指定精度的多项式乘法计算 返回指定精度的多项式
    '''
    l1, r1, l2, r2 = poly1.left, poly1.right, poly2.left, poly2.right
    mul_result = Polynomial([])
    for term1 in poly1.terms:
        for term2 in poly2.terms:
            mul_result += tmult_targeted(term1, term2, target)
    return mul_result  



    
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
        expland_p_term2 = get_leibniz_plus_term(order1, term2)
        return Polynomial([Term(term1.front_coefficient, term1.coefficient, 0)]) * expland_p_term2
    
    elif order1 < 0:
        expland_p_term2 = get_leibniz_plus_term(order1, term2, -default_cal)
        return Polynomial([Term(term1.front_coefficient, term1.coefficient, 0)]) * expland_p_term2
        
def tmult_targeted(term1, term2, target):
    order1 = term1.exponent
    order2 = term2.exponent

    if order1 == 0:
        if order2 < target:
            return Polynomial([Term(0, [], 0)])
        front_coeffi = term1.front_coefficient * term2.front_coefficient
        coeffi = []
        coeffi.extend(copy.deepcopy(term1.coefficient))
        coeffi.extend(term2.coefficient)
        return Polynomial([Term(front_coeffi, coeffi, order2)])
    
    elif term2.is_single():
        if order1 + order2 < target:
            return Polynomial([Term(0, [], 0)])
        return Polynomial([Term(term1.front_coefficient * term2.front_coefficient, term1.coefficient, order1 + order2)])
    
    elif order1 > 0:
        if order2 >= target:
            expland_p_term2 = get_leibniz_plus_term(order1, term2)
        elif order1 + order2 < target:
            return Polynomial([Term(0, [], 0)])
        else:
            expland_p_term2 = get_leibniz_plus_term(order1, term2, order2+order1+1-target)
        return Polynomial([Term(term1.front_coefficient, term1.coefficient, 0)]) * expland_p_term2
    
    elif order1 < 0:
        if order1 + order2 < target:
            return Polynomial([Term(0, [], 0)])
        else:
            expland_p_term2 = get_leibniz_plus_term(order1, term2, order1 + order2 + 1 - target)
        return Polynomial([Term(term1.front_coefficient, term1.coefficient, 0)]) * expland_p_term2

def times(i, j):
    result = 1
    for t in range(i, j+1):
        result *= t
    return result

def get_leibniz_plus_term(order, term, item_nums = -1):
    if item_nums == -1:
        item_nums = order+1
    result_Poly = Polynomial([])
    for i in range(0, item_nums):
        C = times(order-i+1, order) / times(1, i) 
        result_Poly +=  Polynomial([Term(C, [], 0)]) * Polynomial([term]).derivatives(i) * Polynomial([Term(1, [], order - i)])
    return result_Poly

def get_mul_list(poly_list, target):
    if len(poly_list) == 1:
        return poly_list[0].remain_order(target)
    elif len(poly_list) == 2:
        #print(f"target_mul:\n{poly_list[0].print_lex()}\n{poly_list[1].print_lex()}\n计算精度：{target}")
        return target_mul(poly_list[0], poly_list[1], target)
    else:
        pre = get_mul_list(poly_list[:-1], target-poly_list[0].left)
        #print(f"pre:{pre.print_lex()}*\nlist[-1]:{poly_list[-1].print_lex()}\ntarget:{target}\n")
        return target_mul(pre, poly_list[-1], target)


#terms_L = [Term(1, [CoffTerm("u_-1", 1, 0)], 1)]
#terms_L.extend([Term(1, [CoffTerm("u_-2", 1, 0)], 2)])
#terms_L.extend([Term(1, [CoffTerm(f"u_{i}", 1, 0)], -i) for i in range(0, 3)])
#Poly_L = Polynomial(terms_L)

#print(Poly_L)
#print(Poly_L * Poly_L)
#print(target_mul(Poly_L, Poly_L, -4))

#terms_m1 = Term(1, [], -3)
#terms_f = Term(1, [CoffTerm("f", 1, 0)], 0)
#print(Polynomial([terms_m1]) * Polynomial([terms_f]))
#print(get_leibniz_plus_term(-3, terms_f, 7))



terms_L = [Term(1, [], 1)]
terms_L.extend([Term(1, [CoffTerm(f"u_{i}", 1, 0)], -i) for i in range(1, 10)])
Poly_L = Polynomial(terms_L)

print(Poly_L.print_lex())
L_greater_then = 0
L_debug_greater_then = 1

Poly = [[Poly_L for i in range(0, j)] for j in range (1, 4)]
B_s = []
for i, poly_list in enumerate(Poly):
    B = get_mul_list(poly_list, L_greater_then)
    if debug == True:
        for term in B.terms:
            if term.exponent < L_debug_greater_then:
                for coff in term.coefficient:
                    coff.word += "-"

    print(f"B_{i+1} = {B.print_lex()}") 
    B_s.append(B)

#B_1 = Poly_L.remain_order(0)

#Poly_L_2 = target_mul(Poly_L, Poly_L, target = -4)
#B_2 = target_mul(Poly_L, Poly_L, target=0)
#B_3 = target_mul(Poly_L, Poly_L_2, target=0)

#print(f"B_1:{B_1}")
#print(f"B_2:{B_2}")
#print(f"Poly_L_2:{Poly_L_2}")
#print(f"B_3:{B_3}")

#B_3_list = get_mul_list(Poly_list, 0)
#print(f"B_3_list:{B_3_list}")

Lax_Poly_s = []

for i, B in enumerate(B_s):
    one = target_mul(B, Poly_L, -4)
    two = target_mul(Poly_L, B, -4)
    Lax_Poly = one - two
    print(f"---------------B_{i+1}-----------------")
    print(f"B*L:{one.print_lex()}\nL*B:{two.print_lex()}\n")
    print(f"[B_{i+1}, L] = B_{i+1}*L-L*B_{i+1} = {Lax_Poly.print_lex()}")
    Lax_Poly_s.append(Lax_Poly)






