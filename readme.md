```
CoffTerm(word, exponent, derivatives)
    """
        exp
    word, der
    """


Term(front_cofficient, coefficient, exponent)
    """
    front_cofficient   cofficient    exponent
           a           [CT, CT, ..]   exp
          exp    
    a CT p
    """


Polynomial(terms)
    """
    T1 + T2 + T3 + ...
        terms    
    [T1, T2, T3, ...]
    """
```

下面的两个参量需要按计算需求设定，如果算不对，就开高点
default_cal = 5 负Leibnitz 规则计算长度
neg_degree = -5 全环境计算深度

多项式由多项式项们(terms): terms = [Term, Term, Term, ...]

多项式(Term)项: 前系数(front_cofficient) 
                符号系数(coefficient = [CoffTerm, CoffTerm, ...]) 
                指数项系数(exponent)

符号系数(CoffTerm): 符号字符(word)
                    符号字符的系数(exponent)
                    符号字符对x的导数次数(derivatives)

例如:要生成多项式 $\partial = \partial + u_0 + u_1\partial^{-1} + u_2\partial^{-2} + \cdots$

从-2次项截断
terms_L = [Term(1, [], 1)] 一次项
terms_L.extend([Term(1, [CoffTerm(f"u_{i}", 1, 0)], -i) for i in range(0, 3)]) 0到(-2) 次项

Poly_L = Polynomial(terms_L) 从项列表生成多项式

结果解释：
6*u_1,x(2)*u_3*p**-4 即为
$6u_1,_{xx}u_3\partial^{-4}$

多项式中支持的运算
Poly * Poly 
Poly.remain_pos(pos) 保留到pos次
Poly.dell(word) 将word的值设为0后的多项式  对应于u_0=0


已经实现了可变精度计算 
具体使用方法：
get_mul_list(poly_list, target)
poly_list 接受一个由多项式组成的列表，返回列表中多项式依次相乘的结果, 计算顺序为从左往右
target为计算截至项

例如要计算 $B_3$ 则需要计算 $(L*L*L)_{>0}$ 
首先构造L 为 
terms_L = [Term(1, [], 1)] 
terms_L.extend([Term(1, [CoffTerm(f"u_{i}", 1, 0)], -i) for i in range(1, 10)])
Poly_L = Polynomial(terms_L)

然后构造多项式列表 poly_list = [Poly, Poly, Poly]
接着使用get_mul_list得到 B_3 = get_mul_list(poly_list, 0)


已经加入了lex模式显示 调用 .print_lex() 即可以使用lex代码显示