:blush: Hello , I`m yhm.

![](https://img.shields.io/badge/Matrix-Calculate-green) 								![](https://img.shields.io/badge/Matrix-inverse_matrix-r)						 ![](https://img.shields.io/badge/Matrix-Inverse_Matrix-orange) 

​			![](https://img.shields.io/badge/Matrix-Det-purple) 												![](https://img.shields.io/badge/Matrix-Rank-pink)							![](https://img.shields.io/badge/Matrix-vector-blue)

:penguin: 因为某些原因需要用到矩阵运算，因此我写了一个较为简单的**矩阵类**，方便后续的学习。

:herb: 主要包括以下几种操作与运算：

- :traffic_light: <font color = red>*矩阵的初始化*</font>
- :traffic_light: <font color = orange>*单位矩阵的初始化*</font>
- :traffic_light: <font color = yellow>*矩阵的打印操作*</font>
- :traffic_light: <font color = green>*求解转置矩阵*</font>
- :traffic_light: <font color = blue>*矩阵行列式求解*</font>
- :traffic_light: <font color = purple>*矩阵的代数余子式求解*</font>
- :traffic_light: <font color = yellowgreen>*矩阵的求逆运算*</font>
- :traffic_light: <font color = pink>*矩阵的加减乘基础运算*</font>
- :traffic_light: <font color = blackred>*矩阵的初等行变换*</font>

------

:triangular_ruler: <font color = red>**How To Use ? **</font>

**:one:** <font color = yellow>**矩阵求逆**</font>

````c++
yhm_Matrix a(3, 3);
```
    your code to give a  a value
```
auto b = ~a;

b即为a逆矩阵
````

**:two:** <font color = yellow>**矩阵求转置**</font>

```c++
auto b = !a;

b为a的转置矩阵
```

**:three:** <font color = yellow>**矩阵加法**</font>

```c++
auto c = a + b;

c为a与b矩阵相加的结果
```

**:four:** <font color = yellow>**矩阵减法**</font>

```c++
auto c = a - b

c为a与b矩阵相减的结果
```

**:five:** <font color = yellow>**矩阵乘法**</font>

```c++
auto c = a * b;

c为矩阵a与b相乘的运算结果
```

**:six:** <font color = yellow>**矩阵打印**</font>

```c++
a = yhm_Matrix(3, 3, true);

std::cout << a << std::endl;
```









