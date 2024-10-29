#include "yhm_Matrix.h"
#include <stdexcept>
#include <algorithm>

/// @brief 进行矩阵的初始化，初始化矩阵为内部均为(0.0)的矩阵
/// @param rows 矩阵行数
/// @param cols 矩阵列数
yhm_Matrix::yhm_Matrix(size_t rows, size_t cols) : yrow(rows), ycol(cols), yMaxtrix(rows, std::vector(cols, 0.0))
{
    if(this->yrow == 0 || this->ycol == 0)
    {
        throw std::invalid_argument("Matrix rows and columns cannot be initialized to 0!");
    }
}

/// @brief 初始化的矩阵必须为方阵，
/// @param rows 初始化矩阵行数
/// @param cols 初始化矩阵列数
/// @param is_Identy_Matrix 是否初始化为单位阵(True, 是 | False, 否)
yhm_Matrix::yhm_Matrix(size_t rows, size_t cols, bool is_Identy_Matrix)
{
    if(rows != cols)
    {
        throw std::invalid_argument("Non-square matrices cannot be initialized to the identity matrix!");
    }
    double uvalue = is_Identy_Matrix ? 1.0 : 0.0;

    yrow = rows, ycol = cols;
    yMaxtrix = std::vector<std::vector<double>>(rows, std::vector(cols, 0.0));

    for(size_t i = 0; i<yrow; i++)
        yMaxtrix[i][i] += uvalue;
}

/// @brief 重载构造函数，默认为两行两列的矩阵，内部数值为0.0
yhm_Matrix::yhm_Matrix() : yrow(2), ycol(2), yMaxtrix(2, std::vector<double>(2, 0)){}

/// @brief 析构函数,vector将会自动管理内存，在这里无需进行操作
yhm_Matrix::~yhm_Matrix(){}

/// @brief 为这个矩阵的制定位置进行赋值
/// @param x 元素行数位置，(0 -> Matrix.Row)
/// @param y 元素列数位置，(0 -> Matrix.Col)
/// @param value 矩阵指定位置的元素值
void yhm_Matrix::setElement(size_t x, size_t y, double value) 
{ 
    if(x < 0 || x > yrow || y < 0 || y > ycol)
    {
        throw std::invalid_argument("The number of matrix rows or columns is out of bounds");
        return;
    }
    yMaxtrix[x][y] = value; 
}

/// @brief 获取矩阵某一位置的元素值
/// @param x 元素所在行数
/// @param y 元素所在列数
/// @return 返回当前位置的元素的值
double yhm_Matrix::getElement(size_t x, size_t y) const 
{
    if(x < 0 || x > yrow || y < 0 || y > ycol)
    {
        throw std::invalid_argument("The number of matrix rows or columns is out of bounds");
        return INT_MAX;
    }
    return yMaxtrix[x][y];
}

/// @brief 获取行数
/// @return 返回当前矩阵的行数
size_t yhm_Matrix::getRow() const { return yrow; }

/// @brief 获取列数
/// @return 返回当前矩阵的列数
size_t yhm_Matrix::getCol() const { return ycol; }


void yhm_Matrix::yhm_Print_Matrix()
{
    const size_t row = this->yrow;
    const size_t col = this->ycol;
    const auto &Matrix = this->yMaxtrix;

    //打印矩阵
    for(size_t i = 0; i<row; i++)
    {
        for(size_t j = 0; j<col; j++)
            printf("%.3lf\t", yMaxtrix[i][j]);
        std::cout << std::endl;
    }
}

/// @brief 求解当前矩阵的转置矩阵
yhm_Matrix yhm_Matrix::yhm_Transpose_Matrix()
{

    const size_t ori_row = this->yrow;
    const size_t ori_col = this->ycol;
    const size_t row = this->ycol;
    const size_t col = this->yrow;

    if(row == 0 || col == 0)
    {
        throw std::runtime_error("The matrix cannot be transposed, the number of rows or columns is 0");
    }

    yhm_Matrix Result_Matrix = yhm_Matrix(row, col);

    auto &Matrix = Result_Matrix.yMaxtrix;
    auto &ori_Matrix = this->yMaxtrix;

    for(size_t i = 0; i<ori_row; i++)
        for(size_t j = 0; j<ori_col; j++)
            Matrix[j][i] = ori_Matrix[i][j];

    return Result_Matrix;
}

/// @brief 求解当前矩阵的行列式
/// @return 返回当前矩阵的行列式值
double yhm_Matrix::yhm_Cal_Det()
{
    if(this->yrow == 0 || this->ycol == 0 || this->ycol != this->yrow)
    {
        throw std::runtime_error("This matrix is ​​not valid, the number of rows or columns is 0 or it is not a square matrix!");
        return -1;
    }

    if(this->ycol == 1 && this->yrow == 1)
        return this->yMaxtrix[0][0];

    if(this->ycol == 2 && this->yrow == 0)
    {
        double res = 0;
        res = this->yMaxtrix[0][0] * this->yMaxtrix[1][1] - this->yMaxtrix[0][1] * this->yMaxtrix[1][0];
        return res;
    }

    //默认第零行
    int Max_Zero = 0;
    size_t Best_Row = 0;
    for(size_t i = 0; i<yrow; i++)
    {
        int sum_zero = 0;
        for(size_t j = 0; j<ycol; j++)
            if(yMaxtrix[i][j] == 0)
                sum_zero ++;
        Best_Row = Max_Zero < sum_zero ? i : Best_Row;
        Max_Zero = std::max(Max_Zero, sum_zero);
    }

    double res = 0;

    const auto col = this->ycol;
    const auto row = this->yrow;
    const auto &Matrix = this->yMaxtrix;

    for(size_t i = 0; i<row; i++)
    {
        for(size_t j = 0; j<col; j++)
        {
            //剪枝，为0可以直接忽略
            if(yMaxtrix[i][j] == 0) continue;
            auto cosine = getAlgebraic_cofactor(Best_Row, j);
            res += pow(-1, i + j) * yMaxtrix[i][j] * cosine.yhm_Cal_Det();
        }
    }

    return res;
}

/// @brief 返回对应的代数余子式对应的矩阵
/// @param x 某一行
/// @param y 某一列
/// @return 返回一个矩阵
yhm_Matrix yhm_Matrix::getAlgebraic_cofactor(size_t x, size_t y)
{
    size_t col = 0;
    size_t row = 0;

    yhm_Matrix res(yrow - 1, ycol - 1);

    for(size_t i = 0; i<yrow; i++)
    {
        if(i == x)  continue;
        col = 0;
        for(size_t j = 0; j<ycol; j++)
        {
            if(j == y)  continue;
            res.setElement(row, col++, yMaxtrix[i][j]);
        }
        row++;
    }

    return res;
}


/// @brief 求解逆矩阵(使用的是代数余子式求解逆矩阵)
/// @return 返回当前矩阵的逆矩阵
yhm_Matrix yhm_Matrix::getInverse()
{
    // 在这里无需判定矩阵是否行数或者列数为0 因为矩阵无法初始化为0行0列的矩阵 这样的矩阵在这里是不被允许创建的
    // 判定是否为方阵
    if(ycol != yrow)
    {
        throw std::runtime_error("This matrix is ​​not a square matrix and cannot be inverted!");
    }

    // 判定是否为奇异矩阵
    // 行列式为零，则当前的矩阵为奇异矩阵，无法进行求逆操作
    if(this->yhm_Cal_Det() == 0)
    {
        throw std::runtime_error("This matrix is ​​a singular matrix and cannot be inverted!");
    }

    // 求解余子式的转置矩阵进而求得伴随矩阵
    auto res = this->getConfactor_Matrix().yhm_Transpose_Matrix();
    auto rec_det = 1.0/this->yhm_Cal_Det();

    // 在最后计算结果的时候再加上符号
    for(size_t i = 0; i<yrow; i++)
        for(size_t j = 0; j<ycol; j++)
            res.setElement(i, j, pow(-1, i + j) * rec_det * res.getElement(i, j));

    return res;
}

/// @brief 求解当前矩阵的余子式矩阵
/// @return 返回当前矩阵的余子式矩阵
yhm_Matrix yhm_Matrix::getConfactor_Matrix()
{
    // 只可对方阵进行余子式的求解
    if(ycol != yrow)
    {
        throw std::invalid_argument("Unable to solve cofactors for non-square matrices!");
    }

    yhm_Matrix res(yrow, ycol);

    for(size_t i = 0; i<yrow; i++)
        for(size_t j = 0; j<ycol; j++)
            res.setElement(i, j, this->getAlgebraic_cofactor(i, j).yhm_Cal_Det());

    return res;
}

/// @brief 返回两矩阵相加后的矩阵
/// @param b 作为右侧矩阵
/// @return 返回一个同类型的矩阵
yhm_Matrix yhm_Matrix::yhm_Add(yhm_Matrix b)
{
    auto &a = yMaxtrix;
    if(yrow != b.getRow() || ycol != b.getCol())
    {
        throw std::runtime_error("The two matrices cannot be added, and the number of rows or columns is not the same!");
    }

    yhm_Matrix res(yrow, ycol);

    for(size_t i = 0; i<yrow; i++)
        for(size_t j = 0; j<ycol; j++)
            res.setElement(i, j, a[i][j] + b.yMaxtrix[i][j]);

    return res;
}

/// @brief 返回两矩阵相减的结果
/// @param b 作为左侧矩阵
/// @return 返回一个同类型的矩阵
yhm_Matrix yhm_Matrix::yhm_Sub(yhm_Matrix b)
{
    auto &a = yMaxtrix;
    if(yrow != b.getRow() || ycol != b.getCol())
    {
        throw std::runtime_error("The two matrices cannot be added, and the number of rows or columns is not the same!");
    }

    yhm_Matrix res(yrow, ycol);

    for(size_t i = 0; i<yrow; i++)
        for(size_t j = 0; j<ycol; j++)
            res.setElement(i, j, a[i][j] - b.yMaxtrix[i][j]);

    return res;
}

/// @brief 矩阵相乘，所传入的矩阵为右矩阵
/// @param b 右矩阵
/// @return 返回两矩阵相乘的结果
yhm_Matrix yhm_Matrix::yhm_Mul(yhm_Matrix b)
{
    if(ycol != b.getRow())
    {
        throw std::runtime_error("Two matrices cannot be multiplied!");
    }

    auto a = yMaxtrix;
    yhm_Matrix res(yrow, b.getCol());

    for(size_t i = 0; i<yrow; i++)
    {
        for(size_t j = 0; j<b.getCol(); j++)
        {
            double sum = 0;
            for(size_t k = 0; k<ycol; k++)
                sum += a[i][k] * b.yMaxtrix[k][j];
            res.setElement(i, j, sum);
        }
    }

    return res;
}

int yhm_Matrix::yhm_GetRank()
{
    size_t row = yrow, col = ycol;
    yhm_Matrix a(yrow, ycol);
    a.yMaxtrix = yMaxtrix;
    auto &Ma = a.yMaxtrix;

    for(size_t i = 0; i<row - 1; i++)
    {
        if(Ma[i][i] == 0)
        {
            size_t cur = 0;
            for(size_t j = i + 1; j<row; j++)
                if(Ma[j][i] != 0)
                {
                    cur = j;
                    break;
                }
            a.yhm_RowTransform(i, cur);
        }
        if(Ma[i][i] == 0)   continue;

        double factor = Ma[i][i];
        for(size_t j = i + 1; j<row; j++)
        {
            double div = Ma[j][i]/factor;
            for(size_t k = i; k<col; k++)
            {
                Ma[j][k] -= Ma[i][k] * div;
            }
        }
    }

    int rank = 0;
    for(const auto row : Ma)
    {
        for(const auto col : row)    
        {
            if(std::abs(col) > 1e-6)
            {
                rank++;
                break;
            }
        }
    }

    a.yhm_Print_Matrix();

    return rank;
}

void yhm_Matrix::yhm_RowTransform(size_t x, size_t y)
{
    if(x < 0 || x > yrow || y < 0 || y > yrow)
    {
        throw std::invalid_argument("The number of rows or columns exceeds the limit!");
    }

    if(x == y)  return;

    for(size_t j = 0; j<ycol; j++)
        std::swap(yMaxtrix[x][j], yMaxtrix[y][j]);
}

/// @brief 乘号运算符重载
/// @param other 右侧矩阵
/// @return 
yhm_Matrix yhm_Matrix::operator*(const yhm_Matrix& other)
{
    return this->yhm_Mul(other);
}

/// @brief 加号运算符重载
/// @param other 
/// @return 
yhm_Matrix yhm_Matrix::operator+(const yhm_Matrix& other)
{
    return this->yhm_Add(other);
}

/// @brief 减号运算符重载
/// @param other 
/// @return 
yhm_Matrix yhm_Matrix::operator-(const yhm_Matrix& other)
{
    return this->yhm_Sub(other);
}

/// @brief 求解矩阵的逆矩阵
/// @return 
yhm_Matrix yhm_Matrix::operator~()
{
    return this->getInverse();
}

yhm_Matrix yhm_Matrix::operator!()
{
    return this->yhm_Transpose_Matrix();
}