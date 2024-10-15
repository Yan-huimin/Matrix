#include "yhm_Matrix.h"

int main()
{
    const std::string filename = "source.txt";
    std::ifstream  FILE(filename);

    if(!FILE.is_open())
    {
        throw std::runtime_error("Unable to open file : " + filename);
    }

    std::vector<std::vector<double>> cur;

    std::string line;
    while(std::getline(FILE, line))
    {
        std::stringstream ss(line);
        std::vector<double> tmp;

        double value;
        while(ss >> value){ tmp.push_back(value); }
        cur.push_back(tmp);
    }

    size_t row = cur.size(), col = cur[0].size();
    // yhm_Matrix a(row, col);

    yhm_Matrix a(row, col);

    std::cout << "Ori_Matrix" << std::endl;
    std::cout << "********************************Debug*********************************" << std::endl;
    for(size_t i = 0; i<row; i++)
        for(size_t j = 0; j<col; j++)
            a.setElement(i, j, cur[i][j]);
    a.yhm_Print_Matrix();

    // std::cout << "Det_Value" << std::endl;
    // std::cout << "********************************Debug*********************************" << std::endl;
    // auto value = a.yhm_Cal_Det();
    // std::cout << value << std::endl;

    // std::cout << "Inverse_Matrix" << std::endl;
    // std::cout << "********************************Debug*********************************" << std::endl;
    // auto inv_ma = a.getInverse();
    // inv_ma.yhm_Print_Matrix();

    std::cout << "Rank" << std::endl;
    std::cout << "********************************Debug*********************************" << std::endl;
    int rank = a.yhm_GetRank();
    std::cout << rank << std::endl;

    return 0;
}