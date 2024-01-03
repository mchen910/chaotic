#ifndef DATAFRAME_H
#define DATAFRAME_H

#include <tuple>
#include <vector>


namespace DES 
{

template <typename T>
class DataFrame
{

private:
    size_t numRows, numCols;
    std::vector<std::vector<T>> entries;

public:
    DataFrame(size_t rows, size_t cols);
    std::vector<T> getRow(size_t row);
    std::vector<T> getCol(size_t col);

    void addRow(std::vector<T> row);
    void addCol(std::vector<T> col);
    
};


}

#endif