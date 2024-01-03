#ifndef DATAFRAME_H
#define DATAFRAME_H

#include <tuple>
#include <vector>


#define DATAFRAME_EXPAND    0x01
#define DATAFRAME_TRUNCATE  0x02


namespace DES 
{

template <typename T>
class DataFrame
{

private:
    size_t rows, cols;
    std::vector<std::vector<T>> entries;

public:
    DataFrame<T>(size_t rows, size_t cols);
    std::vector<T> getRow(size_t row);
    std::vector<T> getCol(size_t col);

    void addRow(std::vector<T>& row);

    std::vector<T> operator[] (size_t row);
    
};


template <typename T>
void plot(DataFrame<T>& frame);


}

#endif