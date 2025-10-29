#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <ctime>
#include <random>

#include "dataset.h"

void DataSet::read_data(const std::string &dataset_path)
{
    // get the suffix of filename
    std::string::size_type filepos = dataset_path.find_last_of('/') + 1;
    std::string filename = dataset_path.substr(filepos, dataset_path.length() - filepos);
    std::string suffix = filename.substr(filename.find_last_of('.') + 1);

    std::ifstream reader(dataset_path, std::ios::binary | std::ios::ate);

    std::cout << "begin reading " << filename << std::endl;
    size_t filesize = reader.tellg();
    reader.seekg(0, std::ios::beg);

    if (suffix == "fvecs")
    {
        int datasize = 4;
        reader.read((char *)&dim, 4);
        std::cout << "dim = " << dim << std::endl;
        num = filesize / ((dim * datasize) + 4);
        std::cout << "number = " << num << std::endl;
        for (int i = 0; i < num; i++)
        {
            std::vector<float> vector;
            if (i != 0)
                reader.read((char *)&dim, 4);
            for (int j = 0; j < dim; j++)
            {
                float data;
                reader.read((char *)&data, datasize);
                vector.emplace_back(data);
            }
            vectors.emplace_back(vector);
            vector_id.emplace_back(i);
        }
        std::cout << "finish reading!" << std::endl;
    }

    else if (suffix == "bvecs")
    {
        int datasize = 1;
        reader.read((char *)&dim, 4);
        std::cout << "dim = " << dim << std::endl;

        num = filesize / (dim * datasize + 4);
        std::cout << "number = " << num << std::endl;

        for (int i = 0; i < num; i++)
        {
            std::vector<float> float_vector;
            if (i != 0)
                reader.read((char *)&dim, 4);

            for (int j = 0; j < dim; j++)
            {
                unsigned char data;
                reader.read((char *)&data, datasize); 
                float_vector.push_back(static_cast<float>(data)); // 转换为float
            }
            vectors.emplace_back(float_vector); // 存储到总向量集合
            vector_id.emplace_back(i);          // 记录向量ID
        }
        std::cout << "finish reading!" << std::endl;
    }

    reader.close();
}

// attribute = vector_id
// void BaseSet::get_attribute()
// {
//     size_t attribute_size = attribute_num * num;
//     attribute = (float *)malloc(sizeof(float) * attribute_size);
//     for (int i = 0; i < num; i++)
//         for (int j = 0; j < attribute_num; j++)
//             *(attribute + i * attribute_num + j) = i;
// }