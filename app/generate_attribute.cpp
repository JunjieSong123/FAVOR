#include <iostream>
#include <fstream>
#include <random>
#include <string>
#include <iomanip>

#include "dataset.h"
#include <nlohmann/json.hpp>
using json = nlohmann::json;

int main(int argc, char *argv[])
{
    if (argc < 3)
    {
        std::cerr << "Usage: " << argv[0] << " baseset_path" << " attribute_path\n";
        return 1;
    }

    std::string baseset_path(argv[1]);
    std::string attribute_path(argv[2]);

    // std::string base_attribute_path(argv[3]);

    BaseSet baseset;
    baseset.read_data(baseset_path);

    std::ofstream outFile(attribute_path);
    if (!outFile)
    {
        std::cerr << "Error: Unable to create file at " << attribute_path << std::endl;
        return 1;
    }
    std::cout << "begin generating" << std::endl;

    // random generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> priceDist(0.0, 100.0);
    std::uniform_int_distribution<> intDist(0, 7);
    std::uniform_int_distribution<> boolDist(0, 1);

    outFile << baseset.num << "\n" << 1 << "\n";

    // price
    outFile << "price\n";
    for (int i = 0; i < baseset.num; ++i) {
        outFile << std::fixed << std::setprecision(2) << priceDist(gen) << "\n";
    }

    // // color
    // outFile << "color\n";
    // for (int i = 0; i < baseset.num; ++i) {
    //     outFile << intDist(gen) << "\n";
    // }

    // // gender
    // outFile << "gender\n";
    // for (int i = 0; i < baseset.num; ++i) {
    //     outFile << boolDist(gen) << "\n";
    // }

    outFile.close();
    std::cout << "File successfully generated at: " << attribute_path << std::endl;
    return 0;
}