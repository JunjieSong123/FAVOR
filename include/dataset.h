#include <string>
#include <vector>
#include <map>

#include "filter_condition.h"

class DataSet
{
public:
    int num;
    int dim;
    int attribute_num;
    std::vector<std::vector<float>> vectors;
    std::vector<int> vector_id;

    DataSet() {};

    void read_data(const std::string &dataset_path);
};

class BaseSet : public DataSet
{
public:
    // attribute[vector_id][attribute_id]
    float *attribute{nullptr}; // attribute[i][j] = attribute + (i * attribute_num + j)

    std::map<std::string, int> attribute_map;

    void get_attribute(const std::string &attribute_path);
};

class QuerySet : public DataSet
{
public:
    std::vector<FilterCondition> filtering_conditions;

    bool check_condition(float *attribute, const FilterCondition& condition, std::map<std::string, int> attribute_map) const;

    bool check_all_conditions(float *attribute, std::vector<FilterCondition>& conditions, std::map<std::string, int> attribute_map) const;
};