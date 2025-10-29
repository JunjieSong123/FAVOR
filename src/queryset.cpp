#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <ctime>
#include <random>

#include "dataset.h"

bool QuerySet::check_condition(float *attribute, const FilterCondition& condition, std::map<std::string, int> attribute_map) const
{
    auto it = attribute_map.find(condition.attribute_name);
    if (it == attribute_map.end()) return false;
    
    int attribute_id = it->second;
    float value = attribute[attribute_id];
    
    if (condition.op == "IN") {
        return condition.attribute_value.find(value) != condition.attribute_value.end();
    } 
    else if (condition.attribute_value.size() == 1) {
        float ref_value = *condition.attribute_value.begin();
        if (condition.op == "==") return value == ref_value;
        if (condition.op == "!=") return value != ref_value;
        if (condition.op == ">")  return value > ref_value;
        if (condition.op == "<")  return value < ref_value;
        if (condition.op == ">=") return value >= ref_value;
        if (condition.op == "<=") return value <= ref_value;
    }
    return false;
}

bool QuerySet::check_all_conditions(float* attribute, std::vector<FilterCondition>& conditions, std::map<std::string, int> attribute_map) const
{
    for (const auto& cond : conditions) {
        if (!check_condition(attribute, cond, attribute_map)) {
            return false;
        }
    }
    return true;
}