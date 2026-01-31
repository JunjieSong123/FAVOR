#include <functional>
#include <vector>
#include <set>
#include <string>
#include <algorithm>

typedef float attributetype;

enum class ConditionOp
{
    IN,
    EQ,
    NE,
    GT,
    LT,
    GE,
    LE
};

struct ConditionWrapper {
    int attribute_id;
    ConditionOp op;
    attributetype ref_value;
    std::vector<attributetype> in_values;

    explicit ConditionWrapper(const FilterConditionWithId& cond) {
        attribute_id = cond.attribute_id;
        op = getConditionOp(cond.op);
        
        if (op == ConditionOp::IN) {
            in_values.reserve(cond.attribute_value.size());
            for (const auto& val : cond.attribute_value) {
                in_values.push_back(val);
            }
        } else {
            if (!cond.attribute_value.empty()) {
                ref_value = *cond.attribute_value.begin();
            } else {
                ref_value = 0;
            }
        }
    }

    bool check(attributetype value) const {
        switch (op) {
            case ConditionOp::EQ:
                return value == ref_value;
                
            case ConditionOp::NE:
                return value != ref_value;
                
            case ConditionOp::GT:
                return value > ref_value;
                
            case ConditionOp::LT:
                return value < ref_value;
                
            case ConditionOp::GE:
                return value >= ref_value;
                
            case ConditionOp::LE:
                return value <= ref_value;
                
            case ConditionOp::IN:
                for (auto v : in_values) {
                    if (value == v) {
                        return true;
                    }
                }
                return false;
                
            default:
                return true;
        }
    }
    
private:
    ConditionOp getConditionOp(const std::string &opStr) const
    {
        if (opStr == "==") return ConditionOp::EQ;
        if (opStr == "<")  return ConditionOp::LT;
        if (opStr == ">")  return ConditionOp::GT;
        if (opStr == ">=") return ConditionOp::GE;
        if (opStr == "<=") return ConditionOp::LE;
        if (opStr == "!=") return ConditionOp::NE;
        if (opStr == "IN") return ConditionOp::IN;
        return ConditionOp::EQ;
    }
};

class OptimizedFilter
{
private:
    std::vector<ConditionWrapper> compiled_conditions_;

public:
    OptimizedFilter(const std::vector<FilterConditionWithId> &filtering_conditions)
    {
        compiled_conditions_.reserve(filtering_conditions.size());
        for (const auto &cond : filtering_conditions)
        {
            compiled_conditions_.emplace_back(cond);
        }
    }

    bool check(attributetype *attribute) const
    {
        for (const auto &condition : compiled_conditions_)
        {
            attributetype value = attribute[condition.attribute_id];
            if (!condition.check(value))
            {
                return false;
            }
        }
        return true;
    }
};