#include <functional>

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

using FilterFunc = std::function<bool(attributetype*, const FilterConditionWithId&)>;

ConditionOp getConditionOp(const std::string &opStr)
{
    if (opStr == "IN") return ConditionOp::IN;
    if (opStr == "==") return ConditionOp::EQ;
    if (opStr == "!=") return ConditionOp::NE;
    if (opStr == ">")  return ConditionOp::GT;
    if (opStr == "<")  return ConditionOp::LT;
    if (opStr == ">=") return ConditionOp::GE;
    if (opStr == "<=") return ConditionOp::LE;
    return ConditionOp::EQ;
}

FilterFunc getFilterFunction(ConditionOp op)
{
    switch (op)
    {
    case ConditionOp::IN:
        return [](attributetype* attribute, const FilterConditionWithId &condition) {
            attributetype value = attribute[condition.attribute_id];
            return condition.attribute_value.find(value) != condition.attribute_value.end();
        };
    
    case ConditionOp::EQ:
        return [](attributetype* attribute, const FilterConditionWithId &condition) {
            attributetype value = attribute[condition.attribute_id];
            attributetype ref_value = *condition.attribute_value.begin();
            return value == ref_value;
        };
    
    case ConditionOp::NE:
        return [](attributetype* attribute, const FilterConditionWithId &condition) {
            attributetype value = attribute[condition.attribute_id];
            attributetype ref_value = *condition.attribute_value.begin();
            return value != ref_value;
        };
    
    case ConditionOp::GT:
        return [](attributetype* attribute, const FilterConditionWithId &condition) {
            attributetype value = attribute[condition.attribute_id];
            attributetype ref_value = *condition.attribute_value.begin();
            return value > ref_value;
        };
    
    case ConditionOp::LT:
        return [](attributetype* attribute, const FilterConditionWithId &condition) {
            attributetype value = attribute[condition.attribute_id];
            attributetype ref_value = *condition.attribute_value.begin();
            return value < ref_value;
        };
    
    case ConditionOp::GE:
        return [](attributetype* attribute, const FilterConditionWithId &condition) {
            attributetype value = attribute[condition.attribute_id];
            attributetype ref_value = *condition.attribute_value.begin();
            return value >= ref_value;
        };
    
    case ConditionOp::LE:
        return [](attributetype* attribute, const FilterConditionWithId &condition) {
            attributetype value = attribute[condition.attribute_id];
            attributetype ref_value = *condition.attribute_value.begin();
            return value <= ref_value;
        };
    
    default:
        return [](attributetype*, const FilterConditionWithId&) { return true; };
    }
}

bool checkCondition(attributetype *attribute, const FilterConditionWithId &condition, const FilterFunc &filterFunc)
{
    return filterFunc(attribute, condition);
}

bool checkConditions(attributetype *attribute, const std::vector<FilterConditionWithId> &filtering_conditions)
{
    if (filtering_conditions.size() == 1)
    {
        const auto &condition = filtering_conditions[0];
        ConditionOp op = getConditionOp(condition.op);
        FilterFunc filterFunc = getFilterFunction(op);
        return filterFunc(attribute, condition);
    }
    
    for (const auto &cond : filtering_conditions)
    {
        ConditionOp op = getConditionOp(cond.op);
        FilterFunc filterFunc = getFilterFunction(op);
        if (!filterFunc(attribute, cond))
        {
            return false;
        }
    }
    return true;
}

class OptimizedFilter
{
private:
    std::vector<std::pair<FilterConditionWithId, FilterFunc>> compiled_conditions_;

public:
    OptimizedFilter(const std::vector<FilterConditionWithId> &filtering_conditions)
    {
        for (const auto &cond : filtering_conditions)
        {
            ConditionOp op = getConditionOp(cond.op);
            compiled_conditions_.emplace_back(cond, getFilterFunction(op));
        }
    }

    bool check(attributetype *attribute) const
    {
        for (const auto &[condition, filterFunc] : compiled_conditions_)
        {
            if (!filterFunc(attribute, condition))
            {
                return false;
            }
        }
        return true;
    }
};