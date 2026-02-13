#pragma once

#include <functional>
#include <vector>
#include <set>
#include <string>
#include <algorithm>
#include <cstdint>
#include <cstring>

#include "filter_condition.h"

typedef float attributetype;

enum class ConditionOp : uint8_t
{
    IN,
    EQ,
    NE,
    GT,
    LT,
    GE,
    LE
};

struct alignas(16) FastCondition {
    uint32_t attribute_id;
    ConditionOp op;
    attributetype ref_value;
    uint32_t in_values_count;
    const attributetype* in_values_ptr;

    FastCondition() = default;

    explicit FastCondition(const FilterConditionWithId& cond, std::vector<attributetype>& in_vals_storage, size_t& storage_offset) {
        attribute_id = static_cast<uint32_t>(cond.attribute_id);
        op = getConditionOp(cond.op);

        if (op == ConditionOp::IN) {
            in_values_count = static_cast<uint32_t>(cond.attribute_value.size());
            if (in_values_count > 0) {
                in_values_ptr = &in_vals_storage[storage_offset];
                for (const auto& val : cond.attribute_value) {
                    in_vals_storage[storage_offset++] = val;
                }
            } else {
                in_values_ptr = nullptr;
            }
            ref_value = 0;
        } else {
            in_values_count = 0;
            in_values_ptr = nullptr;
            ref_value = cond.attribute_value.empty() ? 0 : *cond.attribute_value.begin();
        }
    }

    inline __attribute__((always_inline)) bool check(const attributetype* attribute) const {
        attributetype value = attribute[attribute_id];

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
            case ConditionOp::IN: {
                if (in_values_count == 0) return false;
                for (uint32_t i = 0; i < in_values_count; ++i) {
                    if (in_values_ptr[i] == value) return true;
                }
                return false;
            }
            default:
                return true;
        }
    }

private:
    static ConditionOp getConditionOp(const std::string &opStr) {
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
    std::vector<FastCondition> conditions_;
    std::vector<attributetype> in_values_storage_;

public:
    OptimizedFilter() = default;

    OptimizedFilter(const std::vector<FilterConditionWithId> &filtering_conditions) {
        size_t total_in_values = 0;
        for (const auto &cond : filtering_conditions) {
            if (cond.op == "IN") {
                total_in_values += cond.attribute_value.size();
            }
        }

        in_values_storage_.reserve(total_in_values);
        conditions_.reserve(filtering_conditions.size());

        size_t storage_offset = 0;
        for (const auto &cond : filtering_conditions) {
            conditions_.emplace_back(cond, in_values_storage_, storage_offset);
        }
    }

    inline __attribute__((always_inline)) bool check(const attributetype* attribute) const {
        for (const auto &condition : conditions_) {
            if (!condition.check(attribute)) {
                return false;
            }
        }
        return true;
    }

    size_t conditionCount() const {
        return conditions_.size();
    }
};
