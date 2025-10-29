#pragma once

#include "hnswlib/hnswlib.h"
#include "hnswlib/hnswalg.h"
#include "filter_condition.h"

using namespace hnswlib;
namespace favor
{
    typedef float attributetype;

    template <typename dist_t>
    class FAVOR : public HierarchicalNSW<dist_t>
    {
    public:
        using HierarchicalNSW<dist_t>::searchKnn;
        using HierarchicalNSW<dist_t>::addPoint;
        static const tableint MAX_LABEL_OPERATION_LOCKS = 65536;
        size_t num_attribute_{0};
        size_t attribute_offset_{0};
        std::mutex delta_mutex;
        dist_t delta_d = 0.0;
        const dist_t LARGE_DIST = 100000.0;

        FAVOR(
            SpaceInterface<dist_t> *s,
            const std::string &location,
            bool nmslib = false,
            size_t max_elements = 0,
            size_t ef = 100,
            bool allow_replace_deleted = false)
            : HierarchicalNSW<dist_t>(max_elements, allow_replace_deleted)
        {
            loadIndex(location, s, max_elements, ef);
        }

        FAVOR(
            SpaceInterface<dist_t> *s,
            size_t max_elements,
            size_t M = 16,
            size_t ef_construction = 200,
            size_t num_attribute = 1,
            size_t random_seed = 100,
            bool allow_replace_deleted = false) : HierarchicalNSW<dist_t>(max_elements, allow_replace_deleted)
        {
            this->max_elements_ = max_elements;
            this->num_deleted_ = 0;
            this->data_size_ = s->get_data_size();
            this->fstdistfunc_ = s->get_dist_func();
            this->dist_func_param_ = s->get_dist_func_param();
            if (M <= 10000)
            {
                this->M_ = M;
            }
            else
            {
                HNSWERR << "warning: M parameter exceeds 10000 which may lead to adverse effects." << std::endl;
                HNSWERR << "         Cap to 10000 will be applied for the rest of the processing." << std::endl;
                this->M_ = 10000;
            }
            this->maxM_ = this->M_;
            this->maxM0_ = this->M_ * 2;
            this->ef_construction_ = std::max(ef_construction, this->M_);
            this->ef_ = 100;
            num_attribute_ = num_attribute;

            this->level_generator_.seed(random_seed);
            this->update_probability_generator_.seed(random_seed + 1);

            this->size_links_level0_ = this->maxM0_ * sizeof(tableint) + sizeof(linklistsizeint);
            this->size_data_per_element_ = this->size_links_level0_ + this->data_size_ + sizeof(labeltype) + sizeof(attributetype) * num_attribute_;
            this->offsetData_ = this->size_links_level0_;
            this->label_offset_ = this->size_links_level0_ + this->data_size_;
            attribute_offset_ = this->size_links_level0_ + this->data_size_ + sizeof(labeltype);
            this->offsetLevel0_ = 0;

            this->data_level0_memory_ = (char *)malloc(this->max_elements_ * this->size_data_per_element_);
            if (this->data_level0_memory_ == nullptr)
                throw std::runtime_error("Not enough memory");

            this->cur_element_count = 0;

            this->visited_list_pool_ = std::unique_ptr<VisitedListPool>(new VisitedListPool(1, max_elements));

            // initializations for special treatment of the first node
            this->enterpoint_node_ = -1;
            this->maxlevel_ = -1;

            this->linkLists_ = (char **)malloc(sizeof(void *) * this->max_elements_);
            if (this->linkLists_ == nullptr)
                throw std::runtime_error("Not enough memory: HierarchicalNSW failed to allocate linklists");
            this->size_links_per_element_ = this->maxM_ * sizeof(tableint) + sizeof(linklistsizeint);
            this->mult_ = 1 / log(1.0 * this->M_);
            this->revSize_ = 1.0 / this->mult_;
        }

        struct CompareByFirst_
        {
            constexpr bool operator()(std::pair<dist_t, tableint> const &a,
                                      std::pair<dist_t, tableint> const &b) const noexcept
            {
                return a.first < b.first;
            }
        };

        bool checkCondition(attributetype *attribute, const FilterConditionWithId &condition) const
        {
            attributetype value = attribute[condition.attribute_id];

            if (condition.op == "IN")
            {
                return condition.attribute_value.find(value) != condition.attribute_value.end();
            }
            else if (condition.attribute_value.size() == 1)
            {
                attributetype ref_value = *condition.attribute_value.begin();
                if (condition.op == "==")
                    return value == ref_value;
                if (condition.op == "!=")
                    return value != ref_value;
                if (condition.op == ">")
                    return value > ref_value;
                if (condition.op == "<")
                    return value < ref_value;
                if (condition.op == ">=")
                    return value >= ref_value;
                if (condition.op == "<=")
                    return value <= ref_value;
            }
            return true;
        }

        bool checkConditions(attributetype *attribute, std::vector<FilterConditionWithId> filtering_conditions) const
        {
            for (const auto &cond : filtering_conditions)
            {
                if (!checkCondition(attribute, cond))
                {
                    return false;
                }
            }
            return true;
        }

        dist_t distFilter(float p) const // selectivity
        {
            // return (1 - p) * (this->ef_ - p) * delta_d / (40 * p);
            return 2200;
        }

        bool stopSearch(dist_t candidate_dist, dist_t lowerBound, int num) const
        {
            bool flag_stop_search = candidate_dist > lowerBound;
            // bool flag_num = static_cast<size_t>(num) > this->ef_ / 2;
            // return flag_stop_search && flag_num;
            return flag_stop_search;
        }

        inline void setAttribute(tableint internal_id, attributetype *attribute) const
        {
            memcpy((this->data_level0_memory_ + internal_id * this->size_data_per_element_ + attribute_offset_), attribute, this->num_attribute_ * sizeof(attributetype));
        }

        inline attributetype *getAttributeByInternalId(tableint internal_id) const
        {
            return (attributetype *)(this->data_level0_memory_ + internal_id * this->size_data_per_element_ + attribute_offset_);
        }

        void saveIndex(const std::string &location)
        {
            std::ofstream output(location, std::ios::binary);
            std::streampos position;

            writeBinaryPOD(output, this->offsetLevel0_);
            writeBinaryPOD(output, this->max_elements_);
            writeBinaryPOD(output, this->cur_element_count);
            writeBinaryPOD(output, this->size_data_per_element_);
            writeBinaryPOD(output, this->label_offset_);
            writeBinaryPOD(output, this->offsetData_);
            writeBinaryPOD(output, this->maxlevel_);
            writeBinaryPOD(output, this->enterpoint_node_);
            writeBinaryPOD(output, this->maxM_);

            writeBinaryPOD(output, this->maxM0_);
            writeBinaryPOD(output, this->M_);
            writeBinaryPOD(output, this->mult_);
            writeBinaryPOD(output, this->ef_construction_);

            writeBinaryPOD(output, num_attribute_);
            writeBinaryPOD(output, attribute_offset_);
            writeBinaryPOD(output, delta_d);

            output.write(this->data_level0_memory_, this->cur_element_count * this->size_data_per_element_);

            for (size_t i = 0; i < this->cur_element_count; i++)
            {
                unsigned int linkListSize = this->element_levels_[i] > 0 ? this->size_links_per_element_ * this->element_levels_[i] : 0;
                writeBinaryPOD(output, linkListSize);
                if (linkListSize)
                    output.write(this->linkLists_[i], linkListSize);
            }
            output.close();
        }

        void loadIndex(const std::string &location, SpaceInterface<dist_t> *s, size_t max_elements_i = 0, size_t ef = 100)
        {
            std::ifstream input(location, std::ios::binary);

            if (!input.is_open())
                throw std::runtime_error("Cannot open file");

            this->clear();
            // get file size:
            input.seekg(0, input.end);
            std::streampos total_filesize = input.tellg();
            input.seekg(0, input.beg);

            readBinaryPOD(input, this->offsetLevel0_);
            readBinaryPOD(input, this->max_elements_);
            readBinaryPOD(input, this->cur_element_count);

            size_t max_elements = max_elements_i;
            if (max_elements < this->cur_element_count)
                max_elements = this->max_elements_;
            this->max_elements_ = max_elements;
            readBinaryPOD(input, this->size_data_per_element_);
            readBinaryPOD(input, this->label_offset_);
            readBinaryPOD(input, this->offsetData_);
            readBinaryPOD(input, this->maxlevel_);
            readBinaryPOD(input, this->enterpoint_node_);

            readBinaryPOD(input, this->maxM_);
            readBinaryPOD(input, this->maxM0_);
            readBinaryPOD(input, this->M_);
            readBinaryPOD(input, this->mult_);
            readBinaryPOD(input, this->ef_construction_);

            readBinaryPOD(input, num_attribute_);
            readBinaryPOD(input, attribute_offset_);
            readBinaryPOD(input, delta_d);

            this->data_size_ = s->get_data_size();
            this->fstdistfunc_ = s->get_dist_func();
            this->dist_func_param_ = s->get_dist_func_param();

            auto pos = input.tellg();

            /// Optional - check if index is ok:
            input.seekg(this->cur_element_count * this->size_data_per_element_, input.cur);
            for (size_t i = 0; i < this->cur_element_count; i++)
            {
                if (input.tellg() < 0 || input.tellg() >= total_filesize)
                {
                    throw std::runtime_error("Index seems to be corrupted or unsupported");
                }

                unsigned int linkListSize;
                readBinaryPOD(input, linkListSize);
                if (linkListSize != 0)
                {
                    input.seekg(linkListSize, input.cur);
                }
            }

            // throw exception if it either corrupted or old index
            if (input.tellg() != total_filesize)
                throw std::runtime_error("Index seems to be corrupted or unsupported");

            input.clear();
            /// Optional check end

            input.seekg(pos, input.beg);

            this->data_level0_memory_ = (char *)malloc(max_elements * this->size_data_per_element_);
            if (this->data_level0_memory_ == nullptr)
                throw std::runtime_error("Not enough memory: loadIndex failed to allocate level0");
            input.read(this->data_level0_memory_, this->cur_element_count * this->size_data_per_element_);

            this->size_links_per_element_ = this->maxM_ * sizeof(tableint) + sizeof(linklistsizeint);

            this->size_links_level0_ = this->maxM0_ * sizeof(tableint) + sizeof(linklistsizeint);
            std::vector<std::mutex>(max_elements).swap(this->link_list_locks_);
            std::vector<std::mutex>(MAX_LABEL_OPERATION_LOCKS).swap(this->label_op_locks_);

            this->visited_list_pool_.reset(new VisitedListPool(1, max_elements));

            this->linkLists_ = (char **)malloc(sizeof(void *) * max_elements);
            if (this->linkLists_ == nullptr)
                throw std::runtime_error("Not enough memory: loadIndex failed to allocate linklists");
            this->element_levels_ = std::vector<int>(max_elements);
            this->revSize_ = 1.0 / this->mult_;
            this->ef_ = ef;
            for (size_t i = 0; i < this->cur_element_count; i++)
            {
                this->label_lookup_[this->getExternalLabel(i)] = i;
                unsigned int linkListSize;
                readBinaryPOD(input, linkListSize);
                if (linkListSize == 0)
                {
                    this->element_levels_[i] = 0;
                    this->linkLists_[i] = nullptr;
                }
                else
                {
                    this->element_levels_[i] = linkListSize / this->size_links_per_element_;
                    this->linkLists_[i] = (char *)malloc(linkListSize);
                    if (this->linkLists_[i] == nullptr)
                        throw std::runtime_error("Not enough memory: loadIndex failed to allocate linklist");
                    input.read(this->linkLists_[i], linkListSize);
                }
            }

            input.close();

            return;
        }

        void getNeighborsByHeuristic2(
            std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst_> &top_candidates,
            const size_t M)
        {
            if (top_candidates.size() < M)
            {
                return;
            }

            std::priority_queue<std::pair<dist_t, tableint>> queue_closest;
            std::vector<std::pair<dist_t, tableint>> return_list;
            while (top_candidates.size() > 0)
            {
                queue_closest.emplace(-top_candidates.top().first, top_candidates.top().second);
                top_candidates.pop();
            }

            while (queue_closest.size())
            {
                if (return_list.size() >= M)
                    break;
                std::pair<dist_t, tableint> curent_pair = queue_closest.top();
                dist_t dist_to_query = -curent_pair.first;
                queue_closest.pop();
                bool good = true;

                for (std::pair<dist_t, tableint> second_pair : return_list)
                {
                    dist_t curdist =
                        this->fstdistfunc_(this->getDataByInternalId(second_pair.second),
                                           this->getDataByInternalId(curent_pair.second),
                                           this->dist_func_param_);
                    if (curdist < dist_to_query)
                    {
                        good = false;
                        break;
                    }
                }
                if (good)
                {
                    return_list.push_back(curent_pair);
                }
            }

            for (std::pair<dist_t, tableint> curent_pair : return_list)
            {
                top_candidates.emplace(-curent_pair.first, curent_pair.second);
            }
        }

        tableint mutuallyConnectNewElement(
            const void *data_point,
            tableint cur_c,
            std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst_> &top_candidates,
            int level)
        {
            size_t Mcurmax = level ? this->maxM_ : this->maxM0_;
            this->getNeighborsByHeuristic2(top_candidates, this->M_);
            if (top_candidates.size() > this->M_)
                throw std::runtime_error("Should be not be more than M_ candidates returned by the heuristic");

            std::vector<tableint> selectedNeighbors;
            selectedNeighbors.reserve(this->M_);
            while (top_candidates.size() > 0)
            {
                selectedNeighbors.push_back(top_candidates.top().second);
                top_candidates.pop();
            }
            tableint next_closest_entry_point = selectedNeighbors.back();

            {
                // lock only during the update
                // because during the addition the lock for cur_c is already acquired
                std::unique_lock<std::mutex> lock(this->link_list_locks_[cur_c], std::defer_lock);
                linklistsizeint *ll_cur;
                if (level == 0)
                    ll_cur = this->get_linklist0(cur_c);
                else
                    ll_cur = this->get_linklist(cur_c, level);

                if (*ll_cur)
                {
                    throw std::runtime_error("The newly inserted element should have blank link list");
                }
                this->setListCount(ll_cur, selectedNeighbors.size());
                tableint *data = (tableint *)(ll_cur + 1);
                for (size_t idx = 0; idx < selectedNeighbors.size(); idx++)
                {
                    if (data[idx])
                        throw std::runtime_error("Possible memory corruption");
                    if (level > this->element_levels_[selectedNeighbors[idx]])
                        throw std::runtime_error("Trying to make a link on a non-existent level");

                    data[idx] = selectedNeighbors[idx];
                }
            }
            for (size_t idx = 0; idx < selectedNeighbors.size(); idx++)
            {
                std::unique_lock<std::mutex> lock(this->link_list_locks_[selectedNeighbors[idx]]);

                linklistsizeint *ll_other;
                if (level == 0)
                    ll_other = this->get_linklist0(selectedNeighbors[idx]);
                else
                    ll_other = this->get_linklist(selectedNeighbors[idx], level);

                size_t sz_link_list_other = this->getListCount(ll_other);

                if (sz_link_list_other > Mcurmax)
                    throw std::runtime_error("Bad value of sz_link_list_other");
                if (selectedNeighbors[idx] == cur_c)
                    throw std::runtime_error("Trying to connect an element to itself");
                if (level > this->element_levels_[selectedNeighbors[idx]])
                    throw std::runtime_error("Trying to make a link on a non-existent level");

                tableint *data = (tableint *)(ll_other + 1);

                bool is_cur_c_present = false;

                // If cur_c is already present in the neighboring connections of `selectedNeighbors[idx]` then no need to modify any connections or run the heuristics.
                if (!is_cur_c_present)
                {
                    if (sz_link_list_other < Mcurmax)
                    {
                        data[sz_link_list_other] = cur_c;
                        this->setListCount(ll_other, sz_link_list_other + 1);
                    }
                    else
                    {
                        // finding the "weakest" element to replace it with the new one
                        dist_t d_max = this->fstdistfunc_(this->getDataByInternalId(cur_c), this->getDataByInternalId(selectedNeighbors[idx]),
                                                          this->dist_func_param_);
                        // Heuristic:
                        std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst_> candidates;
                        candidates.emplace(d_max, cur_c);

                        for (size_t j = 0; j < sz_link_list_other; j++)
                        {
                            candidates.emplace(
                                this->fstdistfunc_(this->getDataByInternalId(data[j]), this->getDataByInternalId(selectedNeighbors[idx]),
                                                   this->dist_func_param_),
                                data[j]);
                        }

                        this->getNeighborsByHeuristic2(candidates, Mcurmax);

                        int indx = 0;
                        while (candidates.size() > 0)
                        {
                            data[indx] = candidates.top().second;
                            candidates.pop();
                            indx++;
                        }

                        this->setListCount(ll_other, indx);
                        // Nearest K:
                        /*int indx = -1;
                        for (int j = 0; j < sz_link_list_other; j++) {
                            dist_t d = fstdistfunc_(getDataByInternalId(data[j]), getDataByInternalId(rez[idx]), dist_func_param_);
                            if (d > d_max) {
                                indx = j;
                                d_max = d;
                            }
                        }
                        if (indx >= 0) {
                            data[indx] = cur_c;
                        } */
                    }
                }
            }

            return next_closest_entry_point;
        }

        std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst_>
        searchBaseLayerFilter(tableint ep_id, const void *data_point, int layer)
        {
            VisitedList *vl = this->visited_list_pool_->getFreeVisitedList();
            vl_type *visited_array = vl->mass;
            vl_type visited_array_tag = vl->curV;

            std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst_> top_candidates;
            std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst_> candidateSet;

            dist_t lowerBound;
            if (!this->isMarkedDeleted(ep_id))
            {
                dist_t dist = this->fstdistfunc_(data_point, this->getDataByInternalId(ep_id), this->dist_func_param_);
                top_candidates.emplace(dist, ep_id);
                lowerBound = dist;
                candidateSet.emplace(-dist, ep_id);
            }
            else
            {
                lowerBound = std::numeric_limits<dist_t>::max();
                candidateSet.emplace(-lowerBound, ep_id);
            }
            visited_array[ep_id] = visited_array_tag;

            while (!candidateSet.empty())
            {
                std::pair<dist_t, tableint> curr_el_pair = candidateSet.top();
                if ((-curr_el_pair.first) > lowerBound && top_candidates.size() == this->ef_construction_)
                {
                    break;
                }
                candidateSet.pop();

                tableint curNodeNum = curr_el_pair.second;

                std::unique_lock<std::mutex> lock(this->link_list_locks_[curNodeNum]);

                int *data; // = (int *)(linkList0_ + curNodeNum * size_links_per_element0_);
                if (layer == 0)
                {
                    data = (int *)this->get_linklist0(curNodeNum);
                }
                else
                {
                    data = (int *)this->get_linklist(curNodeNum, layer);
                    //                    data = (int *) (linkLists_[curNodeNum] + (layer - 1) * size_links_per_element_);
                }
                size_t size = this->getListCount((linklistsizeint *)data);
                tableint *datal = (tableint *)(data + 1);
#ifdef USE_SSE
                _mm_prefetch((char *)(visited_array + *(data + 1)), _MM_HINT_T0);
                _mm_prefetch((char *)(visited_array + *(data + 1) + 64), _MM_HINT_T0);
                _mm_prefetch(this->getDataByInternalId(*datal), _MM_HINT_T0);
                _mm_prefetch(this->getDataByInternalId(*(datal + 1)), _MM_HINT_T0);
#endif

                for (size_t j = 0; j < size; j++)
                {
                    tableint candidate_id = *(datal + j);
//                    if (candidate_id == 0) continue;
#ifdef USE_SSE
                    _mm_prefetch((char *)(visited_array + *(datal + j + 1)), _MM_HINT_T0);
                    _mm_prefetch(this->getDataByInternalId(*(datal + j + 1)), _MM_HINT_T0);
#endif
                    if (visited_array[candidate_id] == visited_array_tag)
                        continue;
                    visited_array[candidate_id] = visited_array_tag;
                    char *currObj1 = (this->getDataByInternalId(candidate_id));

                    dist_t dist1 = this->fstdistfunc_(data_point, currObj1, this->dist_func_param_);
                    if (top_candidates.size() < this->ef_construction_ || lowerBound > dist1)
                    {
                        candidateSet.emplace(-dist1, candidate_id);
#ifdef USE_SSE
                        _mm_prefetch(this->getDataByInternalId(candidateSet.top().second), _MM_HINT_T0);
#endif

                        if (!this->isMarkedDeleted(candidate_id))
                            top_candidates.emplace(dist1, candidate_id);

                        if (top_candidates.size() > this->ef_construction_)
                            top_candidates.pop();

                        if (!top_candidates.empty())
                            lowerBound = top_candidates.top().first;
                    }
                }
            }
            this->visited_list_pool_->releaseVisitedList(vl);

            return top_candidates;
        }

        std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst_>
        searchBaseLayerSTFilter(tableint ep_id,
                                const void *data_point,
                                size_t ef,
                                dist_t e_distance,
                                std::vector<FilterConditionWithId> filtering_conditions) const
        {
            VisitedList *vl = this->visited_list_pool_->getFreeVisitedList();
            vl_type *visited_array = vl->mass;
            vl_type visited_array_tag = vl->curV;

            std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst_> top_candidates;
            std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst_> candidate_set;

            dist_t lowerBound;
            // get the entry point for the top
            char *ep_data = this->getDataByInternalId(ep_id);
            // calculate the distance between entry and query.
            dist_t dist;
            if (checkConditions(getAttributeByInternalId(ep_id), filtering_conditions))
                dist = this->fstdistfunc_(data_point, ep_data, this->dist_func_param_);
            else
                dist = this->fstdistfunc_(data_point, ep_data, this->dist_func_param_) + e_distance;
            lowerBound = dist;
            top_candidates.emplace(dist, ep_id);
            candidate_set.emplace(-dist, ep_id);

            visited_array[ep_id] = visited_array_tag;

            int num_in_range = 0;
            // int length = 0;

            while (!candidate_set.empty())
            {
                // candidate_set.top() is the closest point
                std::pair<dist_t, tableint> current_node_pair = candidate_set.top();
                dist_t candidate_dist = -current_node_pair.first;

                // if part of vectors in candidate are in range, stop searching
                if (stopSearch(candidate_dist, lowerBound, num_in_range))
                {
                    break;
                }

                // pop the closest point
                candidate_set.pop();

                tableint current_node_id = current_node_pair.second;

                // get the closest neighbor
                int *data = (int *)this->get_linklist0(current_node_id);
                size_t size = this->getListCount((linklistsizeint *)data);

#ifdef USE_SSE
                _mm_prefetch((char *)(visited_array + *(data + 1)), _MM_HINT_T0);
                _mm_prefetch((char *)(visited_array + *(data + 1) + 64), _MM_HINT_T0);
                _mm_prefetch(data_level0_memory_ + (*(data + 1)) * size_data_per_element_ + offsetData_, _MM_HINT_T0);
                _mm_prefetch((char *)(data + 2), _MM_HINT_T0);
#endif

                // visit all the neighbors
                for (size_t j = 1; j <= size; j++)
                {
                    int candidate_id = *(data + j);
#ifdef USE_SSE
                    _mm_prefetch((char *)(visited_array + *(data + j + 1)), _MM_HINT_T0);
                    _mm_prefetch(data_level0_memory_ + (*(data + j + 1)) * size_data_per_element_ + offsetData_,
                                 _MM_HINT_T0); ////////////
#endif
                    // check if the point is visited
                    if (!(visited_array[candidate_id] == visited_array_tag))
                    {
                        visited_array[candidate_id] = visited_array_tag; // mark the point which is visited

                        char *currObj1 = (this->getDataByInternalId(candidate_id));
                        dist_t dist1;
                        if (checkConditions(getAttributeByInternalId(candidate_id), filtering_conditions))
                            dist1 = this->fstdistfunc_(data_point, currObj1, this->dist_func_param_);
                        else
                            dist1 = this->fstdistfunc_(data_point, currObj1, this->dist_func_param_) + e_distance;
                        bool flag_consider_candidate;
                        flag_consider_candidate = top_candidates.size() < ef || lowerBound > dist1;

                        if (flag_consider_candidate)
                        {
                            candidate_set.emplace(-dist1, candidate_id);
                            if (checkConditions(getAttributeByInternalId(candidate_id), filtering_conditions))
                                num_in_range++;
#ifdef USE_SSE
                            _mm_prefetch(data_level0_memory_ + candidate_set.top().second * size_data_per_element_ +
                                             offsetLevel0_, ///////////
                                         _MM_HINT_T0);      ////////////////////////
#endif

                            // if (checkConditions(getAttributeByInternalId(candidate_id), filtering_conditions))
                            top_candidates.emplace(dist1, candidate_id);

                            bool flag_remove_extra = false;

                            flag_remove_extra = top_candidates.size() > ef;

                            while (flag_remove_extra)
                            {
                                top_candidates.pop();
                                flag_remove_extra = top_candidates.size() > ef;
                                if (checkConditions(getAttributeByInternalId(candidate_id), filtering_conditions))
                                    num_in_range--;
                            }

                            if (!top_candidates.empty())
                                lowerBound = top_candidates.top().first;
                        }
                    }
                }
            }

            this->visited_list_pool_->releaseVisitedList(vl);
            return top_candidates;
        }

        void addPoint(const void *data_point, labeltype label, attributetype *attribute, bool replace_deleted = false)
        {
            if ((this->allow_replace_deleted_ == false) && (replace_deleted == true))
            {
                throw std::runtime_error("Replacement of deleted elements is disabled in constructor");
            }

            // lock all operations with element by label
            std::unique_lock<std::mutex> lock_label(this->getLabelOpMutex(label));
            if (!replace_deleted)
            {
                addPoint(data_point, label, -1, attribute);
                return;
            }
            // check if there is vacant place
            tableint internal_id_replaced;
            std::unique_lock<std::mutex> lock_deleted_elements(this->deleted_elements_lock);
            bool is_vacant_place = !this->deleted_elements.empty();
            if (is_vacant_place)
            {
                internal_id_replaced = *this->deleted_elements.begin();
                this->deleted_elements.erase(internal_id_replaced);
            }
            lock_deleted_elements.unlock();

            // if there is no vacant place then add or update point
            // else add point to vacant place
            if (!is_vacant_place)
            {
                addPoint(data_point, label, -1, attribute);
            }
            else
            {
                // we assume that there are no concurrent operations on deleted element
                labeltype label_replaced = this->getExternalLabel(internal_id_replaced);
                this->setExternalLabel(internal_id_replaced, label);
                setAttribute(internal_id_replaced, attribute);

                std::unique_lock<std::mutex> lock_table(this->label_lookup_lock);
                this->label_lookup_.erase(label_replaced);
                this->label_lookup_[label] = internal_id_replaced;
                lock_table.unlock();

                this->unmarkDeletedInternal(internal_id_replaced);
                this->updatePoint(data_point, internal_id_replaced, 1.0);
            }
        }

        tableint addPoint(const void *data_point, labeltype label, int level, attributetype *attribute)
        {
            tableint cur_c = 0;
            {
                std::unique_lock<std::mutex> lock_table(this->label_lookup_lock);
                auto search = this->label_lookup_.find(label);
                if (search != this->label_lookup_.end())
                {
                    tableint existingInternalId = search->second;
                    if (this->allow_replace_deleted_)
                    {
                        if (this->isMarkedDeleted(existingInternalId))
                        {
                            throw std::runtime_error("Can't use addPoint to update deleted elements if replacement of deleted elements is enabled.");
                        }
                    }
                    lock_table.unlock();

                    if (this->isMarkedDeleted(existingInternalId))
                    {
                        this->unmarkDeletedInternal(existingInternalId);
                    }
                    this->updatePoint(data_point, existingInternalId, 1.0);

                    return existingInternalId;
                }

                if (this->cur_element_count >= this->max_elements_)
                {
                    throw std::runtime_error("The number of elements exceeds the specified limit");
                }

                cur_c = this->cur_element_count;
                this->cur_element_count++;
                this->label_lookup_[label] = cur_c;
            }

            std::unique_lock<std::mutex> lock_el(this->link_list_locks_[cur_c]);
            int curlevel = this->getRandomLevel(this->mult_);
            if (level > 0)
                curlevel = level;

            this->element_levels_[cur_c] = curlevel;

            std::unique_lock<std::mutex> templock(this->global);
            int maxlevelcopy = this->maxlevel_;
            if (curlevel <= maxlevelcopy)
                templock.unlock();
            tableint currObj = this->enterpoint_node_;

            memset(this->data_level0_memory_ + cur_c * this->size_data_per_element_ + this->offsetLevel0_, 0, this->size_data_per_element_);

            // Initialisation of the data, label and attribute
            memcpy(this->getExternalLabeLp(cur_c), &label, sizeof(labeltype));
            memcpy(this->getDataByInternalId(cur_c), data_point, this->data_size_);
            setAttribute(cur_c, attribute);

            if (curlevel)
            {
                this->linkLists_[cur_c] = (char *)malloc(this->size_links_per_element_ * curlevel + 1);
                if (this->linkLists_[cur_c] == nullptr)
                    throw std::runtime_error("Not enough memory: addPoint failed to allocate linklist");
                memset(this->linkLists_[cur_c], 0, this->size_links_per_element_ * curlevel + 1);
            }

            dist_t local_delta_d = 0;

            if ((signed)currObj != -1)
            {
                if (curlevel < maxlevelcopy)
                {
                    dist_t curdist = this->fstdistfunc_(data_point, this->getDataByInternalId(currObj), this->dist_func_param_);
                    for (int level = maxlevelcopy; level > curlevel; level--)
                    {
                        bool changed = true;
                        while (changed)
                        {
                            changed = false;
                            unsigned int *data;
                            std::unique_lock<std::mutex> lock(this->link_list_locks_[currObj]);
                            data = this->get_linklist(currObj, level);
                            int size = this->getListCount(data);

                            tableint *datal = (tableint *)(data + 1);
                            for (int i = 0; i < size; i++)
                            {
                                tableint cand = datal[i];
                                if (cand < 0 || cand > this->max_elements_)
                                    throw std::runtime_error("cand error");
                                dist_t d = this->fstdistfunc_(data_point, this->getDataByInternalId(cand), this->dist_func_param_);
                                if (d < curdist)
                                {
                                    curdist = d;
                                    currObj = cand;
                                    changed = true;
                                }
                            }
                        }
                    }
                }

                for (int level = std::min(curlevel, maxlevelcopy); level >= 0; level--)
                {
                    if (level > maxlevelcopy || level < 0) // possible?
                        throw std::runtime_error("Level error");
                    std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst_> top_candidates = searchBaseLayerFilter(
                        currObj, data_point, level);

                    // delta_d
                    if (level == 0 && top_candidates.size() == this->ef_construction_)
                    {
                        dist_t rate = this->ef_construction_ - 10;
                        auto temp = top_candidates;
                        std::vector<std::pair<dist_t, tableint>> elements;
                        while (!temp.empty())
                        {
                            elements.push_back(temp.top());
                            temp.pop();
                        }

                        dist_t diff = elements[0].first - elements[elements.size() - 10].first;
                        local_delta_d += diff / (rate * this->max_elements_);
                    }
                    currObj = this->mutuallyConnectNewElement(data_point, cur_c, top_candidates, level);
                }
            }
            else
            {
                // Do nothing for the first element
                this->enterpoint_node_ = 0;
                this->maxlevel_ = curlevel;
            }

            std::lock_guard<std::mutex> lock(delta_mutex);
            delta_d += local_delta_d;

            // Releasing lock for the maximum level
            if (curlevel > maxlevelcopy)
            {
                this->enterpoint_node_ = cur_c;
                this->maxlevel_ = curlevel;
            }
            return cur_c;
        }

        std::priority_queue<std::pair<dist_t, labeltype>>
        searchGraph(const void *query_data, size_t k, float p, std::vector<FilterConditionWithId> filtering_conditions) const
        {
            std::priority_queue<std::pair<dist_t, labeltype>> result;
            if (this->cur_element_count == 0)
                return result;

            dist_t e_distance = distFilter(p);

            tableint currObj = this->enterpoint_node_;
            // calculate the distance between query and enterpoint
            dist_t curdist = this->fstdistfunc_(query_data, this->getDataByInternalId(this->enterpoint_node_), this->dist_func_param_);
            // visit from top to bottom, get the closest point currObj
            for (int level = this->maxlevel_; level > 0; level--)
            {
                bool changed = true;
                while (changed)
                {
                    changed = false;
                    unsigned int *data;

                    // get the neighbor list of currObj at level
                    data = (unsigned int *)this->get_linklist(currObj, level);
                    int size = this->getListCount(data);

                    tableint *datal = (tableint *)(data + 1);
                    // visit the neighbors of currObj
                    for (int i = 0; i < size; i++)
                    {
                        tableint cand = datal[i];
                        if (cand < 0 || cand > this->max_elements_)
                            throw std::runtime_error("cand error");
                        dist_t d = this->fstdistfunc_(query_data, this->getDataByInternalId(cand), this->dist_func_param_);
                        // if the distance is lower the curdist, continue visiting.
                        if (d < curdist)
                        {
                            curdist = d;
                            currObj = cand;
                            changed = true;
                        }
                    }
                }
            }

            std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst_> top_candidates;
            top_candidates = searchBaseLayerSTFilter(currObj, query_data, std::max(this->ef_, k), e_distance, filtering_conditions);
            while (top_candidates.size() > 0)
            {
                std::pair<dist_t, tableint> rez = top_candidates.top();
                attributetype *attribute = getAttributeByInternalId(rez.second);
                if (checkConditions(attribute, filtering_conditions))
                {
                    result.push(std::pair<dist_t, labeltype>(rez.first, this->getExternalLabel(rez.second)));
                }

                top_candidates.pop();
            }
            while (result.size() < k)
                result.push(std::pair<dist_t, labeltype>(1000000, -1));
            while (result.size() > k)
            {
                result.pop();
            }

            return result;
        }

        std::priority_queue<std::pair<dist_t, labeltype>>
        searchBruteForce(const void *query_data, size_t k, std::vector<FilterConditionWithId> filtering_conditions) const
        {
            std::priority_queue<std::pair<dist_t, labeltype>> result;
            for (int i = 0; i < this->max_elements_; i++)
            {
                if (checkConditions(getAttributeByInternalId(i), filtering_conditions))
                {
                    dist_t d = this->fstdistfunc_(query_data, this->getDataByInternalId(i), this->dist_func_param_);
                    result.push(std::pair<dist_t, labeltype>(d, this->getExternalLabel(i)));
                }
            }
            while (result.size() < k)
                result.push(std::pair<dist_t, labeltype>(1000000, -1));
            while (result.size() > k)
            {
                result.pop();
            }
            return result;
        }

        float selectivityEstimator(std::vector<FilterConditionWithId> filtering_conditions) const
        {
            float p;
            size_t count = 0;
            size_t sample_size = std::max(this->max_elements_/10000, static_cast<size_t>(1));
            for(size_t i = 0; i< this->max_elements_; i+=10000)
            {
                if (checkConditions(getAttributeByInternalId(i), filtering_conditions))
                count++;
            }
            p = (float)count / sample_size;
            return p;
        }

        std::priority_queue<std::pair<dist_t, labeltype>>
        searchKnn(const void *query_data, size_t k, std::vector<FilterConditionWithId> filtering_conditions) const
        {
            float p = selectivityEstimator(filtering_conditions);
            if (p > 0.01)
                return searchGraph(query_data, k, p, filtering_conditions);
            else
                return searchBruteForce(query_data, k, filtering_conditions);
        }

        ~FAVOR()
        {
            this->clear();
        }
    };
}