#include "dataset.h"
#include "favor.h"
#include <thread>

template<class Function>
inline void ParallelFor(size_t start, size_t end, size_t numThreads, Function fn) {
    if (numThreads <= 0) {
        numThreads = std::thread::hardware_concurrency();
    }

    if (numThreads == 1) {
        for (size_t id = start; id < end; id++) {
            fn(id, 0);
        }
    } else {
        std::vector<std::thread> threads;
        std::atomic<size_t> current(start);

        std::exception_ptr lastException = nullptr;
        std::mutex lastExceptMutex;

        for (size_t threadId = 0; threadId < numThreads; ++threadId) {
            threads.push_back(std::thread([&, threadId] {
                while (true) {
                    size_t id = current.fetch_add(1);

                    if (id >= end) {
                        break;
                    }

                    try {
                        fn(id, threadId);
                    } catch (...) {
                        std::unique_lock<std::mutex> lastExcepLock(lastExceptMutex);
                        lastException = std::current_exception();
                        current = end;
                        break;
                    }
                }
            }));
        }
        for (auto &thread : threads) {
            thread.join();
        }
        if (lastException) {
            std::rethrow_exception(lastException);
        }
    }
}

int main(int argc, char* argv[]) {
    int num_threads = 32;
    
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " baseset_path" << " attribute_path" << " index_path\n";
        return 1;
    }

    std::string baseset_path(argv[1]);
    std::string attribute_path(argv[2]);
    std::string index_path(argv[3]);

    BaseSet baseset;
    baseset.read_data(baseset_path);
    baseset.get_attribute(attribute_path);

    int dim = baseset.dim;
    int num = baseset.num;
    int M = 64;
    int ef_construction = 200;

    hnswlib::L2Space space(dim);
    favor::FAVOR<float> *alg_hnsw = new favor::FAVOR<float>(&space, num, M, ef_construction, baseset.attribute_num);

    std::cout << "begin building graph" << std::endl;
    ParallelFor(0, num, num_threads, [&](size_t i, size_t threadId) {
        alg_hnsw->addPoint(baseset.vectors.at(i).data(), baseset.vector_id[i], baseset.attribute + i * baseset.attribute_num);
    });
    std::cout << "finish building graph" << std::endl;
    float delta_d = alg_hnsw->delta_d;
    std::cout << "delta d = " << delta_d << std::endl;

    alg_hnsw->saveIndex(index_path);

    std::cout << "save index in " << index_path << std::endl;

    delete alg_hnsw;
    return 0;
}