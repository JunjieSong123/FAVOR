# FAVOR
This repository contains the official implementation of FAVOR, a novel vector search system for hybrid queries that combine approximate nearest neighbor search (ANNS) with complex attribute filtering.

FAVOR efficiently handles arbitrary filtering conditions and maintains stable, high performance across different selectivity levels. It features an integrated architecture, a new HNSW-based search algorithm with an exclusion mechanism, and a dynamic search selector.

## Compilation
```bash
mkdir build && cd build
cmake ..
make
cd ..
```

## Get Attributes
FAVOR uses a simple TXT file format for storing vector attributes.
```bash
./build/app/generate_attribute {baseset_path} {attribute_path}
```

## Ground Truth
```bash
./build/app/generate_groundtruth {baseset_path} {queryset_path} {attribute_path} {topk} {groundtruth_path} {filtering_conditions}
```

## Build Index
```bash
./build/app/build_index {baseset_path} {attribute_path} {index_path}
```

## Search
```
./build/app/search {baseset_path} {queryset_path} {topk} {groundtruth_path} {filtering_conditions} {index_path} {ef} 
```