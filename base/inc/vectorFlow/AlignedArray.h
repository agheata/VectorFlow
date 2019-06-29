#ifndef VECTORFLOW_ALIGNEDARRAY_H
#define VECTORFLOW_ALIGNEDARRAY_H

#include <cstdlib>
#include <VecCore/VecCore>

namespace vectorflow {

/**
 * The AlignedArray class will serve as a container for an
 * aligned array to handle any arbitrary basic data types.
 * The container receives 3 templated parameters:
 * + Data: Type of data to be stored in the array.
 * + Size: Size of the array.
 * + Alignment: SIMD alignment size.
 *
 * Note: 1st prototype, it can change to include new
 *       features to make it a resizable array, etc.
 */

template <typename Data, std::size_t Size, std::size_t Alignment>
class AlignedArray {
  private:
    Data* fArray;

  public:
    AlignedArray() {
      fArray = (Data*) vecCore::AlignedAlloc(Alignment, Size * sizeof(Data));
    }

    virtual ~AlignedArray() {
      vecCore::AlignedFree(fArray);
    }

    Data& operator[](int index) {
      if (!(index >= 0 && index < Size)) {
        printf("\nERROR: Segmentation fault in AlignedArray, index = %d", index);
        exit(-1);
      }
      return fArray[index];
    }
};

} // namespace vectorflow

#endif
