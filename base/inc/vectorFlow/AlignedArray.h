#ifndef VECTORFLOW_ALIGNEDARRAY_H
#define VECTORFLOW_ALIGNEDARRAY_H

#include <VecCore/VecCore>
#include <vectorFlow/VectorTypes.h>

namespace vectorflow {

/**
 * The AlignedArray class will serve as a container for an
 * aligned array to handle any arbitrary basic data types.
 *
 * The container receives a templated parameters:
 * + Data: Type of data to be stored in the array.
 * 
 * The constructor is defined with the size of the array:
 * + fSize: Size of the array.
 *
 * The aligment (VectorTypes.h) is defined by the architecture:
 * + kVecAlignment: SIMD alignment size.
 *
 * Note: 1st prototype, it can change to include new
 *       features to make it a resizable array, etc.
 */

template <typename Data>
class AlignedArray {
  private:
    Data* fArray;
    std::size_t fSize;

  public:
    AlignedArray(std::size_t size) {
      fArray = (Data*) vecCore::AlignedAlloc(kVecAlignment, size * sizeof(Data));
      fSize  = size;
    }

    virtual ~AlignedArray() {
      vecCore::AlignedFree(fArray);
    }

    Data& operator[](int index) {
      if (!(index >= 0 && index < fSize)) {
        printf("\nERROR: Segmentation fault in AlignedArray, index = %d", index);
        exit(-1);
      }
      return fArray[index];
    }

    Data* Head() { return fArray; }

    std::size_t Size() { return fSize; }
};

} // namespace vectorflow

#endif
