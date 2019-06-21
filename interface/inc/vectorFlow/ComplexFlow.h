#ifndef VECTORFLOW_COMPLEXFLOW_H
#define VECTORFLOW_COMPLEXFLOW_H

#include <vector>
#include <array>
#include <vectorFlow/Flow.h>

namespace vectorflow {

/// Class representing a complex workflow where work tasks can dispatch data from one to another.
/**
The complex workflow allows registering a set of work tasks, creating a buffer (basket) of input states
for each of them. Each task can execute in scalar or vector mode only on the data available on its input basket.
After execution, the data will be dispatched to another task basket, not necessary the one following in the
sequence of work.
*/ 
template <typename Data, typename DataContainer, std::size_t NSeq>
class ComplexFlow : public Flow<Data, DataContainer, NSeq>
{
  using Base_t = Flow<Data, DataContainer, NSeq>;
  using typename Base_t::size_t;
  using typename Base_t::Work_t;
  using typename Base_t::WorkSeq_t;

  using typename Base_t::BasketSeq_t;

private:
  BasketSeq_t fBaskets;             ///< Sequence of baskets for each stage in the workflow
  DataContainer fOutput;            ///< Output container for the complex flow


public:
  virtual ~ComplexFlow() {}

  /// Add data the workflow
  void AddData(Data *data) { fBaskets[0].push_back(data); }
  void AddData(DataContainer const &vdata)
  {
    std::copy(vdata.begin(), vdata.end(), std::back_inserter(fBaskets[0]));
  }

  void Clear(size_t stage) { fBaskets[stage].clear(); }

  /// Getters for input data
  size_t GetNinput(size_t stage) const { return fBaskets[stage].size(); }
  size_t GetNstates() const
  { 
    size_t nstates = 0;
    for (auto &basket : fBaskets)
      nstates += basket.size();
    return nstates;
  }

  DataContainer &InputData(size_t stage) { return fBaskets[stage]; }

  /// Execute a given stage in scalar mode
  void Execute_s(size_t stage) { Base_t::fWorkSeq[stage]->ExecuteLoop(fBaskets[stage]); fBaskets[stage].clear(); }

  /// Execute the full sequence once in scalar mode
  void Execute_s() { for (size_t stage = 0; stage < NSeq; ++stage) Execute_s(stage); }

  /// Execute a given stage in vector mode
  void Execute_v(size_t stage) { Base_t::fWorkSeq[stage]->Execute(fBaskets[stage]); fBaskets[stage].clear(); }
  
  /// Execute all stages in vector mode
  void Execute_v() { for (size_t stage = 0; stage < NSeq; ++stage) Execute_v(stage); }

  /// Execute the full pipeline according the vector mode flags
  void Execute()
  {
    for (size_t stage = 0; stage < NSeq; ++stage) {
      if (Base_t::fVectorMode[stage]) Execute_v(stage);
      else                            Execute_s(stage);
    }
  }

};

} // namespace vectorflow

#endif

