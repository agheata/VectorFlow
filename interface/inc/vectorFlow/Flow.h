#ifndef VECTORFLOW_FLOW_H
#define VECTORFLOW_FLOW_H

#include <vectorFlow/Work.h>

namespace vectorflow {

/// Class providing common functionality for different types of worflows.
/**
  A flow represents a sequence of work tasks that can be executed in scalar or vector mode.
The class does not provide execution API, one should use explicitely derived classes for that.
*/
template <typename Data, typename DataContainer, std::size_t NSeq>
class Flow {
public:
  using size_t      = std::size_t;
  using Work_t      = Work<Data, DataContainer>;
  using WorkSeq_t   = std::array<Work_t*, NSeq>;
  using BasketSeq_t = std::array<DataContainer, NSeq>;

protected:
  WorkSeq_t fWorkSeq;               ///< Work sequence
  bool fVectorMode[NSeq] = {false}; ///< Stages to be executed in vector mode

public:
  virtual ~Flow() {}

  /// Add a task for a given stage in the pipeline
  void    AddWork(Work_t *task, int stage) { fWorkSeq[stage] = task; }

  /// Mark a stage to be executed in vector mode
  void    SetVectorMode(size_t stage, bool flag) { fVectorMode[stage] = flag; }
  bool    GetVectorMode(size_t stage) const { return fVectorMode[stage]; }

  /// Getter for number of tasks
  size_t  GetNtasks() const { return NSeq; }

  /// Getter for the work at a given stage
  Work_t &GetWork(size_t stage) { return fWorkSeq[stage]; }

};

} // namespace vectorflow

#endif
