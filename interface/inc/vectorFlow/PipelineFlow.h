#ifndef VECTORFLOW_PIPELINEFLOW_H
#define VECTORFLOW_PIPELINEFLOW_H

#include <vector>
#include <array>
#include <vectorFlow/Work.h>

namespace vectorflow {

/// Class representing a pipeline flow of tasks:
/**
   The pipeline flow is a sequence of work tasks: A -> B -> C -> ...
   The tasks have to derive from the Work class). All tasks will run on the same input container of data.
   The data has to be added to the input basket before executing tasks in the sequence.
   
   The execution of a given work task can be done in a loop over the input data, calling the scalar interface
   for executing on each data element, or on the full container, passing to the task vector interface.

   One can trigger manual execution in scalar or vector mode for individual tasks in the pipeline, or can execute
   the full pipeline according a set of switches defined by the user (scalar or vector mode per stage)

   Connection between 2 pipelines has to be done manually.

   DataContainer is a container holding Data* as elements, and should be dynamically exetensible (e.g. std::vector
   but not std::array)
*/
template <typename Data, typename DataContainer, std::size_t NSeq>
class PipelineFlow() : public Flow<Data, DataContainer, NSeq>
{
  using size_t = std::size_t;
  using Work_t = Work<Data, DataContainer>;
  using WorkSeq_t = std::array<Work, NSeq>;

private:
  DataContainer fBasket;               ///< The container for input data

public:
  virtual ~PipelineFlow() {}

  /// Add data to a given stage
  void AddData(Data *data) { fBaskets.push_back(data); }
  void AddData(DataContainer const &vdata)
  {
    std::copy(vdata.begin(), vdata.end(), std::back_inserter(fBasket));
  }

  /// Clear the input data
  void Clear() { fBasket.clear(); }

  /// Getters for input data
  size_t GetNinput() const { return fBasket.size(); }
  DataContainer &InputData() { return fBasket; }

  /// Execute a given stage in scalar mode
  void Execute_s(size_t stage) { fWorkSeq[stage].ExecuteLoop(fBasket); }

  /// Execute the full pipeline in scalar mode
  void Execute_s() { for (size_t stage = 0; stage < NSeq; ++stage) Execute_s(stage); }

  /// Execute a given stage in vector mode
  void Execute_v(size_t stage) { fWorkSeq[stage].Execute(fBasket); }
  
  /// Execute all stages in vector mode
  void Execute_v() { for (size_t stage = 0; stage < NSeq; ++stage) Execute_v(stage); }

  /// Execute the full pipeline according the vector mode flags
  void Execute()
  {
    for (size_t stage = 0; stage < NSeq; ++stage) {
      if (fVectorMode[stage]) Execute_v(stage);
      else                    Execute_s(stage);
    }
  }

};

} // namespace vectorflow

#endif

